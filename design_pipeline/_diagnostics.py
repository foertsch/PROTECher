"""Diagnostics: the checks that need no answer key, so they run on the target.

Reached as `design_pipeline.metrics.diagnostics` - the file is private because
that is the only name it should be imported under, and `metrics.__getattr__`
loads it on first use so that `import design_pipeline` still costs nothing.

**None of these measure whether a prediction is correct.** Every function here
is computable from a campaign's own output, which is exactly why none of them
can speak to accuracy: there is nothing in a run that knows the answer. They
detect three other things, all of which are worth catching and none of which
is quality:

* **Bugs** - a chain assignment that produced a four-residue interface, a sign
  convention flipped so that the worst mutations sort to the top, a scorer
  whose argmax is never the wild type on a matured antibody.
* **Redundancy** - two score columns correlating at 0.8 are one opinion with
  two names, and an ensemble built from them is thinner than its column count
  suggests.
* **Known-bad regimes** - an ensemble that split across two epitopes, scores
  that move when the seed moves, a shortlist that is a hundred bets on three
  positions.

Only the Category 1 metrics in `design_pipeline.metrics`, run in the harness
or the validation notebook against measured affinities or deposited
structures, say anything about accuracy. Any report built out of this module
carries `CAVEAT` verbatim, and there is deliberately no function here that
blends a diagnostic with a Category 1 metric into a single score: such a
number would be read as "how good are the predictions", which is the one
question these cannot answer.

The functions are grouped by the pipeline step they diagnose, because that is
the order in which things go wrong and the order in which you look:

1. structure - `seed_agreement`, `decoy_gap`
2. enumerate - `variant_count`, `position_buckets`, `interface_size`
3. score - `inter_model_spearman`, `beneficial_fraction`, `wt_argmax_rate`,
   `positional_entropy`, `seed_rank_stability`
4. gate - `pass_rate`, `gate_score_correlation`
5. pairs - `epistasis`, `epistasis_magnitude`, `distance_vs_epistasis`
6. shortlist - `position_diversity`, `chain_balance`,
   `edit_distance_histogram`, `final_rank_agreement`

Every threshold here is a named constant and a keyword argument, and every one
of them is a loose first guess. They are written down so that a report can say
"redundant" instead of "0.81", not because the numbers are settled; the
harness settles numbers with data, and these have not been through it.

Anything that takes a table takes either a pandas frame or one of this
package's own tables (`VariantTable`, `ScoreTable`) and reads `.frame` off it.
"""

from __future__ import annotations

import math
from pathlib import Path
from typing import TYPE_CHECKING, Any, Iterable, Mapping, Sequence

import numpy as np

from .metrics import (
    INTERFACE_CUTOFF,
    Scalar,
    interface_residues,
    interface_rmsd,
    nanmean,
    score_columns,
    scorer_name,
    spearman,
    split_chains,
)
from .structures import ca_coordinates, chain_sequence
from .structures import load as load_structure

if TYPE_CHECKING:
    import pandas as pd

CAVEAT = (
    "These are diagnostics, not quality measures: they are computed from the "
    "run's own output and so cannot say whether any prediction is correct. "
    "Only the ground-truth metrics - run in the harness against measured "
    "affinities, or against a deposited structure - speak to accuracy."
)
"""Goes verbatim into any report these feed. The failure mode this guards
against is a reader seeing a table of green numbers and concluding the
predictions are good, when every number in it was computed without reference
to a single measurement."""


# --------------------------------------------------------------------------
# Thresholds. All loose first guesses, all overridable per call.
# --------------------------------------------------------------------------

SEED_CLUSTER_THRESHOLD = 5.0
"""Interface RMSD below which two seeds are called the same answer."""

REDUNDANT_SPEARMAN = 0.8
"""Above this, two score columns are one opinion with two names."""

NOISE_SPEARMAN = 0.1
"""Below this, two columns share nothing, which means at least one of them is
noise - the pair cannot say which."""

FULL_SCAN_ALPHABET = 19
"""Mutations per position in a full scan: the 20 amino acids less the wild
type, which `enumerate_variants` does not emit."""

SECOND_SHELL_CUTOFF = 10.0
"""Within this of the antigen but not in contact with it: the layer that moves
the contact residues without touching anything itself."""

BURIAL_RADIUS = 10.0
BURIED_NEIGHBOURS = 20
"""Burial as a CA contact number rather than relative solvent accessibility.
Cruder, and it is the only one of the two that survives the CA-only files this
package writes for fixtures and that a backbone-only prediction gives you. The
area in `interface_size` is the place that needs real side chains."""

TYPICAL_BSA_PER_SIDE = (700.0, 900.0)
"""Angstroms squared buried per side at an antibody-antigen interface, from the
literature. Reported against, never enforced: a measured interface outside it
is a thing to look at, not an error."""

TYPICAL_INTERFACE_RESIDUES = (15, 25)
"""Residues per side at an antibody-antigen interface, likewise."""

MIN_PLAUSIBLE_INTERFACE = 10
"""Fewer contact residues than this on either side and something is wrong with
the chain assignment or the pose rather than unusual about the epitope. This is
the cheapest place in the whole pipeline to catch that."""

PLAUSIBLE_BENEFICIAL = 0.25
"""Most mutations to a matured antibody are neutral or deleterious. A column
calling appreciably more than this share improving is mis-signed or
miscalibrated, and that is worth knowing before it inverts the shortlist."""

EXPECTED_WT_ARGMAX = 0.6
"""A matrix scorer should mostly prefer the residue that is already there, on a
molecule that affinity maturation already optimised. Materially below this and
the run is invalid rather than merely weak."""

TOP_MUTATIONS = 100
"""The shortlist length `positional_entropy` is asked about."""


def _frame(table: Any) -> "pd.DataFrame":
    """The DataFrame behind either a frame or one of this package's tables."""
    frame = getattr(table, "frame", None)
    return table if frame is None else frame


def _columns(frame: "pd.DataFrame", columns: Sequence[str] | None) -> list[str]:
    """The score columns to work over: the ones named, or all of them.

    Named columns are taken as given - a caller passing `["score_esm_if"]` is
    being deliberate - but checked, because a typo would otherwise read as
    "this scorer had no opinion" rather than as a mistake.
    """
    if columns is None:
        return score_columns(frame)
    missing = [c for c in columns if c not in frame.columns]
    if missing:
        raise KeyError(
            f"no column(s) {', '.join(missing)}; present: "
            f"{', '.join(sorted(frame.columns))}"
        )
    return list(columns)


# --------------------------------------------------------------------------
# Step 1: the structure
# --------------------------------------------------------------------------


def seed_agreement(
    paths: Sequence[Path | str],
    chains: Sequence[str] | None = None,
    partners: Sequence[str] | None = None,
    threshold: float = SEED_CLUSTER_THRESHOLD,
) -> tuple[float, int]:
    """(fraction of seeds in the largest cluster, number of clusters).

    Single-linkage clustering of the ensemble by interface RMSD. Eighteen of
    twenty seeds agreeing on one epitope is stronger evidence than any single
    ipTM, because the documented bottleneck in antibody-antigen prediction is
    identifying the correct model among those generated rather than generating
    one at all - and disagreement is the one signal here that does not come
    from the model's own self-assessment.

    It still says nothing about correctness: twenty seeds can agree on the
    wrong epitope, and published decoy controls show they often do. What it
    rules out is the case where the predictor had no answer and the top pick
    was a coin toss.

    One path is (1.0, 1) by definition and says nothing; ask for seeds.
    """
    paths = [Path(p) for p in paths]
    if not paths:
        raise ValueError("seed_agreement needs at least one prediction")
    count = len(paths)
    if count == 1:
        return 1.0, 1

    parent = list(range(count))

    def find(i: int) -> int:
        while parent[i] != i:
            parent[i] = parent[parent[i]]
            i = parent[i]
        return i

    for i in range(count):
        for j in range(i + 1, count):
            distance = interface_rmsd(paths[j], paths[i], chains, partners)
            if distance == distance and distance <= threshold:
                a, b = find(i), find(j)
                if a != b:
                    parent[a] = b

    sizes: dict[int, int] = {}
    for i in range(count):
        root = find(i)
        sizes[root] = sizes.get(root, 0) + 1
    return max(sizes.values()) / float(count), len(sizes)


def decoy_gap(real_iptm: float | None, decoy_iptm: float | None) -> float:
    """real - decoy. The control that says whether confidence means anything.

    These predictors assign confident scores to deliberately shuffled,
    non-cognate antibody-antigen pairs, so an ipTM of 0.8 is only evidence if a
    mismatched antigen scores materially lower on the same settings. A gap near
    zero means confidence is uninformative *for this system* - not that the
    prediction is wrong, which is a different and less actionable thing.
    """
    if real_iptm is None or decoy_iptm is None:
        return float("nan")
    return float(real_iptm) - float(decoy_iptm)


# --------------------------------------------------------------------------
# Step 2: enumerate
# --------------------------------------------------------------------------


def variant_count(
    variants: Any,
    sequence_length: int | None = None,
    *,
    alphabet: int = FULL_SCAN_ALPHABET,
    strict: bool = False,
) -> dict[str, Scalar]:
    """Candidate and row counts, against what a full scan should have produced.

    A full scan of a chain is `len(chain) * 19` candidates, exactly. Anything
    less means positions were dropped somewhere - an X the reader could not
    resolve, a `positions` argument left in from a smaller test run - and the
    consequence is a shortlist chosen from a set nobody intended, which is
    invisible downstream because every later step is perfectly happy with a
    short table.

    `complete` is that comparison, and is None when no `sequence_length` was
    given to compare against. `strict=True` raises on a mismatch instead, for
    a caller that wants the run to stop rather than the number recorded.
    """
    frame = _frame(variants)
    candidates = int(frame["variant_id"].nunique())
    positions = int(frame[["chain", "position"]].drop_duplicates().shape[0])
    expected = None if sequence_length is None else int(sequence_length) * int(alphabet)
    complete = None if expected is None else candidates == expected
    if strict and complete is False:
        raise ValueError(
            f"a full scan of {sequence_length} residues is {expected} "
            f"candidates and this table has {candidates}, over {positions} "
            f"position(s). Check `positions=`, and check for X residues the "
            f"structure reader could not resolve."
        )
    return {
        "candidates": candidates,
        "rows": int(len(frame)),
        "positions": positions,
        "expected": expected,
        "complete": complete,
    }


def position_buckets(
    structure: Path | str,
    chain: str,
    *,
    antibody_chains: Sequence[str] | None = None,
    partners: Sequence[str] | None = None,
    cdr_positions: Iterable[int] | None = None,
    cutoff: float = INTERFACE_CUTOFF,
    second_shell_cutoff: float = SECOND_SHELL_CUTOFF,
    buried_neighbours: int = BURIED_NEIGHBOURS,
) -> dict[str, Scalar]:
    """Every position of `chain`, counted by what kind of position it is.

    Six classes, as a partition in this precedence order: `interface` (within
    `cutoff` of the antigen), `second_shell` (within `second_shell_cutoff` of
    it but not in contact), `cdr_non_contact`, `vh_vl_interface`,
    `buried_core`, `surface_framework`. The precedence is a choice and the
    classes genuinely overlap - a CDR residue eight angstroms off the antigen
    is both second shell and CDR - so it is written down here rather than left
    to be inferred from the numbers.

    What it is for: the mix tells you whether the scan you are about to pay for
    is aimed at the interface at all. A scan whose interface bucket holds three
    positions is not an affinity-maturation scan, whatever the chain is called,
    and a buried core bucket full of mutations is a stability experiment.

    `cdr_positions` is passed in rather than derived, because deriving it needs
    IMGT labels and therefore ANARCI, and nothing in this module is allowed to
    reach for a model or a tool. `numbering.imgt_labels` plus `denovo.IMGT_CDR`
    is where a caller gets them; without them `cdr_known` is False and that
    bucket is empty rather than quietly folded into the others.
    """
    structure = Path(structure)
    side, antigen = split_chains(structure, antibody_chains, partners)
    if chain not in side:
        raise KeyError(
            f"chain {chain!r} is not on the antibody side ({', '.join(side)}); "
            f"pass antibody_chains= to say which chains are the antibody"
        )
    length = len(chain_sequence(structure, chain))
    partner_side = [c for c in side if c != chain]

    def positions_of(target: Sequence[str], radius: float) -> set[int]:
        if not target:
            return set()
        return {
            position
            for found_chain, position in interface_residues(
                structure, [chain], target, radius
            )
            if found_chain == chain
        }

    contact = positions_of(antigen, cutoff)
    near = positions_of(antigen, second_shell_cutoff) - contact
    vh_vl = positions_of(partner_side, cutoff)
    cdr = {int(p) for p in cdr_positions} if cdr_positions is not None else set()
    buried = {
        position
        for position, neighbours in enumerate(
            _contact_numbers(structure, chain, BURIAL_RADIUS), start=1
        )
        if neighbours >= buried_neighbours
    }

    counts = dict.fromkeys(
        (
            "interface",
            "second_shell",
            "cdr_non_contact",
            "vh_vl_interface",
            "buried_core",
            "surface_framework",
        ),
        0,
    )
    for position in range(1, length + 1):
        if position in contact:
            counts["interface"] += 1
        elif position in near:
            counts["second_shell"] += 1
        elif position in cdr:
            counts["cdr_non_contact"] += 1
        elif position in vh_vl:
            counts["vh_vl_interface"] += 1
        elif position in buried:
            counts["buried_core"] += 1
        else:
            counts["surface_framework"] += 1
    return {
        **counts,
        "positions": length,
        "cdr_known": cdr_positions is not None,
    }


def _contact_numbers(structure: Path, chain: str, radius: float) -> np.ndarray:
    """CA neighbours within `radius` for each position of `chain`, itself excluded.

    A residue whose CA is not modelled gets 0 - "cannot tell" reads as not
    buried, which is the direction that does not invent a buried position out
    of a gap in the file.
    """
    import biotite.structure as struc

    atoms = load_structure(structure)
    atoms = atoms[struc.filter_amino_acids(atoms)]
    everything = atoms[atoms.atom_name == "CA"].coord
    own = ca_coordinates(structure, chain)
    counts = np.zeros(len(own), dtype=int)
    known = np.isfinite(own).all(axis=1)
    if not known.any() or not len(everything):
        return counts
    distances = np.linalg.norm(
        everything[None, :, :] - own[known][:, None, :], axis=-1
    )
    # Minus one for the residue's own CA, which is in `everything`.
    counts[known] = (distances <= radius).sum(axis=1) - 1
    return counts


def interface_size(
    structure: Path | str,
    chains: Sequence[str] | None = None,
    partners: Sequence[str] | None = None,
    cutoff: float = INTERFACE_CUTOFF,
    *,
    min_residues: int = MIN_PLAUSIBLE_INTERFACE,
) -> dict[str, Scalar]:
    """How big the interface is: residues each side, and area buried.

    An antibody-antigen interface buries roughly 700-900 A^2 per side across
    15-25 residues (`TYPICAL_BSA_PER_SIDE`, `TYPICAL_INTERFACE_RESIDUES`). A
    computed interface of four residues does not mean an unusual epitope; it
    means the chain assignment is wrong, or the pose is, and every
    structure-based scorer downstream will run happily on it and return
    numbers. This is the cheapest place in the pipeline to catch that, which is
    the entire reason the function exists.

    `plausible` is the residue-count check only. The buried area is reported
    but not judged: it needs real side chains to mean anything, and half the
    structures this package handles are backbone-only predictions or the
    CA-only files the test suite writes.
    """
    structure = Path(structure)
    side, other = split_chains(structure, chains, partners)
    near = interface_residues(structure, side, other, cutoff)
    far = interface_residues(structure, other, side, cutoff)
    buried = _buried_area(structure, side, other)
    return {
        "n_residues_side": len(near),
        "n_residues_partner": len(far),
        "n_residues": len(near) + len(far),
        "buried_area": buried,
        "buried_area_per_side": buried / 2.0,
        "side": ",".join(side),
        "partner": ",".join(other),
        "plausible": min(len(near), len(far)) >= int(min_residues),
    }


def _buried_area(
    structure: Path, side: Sequence[str], other: Sequence[str]
) -> float:
    """SASA(side) + SASA(partner) - SASA(complex), in angstroms squared.

    One radius for every atom rather than biotite's per-element table: these
    structures are routinely backbone-only, where a per-element radius is
    precision about coordinates that are not there.
    """
    import biotite.structure as struc

    atoms = load_structure(structure)
    atoms = atoms[struc.filter_amino_acids(atoms)]

    def area(chains: Sequence[str]) -> float:
        subset = atoms[np.isin(atoms.chain_id, list(chains))]
        if subset.array_length() == 0:
            return 0.0
        return float(np.nansum(struc.sasa(subset, vdw_radii="Single")))

    return area(side) + area(other) - area([*side, *other])


# --------------------------------------------------------------------------
# Step 3: score. The richest group, and the one that earns its keep.
# --------------------------------------------------------------------------


def inter_model_spearman(
    frame: Any, columns: Sequence[str] | None = None
) -> dict[tuple[str, str], float]:
    """Rank correlation between every pair of score columns.

    The single most informative diagnostic in the pipeline, because it is the
    only one that tests the premise the whole ensemble rests on - that the
    columns are independent opinions. Near `REDUNDANT_SPEARMAN` and two columns
    are one opinion under two names, so agreement between them is not
    corroboration and the ensemble is thinner than its width suggests. Near
    `NOISE_SPEARMAN` and the two share nothing, which means at least one of
    them is noise; it cannot say which, and that is what the harness is for.

    Keyed by the scorer names with the `score_` prefix stripped, each pair
    sorted, so the result is stable across column orderings.
    """
    frame = _frame(frame)
    chosen = _columns(frame, columns)
    pairwise: dict[tuple[str, str], float] = {}
    for index, left in enumerate(chosen):
        for right in chosen[index + 1 :]:
            key = tuple(sorted((scorer_name(left), scorer_name(right))))
            pairwise[key] = spearman(frame[left], frame[right])
    return pairwise


def beneficial_fraction(
    frame: Any,
    columns: Sequence[str] | None = None,
    *,
    threshold: float = 0.0,
) -> dict[str, float]:
    """Share of candidates each column calls improving.

    Score columns in this package are oriented so that higher is better
    binding, which makes "improving" the share above `threshold`. Most
    mutations to an already-matured antibody are neutral or deleterious, so a
    column reporting 40% beneficial is mis-signed or miscalibrated rather than
    optimistic - and a flipped sign convention caught here is a flipped sign
    convention that does not invert the shortlist. `PLAUSIBLE_BENEFICIAL` is
    the loose upper bound to read it against.

    NaN for a column with nothing finite in it, which is a scorer that
    declined, not a scorer that found nothing good.
    """
    frame = _frame(frame)
    shares: dict[str, float] = {}
    for column in _columns(frame, columns):
        values = np.asarray(frame[column], dtype=float)
        usable = values[np.isfinite(values)]
        shares[scorer_name(column)] = (
            float((usable > threshold).mean()) if len(usable) else float("nan")
        )
    return shares


def wt_argmax_rate(frame: Any, column: str) -> dict[str, Scalar]:
    """How often a matrix scorer's best residue at a position is the wild type.

    On a molecule that affinity maturation has already optimised, a model that
    understands it should mostly prefer the residue that is there. A low rate
    does not weaken the run, it invalidates it: the model is not reading this
    antibody, and every number in the column is then about something else.
    Read it against `EXPECTED_WT_ARGMAX`.

    Needs the wild type to be among the scored residues, and the ordinary
    variant table does not have it - `enumerate_variants` skips wt-to-wt
    because it is not a variant. So this takes a *matrix* run: the scorer over
    the full alphabet including the wild type. Positions with no wild-type row
    are counted in `skipped` and left out of the rate rather than guessed at,
    and the rate is NaN when that is all of them.
    """
    frame = _frame(frame)
    if column not in frame.columns:
        raise KeyError(f"no column {column!r} in this table")
    missing = [c for c in ("chain", "position", "wt", "mut") if c not in frame.columns]
    if missing:
        raise KeyError(
            f"wt_argmax_rate reads a variant table and this one is missing "
            f"{', '.join(missing)}"
        )

    wins = 0
    scored = 0
    skipped = 0
    for _, group in frame.groupby(["chain", "position"], sort=True):
        usable = group[np.isfinite(np.asarray(group[column], dtype=float))]
        wild = set(usable["wt"].astype(str))
        if not len(usable) or not wild & set(usable["mut"].astype(str)):
            skipped += 1
            continue
        scored += 1
        best = usable.loc[usable[column].idxmax()]
        if str(best["mut"]) == str(best["wt"]):
            wins += 1
    return {
        "scorer": scorer_name(column),
        "positions": scored + skipped,
        "scored": scored,
        "skipped": skipped,
        "rate": wins / float(scored) if scored else float("nan"),
    }


def positional_entropy(
    frame: Any, column: str, *, top: int = TOP_MUTATIONS
) -> dict[str, Scalar]:
    """How spread the top `top` mutations are across positions.

    A hundred mutations clustered on three positions is three bets, not a
    hundred, and that is the number the submission's diversity constraint
    actually cares about. `effective_positions` is the one to read: it is
    2**entropy, so a perfectly even spread over nine positions reports nine,
    and a shortlist that is 90% one position reports close to one whatever its
    length.

    Entropy is over the distribution of the top rows across (chain, position),
    in bits. `normalized` divides by the most it could have been given how many
    rows were used, so it is comparable between a top-10 and a top-100.
    """
    frame = _frame(frame)
    if column not in frame.columns:
        raise KeyError(f"no column {column!r} in this table")
    usable = frame[np.isfinite(np.asarray(frame[column], dtype=float))]
    chosen = usable.nlargest(int(top), column, keep="all").head(int(top))
    if not len(chosen):
        return {
            "n": 0,
            "positions": 0,
            "entropy": float("nan"),
            "normalized": float("nan"),
            "effective_positions": float("nan"),
            "top_position_share": float("nan"),
        }

    counts = chosen.groupby(["chain", "position"], sort=False).size().to_numpy()
    shares = counts / float(counts.sum())
    entropy = float(-(shares * np.log2(shares)).sum())
    ceiling = math.log2(len(chosen)) if len(chosen) > 1 else 0.0
    return {
        "n": int(len(chosen)),
        "positions": int(len(counts)),
        "entropy": entropy,
        "normalized": entropy / ceiling if ceiling else float("nan"),
        "effective_positions": float(2.0**entropy),
        "top_position_share": float(shares.max()),
    }


def seed_rank_stability(
    per_seed: Mapping[str, Any] | Sequence[Any], column: str
) -> dict[str, Scalar]:
    """Rank correlation of one score column across structures from other seeds.

    The link back to step 1, and the question step 1 cannot answer on its own.
    Score the same variants against each seed's structure: a column that
    reorders when the structure changes is measuring structure noise, not
    mutation effects, and no amount of seed agreement upstream makes it mean
    anything. This is what says whether the structure was good enough to score
    on at all.

    Takes the per-seed score tables, joined here on `variant_id` so that only
    the candidates every seed answered are compared. `mean_spearman` is what to
    report and `min_spearman` is what to worry about - one bad seed is enough
    to tell you the column is reading coordinates.
    """
    frames = list(per_seed.values()) if isinstance(per_seed, Mapping) else list(per_seed)
    labels = (
        [str(k) for k in per_seed]
        if isinstance(per_seed, Mapping)
        else [str(i) for i in range(len(frames))]
    )
    if len(frames) < 2:
        raise ValueError(
            "seed_rank_stability compares seeds and needs at least two score "
            "tables; one seed cannot disagree with itself"
        )

    series = []
    for label, table in zip(labels, frames):
        frame = _frame(table)
        if column not in frame.columns:
            raise KeyError(f"seed {label!r} has no column {column!r}")
        unique = frame.drop_duplicates("variant_id").set_index("variant_id")[column]
        series.append(unique.rename(label))

    import pandas as pd

    joined = pd.concat(series, axis=1, join="inner").dropna()
    correlations = [
        spearman(joined[a], joined[b])
        for index, a in enumerate(joined.columns)
        for b in joined.columns[index + 1 :]
    ]
    return {
        "seeds": len(frames),
        "shared_candidates": int(len(joined)),
        "mean_spearman": nanmean(correlations),
        "min_spearman": min((c for c in correlations if c == c), default=float("nan")),
        "scorer": scorer_name(column),
    }


# --------------------------------------------------------------------------
# Step 4: the gate
# --------------------------------------------------------------------------


def check_names(frame: Any) -> list[str]:
    """The developability checks a gated table carries, from its `pass_` columns."""
    return sorted(
        column[len("pass_") :] for column in _frame(frame).columns
        if column.startswith("pass_")
    )


def pass_rate(frame: Any, checks: Sequence[str] | None = None) -> dict[str, float]:
    """Share of candidates each check passed.

    Per check, never as one survival rate: the point is to find the check that
    is misbehaving, and a pooled number hides it. A liabilities check dropping
    60% of a scan means the regex is wrong or the parent antibody is unusual,
    and those are both worth a look before the survivors go forward - a gate
    this package applies by dropping, so what it removes is gone.
    """
    frame = _frame(frame)
    names = list(checks) if checks is not None else check_names(frame)
    rates: dict[str, float] = {}
    for check in names:
        column = f"pass_{check}"
        if column not in frame.columns:
            raise KeyError(
                f"no column {column!r}; this table carries "
                f"{', '.join(check_names(frame)) or 'no checks'}"
            )
        verdicts = frame.drop_duplicates("variant_id")[column].astype(bool)
        rates[check] = float(verdicts.mean()) if len(verdicts) else float("nan")
    return rates


def gate_score_correlation(
    frame: Any,
    columns: Sequence[str] | None = None,
    checks: Sequence[str] | None = None,
) -> dict[tuple[str, str], float]:
    """Does each check preferentially kill the high-scoring variants?

    Rank correlation between a check's verdict (1 for pass) and a score column,
    keyed by (check, scorer). Negative means the gate is removing the
    candidates the scorer liked most - a real tension between binding and
    developability, and one worth surfacing rather than silently absorbing,
    because the shortlist that comes out the other side is then the
    second-best-binding set and nobody said so.

    NaN where a check passed or failed everything: there is no correlation to
    compute against a constant, and zero would read as "no tension".
    """
    frame = _frame(frame)
    chosen = _columns(frame, columns)
    names = list(checks) if checks is not None else check_names(frame)
    candidates = frame.drop_duplicates("variant_id")
    out: dict[tuple[str, str], float] = {}
    for check in names:
        column = f"pass_{check}"
        if column not in frame.columns:
            raise KeyError(f"no column {column!r} in this table")
        verdict = candidates[column].astype(bool).astype(float)
        for score in chosen:
            out[(check, scorer_name(score))] = spearman(candidates[score], verdict)
    return out


# --------------------------------------------------------------------------
# Step 5: pairs
# --------------------------------------------------------------------------


def epistasis(pairs: Any, singles: Any, column: str) -> "pd.DataFrame":
    """One row per pair: its score, its singles' sum, and the difference.

    score(AB) - [score(A) + score(B)], which is directly measurable here
    because this package scores a pair as its own candidate rather than adding
    its singles together. Both summaries below are built on this, so the join -
    a pair's `variant_id` is its singles' ids joined with '+' - lives in one
    place.

    Pairs whose singles are not both in the table are dropped: the alternative
    is a sum over one mutation, which would read as an epistasis equal to the
    missing single's score.
    """
    import pandas as pd

    pair_frame = _frame(pairs).drop_duplicates("variant_id")
    single_frame = _frame(singles).drop_duplicates("variant_id")
    for name, table in (("pairs", pair_frame), ("singles", single_frame)):
        if column not in table.columns:
            raise KeyError(f"the {name} table has no column {column!r}")
    lookup = single_frame.set_index("variant_id")[column].to_dict()

    rows: list[dict[str, Any]] = []
    for row in pair_frame.itertuples():
        parts = str(row.variant_id).split("+")
        if len(parts) < 2:
            continue
        values = [lookup.get(part) for part in parts]
        if any(v is None or v != v for v in values):
            continue
        observed = getattr(row, column)
        if observed != observed:
            continue
        additive = float(sum(values))
        rows.append(
            {
                "variant_id": str(row.variant_id),
                "observed": float(observed),
                "additive": additive,
                "epistasis": float(observed) - additive,
                "distance": float(getattr(row, "distance", float("nan"))),
            }
        )
    return pd.DataFrame(
        rows, columns=["variant_id", "observed", "additive", "epistasis", "distance"]
    )


def epistasis_magnitude(
    pairs: Any, singles: Any, column: str, *, tolerance: float | None = None
) -> dict[str, Scalar]:
    """The distribution of score(AB) - [score(A) + score(B)] over the pairs.

    What it is really testing is whether re-scoring the pairs bought anything.
    If the epistasis is near zero across every pair, the scorer is additive by
    construction - it is reading each mutation independently and summing - and
    the second, expensive `score_variants` call could have been a sum. If it is
    not near zero, the pairs genuinely needed scoring, which is the assumption
    this pipeline is built on and is worth confirming on our own data rather
    than on the literature's.

    `tolerance` is what counts as near zero. Score units differ per column, so
    it defaults to a tenth of the spread of the pair scores themselves rather
    than to an absolute number that would mean something different for every
    scorer.
    """
    table = epistasis(pairs, singles, column)
    values = np.asarray(table["epistasis"], dtype=float)
    values = values[np.isfinite(values)]
    if not len(values):
        return {
            "pairs": 0,
            "mean": float("nan"),
            "median": float("nan"),
            "std": float("nan"),
            "max_abs": float("nan"),
            "tolerance": float("nan"),
            "near_additive_fraction": float("nan"),
        }
    observed = np.asarray(table["observed"], dtype=float)
    spread = float(np.nanstd(observed[np.isfinite(observed)]))
    cut = float(tolerance) if tolerance is not None else 0.1 * spread
    return {
        "pairs": int(len(values)),
        "mean": float(values.mean()),
        "median": float(np.median(values)),
        "std": float(values.std()),
        "max_abs": float(np.abs(values).max()),
        "tolerance": cut,
        "near_additive_fraction": (
            float((np.abs(values) <= cut).mean()) if cut > 0 else float("nan")
        ),
    }


def distance_vs_epistasis(pairs: Any, singles: Any, column: str) -> dict[str, Scalar]:
    """Does distance between the two mutations predict how additive they are?

    `select_pairs` keeps pairs at least 10 A apart on the premise that distant
    mutations do not interact, so their effects add and a pair is a safe
    combination of two singles. That is a reasonable premise and it is also
    testable on our own output, which is what this does: a negative
    `spearman_abs` - bigger distance, smaller departure from additivity - is
    the rule working. Near zero means the 10 A cutoff is not buying what it is
    there to buy, and the pair-selection rule is then a convention rather than
    a constraint.

    NaN when the pairs table carries no distances, which is what a table built
    without a structure looks like.
    """
    table = epistasis(pairs, singles, column)
    distance = np.asarray(table["distance"], dtype=float)
    values = np.asarray(table["epistasis"], dtype=float)
    usable = np.isfinite(distance) & np.isfinite(values)
    return {
        "pairs": int(usable.sum()),
        "spearman_abs": spearman(distance[usable], np.abs(values[usable])),
        "spearman_signed": spearman(distance[usable], values[usable]),
    }


# --------------------------------------------------------------------------
# Step 6: the shortlist. Composition, not quality.
# --------------------------------------------------------------------------


def position_diversity(frame: Any) -> dict[str, Scalar]:
    """How many distinct positions the shortlist touches, and how evenly.

    A submission check rather than a quality one: nothing here says the
    positions chosen are the right positions. What it says is whether a
    shortlist of fifty candidates is fifty bets or five, which is a property of
    the list that no scorer reports and that the diversity constraint is about.
    """
    frame = _frame(frame)
    candidates = int(frame["variant_id"].nunique())
    per_position = frame.drop_duplicates(["variant_id", "chain", "position"])
    counts = per_position.groupby(["chain", "position"], sort=False).size()
    return {
        "candidates": candidates,
        "positions": int(len(counts)),
        "candidates_per_position": (
            float(counts.mean()) if len(counts) else float("nan")
        ),
        "max_per_position": int(counts.max()) if len(counts) else 0,
        "largest_position_share": (
            float(counts.max() / candidates) if candidates else float("nan")
        ),
    }


def chain_balance(frame: Any) -> dict[str, float]:
    """Share of the shortlist's mutations on each chain.

    Heavy-chain CDRs carry most of the affinity in most antibodies, so a
    shortlist skewed to H is expected rather than wrong. It is worth seeing all
    the same: a shortlist that is 100% one chain has left half the paratope
    unexplored, and a shortlist that is 100% framework has probably mistaken
    which chain is which.
    """
    frame = _frame(frame)
    rows = frame.drop_duplicates(["variant_id", "chain", "position"])
    counts = rows.groupby("chain", sort=True).size()
    total = float(counts.sum())
    return {str(chain): float(n / total) for chain, n in counts.items()} if total else {}


def edit_distance_histogram(frame: Any) -> dict[int, int]:
    """How many candidates carry one mutation, two, three ...

    Mutations per candidate, which for this package's tables is rows sharing a
    `variant_id`. Flat and small, and it catches the thing a glance at a
    shortlist does not: a table that was meant to be the paired set and is
    still all singles, or a `select_pairs` run whose output kept its singles.
    """
    frame = _frame(frame)
    sizes = frame.groupby("variant_id", sort=False).size()
    histogram = sizes.value_counts()
    return {int(k): int(v) for k, v in sorted(histogram.items())}


def final_rank_agreement(
    frame: Any, columns: Sequence[str] | None = None
) -> dict[str, Scalar]:
    """How much the scorers agree about the order *of the shortlist*.

    `inter_model_spearman` over the survivors rather than over the whole scan,
    which is a different question and the one that matters at the end: the
    columns can disagree wildly across 4,598 candidates and still agree about
    the forty that got through, or the reverse. Low agreement here means the
    composition of the submission is an artifact of which column it was sorted
    on, and that is worth knowing before it is sorted.
    """
    pairwise = inter_model_spearman(frame, columns)
    values = list(pairwise.values())
    return {
        "scorers": len(_columns(_frame(frame), columns)),
        "pairs": len(values),
        "mean_spearman": nanmean(values),
        "min_spearman": min((v for v in values if v == v), default=float("nan")),
        "max_spearman": max((v for v in values if v == v), default=float("nan")),
    }


# --------------------------------------------------------------------------
# Into the schema
# --------------------------------------------------------------------------


def flatten(values: Mapping[Any, Any], prefix: str = "") -> dict[str, Scalar]:
    """Any of the above, as the flat scalar map `metrics.save` takes.

    Several functions here return a nested dict or one keyed by a tuple of
    scorer names, because that is the shape the analysis wants. `metrics.json`
    is deliberately flat, so the translation happens once, here, rather than in
    every caller that wants to record a diagnostic.

        metrics.save(run_dir, diagnostics.flatten(rates, "pass_rate"))
    """
    flat: dict[str, Scalar] = {}
    for key, value in values.items():
        name = "_".join(str(part) for part in key) if isinstance(key, tuple) else str(key)
        full = f"{prefix}_{name}" if prefix else name
        if isinstance(value, Mapping):
            flat.update(flatten(value, full))
        elif isinstance(value, (float, int, str, bool, type(None))):
            flat[full] = value
        else:
            flat[full] = str(value)
    return flat
