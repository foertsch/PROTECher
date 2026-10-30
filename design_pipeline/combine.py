"""Combine single substitutions into variants of 5 to 10, and choose designs among the scored ones.

`select_pairs` stops at two. Track A needs at least five substitutions per design, so this module takes the gated
singles to sets of five to ten, and then picks the designs that go forward. It is three plain functions and a
loop, like the rest of `variants`: nothing here runs a model, so nothing here is an adapter.

    sets = propose_sets(survivors, by="score_antifold", regions=regions)     # singles -> candidate sets
    scored = score_variants(complex_, sets, models=["antifold_seq"], ...)    # the caller scores them whole
    designs = pick_designs(scored, by="score_antifold_seq", regions=regions, quota={5: 6, 6: 4, 7: 3, 9: 2})
    table = frontier(scored, by="score_antifold_seq")                        # best score at each mutation count

Scoring is not done here, on purpose. A set's score is not a function of its singles' scores, so nothing in this
module carries a single's score forward or adds two together; `propose_sets` drops every `score_` column, as
`select_pairs` does, and the sets are scored as whole candidates by `score_variants`. Of the scorers that can do
that, `antifold_seq` is the one that scores a candidate with several mutations at all. Six sizes at 300 sets each
come to 1,800 candidates, under the 2,000 that `score_variants` accepts for a per-variant scorer without a cap.

Two kinds of rule, in two places. A rule about one set (positions distinct, mutations not on top of each other,
no more than so many in one region) is enforced while the set is drawn. A rule about the designs taken together
(each mutation in so many designs, designs that differ enough, a quota per size) can only be enforced while they are
chosen, so it lives in `pick_designs`. Region labels come in as a plain mapping from (chain, position) to a label,
because the classes the plan uses (paratope, CDR-H3, second shell, framework) come from the position map and not
from anything in this package.
"""

from __future__ import annotations

import warnings
from collections import Counter
from pathlib import Path
from typing import TYPE_CHECKING, Callable, Iterable, Mapping, Sequence

import numpy as np

from .errors import MissingDependencyError
from .structures import chain_sequence
from .variants import (
    Mutation,
    VariantTable,
    _coordinates,
    _distance,
    _imgt_lookup,
    _rows_for,
    make_id,
    parse_id,
    write_table,
)

if TYPE_CHECKING:
    import pandas as pd

Regions = Mapping[tuple[str, int], str]
"""(chain, position) -> a label such as "H3" or "framework". A position with no entry has no label and is never limited."""

ATTEMPTS_PER_SET = 20
"""Random draws tried for each set wanted before a size is given up on. A draw that cannot reach the size, because
the pool is too crowded or a limit is hit, counts as an attempt."""


def propose_sets(
    variants: VariantTable,
    sizes: Iterable[int] = range(5, 11),
    per_size: int = 300,
    output_dir: Path | str | None = None,
    *,
    by: str,
    pool: int = 30,
    min_distance: float = 4.0,
    regions: Regions | None = None,
    limits: Mapping[str, int] | None = None,
    seed: int = 0,
    structure: Path | str | None = None,
) -> VariantTable:
    """Draw candidate sets of several substitutions from the best singles. Returns them as a variant table.

    `by` names the score column that ranks the singles, so the ranking is the caller's and this function never
    combines scores; the best `pool` of them are the shortlist. The pool sets how often a substitution can recur
    across the designs picked later. Fifteen designs of about six substitutions are some 90 slots, and the plan wants
    each substitution in two to four designs, so a pool of about 30 suits it better than the shortlist's full 60. In
    a dry run on 1N8Z with the offline scorer (mechanics only, not science), 15 designs came from a pool of 30 with
    one substitution of 29 in fewer than two designs, against 22 of 51 from a pool of 60. Too small a pool, 25
    there, leaves too few designs that differ enough. Every set is a random greedy draw from the shortlist: the
    singles are shuffled and taken in turn when
    they fit, until the set has its size. A draw that cannot reach the size is dropped, and a set seen twice is
    kept once, so the number returned per size can fall short of `per_size`; `config["found"]` says by how much.

    A set obeys these, checked as it grows:

    * one substitution per position;
    * no two substitutions closer than `min_distance` angstroms, CA to CA, read off the structure the table
      carries (residues along the chain without one, as in `select_pairs`). The plan's "within 4 Å" rule is
      written for singles scored apart; the sets here are scored whole, so what is left is the physical one, two
      neighbours changed at once, which 4 Å CA to CA excludes only for adjacent residues. Raise it to be stricter;
    * at most `limits[label]` substitutions whose position has that label in `regions`, for example
      `limits={"H3": 1}`. Positions with no label are not limited.

    The same seed gives the same sets. Score columns are dropped, because a set's score is not a function of its
    singles'; score the sets with `score_variants`.
    """
    import pandas as pd

    frame = variants.frame
    sizes = sorted({int(k) for k in sizes})
    if not sizes or min(sizes) < 2:
        raise ValueError(f"sizes must be whole numbers of at least 2, got {sizes}")
    per_candidate = frame.groupby("variant_id")["position"].size()
    multi = per_candidate[per_candidate > 1]
    if len(multi):
        raise ValueError(
            f"propose_sets takes single mutations and {len(multi)} of these candidates already have more than "
            f"one (e.g. {multi.index[0]!r}). Pass the singles table, or the gated singles."
        )
    if by not in frame.columns:
        raise KeyError(f"no column {by!r} to rank on; present: {', '.join(sorted(frame.columns))}")
    ranked = frame[frame[by].notna()].sort_values(by, ascending=False, kind="stable").head(pool)
    if len(ranked) < max(sizes):
        raise ValueError(
            f"only {len(ranked)} singles have a {by!r} score, fewer than the {max(sizes)} a set of that size needs. "
            f"Widen the singles, or lower the largest size."
        )
    limits = dict(limits or {})
    if limits and regions is None:
        raise ValueError("limits needs regions: the labels the limits are counted over")
    regions = dict(regions or {})

    target = Path(structure) if structure is not None else variants.structure
    singles = [
        (str(r.variant_id), str(r.chain), int(r.position), str(r.wt), str(r.mut)) for r in ranked.itertuples()
    ]
    labels = [regions.get((s[1], s[2])) for s in singles]
    coords = _coordinates(target, {s[1] for s in singles})
    n = len(singles)
    compatible = np.zeros((n, n), dtype=bool)
    for i in range(n):
        for j in range(i):
            if (singles[i][1], singles[i][2]) == (singles[j][1], singles[j][2]):
                continue  # two substitutions of one position are not a set
            distance = _distance(coords, singles[i], singles[j])
            ok = distance is not None and distance >= min_distance
            compatible[i, j] = compatible[j, i] = ok

    rng = np.random.default_rng(seed)
    rows: list[dict] = []
    found: dict[int, int] = {}
    label_lookup = _imgt_lookup(ranked)
    for size in sizes:
        drawn: dict[str, list[int]] = {}
        attempts = 0
        while len(drawn) < per_size and attempts < per_size * ATTEMPTS_PER_SET:
            attempts += 1
            chosen: list[int] = []
            in_label: Counter[str] = Counter()
            for i in rng.permutation(n):
                if not all(compatible[i, j] for j in chosen):
                    continue
                label = labels[i]
                if label in limits and in_label[label] >= limits[label]:
                    continue
                chosen.append(int(i))
                if label is not None:
                    in_label[label] += 1
                if len(chosen) == size:
                    break
            if len(chosen) == size:
                mutations = [(singles[i][1], singles[i][2], singles[i][3], singles[i][4]) for i in chosen]
                drawn.setdefault(make_id(mutations), chosen)
        found[size] = len(drawn)
        for chosen in drawn.values():
            mutations = [(singles[i][1], singles[i][2], singles[i][3], singles[i][4]) for i in chosen]
            for row in _rows_for(mutations, label_lookup):
                rows.append({**row, "n_mutations": size})

    if not rows:
        raise ValueError(
            f"no set of sizes {sizes} could be drawn from these {n} singles under min_distance={min_distance} "
            f"and limits={limits}. Relax a limit, lower min_distance, or widen the singles."
        )
    out = Path(output_dir) if output_dir is not None else variants.output_dir / "sets"
    return write_table(
        pd.DataFrame(rows),
        out / "sets.parquet",
        output_dir=variants.output_dir,
        structure=target,
        config={
            "sizes": sizes,
            "per_size": per_size,
            "found": found,
            "pool": n,
            "by": by,
            "min_distance": min_distance,
            "limits": limits,
            "seed": seed,
            "from": str(variants.path),
            "distance_units": "angstrom" if target is not None else "residues",
        },
    )


def pick_designs(
    scored: VariantTable,
    *,
    by: str,
    n: int = 15,
    min_differences: int = 3,
    max_uses: int = 4,
    min_uses: int = 2,
    regions: Regions | None = None,
    use_limits: Mapping[str, int] | None = None,
    cdr_labels: Iterable[str] = (),
    quota: Mapping[int, int] | None = None,
    max_pi_shift: float | None = 1.0,
    sequences: Mapping[str, str] | None = None,
    pi_function: Callable[[str], float] | None = None,
    output_dir: Path | str | None = None,
) -> VariantTable:
    """Choose `n` designs from scored sets, best first, under the rules that hold across designs.

    Candidates are taken in order of `by`, highest first. Ties go to fewer substitutions in `cdr_labels`, then to
    fewer substitutions overall, then to the id, so the same table always gives the same designs. A candidate is
    kept when all of these hold:

    * it differs from every design already kept at `min_differences` positions or more, counting a position as
      different when the two designs put different residues there, the parent's included;
    * none of its substitutions is already in `max_uses` designs, and one whose label is in `use_limits` is
      in no more designs than that, for example `use_limits={"H3": 2}` for the plan's CDR-H3 rule;
    * its size still has room in `quota`, a mapping from number of substitutions to designs. Without a quota the
      best scores win outright, which favours the smallest sets, because every extra substitution costs score.
      With one, a size not named gets no designs;
    * its isoelectric point is within `max_pi_shift` of the parent's. The point is computed on the chains joined
      in name order, from `sequences` (chain to sequence) or else from the structure the table carries. For a Fab
      structure that includes the constant domains and so understates the shift of an scFv; pass the scFv's own
      chains in `sequences`. `pi_function` replaces the calculation; the default is Biopython's, imported when
      first needed. Pass `max_pi_shift=None` to skip the rule.

    The plan also wants each substitution in at least `min_uses` designs. A greedy pass cannot guarantee a lower
    bound, so it is reported rather than enforced: `config["underused"]` lists the substitutions in fewer designs,
    and a warning says so, as it does when fewer than `n` designs could be chosen. Scores of the chosen sets are
    kept as they are, because each is the set's own, scored whole.
    """
    import pandas as pd

    frame = scored.frame
    if by not in frame.columns:
        raise KeyError(f"no column {by!r} to rank on; present: {', '.join(sorted(frame.columns))}")
    regions = dict(regions or {})
    use_limits = dict(use_limits or {})
    cdr_labels = set(cdr_labels)
    if use_limits and not regions:
        raise ValueError("use_limits needs regions: the labels the limits are counted over")

    values = frame.groupby("variant_id", sort=False)[by].first()
    candidates = []
    for variant_id, score in values.items():
        if pd.isna(score):
            continue  # a candidate the scorer could not score cannot be ranked
        mutations = parse_id(str(variant_id))
        cdr_changes = sum(regions.get((m[0], m[1])) in cdr_labels for m in mutations)
        candidates.append((-float(score), cdr_changes, len(mutations), str(variant_id), mutations))
    candidates.sort(key=lambda c: c[:4])
    unscored = len(values) - len(candidates)

    shift_of: Callable[[Sequence[Mutation]], float] | None = None
    if max_pi_shift is not None:
        shift_of = _pi_shift_function(scored, sequences, pi_function)

    chosen: list[tuple[str, Sequence[Mutation]]] = []
    uses: Counter[Mutation] = Counter()
    by_size: Counter[int] = Counter()
    for _negative, _cdr, size, variant_id, mutations in candidates:
        if len(chosen) >= n:
            break
        if quota is not None and by_size[size] >= quota.get(size, 0):
            continue
        if any(_hamming(mutations, other) < min_differences for _, other in chosen):
            continue
        if any(uses[m] >= min(max_uses, use_limits.get(regions.get((m[0], m[1])), max_uses)) for m in mutations):
            continue
        if shift_of is not None and shift_of(mutations) >= max_pi_shift:
            continue
        chosen.append((variant_id, mutations))
        by_size[size] += 1
        uses.update(mutations)

    underused = sorted(m for m, count in uses.items() if count < min_uses)
    if len(chosen) < n:
        warnings.warn(
            f"chose {len(chosen)} of {n} designs: the rules left no more candidates among {len(candidates)} scored. "
            f"Propose more sets, relax min_differences or max_uses, or widen the quota.",
            RuntimeWarning,
            stacklevel=2,
        )
    if underused:
        warnings.warn(
            f"{len(underused)} substitutions are in fewer than {min_uses} designs (e.g. {make_id([underused[0]])}); "
            f"see config['underused'].",
            RuntimeWarning,
            stacklevel=2,
        )
    if not chosen:
        raise ValueError(f"no design could be chosen from {len(candidates)} scored candidates ({unscored} unscored)")

    rank = {variant_id: i + 1 for i, (variant_id, _) in enumerate(chosen)}
    picked = frame[frame["variant_id"].isin(rank)].copy()
    picked["design"] = picked["variant_id"].map(rank)
    picked["n_mutations"] = picked["variant_id"].map(picked.groupby("variant_id")["position"].size())
    picked = picked.sort_values(["design"], kind="stable")
    out = Path(output_dir) if output_dir is not None else scored.output_dir / "designs"
    return write_table(
        picked,
        out / "designs.parquet",
        output_dir=scored.output_dir,
        structure=scored.structure,
        config={
            "n": n,
            "chosen": len(chosen),
            "by": by,
            "min_differences": min_differences,
            "max_uses": max_uses,
            "min_uses": min_uses,
            "use_limits": use_limits,
            "quota": dict(quota) if quota is not None else None,
            "max_pi_shift": max_pi_shift,
            "unscored": unscored,
            "underused": [make_id([m]) for m in underused],
            "uses": {make_id([m]): c for m, c in sorted(uses.items())},
            "from": str(scored.path),
        },
    )


def frontier(scored: VariantTable, by: str, keep: int = 1) -> "pd.DataFrame":
    """The best `keep` candidates at each number of substitutions, for the jury's score-against-count plot.

    One row per kept candidate, with `n_mutations`, `variant_id` and the `by` column, ordered by size and then by
    score. This is the plan's "best-scoring sets for each mutation count"; reading where the curve flattens is how
    the mutation-count quota is set. A candidate with no score is left out.
    """
    frame = scored.frame
    if by not in frame.columns:
        raise KeyError(f"no column {by!r} to rank on; present: {', '.join(sorted(frame.columns))}")
    per_candidate = frame.groupby("variant_id", sort=False).agg(n_mutations=("position", "size"), **{by: (by, "first")})
    per_candidate = per_candidate.dropna(subset=[by]).reset_index()
    ordered = per_candidate.sort_values(["n_mutations", by], ascending=[True, False], kind="stable")
    return ordered.groupby("n_mutations", sort=True).head(keep).reset_index(drop=True)[["n_mutations", "variant_id", by]]


def _hamming(left: Sequence[Mutation], right: Sequence[Mutation]) -> int:
    """Positions at which two candidates put different residues, with the parent's residue where one has no change."""
    a = {(c, p): m for c, p, _wt, m in left}
    b = {(c, p): m for c, p, _wt, m in right}
    return sum(a.get(position) != b.get(position) for position in a.keys() | b.keys())


def _pi_shift_function(
    scored: VariantTable, sequences: Mapping[str, str] | None, pi_function: Callable[[str], float] | None
) -> Callable[[Sequence[Mutation]], float]:
    """A function from a candidate's mutations to the absolute shift of the isoelectric point from the parent's."""
    calculate = pi_function or _biopython_pi
    chains = dict(sequences) if sequences is not None else {}
    if not chains:
        if scored.structure is None:
            raise ValueError(
                "the isoelectric point rule needs the chain sequences: pass sequences={chain: sequence}, "
                "or a table that carries its structure. max_pi_shift=None skips the rule."
            )
        wanted = sorted(set(scored.frame["chain"].astype(str)))
        chains = {chain: chain_sequence(scored.structure, chain) for chain in wanted}
    order = sorted(chains)
    parent = calculate("".join(chains[c] for c in order))

    def shift(mutations: Sequence[Mutation]) -> float:
        residues = {c: list(chains[c]) for c in order}
        for chain, position, wt, mut in mutations:
            if chain not in residues or not 0 < position <= len(residues[chain]) or residues[chain][position - 1] != wt:
                raise ValueError(f"{make_id([(chain, position, wt, mut)])} does not match the sequence of chain {chain}")
            residues[chain][position - 1] = mut
        return abs(calculate("".join("".join(residues[c]) for c in order)) - parent)

    return shift


def _biopython_pi(sequence: str) -> float:
    try:
        from Bio.SeqUtils.ProtParam import ProteinAnalysis
    except ModuleNotFoundError as exc:
        raise MissingDependencyError(
            "the isoelectric point rule uses Biopython, which is not installed here: `pip install biopython`, "
            "or pass max_pi_shift=None to skip the rule, or pi_function= to bring your own."
        ) from exc
    return float(ProteinAnalysis(sequence).isoelectric_point())
