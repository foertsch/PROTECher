"""Everything that computes a number. Pure functions, one vocabulary.

No adapters, no subprocesses, no registry, no cache. The notebook and the
validation harness both call in here rather than each growing their own
statistics, which is the only way the two stay comparable: a Spearman in the
harness report and a Spearman in a notebook plot have to be the same Spearman.

The functions fall into two categories and the split is the point.

**Category 1 - ground truth.** Everything in this module. Each one needs a
known answer - a deposited structure, a measured affinity - so they are
available in the validation harness and the validation notebook and nowhere
else. They are the only things in this package that speak to whether a
prediction is *correct*.

* Structure comparison: `dockq`, `interface_rmsd`, `epitope_recall`,
  `interface_residues`. How close a predicted complex is to a known one.
* Ranking: `spearman`, `precision_at_k`, `top_k_recall`,
  `magnitude_agreement`. How well a column orders candidates against measured
  affinity - and whether its magnitudes carry information or only its order
  does, which decides whether candidates can be picked by threshold or only by
  rank.

**Category 2 - diagnostics**, in `design_pipeline.metrics.diagnostics`. They
need no answer key, so they run on the competition target. They catch bugs,
redundancy and known-bad regimes; they cannot say whether a prediction is
right, and that module's docstring says so at length because it is the thing
most likely to be forgotten.

The two are never blended. There is no summary score here that mixes them,
because a number combining "the seeds agreed with each other" and "the
predictions were correct" would be read as the second.

Plus one schema. Every step that produces numbers writes `metrics.json` into
its cache directory as a flat map, and `load` reads it back; nothing parses
another step's output by hand.

Positions are 1-based into the chain as given in the file, the same convention
the variant table uses, so an interface residue and a variant can be compared
without a second numbering scheme. Comparing two *different* files - a
prediction against a deposited structure - goes through `align_positions`,
because a predictor numbers its output 1..N while the PDB uses author numbering
with gaps and insertion codes.

`numpy` and `pandas` are core. DockQ - and biopython, which it brings - is an
extra (`design-pipeline[metrics]`), imported inside the one function that needs
it so that `import design_pipeline` stays light. `diagnostics` is reached
through this module's `__getattr__` for the same reason: it speaks pandas, and
nothing here that reads a structure needs pandas loaded to do it.
"""

from __future__ import annotations

import json
from pathlib import Path
from typing import Any, Iterable, Mapping, Sequence

import numpy as np

from .errors import MissingDependencyError
# `load` is imported under another name: this module defines its own load()
# for metrics.json, and the structure reader would be shadowed by it.
from .structures import (
    ca_coordinates,
    chain_ids,
    chain_sequence,
    residue_starts,
)
from .structures import load as load_structure

FILENAME = "metrics.json"

# DockQ's published quality bands. Named rather than inlined because a table in
# a notebook that says "0.49" means nothing to a reader and "medium" does.
ACCEPTABLE = 0.23
MEDIUM = 0.49
HIGH = 0.80

INTERFACE_CUTOFF = 5.0
"""Heavy-atom distance defining a contact, in angstroms. 5 A is the usual
choice for an epitope definition."""

IRMSD_CUTOFF = 10.0
"""Wider cutoff for interface RMSD, matching DockQ's own definition: the
superposition needs the interface's whole neighbourhood, not only its contacts."""

MIN_PAIRS = 3
"""Below this, a rank correlation is noise with a number attached, and a
superposition is underdetermined."""

MAGNITUDE_GAP = 0.2
"""How far Spearman has to exceed Pearson before a column is called
rank-only. A loose first guess, like every threshold in this package: it is
here so a report can say "rank-only" rather than print two correlations and
leave the reader to compare them."""

SCORE_PREFIX = "score_"
RAW_SUFFIX = "_raw"
"""A scorer's column is `score_<model>`; `score_<model>_raw` is the same
opinion before orientation, kept for the record and never correlated against
anything, since its sign convention is the model's rather than ours."""


def score_columns(frame: Any) -> list[str]:
    """The score columns of a frame, in order, raw ones excluded.

    One definition, used by the harness report and by every diagnostic that
    sweeps the columns. It was written out at three call sites before this
    existed, and the third one disagreed with the other two about `_raw`.
    """
    return [
        column
        for column in frame.columns
        if column.startswith(SCORE_PREFIX) and not column.endswith(RAW_SUFFIX)
    ]


def scorer_name(column: str) -> str:
    """`score_esm_if` -> `esm_if`. What a report column is headed with."""
    return column[len(SCORE_PREFIX) :] if column.startswith(SCORE_PREFIX) else column


# --------------------------------------------------------------------------
# Category 1: ground truth - residue correspondence
# --------------------------------------------------------------------------


def align_positions(a: str, b: str) -> dict[int, int]:
    """{position in `a`: position in `b`} for the aligned pairs, 1-based.

    Needleman-Wunsch with unit scores, written out rather than pulled from
    biopython: the job is mapping a deposited chain's author numbering onto a
    prediction's 1..N, the sequences are near-identical, and thirty lines that
    the test suite can check beat a dependency that it cannot.

    Gapped positions are simply absent from the map, so a caller that looks up
    a residue the other structure does not model gets a miss rather than a
    silently wrong neighbour.
    """
    if a == b:
        return {i: i for i in range(1, len(a) + 1)}
    n, m = len(a), len(b)
    if not n or not m:
        return {}

    gap = -1.0
    score = np.zeros((n + 1, m + 1))
    score[:, 0] = np.arange(n + 1) * gap
    score[0, :] = np.arange(m + 1) * gap
    # 0 = diagonal, 1 = up (gap in b), 2 = left (gap in a)
    trace = np.zeros((n + 1, m + 1), dtype=np.int8)
    trace[1:, 0] = 1
    trace[0, 1:] = 2
    for i in range(1, n + 1):
        row = score[i - 1]
        for j in range(1, m + 1):
            diagonal = row[j - 1] + (1.0 if a[i - 1] == b[j - 1] else -1.0)
            up = row[j] + gap
            left = score[i][j - 1] + gap
            best = max(diagonal, up, left)
            score[i][j] = best
            trace[i][j] = 0 if best == diagonal else (1 if best == up else 2)

    mapping: dict[int, int] = {}
    i, j = n, m
    while i > 0 and j > 0:
        move = trace[i][j]
        if move == 0:
            mapping[i] = j
            i -= 1
            j -= 1
        elif move == 1:
            i -= 1
        else:
            j -= 1
    return dict(sorted(mapping.items()))


def sequence_identity(a: str, b: str) -> float:
    """Fraction of aligned positions that match, over the shorter sequence.

    Used to deduplicate a benchmark set by antigen: without it a SAbDab
    selection fills up with near-identical SARS-CoV-2 RBD complexes and the
    sample size is a fiction.
    """
    if not a or not b:
        return 0.0
    mapping = align_positions(a, b)
    matches = sum(1 for i, j in mapping.items() if a[i - 1] == b[j - 1])
    return matches / float(min(len(a), len(b)))


# --------------------------------------------------------------------------
# Category 1: ground truth - structure comparison
# --------------------------------------------------------------------------


def split_chains(
    structure: Path, chains: Sequence[str] | None, partners: Sequence[str] | None
) -> tuple[list[str], list[str]]:
    """(side, other side). Defaults to the antibody against everything else.

    The default reads H and L as the antibody when they are present, because
    that is what this package's own predictions are labelled with, and falls
    back to the first chain against the rest. Both guesses are worth naming
    explicitly in anything that matters, which is why both are arguments.
    """
    present = list(chain_ids(structure))
    if chains:
        side = [c for c in chains if c in present]
        if not side:
            raise KeyError(
                f"{structure.name} has none of chain(s) {', '.join(chains)}; "
                f"it has {', '.join(present)}"
            )
    else:
        paired = [c for c in ("H", "L") if c in present]
        side = paired or present[:1]
    other = (
        [c for c in partners if c in present]
        if partners
        else [c for c in present if c not in side]
    )
    if not other:
        raise ValueError(
            f"{structure.name} has no chain outside {', '.join(side)}, so it "
            f"has no interface to measure"
        )
    return side, other


def interface_residues(
    structure: Path | str,
    chains: Sequence[str] | None = None,
    partners: Sequence[str] | None = None,
    cutoff: float = INTERFACE_CUTOFF,
) -> set[tuple[str, int]]:
    """{(chain, position)} of residues in `chains` within `cutoff` of `partners`.

    Positions are 1-based into the chain as given. Call it twice to get both
    sides: `interface_residues(p, ["H", "L"], ["A"])` is the paratope and
    `interface_residues(p, ["A"], ["H", "L"])` the epitope.
    """
    import biotite.structure as struc

    structure = Path(structure)
    side, other = split_chains(structure, chains, partners)
    atoms = load_structure(structure)
    atoms = atoms[struc.filter_amino_acids(atoms)]
    partner_atoms = atoms[np.isin(atoms.chain_id, other)]
    if partner_atoms.array_length() == 0:
        return set()
    # Neighbour search rather than a full distance matrix: a Fab against an
    # antigen is a few thousand atoms each way, and this stays linear.
    cells = struc.CellList(partner_atoms, cell_size=max(cutoff, 1.0))

    found: set[tuple[str, int]] = set()
    for chain in side:
        sub, starts = residue_starts(atoms, chain)
        contacts = cells.get_atoms(sub.coord, radius=cutoff)
        touching = (
            (contacts != -1).any(axis=1)
            if contacts.ndim == 2 and contacts.shape[1]
            else np.zeros(sub.array_length(), dtype=bool)
        )
        bounds = list(starts) + [sub.array_length()]
        for index, (begin, end) in enumerate(zip(bounds, bounds[1:]), start=1):
            if touching[begin:end].any():
                found.add((chain, index))
    return found


def _paired_coordinates(
    predicted: Path,
    reference: Path,
    residues: Iterable[tuple[str, int]],
) -> tuple[np.ndarray, np.ndarray]:
    """CA coordinates of `residues`, in both files, for the ones both model.

    `residues` are reference positions. Each chain's sequences are aligned, so
    a prediction numbered 1..N lines up with a deposited chain that starts at
    residue 3 and skips a disordered loop.
    """
    by_chain: dict[str, list[int]] = {}
    for chain, position in residues:
        by_chain.setdefault(chain, []).append(position)

    reference_points: list[np.ndarray] = []
    predicted_points: list[np.ndarray] = []
    for chain, positions in sorted(by_chain.items()):
        mapping = align_positions(
            chain_sequence(reference, chain), chain_sequence(predicted, chain)
        )
        reference_coords = ca_coordinates(reference, chain)
        predicted_coords = ca_coordinates(predicted, chain)
        for position in sorted(positions):
            target = mapping.get(position)
            if target is None:
                continue
            a = reference_coords[position - 1]
            b = predicted_coords[target - 1]
            if np.isfinite(a).all() and np.isfinite(b).all():
                reference_points.append(a)
                predicted_points.append(b)
    if not reference_points:
        return np.empty((0, 3)), np.empty((0, 3))
    return np.asarray(reference_points), np.asarray(predicted_points)


def interface_rmsd(
    predicted: Path | str,
    reference: Path | str,
    chains: Sequence[str] | None = None,
    partners: Sequence[str] | None = None,
    cutoff: float = IRMSD_CUTOFF,
) -> float:
    """CA RMSD over the reference interface, after optimal superposition on it.

    The interface is defined on the *reference*, so a prediction that puts the
    paratope somewhere else is measured against where it should have been
    rather than against wherever it happens to touch.

    CA atoms rather than the full backbone DockQ uses: it works on the CA-only
    files this package writes for fixtures, and `dockq` reports the
    authoritative iRMSD when that package is installed. NaN when fewer than
    three residues correspond.
    """
    import biotite.structure as struc

    predicted, reference = Path(predicted), Path(reference)
    side, other = split_chains(reference, chains, partners)
    residues = interface_residues(
        reference, side, other, cutoff
    ) | interface_residues(reference, other, side, cutoff)
    fixed, mobile = _paired_coordinates(predicted, reference, residues)
    if len(fixed) < MIN_PAIRS:
        return float("nan")
    fitted = struc.superimpose(fixed, mobile).apply(mobile)
    return float(struc.rmsd(fixed, fitted))


def epitope_recall(
    predicted: Path | str,
    reference: Path | str,
    antigen_chains: Sequence[str] | None = None,
    antibody_chains: Sequence[str] | None = None,
    cutoff: float = INTERFACE_CUTOFF,
) -> float:
    """Fraction of the true epitope that the prediction also puts at the interface.

    The one that matters downstream. A prediction can score passably on DockQ
    while sitting on the wrong patch of the antigen, and every structure-based
    scorer downstream would then be reading out a fiction with no sign that
    anything is wrong. NaN when the reference has no epitope to recall.
    """
    predicted, reference = Path(predicted), Path(reference)
    if antigen_chains is None:
        side, other = split_chains(reference, antibody_chains, None)
        antibody_chains, antigen_chains = side, other
    true = interface_residues(reference, antigen_chains, antibody_chains, cutoff)
    if not true:
        return float("nan")
    found = interface_residues(predicted, antigen_chains, antibody_chains, cutoff)

    # Predicted positions translated into reference positions, so the two sets
    # are comparable even though the files are numbered differently.
    translated: set[tuple[str, int]] = set()
    for chain in {c for c, _ in found}:
        reverse = {
            target: source
            for source, target in align_positions(
                chain_sequence(reference, chain), chain_sequence(predicted, chain)
            ).items()
        }
        for other_chain, position in found:
            if other_chain != chain:
                continue
            mapped = reverse.get(position)
            if mapped is not None:
                translated.add((chain, mapped))
    return len(true & translated) / float(len(true))


def quality(score: float) -> str:
    """DockQ's band for a score: incorrect, acceptable, medium or high."""
    if score != score:
        return "unknown"
    if score >= HIGH:
        return "high"
    if score >= MEDIUM:
        return "medium"
    if score >= ACCEPTABLE:
        return "acceptable"
    return "incorrect"


def summarize_dockq(interfaces: Mapping[str, Mapping[str, Any]]) -> dict[str, Any]:
    """DockQ's per-interface output, reduced to one flat row.

    Separate from `dockq` so the parsing is testable without the package
    installed - the component names are DockQ's and they have moved between
    versions, which is exactly the kind of thing worth a test.

    The best-scoring interface is the one reported: a Fab against an antigen
    has an H-antigen and an L-antigen interface, and the question is whether
    the complex is right, not which chain carries more of the contact.
    """
    if not interfaces:
        return {"dockq": float("nan"), "quality": "unknown", "interfaces": 0}

    def component(row: Mapping[str, Any], *names: str) -> float:
        for name in names:
            if name in row and row[name] is not None:
                return float(row[name])
        return float("nan")

    best_name, best = max(
        interfaces.items(), key=lambda item: component(item[1], "DockQ", "dockq")
    )
    score = component(best, "DockQ", "dockq")
    return {
        "dockq": score,
        "quality": quality(score),
        "fnat": component(best, "fnat", "F1"),
        "irmsd": component(best, "iRMSD", "irmsd", "iRMS"),
        "lrmsd": component(best, "LRMSD", "lrmsd", "LRMS"),
        "interfaces": len(interfaces),
        "interface": str(best_name),
    }


def dockq(
    predicted: Path | str,
    reference: Path | str,
    chain_map: Mapping[str, str] | None = None,
) -> dict[str, Any]:
    """DockQ for a predicted complex against a deposited one, with components.

    Wraps the reference implementation rather than reimplementing it: DockQ is
    a composite of fnat, iRMSD and LRMSD with published weights and a decade of
    convention behind the exact definitions, and a second implementation of it
    would only be a second set of bugs. The components come back alongside the
    score because the score alone cannot tell a near-miss from a wrong epitope.

    `chain_map` maps predicted chain ids to reference ones, which is usually
    needed: the deposited structure calls the heavy chain B and a prediction
    from this package calls it H.
    """
    predicted, reference = Path(predicted), Path(reference)
    try:
        from DockQ.DockQ import load_PDB, run_on_all_native_interfaces
    except ModuleNotFoundError as exc:
        raise MissingDependencyError(
            "DockQ is not installed. `uv pip install "
            "'design-pipeline[metrics]'`. Without it, `interface_rmsd` and "
            "`epitope_recall` still work and need nothing extra."
        ) from exc

    model = load_PDB(str(predicted))
    native = load_PDB(str(reference))
    interfaces, total = run_on_all_native_interfaces(
        model, native, chain_map=dict(chain_map) if chain_map else None
    )
    summary = summarize_dockq(interfaces)
    summary["dockq_total"] = float(total) if total is not None else float("nan")
    return summary


# --------------------------------------------------------------------------
# Category 1: ground truth - ranking against what was measured
# --------------------------------------------------------------------------


def rankdata(values: np.ndarray) -> np.ndarray:
    """1-based ranks, ties averaged - the convention Spearman assumes.

    Ties are the normal case here: a matrix scorer gives the same number to
    every variant at a position it has no opinion about, and a measured set
    often reports a detection limit as one repeated value.
    """
    values = np.asarray(values, dtype=float)
    order = np.argsort(values, kind="stable")
    ranks = np.empty(len(values), dtype=float)
    sorted_values = values[order]
    start = 0
    while start < len(values):
        stop = start
        while stop + 1 < len(values) and sorted_values[stop + 1] == sorted_values[start]:
            stop += 1
        ranks[order[start : stop + 1]] = 0.5 * (start + stop) + 1.0
        start = stop + 1
    return ranks


def spearman(predicted: np.ndarray, measured: np.ndarray) -> float:
    """Rank correlation, NaN when it cannot be computed.

    NaN rather than 0.0 for too few points or a constant column: zero would
    read as "this scorer has no skill", which is a claim, where NaN is the
    absence of one.
    """
    predicted = np.asarray(predicted, dtype=float)
    measured = np.asarray(measured, dtype=float)
    usable = np.isfinite(predicted) & np.isfinite(measured)
    if usable.sum() < MIN_PAIRS:
        return float("nan")
    a = rankdata(predicted[usable])
    b = rankdata(measured[usable])
    if a.std() == 0 or b.std() == 0:
        return float("nan")
    return float(np.corrcoef(a, b)[0, 1])


def precision_at_k(predicted: np.ndarray, measured: np.ndarray, k: int = 10) -> float:
    """Overlap between the scorer's top k and the measured top k, over k.

    The question a campaign actually asks: of the handful we would order, how
    many are among the ones that turned out best. Defined as set overlap rather
    than "how many beat the wild type", because that needs a wild-type
    measurement the published sets do not always include - a definition this
    module skips rather than guesses at.
    """
    predicted, measured = _usable(predicted, measured)
    if len(predicted) < k:
        return float("nan")
    top_predicted = set(np.argsort(-predicted, kind="stable")[:k].tolist())
    top_measured = set(np.argsort(-measured, kind="stable")[:k].tolist())
    return len(top_predicted & top_measured) / float(k)


def top_k_recall(
    predicted: np.ndarray, measured: np.ndarray, k: int = 20, n_true: int = 10
) -> float:
    """How much of the measured top `n_true` a shortlist of `k` would capture.

    `precision_at_k` with k == n_true is the same number; the two differ when
    the shortlist is larger than the set of things that actually worked, which
    is the realistic case - you order 20 variants hoping to catch the 5 good
    ones. NaN when there are fewer than `n_true` measurements.
    """
    predicted, measured = _usable(predicted, measured)
    if len(predicted) < max(k, n_true) or n_true <= 0:
        return float("nan")
    shortlist = set(np.argsort(-predicted, kind="stable")[:k].tolist())
    best = set(np.argsort(-measured, kind="stable")[:n_true].tolist())
    return len(shortlist & best) / float(n_true)


def magnitude_agreement(
    predicted: np.ndarray,
    measured: np.ndarray,
    *,
    gap: float = MAGNITUDE_GAP,
) -> dict[str, float]:
    """Does score *magnitude* predict effect size, or only score order?

    Spearman answers whether a scorer ranks candidates correctly. It cannot
    answer the question that decides how the shortlist is cut: if a mutation
    scoring twice as high does not bind twice as tightly, then only the
    ordering carries information and candidates have to be picked by rank
    rather than by a threshold on the score. A threshold on a column whose
    magnitudes are meaningless is an arbitrary cut dressed up as a criterion.

    So this returns both correlations plus the linear fit of measured on
    predicted. `rank_only` is True when the rank correlation is respectable and
    the linear one is materially worse by `gap` - the signature of a column
    that orders well and scales badly. It is a flag for a human to look at, not
    a verdict: `gap` is a first guess and the fit is a straight line through
    data that has no reason to be straight.

    Only answerable against measured affinities, so it is Category 1 and lives
    with the rest of the ranking metrics: the harness is the only place in this
    package that can call it.
    """
    predicted, measured = _usable(predicted, measured)
    if len(predicted) < MIN_PAIRS or predicted.std() == 0 or measured.std() == 0:
        return {
            "n": float(len(predicted)),
            "pearson": float("nan"),
            "spearman": float("nan"),
            "slope": float("nan"),
            "intercept": float("nan"),
            "r2": float("nan"),
            "rank_only": False,
        }
    pearson = float(np.corrcoef(predicted, measured)[0, 1])
    rho = spearman(predicted, measured)
    slope, intercept = np.polyfit(predicted, measured, 1)
    return {
        "n": float(len(predicted)),
        "pearson": pearson,
        "spearman": rho,
        "slope": float(slope),
        "intercept": float(intercept),
        "r2": pearson**2,
        "rank_only": bool(
            rho == rho and abs(rho) >= 0.3 and abs(rho) - abs(pearson) > gap
        ),
    }


def nanmean(values: Sequence[float]) -> float:
    """Mean of the values that are not NaN, NaN when none are.

    Here rather than in the harness because the harness is not allowed its own
    statistics, and because the choice it encodes - average each scorer's
    per-system correlations rather than pool the rows - is a decision about
    what the number means, not a formatting detail.
    """
    usable = [float(v) for v in values if v == v]
    return sum(usable) / len(usable) if usable else float("nan")


def _usable(predicted: np.ndarray, measured: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    predicted = np.asarray(predicted, dtype=float)
    measured = np.asarray(measured, dtype=float)
    keep = np.isfinite(predicted) & np.isfinite(measured)
    return predicted[keep], measured[keep]


# --------------------------------------------------------------------------
# The schema
# --------------------------------------------------------------------------

Scalar = float | int | str | bool | None


def save(run_dir: Path | str, metrics: Mapping[str, Scalar]) -> Path:
    """Write `metrics.json` into a run directory. Flat map, scalars only.

    Flatness is the point and is enforced: the whole reason this file exists is
    that a notebook and the harness can both read any step's numbers without
    knowing which step wrote them, and a nested payload would put that
    knowledge back into every reader.
    """
    run_dir = Path(run_dir)
    bad = {
        key: type(value).__name__
        for key, value in metrics.items()
        if not isinstance(value, (float, int, str, bool, type(None)))
    }
    if bad:
        raise TypeError(
            f"metrics.json holds a flat map of scalars; these are not: {bad}. "
            f"Flatten the names instead, e.g. 'confidence_iptm'."
        )
    run_dir.mkdir(parents=True, exist_ok=True)
    path = run_dir / FILENAME
    path.write_text(json.dumps(dict(metrics), indent=1, sort_keys=True))
    return path


def load(run_dir: Path | str) -> dict[str, Scalar]:
    """Read `metrics.json` back. Raises if the step never wrote one."""
    path = Path(run_dir) / FILENAME
    if not path.is_file():
        raise FileNotFoundError(
            f"no {FILENAME} in {run_dir}. Steps that produce numbers write one; "
            f"use load_all() to sweep a tree and skip the ones that do not."
        )
    return json.loads(path.read_text())


def load_all(root: Path | str) -> dict[Path, dict[str, Scalar]]:
    """Every `metrics.json` under `root`, keyed by its run directory.

    How the notebook collects a campaign without knowing its layout.
    """
    found: dict[Path, dict[str, Scalar]] = {}
    for path in sorted(Path(root).rglob(FILENAME)):
        try:
            found[path.parent] = json.loads(path.read_text())
        except json.JSONDecodeError:
            continue  # a half-written file is not a result
    return found


# --------------------------------------------------------------------------
# Category 2 lives next door
# --------------------------------------------------------------------------


def __getattr__(name: str):
    """`metrics.diagnostics` - Category 2, loaded on first use.

    A deliberate seam rather than a plain import at the top of the file. Two
    reasons, in order: the categories stay separate, so nothing in Category 1
    can quietly come to depend on a diagnostic and no caller can reach one
    without naming `diagnostics` out loud; and the diagnostics speak pandas,
    which the structure comparisons above do not need and which `import
    design_pipeline` should not pay for.

    `from design_pipeline.metrics import diagnostics` works, because that is
    what PEP 562 is for.
    """
    if name == "diagnostics":
        from . import _diagnostics

        return _diagnostics
    raise AttributeError(f"module {__name__!r} has no attribute {name!r}")
