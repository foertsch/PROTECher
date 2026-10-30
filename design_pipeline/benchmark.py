"""Choosing targets for the structure-prediction validation, and scoring them.

The notebook that uses this asks one question: does `get_structure` recover a
known antibody-antigen complex from sequence alone, and do the confidence
metrics predict when it does not. Everything here is the part of that worth
testing - target selection, validation and the per-target protocol - so the
notebook stays a thin driver over functions with unit tests.

Two tiers, kept apart on purpose.

**Smoke tier.** Three complexes everybody uses. They are in every predictor's
training set, so a good result on them says the plumbing works and nothing
else. Reported in their own table, labelled as memorized, and never averaged
with the honest tier.

**Honest tier.** Selected from SAbDab, filtered to complexes deposited after
the predictor's training cutoff. Raw RCSB attribute queries detect antibodies
badly - they find lysozyme complexes and miss nanobodies - so the
purpose-built resource is used instead. Deduplicated by antigen, because
without that the set fills with near-identical SARS-CoV-2 RBD complexes and
the sample size is a fiction.

Nothing here falls back to a hardcoded honest tier. A silently stale benchmark
is worse than no benchmark, so an unreachable SAbDab is an error that names the
URL.
"""

from __future__ import annotations

import io
import json
import re
from dataclasses import dataclass, replace
from datetime import date, datetime
from pathlib import Path
from typing import TYPE_CHECKING, Any, Callable, Mapping, Sequence
from urllib.error import HTTPError, URLError
from urllib.request import urlopen

from . import metrics
from .errors import MissingDependencyError, PipelineError

if TYPE_CHECKING:
    import pandas as pd

RCSB_ENTRY = "https://data.rcsb.org/rest/v1/core/entry/{pdb_id}"
RCSB_FASTA = "https://www.rcsb.org/fasta/entry/{pdb_id}"
SABDAB_SUMMARY = (
    "https://opig.stats.ox.ac.uk/webapps/sabdab-sabpred/sabdab/summary/all/"
)
"""The summary TSV. Reachable from the Downloads section of
https://opig.stats.ox.ac.uk/webapps/sabdab-sabpred/sabdab/ ; this is the direct
link that page offers."""

TIMEOUT = 60

MAX_RESOLUTION = 3.0
MIN_PROTEIN_ENTITIES = 3
"""Two variable domains and a protein antigen. Fewer means it is a nanobody,
an unbound Fab, or a complex with something that is not a protein."""

IDENTITY_THRESHOLD = 0.70
XRAY = "X-RAY DIFFRACTION"

# AlphaFold3's cutoff is published in the paper's methods. The others are not
# recorded here on purpose: a wrong cutoff silently puts training data into the
# honest tier, and a date invented by this package would be indistinguishable
# from one that was checked. `training_cutoff` returns verified=False for them
# so the notebook has to print where its number came from.
TRAINING_CUTOFFS: Mapping[str, tuple[str, str]] = {
    "alphafold3": ("2021-09-30", "AlphaFold3 paper, Methods: PDB training cutoff"),
}

RELEASE_MARGIN_DAYS = 183
"""Six months. The fallback when a predictor documents no cutoff: its release
date minus this. Crude, and deliberately conservative - a stricter cutoff only
shrinks the pool, while a loose one invalidates the whole exercise."""


class BenchmarkError(PipelineError):
    """A benchmark set could not be built, with the reason and the URL."""


@dataclass(frozen=True)
class Target:
    """One complex to predict and score against."""

    pdb_id: str
    heavy: str
    light: str
    antigen: str
    tier: str = "honest"
    date: str = ""
    resolution: float | None = None
    note: str = ""

    @property
    def chains(self) -> tuple[str, str, str]:
        return (self.heavy, self.light, self.antigen)

    @property
    def mapping(self) -> dict[str, str]:
        """Deposited chain id -> the id this package predicts under.

        Relabelling the reference rather than the prediction: every comparison
        in `metrics` matches chains by id across the two files, and the
        prediction is the thing being scored.
        """
        return {self.heavy: "H", self.light: "L", self.antigen: "A"}


SMOKE_TIER: tuple[Target, ...] = (
    Target(
        pdb_id="1VFB",
        heavy="A",
        light="B",
        antigen="C",
        tier="smoke",
        note="D1.3 Fv against hen egg lysozyme",
    ),
    Target(
        pdb_id="1N8Z",
        heavy="B",
        light="A",
        antigen="C",
        tier="smoke",
        note="trastuzumab Fab against HER2 domain IV; the heavy chain is B",
    ),
    Target(
        pdb_id="3HFM",
        heavy="H",
        light="L",
        antigen="Y",
        tier="smoke",
        note="HyHEL-10 Fab against hen egg lysozyme",
    ),
)
"""The three the spec names. The chain assignments here are what the
literature reports and are **not** to be trusted blindly - `check_target`
verifies every one of them against RCSB before use, and a target that fails is
printed and dropped rather than predicted."""


# --------------------------------------------------------------------------
# Fetching
# --------------------------------------------------------------------------


def _fetch(url: str, timeout: int = TIMEOUT) -> str:
    try:
        with urlopen(url, timeout=timeout) as response:  # noqa: S310 - fixed scheme
            return response.read().decode("utf-8", errors="replace")
    except HTTPError as exc:
        raise BenchmarkError(f"{url} returned HTTP {exc.code}") from exc
    except URLError as exc:
        raise BenchmarkError(f"could not reach {url}: {exc.reason}") from exc


def fetch_entry(pdb_id: str, timeout: int = TIMEOUT) -> dict[str, Any]:
    """RCSB's entry record: resolution, method, entity counts, dates."""
    return json.loads(_fetch(RCSB_ENTRY.format(pdb_id=pdb_id.upper()), timeout))


def parse_rcsb_fasta(text: str) -> dict[str, str]:
    """{auth chain id: sequence} from RCSB's entry FASTA.

    Headers look like `>1N8Z_1|Chains A, B|Immunoglobulin|Mus musculus`, and
    one entity can cover several chains, so each id gets its own entry. Pure,
    because this parse is the part that breaks when RCSB changes its header.
    """
    sequences: dict[str, str] = {}
    for block in text.split(">"):
        if not block.strip():
            continue
        header, _, body = block.partition("\n")
        fields = header.split("|")
        if len(fields) < 2:
            continue
        labels = re.sub(r"^Chains?\s+", "", fields[1].strip())
        sequence = "".join(body.split())
        for label in labels.split(","):
            # "A[auth H]" is how RCSB writes a renamed asym id; the author id
            # in brackets is the one a PDB file uses and the one we need.
            label = label.strip()
            authored = re.search(r"\[auth\s+([^\]]+)\]", label)
            chain = (authored.group(1) if authored else label).strip()
            if chain and sequence:
                sequences[chain] = sequence
    return sequences


def fetch_sequences(pdb_id: str, timeout: int = TIMEOUT) -> dict[str, str]:
    """Every chain's sequence, for chain verification and antigen clustering."""
    return parse_rcsb_fasta(_fetch(RCSB_FASTA.format(pdb_id=pdb_id.upper()), timeout))


# --------------------------------------------------------------------------
# Validation
# --------------------------------------------------------------------------


def check_target(
    target: Target,
    entry: Mapping[str, Any],
    sequences: Mapping[str, str] | None = None,
    *,
    max_resolution: float = MAX_RESOLUTION,
    min_protein_entities: int = MIN_PROTEIN_ENTITIES,
) -> list[str]:
    """Every reason this target is unusable. Empty means it passed.

    Returns reasons rather than raising, so a caller can print them and drop
    the entry - which beats proceeding and discovering at the DockQ step that
    the "antigen" chain was a sulfate ion.
    """
    failures: list[str] = []
    info = entry.get("rcsb_entry_info", {}) or {}

    resolutions = info.get("resolution_combined") or []
    if not resolutions:
        failures.append("no resolution reported (not a crystal structure?)")
    elif min(float(r) for r in resolutions) > max_resolution:
        failures.append(
            f"resolution {min(float(r) for r in resolutions):.2f} A is worse "
            f"than {max_resolution} A"
        )

    entities = info.get("polymer_entity_count_protein")
    if entities is None:
        failures.append("no protein entity count reported")
    elif int(entities) < min_protein_entities:
        failures.append(
            f"{entities} protein polymer entities, need {min_protein_entities}: "
            f"both variable domains and a protein antigen have to be resolved"
        )

    if sequences is not None:
        missing = [c for c in target.chains if c not in sequences]
        if missing:
            failures.append(
                f"chain(s) {', '.join(missing)} are not in the entry; it has "
                f"{', '.join(sorted(sequences))}"
            )
        else:
            antigen = sequences[target.antigen]
            if antigen in (sequences[target.heavy], sequences[target.light]):
                failures.append(
                    f"chain {target.antigen} has the same sequence as an "
                    f"antibody chain, so it is not the antigen"
                )
            if len(antigen) < 25:
                failures.append(
                    f"chain {target.antigen} is {len(antigen)} residues, which "
                    f"is a peptide rather than a protein antigen"
                )
    return failures


def validate_smoke_tier(
    targets: Sequence[Target] = SMOKE_TIER,
    *,
    max_resolution: float = MAX_RESOLUTION,
    fetch: Callable[[str], dict[str, Any]] = fetch_entry,
    sequences: Callable[[str], dict[str, str]] | None = fetch_sequences,
) -> tuple[list[Target], list[tuple[Target, list[str]]]]:
    """(kept, dropped with reasons). Checks the hardcoded list before using it."""
    kept: list[Target] = []
    dropped: list[tuple[Target, list[str]]] = []
    for target in targets:
        try:
            entry = fetch(target.pdb_id)
            seqs = sequences(target.pdb_id) if sequences is not None else None
        except BenchmarkError as exc:
            dropped.append((target, [str(exc)]))
            continue
        failures = check_target(
            target, entry, seqs, max_resolution=max_resolution
        )
        if failures:
            dropped.append((target, failures))
            continue
        info = entry.get("rcsb_entry_info", {}) or {}
        released = (entry.get("rcsb_accession_info", {}) or {}).get(
            "initial_release_date", ""
        )
        kept.append(
            replace(
                target,
                resolution=min(
                    float(r) for r in info.get("resolution_combined") or [0.0]
                )
                or None,
                date=str(released)[:10],
            )
        )
    return kept, dropped


# --------------------------------------------------------------------------
# Training cutoffs
# --------------------------------------------------------------------------


def training_cutoff(model_name: str) -> tuple[str | None, str, bool]:
    """(cutoff, where it came from, verified).

    Only AlphaFold3's is recorded, from its paper. For Boltz-2, Chai-1 and
    anything else this returns None and says to read the model card: a date
    guessed here would be indistinguishable from one that was checked, and
    getting it wrong puts training data into the honest tier, which quietly
    turns the whole benchmark into a memorization test.
    """
    known = TRAINING_CUTOFFS.get(model_name)
    if known is not None:
        return known[0], known[1], True
    return (
        None,
        f"no cutoff recorded for {model_name!r}: read it from the model card "
        f"or repository, or use cutoff_from_release(<release date>)",
        False,
    )


def cutoff_from_release(release_date: str, margin_days: int = RELEASE_MARGIN_DAYS) -> str:
    """The documented fallback: release date minus six months.

    For a predictor whose training cutoff is not published. Conservative by
    construction - too strict only shrinks the pool, too loose invalidates the
    result - and the notebook prints which rule produced its cutoff.
    """
    from datetime import timedelta

    parsed = parse_date(release_date)
    if parsed is None:
        raise ValueError(f"{release_date!r} is not a date this can read")
    return (parsed - timedelta(days=margin_days)).isoformat()


def parse_date(value: Any) -> date | None:
    """ISO, or SAbDab's MM/DD/YY. None when it is neither.

    SAbDab has written its `date` column both ways over the years, so this
    tries both rather than assuming the one that happens to be current.
    """
    text = str(value or "").strip()
    if not text:
        return None
    for pattern in ("%Y-%m-%d", "%m/%d/%y", "%m/%d/%Y", "%d/%m/%Y"):
        try:
            return datetime.strptime(text[:10], pattern).date()
        except ValueError:
            continue
    return None


# --------------------------------------------------------------------------
# SAbDab selection
# --------------------------------------------------------------------------

SABDAB_COLUMNS = (
    "pdb",
    "Hchain",
    "Lchain",
    "antigen_chain",
    "antigen_type",
    "date",
    "resolution",
    "method",
)


def read_sabdab(source: str | Path = SABDAB_SUMMARY) -> "pd.DataFrame":
    """The SAbDab summary TSV, from a URL or a local copy.

    A local path is the reproducible option: the file grows, so a notebook rerun
    six months later selects a different set unless it reads the copy it kept.
    """
    import pandas as pd

    text = (
        Path(source).read_text(encoding="utf-8", errors="replace")
        if not str(source).startswith("http")
        else _fetch(str(source))
    )
    frame = pd.read_csv(io.StringIO(text), sep="\t", dtype=str)
    missing = [c for c in SABDAB_COLUMNS if c not in frame.columns]
    if missing:
        raise BenchmarkError(
            f"{source} is missing column(s) {missing}; expected the SAbDab "
            f"summary file from {SABDAB_SUMMARY}"
        )
    return frame


def cluster_by_identity(
    sequences: Mapping[str, str], threshold: float = IDENTITY_THRESHOLD
) -> list[list[str]]:
    """Single-linkage clusters of keys whose sequences are at least `threshold` identical.

    Quadratic and fine: it runs over the few hundred rows that survive the
    other filters, not over SAbDab.
    """
    keys = list(sequences)
    parent = {key: key for key in keys}

    def find(key: str) -> str:
        while parent[key] != key:
            parent[key] = parent[parent[key]]
            key = parent[key]
        return key

    for i, a in enumerate(keys):
        for b in keys[i + 1 :]:
            if metrics.sequence_identity(sequences[a], sequences[b]) >= threshold:
                ra, rb = find(a), find(b)
                if ra != rb:
                    parent[ra] = rb

    grouped: dict[str, list[str]] = {}
    for key in keys:
        grouped.setdefault(find(key), []).append(key)
    return list(grouped.values())


def select_sabdab(
    frame: "pd.DataFrame",
    cutoff: str,
    *,
    max_resolution: float = MAX_RESOLUTION,
    identity: float = IDENTITY_THRESHOLD,
    limit: int = 6,
    sequences: Callable[[str], dict[str, str]] | None = None,
) -> tuple[list[Target], list[str]]:
    """(targets, a log line per filter step).

    The filters, in the order the surviving count is reported:

    1. a protein antigen - drops peptides, haptens, nucleic acids and unbound
    2. both a heavy and a light chain - drops nanobodies and single domains
    3. an antigen chain distinct from both of them
    4. X-ray, at or better than `max_resolution`
    5. deposited after `cutoff`
    6. one entry per antigen cluster at `identity`, best resolution kept

    `sequences` fetches a candidate's chain sequences for step 6. Without it
    the deduplication falls back to SAbDab's `antigen_name`, which is coarser:
    two names for one protein survive as two entries.
    """
    log: list[str] = [f"{len(frame)} rows in the summary file"]
    working = frame.copy()

    working = working[
        working["antigen_type"].fillna("").str.strip().str.lower() == "protein"
    ]
    log.append(f"{len(working)} with a protein antigen")

    for column in ("Hchain", "Lchain", "antigen_chain"):
        working = working[_present(working[column])]
    log.append(f"{len(working)} with a heavy chain, a light chain and an antigen chain")

    antigen_first = working["antigen_chain"].map(_first_chain)
    working = working[
        (antigen_first != working["Hchain"]) & (antigen_first != working["Lchain"])
    ]
    working = working.assign(antigen_first=antigen_first)
    log.append(f"{len(working)} whose antigen chain is not an antibody chain")

    resolution = working["resolution"].map(_number)
    working = working[
        (working["method"].fillna("").str.upper().str.strip() == XRAY)
        & resolution.notna()
        & (resolution <= max_resolution)
    ]
    working = working.assign(resolution_value=working["resolution"].map(_number))
    log.append(f"{len(working)} solved by X-ray at {max_resolution} A or better")

    limit_date = parse_date(cutoff)
    if limit_date is None:
        raise ValueError(f"{cutoff!r} is not a date")
    parsed = working["date"].map(parse_date)
    working = working[
        parsed.notna() & parsed.map(lambda d: d is not None and d > limit_date)
    ]
    working = working.assign(date_value=working["date"].map(parse_date))
    log.append(f"{len(working)} deposited after {limit_date.isoformat()}")

    # One per antigen. Without this the set is mostly one antigen: SAbDab holds
    # hundreds of near-identical SARS-CoV-2 RBD complexes, and six of those are
    # one data point wearing six hats.
    working = working.sort_values(["date_value", "resolution_value"], ascending=[False, True])
    kept, dedupe_log = _deduplicate(working, identity, sequences)
    log.append(dedupe_log)

    targets = [
        Target(
            pdb_id=str(row.pdb).upper(),
            heavy=str(row.Hchain),
            light=str(row.Lchain),
            antigen=str(row.antigen_first),
            tier="honest",
            date=row.date_value.isoformat() if row.date_value else "",
            resolution=row.resolution_value,
        )
        for row in kept.head(limit).itertuples()
    ]
    log.append(f"{len(targets)} selected, most recent first")
    return targets, log


def _deduplicate(
    working: "pd.DataFrame",
    identity: float,
    sequences: Callable[[str], dict[str, str]] | None,
) -> tuple["pd.DataFrame", str]:
    if working.empty:
        return working, "0 after antigen deduplication"
    if sequences is None:
        if "antigen_name" in working.columns:
            kept = working.drop_duplicates("antigen_name")
            return kept, (
                f"{len(kept)} after deduplicating by antigen_name (no sequence "
                f"fetcher given, so this is coarser than {identity:.0%} identity)"
            )
        return working, "no antigen deduplication possible: no sequences, no antigen_name"

    antigens: dict[str, str] = {}
    for row in working.itertuples():
        try:
            chains = sequences(str(row.pdb))
        except BenchmarkError:
            continue
        sequence = chains.get(str(row.antigen_first))
        if sequence:
            antigens[str(row.pdb)] = sequence
    clusters = cluster_by_identity(antigens, identity)
    # The frame is already sorted by date then resolution, so the first member
    # of each cluster is the one to keep.
    order = {str(row.pdb): i for i, row in enumerate(working.itertuples())}
    representatives = {min(cluster, key=lambda p: order.get(p, 1 << 30)) for cluster in clusters}
    kept = working[working["pdb"].astype(str).isin(representatives)]
    return kept, (
        f"{len(kept)} after clustering antigens at {identity:.0%} identity "
        f"({len(clusters)} cluster(s) from {len(antigens)} fetched)"
    )


def _present(column: "pd.Series") -> "pd.Series":
    text = column.fillna("").astype(str).str.strip()
    return (text != "") & (text.str.upper() != "NA") & (text.str.upper() != "NAN")


def _first_chain(value: Any) -> str:
    """SAbDab writes a multi-chain antigen as 'C | D'. The first is enough:
    predicting one antigen chain is the protocol here."""
    return str(value or "").split("|")[0].strip()


def _number(value: Any) -> float | None:
    try:
        return float(str(value).strip())
    except (TypeError, ValueError):
        return None


# --------------------------------------------------------------------------
# The per-target protocol
# --------------------------------------------------------------------------



def reference_structure(
    target: Target,
    *,
    output_dir: Path | str,
    reference_source: str | Path | None = None,
    cache: bool = True,
) -> tuple[Path, tuple[str, str, str]]:
    """(normalized reference, (heavy, light, antigen) sequences).

    The ground truth, relabelled to H/L/A so a prediction can be compared to
    it chain by chain, and the three sequences read back off it - which are
    exactly what the predictor is given. Reading them from the structure rather
    than from the FASTA matters: the prediction has to be made from the
    residues that were actually *resolved*, or a disordered loop in the
    reference turns into a length mismatch at the comparison.

    Separate from `evaluate_target` because a caller needs the antigen
    sequences before the loop - the natural decoy for one target is another
    target's antigen, a non-cognate pair of real proteins rather than a
    shuffled sequence that is not a protein at all.
    """
    from .api import get_structure
    from .structures import chain_sequence, relabel_chains

    root = Path(output_dir) / target.pdb_id.upper()
    reference = get_structure(
        root / "reference",
        model_name="pdb_download",
        pdb_id=None if reference_source else target.pdb_id,
        source=str(reference_source) if reference_source else None,
        chains=list(target.chains),
        cache=cache,
        confidence=False,
    )
    normalized = relabel_chains(
        reference.path, target.mapping, root / "reference" / "normalized.pdb"
    )
    return normalized, (
        chain_sequence(normalized, "H"),
        chain_sequence(normalized, "L"),
        chain_sequence(normalized, "A"),
    )


def evaluate_target(
    target: Target,
    *,
    output_dir: Path | str,
    model_name: str = "boltz2",
    num_seeds: int = 15,
    decoy_antigen_seq: str | None = None,
    reference_source: str | Path | None = None,
    cache: bool = True,
    with_dockq: bool = True,
    **kwargs: Any,
) -> dict[str, Any]:
    """Predict one complex from sequence and score it against the deposited one.

    The protocol: take the deposited structure, read the three chain sequences
    off it, throw the coordinates away, predict, compare. The reference is
    relabelled to H/L/A so both files speak the same chain ids.

    Returns one flat row - the same shape `metrics.save` takes - and writes it
    as `metrics.json` under the target's directory.
    """
    from .api import get_structure

    root = Path(output_dir) / target.pdb_id.upper()
    normalized, (heavy, light, antigen) = reference_structure(
        target, output_dir=output_dir, reference_source=reference_source, cache=cache
    )

    predicted = get_structure(
        root / "predicted",
        model_name=model_name,
        antibody_seq=f"{heavy}:{light}",
        antigen_seq=antigen,
        chain_ids=("H", "L", "A"),
        num_seeds=num_seeds,
        decoy_antigen_seq=decoy_antigen_seq,
        cache=cache,
        **kwargs,
    )

    row: dict[str, Any] = {
        "pdb_id": target.pdb_id.upper(),
        "tier": target.tier,
        "date": target.date,
        "resolution": target.resolution,
        "model_name": model_name,
        "n_seeds": len(predicted.paths),
        "interface_rmsd": metrics.interface_rmsd(predicted.path, normalized),
        "epitope_recall": metrics.epitope_recall(
            predicted.path, normalized, ["A"], ["H", "L"]
        ),
    }
    if with_dockq:
        try:
            row.update(
                {
                    f"dockq_{key}" if key != "dockq" else key: value
                    for key, value in metrics.dockq(predicted.path, normalized).items()
                }
            )
        except MissingDependencyError as exc:
            row["dockq"] = float("nan")
            row["dockq_note"] = str(exc)[:120]
    row.update(
        {
            f"confidence_{key}": value
            for key, value in (predicted.confidence or {}).items()
            if isinstance(value, (float, int, str, bool, type(None)))
        }
    )
    metrics.save(root, {k: v for k, v in row.items() if _scalar(v)})
    return row


def _scalar(value: Any) -> bool:
    return isinstance(value, (float, int, str, bool, type(None)))


def fixture_targets(
    output_dir: Path | str, count: int = 3
) -> list[tuple[Target, Path]]:
    """Synthetic (target, reference) pairs, for a run with no GPU and no network.

    The validation notebook needs to execute in CI, where there is no predictor
    and no RCSB. These stand in: deposited-looking complexes with non-H/L/A
    chain labels, a real interface, and a different antigen in each one so the
    ensemble and the decoy control come out differently per target.

    The numbers they produce are meaningless and the notebook says so. What
    they exercise - the selection, the relabelling, the prediction, every
    metric, the calibration plot - is the real code path.
    """
    from .adapters._pdb import (
        ANTIBODY_SPACING,
        DOCKED,
        FIXTURE_CHAINS,
        write_ca_pdb,
    )
    from .structures import relabel_chains

    root = Path(output_dir)
    antigen = FIXTURE_CHAINS["A"]
    made: list[tuple[Target, Path]] = []
    for index in range(count):
        # A different antigen per target: shorter by a few residues each time,
        # which changes its hash and so its ipTM, its poses and its decoy gap.
        chains = {
            "H": FIXTURE_CHAINS["H"],
            "L": FIXTURE_CHAINS["L"],
            "A": antigen[: len(antigen) - 4 * index],
        }
        origins = {
            "H": (0.0, 0.0, 0.0),
            "L": (ANTIBODY_SPACING, 0.0, 0.0),
            "A": DOCKED,
        }
        raw = write_ca_pdb(root / f"fixture_{index}" / "raw.pdb", chains, origins=origins)
        # Relabelled the way a deposited entry is, so the notebook exercises
        # the chain mapping rather than skipping it.
        deposited = relabel_chains(
            raw, {"H": "B", "L": "A", "A": "C"}, root / f"fixture_{index}" / "deposited.pdb"
        )
        made.append(
            (
                Target(
                    pdb_id=f"FIX{index}",
                    heavy="B",
                    light="A",
                    antigen="C",
                    tier="fixture",
                    date="1970-01-01",
                    note="synthetic: no GPU, no network, meaningless numbers",
                ),
                deposited,
            )
        )
    return made
