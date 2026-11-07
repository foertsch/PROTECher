"""Gate checks: a boolean and a reason per candidate.

Same Protocol as everything else, same candidates CSV in, a `checks.csv` out.
A failure is dropped by `gate`, never traded off against binding score - which
is the point of having them separate from the scorers.

Gate only on confident failures. The developability predictors are weak, so a
strict threshold on a marginal model score throws away good candidates; the
reliable signal is a hard sequence liability. The thresholds below are
deliberately loose, they are first guesses rather than calibrated numbers, and
every one of them is overridable per call
(`gate(..., threshold=2.0)`) so the harness can settle them with data.
"""

from __future__ import annotations

import csv
import os
import re
from dataclasses import dataclass, replace
from pathlib import Path
from typing import Any, Mapping, Sequence

from ..containers import ContainerAdapter, ProcessAdapter, make_id
from ..designs import Design
from ..errors import MissingRuntimeError
from ..routes import require
from ..structures import chain_residue_ids, chain_sequence
from ._tabular import CANDIDATES, CHECKS, read_candidates, read_tool_table, write_checks

# Each entry is (label, pattern). The pattern is wrapped in a lookahead when it
# is searched, so overlapping occurrences are all found: NNS holds both an
# N-glycosylation sequon and a deamidation motif.
LIABILITIES: tuple[tuple[str, str], ...] = (
    # N-x-S/T with x not proline: an N-linked glycosylation sequon. The most
    # common own-goal in an affinity maturation campaign.
    ("n_glycosylation", r"N[^P][ST]"),
    ("deamidation", r"N[GST]"),
    ("isomerization", r"D[GSD]"),
    ("fragmentation", r"DP"),
)


def new_motifs(wild: str, mutant: str) -> list[str]:
    """Motifs present in `mutant` and not in `wild`, as "label@position".

    Occurrence-wise rather than count-wise: a motif that moves along the chain
    is a new liability at its new position, and one that was already there is
    not this variant's fault.
    """
    found: list[str] = []
    for label, pattern in LIABILITIES:
        lookahead = f"(?={pattern})"
        before = {m.start() for m in re.finditer(lookahead, wild)}
        after = {m.start() for m in re.finditer(lookahead, mutant)}
        found += [f"{label}@{start + 1}" for start in sorted(after - before)]
    return found


@dataclass(frozen=True)
class LiabilityCheck:
    """Pure sequence regex. Zero dependencies, milliseconds, no structure needed
    beyond the chain sequences the complex already holds.

    Cheap enough to run on the full scan, which is why it is the default check:
    it catches the mistakes that no amount of binding score makes up for.
    """

    name: str = "liabilities"
    per_variant: bool = False
    requires_structure: bool = False
    """Reads the chain sequences and never the coordinates, so it is just as
    valid on a complex whose docking is wrong."""

    def run(
        self,
        inputs: Sequence[Design],
        output_dir: Path,
        *,
        candidates: Path | str | None = None,
        extra_args: Sequence[str] | None = None,
        **kwargs: Any,
    ) -> list[Design]:
        if candidates is None:
            raise ValueError(f"{self.name} needs candidates=<path to {CANDIDATES}>")
        if not inputs:
            raise ValueError(f"{self.name} needs the complex as its input")
        structure = inputs[0].path
        wanted = read_candidates(Path(candidates))
        wild: dict[str, str] = {}

        rows: list[tuple[str, bool, str]] = []
        for variant_id, mutations in wanted.items():
            reasons: list[str] = []
            for chain in sorted({chain for chain, _, _, _ in mutations}):
                if chain not in wild:
                    wild[chain] = chain_sequence(structure, chain)
                residues = list(wild[chain])
                for mut_chain, position, _wt, mut in mutations:
                    if mut_chain == chain:
                        residues[position - 1] = mut
                reasons += [
                    f"{chain}:{motif}"
                    for motif in new_motifs(wild[chain], "".join(residues))
                ]
            # An introduced cysteine has nothing to pair with, so it is a free
            # thiol: aggregation, and a heterogeneous product. Whether an
            # existing cysteine is in a disulfide is a structure question this
            # check deliberately does not ask.
            reasons += [
                f"{chain}:free_cysteine@{position}"
                for chain, position, _wt, mut in mutations
                if mut == "C"
            ]
            rows.append((variant_id, not reasons, ", ".join(sorted(reasons))))

        path = write_checks(output_dir / CHECKS, rows)
        return [
            Design(
                id=make_id(path, output_dir),
                path=path,
                parent_id=inputs[0].id,
                producer=self.name,
            )
        ]


def threshold_rows(
    name: str,
    wanted: Mapping[str, object],
    values: Mapping[str, float],
    column: str,
    limit: float,
    fail_above: bool,
) -> list[tuple[str, bool, str]]:
    """One (variant_id, passed, reason) per candidate, in candidate order.

    A candidate with no value survives. Gate only on confident failures, and
    "the tool said nothing about this one" is not one of those - the reason
    records that nothing was measured rather than dropping it.
    """
    rows: list[tuple[str, bool, str]] = []
    for variant_id in wanted:
        value = values.get(variant_id)
        if value is None:
            rows.append((variant_id, True, f"{name}: no value"))
            continue
        failed = value > limit if fail_above else value < limit
        side = ">" if fail_above else "<"
        rows.append(
            (
                variant_id,
                not failed,
                f"{name}: {column}={value:.3f} {side} {limit:g}" if failed else "",
            )
        )
    return rows


@dataclass(frozen=True)
class ContainerCheck:
    """A developability container, reduced to a threshold on one column.

    `fail_above` says which side of `threshold` is the failure. The reason
    string carries the number, so a gate that drops too much can be re-tuned
    from the annotated table without re-running the tool.
    """

    name: str
    value_column: str
    threshold: float
    tool: ContainerAdapter | ProcessAdapter | None = None
    """None for a check this package cannot invoke: no published image and no
    verified command line. Not a guess at one - `run` raises out of `routes`."""
    fail_above: bool = True
    per_variant: bool = False
    requires_structure: bool = True
    id_column: str = "variant_id"

    def run(
        self,
        inputs: Sequence[Design],
        output_dir: Path,
        *,
        candidates: Path | str | None = None,
        extra_args: Sequence[str] | None = None,
        threshold: float | None = None,
        id_column: str | None = None,
        value_column: str | None = None,
        **kwargs: Any,
    ) -> list[Design]:
        if candidates is None:
            raise ValueError(f"{self.name} needs candidates=<path to {CANDIDATES}>")
        limit = self.threshold if threshold is None else float(threshold)
        column = value_column or self.value_column
        if self.tool is None:
            require(self.name)
            raise MissingRuntimeError(
                f"{self.name} is installed, but this package has no verified "
                f"way to invoke it. Set REGISTRY[{self.name!r}]; "
                f"`python -m design_pipeline.preflight` prints the table."
            )
        produced = self.tool.run(
            inputs,
            output_dir,
            candidates=str(candidates),
            extra_args=extra_args,
            **kwargs,
        )
        values = read_tool_table(
            [d.path for d in produced], id_column or self.id_column, column
        )
        wanted = read_candidates(Path(candidates))

        rows = threshold_rows(
            self.name, wanted, values, column, limit, self.fail_above
        )
        path = write_checks(output_dir / CHECKS, rows)
        return [
            Design(
                id=make_id(path, output_dir),
                path=path,
                parent_id=inputs[0].id if inputs else None,
                producer=self.name,
            )
        ]


def check_container(image: str, *, gpu: bool = True) -> ContainerAdapter:
    """A gate-check image, for one you have built yourself.

    Public for the same reason as `_scorers.matrix_container`: this package
    ships no image names, so a check that needs a container takes it from you.
    The flags are this package's convention, not any tool's.
    """
    return ContainerAdapter(
        image=image,
        input_template=("--pdb", "{path}"),
        output_template=("--out_dir", "{path}"),
        aliases={"candidates": "mutations_csv"},
        path_kwargs=("candidates",),
        output_globs=("*.csv",),
        gpu=gpu,
    )


INFERENCE_SCRIPT = ("analysis", "custom_inference.py")
"""Where ThermoMPNN's prediction script sits inside the clone."""

# ThermoMPNN ships no console script - its pyproject sets `package = false` -
# so there is nothing to put on PATH and nothing for `-m` to import. The
# adapter runs the script out of the clone with the interpreter of whichever
# environment it was installed into, uv or mamba. The script path is
# per-install, so `command` is filled in at call time from THERMOMPNN_SRC
# rather than frozen here.
thermompnn_process = ProcessAdapter(
    tool="thermompnn",
    executable="python",
    env_var="DESIGN_PIPELINE_THERMOMPNN",
    input_template=("--pdb", "{path}"),
    output_template=("--out_dir", "{path}"),
    output_globs=("ThermoMPNN_inference_*.csv",),
)


def thermompnn_source() -> Path:
    """The ThermoMPNN clone, from THERMOMPNN_SRC."""
    found = os.environ.get("THERMOMPNN_SRC")
    if not found:
        raise MissingRuntimeError(
            "thermompnn needs THERMOMPNN_SRC (or src=): the ThermoMPNN clone, "
            f"the folder holding {'/'.join(INFERENCE_SCRIPT)}. The "
            "interpreter in $DESIGN_PIPELINE_THERMOMPNN cannot import it - the "
            "repo sets `package = false`, so there is nothing on its path to "
            "find. `python -m design_pipeline.preflight` prints the table."
        )
    return Path(found)


def thermompnn_index(labels: Sequence[str]) -> dict[int, int]:
    """Chain position (1-based, file order) -> ThermoMPNN's sequence index.

    ThermoMPNN does not speak the file's residue numbers, and its output is not
    in chain positions either. `alt_parse_PDB` keys every residue by
    `int(resSeq) - 1`, walks `range(min, max + 1)` and emits a '-' for each
    number in that span with no residue, so the `position` column is an index
    into that gap-padded string: the author number less the lowest author
    number present, shifted along by every insertion code before it.

    Reconstructed from the labels rather than applied as one offset, because an
    antibody chain has both - gaps where a loop was not modelled, and 100A/100B
    insertions in CDR H3 - and an offset alone would be wrong after the first
    of either.
    """
    parsed: list[tuple[int, str]] = []
    for label in labels:
        text = str(label).strip()
        code = text[-1] if text and text[-1].isalpha() else ""
        parsed.append((int(text[: len(text) - len(code)] if code else text), code))

    by_number: dict[int, list[str]] = {}
    for number, code in parsed:
        by_number.setdefault(number, []).append(code)

    index_of: dict[tuple[int, str], int] = {}
    cursor = 0
    for number in range(min(by_number), max(by_number) + 1):
        if number not in by_number:
            cursor += 1  # the '-' standing in for a number with no residue
            continue
        for code in sorted(by_number[number]):
            index_of[(number, code)] = cursor
            cursor += 1
    return {position + 1: index_of[key] for position, key in enumerate(parsed)}


def read_ssm_table(
    paths: Sequence[Path], value_column: str = "ddG_pred"
) -> dict[tuple[int, str], tuple[float, str]]:
    """ThermoMPNN's own CSV as (index, mutant) -> (value, wild type).

    The wild-type letter comes along because it is the only thing in that file
    that can tell us the numbering lined up, which is what the caller checks.
    """
    found: dict[tuple[int, str], tuple[float, str]] = {}
    needed = {value_column, "position", "wildtype", "mutation"}
    for path in paths:
        if path.suffix.lower() != ".csv":
            continue
        with path.open(newline="", encoding="utf-8", errors="replace") as handle:
            reader = csv.DictReader(handle)
            if not needed <= set(reader.fieldnames or ()):
                continue
            for row in reader:
                try:
                    found[(int(row["position"]), row["mutation"])] = (
                        float(row[value_column]),
                        row["wildtype"],
                    )
                except (TypeError, ValueError):
                    continue
    if not found:
        names = ", ".join(str(p.name) for p in paths) or "nothing"
        raise ValueError(
            f"no rows carrying {value_column!r}, 'position', 'wildtype' and "
            f"'mutation' in {names}. Pass value_column= if this version of "
            f"ThermoMPNN names its prediction column differently."
        )
    return found


@dataclass(frozen=True)
class ThermoMPNNCheck:
    """Predicted folding ddG, destabilizing positive, from site saturation.

    Its own contract rather than `ContainerCheck`'s, because the tool does not
    fit that one. `analysis/custom_inference.py` takes a PDB and a single
    chain, mutates every position of that chain to all twenty residues, and
    writes a table keyed by its own sequence index. So: one call per chain
    rather than one per candidate, the candidates joined onto that table
    afterwards, and no value at all for a multi-mutation candidate, since the
    script cannot be asked about a combination.

    The join is checked rather than trusted. ThermoMPNN reports the wild-type
    residue it thinks is at each index, and a disagreement with the candidate
    table means the numbering did not line up, which `run` raises on instead of
    passing or failing candidates on values that belong to other positions.

    1.5 kcal/mol is loose for the same reason as the rest of the gate: a
    marginal prediction is not a confident failure.
    """

    name: str = "thermompnn"
    value_column: str = "ddG_pred"
    threshold: float = 1.5
    tool: ProcessAdapter | None = None
    fail_above: bool = True
    per_variant: bool = True
    requires_structure: bool = True
    id_column: str = "variant_id"

    def run(
        self,
        inputs: Sequence[Design],
        output_dir: Path,
        *,
        candidates: Path | str | None = None,
        src: Path | str | None = None,
        extra_args: Sequence[str] | None = None,
        threshold: float | None = None,
        value_column: str | None = None,
        **kwargs: Any,
    ) -> list[Design]:
        if candidates is None:
            raise ValueError(f"{self.name} needs candidates=<path to {CANDIDATES}>")
        # Asked before anything is read, so a missing install arrives as its
        # install route rather than as a parse error on someone's structure.
        require(self.name)
        if self.tool is None:
            raise MissingRuntimeError(
                f"{self.name} is installed, but this package has no verified "
                f"way to invoke it. Set REGISTRY[{self.name!r}]; "
                f"`python -m design_pipeline.preflight` prints the table."
            )
        limit = self.threshold if threshold is None else float(threshold)
        column = value_column or self.value_column
        script = (Path(src) if src is not None else thermompnn_source()).joinpath(
            *INFERENCE_SCRIPT
        )
        if not script.is_file():
            raise MissingRuntimeError(
                f"{script} does not exist. THERMOMPNN_SRC (or src=) has to be "
                f"the clone itself, the folder holding "
                f"{'/'.join(INFERENCE_SCRIPT)}."
            )
        if not inputs:
            raise ValueError(
                f"{self.name} needs the structure its candidates are indexed "
                f"against, as inputs[0]"
            )

        structure = inputs[0].path
        wanted = read_candidates(Path(candidates))
        # One mutation only: site saturation says nothing about a combination,
        # and summing two of its numbers would be inventing an answer.
        singles = {
            variant_id: mutations[0]
            for variant_id, mutations in wanted.items()
            if len(mutations) == 1
        }

        runner = replace(self.tool, command=(str(script),))
        values: dict[str, float] = {}
        for chain in sorted({mutation[0] for mutation in singles.values()}):
            # Per chain, because the script takes one --chain and writes the
            # same filename every time: a shared directory would look like a
            # run that produced nothing the second time around.
            chain_dir = output_dir / f"chain_{chain}"
            chain_dir.mkdir(parents=True, exist_ok=True)
            produced = runner.run(
                inputs[:1],
                chain_dir,
                extra_args=extra_args,
                chain=chain,
                **kwargs,
            )
            table = read_ssm_table([d.path for d in produced], column)
            index_of = thermompnn_index(chain_residue_ids(structure, chain))
            for variant_id, (its_chain, position, wt, mut) in singles.items():
                if its_chain != chain:
                    continue
                answer = table.get((index_of.get(position, -1), mut))
                if answer is None:
                    continue
                value, wildtype = answer
                if wildtype != wt:
                    raise ValueError(
                        f"{self.name} disagrees about position {position} of "
                        f"chain {chain}: the candidates say the wild type is "
                        f"{wt!r}, ThermoMPNN's index "
                        f"{index_of.get(position)} says {wildtype!r}. Either "
                        f"the structure handed to the gate is not the one the "
                        f"candidates were built against, or this version "
                        f"numbers its output differently."
                    )
                values[variant_id] = value

        rows = threshold_rows(self.name, wanted, values, column, limit, self.fail_above)
        path = write_checks(output_dir / CHECKS, rows)
        return [
            Design(
                id=make_id(path, output_dir),
                path=path,
                parent_id=inputs[0].id,
                producer=self.name,
            )
        ]


liabilities = LiabilityCheck()

thermompnn = ThermoMPNNCheck(tool=thermompnn_process)
# Aggregation propensity, structure-based, so one model per candidate.
aggrescan3d = ContainerCheck(
    name="aggrescan3d",
    tool=None,
    value_column="a3d_score",
    threshold=1.0,
    fail_above=True,
    per_variant=True,
)
# Predicted solubility, higher is better, hence fail_above=False.
netsolp = ContainerCheck(
    name="netsolp",
    tool=None,
    value_column="solubility",
    threshold=0.3,
    fail_above=False,
)


@dataclass(frozen=True)
class DummyCheck:
    """Offline fixture: fails the candidates whose id hashes into `fail_rate`.

    Deterministic, so a gated table is reproducible, and enough to exercise the
    fan-out, the join and the dropping.
    """

    name: str = "dummy_check"
    per_variant: bool = False
    requires_structure: bool = True
    fail_rate: float = 0.25

    def run(
        self,
        inputs: Sequence[Design],
        output_dir: Path,
        *,
        candidates: Path | str | None = None,
        extra_args: Sequence[str] | None = None,
        **kwargs: Any,
    ) -> list[Design]:
        import hashlib

        if candidates is None:
            raise ValueError(f"{self.name} needs candidates=<path to {CANDIDATES}>")
        rows: list[tuple[str, bool, str]] = []
        for variant_id in read_candidates(Path(candidates)):
            digest = hashlib.sha256(f"{self.name}:{variant_id}".encode()).digest()
            unit = int.from_bytes(digest[:8], "big") / 2**64
            failed = unit < self.fail_rate
            rows.append(
                (variant_id, not failed, f"{self.name}: {unit:.3f}" if failed else "")
            )
        path = write_checks(output_dir / CHECKS, rows)
        return [
            Design(
                id=make_id(path, output_dir),
                path=path,
                parent_id=inputs[0].id if inputs else None,
                producer=self.name,
            )
        ]


dummy_check = DummyCheck()

CHECKS_AVAILABLE: tuple[str, ...] = (
    "liabilities",
    "thermompnn",
    "aggrescan3d",
    "netsolp",
    "dummy_check",
)
DEFAULT_CHECKS: tuple[str, ...] = ("liabilities",)
"""What `gate` runs when nothing is named: the free one that catches the most
common own-goal. The container checks are opt-in: name them explicitly."""
