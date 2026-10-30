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

import re
from dataclasses import dataclass
from pathlib import Path
from typing import Any, Sequence

from ..containers import ContainerAdapter, ProcessAdapter, make_id
from ..designs import Design
from ..errors import MissingRuntimeError
from ..routes import require
from ..structures import chain_sequence
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

        rows: list[tuple[str, bool, str]] = []
        for variant_id in wanted:
            value = values.get(variant_id)
            if value is None:
                # No answer is not a confident failure. The candidate survives
                # and the reason records that nothing was measured.
                rows.append((variant_id, True, f"{self.name}: no value"))
                continue
            failed = value > limit if self.fail_above else value < limit
            side = ">" if self.fail_above else "<"
            rows.append(
                (
                    variant_id,
                    not failed,
                    f"{self.name}: {column}={value:.3f} {side} {limit:g}"
                    if failed
                    else "",
                )
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


# ThermoMPNN installs from its repository into its own mamba environment - its
# environment.yaml, which its README warns can resolve a CPU-only torch - so
# the adapter runs that environment's interpreter rather than anything on PATH.
thermompnn_process = ProcessAdapter(
    tool="thermompnn",
    executable="thermompnn",
    env_var="DESIGN_PIPELINE_THERMOMPNN",
    input_template=("--pdb", "{path}"),
    output_template=("--out_dir", "{path}"),
    aliases={"candidates": "mutations_csv"},
    path_kwargs=("candidates",),
    output_globs=("*.csv",),
)


liabilities = LiabilityCheck()

# Predicted change in folding stability, destabilizing positive. 1.5 kcal/mol
# is loose on purpose: a marginal prediction is not a confident failure.
thermompnn = ContainerCheck(
    name="thermompnn",
    tool=thermompnn_process,
    value_column="ddg",
    threshold=1.5,
    fail_above=True,
    per_variant=True,
)
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
