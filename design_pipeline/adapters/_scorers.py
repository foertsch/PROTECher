"""Scorers: one value per candidate, written to a CSV, joined on `variant_id`.

They satisfy the one existing adapter Protocol. There is no second Protocol for
scoring and there is no scoring base class; what a scorer adds is two pieces of
data on the instance:

* `per_variant` - its cost class. A matrix scorer runs once on the complex and
  reads every single mutation off one forward pass, so the whole 4,598-row
  table costs about as much as one candidate. A per-variant scorer builds and
  relaxes a mutant structure per candidate, seconds to minutes each, so
  `score_variants` refuses to start one on an uncapped table.
* `sign` - which way the tool's numbers point. Normalization happens here, at
  the adapter boundary: `score` always means higher-is-better-binding, and the
  tool's own number is kept in `raw` so a convention we got backwards is
  recoverable without re-running anything.

On pairs. Every scorer here is handed the candidate's *whole* mutation set and
answers for the mutant as a whole. A probability matrix cannot do that for a
double mutant - the lookup is exact only for a single - so a matrix scorer
either conditions on the full mutant sequence (what the inverse-folding models
do, one pass per multi-mutation candidate) or returns no score. What it must
never do is add two single scores together, so nothing here does.
"""

from __future__ import annotations

import hashlib
import math
import os
from dataclasses import dataclass
from pathlib import Path
from typing import Any, Mapping, Sequence

from ..containers import ContainerAdapter, ProcessAdapter, make_id
from ..designs import Design
from ..errors import MissingDependencyError, MissingRuntimeError, RunFailedError
from ..routes import require
from ._tabular import (
    CANDIDATES,
    SCORES,
    read_candidates,
    read_tool_table,
    write_scores,
)

MATRIX = "matrix"
PER_VARIANT = "per_variant"


@dataclass(frozen=True)
class ContainerScorer:
    """An external scoring tool plus the three facts needed to read it back.

    The tool itself is an ordinary `ContainerAdapter` or `ProcessAdapter` -
    same argument assembly either way - configured with the candidates CSV as a
    path kwarg so it is mounted and rewritten like any other input file. This
    class only normalizes the output.

    `tool` is None for a scorer this package cannot yet invoke: no published
    image and no verified command line. It is not a guess at one. `run` then
    raises out of `routes`, naming the install route if there is one and saying
    there is none if there is not.
    """

    name: str
    tool: ContainerAdapter | ProcessAdapter | None = None
    sign: float = 1.0
    """+1 when the tool already scores higher-is-better-binding, -1 for a
    ddG-style number where negative is the improvement."""
    per_variant: bool = False
    requires_structure: bool = True
    """Every container here reads coordinates. False only for the scorers whose
    numbers would be unchanged by a wrong structure - see _plm.py."""
    id_column: str = "variant_id"
    value_column: str = "score"

    def run(
        self,
        inputs: Sequence[Design],
        output_dir: Path,
        *,
        candidates: Path | str | None = None,
        extra_args: Sequence[str] | None = None,
        id_column: str | None = None,
        value_column: str | None = None,
        **kwargs: Any,
    ) -> list[Design]:
        if candidates is None:
            raise ValueError(f"{self.name} needs candidates=<path to {CANDIDATES}>")
        if self.tool is None:
            # Raises when the tool is absent, naming the route. If it is
            # present, the gap is ours rather than the environment's, and
            # saying which is the whole point of this branch.
            require(self.name)
            raise MissingRuntimeError(
                f"{self.name} is installed, but this package has no verified "
                f"way to invoke it - it exposes a Python API rather than a CLI, "
                f"or its result format has not been checked. Set "
                f"REGISTRY[{self.name!r}] to an adapter that runs it; "
                f"`python -m design_pipeline.preflight` prints the table."
            )
        produced = self.tool.run(
            inputs,
            output_dir,
            candidates=str(candidates),
            extra_args=extra_args,
            **kwargs,
        )
        raw = read_tool_table(
            [d.path for d in produced],
            id_column or self.id_column,
            value_column or self.value_column,
        )
        wanted = read_candidates(Path(candidates))
        unknown = sorted(set(raw) - set(wanted))
        if unknown:
            raise RunFailedError(
                f"{self.name} scored {len(unknown)} ids that were not asked "
                f"for (e.g. {unknown[0]!r}); the join key convention has "
                f"drifted. Expected the ids in {CANDIDATES}."
            )
        path = write_scores(
            output_dir / SCORES,
            (
                (variant_id, _signed(raw.get(variant_id), self.sign), raw.get(variant_id))
                for variant_id in wanted
            ),
        )
        return [
            Design(
                id=make_id(path, output_dir),
                path=path,
                parent_id=inputs[0].id if inputs else None,
                producer=self.name,
            )
        ]


def _signed(value: float | None, sign: float) -> float | None:
    return None if value is None else sign * value


@dataclass(frozen=True)
class ConsensusScorer:
    """Positional consensus from a cached OAS frequency table. No model, no GPU.

    The score is log(f_mut / f_wt) at the position: how much more often the
    repertoire puts the replacement there than the residue that is there now.
    Weak, free, and uncorrelated with the structure-based models, which is the
    whole reason to have it in the default set.

    A double mutant gets no score. Positional frequencies hold one marginal per
    position and combining two of them is multiplying marginals, which is the
    additivity shortcut by another name; an empty cell is the honest answer.

    The table is not shipped: OAS is a download. Point `frequencies` at a CSV
    of `imgt,aa,frequency` (or `position,aa,frequency` when ANARCI is not
    available), or set DESIGN_PIPELINE_OAS_TABLE.
    """

    name: str = "consensus"
    per_variant: bool = False
    requires_structure: bool = False
    """A positional prior. It never opens the structure, so it survives a
    prediction that put the paratope on the wrong patch."""
    frequencies: Path | None = None
    pseudocount: float = 1e-4
    """Keeps log() finite for a residue the repertoire never uses at a
    position, which is information rather than an error."""

    def run(
        self,
        inputs: Sequence[Design],
        output_dir: Path,
        *,
        candidates: Path | str | None = None,
        frequencies: Path | str | None = None,
        extra_args: Sequence[str] | None = None,
        **kwargs: Any,
    ) -> list[Design]:
        if candidates is None:
            raise ValueError(f"{self.name} needs candidates=<path to {CANDIDATES}>")
        table = _frequency_table(frequencies or self.frequencies)
        wanted = read_candidates(Path(candidates))

        rows: list[tuple[str, float | None, float | None]] = []
        for variant_id, mutations in wanted.items():
            if len(mutations) != 1:
                rows.append((variant_id, None, None))
                continue
            _chain, position, wt, mut = mutations[0]
            key = str(position)
            wt_frequency = table.get((key, wt))
            mut_frequency = table.get((key, mut))
            if wt_frequency is None or mut_frequency is None:
                rows.append((variant_id, None, None))
                continue
            score = math.log(
                (mut_frequency + self.pseudocount) / (wt_frequency + self.pseudocount)
            )
            # Already higher-is-better, so raw and normalized agree.
            rows.append((variant_id, score, score))

        path = write_scores(output_dir / SCORES, rows)
        return [
            Design(
                id=make_id(path, output_dir),
                path=path,
                parent_id=inputs[0].id if inputs else None,
                producer=self.name,
            )
        ]


def _frequency_table(source: Path | str | None) -> Mapping[tuple[str, str], float]:
    """(position-or-imgt label, residue) -> frequency, from a CSV.

    Keyed on whichever of `imgt` and `position` the file has. IMGT is the one
    that transfers between antibodies; `position` is the fallback for a run
    without ANARCI, and then the table only means anything for the antibody it
    was built from.
    """
    import csv

    configured = str(source or os.environ.get("DESIGN_PIPELINE_OAS_TABLE", "")).strip()
    if not configured:
        raise MissingDependencyError(
            "the consensus scorer needs an OAS frequency table and none is "
            "configured. Pass frequencies=<path>, or set "
            "DESIGN_PIPELINE_OAS_TABLE. The file is a CSV of imgt,aa,frequency "
            "(or position,aa,frequency), built once from the Observed Antibody "
            "Space download at https://opig.stats.ox.ac.uk/webapps/oas/. "
            "Drop 'consensus' from models=[...] to run without it."
        )
    path = Path(configured)
    if not path.is_file():
        raise MissingDependencyError(
            f"the consensus frequency table {path} does not exist. Build it "
            f"from OAS, or drop 'consensus' from models=[...]."
        )
    table: dict[tuple[str, str], float] = {}
    with path.open(newline="", encoding="utf-8") as handle:
        reader = csv.DictReader(handle)
        fields = reader.fieldnames or []
        key_column = "imgt" if "imgt" in fields else "position"
        if key_column not in fields or "aa" not in fields or "frequency" not in fields:
            raise ValueError(
                f"{path} needs columns (imgt|position),aa,frequency; found "
                f"{', '.join(fields) or 'none'}"
            )
        for row in reader:
            try:
                table[(str(row[key_column]), str(row["aa"]))] = float(row["frequency"])
            except (TypeError, ValueError):
                continue
    return table


@dataclass(frozen=True)
class DummyScorer:
    """Deterministic pseudo-random scores keyed on `variant_id`.

    So the whole pipeline runs end to end on a laptop with no GPU, no container
    runtime and no network. Keyed on the full id, which means a pair gets its
    own independent number rather than anything derived from its halves - the
    fixture obeys the same no-additivity rule as the real scorers.

    `per_variant` is a constructor argument because the cap that
    `score_variants` enforces on expensive scorers has to be testable offline.
    """

    name: str = "dummy_scorer"
    per_variant: bool = False
    requires_structure: bool = True
    """True by default because the fixture stands in for the structure-based
    scorers; `dummy_sequence_scorer` is the structure-free counterpart."""
    spread: float = 5.0
    """Scores land in [-spread, +spread], the rough range of a ddG model."""

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
        wanted = read_candidates(Path(candidates))
        rows = []
        for variant_id in wanted:
            digest = hashlib.sha256(f"{self.name}:{variant_id}".encode()).digest()
            unit = int.from_bytes(digest[:8], "big") / 2**64
            raw = self.spread * (2.0 * unit - 1.0)
            rows.append((variant_id, raw, raw))
        path = write_scores(output_dir / SCORES, rows)
        return [
            Design(
                id=make_id(path, output_dir),
                path=path,
                parent_id=inputs[0].id if inputs else None,
                producer=self.name,
            )
        ]


def matrix_container(
    image: str,
    *,
    mutations_flag: str = "mutations_csv",
    gpu: bool = True,
    fixed_args: tuple[str, ...] = (),
    command: tuple[str, ...] = (),
) -> ContainerAdapter:
    """A scoring image: structure in, candidates CSV in, CSV out.

    Public because this package no longer ships image names - every scorer
    below that needs a container has `tool=None` - so building one is the
    caller's job:

        REGISTRY["gearbind"] = replace(
            _scorers.gearbind,
            tool=_scorers.matrix_container("ourlab/gearbind:1.0"),
        )

    `per_input=True` with exactly one input - the structure - so the container
    layer mounts the structure's directory, the candidates CSV's directory and
    the run directory, and nothing else. The flags below are this package's
    convention, not any tool's: override them for the image you built.
    """
    return ContainerAdapter(
        image=image,
        command=command,
        input_template=("--pdb", "{path}"),
        output_template=("--out_dir", "{path}"),
        fixed_args=fixed_args,
        aliases={"candidates": mutations_flag},
        path_kwargs=("candidates",),
        output_globs=("*.csv",),
        gpu=gpu,
    )


# AntiFold cannot live in Env A. `pip install antifold` pins biotite==0.38,
# torch==2.3.1 and numpy==1.26.*, and this package needs biotite>=2.0; pip calls
# that combination ResolutionImpossible. So it runs out of its own environment
# and `executable=` (or DESIGN_PIPELINE_ANTIFOLD) is how you point at it.
#
# Its own command line (`python -m antifold.main --pdb_file --out_dir`) cannot
# be called as a scorer. It takes no candidates CSV, so the `--candidates` this
# adapter always forwards is an argparse error, and what it writes is a
# per-residue matrix, not the `variant_id,score` table `read_tool_table` reads.
# `antifold_score` (a module at the repo root, run inside the AntiFold
# environment) closes both gaps: it takes the candidates CSV, calls AntiFold's
# own functions and writes `scores.csv`. It is started through `-c` with the
# repo root put on the path, so it is found whatever the working directory.
#
# Needs, besides DESIGN_PIPELINE_ANTIFOLD: ANTIFOLD_SRC, the AntiFold source
# folder that holds models/model.pt.
_ANTIFOLD_SCORE = (
    "import runpy, sys; "
    f"sys.path.insert(0, {str(Path(__file__).resolve().parents[2])!r}); "
    "runpy.run_module('antifold_score', run_name='__main__', alter_sys=True)"
)


antifold_process = ProcessAdapter(
    tool="antifold",
    executable="python",
    command=("-c", _ANTIFOLD_SCORE),
    env_var="DESIGN_PIPELINE_ANTIFOLD",
    input_template=("--pdb", "{path}"),
    output_template=("--out_dir", "{path}"),
    fixed_args=("--route", "table"),
    aliases={"candidates": "mutations_csv"},
    path_kwargs=("candidates",),
    output_globs=("scores.csv",),
)


# Matrix scorers. One forward pass over the complex, then a lookup per single
# mutation. The lookup is exact for a single only, so a candidate with several
# mutations gets no score here (an empty cell and a note, never a sum); score
# those with `antifold_seq`, which conditions on the whole mutant.
antifold = ContainerScorer(
    name="antifold",
    tool=antifold_process,
    # AntiFold reports a log-likelihood ratio against the wild type: already
    # higher-is-better.
    sign=1.0,
)

# The whole-sequence counterpart of `antifold`. The matrix scorer above is exact
# for a single mutation and has no answer for several; this one asks AntiFold's
# own score_sequence_in_complex how likely the mutant chain is given the
# backbone, one forward pass per mutated chain, so a multi-mutation candidate is
# scored as a whole and never as a sum. It disagrees with the lookup on a single
# mutation (20 trastuzumab singles: correlation 0.81), so it is its own column
# and not a replacement. The script runs in AntiFold's environment, by path, and
# imports nothing from this package.
#
# `per_variant` because it is one forward pass per candidate, about 0.7 s on a
# CPU (measured on 1N8Z), so the uncapped-table guard in score_variants applies:
# a scan of all 4,598 singles is about an hour. Needs AntiFold 0.4.0 or later;
# chain ids and, to reuse a download, `weights=` go in as keyword arguments.
ANTIFOLD_SEQ_SCRIPT = Path(__file__).with_name("_antifold_seq.py")
antifold_seq_process = ProcessAdapter(
    tool="antifold_seq",
    executable="python",
    command=(str(ANTIFOLD_SEQ_SCRIPT),),
    env_var="DESIGN_PIPELINE_ANTIFOLD",
    input_template=("--pdb", "{path}"),
    output_template=("--out_dir", "{path}"),
    aliases={"candidates": "mutations_csv"},
    path_kwargs=("candidates",),
    output_globs=("antifold_seq.csv",),
)
antifold_seq = ContainerScorer(
    name="antifold_seq",
    tool=antifold_seq_process,
    # A log-likelihood ratio against the wild type: already higher-is-better.
    sign=1.0,
    per_variant=True,
)
esm_if = ContainerScorer(
    name="esm_if",
    # ESM-IF is also pip-installable, but only alongside torch; the container
    # keeps torch out of this environment. The in-process variant is skipped.
    # fair-esm installs, but ESM-IF is a Python API with no command line, so
    # there is nothing to shell out to. Wiring one would mean inventing it.
    tool=None,
    sign=1.0,
)
abmpnn = ContainerScorer(
    name="abmpnn",
    tool=None,
    sign=1.0,
)
proteinmpnn_score = ContainerScorer(
    name="proteinmpnn_score",
    tool=None,
    sign=1.0,
)
consensus = ConsensusScorer()

# Per-variant scorers. A mutant structure built and relaxed per candidate.
gearbind = ContainerScorer(
    name="gearbind",
    # The one genuine environment conflict in the set, so this is the one
    # scorer that should be a container. Build it and set the image.
    tool=None,
    # A predicted ddG of binding: negative is the improvement.
    sign=-1.0,
    per_variant=True,
    value_column="ddg",
)
rosetta_ddg = ContainerScorer(
    name="rosetta_ddg",
    # Licence-gated and deliberately not automated.
    tool=None,
    sign=-1.0,
    per_variant=True,
    value_column="ddg",
)
foldx = ContainerScorer(
    name="foldx",
    # FoldX weights and binary are licensed; the image is the user's own build.
    tool=None,
    sign=-1.0,
    per_variant=True,
    value_column="ddg",
)

dummy_scorer = DummyScorer()
dummy_per_variant = DummyScorer(name="dummy_per_variant", per_variant=True)
# So the structure-free subset, and the weighting that depends on it, are
# testable with no GPU and no model weights.
dummy_sequence_scorer = DummyScorer(
    name="dummy_sequence_scorer", requires_structure=False
)

MATRIX_SCORERS: tuple[str, ...] = (
    "antifold",
    "esm_if",
    "abmpnn",
    "proteinmpnn_score",
    "consensus",
    "dummy_scorer",
    "dummy_sequence_scorer",
)
PER_VARIANT_SCORERS: tuple[str, ...] = (
    "gearbind",
    "rosetta_ddg",
    "foldx",
    "dummy_per_variant",
)
DEFAULT_SCORERS: tuple[str, ...] = ("ablang2", "consensus")
"""What `score_variants` runs when nothing is named.

Two columns that disagree by construction - sequence-based and
frequency-based - and both install into one environment:

* `ablang2`    a paired-antibody language model: the only one that knows what
               an antibody is allowed to look like at a framework position.
* `consensus`  positional repertoire frequencies: no model, free, and
               uncorrelated with the above.

`antifold` was the third column until 2026-10-31 and is not here any more
(`antifold_seq`, its per-variant counterpart, was never a default). Both
adapters and both routes stay: naming either in `models=` runs it exactly as
before, and the structure-based column is the one to put back first if a panel
shows it ranking. What took it out of the defaults is that on the ADAPT ddG
panel it does not rank. With HER2 in AntiFold's input - the condition closest
to a binding campaign - the matrix lookup over the 11 singles correlates with
ddG at Spearman +0.045 (R^2 of ddG on score 0.002), against +0.291 (R^2 0.066)
on the Fab alone, and `antifold_seq` over the combinations, held at a fixed
mutation count so that counting mutations cannot help, comes to +0.357 over the
7 doubles (+0.179 without the antigen) and 0.000 over the 4 triples either way.
Both runs are in `notebooks/03_adapt_calibration.ipynb` and
`notebooks/03_adapt_calibration_antigen.ipynb`, and CHANGELOG.md reads them.
A default that costs a second environment and a GPU has to earn more than
that.

`esm_if` is deliberately not here either. It is a Python API with no command
line, `torch-scatter` needs a compiler and `--no-build-isolation` to install,
and published benchmarks put it below AntiFold on antibody affinity ranking -
so the install cost buys a column we would trust less. The adapter and its
route stay, for anyone who wants it.
"""
