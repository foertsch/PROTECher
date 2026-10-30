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

from ..containers import ContainerAdapter, make_id
from ..designs import Design
from ..errors import MissingDependencyError, RunFailedError
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
    """A scoring container plus the three facts needed to read it back.

    The container itself is an ordinary `ContainerAdapter` - same runtime
    detection, same mounts, same GPU flags - configured with the candidates CSV
    as a path kwarg so it is bind-mounted and rewritten like any other input
    file. This class only normalizes the output.
    """

    name: str
    container: ContainerAdapter
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
        produced = self.container.run(
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


def _matrix_container(
    image: str,
    *,
    mutations_flag: str = "mutations_csv",
    gpu: bool = True,
    fixed_args: tuple[str, ...] = (),
    command: tuple[str, ...] = (),
) -> ContainerAdapter:
    """A scoring image: structure in, candidates CSV in, CSV out.

    `per_input=True` with exactly one input - the structure - so the container
    layer mounts the structure's directory, the candidates CSV's directory and
    the run directory, and nothing else.
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


# Matrix scorers. One forward pass over the complex, then a lookup per single
# mutation; a multi-mutation candidate is re-run conditioned on the full mutant
# sequence, which is the tool's business and the reason the candidates CSV
# carries every mutation of every id.
antifold = ContainerScorer(
    name="antifold",
    container=_matrix_container("biolib/antifold"),
    # AntiFold reports a log-likelihood ratio against the wild type: already
    # higher-is-better.
    sign=1.0,
)
esm_if = ContainerScorer(
    name="esm_if",
    # ESM-IF is also pip-installable, but only alongside torch; the container
    # keeps torch out of this environment. The in-process variant is skipped.
    container=_matrix_container("facebookresearch/esm-if1"),
    sign=1.0,
)
abmpnn = ContainerScorer(
    name="abmpnn",
    container=_matrix_container("oxpig/abmpnn"),
    sign=1.0,
)
proteinmpnn_score = ContainerScorer(
    name="proteinmpnn_score",
    container=_matrix_container("dauparas/proteinmpnn"),
    sign=1.0,
)
consensus = ConsensusScorer()

# Per-variant scorers. A mutant structure built and relaxed per candidate.
gearbind = ContainerScorer(
    name="gearbind",
    container=_matrix_container("deepgraphlearning/gearbind"),
    # A predicted ddG of binding: negative is the improvement.
    sign=-1.0,
    per_variant=True,
    value_column="ddg",
)
rosetta_ddg = ContainerScorer(
    name="rosetta_ddg",
    container=_matrix_container(
        "rosettacommons/rosetta",
        gpu=False,
        command=("rosetta_scripts.linuxgccrelease",),
    ),
    sign=-1.0,
    per_variant=True,
    value_column="ddg",
)
foldx = ContainerScorer(
    name="foldx",
    # FoldX weights and binary are licensed; the image is the user's own build.
    container=_matrix_container("local/foldx", gpu=False),
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
DEFAULT_SCORERS: tuple[str, ...] = ("antifold", "esm_if", "consensus")
"""What `score_variants` runs when nothing is named: the cheap set, AntiFold
first because it is the antibody-specific one."""
