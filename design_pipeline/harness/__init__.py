"""The validation harness: the ordinary pipeline in a loop, scored.

    from design_pipeline.harness import run_systems, report

    results = run_systems(
        systems=["1n8z_her2", "4fqi_flu"],
        models=["antifold", "esm_if", "gearbind", "consensus"],
        output_dir="runs/harness",
    )
    print(report(results))
    magnitude_plot(results, "runs/harness/magnitude.png")

This is the single most important thing in the package. A campaign gets no
feedback while it is running, so the only way to know which scorer to trust is
to have measured it beforehand on systems where the answers are known.

It is `get_structure`, `enumerate_variants` and `score_variants` in a for
loop, with a correlation bolted on the end. It is not a second code path, and
if it ever needs to become one, the abstraction underneath it is wrong.

`report` breaks everything out per scorer, never as an ensemble. An
ensemble-level number is precisely what hides the thing you came here to learn:
which single column to rank on.

Every statistic it prints comes from `design_pipeline.metrics`, which the
validation notebook also calls. The harness deliberately contains none of its
own: a Spearman here and a Spearman in a plot have to be the same Spearman.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from pathlib import Path
from typing import TYPE_CHECKING, Any, Mapping, Sequence

import numpy as np

from ..adapters._scorers import DEFAULT_SCORERS
from .. import metrics as _metrics
from ..api import get_structure, score_variants
from ..variants import enumerate_variants
from ..metrics import (
    magnitude_agreement,
    nanmean,
    precision_at_k,
    score_columns,
    scorer_name,
    spearman,
    top_k_recall,
)
from .systems import (
    FIXTURE_SYSTEMS,
    GROUND_TRUTH_ENV,
    SYSTEMS,
    System,
    get_system,
)

if TYPE_CHECKING:
    import pandas as pd

__all__ = [
    "System",
    "SYSTEMS",
    "FIXTURE_SYSTEMS",
    "GROUND_TRUTH_ENV",
    "HarnessResults",
    "run_systems",
    "report",
    "magnitude_plot",
    "spearman",
    "precision_at_k",
    "top_k_recall",
    "magnitude_agreement",
]

MEASURED = "measured_affinity"
TRUTH_COLUMNS = ("variant_id", MEASURED)


@dataclass(frozen=True)
class HarnessResults:
    """One row per (system, candidate): the score columns plus what was measured.

    `oriented` is the measured column flipped where a system reports a KD or a
    ddG, so that a positive correlation always means the scorer agreed with
    better binding, whatever the assay wrote down.
    """

    path: Path
    systems: tuple[str, ...]
    models: tuple[str, ...]
    rows: int
    config: Mapping[str, Any] = field(default_factory=dict)

    @property
    def frame(self) -> "pd.DataFrame":
        import pandas as pd

        return pd.read_parquet(self.path)

    def __len__(self) -> int:
        return self.rows


def run_systems(
    systems: Sequence[str] | None = None,
    models: Sequence[str] = DEFAULT_SCORERS,
    output_dir: Path | str = Path("runs") / "harness",
    *,
    max_variants: int | None = None,
    cache: bool = True,
    **kwargs: Any,
) -> HarnessResults:
    """Score every measured variant of every system with every model.

    `systems` defaults to the shipped fixtures, so an unconfigured call runs
    offline and proves the loop works rather than failing on a download.

    `enumerate_variants(only=...)` restricts each table to exactly the
    candidates that were measured - which also validates them: every
    wild-type residue is checked against the structure, so a set numbered in
    IMGT rather than in chain positions is an error here instead of a
    correlation computed against the wrong residues.
    """
    import pandas as pd

    names = list(systems) if systems else list(FIXTURE_SYSTEMS)
    root = Path(output_dir)
    frames: list[pd.DataFrame] = []

    for name in names:
        system = get_system(name)
        truth = _read_truth(system.truth_path())
        complex_ = get_structure(
            root / name,
            model_name=system.model_name,
            pdb_id=system.pdb_id,
            source=system.source_path(),
            chains=list(system.chains) or None,
            cache=cache,
        )
        variants = enumerate_variants(
            complex_,
            chain=system.chain,
            only=truth["variant_id"],
        )
        scores = score_variants(
            complex_,
            variants,
            root / name / "scores",
            models=models,
            max_variants=max_variants,
            cache=cache,
            **kwargs,
        )
        # One row per candidate: the table has one row per mutation, and a
        # double mutant would otherwise be counted twice in the correlation.
        frame = scores.frame.drop_duplicates("variant_id")
        frame = frame.merge(truth, on="variant_id", how="inner")
        frame.insert(0, "system", name)
        frame["oriented"] = (
            frame[MEASURED] if system.higher_is_better else -frame[MEASURED]
        )
        # One flat metrics.json per system, in the same schema every other step
        # writes, so the notebook reads harness output the same way it reads a
        # campaign's.
        _metrics.save(
            root / name,
            {
                "system": name,
                "measured": len(frame),
                **{
                    f"spearman_{scorer_name(column)}": spearman(
                        frame[column], frame["oriented"]
                    )
                    for column in score_columns(frame)
                },
            },
        )
        frames.append(frame)

    combined = pd.concat(frames, ignore_index=True)
    root.mkdir(parents=True, exist_ok=True)
    path = root / "results.parquet"
    combined.to_parquet(path, index=False)
    return HarnessResults(
        path=path,
        systems=tuple(names),
        models=tuple(models),
        rows=len(combined),
        config={"max_variants": max_variants, "output_dir": str(root)},
    )


def _read_truth(path: Path) -> "pd.DataFrame":
    import pandas as pd

    truth = pd.read_csv(path)
    missing = [c for c in TRUTH_COLUMNS if c not in truth.columns]
    if missing:
        raise ValueError(
            f"{path} is not a ground-truth table: missing {missing}; expected "
            f"{','.join(TRUTH_COLUMNS)}"
        )
    return truth[list(TRUTH_COLUMNS)].drop_duplicates("variant_id")


def report(results: HarnessResults, k: int = 10) -> str:
    """Spearman and precision@k, per scorer and per system, then magnitudes.

    Per scorer, because choosing which column to trust is the entire point. The
    `(mean)` row averages each scorer's per-system correlations rather than
    pooling the rows: two assays do not share a scale, so a pooled correlation
    would mostly measure the gap between the systems.

    The second table answers a different question from the first, and one that
    only the harness can answer: whether score *magnitude* tracks effect size
    or whether only the order does. A column with a good Spearman and a poor
    r^2 ranks candidates correctly and scales badly, so its numbers cannot
    carry a threshold - candidates have to be taken by rank. `magnitude_plot`
    draws the same thing, which is the honest way to read it, because a
    straight-line fit through a scatter that is not straight still returns a
    slope.

    A blank cell is an honest NaN - too few measured variants, or a column the
    scorer declined to answer, as a positional model does for a double mutant.
    """
    frame = results.frame
    columns = score_columns(frame)
    if not columns:
        return "no score columns in these results"

    lines = [
        f"{len(frame)} measured variants, "
        f"{frame['system'].nunique()} system(s), k={k}",
        "",
        f"{'scorer':<20}{'system':<16}{'n':>6}{'spearman':>10}"
        f"{'p@' + str(k):>8}{'r@' + str(2 * k):>8}",
        "-" * 68,
    ]
    for column in columns:
        scorer = scorer_name(column)
        correlations: list[float] = []
        precisions: list[float] = []
        recalls: list[float] = []
        for system, group in frame.groupby("system", sort=True):
            rho = spearman(group[column], group["oriented"])
            precision = precision_at_k(group[column], group["oriented"], k=k)
            # A shortlist twice the length of the set of things that worked:
            # what you would actually order, against what actually bound.
            recall = top_k_recall(group[column], group["oriented"], k=2 * k, n_true=k)
            answered = int(group[column].notna().sum())
            lines.append(
                f"{scorer:<20}{system:<16}{answered:>6}"
                f"{_cell(rho):>10}{_cell(precision):>8}{_cell(recall):>8}"
            )
            correlations.append(rho)
            precisions.append(precision)
            recalls.append(recall)
        if len(correlations) > 1:
            lines.append(
                f"{scorer:<20}{'(mean)':<16}{'':>6}"
                f"{_cell(nanmean(correlations)):>10}"
                f"{_cell(nanmean(precisions)):>8}"
                f"{_cell(nanmean(recalls)):>8}"
            )
    return "\n".join([*lines, "", *_magnitude_lines(frame, columns)])


def _magnitude_lines(frame: "pd.DataFrame", columns: Sequence[str]) -> list[str]:
    """The magnitude table: does score size predict effect size, or only order?

    Per system as well as per scorer, for the same reason the first table is:
    the fit is in the assay's own units, and two assays pooled would be fitting
    a line through two clouds.
    """
    lines = [
        "score magnitude against measured effect size (rank-only = orders "
        "well, scales badly)",
        f"{'scorer':<20}{'system':<16}{'n':>6}{'pearson':>9}{'r2':>7}"
        f"{'slope':>10}{'rank-only':>11}",
        "-" * 79,
    ]
    for column in columns:
        for system, group in frame.groupby("system", sort=True):
            fit = magnitude_agreement(group[column], group["oriented"])
            lines.append(
                f"{scorer_name(column):<20}{system:<16}{int(fit['n']):>6}"
                f"{_cell(fit['pearson']):>9}{_cell(fit['r2']):>7}"
                f"{_cell(fit['slope']):>10}"
                f"{('yes' if fit['rank_only'] else ''):>11}"
            )
    return lines


def magnitude_plot(
    results: HarnessResults,
    path: Path | str | None = None,
    *,
    columns: Sequence[str] | None = None,
):
    """Scatter each scorer's column against the measured effect, with its fit.

    The step-1 notebook asks whether confidence predicts structural accuracy.
    This is the same question one step later and the one that decides how the
    shortlist is cut: whether score magnitude predicts effect size, not just
    rank order. If a mutation scoring twice as high does not bind twice as
    tightly, only the ordering is usable and candidates are picked by rank
    rather than by a threshold on the score.

    It is a plot and not only a number on purpose. `magnitude_agreement`
    reports a slope and an r^2 for every column, and a straight line through a
    scatter that is flat at one end and vertical at the other reports a
    perfectly respectable slope. The shape is the answer; the fit is a summary
    of it.

    One panel per score column, one marker colour per system, since the
    measured values are in each assay's own units and must not be pooled.
    Returns the matplotlib figure, and writes it to `path` when given.
    Matplotlib is the `notebooks` extra, imported here so the rest of the
    harness runs without it.
    """
    import matplotlib

    matplotlib.use("Agg")  # a harness run is a script, not a session
    import matplotlib.pyplot as plt

    frame = results.frame
    chosen = list(columns) if columns is not None else score_columns(frame)
    if not chosen:
        raise ValueError("no score columns in these results to plot")

    width = min(4, len(chosen))
    height = -(-len(chosen) // width)
    figure, axes = plt.subplots(
        height, width, figsize=(4.0 * width, 3.4 * height), squeeze=False
    )
    flat = [ax for row in axes for ax in row]
    for ax, column in zip(flat, chosen):
        fits: list[str] = []
        for system, group in frame.groupby("system", sort=True):
            points = ax.scatter(
                group[column], group["oriented"], s=14, alpha=0.7, label=system
            )
            # One fit per system, in that system's colour. Not one line through
            # everything: the measured values are in each assay's own units, so
            # a pooled line would mostly be fitting the gap between two clouds.
            fit = magnitude_agreement(group[column], group["oriented"])
            if fit["slope"] == fit["slope"]:
                x = np.linspace(group[column].min(), group[column].max(), 2)
                ax.plot(
                    x,
                    fit["slope"] * x + fit["intercept"],
                    lw=1,
                    color=points.get_facecolor()[0],
                )
            fits.append(
                f"{system} r2={_cell(fit['r2'])}"
                f"{' rank-only' if fit['rank_only'] else ''}"
            )
        ax.set_title(f"{scorer_name(column)}\n" + "  ".join(fits), fontsize=8)
        ax.set_xlabel("score")
        ax.set_ylabel("measured (oriented)")
        if frame["system"].nunique() > 1:
            ax.legend(fontsize=7)
    for ax in flat[len(chosen) :]:
        ax.set_axis_off()
    figure.tight_layout()
    if path is not None:
        path = Path(path)
        path.parent.mkdir(parents=True, exist_ok=True)
        figure.savefig(path, dpi=150)
    return figure


def _cell(value: float) -> str:
    return "" if value != value else f"{value:.3f}"
