"""The CSV files a per-candidate adapter reads and writes.

Scorers and gate checks both take the same two inputs - a structure and a list
of candidates - and both answer one value per candidate. That is why neither
needs its own Protocol: they satisfy the existing `run(inputs, output_dir,
**kwargs) -> list[Design]`, and what differs is only the small table they
write into their run directory.

    candidates.csv   variant_id,chain,position,wt,mut   (one row per mutation)
    scores.csv       variant_id,score[,raw]             (one row per candidate)
    checks.csv       variant_id,passed,reason           (one row per candidate)

CSV rather than parquet on purpose: a container has to be able to read the
input with nothing installed, and `csv` is in the standard library.
"""

from __future__ import annotations

import csv
import math
from pathlib import Path
from typing import Iterable, Mapping, Sequence

CANDIDATES = "candidates.csv"
SCORES = "scores.csv"
CHECKS = "checks.csv"

Mutation = tuple[str, int, str, str]
"""(chain, position, wt, mut), as in design_pipeline.variants."""


def write_candidates(path: Path, rows: Iterable[Mapping[str, object]]) -> Path:
    """The candidate set handed to one adapter, in variant-table column order."""
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.writer(handle)
        writer.writerow(["variant_id", "chain", "position", "wt", "mut"])
        for row in rows:
            writer.writerow(
                [
                    row["variant_id"],
                    row["chain"],
                    int(row["position"]),  # type: ignore[arg-type]
                    row["wt"],
                    row["mut"],
                ]
            )
    return path


def read_candidates(path: Path) -> dict[str, tuple[Mutation, ...]]:
    """variant_id -> every mutation it names, in file order.

    Grouping here rather than in each adapter is what keeps a pair an ordinary
    candidate: an adapter sees one id with two mutations and scores the mutant,
    with no opportunity to score the two halves and add them up.
    """
    grouped: dict[str, list[Mutation]] = {}
    with Path(path).open(newline="", encoding="utf-8") as handle:
        for row in csv.DictReader(handle):
            grouped.setdefault(row["variant_id"], []).append(
                (row["chain"], int(row["position"]), row["wt"], row["mut"])
            )
    return {key: tuple(value) for key, value in grouped.items()}


def write_scores(
    path: Path, rows: Iterable[tuple[str, float | None, float | None]]
) -> Path:
    """variant_id, normalized score, raw score.

    `score` is normalized at the adapter boundary so higher always means better
    binding, whatever the tool's own convention. `raw` is what the tool said,
    kept so a flipped convention is recoverable without re-running anything. An
    empty cell is "no score for this candidate", which is a legitimate answer -
    a positional model cannot score a double mutant - and reads back as NaN.
    """
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.writer(handle)
        writer.writerow(["variant_id", "score", "raw"])
        for variant_id, score, raw in rows:
            writer.writerow([variant_id, _cell(score), _cell(raw)])
    return path


def write_checks(path: Path, rows: Iterable[tuple[str, bool, str]]) -> Path:
    """variant_id, passed, reason. The reason is for the failures."""
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.writer(handle)
        writer.writerow(["variant_id", "passed", "reason"])
        for variant_id, passed, reason in rows:
            writer.writerow([variant_id, "true" if passed else "false", reason])
    return path


def _cell(value: float | None) -> str:
    if value is None or (isinstance(value, float) and math.isnan(value)):
        return ""
    return repr(float(value))


def read_tool_table(
    paths: Sequence[Path], id_column: str, value_column: str
) -> dict[str, float]:
    """A tool's own CSV, reduced to id -> value.

    The column names are the tool's, passed in by its adapter and overridable
    per call, because they move between versions. A row whose value does not
    parse as a number is skipped rather than guessed at.
    """
    found: dict[str, float] = {}
    for path in paths:
        if path.suffix.lower() != ".csv":
            continue
        with path.open(newline="", encoding="utf-8", errors="replace") as handle:
            reader = csv.DictReader(handle)
            if reader.fieldnames is None:
                continue
            if id_column not in reader.fieldnames or value_column not in reader.fieldnames:
                continue
            for row in reader:
                try:
                    found[row[id_column]] = float(row[value_column])
                except (TypeError, ValueError):
                    continue
    if not found:
        names = ", ".join(str(p.name) for p in paths) or "nothing"
        raise ValueError(
            f"no rows with columns {id_column!r} and {value_column!r} in {names}. "
            f"Pass id_column= and value_column= if this version of the tool "
            f"names them differently."
        )
    return found
