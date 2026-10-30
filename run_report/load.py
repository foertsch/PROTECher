"""Read the design records that design_pipeline writes (the DONE manifests) into one table."""

from __future__ import annotations

import json
from pathlib import Path
from typing import Iterable, Iterator

import pandas as pd

_SCALARS = (str, int, float, bool)


def find_manifests(roots: Iterable[Path | str]) -> Iterator[Path]:
    for root in map(Path, roots):
        if root.is_file() and root.name == "DONE":
            yield root
        elif root.is_dir():
            yield from sorted(root.rglob("DONE"))


def load_designs(roots: Iterable[Path | str]) -> pd.DataFrame:
    """One row per design across every manifest under `roots`.

    A `designs_meta.csv` in a root directory adds columns the models cannot know, such as the
    mutation set and the design group. It needs a `design` column holding the id of the sequence
    that was folded (the `parent_id` of a prediction), and every other column becomes
    `meta_<column>`.

    Scores become columns under their own names, scalar meta values become
    `meta_<key>` columns, and `step` is the name of the directory the step wrote
    into (`predictions`, `sequences`, ...). A design that appears in several
    manifests (a cache hit copied forward) is kept once, last one wins.
    """
    rows: dict[str, dict] = {}
    for marker in find_manifests(roots):
        try:
            records = json.loads(marker.read_text())
        except (OSError, json.JSONDecodeError):
            continue
        if not isinstance(records, list):
            continue
        for record in records:
            if not isinstance(record, dict) or "id" not in record:
                continue
            row: dict = {
                "id": record["id"],
                "parent_id": record.get("parent_id"),
                "producer": record.get("producer"),
                "step": marker.parent.parent.name,
            }
            for key, value in (record.get("scores") or {}).items():
                if isinstance(value, (int, float)) and not isinstance(value, bool):
                    row[key] = float(value)
            for key, value in (record.get("meta") or {}).items():
                if isinstance(value, _SCALARS):
                    row[f"meta_{key}"] = value
            rows[record["id"]] = row
    return _attach_sidecars(pd.DataFrame(list(rows.values())), roots)


def _attach_sidecars(df: pd.DataFrame, roots: Iterable[Path | str]) -> pd.DataFrame:
    if df.empty or "parent_id" not in df:
        return df
    for root in map(Path, roots):
        sidecar = root / "designs_meta.csv"
        if not sidecar.is_file():
            continue
        side = pd.read_csv(sidecar)
        if "design" not in side:
            continue
        side = side.drop_duplicates("design").set_index("design")
        for column in side.columns:
            name = f"meta_{column}"
            values = df["parent_id"].map(side[column])
            df[name] = values if name not in df else df[name].where(df[name].notna(), values)
    return df
