#!/usr/bin/env python3
"""Join real model, filtering and developability results by variant_id.

The scores table is the anchor. Optional filter/developability tables are
joined with an OUTER join so IDs found in any supplied file are retained.
Non-key columns from the secondary tables receive prefixes to avoid collisions.
Missing measurements remain missing; this script never fills them with zero.
"""

from __future__ import annotations

import argparse
from pathlib import Path
from typing import Any

import pandas as pd


def reject_dummy_path(path: Path, label: str) -> None:
    if any("dummy" in part.lower() or "test" in part.lower() for part in path.parts):
        raise ValueError(
            f"Refusing {label} path containing 'dummy' or 'test': {path}"
        )


def read_table(path: Path, label: str) -> pd.DataFrame:
    if not path.is_file():
        raise FileNotFoundError(f"{label} file not found: {path}")
    reject_dummy_path(path, label)
    if path.suffix.lower() == ".csv":
        frame = pd.read_csv(path)
    elif path.suffix.lower() == ".parquet":
        frame = pd.read_parquet(path)
    else:
        raise ValueError(f"{label} input must be .csv or .parquet")

    dummy_columns = [c for c in frame.columns if "dummy" in str(c).lower()]
    if dummy_columns:
        raise ValueError(
            f"{label} file contains dummy/test columns: {dummy_columns}"
        )
    if "variant_id" not in frame.columns:
        raise ValueError(f"{label} file must contain a variant_id column")
    if frame["variant_id"].isna().any():
        raise ValueError(f"{label} file contains empty variant_id values")
    frame["variant_id"] = frame["variant_id"].astype(str)
    if frame["variant_id"].duplicated().any():
        duplicates = frame.loc[
            frame["variant_id"].duplicated(keep=False), "variant_id"
        ].head(10).tolist()
        raise ValueError(
            f"{label} file must have one row per variant_id. Examples of duplicates: "
            + ", ".join(duplicates)
        )
    if frame["variant_id"].str.contains("dummy|test", case=False, regex=True).any():
        raise ValueError(f"{label} file appears to contain dummy/test variant IDs")
    return frame


def prefix_non_key(frame: pd.DataFrame, prefix: str) -> pd.DataFrame:
    return frame.rename(
        columns={col: f"{prefix}__{col}" for col in frame.columns if col != "variant_id"}
    )


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--scores", required=True, type=Path,
                        help="One-row-per-variant output from inspect_scores.py")
    parser.add_argument("--filter-decisions", type=Path, default=None,
                        help="Optional filter_decisions.csv from filter_variants.py")
    parser.add_argument("--developability", type=Path, default=None,
                        help="Optional Iason CSV/Parquet with stability/expression metrics")
    parser.add_argument("--output-dir", type=Path,
                        default=Path("runs/her2_v1/reports/merged"))
    args = parser.parse_args()

    reject_dummy_path(args.output_dir, "output")
    scores = read_table(args.scores, "scores")
    if not any(
        str(c).startswith("score_") and not str(c).endswith("_raw")
        for c in scores.columns
    ):
        raise ValueError("Scores input has no normalized score_ columns")

    merged = scores.copy()
    merged["present_in_scores"] = True
    summary_rows: list[dict[str, Any]] = [{
        "source": "scores",
        "file": str(args.scores),
        "provided": True,
        "rows": len(scores),
        "unique_variant_ids": scores["variant_id"].nunique(),
    }]

    if args.filter_decisions is not None:
        filt = read_table(args.filter_decisions, "filter decisions")
        filt = prefix_non_key(filt, "filter")
        filt["present_in_filter"] = True
        merged = merged.merge(filt, on="variant_id", how="outer", validate="one_to_one")
        summary_rows.append({
            "source": "filter_decisions",
            "file": str(args.filter_decisions),
            "provided": True,
            "rows": len(filt),
            "unique_variant_ids": filt["variant_id"].nunique(),
        })
    else:
        merged["present_in_filter"] = False
        summary_rows.append({
            "source": "filter_decisions", "file": "", "provided": False,
            "rows": 0, "unique_variant_ids": 0,
        })

    if args.developability is not None:
        dev = read_table(args.developability, "developability")
        dev = prefix_non_key(dev, "developability")
        dev["present_in_developability"] = True
        merged = merged.merge(dev, on="variant_id", how="outer", validate="one_to_one")
        summary_rows.append({
            "source": "developability",
            "file": str(args.developability),
            "provided": True,
            "rows": len(dev),
            "unique_variant_ids": dev["variant_id"].nunique(),
        })
    else:
        merged["present_in_developability"] = False
        summary_rows.append({
            "source": "developability", "file": "", "provided": False,
            "rows": 0, "unique_variant_ids": 0,
        })

    for flag in ["present_in_scores", "present_in_filter", "present_in_developability"]:
        if flag not in merged.columns:
            merged[flag] = False
        merged[flag] = merged[flag].fillna(False).astype(bool)

    # Do not silently drop unmatched IDs. Keep them in the main table and also
    # write a focused file so the team can resolve input mismatches.
    sources_provided = ["present_in_scores"]
    if args.filter_decisions is not None:
        sources_provided.append("present_in_filter")
    if args.developability is not None:
        sources_provided.append("present_in_developability")
    unmatched = merged.loc[~merged[sources_provided].all(axis=1)].copy()

    args.output_dir.mkdir(parents=True, exist_ok=True)
    merged_path = args.output_dir / "all_candidates_merged.csv"
    unmatched_path = args.output_dir / "unmatched_variant_ids.csv"
    summary_path = args.output_dir / "merge_summary.csv"

    merged.to_csv(merged_path, index=False)
    unmatched[["variant_id"] + sources_provided].to_csv(unmatched_path, index=False)
    summary = pd.DataFrame(summary_rows)
    summary = pd.concat([
        summary,
        pd.DataFrame([
            {"source": "merged_union", "file": str(merged_path), "provided": True,
             "rows": len(merged), "unique_variant_ids": merged["variant_id"].nunique()},
            {"source": "unmatched_union_ids", "file": str(unmatched_path), "provided": True,
             "rows": len(unmatched), "unique_variant_ids": unmatched["variant_id"].nunique()},
        ])
    ], ignore_index=True)
    summary.to_csv(summary_path, index=False)

    print(f"Unique candidate IDs in merged union: {merged['variant_id'].nunique()}")
    print(f"IDs needing source reconciliation: {len(unmatched)}")
    print("Reports saved:")
    print(f"  {merged_path}")
    print(f"  {unmatched_path}")
    print(f"  {summary_path}")
    print("Missing values were preserved; no score or developability value was imputed.")


if __name__ == "__main__":
    main()
