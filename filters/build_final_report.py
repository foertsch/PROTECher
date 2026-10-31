#!/usr/bin/env python3
"""Build a transparent candidate report for team review.

Input: all_candidates_merged.csv from merge_developability.py.
Outputs: all_candidates.csv, final_shortlist.csv, selection_summary.csv.

No composite score is invented, and the script does not autonomously select the
final 15 designs. PASS/PROVISIONAL_PASS candidates are placed in a review queue;
the team makes the final design selection using the agreed Track A design arms
and the complete evidence.
"""

from __future__ import annotations

import argparse
from pathlib import Path

import pandas as pd


def reject_dummy_path(path: Path, label: str) -> None:
    if any("dummy" in part.lower() or "test" in part.lower() for part in path.parts):
        raise ValueError(
            f"Refusing {label} path containing 'dummy' or 'test': {path}"
        )


def read_table(path: Path) -> pd.DataFrame:
    if not path.is_file():
        raise FileNotFoundError(f"Input file not found: {path}")
    reject_dummy_path(path, "input")
    if path.suffix.lower() == ".csv":
        data = pd.read_csv(path)
    elif path.suffix.lower() == ".parquet":
        data = pd.read_parquet(path)
    else:
        raise ValueError("Input must be .csv or .parquet")
    dummy_cols = [c for c in data.columns if "dummy" in str(c).lower()]
    if dummy_cols:
        raise ValueError(
            "Dummy/test columns are not accepted: " + ", ".join(dummy_cols)
        )
    if "variant_id" not in data.columns:
        raise ValueError("Merged input must contain variant_id")
    if data["variant_id"].isna().any() or data["variant_id"].duplicated().any():
        raise ValueError("variant_id must be present and unique")
    return data


def choose_filter_status_column(data: pd.DataFrame, configured: str | None) -> str | None:
    if configured:
        return configured if configured in data.columns else None
    for candidate in ("filter__filter_status", "filter_status"):
        if candidate in data.columns:
            return candidate
    return None


def build_summary(data: pd.DataFrame, shortlist: pd.DataFrame, status_col: str | None) -> pd.DataFrame:
    rows = [
        {"metric": "total_unique_candidates", "value": int(data["variant_id"].nunique())},
        {"metric": "candidates_in_review_shortlist", "value": int(len(shortlist))},
        {"metric": "final_15_selected_by_script", "value": 0},
        {"metric": "note", "value": "Final 15 are selected by the team; this script only creates a review shortlist."},
    ]

    if status_col is None:
        rows.append({
            "metric": "filter_status_column",
            "value": "NOT FOUND - shortlist is empty until filter decisions are merged",
        })
    else:
        counts = data[status_col].fillna("MISSING").astype(str).value_counts()
        for status, count in counts.items():
            rows.append({"metric": f"filter_status::{status}", "value": int(count)})
        rows.append({"metric": "filter_status_column", "value": status_col})

    if "present_in_developability" in data.columns:
        complete = data["present_in_developability"].fillna(False).astype(bool)
        rows.append({"metric": "candidates_with_developability_record", "value": int(complete.sum())})
        rows.append({"metric": "candidates_missing_developability_record", "value": int((~complete).sum())})
    else:
        dev_cols = [c for c in data.columns if str(c).startswith("developability__")]
        if dev_cols:
            complete = data[dev_cols].notna().any(axis=1)
            rows.append({"metric": "candidates_with_developability_values", "value": int(complete.sum())})
            rows.append({"metric": "candidates_missing_developability_values", "value": int((~complete).sum())})
        else:
            rows.append({"metric": "developability_data", "value": "NOT PROVIDED"})

    for flag in ("present_in_scores", "present_in_filter", "present_in_developability"):
        if flag in data.columns:
            rows.append({
                "metric": f"records_{flag}",
                "value": int(data[flag].fillna(False).astype(bool).sum()),
            })

    return pd.DataFrame(rows)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input", required=True, type=Path,
                        help="all_candidates_merged.csv from merge_developability.py")
    parser.add_argument("--output-dir", type=Path,
                        default=Path("runs/her2_v1/reports/final"))
    parser.add_argument("--filter-status-column", default=None,
                        help="Optional explicit filter status column name")
    parser.add_argument("--exclude-provisional", action="store_true",
                        help="Exclude PROVISIONAL_PASS candidates from review shortlist")
    args = parser.parse_args()

    reject_dummy_path(args.output_dir, "output")
    data = read_table(args.input)
    status_col = choose_filter_status_column(data, args.filter_status_column)

    all_path = args.output_dir / "all_candidates.csv"
    shortlist_path = args.output_dir / "final_shortlist.csv"
    summary_path = args.output_dir / "selection_summary.csv"
    args.output_dir.mkdir(parents=True, exist_ok=True)

    # Preserve every candidate in the full report, including missing data.
    data.to_csv(all_path, index=False)

    if status_col is None:
        shortlist = data.iloc[0:0].copy()
        shortlist["review_status"] = pd.Series(dtype="object")
    else:
        statuses = data[status_col].fillna("").astype(str).str.upper()
        keep = statuses.eq("PASS")
        if not args.exclude_provisional:
            keep |= statuses.eq("PROVISIONAL_PASS")
        shortlist = data.loc[keep].copy()
        shortlist["review_status"] = statuses.loc[keep].map({
            "PASS": "PASS - eligible for team review",
            "PROVISIONAL_PASS": "PROVISIONAL - calibration remains incomplete",
        })

    shortlist["not_finally_selected"] = True
    shortlist["selection_note"] = (
        "Review queue only. The team must confirm all Track A requirements, "
        "developability evidence, structure checks, mutation-count/spread rules, "
        "and choose the final 15 designs."
    )
    shortlist.to_csv(shortlist_path, index=False)

    summary = build_summary(data, shortlist, status_col)
    summary.to_csv(summary_path, index=False)

    print(f"All candidates retained: {len(data)}")
    print(f"Candidates in team-review shortlist: {len(shortlist)}")
    if status_col is None:
        print("WARNING: no filter status column found; shortlist is empty.")
    else:
        print(f"Filter status column: {status_col}")
    print("Reports saved:")
    print(f"  {all_path}")
    print(f"  {shortlist_path}")
    print(f"  {summary_path}")
    print("No composite score was calculated; final 15 remain a team decision.")


if __name__ == "__main__":
    main()
