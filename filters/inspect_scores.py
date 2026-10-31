#!/usr/bin/env python3
"""Inspect and compare real PROTECher model scores.

Reads a combined score table produced by PROTECher's score_variants() or an
externally prepared table that follows the same schema. It writes one row per
variant_id, model rankings, site-specific rankings for single substitutions,
Spearman correlations, and a per-model coverage summary.

This script intentionally refuses inputs/outputs marked as dummy or test data.
It does not run models, predict affinity, or filter candidates.
"""

from __future__ import annotations

import argparse
from pathlib import Path
from typing import Any

import pandas as pd


def reject_dummy_path(path: Path, label: str) -> None:
    """Prevent accidental use or creation of dummy/test reports."""
    parts = [part.lower() for part in path.parts]
    if any("dummy" in part or "test" in part for part in parts):
        raise ValueError(
            f"Refusing {label} path containing 'dummy' or 'test': {path}. "
            "Use an actual non-dummy model output."
        )


def read_table(path: Path) -> pd.DataFrame:
    if not path.is_file():
        raise FileNotFoundError(f"Input file not found: {path}")
    reject_dummy_path(path, "input")
    suffix = path.suffix.lower()
    if suffix == ".parquet":
        frame = pd.read_parquet(path)
    elif suffix == ".csv":
        frame = pd.read_csv(path)
    else:
        raise ValueError("Input must be .parquet or .csv")
    return frame


def reject_dummy_content(data: pd.DataFrame) -> None:
    """Catch obvious dummy scorer outputs even if they were renamed."""
    dummy_columns = [
        col for col in data.columns
        if "dummy" in str(col).lower() or "test_scorer" in str(col).lower()
    ]
    if dummy_columns:
        raise ValueError(
            "Input contains dummy/test model columns: "
            + ", ".join(dummy_columns)
            + ". This script is configured for real results only."
        )


def make_one_row_per_variant(
    data: pd.DataFrame,
) -> tuple[pd.DataFrame, list[str]]:
    """Collapse mutation-level rows to one row per variant_id."""
    required = {"variant_id", "chain", "position", "wt", "mut"}
    missing = sorted(required - set(data.columns))
    if missing:
        raise ValueError(f"Input is missing required columns: {missing}")
    if data["variant_id"].isna().any():
        raise ValueError("Input has empty variant_id values")

    score_columns = [
        col for col in data.columns
        if str(col).startswith("score_") and not str(col).endswith("_raw")
    ]
    if not score_columns:
        raise ValueError(
            "No normalized score_ columns found. The input must contain real "
            "model scores produced by score_variants() or an equivalent export."
        )

    rows: list[dict[str, Any]] = []
    for variant_id, group in data.groupby("variant_id", sort=False, dropna=False):
        row: dict[str, Any] = {"variant_id": variant_id}
        mutations = group[["chain", "position", "wt", "mut"]].drop_duplicates()
        mutations = mutations.sort_values(["chain", "position", "wt", "mut"])
        row["mutation_count"] = len(mutations)
        row["mutations"] = ";".join(
            f"{m.chain}:{m.wt}{int(m.position)}{m.mut}"
            for m in mutations.itertuples(index=False)
        )

        if len(mutations) == 1:
            mutation = mutations.iloc[0]
            row["chain"] = str(mutation["chain"])
            row["position"] = int(mutation["position"])
            row["wt"] = str(mutation["wt"])
            row["mut"] = str(mutation["mut"])

        # Keep metadata that has one consistent value per candidate. The
        # residue-level columns above are reconstructed from the mutation list.
        excluded = {"variant_id", "chain", "position", "wt", "mut", *score_columns}
        for col in data.columns:
            if col in excluded or col.endswith("_raw"):
                continue
            values = group[col].dropna().astype(str).unique()
            if len(values) == 1:
                row[col] = group[col].dropna().iloc[0]

        # A model's score must be the same on every mutation row of one design.
        for col in score_columns:
            values = pd.to_numeric(group[col], errors="coerce").dropna().unique()
            if len(values) > 1:
                raise ValueError(
                    f"{variant_id}: conflicting values in {col}. Check the "
                    "score_variants() join and source export."
                )
            row[col] = float(values[0]) if len(values) == 1 else float("nan")

        rows.append(row)

    candidates = pd.DataFrame(rows)
    reject_dummy_content(candidates)
    return candidates, score_columns


def build_reports(
    data: pd.DataFrame,
) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    candidates, score_columns = make_one_row_per_variant(data)

    for col in score_columns:
        model = col.removeprefix("score_")
        candidates[f"rank_{model}"] = candidates[col].rank(
            ascending=False, method="min"
        )

    # For single substitutions, rank alternatives only against mutations at
    # the same chain and residue position. Multi-mutants have no site rank.
    single_mask = candidates["mutation_count"] == 1
    for col in score_columns:
        model = col.removeprefix("score_")
        rank_col = f"rank_at_site_{model}"
        candidates[rank_col] = float("nan")
        singles = candidates.loc[single_mask]
        if not singles.empty:
            candidates.loc[singles.index, rank_col] = (
                singles.groupby(["chain", "position"])[col]
                .rank(ascending=False, method="min")
            )

    candidates["models_scored"] = candidates[score_columns].notna().sum(axis=1)
    correlations = candidates[score_columns].corr(
        method="spearman", min_periods=3
    )

    summary = pd.DataFrame(
        [
            {
                "model": col.removeprefix("score_"),
                "score_column": col,
                "n_scored": int(candidates[col].notna().sum()),
                "n_missing": int(candidates[col].isna().sum()),
                "mean_score": candidates[col].mean(),
                "median_score": candidates[col].median(),
                "minimum_score": candidates[col].min(),
                "maximum_score": candidates[col].max(),
            }
            for col in score_columns
        ]
    )
    return candidates, correlations, summary


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input", required=True, type=Path)
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=Path("runs/her2_v1/reports/score_comparison"),
    )
    args = parser.parse_args()

    reject_dummy_path(args.output_dir, "output")
    raw = read_table(args.input)
    reject_dummy_content(raw)
    candidates, correlations, summary = build_reports(raw)

    args.output_dir.mkdir(parents=True, exist_ok=True)
    candidates_path = args.output_dir / "all_candidates_compared.csv"
    correlation_path = args.output_dir / "model_spearman_correlation.csv"
    summary_path = args.output_dir / "model_summary.csv"

    candidates.to_csv(candidates_path, index=False)
    correlations.to_csv(correlation_path)
    summary.to_csv(summary_path, index=False)

    print(f"Input rows: {len(raw)}")
    print(f"Unique variants: {len(candidates)}")
    print(f"Real score columns found: {len(summary)}")
    for row in summary.itertuples(index=False):
        print(f"  {row.score_column}: {row.n_scored} scores present; {row.n_missing} missing")
    print("\nReports saved:")
    print(f"  {candidates_path}")
    print(f"  {correlation_path}")
    print(f"  {summary_path}")
    print("\nReporting only: no candidate has been selected or filtered.")


if __name__ == "__main__":
    main()
