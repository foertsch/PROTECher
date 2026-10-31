
#!/usr/bin/env python3
"""Inspect and compare PROTECher model scores.

This script reads an existing combined score table and produces:
1. One row per variant, with model scores and rankings.
2. Rankings among alternatives at the same position.
3. A Spearman correlation matrix between model scores.
4. A summary of available and missing scores.

It does NOT run protein models or apply the final Track A filter.
"""

from __future__ import annotations

import argparse
from pathlib import Path

import pandas as pd


def load_results(path: Path) -> pd.DataFrame:
    """Read a PROTECher Parquet or CSV score table."""
    if not path.exists():
        raise FileNotFoundError(f"Cannot find input file: {path}")

    suffix = path.suffix.lower()

    if suffix == ".parquet":
        return pd.read_parquet(path)

    if suffix == ".csv":
        return pd.read_csv(path)

    raise ValueError("Input must end in .parquet or .csv")


def make_one_row_per_variant(data: pd.DataFrame):
    """Convert mutation-level rows into one row per candidate variant."""

    required = {"variant_id", "chain", "position", "wt", "mut"}
    missing = sorted(required - set(data.columns))

    if missing:
        raise ValueError(f"Input is missing required columns: {missing}")

    score_columns = [
        col for col in data.columns
        if col.startswith("score_") and not col.endswith("_raw")
    ]

    if not score_columns:
        raise ValueError(
            "No normalized score_ columns found. Have the models been run?"
        )

    rows = []

    for variant_id, group in data.groupby(
        "variant_id", sort=False, dropna=False
    ):
        row = {"variant_id": variant_id}

        mutations = group[
            ["chain", "position", "wt", "mut"]
        ].drop_duplicates()

        mutations = mutations.sort_values(
            ["chain", "position", "wt", "mut"]
        )

        row["mutation_count"] = len(mutations)

        row["mutations"] = ";".join(
            f"{m.chain}:{m.wt}{int(m.position)}{m.mut}"
            for m in mutations.itertuples(index=False)
        )

        if len(mutations) == 1:
            m = mutations.iloc[0]

            row["chain"] = str(m["chain"])
            row["position"] = int(m["position"])
            row["wt"] = str(m["wt"])
            row["mut"] = str(m["mut"])

        for col in score_columns:
            values = pd.to_numeric(
                group[col], errors="coerce"
            ).dropna().unique()

            if len(values) > 1:
                raise ValueError(
                    f"{variant_id}: {col} has conflicting scores"
                )

            row[col] = (
                float(values[0]) if len(values) == 1 else float("nan")
            )

        rows.append(row)

    candidates = pd.DataFrame(rows)

    return candidates, score_columns


def build_reports(data: pd.DataFrame):
    """Calculate rankings, model agreement and score summaries."""

    candidates, score_columns = make_one_row_per_variant(data)

    # Rank all candidates for each model.
    for col in score_columns:
        model = col.removeprefix("score_")

        candidates[f"rank_{model}"] = candidates[col].rank(
            ascending=False,
            method="min"
        )

    # Rank single mutations against alternatives at the same site.
    single_mask = candidates["mutation_count"] == 1

    for col in score_columns:
        model = col.removeprefix("score_")
        rank_col = f"rank_at_site_{model}"

        candidates[rank_col] = float("nan")

        singles = candidates.loc[single_mask].copy()

        if not singles.empty:
            candidates.loc[singles.index, rank_col] = (
                singles.groupby(["chain", "position"])[col]
                .rank(ascending=False, method="min")
            )

    candidates["models_scored"] = (
        candidates[score_columns].notna().sum(axis=1)
    )

    # Spearman correlation between model scores.
    correlations = candidates[score_columns].corr(
        method="spearman",
        min_periods=3
    )

    # Summary for each model.
    summary = pd.DataFrame([
        {
            "model": col.removeprefix("score_"),
            "n_scored": int(candidates[col].notna().sum()),
            "n_missing": int(candidates[col].isna().sum()),
            "mean_score": candidates[col].mean(),
            "median_score": candidates[col].median(),
            "minimum_score": candidates[col].min(),
            "maximum_score": candidates[col].max(),
        }
        for col in score_columns
    ])

    return candidates, correlations, summary


def main():
    parser = argparse.ArgumentParser(
        description="Compare PROTECher model scores."
    )

    parser.add_argument(
        "--input",
        required=True,
        type=Path,
        help="Input score table (.parquet or .csv)"
    )

    parser.add_argument(
        "--output-dir",
        type=Path,
        default=Path("runs/score_report")
    )

    args = parser.parse_args()

    raw = load_results(args.input)

    candidates, correlations, summary = build_reports(raw)

    args.output_dir.mkdir(parents=True, exist_ok=True)

    candidates_path = (
        args.output_dir / "all_candidates_compared.csv"
    )

    correlation_path = (
        args.output_dir / "model_spearman_correlation.csv"
    )

    summary_path = (
        args.output_dir / "model_summary.csv"
    )

    candidates.to_csv(candidates_path, index=False)
    correlations.to_csv(correlation_path)
    summary.to_csv(summary_path, index=False)

    score_columns = [
        c for c in candidates.columns
        if c.startswith("score_") and not c.endswith("_raw")
    ]

    print(f"Input rows: {len(raw)}")
    print(f"Unique variants: {len(candidates)}")
    print(f"Score columns found: {len(score_columns)}")

    for col in score_columns:
        print(
            f"  {col}: "
            f"{candidates[col].notna().sum()} scores present"
        )

    print("\nReports saved:")
    print(f"  {candidates_path}")
    print(f"  {correlation_path}")
    print(f"  {summary_path}")

    print(
        "\nNote: comparison/reporting only. "
        "No candidate has been filtered."
    )


if __name__ == "__main__":
    main()
