
#!/usr/bin/env python3
"""Apply the PROTECher Track A Stage 3 single-mutation consensus screen.

Input: the one-row-per-variant report created by inspect_scores.py.
Output: candidate-level filtering decisions and a model summary.

This script does not run protein models.
It does not perform the later liability, energy, stability,
multi-mutation design or structural-veto stages.
"""

from __future__ import annotations

import argparse
import math
import tomllib
from pathlib import Path
from typing import Any

import pandas as pd


def read_table(path: Path) -> pd.DataFrame:
    """Read a CSV or Parquet report."""

    if not path.is_file():
        raise FileNotFoundError(
            f"Input file does not exist: {path}"
        )

    if path.suffix.lower() == ".csv":
        return pd.read_csv(path)

    if path.suffix.lower() == ".parquet":
        return pd.read_parquet(path)

    raise ValueError(
        "Input must be a .csv or .parquet file"
    )


def read_config(path: Path) -> dict[str, Any]:
    """Read filtering settings from TOML."""

    if not path.is_file():
        raise FileNotFoundError(
            f"Config file does not exist: {path}"
        )

    with path.open("rb") as handle:
        return tomllib.load(handle)


def model_is_active(
    model: dict[str, Any],
    chance_cutoff: float
) -> tuple[bool, str, bool]:
    """Exclude models with measured AUROC <= chance."""

    if not bool(model.get("enabled", True)):
        return False, "disabled in config", False

    raw_auc = model.get("auroc", "not_measured")

    try:
        auc = float(raw_auc)
    except (TypeError, ValueError):
        return (
            True,
            "active provisionally; AUROC not measured",
            True
        )

    if not math.isfinite(auc) or auc < 0 or auc > 1:
        return (
            True,
            "active provisionally; AUROC unavailable",
            True
        )

    if auc <= chance_cutoff:
        return (
            False,
            f"excluded: AUROC {auc:.3f} <= {chance_cutoff:.3f}",
            False
        )

    return (
        True,
        f"active; AUROC {auc:.3f}",
        False
    )


def numeric_value(value: Any) -> float | None:
    """Convert a cell to a number."""

    parsed = pd.to_numeric(
        pd.Series([value]),
        errors="coerce"
    ).iloc[0]

    if pd.isna(parsed):
        return None

    return float(parsed)


def summarize_models(
    groups: list[tuple[str, list[dict[str, Any]]]],
    data: pd.DataFrame,
    chance_cutoff: float
):
    """Select active models and create a model summary."""

    active_groups = {}
    summary_rows = []
    pending_models = []

    for group_name, models in groups:

        active = []

        for model in models:

            use_model, reason, pending = model_is_active(
                model,
                chance_cutoff
            )

            if use_model and group_name != "report_only":
                active.append(model)

            if (
                pending
                and use_model
                and group_name != "report_only"
            ):
                pending_models.append(model["name"])

            score_col = model["score_column"]

            site_rank_col = model.get(
                "site_rank_column", ""
            )

            has_score = score_col in data.columns

            has_rank = (
                not site_rank_col
                or site_rank_col in data.columns
            )

            scored_n = (
                int(data[score_col].notna().sum())
                if has_score else 0
            )

            if group_name == "report_only":
                reason = "report only; not counted in consensus"
                use_for_consensus = False
            else:
                use_for_consensus = use_model

            summary_rows.append({
                "group": group_name,
                "model": model["name"],
                "score_column": score_col,
                "site_rank_column": site_rank_col,
                "auroc": model.get(
                    "auroc", "not_measured"
                ),
                "active_for_consensus": use_for_consensus,
                "score_column_present": has_score,
                "site_rank_column_present": has_rank,
                "candidates_with_score": scored_n,
                "candidates_missing_score": (
                    len(data) - scored_n
                    if has_score else len(data)
                ),
                "reason": reason,
            })

        active_groups[group_name] = active

    return (
        active_groups,
        pd.DataFrame(summary_rows),
        pending_models
    )


def count_positive_votes(
    row: pd.Series,
    models: list[dict[str, Any]],
    threshold: float
):
    """Count positive model scores."""

    positives = 0
    missing = []

    for model in models:

        col = model["score_column"]

        value = numeric_value(row.get(col))

        if value is None:
            missing.append(model["name"])

        elif value > threshold:
            positives += 1

    return positives, missing


def main():

    parser = argparse.ArgumentParser(
        description=__doc__
    )

    parser.add_argument(
        "--input",
        required=True,
        type=Path,
        help="Single-mutation report from inspect_scores.py"
    )

    parser.add_argument(
        "--config",
        type=Path,
        default=Path("filters/filtering_config.toml")
    )

    parser.add_argument(
        "--output-dir",
        type=Path,
        default=Path("runs/filter_report")
    )

    args = parser.parse_args()

    # 1. Load configuration

    config = read_config(args.config)

    rules = config.get("filter", {})

    score_threshold = float(
        rules.get("positive_score_threshold", 0.0)
    )

    chance_cutoff = float(
        rules.get("chance_auroc_cutoff", 0.5)
    )

    min_sequence = int(
        rules.get("minimum_sequence_positive", 2)
    )

    min_inverse = int(
        rules.get("minimum_inverse_positive", 2)
    )

    require_top_at_site = bool(
        rules.get("require_top_at_site", True)
    )

    sequence_models = config.get(
        "sequence_models", []
    )

    inverse_models = config.get(
        "inverse_folding_models", []
    )

    report_only_models = config.get(
        "report_only_models", []
    )

    if not sequence_models or not inverse_models:
        raise ValueError(
            "Config must define sequence and inverse-folding models"
        )

    # 2. Load candidate results

    data = read_table(args.input)

    required = {
        "variant_id",
        "mutation_count"
    }

    missing = sorted(required - set(data.columns))

    if missing:
        raise ValueError(
            "Input must be the one-row-per-variant "
            f"report from inspect_scores.py. Missing: {missing}"
        )

    if data["variant_id"].duplicated().any():
        raise ValueError(
            "Duplicate variant_id found. Run inspect_scores.py first."
        )

    # Stage 3 applies to single substitutions only.

    mutation_counts = pd.to_numeric(
        data["mutation_count"],
        errors="coerce"
    )

    if (
        mutation_counts.isna().any()
        or not (mutation_counts == 1).all()
    ):
        raise ValueError(
            "This script filters single substitutions only. "
            "Use score_report_singles, not score_report_pairs."
        )

    # 3. Prepare model groups

    groups = [
        ("sequence", sequence_models),
        ("inverse_folding", inverse_models),
        ("report_only", report_only_models),
    ]

    active, model_summary, pending_models = summarize_models(
        groups,
        data,
        chance_cutoff
    )

    seq_active = active["sequence"]
    inv_active = active["inverse_folding"]

    if (
        len(seq_active) < min_sequence
        or len(inv_active) < min_inverse
    ):
        raise ValueError(
            "Too few active models remain. "
            f"Sequence: {len(seq_active)}, "
            f"inverse folding: {len(inv_active)}."
        )

    # 4. Check whether all required columns exist

    missing_score_columns = [
        model["score_column"]
        for model in seq_active + inv_active
        if model["score_column"] not in data.columns
    ]

    missing_rank_columns = [
        model.get("site_rank_column", "")
        for model in seq_active
        if require_top_at_site
        and model.get("site_rank_column", "")
        not in data.columns
    ]

    # 5. Apply consensus rules to each candidate

    decisions = []

    for _, row in data.iterrows():

        seq_votes, missing_seq = count_positive_votes(
            row,
            seq_active,
            score_threshold
        )

        inv_votes, missing_inv = count_positive_votes(
            row,
            inv_active,
            score_threshold
        )

        top_models = []
        missing_ranks = []

        if require_top_at_site:

            for model in seq_active:

                rank_col = model.get(
                    "site_rank_column", ""
                )

                if not rank_col or rank_col not in data.columns:
                    missing_ranks.append(model["name"])
                    continue

                rank = numeric_value(row.get(rank_col))

                if rank is None:
                    missing_ranks.append(model["name"])

                elif rank == 1:
                    top_models.append(model["name"])

        # Missing evidence is INCOMPLETE, not a negative vote.

        incomplete = bool(
            missing_score_columns
            or missing_rank_columns
            or missing_seq
            or missing_inv
            or missing_ranks
        )

        if incomplete:

            status = "INCOMPLETE"

            reasons = []

            if missing_score_columns:
                reasons.append(
                    "missing score columns: "
                    + ", ".join(missing_score_columns)
                )

            if missing_rank_columns:
                reasons.append(
                    "missing site-rank columns: "
                    + ", ".join(missing_rank_columns)
                )

            if missing_seq:
                reasons.append(
                    "missing sequence scores: "
                    + ", ".join(missing_seq)
                )

            if missing_inv:
                reasons.append(
                    "missing inverse-folding scores: "
                    + ", ".join(missing_inv)
                )

            if missing_ranks:
                reasons.append(
                    "missing site ranks: "
                    + ", ".join(missing_ranks)
                )

            reason = "; ".join(reasons)

        else:

            seq_pass = seq_votes >= min_sequence

            site_pass = (
                len(top_models) >= 1
                if require_top_at_site
                else True
            )

            inv_pass = inv_votes >= min_inverse

            passed = seq_pass and site_pass and inv_pass

            if pending_models:
                status = (
                    "PROVISIONAL_PASS"
                    if passed else "PROVISIONAL_FAIL"
                )
            else:
                status = "PASS" if passed else "FAIL"

            reason = (
                f"sequence positive votes={seq_votes}/{len(seq_active)} "
                f"(minimum {min_sequence}); "
                f"top-at-site models={len(top_models)}; "
                f"inverse-folding positive votes={inv_votes}/{len(inv_active)} "
                f"(minimum {min_inverse})"
            )

            if pending_models:
                reason += (
                    "; AUROC calibration pending: "
                    + ", ".join(pending_models)
                )

        decisions.append({
            "variant_id": row["variant_id"],
            "mutations": row.get("mutations", ""),
            "sequence_positive_votes": seq_votes,
            "sequence_models_active": len(seq_active),
            "top_at_site_models": ";".join(top_models),
            "inverse_positive_votes": inv_votes,
            "inverse_models_active": len(inv_active),
            "filter_status": status,
            "filter_reason": reason,
        })

    # 6. Save reports

    results = pd.DataFrame(decisions)

    args.output_dir.mkdir(
        parents=True,
        exist_ok=True
    )

    decisions_path = (
        args.output_dir / "filter_decisions.csv"
    )

    eligible_path = (
        args.output_dir / "eligible_variants.csv"
    )

    summary_path = (
        args.output_dir / "filter_summary.csv"
    )

    results.to_csv(
        decisions_path,
        index=False
    )

    results[
        results["filter_status"].isin([
            "PASS",
            "PROVISIONAL_PASS"
        ])
    ].to_csv(
        eligible_path,
        index=False
    )

    model_summary.to_csv(
        summary_path,
        index=False
    )

    # 7. Print a short report

    print(f"Input candidates: {len(data)}")

    print("\nFilter status counts:")

    print(
        results["filter_status"]
        .value_counts(dropna=False)
        .to_string()
    )

    print("\nReports saved:")

    print(f"  {decisions_path}")
    print(f"  {eligible_path}")
    print(f"  {summary_path}")

    print(
        "\nStage 3 single-mutation consensus only."
    )

    print(
        "No liabilities, energy, stability or structure vetoes applied."
    )


if __name__ == "__main__":
    main()
