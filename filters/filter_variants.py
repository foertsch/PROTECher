#!/usr/bin/env python3
"""Apply the PROTECher Track A Stage 3 single-substitution consensus screen.

Input is the one-row-per-variant CSV produced by inspect_scores.py.
Expected model groups and columns are listed in filtering_config.toml.
Missing columns/scores are labelled INCOMPLETE, never converted to negative
votes. This script refuses dummy/test inputs and outputs.

This implements Stage 3 only. It does not run model scorers, inspect structural
liabilities, calculate FoldX/Rosetta energy, test stability/expression, build
5-10 mutation designs, or apply the Stage 6 structure veto.
"""

from __future__ import annotations

import argparse
import math
import tomllib
from pathlib import Path
from typing import Any

import pandas as pd


def reject_dummy_path(path: Path, label: str) -> None:
    if any("dummy" in part.lower() or "test" in part.lower() for part in path.parts):
        raise ValueError(
            f"Refusing {label} path containing 'dummy' or 'test': {path}. "
            "Use actual, non-dummy model results."
        )


def read_table(path: Path) -> pd.DataFrame:
    if not path.is_file():
        raise FileNotFoundError(f"Input file not found: {path}")
    reject_dummy_path(path, "input")
    if path.suffix.lower() == ".csv":
        frame = pd.read_csv(path)
    elif path.suffix.lower() == ".parquet":
        frame = pd.read_parquet(path)
    else:
        raise ValueError("Input must be .csv or .parquet")
    dummy_columns = [c for c in frame.columns if "dummy" in str(c).lower()]
    if dummy_columns:
        raise ValueError(
            "Dummy score columns are not accepted: " + ", ".join(dummy_columns)
        )
    return frame


def read_config(path: Path) -> dict[str, Any]:
    if not path.is_file():
        raise FileNotFoundError(f"Config file not found: {path}")
    with path.open("rb") as handle:
        return tomllib.load(handle)


def to_number(value: Any) -> float | None:
    parsed = pd.to_numeric(pd.Series([value]), errors="coerce").iloc[0]
    if pd.isna(parsed):
        return None
    return float(parsed)


def model_calibration(model: dict[str, Any], cutoff: float) -> tuple[bool, str, bool]:
    """Return (included, explanation, calibration_pending)."""
    if not model.get("enabled", True):
        return False, "disabled in configuration", False

    raw_auc = model.get("auroc", "not_measured")
    try:
        auc = float(raw_auc)
    except (TypeError, ValueError):
        return True, "AUROC not measured; provisionally included", True

    if not math.isfinite(auc) or not 0 <= auc <= 1:
        return True, "AUROC unavailable/invalid; provisionally included", True
    if auc <= cutoff:
        return False, f"excluded: AUROC {auc:.3f} <= chance cutoff {cutoff:.3f}", False
    return True, f"included: AUROC {auc:.3f} > chance cutoff {cutoff:.3f}", False


def get_active_models(
    models: list[dict[str, Any]], data: pd.DataFrame, group: str, cutoff: float
) -> tuple[list[dict[str, Any]], list[dict[str, Any]], list[str]]:
    active: list[dict[str, Any]] = []
    summary: list[dict[str, Any]] = []
    pending: list[str] = []

    for model in models:
        include, reason, calibration_pending = model_calibration(model, cutoff)
        if include:
            active.append(model)
        if include and calibration_pending:
            pending.append(model["name"])

        score_col = model["score_column"]
        rank_col = model.get("site_rank_column", "")
        score_present = score_col in data.columns
        rank_present = not rank_col or rank_col in data.columns
        n_scored = int(data[score_col].notna().sum()) if score_present else 0

        summary.append({
            "group": group,
            "model": model["name"],
            "score_column": score_col,
            "site_rank_column": rank_col,
            "configured_enabled": bool(model.get("enabled", True)),
            "included_in_consensus": include,
            "auroc": model.get("auroc", "not_measured"),
            "repo_status": model.get("repo_status", "not specified"),
            "score_column_present": score_present,
            "site_rank_column_present": rank_present,
            "candidates_with_score": n_scored,
            "candidates_missing_score": len(data) - n_scored if score_present else len(data),
            "calibration_pending": calibration_pending,
            "reason": reason,
        })

    return active, summary, pending


def count_votes(
    row: pd.Series, models: list[dict[str, Any]], threshold: float
) -> tuple[int, list[str]]:
    positives = 0
    missing: list[str] = []
    for model in models:
        value = to_number(row.get(model["score_column"]))
        if value is None:
            missing.append(model["name"])
        elif value > float(model.get("positive_score_threshold", threshold)):
            positives += 1
    return positives, missing


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input", required=True, type=Path,
                        help="Real, one-row-per-variant report from inspect_scores.py")
    parser.add_argument("--config", type=Path,
                        default=Path("filters/filtering_config.toml"))
    parser.add_argument("--output-dir", type=Path,
                        default=Path("runs/her2_v1/reports/filtering"))
    args = parser.parse_args()

    reject_dummy_path(args.output_dir, "output")
    data = read_table(args.input)
    config = read_config(args.config)
    rules = config.get("filter", {})

    required = {"variant_id", "mutation_count"}
    missing_required = sorted(required - set(data.columns))
    if missing_required:
        raise ValueError(
            "Input must be generated by inspect_scores.py; missing columns: "
            + ", ".join(missing_required)
        )
    if data["variant_id"].isna().any() or data["variant_id"].duplicated().any():
        raise ValueError("variant_id must be present and unique (one row per candidate)")

    mutation_counts = pd.to_numeric(data["mutation_count"], errors="coerce")
    if mutation_counts.isna().any() or not (mutation_counts == 1).all():
        raise ValueError(
            "Track A Stage 3 is a single-substitution scan. Supply the singles "
            "report, not the pairs/multi-mutant report."
        )

    threshold = float(rules.get("positive_score_threshold", 0.0))
    chance_cutoff = float(rules.get("chance_auroc_cutoff", 0.5))
    min_sequence = int(rules.get("minimum_sequence_positive", 2))
    min_inverse = int(rules.get("minimum_inverse_positive", 2))
    require_top = bool(rules.get("require_top_at_site", True))

    seq_models = config.get("sequence_models", [])
    inv_models = config.get("inverse_folding_models", [])
    report_only = config.get("report_only_models", [])
    if not seq_models or not inv_models:
        raise ValueError("Config must define sequence_models and inverse_folding_models")

    seq_active, seq_summary, seq_pending = get_active_models(
        seq_models, data, "sequence", chance_cutoff
    )
    inv_active, inv_summary, inv_pending = get_active_models(
        inv_models, data, "inverse_folding", chance_cutoff
    )
    _, report_summary, _ = get_active_models(
        report_only, data, "report_only", chance_cutoff
    ) if report_only else ([], [], [])

    model_summary = pd.DataFrame(seq_summary + inv_summary + report_summary)

    # Missing whole columns are reported before candidate-level vote logic.
    missing_score_columns = sorted({
        model["score_column"]
        for model in seq_active + inv_active
        if model["score_column"] not in data.columns
    })
    missing_rank_columns = sorted({
        model.get("site_rank_column", "")
        for model in seq_active
        if require_top and model.get("site_rank_column")
        and model["site_rank_column"] not in data.columns
    })

    # If calibration or missing columns leave too few active models, we still
    # produce an INCOMPLETE report rather than stopping or inventing votes.
    group_incomplete = (
        len(seq_active) < min_sequence
        or len(inv_active) < min_inverse
        or bool(missing_score_columns)
        or bool(missing_rank_columns)
    )

    decisions: list[dict[str, Any]] = []
    for _, row in data.iterrows():
        seq_votes, missing_seq = count_votes(row, seq_active, threshold)
        inv_votes, missing_inv = count_votes(row, inv_active, threshold)

        top_models: list[str] = []
        missing_ranks: list[str] = []
        if require_top:
            for model in seq_active:
                rank_col = model.get("site_rank_column", "")
                if not rank_col or rank_col not in data.columns:
                    missing_ranks.append(model["name"])
                    continue
                rank = to_number(row.get(rank_col))
                if rank is None:
                    missing_ranks.append(model["name"])
                elif rank == 1:
                    top_models.append(model["name"])

        # Any missing required score/rank or insufficient model group means
        # the candidate cannot be assessed; it is not labelled a biological FAIL.
        incomplete = (
            group_incomplete or bool(missing_seq) or bool(missing_inv)
            or bool(missing_ranks)
        )

        if incomplete:
            status = "INCOMPLETE"
            reasons = []
            if len(seq_active) < min_sequence:
                reasons.append(
                    f"only {len(seq_active)} active sequence models; need at least {min_sequence}"
                )
            if len(inv_active) < min_inverse:
                reasons.append(
                    f"only {len(inv_active)} active inverse-folding models; need at least {min_inverse}"
                )
            if missing_score_columns:
                reasons.append("missing score columns: " + ", ".join(missing_score_columns))
            if missing_rank_columns:
                reasons.append("missing site-rank columns: " + ", ".join(missing_rank_columns))
            if missing_seq:
                reasons.append("missing candidate sequence scores: " + ", ".join(missing_seq))
            if missing_inv:
                reasons.append("missing candidate inverse-folding scores: " + ", ".join(missing_inv))
            if missing_ranks:
                reasons.append("missing candidate site ranks: " + ", ".join(missing_ranks))
            reason = "; ".join(reasons)
        else:
            seq_pass = seq_votes >= min_sequence
            site_pass = bool(top_models) if require_top else True
            inv_pass = inv_votes >= min_inverse
            passed = seq_pass and site_pass and inv_pass
            calibration_pending = bool(seq_pending or inv_pending)
            if calibration_pending:
                status = "PROVISIONAL_PASS" if passed else "PROVISIONAL_FAIL"
            else:
                status = "PASS" if passed else "FAIL"
            reason = (
                f"sequence positive votes {seq_votes}/{len(seq_active)} "
                f"(minimum {min_sequence}); top-at-site models {len(top_models)}; "
                f"inverse-folding positive votes {inv_votes}/{len(inv_active)} "
                f"(minimum {min_inverse})"
            )
            if calibration_pending:
                reason += "; AUROC calibration pending for: " + ", ".join(
                    seq_pending + inv_pending
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

    result = pd.DataFrame(decisions)
    args.output_dir.mkdir(parents=True, exist_ok=True)
    decisions_path = args.output_dir / "filter_decisions.csv"
    eligible_path = args.output_dir / "eligible_variants.csv"
    summary_path = args.output_dir / "filter_model_summary.csv"

    result.to_csv(decisions_path, index=False)
    result[result["filter_status"].isin(["PASS", "PROVISIONAL_PASS"])].to_csv(
        eligible_path, index=False
    )
    model_summary.to_csv(summary_path, index=False)

    print(f"Candidates evaluated: {len(data)}")
    print("Filter status counts:")
    print(result["filter_status"].value_counts(dropna=False).to_string())
    print("\nModel availability:")
    if not model_summary.empty:
        print(model_summary[["model", "score_column_present", "candidates_with_score", "reason"]].to_string(index=False))
    print("\nReports saved:")
    print(f"  {decisions_path}")
    print(f"  {eligible_path}")
    print(f"  {summary_path}")
    print("\nStage 3 single-substitution consensus only; no Stage 4-6 filters applied.")


if __name__ == "__main__":
    main()
