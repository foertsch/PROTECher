"""Grade designs against the thresholds, compare with the parent, and summarise a run."""

from __future__ import annotations

import math
from datetime import datetime, timezone
from typing import Any

import pandas as pd

from .config import Config, Metric


def grade(value: Any, metric: Metric) -> str:
    """"ok" at or beyond the target, "fail" worse than the gate, "warn" in between."""
    if value is None or (isinstance(value, float) and math.isnan(value)):
        return "n/a"
    reaches = (lambda v, t: v >= t) if metric.better == "higher" else (lambda v, t: v <= t)
    if metric.target is not None and reaches(value, metric.target):
        return "ok"
    if metric.gate is not None and not reaches(value, metric.gate):
        return "fail"
    return "warn" if metric.target is not None else "ok"


def present_metrics(df: pd.DataFrame, config: Config) -> list[Metric]:
    return [m for m in config.metrics if m.key in df.columns and df[m.key].notna().any()]


def evaluate(df: pd.DataFrame, config: Config) -> pd.DataFrame:
    """Add a grade column per metric, `pass` (fails no gate) and `n_warn`."""
    out = df.copy()
    failed = pd.Series(False, index=out.index)
    warned = pd.Series(0, index=out.index)
    for metric in present_metrics(out, config):
        grades = out[metric.key].map(lambda v, m=metric: grade(v, m))
        out[f"grade_{metric.key}"] = grades
        failed |= grades.eq("fail")
        warned += grades.eq("warn").astype(int)
    out["pass"] = ~failed
    out["n_warn"] = warned
    return out


def find_parent(df: pd.DataFrame, substring: str | None) -> pd.Series | None:
    """The parent design: every design whose id contains `substring`.

    With several matches (the parent folded with several seeds) the numeric columns are
    averaged, so deltas and the verdict use the mean parent, not one seed. `id` is the first
    match and `matched_ids` joins all of them with "|".
    """
    if not substring or df.empty:
        return None
    hits = df[df["id"].str.contains(substring, regex=False)].sort_values("id")
    if hits.empty:
        return None
    parent = hits.iloc[0].copy()
    if len(hits) > 1:
        for column, value in hits.select_dtypes("number").mean().items():
            parent[column] = value
    parent["matched_ids"] = "|".join(hits["id"])
    return parent


def add_deltas(df: pd.DataFrame, config: Config, parent: pd.Series | None) -> pd.DataFrame:
    """Add `d_<metric>` = design minus parent for every metric the parent has."""
    out = df.copy()
    if parent is None:
        return out
    for metric in present_metrics(out, config):
        value = parent.get(metric.key)
        if value is not None and not pd.isna(value):
            out[f"d_{metric.key}"] = out[metric.key] - value
    return out


def _best(series: pd.Series, better: str) -> float:
    return float(series.max() if better == "higher" else series.min())


def summarise(
    scored: pd.DataFrame,
    config: Config,
    *,
    run_id: str,
    label: str = "",
    code_commit: str = "",
    gpu_job_id: str = "",
    parent: pd.Series | None = None,
) -> dict[str, Any]:
    """The numbers the leaderboard and the next run's comparison need, as plain JSON."""
    metrics: dict[str, dict] = {}
    for metric in present_metrics(scored, config):
        values = scored[metric.key].dropna()
        grades = scored[f"grade_{metric.key}"]
        entry = {
            "label": metric.label,
            "better": metric.better,
            "target": metric.target,
            "gate": metric.gate,
            "digits": metric.digits,
            "n": int(len(values)),
            "median": float(values.median()),
            "mean": float(values.mean()),
            "min": float(values.min()),
            "max": float(values.max()),
            "best": _best(values, metric.better),
            "ok": int(grades.eq("ok").sum()),
            "warn": int(grades.eq("warn").sum()),
            "fail": int(grades.eq("fail").sum()),
        }
        if parent is not None and not pd.isna(parent.get(metric.key)):
            entry["parent"] = float(parent[metric.key])
        metrics[metric.key] = entry
    n = int(len(scored))
    n_pass = int(scored["pass"].sum())
    seeds = sorted({int(s) for s in scored.get("meta_seed", pd.Series(dtype=float)).dropna()})
    groups = scored[config.group].value_counts().to_dict() if config.group in scored else {}
    primary = config.primary
    return {
        "run_id": run_id,
        "label": label,
        "created": datetime.now(timezone.utc).strftime("%Y-%m-%d %H:%M UTC"),
        "code_commit": code_commit,
        "gpu_job_id": gpu_job_id,
        "n_designs": n,
        "n_pass": n_pass,
        "pass_rate": (n_pass / n) if n else 0.0,
        "n_with_warnings": int((scored["n_warn"] > 0).sum()),
        "producers": sorted(scored["producer"].dropna().unique().tolist()),
        "seeds": seeds,
        "groups": {str(k): int(v) for k, v in groups.items()},
        "primary": {"key": primary.key, "label": primary.label, "better": primary.better},
        "parent_id": None if parent is None else str(parent["id"]),
        "metrics": metrics,
    }


def verdict(summary: dict, history: list[dict]) -> str:
    """One plain line a reader can judge the run by."""
    parts = [f"{summary['n_pass']}/{summary['n_designs']} designs pass all gates"]
    stats = summary["metrics"].get(summary["primary"]["key"])
    if stats:
        d = stats["digits"]
        parts.append(f"best {stats['label']} {stats['best']:.{d}f}")
        if "parent" in stats:
            parts.append(f"parent {stats['parent']:.{d}f}")
    if summary["n_with_warnings"]:
        parts.append(f"{summary['n_with_warnings']} with warnings")
    if history:
        best = max(history, key=lambda h: (h["pass_rate"], h["n_designs"]))
        parts.append(f"earlier best run {best['run_id']}: {best['n_pass']}/{best['n_designs']}")
    return " · ".join(parts)
