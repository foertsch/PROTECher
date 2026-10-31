"""
rank_ipsae.py -- aggregate per-fold ipSAE rows into per-design rankings.

Consumes the tidy table of IpsaeRow records (one row per design x seed) written
by ipsae_score.score_fold, plus the PARENT scored the SAME way (same construct,
cutoffs, ipsae_commit). This is the "dataset of all rankings" layer.

Decision logic (consolidated-plan option c): pose-gate -> aggregate over seeds
-> veto vs parent -> rank. Nothing here re-folds, so changing a rule is cheap.
"""
from __future__ import annotations

import pandas as pd


def aggregate(rows: pd.DataFrame, min_pose_seeds: int = 3) -> pd.DataFrame:
    """One row per design. Robust ipSAE = median over pose-passing seeds.

    If no external pose gate was supplied (pose_pass all None), fall back to
    "fold has contacts" as the gate. Prefer feeding site_guard's verdict in.
    """
    rows = rows.copy()
    if "pose_pass" not in rows or rows["pose_pass"].isna().all():
        rows["pose_pass"] = rows["n_contacts"] > 0

    passed = rows[rows["pose_pass"].astype(bool)]
    agg = passed.groupby("design_id").agg(
        ipsae_robust=("ipsae_min", "median"),   # robust central estimate
        ipsae_best=("ipsae_min", "max"),         # best single seed (optimistic)
        n_pose_seeds=("seed", "nunique"),
    )
    pose_frac = (rows.groupby("design_id")["pose_pass"]
                     .mean().rename("pose_pass_frac"))
    out = agg.join(pose_frac)
    out = out[out["n_pose_seeds"] >= min_pose_seeds]   # drop one-lucky-seed designs
    return out.reset_index()


def rank(designs: pd.DataFrame,
         parent_ipsae_robust: float,
         veto_margin: float = 0.0) -> pd.DataFrame:
    """Parent veto + rank survivors by robust ipSAE (descending).

    veto_margin loosens the veto to absorb seed noise (parent spread ~0.007
    on ipTM; set a few hundredths if you want to keep borderline designs).
    """
    d = designs.copy()
    d["veto_pass"] = d["ipsae_robust"] >= (parent_ipsae_robust - veto_margin)
    d = d.sort_values("ipsae_robust", ascending=False).reset_index(drop=True)
    d["rank"] = d.index + 1
    return d


def parent_robust(rows: pd.DataFrame, parent_id: str = "parent",
                  min_pose_seeds: int = 3) -> float:
    """Convenience: robust ipSAE of the parent from the same table."""
    agg = aggregate(rows[rows["design_id"] == parent_id],
                    min_pose_seeds=min_pose_seeds)
    if agg.empty:
        raise ValueError(f"parent '{parent_id}' has < {min_pose_seeds} pose-passing seeds")
    return float(agg["ipsae_robust"].iloc[0])
