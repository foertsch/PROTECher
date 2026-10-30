"""Compare every report under a reports directory in one Markdown file."""

from __future__ import annotations

import json
from collections import Counter
from pathlib import Path

import pandas as pd


def load_summaries(reports_dir: Path | str, exclude: str | None = None) -> list[dict]:
    summaries = []
    for path in sorted(Path(reports_dir).glob("*/summary.json")):
        try:
            summary = json.loads(path.read_text())
        except (OSError, json.JSONDecodeError):
            continue
        if summary.get("run_id") != exclude:
            summaries.append(summary)
    return sorted(summaries, key=lambda s: (s["created"], s["run_id"]))


def _bold_best(cells: list[tuple[float | None, int]], better: str) -> list[str]:
    values = [v for v, _ in cells if v is not None]
    best = None if not values else (max(values) if better == "higher" else min(values))
    out = []
    for value, digits in cells:
        if value is None:
            out.append("-")
        else:
            text = f"{value:.{digits}f}"
            out.append(f"**{text}**" if value == best and len(values) > 1 else text)
    return out


def build_leaderboard(reports_dir: Path | str, top: int = 10) -> str:
    reports_dir = Path(reports_dir)
    summaries = load_summaries(reports_dir)
    if not summaries:
        return "# Leaderboard\n\nNo reports yet.\n"

    counts = Counter(k for s in summaries for k in s["metrics"])
    keys = [k for k, _ in sorted(counts.items(), key=lambda kv: (-kv[1], kv[0]))][:5]
    info = {k: next(s["metrics"][k] for s in summaries if k in s["metrics"]) for k in keys}

    columns = {
        k: _bold_best(
            [
                (s["metrics"][k]["median"] if k in s["metrics"] else None, info[k]["digits"])
                for s in summaries
            ],
            info[k]["better"],
        )
        for k in keys
    }
    header = ["Run", "Created", "Models", "Designs", "Pass"] + [f"{info[k]['label']} (median)" for k in keys]
    rows = []
    for i, s in enumerate(summaries):
        rows.append(
            [
                f"[{s['run_id']}](reports/{s['run_id']}/report.md)",
                s["created"],
                ", ".join(s["producers"]),
                str(s["n_designs"]),
                f"{s['n_pass']} ({s['pass_rate']:.0%})",
            ]
            + [columns[k][i] for k in keys]
        )
    lines = ["# Leaderboard", "", "Bold is the best median in a column. One row per run, oldest first.", ""]
    lines += ["| " + " | ".join(header) + " |", "|" + "|".join("---" for _ in header) + "|"]
    lines += ["| " + " | ".join(r) + " |" for r in rows]
    lines += ["", _best_designs(reports_dir, summaries, top)]
    return "\n".join(lines) + "\n"


def _best_designs(reports_dir: Path, summaries: list[dict], top: int) -> str:
    primary = Counter(s["primary"]["key"] for s in summaries).most_common(1)[0][0]
    better = next(s["primary"]["better"] for s in summaries if s["primary"]["key"] == primary)
    frames = []
    for s in summaries:
        csv = reports_dir / s["run_id"] / "designs.csv"
        if csv.is_file():
            frame = pd.read_csv(csv)
            frame.insert(0, "run", s["run_id"])
            frames.append(frame)
    if not frames:
        return ""
    data = pd.concat(frames, ignore_index=True)
    if primary not in data or "pass" not in data:
        return ""
    data = data[data["pass"].astype(bool)].dropna(subset=[primary])
    data = data.sort_values(primary, ascending=(better == "lower")).head(top)
    label = next(s["metrics"][primary]["label"] for s in summaries if primary in s["metrics"])
    rows = [
        f"| {r['run']} | {str(r['id']).replace('|', '/')} | {r['producer']} | {r[primary]:.3g} |"
        for _, r in data.iterrows()
    ]
    return "\n".join(
        [f"## Best passing designs across all runs, by {label}", "", f"| Run | Design | Model | {label} |", "|---|---|---|---|"]
        + rows
    )
