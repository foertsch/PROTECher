"""Write report.md, designs.csv, summary.json and two plots for one run."""

from __future__ import annotations

import json
from pathlib import Path

import pandas as pd

from .config import Config, Metric
from .scorecard import present_metrics, verdict

MARK = {"ok": "✓", "warn": "!", "fail": "✗", "n/a": ""}


def _fmt(value, digits: int) -> str:
    return "-" if value is None or pd.isna(value) else f"{value:.{digits}f}"


def _cell(value, grade: str, digits: int) -> str:
    return f"{_fmt(value, digits)} {MARK.get(grade, '')}".strip()


def _short(text: str, width: int = 44) -> str:
    text = str(text).replace("|", "/")
    return text if len(text) <= width else text[: width - 1] + "…"


def _table(header: list[str], rows: list[list[str]]) -> str:
    lines = ["| " + " | ".join(header) + " |", "|" + "|".join("---" for _ in header) + "|"]
    lines += ["| " + " | ".join(row) + " |" for row in rows]
    return "\n".join(lines)


def _sorted(scored: pd.DataFrame, config: Config, metrics: list[Metric]) -> pd.DataFrame:
    primary = next((m for m in metrics if m.key == config.x), metrics[0])
    return scored.sort_values(
        ["pass", primary.key, "id"],
        ascending=[False, primary.better == "lower", True],
        na_position="last",
    )


def scorecard_table(summary: dict, history: list[dict]) -> str:
    rows = []
    for key, s in summary["metrics"].items():
        d = s["digits"]
        earlier = [h["metrics"][key] for h in history if key in h["metrics"]]
        if earlier:
            pick = max if s["better"] == "higher" else min
            best = pick(earlier, key=lambda e: e["median"])
            earlier_cell = _fmt(best["median"], d)
        else:
            earlier_cell = "-"
        rows.append(
            [
                s["label"],
                _fmt(s["median"], d),
                f"{_fmt(s['min'], d)} to {_fmt(s['max'], d)}",
                _fmt(s.get("parent"), d),
                earlier_cell,
                f"{s['ok']} / {s['warn']} / {s['fail']}",
            ]
        )
    header = ["Metric", "Median", "Range", "Parent", "Best earlier median", "ok / warn / fail"]
    return _table(header, rows)


def _blank(value) -> bool:
    return value is None or (isinstance(value, float) and pd.isna(value))


def designs_table(
    scored: pd.DataFrame, config: Config, metrics: list[Metric], limit: int = 30, parent_ids: frozenset[str] = frozenset()
) -> str:
    ordered = _sorted(scored, config, metrics).head(limit)
    primary = next((m for m in metrics if m.key == config.x), metrics[0])
    has_delta = f"d_{primary.key}" in ordered
    has_group = config.group in ordered
    has_mut = "meta_mutations" in ordered
    header = ["Design", "Model"]
    header += ["Group"] if has_group else []
    header += ["Mutations"] if has_mut else []
    header += [m.label for m in metrics]
    header += [f"Δ {primary.label}"] if has_delta else []
    header += ["Result"]
    rows = []
    for _, r in ordered.iterrows():
        row = [_short(r["id"]) + (" (parent)" if r["id"] in parent_ids else ""), str(r["producer"])]
        row += ["-" if _blank(r[config.group]) else str(r[config.group])] if has_group else []
        row += ["-" if _blank(r["meta_mutations"]) else _short(r["meta_mutations"], 60)] if has_mut else []
        for m in metrics:
            row.append(_cell(r.get(m.key), r.get(f"grade_{m.key}", "n/a"), m.digits))
        if has_delta:
            delta = r.get(f"d_{primary.key}")
            row.append("-" if _blank(delta) else f"{delta:+.{primary.digits}f}")
        row.append("pass" if r["pass"] else "FAIL")
        rows.append(row)
    shown = f"Showing {len(ordered)} of {len(scored)} designs, passing first, best {config.primary.label} first."
    return _table(header, rows) + "\n\n" + shown


def history_table(history: list[dict], summary: dict) -> str:
    key = summary["primary"]["key"]
    rows = []
    for h in history + [summary]:
        s = h["metrics"].get(key)
        d = s["digits"] if s else 2
        name = h["run_id"] + (" (this run)" if h is summary else "")
        rows.append(
            [
                name,
                h["created"],
                ", ".join(h["producers"]),
                str(h["n_designs"]),
                f"{h['n_pass']} ({h['pass_rate']:.0%})",
                _fmt(s["median"], d) if s else "-",
                _fmt(s["best"], d) if s else "-",
            ]
        )
    label = summary["primary"]["label"]
    header = ["Run", "Created", "Models", "Designs", "Pass", f"{label} median", f"{label} best"]
    return _table(header, rows)


def _pyplot():
    """pyplot with the Agg backend, or None when matplotlib is not installed.

    matplotlib is the `notebooks` extra, so a `[dev]` install has none. The
    report is built out of the tables and the figures are an addition to it, so
    a missing matplotlib leaves the image out - the same answer `_scatter` and
    `_compare` already give when the columns they need are absent. It used to
    raise ModuleNotFoundError here instead, which failed five tests on exactly
    the install INSTALL.md section 2 tells you to make.
    """
    try:
        import matplotlib
    except ImportError:
        return None

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    return plt


def _scatter(scored: pd.DataFrame, config: Config, parent: pd.Series | None, path: Path) -> bool:
    plt = _pyplot()
    if plt is None:
        return False

    x, y = config.x, config.y
    if x not in scored or y not in scored or x == y:
        return False
    data = scored.dropna(subset=[x, y])
    if data.empty:
        return False
    by = config.group if config.group in data else "producer"
    fig, ax = plt.subplots(figsize=(6, 4.5), dpi=100)
    others = data if parent is None else data[data["id"] != parent["id"]]
    for name, part in others.groupby(by):
        ax.scatter(part[x], part[y], s=36, alpha=0.8, label=str(name))
    if parent is not None and not pd.isna(parent.get(x)) and not pd.isna(parent.get(y)):
        ax.scatter([parent[x]], [parent[y]], marker="*", s=220, c="black", label="parent", zorder=5)
    for metric, line in ((config.by_key(x), ax.axvline), (config.by_key(y), ax.axhline)):
        for value, style in ((metric.target, "--"), (metric.gate, ":")):
            if value is not None:
                line(value, color="grey", linestyle=style, linewidth=1)
    ax.set_xlabel(config.by_key(x).label)
    ax.set_ylabel(config.by_key(y).label)
    ax.legend(fontsize=8, loc="best")
    ax.set_title("Designs (dashed: target, dotted: gate)", fontsize=10)
    fig.tight_layout()
    fig.savefig(path)
    plt.close(fig)
    return True


def _compare(scored: pd.DataFrame, config: Config, history: list[dict], summary: dict, path: Path) -> bool:
    plt = _pyplot()
    if plt is None:
        return False

    key = config.x
    if key not in scored:
        return False
    metric = config.by_key(key)
    fig, ax = plt.subplots(figsize=(6, 3.8), dpi=100)
    labels = []
    for i, h in enumerate(history):
        s = h["metrics"].get(key)
        if not s:
            continue
        ax.vlines(i, s["min"], s["max"], color="grey", linewidth=3)
        ax.scatter([i], [s["median"]], color="black", zorder=3)
        labels.append(h["run_id"])
    pos = len(labels)
    values = scored[key].dropna()
    offsets = [((j * 7919) % 100 - 50) / 400 for j in range(len(values))]
    ax.scatter([pos + o for o in offsets], values, s=22, alpha=0.8, color="tab:blue", zorder=3)
    ax.scatter([pos], [values.median()], marker="D", color="black", zorder=4)
    labels.append(summary["run_id"] + "\n(this run)")
    for value, style in ((metric.target, "--"), (metric.gate, ":")):
        if value is not None:
            ax.axhline(value, color="grey", linestyle=style, linewidth=1)
    ax.set_xticks(range(len(labels)))
    ax.set_xticklabels(labels, fontsize=8, rotation=20, ha="right")
    ax.set_ylabel(metric.label)
    ax.set_title(f"{metric.label}: range and median per run, this run's designs in blue", fontsize=9)
    fig.tight_layout()
    fig.savefig(path)
    plt.close(fig)
    return True


def write_report(
    out_dir: Path,
    scored: pd.DataFrame,
    summary: dict,
    config: Config,
    history: list[dict],
    parent: pd.Series | None,
    *,
    max_rows: int = 500,
) -> str:
    """Write the report files and return the one-line verdict."""
    out_dir.mkdir(parents=True, exist_ok=True)
    metrics = present_metrics(scored, config)
    line = verdict(summary, history)
    summary = {**summary, "verdict": line}

    files = ["report.md", "designs.csv", "summary.json"]
    if _scatter(scored, config, parent, out_dir / "scatter.png"):
        files.append("scatter.png")
    if _compare(scored, config, history, summary, out_dir / "compare.png"):
        files.append("compare.png")
    summary["files"] = files

    _sorted(scored, config, metrics).head(max_rows).to_csv(out_dir / "designs.csv", index=False, float_format="%.5g")
    # encoding spelled out on every write in this package: the report carries
    # U+2717 and U+00B7, and `write_text` without it uses the platform default,
    # which is cp1252 on Windows - building a report there raised
    # UnicodeEncodeError and took five tests with it. Linux never saw it because
    # PEP 538 coerces the C locale to UTF-8, so this is a Windows-only break on
    # exactly the machine INSTALL.md section 2 tells you to run the suite on.
    (out_dir / "summary.json").write_text(json.dumps(summary, indent=1), encoding="utf-8")

    parts = [
        f"# Run {summary['run_id']}" + (f": {summary['label']}" if summary["label"] else ""),
        f"**{line}**",
        f"Created {summary['created']} · models {', '.join(summary['producers'])} · "
        f"seeds {', '.join(map(str, summary['seeds'])) or 'n/a'} · code `{summary['code_commit'] or 'unknown'}`"
        + (f" · GPU job {summary['gpu_job_id']}" if summary["gpu_job_id"] else ""),
        "ok = reaches the target, ! = passes the gate but misses the target, ✗ = fails the gate. "
        "Thresholds are in `thresholds.toml`.",
        "## Scorecard",
        scorecard_table(summary, history),
        "## Designs",
        designs_table(
            scored, config, metrics,
            parent_ids=frozenset() if parent is None else frozenset(str(parent["matched_ids"]).split("|")),
        ),
    ]
    if "scatter.png" in files:
        parts += ["## Two metrics at once", "![scatter](scatter.png)"]
    parts += ["## Compared with earlier runs", history_table(history, summary)]
    if "compare.png" in files:
        parts.append("![compare](compare.png)")
    if summary["groups"]:
        parts += ["## Groups", ", ".join(f"{k}: {v}" for k, v in summary["groups"].items())]
    parts += [
        "## Files",
        "`designs.csv` (one row per design with every score and grade), `summary.json` (numbers for the leaderboard).",
    ]
    (out_dir / "report.md").write_text("\n\n".join(parts) + "\n", encoding="utf-8")
    return line
