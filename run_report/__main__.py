"""Command line: `build` one run, or rebuild the `leaderboard`."""

from __future__ import annotations

import argparse
import os
import subprocess
import sys
from pathlib import Path

from .config import load_config
from .leaderboard import build_leaderboard, load_summaries
from .load import load_designs
from .render import write_report
from .scorecard import add_deltas, evaluate, find_parent, present_metrics, summarise

DEFAULT_THRESHOLDS = Path(__file__).resolve().parent.parent / "thresholds.toml"


def _commit(code_dir: Path | None) -> str:
    if code_dir is None:
        return ""
    try:
        out = subprocess.run(
            ["git", "-C", str(code_dir), "rev-parse", "--short", "HEAD"],
            capture_output=True, text=True, timeout=20, check=True,
        )
        return out.stdout.strip()
    except (OSError, subprocess.SubprocessError):
        return ""


def build(args: argparse.Namespace) -> int:
    config = load_config(args.thresholds)
    designs = load_designs(args.run_dirs)
    if designs.empty:
        print(f"no design records found under {[str(p) for p in args.run_dirs]}", file=sys.stderr)
        return 2
    keys = [m.key for m in config.metrics if m.key in designs.columns]
    if not keys:
        found = sorted(c for c in designs.columns if c not in ("id", "parent_id", "producer", "step", "path"))
        print(f"none of the thresholds-file metrics are in the data; found: {found}", file=sys.stderr)
        return 2
    scored = designs[designs[keys].notna().any(axis=1)].reset_index(drop=True)
    scored = evaluate(scored, config)
    parent = find_parent(scored, args.parent)
    scored = add_deltas(scored, config, parent)
    history = load_summaries(args.history, exclude=args.run_id) if args.history else []
    summary = summarise(
        scored, config,
        run_id=args.run_id, label=args.label or "",
        code_commit=_commit(args.code_dir), gpu_job_id=os.environ.get("GPU_JOB_ID", ""),
        parent=parent,
    )
    line = write_report(args.out, scored, summary, config, history, parent, max_rows=args.max_rows)
    print(f"wrote {args.out}/report.md ({len(scored)} designs, {len(present_metrics(scored, config))} metrics)")
    print(f"VERDICT: {line}")
    return 0


def leaderboard(args: argparse.Namespace) -> int:
    Path(args.out).write_text(build_leaderboard(args.reports_dir, top=args.top))
    print(f"wrote {args.out} from {len(load_summaries(args.reports_dir))} reports")
    return 0


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(prog="run_report", description=__doc__)
    sub = parser.add_subparsers(dest="command", required=True)

    b = sub.add_parser("build", help="write one run's report")
    b.add_argument("run_dirs", nargs="+", type=Path, help="run directories holding DONE manifests")
    b.add_argument("--run-id", required=True)
    b.add_argument("--out", required=True, type=Path)
    b.add_argument("--label", help="short description shown in the title")
    b.add_argument("--thresholds", type=Path, default=DEFAULT_THRESHOLDS)
    b.add_argument("--history", type=Path, help="reports directory with earlier runs, for comparison")
    b.add_argument("--parent", help="substring of the design id of the parent (wild type) design")
    b.add_argument("--code-dir", type=Path, help="git checkout the run used, to record its commit")
    b.add_argument("--max-rows", type=int, default=500, help="rows kept in designs.csv")
    b.set_defaults(func=build)

    l = sub.add_parser("leaderboard", help="compare all reports in a directory")
    l.add_argument("reports_dir", type=Path)
    l.add_argument("--out", default="LEADERBOARD.md", type=Path)
    l.add_argument("--top", type=int, default=10)
    l.set_defaults(func=leaderboard)

    args = parser.parse_args(argv)
    return args.func(args)


if __name__ == "__main__":
    raise SystemExit(main())
