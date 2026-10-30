"""One table: every adapter, its step, whether it resolves here, and how.

    python -m design_pipeline.preflight

Run it before the event, in each environment, and read the `env` column: a row
marked A that does not resolve inside Env A is a problem to fix now rather than
at 2am. Nothing is pulled, imported for real, or executed - the probes are
`importlib.util.find_spec` and `shutil.which`, so this is fast and safe to run
on a login node.

Exit status is 1 when any adapter in the default campaign - the five steps with
their default models - does not resolve, so CI can gate on it. An unresolved
optional model is reported and does not fail the run.
"""

from __future__ import annotations

import argparse
import platform
import sys
from typing import Iterable, Sequence

from .adapters._gates import DEFAULT_CHECKS
from .adapters._scorers import DEFAULT_SCORERS
from .registry import REGISTRY
from .routes import EXTRAS, VERIFIED_ON, resolve, route_for

DEFAULT_CAMPAIGN: tuple[str, ...] = (
    "pdb_download",
    *DEFAULT_SCORERS,
    *DEFAULT_CHECKS,
)
"""What a campaign touches with no model named: step 1 plus the default scorer
and check sets. These are the rows that decide the exit status."""

COLUMNS = ("adapter", "step", "env", "resolved", "how")


def rows(names: Iterable[str]) -> list[tuple[str, ...]]:
    out: list[tuple[str, ...]] = []
    for name in names:
        route = route_for(name)
        ok, how = resolve(name)
        step = route.step if route else "?"
        env = route.env if route else "?"
        mark = "yes" if ok else "NO"
        if route is not None and not route.verified and route.how:
            # A route nobody has checked is not the same as a checked one, and
            # the table is the only place that distinction is visible.
            mark += " (route unchecked)"
        out.append((name, step, env, mark, how))
    return out


def render(table: Sequence[tuple[str, ...]], header: Sequence[str] = COLUMNS) -> str:
    widths = [
        max(len(str(row[i])) for row in (header, *table)) for i in range(len(header))
    ]
    lines = ["  ".join(str(h).ljust(w) for h, w in zip(header, widths)).rstrip()]
    lines.append("  ".join("-" * w for w in widths))
    for row in table:
        lines.append("  ".join(str(c).ljust(w) for c, w in zip(row, widths)).rstrip())
    return "\n".join(lines)


def notes(names: Iterable[str]) -> list[str]:
    out = []
    for name in names:
        route = route_for(name)
        if route is not None and route.note:
            out.append(f"{name}: {route.note}")
    return out


def main(argv: Sequence[str] | None = None) -> int:
    parser = argparse.ArgumentParser(
        prog="python -m design_pipeline.preflight", description=__doc__
    )
    parser.add_argument(
        "--all",
        action="store_true",
        help="include the dummy fixtures, which always resolve",
    )
    parser.add_argument(
        "--notes", action="store_true", help="print the per-tool notes as well"
    )
    args = parser.parse_args(argv)

    names = [n for n in REGISTRY if args.all or not n.startswith("dummy_")]
    extras = list(EXTRAS)

    print(f"design_pipeline preflight - routes last verified {VERIFIED_ON}")
    print(f"python {platform.python_version()} on {sys.platform}, {platform.machine()}")
    print()
    print(render(rows(names)))
    print()
    print("not adapters, but part of an environment:")
    print(render(rows(extras)))

    unresolved = [n for n in DEFAULT_CAMPAIGN if not resolve(n)[0]]
    print()
    if unresolved:
        print(
            f"default campaign NOT ready here: {', '.join(unresolved)} "
            f"{'does' if len(unresolved) == 1 else 'do'} not resolve."
        )
        print("Each row's `how` says what was probed; INSTALL.md has the routes.")
    else:
        print("default campaign resolves here.")

    if args.notes:
        print()
        for line in notes([*names, *extras]):
            print(f"  {line}")

    return 1 if unresolved else 0


if __name__ == "__main__":  # pragma: no cover - entry point
    raise SystemExit(main())
