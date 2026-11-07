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
import os
import platform
import re
import sys
from datetime import date
from pathlib import PurePath
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

FOLDX_STAMP = re.compile(r"(?<!\d)(20\d{2})(0[1-9]|1[0-2])(0[1-9]|[12]\d|3[01])(?!\d)")
"""The expiry date an academic FoldX build carries in its filename, as in
`foldx_20261231`. A convention, not something the binary reports."""

EXPIRY_WARNING_DAYS = 30


def foldx_expiry(binary: str | None = None) -> date | None:
    """The licence expiry in `FOLDX_BINARY`'s filename, or None.

    None means "nothing to say", not "fine": a renamed binary carries no date
    and must not make preflight complain, and must not make it fail either.
    """
    value = (binary if binary is not None else os.environ.get("FOLDX_BINARY", "")).strip()
    if not value:
        return None
    match = FOLDX_STAMP.search(PurePath(value).name)
    if match is None:
        return None
    try:
        return date(*(int(part) for part in match.groups()))
    except ValueError:  # pragma: no cover - the regex already bounds these
        return None


def foldx_notice(today: date | None = None, binary: str | None = None) -> str | None:
    """A warning when the FoldX licence is within 30 days of expiry, or past it.

    A warning rather than a failure, deliberately. The date is read off a
    filename, so a false positive must not stop a campaign - and an academic
    FoldX licence is free to renew but takes a few days, which is exactly the
    kind of thing worth hearing about before the run rather than during it.
    """
    expiry = foldx_expiry(binary)
    if expiry is None:
        return None
    remaining = (expiry - (today or date.today())).days
    if remaining < 0:
        return (
            f"FoldX licence expired {-remaining} day(s) ago ({expiry.isoformat()}, "
            f"read from the FOLDX_BINARY filename). Renew it at "
            f"foldxsuite.crg.eu; `foldx` scores will fail until you do."
        )
    if remaining <= EXPIRY_WARNING_DAYS:
        return (
            f"FoldX licence expires in {remaining} day(s) ({expiry.isoformat()}, "
            f"read from the FOLDX_BINARY filename). An academic renewal is free "
            f"but takes a few days - start it before the next campaign."
        )
    return None


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
        # The tool resolves, but its command may still read a variable that is
        # not set. Reported here, not counted against `resolved`: a keyword
        # argument can stand in for it, so the run is where that is decided.
        unset = [v for v in (route.needs_env if route else ()) if not os.environ.get(v, "").strip()]
        if unset:
            how += f" ({', '.join(unset)} not set)"
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

    notice = foldx_notice()
    if notice is not None:
        print()
        print(f"warning: {notice}")

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
