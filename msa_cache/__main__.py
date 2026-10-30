"""Command line: build the parent's cache, import an existing Boltz run's alignments, or show a cache."""

from __future__ import annotations

import argparse
import json
from pathlib import Path

from .cache import MsaCache, build_parent


def _chains(fasta: Path) -> list[str]:
    records = [r for r in fasta.read_text().split(">") if r.strip()]
    if len(records) != 1:
        raise SystemExit(f"{fasta} must hold exactly one record, with chains separated by ':'")
    body = records[0].partition("\n")[2].replace("\n", "")
    return [c.strip() for c in body.split(":") if c.strip()]


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(prog="msa_cache", description=__doc__)
    sub = parser.add_subparsers(dest="command", required=True)

    b = sub.add_parser("build", help="fold the parent once with the public server and keep its alignments (GPU)")
    b.add_argument("parent_fasta", type=Path, help="one record, chains separated by ':'")
    b.add_argument("--cache", required=True, type=Path)
    b.add_argument("--work-dir", type=Path, help="where the parent fold goes (default: <cache>/parent_fold)")
    b.add_argument("--seed", type=int, default=0)

    i = sub.add_parser("import", help="keep the alignments of an existing Boltz run for the parent")
    i.add_argument("msa_dir", type=Path, help="the run's boltz_results_<name>/msa folder")
    i.add_argument("parent_fasta", type=Path)
    i.add_argument("--cache", required=True, type=Path)

    s = sub.add_parser("info", help="show what a cache holds")
    s.add_argument("--cache", required=True, type=Path)

    args = parser.parse_args(argv)
    cache = MsaCache(args.cache)
    if args.command == "build":
        build_parent(cache, _chains(args.parent_fasta), args.work_dir or args.cache / "parent_fold", seed=args.seed)
    elif args.command == "import":
        cache.import_boltz_msas(args.msa_dir, _chains(args.parent_fasta))
    print(json.dumps(cache.info(), indent=1))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
