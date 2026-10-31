"""Turn chosen designs into scFv sequences and a submission file.

The last step of the Track A chain. `pick_designs` returns designs as crystal-numbered substitutions (`B:S97A+A:Q27E`:
heavy chain B, light chain A, 1-based along the chain as given). The sequence that is folded and submitted is the scFv.
This script applies each design to the parent through `design_pipeline.scfv`, checks the result, and writes the files the
next steps read. It checks everything before it writes anything, and says every problem it found, not only the first.

    python examples/track_a_export.py --parent parent.csv --designs designs.parquet --structure 1N8Z.cif --out submission/
    python examples/track_a_export.py --parent parent.csv --ids B:S97A+B:G103Y+B:F104Y+A:Q27E+A:H91R --out submission/

`--parent` takes a FASTA file, a `name,sequence` file (one line, like the parent the challenge gives) or a bare sequence.
`--designs` takes the table `pick_designs` wrote (a parquet or a csv with a `variant_id` column, ordered by its `design`
column when it has one), or give ids with `--ids`. `--structure` is the crystal; with it the parent is checked against
the crystal's variable domains, without it that check is skipped and the output says so.

Checked per design: every wild-type residue matches the parent at its mapped position, every position is in a variable
domain (a constant-domain position has no scFv counterpart), the design has at least `--min-substitutions` substitutions
(5 in the Track A rules) counted on the finished sequence, the linker is untouched, no two designs give the same
sequence, and the sequence carries no Twin-Strep tag (the lab appends its own).

Written to `--out`:
    submission.csv           name,sequence for each design, no header unless `--header`. The parent is not in it.
    scfv_for_folding.fasta   the parent first, then each design: the form `examples/track_a_fold.py --scfv` reads.
    designs_report.csv       name, crystal ids, substitution counts per chain, the substitutions in scFv numbering.

Assumption to confirm against the challenge's own instructions: that a submission is `name,sequence` lines like the
parent line it gives. Change `write_submission` if the format differs.
"""

from __future__ import annotations

import argparse
import csv
import sys
from pathlib import Path
from typing import NamedTuple, Sequence

from design_pipeline import scfv
from design_pipeline.variants import make_id, parse_id

AMINO_ACIDS = frozenset("ACDEFGHIKLMNPQRSTVWY")
STREP_MOTIF = "WSHPQFEK"
"""Adaptyv appends the Twin-Strep tag itself. A sequence that already carries it would get it twice."""


class Design(NamedTuple):
    name: str
    variant_id: str
    sequence: str
    substitutions: list[tuple[int, str, str]]  # (scFv position, parent residue, design residue)


def read_parent(source: str) -> tuple[str, str]:
    """(name, sequence) of the parent from a FASTA file, a `name,sequence` file or line, or a bare sequence."""
    try:
        path = Path(source) if Path(source).is_file() else None
    except OSError:  # a bare sequence can be too long to be a file name
        path = None
    text = path.read_text() if path is not None else source
    lines = [line.strip() for line in text.splitlines() if line.strip()]
    if lines and lines[0].replace(" ", "").lower() in ("name,sequence", "id,sequence"):
        lines = lines[1:]  # a header line
    if not lines:
        raise ValueError("the parent is empty")
    if lines[0].startswith(">"):
        name = (lines[0][1:].split() or ["parent"])[0]
        return name, _clean(name, "".join(lines[1:]))
    if "," in lines[0]:
        name, _, sequence = lines[0].partition(",")
        return name.strip() or "parent", _clean(name, sequence)
    return "parent", _clean("parent", "".join(lines))


def _clean(name: str, sequence: str) -> str:
    seq = "".join(sequence.split()).upper()
    if not seq:
        raise ValueError(f"{name}: empty sequence")
    bad = [(i, aa) for i, aa in enumerate(seq, start=1) if aa not in AMINO_ACIDS]
    if bad:
        shown = ", ".join(f"{aa} at {i}" for i, aa in bad[:5])
        raise ValueError(f"{name}: not a standard amino acid: {shown}")
    return seq


def read_designs(path: str) -> list[str]:
    """The variant ids of a design table, in `design` order when the table has that column, else in file order."""
    import pandas as pd

    table = Path(path)
    frame = pd.read_parquet(table) if table.suffix.lower() == ".parquet" else pd.read_csv(table)
    if "variant_id" not in frame.columns:
        raise ValueError(f"{path} has no variant_id column; it has {', '.join(map(str, frame.columns))}")
    if "design" in frame.columns:
        frame = frame.sort_values("design", kind="stable")
    return list(dict.fromkeys(str(v) for v in frame["variant_id"]))


def build(
    layout: scfv.ScfvLayout,
    ids: Sequence[str],
    *,
    min_substitutions: int = 5,
    max_substitutions: int | None = None,
    prefix: str = "design",
) -> list[Design]:
    """The designs as scFv sequences, or one ValueError that lists every problem found."""
    problems: list[str] = []
    designs: list[Design] = []
    first_with: dict[str, str] = {}
    for number, variant_id in enumerate(ids, start=1):
        name = f"{prefix}_{number:02d}"
        try:
            mutations = parse_id(variant_id)
            sequence = scfv.apply_mutations(layout, mutations)
        except ValueError as exc:
            problems.append(f"{name} ({variant_id}): {exc}")
            continue
        substitutions = scfv.substitutions(layout.scfv, sequence)
        count = len(substitutions)
        if count < min_substitutions:
            problems.append(f"{name} ({variant_id}): {count} substitutions, the rule is at least {min_substitutions}")
        if max_substitutions is not None and count > max_substitutions:
            problems.append(f"{name} ({variant_id}): {count} substitutions, more than the {max_substitutions} allowed")
        if sorted(scfv.crystal_mutations(layout, sequence)) != sorted(mutations):
            problems.append(f"{name} ({variant_id}): does not map back to the same crystal mutations; this is a bug in the mapping")
        if STREP_MOTIF in sequence:
            problems.append(f"{name} ({variant_id}): the sequence carries {STREP_MOTIF}; the lab appends the Twin-Strep tag itself")
        if sequence in first_with:
            problems.append(f"{name} ({variant_id}): the same sequence as {first_with[sequence]}")
        first_with.setdefault(sequence, name)
        designs.append(Design(name, make_id(mutations), sequence, substitutions))
    if STREP_MOTIF in layout.scfv:
        problems.append(f"the parent carries {STREP_MOTIF}; the lab appends the Twin-Strep tag itself")
    if not ids:
        problems.append("no designs were given")
    if problems:
        raise ValueError(f"{len(problems)} problem(s):\n  " + "\n  ".join(problems))
    return designs


def write_submission(path: Path, designs: Sequence[Design], header: bool = False) -> Path:
    with path.open("w", newline="") as handle:
        writer = csv.writer(handle)
        if header:
            writer.writerow(["name", "sequence"])
        for design in designs:
            writer.writerow([design.name, design.sequence])
    return path


def write_fasta(path: Path, parent: tuple[str, str], designs: Sequence[Design]) -> Path:
    records = [parent] + [(d.name, d.sequence) for d in designs]
    path.write_text("".join(f">{name}\n{sequence}\n" for name, sequence in records))
    return path


def write_report(path: Path, layout: scfv.ScfvLayout, designs: Sequence[Design]) -> Path:
    chains = sorted(layout.offsets)
    with path.open("w", newline="") as handle:
        writer = csv.writer(handle)
        writer.writerow(["name", "variant_id", "n_substitutions", *[f"n_in_chain_{c}" for c in chains], "scfv_substitutions"])
        for design in designs:
            per_chain = [scfv.chain_position(layout, position)[0] for position, _wt, _mut in design.substitutions]
            writer.writerow(
                [
                    design.name,
                    design.variant_id,
                    len(design.substitutions),
                    *[per_chain.count(c) for c in chains],
                    ";".join(f"{wt}{position}{mut}" for position, wt, mut in design.substitutions),
                ]
            )
    return path


def main(argv: Sequence[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    parser.add_argument("--parent", required=True, help="FASTA file, name,sequence file, or a bare sequence")
    source = parser.add_mutually_exclusive_group(required=True)
    source.add_argument("--designs", help="the design table pick_designs wrote (parquet or csv with a variant_id column)")
    source.add_argument("--ids", nargs="+", help="variant ids such as B:S97A+A:Q27E, one per design")
    parser.add_argument("--structure", help="the crystal, to check the parent against its variable domains")
    parser.add_argument("--heavy", default="B", help="the crystal's heavy chain id (default B, as in 1N8Z)")
    parser.add_argument("--light", default="A", help="the crystal's light chain id (default A, as in 1N8Z)")
    parser.add_argument("--light-first", action="store_true", help="the scFv is VL-linker-VH (default VH first)")
    parser.add_argument("--min-substitutions", type=int, default=5)
    parser.add_argument("--max-substitutions", type=int)
    parser.add_argument("--name-prefix", default="design")
    parser.add_argument("--header", action="store_true", help="write a name,sequence header line in submission.csv")
    parser.add_argument("--out", required=True, help="directory for the three output files")
    args = parser.parse_args(argv)

    try:
        parent_name, parent = read_parent(args.parent)
        first, second = (args.light, args.heavy) if args.light_first else (args.heavy, args.light)
        layout = scfv.layout_from_parent(parent, first=first, second=second)
        if args.structure:
            scfv.verify_parent(layout, args.structure)
            print(f"parent {parent_name} ({len(parent)} residues) matches the variable domains of {args.structure}")
        else:
            print(f"parent {parent_name} ({len(parent)} residues): NOT checked against a crystal, pass --structure to check it")
        ids = read_designs(args.designs) if args.designs else list(args.ids)
        designs = build(layout, ids, min_substitutions=args.min_substitutions, max_substitutions=args.max_substitutions, prefix=args.name_prefix)
    except ValueError as exc:
        print(f"not written: {exc}", file=sys.stderr)
        return 1

    out = Path(args.out)
    out.mkdir(parents=True, exist_ok=True)
    files = [
        write_submission(out / "submission.csv", designs, header=args.header),
        write_fasta(out / "scfv_for_folding.fasta", (parent_name, parent), designs),
        write_report(out / "designs_report.csv", layout, designs),
    ]
    counts = sorted(len(d.substitutions) for d in designs)
    print(f"{len(designs)} designs, {counts[0]} to {counts[-1]} substitutions each, all in the variable domains, linker untouched")
    for file in files:
        print(f"  wrote {file}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
