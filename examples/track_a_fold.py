"""Fold announced scFv variants with their antigen: Boltz-2, several seeds.

The Boltz-2 path of docs/track_a_filtering.md. Stage 2 folds the announced scFv
with its antigen to get a structure for the position map. Stage 6 folds the
parent and the finalists as the structure veto, and the veto needs several
seeds to measure seed noise. Both are this script with different input.

The antigen is trimmed first. A whole extracellular domain is several hundred
residues, and the 144 s we measured was for a 242-residue scFv with a
128-residue HER2 fragment, so folding about 100 complexes is only affordable
with a fragment. Trim with `--antigen-region` (1-based, inclusive, numbered as
in the sequence you pass) or with `--epitope` plus `--pad`. With neither, the
whole antigen is folded and a long one only earns a warning.

    python examples/track_a_fold.py --scfv scfv.fasta --antigen antigen.fasta \\
        --antigen-region 480-607 --seeds 5 --out runs/track_a

    python examples/track_a_fold.py --scfv scfv.fasta --antigen antigen.fasta --dry-run
    python examples/track_a_fold.py --scfv scfv.fasta --antigen antigen.fasta --dummy

`--scfv` and `--antigen` take a FASTA file or a bare sequence. The scFv is
chain A and the antigen chain B, the order `pae_interaction` assumes. Give the
parent first and the variants after it, one FASTA record each.

Folding needs the cluster: Boltz-2 wants a CUDA GPU. `--dry-run` writes the
inputs and stops, `--dummy` runs the fixture model offline.

If the cluster is down or full, the folds run on Modal instead: `--dry-run` leaves
`<out>/inputs/*.fa` in the form `cloud/modal_fold.py::fold_dir` reads, trimmed antigen
included, and `fold_dir --seeds 5` folds them in parallel. See cloud/README.md.
"""

from __future__ import annotations

import argparse
import csv
import re
import statistics
import sys
from pathlib import Path
from typing import Sequence

AMINO_ACIDS = frozenset("ACDEFGHIKLMNPQRSTVWY")
STREP_MOTIF = "WSHPQFEK"
"""Adaptyv appends the Twin-Strep tag itself. A sequence that already carries it
would fold with a tag the lab then adds a second time."""
LONG_ANTIGEN = 300
"""Folded antigen length above which fold time stops matching the 144 s we
measured. A warning, not a limit."""


def read_sequences(source: str) -> list[tuple[str, str]]:
    """(name, sequence) pairs from a FASTA file, or one bare sequence."""
    if not _is_file(source):
        return [("seq", clean("seq", source))]
    path = Path(source)
    records: list[tuple[str, str]] = []
    name: str | None = None
    chunks: list[str] = []
    for line in path.read_text().splitlines():
        if line.startswith(">"):
            if name is not None:
                records.append((name, "".join(chunks)))
            name = (line[1:].split() or [f"seq{len(records) + 1}"])[0]
            chunks = []
        elif line.strip():
            if name is None:
                raise ValueError(f"{path}: sequence before the first '>' header")
            chunks.append(line.strip())
    if name is not None:
        records.append((name, "".join(chunks)))
    if not records:
        raise ValueError(f"{path} holds no FASTA records")
    return [(name, clean(name, seq)) for name, seq in records]


def _is_file(text: str) -> bool:
    try:
        return Path(text).is_file()
    except OSError:  # a bare sequence is too long to be a file name
        return False


def clean(name: str, sequence: str) -> str:
    """Uppercase, whitespace removed, 20 standard amino acids only.

    Boltz would fail on anything else deep inside the run, so this says which
    residue and where while the sequence is still on the screen.
    """
    seq = "".join(sequence.split()).upper()
    if not seq:
        raise ValueError(f"{name}: empty sequence")
    bad = [(i, aa) for i, aa in enumerate(seq, start=1) if aa not in AMINO_ACIDS]
    if bad:
        shown = ", ".join(f"{aa} at {i}" for i, aa in bad[:5])
        more = " and more" if len(bad) > 5 else ""
        raise ValueError(f"{name}: not a standard amino acid: {shown}{more}")
    return seq


def parse_region(text: str) -> tuple[int, int]:
    match = re.fullmatch(r"\s*(\d+)\s*-\s*(\d+)\s*", text)
    if match is None:
        raise ValueError(f"{text!r} is not START-END, for example 480-607")
    start, end = int(match[1]), int(match[2])
    if start < 1 or end < start:
        raise ValueError(f"{text!r}: residues count from 1 and END must not precede START")
    return start, end


def parse_residues(text: str) -> list[int]:
    """'557,558,560-565' as a list of residue numbers."""
    residues: list[int] = []
    for part in text.split(","):
        part = part.strip()
        if not part:
            continue
        if "-" in part:
            start, end = parse_region(part)
            residues.extend(range(start, end + 1))
        else:
            residues.append(int(part))
    if not residues or min(residues) < 1:
        raise ValueError(f"{text!r} names no residues, or one below 1")
    return residues


def trim_antigen(
    sequence: str,
    region: tuple[int, int] | None = None,
    epitope: Sequence[int] | None = None,
    pad: int = 20,
) -> tuple[str, int, int]:
    """The antigen fragment to fold, with its 1-based start and end.

    A region is taken as given. Epitope residues give the span they cover plus
    `pad` on each side, clipped to the sequence. Neither means the whole antigen.
    """
    if region is not None and epitope:
        raise ValueError("give an antigen region or epitope residues, not both")
    length = len(sequence)
    if region is None and not epitope:
        return sequence, 1, length
    if region is not None:
        start, end = region
    else:
        assert epitope is not None
        if max(epitope) > length:
            raise ValueError(
                f"epitope residue {max(epitope)} is past the antigen's {length} residues"
            )
        start, end = max(1, min(epitope) - pad), min(length, max(epitope) + pad)
    if end > length:
        raise ValueError(f"antigen region {start}-{end} runs past the antigen's {length} residues")
    return sequence[start - 1 : end], start, end


def complex_fasta(name: str, scfv: str, antigen: str) -> str:
    """One record, chains joined by ':', the form design_sequences writes."""
    return f">{name}\n{scfv}:{antigen}\n"


def warnings_for(scfvs: Sequence[tuple[str, str]], antigen_length: int) -> list[str]:
    found = [
        f"{name} contains {STREP_MOTIF}, part of the Twin-Strep tag Adaptyv adds itself"
        for name, seq in scfvs
        if STREP_MOTIF in seq
    ]
    if antigen_length > LONG_ANTIGEN:
        found.append(
            f"the antigen to fold has {antigen_length} residues; fold time grows with "
            f"length, so trim it with --antigen-region or --epitope"
        )
    return found


def write_inputs(
    out: Path, scfvs: Sequence[tuple[str, str]], antigen: str, antigen_span: tuple[int, int]
) -> list[tuple[str, Path]]:
    """One two-chain FASTA per scFv in `inputs/`, and the trimmed antigen for the record beside it.

    The record stays out of `inputs/` so that folder holds only complexes, which is what a batch fold reads
    (`cloud/modal_fold.py::fold_dir` takes every `*.fa` in it as a design).
    """
    inputs = out / "inputs"
    inputs.mkdir(parents=True, exist_ok=True)
    start, end = antigen_span
    (out / "antigen_trimmed.fa").write_text(f">antigen_{start}-{end}\n{antigen}\n")
    written: list[tuple[str, Path]] = []
    seen: set[str] = set()
    for name, scfv in scfvs:
        safe = re.sub(r"[^A-Za-z0-9_.-]", "_", name)
        if safe in seen:
            raise ValueError(f"two scFv records are named {safe!r}; give each a unique name")
        seen.add(safe)
        path = inputs / f"{safe}.fa"
        path.write_text(complex_fasta(safe, scfv, antigen))
        written.append((safe, path))
    return written


def run_folds(
    inputs: Sequence[tuple[str, Path]], out: Path, model: str, seeds: Sequence[int]
) -> list:
    """Fold every input once per seed with the named adapter, and tag each prediction with its seed.

    The adapter is called directly through `get_adapter`, which the pipeline exports on both its
    older and its current form, instead of through a pipeline step: a step records only the top
    prediction's scores, and the veto needs the ipTM of every seed. There is no cache, so a rerun
    folds again.
    """
    from dataclasses import replace

    from design_pipeline import Design, get_adapter

    adapter = get_adapter(model)
    designs = [
        Design(id=name, path=path, parent_id=None, producer="track_a_fold")
        for name, path in inputs
    ]
    predictions: list = []
    for seed in seeds:
        run_dir = out / "predictions" / f"seed{seed}"
        run_dir.mkdir(parents=True, exist_ok=True)
        for prediction in adapter.run(designs, run_dir, seed=seed):
            predictions.append(replace(prediction, meta={**prediction.meta, "seed": seed}))
    return predictions


def summarise(predictions: Sequence) -> list[dict]:
    """One row per prediction: the complex it came from, the seed and the scores."""
    return [
        {
            "complex": p.parent_id,
            "seed": p.meta.get("seed"),
            "iptm": p.scores.get("iptm"),
            "ptm": p.scores.get("ptm"),
            "plddt": p.scores.get("plddt"),
            "pdb": str(p.path),
        }
        for p in predictions
    ]


def write_summary(rows: Sequence[dict], path: Path) -> None:
    fields = ["complex", "seed", "iptm", "ptm", "plddt", "pdb"]
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields)
        writer.writeheader()
        writer.writerows(rows)


def seed_spread(rows: Sequence[dict]) -> list[tuple[str, int, float, float | None]]:
    """(complex, seeds, mean ipTM, SD of ipTM) per complex. The SD is the seed
    noise the veto compares against, and needs at least two seeds."""
    by_complex: dict[str, list[float]] = {}
    for row in rows:
        if row["iptm"] is not None:
            by_complex.setdefault(row["complex"], []).append(row["iptm"])
    return [
        (
            name,
            len(values),
            statistics.mean(values),
            statistics.stdev(values) if len(values) > 1 else None,
        )
        for name, values in by_complex.items()
    ]


def main(argv: Sequence[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    parser.add_argument("--scfv", required=True, help="FASTA file or bare sequence; parent first")
    parser.add_argument("--antigen", required=True, help="FASTA file or bare sequence")
    parser.add_argument("--antigen-region", help="START-END to keep, 1-based and inclusive")
    parser.add_argument("--epitope", help="epitope residues to keep, e.g. 557,558,560-565")
    parser.add_argument("--pad", type=int, default=20, help="residues kept around the epitope")
    parser.add_argument("--seeds", type=int, default=5, help="seeds 0..N-1 (default 5)")
    parser.add_argument("--out", type=Path, default=Path("runs/track_a"))
    parser.add_argument("--model", default="boltz")
    parser.add_argument("--dummy", action="store_true", help="fixture model, runs offline")
    parser.add_argument("--dry-run", action="store_true", help="write the inputs and stop")
    args = parser.parse_args(argv)
    if args.seeds < 1:
        parser.error("--seeds must be at least 1")

    try:
        scfvs = read_sequences(args.scfv)
        antigens = read_sequences(args.antigen)
        if len(antigens) != 1:
            raise ValueError(f"--antigen holds {len(antigens)} records; give exactly one")
        region = parse_region(args.antigen_region) if args.antigen_region else None
        epitope = parse_residues(args.epitope) if args.epitope else None
        antigen, start, end = trim_antigen(antigens[0][1], region, epitope, args.pad)
        extra = warnings_for(scfvs, len(antigen))
        inputs = write_inputs(args.out, scfvs, antigen, (start, end))
    except ValueError as exc:
        parser.error(str(exc))

    for message in extra:
        print(f"warning: {message}", file=sys.stderr)
    model = "dummy_boltz" if args.dummy else args.model
    seeds = list(range(args.seeds))
    print(
        f"{len(inputs)} complex(es), antigen residues {start}-{end} of "
        f"{len(antigens[0][1])} ({len(antigen)} kept), model {model}, seeds {seeds}"
    )
    if args.dry_run:
        print(f"dry run: inputs written to {args.out / 'inputs'}, nothing folded")
        return 0

    rows = summarise(run_folds(inputs, args.out, model, seeds))
    write_summary(rows, args.out / "summary.csv")
    for name, n, mean, sd in seed_spread(rows):
        spread = f"{sd:.3f}" if sd is not None else "n/a (one seed)"
        print(f"{name}: ipTM mean {mean:.3f}, SD {spread} over {n} seed(s)")
    print(f"wrote {args.out / 'summary.csv'}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
