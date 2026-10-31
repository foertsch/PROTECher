"""Alanine scan with FoldX: the binding energy change of every antibody position set to alanine.

Stage 4 of docs/track_a_filtering.md. Alanine removes the side chain beyond its first carbon, so a
large loss of binding energy marks a position whose side chain carries binding (a hotspot), and the
framework and second-shell arms stay away from those. One run per position on the parent complex.

Per position this runs FoldX BuildModel (the alanine mutant and a wild-type copy built the same way,
`--runs` times) and then AnalyseComplex on both with the antibody and the antigen as the two groups:

    ddg_bind = interaction energy (mutant) - interaction energy (wild type), kcal/mol, positive
               means the alanine weakens binding
    ddg_fold = the total-energy change BuildModel reports, positive means the alanine destabilises
               the fold, which matters for expression and not for binding

Prepare the structure once with `foldx --command=RepairPDB`. The antigen should be trimmed to the
region near the antibody, as FoldX time grows with the number of residues. FoldX needs its
`molecules` folder, so it is run with its own folder as the working directory. It is CPU only.

    python3 examples/foldx_alanine_scan.py --foldx ~/tools/foldx/foldx_20261231 \\
        --pdb repaired.pdb --scan-chains A,B --groups AB,C --runs 3 --workers 12 --out scan_out

`--limit 3` scans the first three positions only, as a smoke test. Needs only the standard library.
FoldX is licensed per user: do not copy the binary or share its output without reading the licence.
"""

from __future__ import annotations

import argparse
import csv
import shutil
import statistics
import subprocess
import sys
import tempfile
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path
from typing import Sequence

THREE_TO_ONE = {
    "ALA": "A", "ARG": "R", "ASN": "N", "ASP": "D", "CYS": "C", "GLN": "Q", "GLU": "E", "GLY": "G",
    "HIS": "H", "ILE": "I", "LEU": "L", "LYS": "K", "MET": "M", "PHE": "F", "PRO": "P", "SER": "S",
    "THR": "T", "TRP": "W", "TYR": "Y", "VAL": "V",
}
FIELDS = ["chain", "resseq", "wt", "ddg_bind", "ddg_bind_sd", "ddg_fold", "n_runs", "status"]


def read_positions(pdb: Path, chains: set[str]) -> list[tuple[str, int, str]]:
    """(chain, residue number, one-letter wild type) from the CA records, in file order."""
    found: list[tuple[str, int, str]] = []
    for line in pdb.read_text().splitlines():
        if not line.startswith("ATOM") or line[12:16].strip() != "CA" or line[21] not in chains:
            continue
        aa = THREE_TO_ONE.get(line[17:20].strip())
        if aa is not None:
            found.append((line[21], int(line[22:26]), aa))
    return found


def mutation_line(wt: str, chain: str, number: int, mutant: str = "A") -> str:
    """FoldX's own notation: wild type, chain, number, mutant, then a semicolon."""
    return f"{wt}{chain}{number}{mutant};"


def parse_interaction_energy(path: Path) -> float:
    """The 'Interaction Energy' column of the last data line of an AnalyseComplex file."""
    lines = [line for line in path.read_text().splitlines() if line.strip()]
    header = next(line for line in lines if line.startswith("Pdb\t")).split("\t")
    return float(lines[-1].split("\t")[header.index("Interaction Energy")])


def parse_total_energies(path: Path) -> list[float]:
    """The 'total energy' of each model row in a BuildModel Dif file."""
    lines = [line for line in path.read_text().splitlines() if line.strip()]
    start = next(i for i, line in enumerate(lines) if line.startswith("Pdb\t"))
    return [float(line.split("\t")[1]) for line in lines[start + 1 :]]


def run_foldx(foldx: Path, args: Sequence[str], expect: Path) -> None:
    """FoldX exits 0 even when it fails, so success is the output file being there."""
    done = subprocess.run(
        [str(foldx), *args], cwd=foldx.parent, capture_output=True, text=True, check=False
    )
    if done.returncode != 0 or not expect.is_file():
        tail = "\n".join((done.stdout + done.stderr).strip().splitlines()[-4:])
        command = args[0].removeprefix("--command=")
        raise RuntimeError(f"foldx {command} did not write {expect.name}:\n{tail}")


def scan_position(job: tuple) -> dict:
    foldx, pdb, chain, number, wt, runs, groups, scratch = job
    row = {"chain": chain, "resseq": number, "wt": wt}
    if wt == "A":
        return {**row, "status": "wt_is_ala"}
    work = Path(tempfile.mkdtemp(dir=scratch, prefix=f"{chain}{number}_")).resolve()
    try:
        (work / "in").mkdir()
        (work / "bm").mkdir()
        shutil.copy(pdb, work / "in" / pdb.name)
        (work / "individual_list.txt").write_text(mutation_line(wt, chain, number) + "\n")
        stem = pdb.stem
        dif = work / "bm" / f"Dif_{stem}.fxout"
        run_foldx(
            foldx,
            ["--command=BuildModel", f"--pdb={pdb.name}", f"--pdb-dir={work / 'in'}",
             f"--mutant-file={work / 'individual_list.txt'}", f"--numberOfRuns={runs}",
             f"--output-dir={work / 'bm'}"],
            dif,
        )
        models = sorted(p for p in (work / "bm").glob("*.pdb") if not p.name.startswith("WT_"))
        shifts: list[float] = []
        for k, model in enumerate(models):
            energies = []
            for name in (model.name, f"WT_{model.name}"):
                out = work / f"ac_{k}_{name[:-4]}"
                out.mkdir()
                result = out / f"Interaction_{name[:-4]}_AC.fxout"
                run_foldx(
                    foldx,
                    ["--command=AnalyseComplex", f"--pdb={name}", f"--pdb-dir={work / 'bm'}",
                     f"--analyseComplexChains={groups}", f"--output-dir={out}"],
                    result,
                )
                energies.append(parse_interaction_energy(result))
            shifts.append(energies[0] - energies[1])
        return {
            **row,
            "ddg_bind": statistics.mean(shifts),
            "ddg_bind_sd": statistics.stdev(shifts) if len(shifts) > 1 else "",
            "ddg_fold": statistics.mean(parse_total_energies(dif)),
            "n_runs": len(shifts),
            "status": "ok",
        }
    except Exception as exc:  # one bad position must not lose the other 200
        return {**row, "status": f"failed: {str(exc).splitlines()[0][:120]}"}
    finally:
        shutil.rmtree(work, ignore_errors=True)


def main(argv: Sequence[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    parser.add_argument("--foldx", required=True, type=Path, help="the FoldX executable")
    parser.add_argument("--pdb", required=True, type=Path, help="a RepairPDB output")
    parser.add_argument("--scan-chains", required=True, help="chains to scan, e.g. A,B")
    parser.add_argument("--groups", required=True, help="the two complex groups, e.g. AB,C")
    parser.add_argument("--runs", type=int, default=3, help="BuildModel runs per position")
    parser.add_argument("--workers", type=int, default=4)
    parser.add_argument("--limit", type=int, help="scan only the first N positions")
    parser.add_argument("--out", type=Path, default=Path("foldx_scan"))
    args = parser.parse_args(argv)
    if args.runs < 1 or args.workers < 1:
        parser.error("--runs and --workers must be at least 1")
    if not args.foldx.is_file() or not args.pdb.is_file():
        parser.error("--foldx and --pdb must be existing files")

    positions = read_positions(args.pdb, set(args.scan_chains.split(",")))
    if not positions:
        parser.error(f"no residues found on chains {args.scan_chains}")
    if args.limit:
        positions = positions[: args.limit]
    args.out.mkdir(parents=True, exist_ok=True)
    scratch = Path(tempfile.mkdtemp(dir=args.out, prefix="work_"))
    jobs = [
        (args.foldx.resolve(), args.pdb.resolve(), chain, number, wt, args.runs, args.groups, scratch)
        for chain, number, wt in positions
    ]
    print(f"scanning {len(jobs)} positions, {args.runs} run(s) each, {args.workers} workers")
    # Threads are enough: each worker spends its time waiting on a FoldX subprocess.
    with ThreadPoolExecutor(max_workers=args.workers) as pool:
        rows = list(pool.map(scan_position, jobs))
    shutil.rmtree(scratch, ignore_errors=True)

    path = args.out / "ddg_alanine_scan.csv"
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=FIELDS)
        writer.writeheader()
        writer.writerows(rows)
    done = [r for r in rows if r["status"] == "ok"]
    failed = [r for r in rows if r["status"].startswith("failed")]
    print(f"{len(done)} scanned, {len(failed)} failed, wrote {path}")
    for r in sorted(done, key=lambda r: -r["ddg_bind"])[:10]:
        print(f"  {r['wt']}{r['chain']}{r['resseq']}A  ddg_bind {r['ddg_bind']:+.2f}  ddg_fold {r['ddg_fold']:+.2f}")
    return 0 if not failed else 1


if __name__ == "__main__":
    sys.exit(main())
