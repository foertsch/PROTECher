"""FoldX binding ddG for each candidate in a variants table: the script `foldx_process` in `_scorers.py` runs.

Standard library only, so it runs under any Python 3, and it imports nothing from this package. Per candidate it
runs FoldX BuildModel (the mutant, with every mutation of the candidate in one model, and a wild-type copy built the
same way, `FOLDX_RUNS` times) and then AnalyseComplex on both with the antibody and the antigen as the two groups:

    ddg = interaction energy (mutant) - interaction energy (wild type), kcal/mol, positive means weaker binding

The same calls as `examples/foldx_alanine_scan.py`, so a number here is comparable with that scan's `ddg_bind`.
FoldX is licensed per user and is not part of this repository: its binary, and a RepairPDB output to run on, come
from the environment.

    FOLDX_BINARY               the FoldX executable (required); it runs with its own folder as the working directory
    DESIGN_PIPELINE_FOLDX_PDB  the structure to build on, replacing --pdb. Use the trimmed, repaired Fv plus antigen
                               domain the alanine scan ran on: FoldX time grows with the residue count. Residue numbers
                               must be those of --pdb.
    FOLDX_GROUPS               the two complex groups as chain letters, default AB,C
    FOLDX_RUNS                 BuildModel runs per candidate, default 3
    FOLDX_WORKERS              FoldX processes in parallel, default 4

Writes foldx.csv (variant_id, ddg, ddg_fold, n_runs). A candidate FoldX fails on is left out and named on stderr, so
the scorer reads it as no score and never as zero.
"""

from __future__ import annotations

import argparse
import csv
import os
import shutil
import statistics
import subprocess
import sys
import tempfile
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path
from typing import Sequence

Mutation = tuple[str, int, str, str]    # chain, position, wild type, mutant


def read_candidates(path: Path) -> dict[str, list[Mutation]]:
    """variant_id -> its mutations, from the one-row-per-mutation table the pipeline writes."""
    grouped: dict[str, list[Mutation]] = {}
    with path.open(newline="", encoding="utf-8") as handle:
        for row in csv.DictReader(handle):
            grouped.setdefault(row["variant_id"], []).append((row["chain"], int(row["position"]), row["wt"], row["mut"]))
    return grouped


def mutation_line(mutations: Sequence[Mutation]) -> str:
    """FoldX's notation, one line for the whole candidate: wild type, chain, number, mutant, comma-separated, then a semicolon."""
    return ",".join(f"{wt}{chain}{position}{mut}" for chain, position, wt, mut in mutations) + ";"


def parse_interaction_energy(path: Path) -> float:
    """The 'Interaction Energy' column of the last data line of an AnalyseComplex file."""
    lines = [line for line in path.read_text().splitlines() if line.strip()]
    header = next(line for line in lines if line.startswith("Pdb\t")).split("\t")
    return float(lines[-1].split("\t")[header.index("Interaction Energy")])


def parse_total_energies(path: Path) -> list[float]:
    """The 'total energy' of each model row in a BuildModel Dif file."""
    lines = [line for line in path.read_text().splitlines() if line.strip()]
    start = next(i for i, line in enumerate(lines) if line.startswith("Pdb\t"))
    return [float(line.split("\t")[1]) for line in lines[start + 1:]]


def run_foldx(foldx: Path, args: Sequence[str], expect: Path) -> None:
    """FoldX exits 0 even when it fails, so success is the output file being there."""
    done = subprocess.run([str(foldx), *args], cwd=foldx.parent, capture_output=True, text=True, check=False)
    if done.returncode != 0 or not expect.is_file():
        tail = "\n".join((done.stdout + done.stderr).strip().splitlines()[-4:])
        raise RuntimeError(f"foldx {args[0].removeprefix('--command=')} did not write {expect.name}:\n{tail}")


def score_candidate(job: tuple) -> dict:
    foldx, pdb, variant_id, mutations, runs, groups, scratch = job
    work = Path(tempfile.mkdtemp(dir=scratch, prefix="v_")).resolve()
    try:
        (work / "in").mkdir()
        (work / "bm").mkdir()
        shutil.copy(pdb, work / "in" / pdb.name)
        (work / "individual_list.txt").write_text(mutation_line(mutations) + "\n")
        dif = work / "bm" / f"Dif_{pdb.stem}.fxout"
        run_foldx(
            foldx,
            ["--command=BuildModel", f"--pdb={pdb.name}", f"--pdb-dir={work / 'in'}",
             f"--mutant-file={work / 'individual_list.txt'}", f"--numberOfRuns={runs}", f"--output-dir={work / 'bm'}"],
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
        return {"variant_id": variant_id, "ddg": statistics.mean(shifts), "ddg_fold": statistics.mean(parse_total_energies(dif)), "n_runs": len(shifts)}
    except Exception as exc:    # one failing candidate must not lose the others
        return {"variant_id": variant_id, "error": str(exc).splitlines()[0][:160]}
    finally:
        shutil.rmtree(work, ignore_errors=True)


def main(argv: Sequence[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    parser.add_argument("--pdb", type=Path, required=True)
    parser.add_argument("--out_dir", type=Path, required=True)
    parser.add_argument("--mutations_csv", type=Path, required=True, help="variant_id,chain,position,wt,mut, one row per mutation")
    args = parser.parse_args(argv)

    binary = os.environ.get("FOLDX_BINARY", "")
    if not binary or not Path(binary).is_file():
        parser.error("set FOLDX_BINARY to the FoldX executable")
    pdb = Path(os.environ.get("DESIGN_PIPELINE_FOLDX_PDB") or args.pdb)
    if not pdb.is_file():
        parser.error(f"structure {pdb} does not exist")
    groups = os.environ.get("FOLDX_GROUPS", "AB,C")
    runs, workers = int(os.environ.get("FOLDX_RUNS", "3")), int(os.environ.get("FOLDX_WORKERS", "4"))
    if runs < 1 or workers < 1:
        parser.error("FOLDX_RUNS and FOLDX_WORKERS must be at least 1")

    candidates = read_candidates(args.mutations_csv)
    args.out_dir.mkdir(parents=True, exist_ok=True)
    scratch = Path(tempfile.mkdtemp(dir=args.out_dir, prefix="work_"))
    jobs = [(Path(binary).resolve(), pdb.resolve(), vid, muts, runs, groups, scratch) for vid, muts in candidates.items()]
    print(f"foldx: {len(jobs)} candidates on {pdb.name}, {runs} run(s) each, {workers} workers", flush=True)
    with ThreadPoolExecutor(max_workers=workers) as pool:
        rows = list(pool.map(score_candidate, jobs))
    shutil.rmtree(scratch, ignore_errors=True)

    done = [r for r in rows if "ddg" in r]
    failed = [r for r in rows if "error" in r]
    with (args.out_dir / "foldx.csv").open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=["variant_id", "ddg", "ddg_fold", "n_runs"])
        writer.writeheader()
        writer.writerows(done)
    print(f"foldx: {len(done)} scored, {len(failed)} failed", flush=True)
    for row in failed[:10]:
        print(f"  failed {row['variant_id']}: {row['error']}", file=sys.stderr)
    return 0 if done else 1


if __name__ == "__main__":
    sys.exit(main())
