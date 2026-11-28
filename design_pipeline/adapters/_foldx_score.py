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

Writes foldx.csv (variant_id, ddg, ddg_sd, ddg_folding, n_runs). A candidate FoldX fails on is left out and named on
stderr, so the scorer reads it as no score and never as zero.

    ddg           mean binding ddG over the FOLDX_RUNS models, kcal/mol, positive means weaker binding
    ddg_sd        sample SD of that mean's runs - FoldX's own precision on this candidate. Empty at FOLDX_RUNS=1,
                  where there is no spread to report. Without it a ranking cannot say whether two candidates
                  differ or whether BuildModel simply landed in two places.
    ddg_folding   mean TOTAL energy change BuildModel reports: the folding free energy, a different quantity from
                  ddg and not an aggregate of it. Called ddg_fold until 2026-11-24; a table under the old header is
                  rewritten in place on the next run.
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
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path
from typing import Sequence

Mutation = tuple[str, int, str, str]    # chain, position, wild type, mutant

FIELDS = ["variant_id", "ddg", "ddg_sd", "ddg_folding", "n_runs"]
RENAMED = {"ddg_fold": "ddg_folding"}   # header changes a table written by an older version has to survive


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
        # The spread of `shifts` is the only measure of FoldX's precision on this candidate, and it was computed and
        # thrown away here: two candidates 0.001 kcal/mol apart used to be indistinguishable from two 1.0 apart. It
        # needs two runs; at FOLDX_RUNS=1 the column is empty rather than 0, which would read as perfect agreement.
        return {"variant_id": variant_id, "ddg": statistics.mean(shifts),
                "ddg_sd": statistics.stdev(shifts) if len(shifts) >= 2 else "",
                "ddg_folding": statistics.mean(parse_total_energies(dif)), "n_runs": len(shifts)}
    except Exception as exc:    # one failing candidate must not lose the others
        return {"variant_id": variant_id, "error": str(exc).splitlines()[0][:160]}
    finally:
        shutil.rmtree(work, ignore_errors=True)


def read_scored(out_csv: Path) -> tuple[dict[str, dict], bool]:
    """The candidates already in foldx.csv, and whether the file had to be rewritten to hold them.

    A table written before ddg_sd existed, or under ddg_fold's old name, carries a header that an append of today's
    FIELDS would write five columns under - rows that no longer line up with their own header, and nothing would say
    so until whoever reads the file next gets the wrong column. So the file is brought up to FIELDS first, keeping
    every row. An older run's per-run spread cannot be recovered, so ddg_sd stays EMPTY for those rows rather than 0,
    which would claim the three runs agreed perfectly.
    """
    if not out_csv.is_file():
        return {}, False
    with out_csv.open(newline="") as handle:
        reader = csv.DictReader(handle)
        present = list(reader.fieldnames or ())
        scored = {row["variant_id"]: row for row in reader if row.get("variant_id")}
    if not present or present == FIELDS:
        return scored, False
    with out_csv.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=FIELDS)
        writer.writeheader()
        for row in scored.values():
            under_new_names = {RENAMED.get(key, key): value for key, value in row.items()}
            writer.writerow({key: under_new_names.get(key, "") for key in FIELDS})
    return scored, True


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
    # "AB,C" is 1N8Z's own layout (A light, B heavy, C HER2) and nothing below
    # would notice if the structure disagreed: AnalyseComplex would report the
    # interaction energy of whichever chains those letters do name, or of
    # nothing, and the ddG would be a real number for the wrong interface. So
    # check the letters against the file rather than trusting the default.
    in_pdb = {
        line[21]
        for line in pdb.read_text(encoding="utf-8", errors="replace").splitlines()
        if line.startswith(("ATOM", "HETATM")) and len(line) > 21
    }
    named = {c for group in groups.split(",") for c in group.strip() if c.strip()}
    if not named or len(groups.split(",")) != 2:
        parser.error(
            f"FOLDX_GROUPS is {groups!r}; it has to be two comma-separated "
            f"groups of chain letters, as in AB,C"
        )
    unknown = sorted(named - in_pdb)
    if unknown:
        parser.error(
            f"FOLDX_GROUPS={groups} names chain(s) {', '.join(unknown)}, which "
            f"{pdb} does not have (it has {', '.join(sorted(in_pdb))}). Every "
            f"ddG would be the wrong interface."
        )
    runs, workers = int(os.environ.get("FOLDX_RUNS", "3")), int(os.environ.get("FOLDX_WORKERS", "4"))
    if runs < 1 or workers < 1:
        parser.error("FOLDX_RUNS and FOLDX_WORKERS must be at least 1")

    candidates = read_candidates(args.mutations_csv)
    args.out_dir.mkdir(parents=True, exist_ok=True)
    scratch = Path(tempfile.mkdtemp(dir=args.out_dir, prefix="work_"))
    jobs = [(Path(binary).resolve(), pdb.resolve(), vid, muts, runs, groups, scratch) for vid, muts in candidates.items()]
    # Resumable, and incremental on purpose. This used to collect every row with
    # pool.map and write foldx.csv once at the end, so an interrupted run - a
    # stopped VM, a dropped connection, Ctrl-C - threw away every candidate it
    # had already paid for. At 7 FoldX invocations per candidate and several
    # minutes each, that is hours. Rows are now appended as they finish and
    # candidates already in foldx.csv are skipped, so a re-run continues instead
    # of restarting.
    fields = FIELDS
    out_csv = args.out_dir / "foldx.csv"
    scored_already, rewritten = read_scored(out_csv)
    if rewritten:
        print(
            f"foldx: rewrote {out_csv.name} under the current header "
            f"({', '.join(fields)}); {len(scored_already)} row(s) kept, ddg_sd empty for them",
            flush=True,
        )
    pending = [job for job in jobs if job[2] not in scored_already]
    print(
        f"foldx: {len(pending)} candidates on {pdb.name}, {runs} run(s) each, "
        f"{workers} workers"
        + (f"; {len(scored_already)} already scored, skipping" if scored_already else ""),
        flush=True,
    )

    failed: list[dict] = []
    completed = 0
    handle = out_csv.open("a" if scored_already else "w", newline="")
    try:
        writer = csv.DictWriter(handle, fieldnames=fields)
        if not scored_already:
            writer.writeheader()
            handle.flush()
        with ThreadPoolExecutor(max_workers=workers) as pool:
            futures = {pool.submit(score_candidate, job): job[2] for job in pending}
            for finished in as_completed(futures):
                row = finished.result()
                if "ddg" in row:
                    writer.writerow({key: row[key] for key in fields})
                    handle.flush()  # the point of the exercise: survive a kill
                    completed += 1
                else:
                    failed.append(row)
                # A progress line, because an 800-candidate block was otherwise
                # silent for hours with no way to tell how far it had got.
                print(
                    f"foldx: {completed + len(failed)}/{len(pending)} "
                    f"({completed} scored, {len(failed)} failed)",
                    flush=True,
                )
    finally:
        handle.close()
    shutil.rmtree(scratch, ignore_errors=True)

    print(f"foldx: {completed} scored this run, {len(failed)} failed, {len(scored_already)} reused", flush=True)
    for row in failed[:10]:
        print(f"  failed {row['variant_id']}: {row['error']}", file=sys.stderr)
    return 0 if (completed or scored_already) else 1


if __name__ == "__main__":
    sys.exit(main())
