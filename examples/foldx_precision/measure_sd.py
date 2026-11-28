"""Per-run FoldX ddG for a sample of the 432, keeping the three shifts the
production wrapper averages away.

Imports the adapter's own `run_foldx` and `parse_interaction_energy` so the
numbers are produced by exactly the calls that made foldx.csv. Writes its own
CSV; it does not touch foldx.csv or the adapter's schema.

Sample: ranks 1-30 of ranking_qualified_v1_432.csv, where the top-15 boundary
question lives, plus 30 evenly spaced over ranks 31-432 for a median SD that
is not just a property of the best designs.
"""
import csv
import os
import shutil
import statistics
import sys
import tempfile
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path

REPO = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(REPO))

from design_pipeline.adapters._foldx_score import (  # noqa: E402
    mutation_line,
    parse_interaction_energy,
    read_candidates,
    run_foldx,
)

RANKING = REPO / "runs/rerank/ranking_qualified_v1_432.csv"
MUTATIONS = REPO / "runs/track_a_foldx/candidates_all432.csv"
OUT = REPO / "runs/track_a_foldx/foldx_sd_sample.csv"

RUNS = int(os.environ.get("FOLDX_RUNS", "3"))
WORKERS = int(os.environ.get("FOLDX_WORKERS", "8"))
GROUPS = os.environ.get("FOLDX_GROUPS", "AB,C")
BINARY = Path(os.environ["FOLDX_BINARY"]).resolve()
PDB = Path(os.environ["DESIGN_PIPELINE_FOLDX_PDB"]).resolve()


def shifts_for(job):
    """The per-run binding ddGs of one candidate. Mirrors the adapter's
    score_candidate, except that it returns the list instead of its mean."""
    variant_id, mutations, scratch = job
    work = Path(tempfile.mkdtemp(dir=scratch, prefix="v_")).resolve()
    try:
        (work / "in").mkdir()
        (work / "bm").mkdir()
        shutil.copy(PDB, work / "in" / PDB.name)
        (work / "individual_list.txt").write_text(mutation_line(mutations) + "\n")
        dif = work / "bm" / f"Dif_{PDB.stem}.fxout"
        run_foldx(
            BINARY,
            ["--command=BuildModel", f"--pdb={PDB.name}", f"--pdb-dir={work / 'in'}",
             f"--mutant-file={work / 'individual_list.txt'}", f"--numberOfRuns={RUNS}",
             f"--output-dir={work / 'bm'}"],
            dif,
        )
        models = sorted(p for p in (work / "bm").glob("*.pdb") if not p.name.startswith("WT_"))
        shifts = []
        for k, model in enumerate(models):
            energies = []
            for name in (model.name, f"WT_{model.name}"):
                out = work / f"ac_{k}_{name[:-4]}"
                out.mkdir()
                result = out / f"Interaction_{name[:-4]}_AC.fxout"
                run_foldx(
                    BINARY,
                    ["--command=AnalyseComplex", f"--pdb={name}", f"--pdb-dir={work / 'bm'}",
                     f"--analyseComplexChains={GROUPS}", f"--output-dir={out}"],
                    result,
                )
                energies.append(parse_interaction_energy(result))
            shifts.append(energies[0] - energies[1])
        return {"variant_id": variant_id, "shifts": shifts}
    except Exception as exc:
        return {"variant_id": variant_id, "error": str(exc).splitlines()[0][:160]}
    finally:
        shutil.rmtree(work, ignore_errors=True)


def main() -> int:
    with RANKING.open(newline="") as handle:
        ranked = list(csv.DictReader(handle))
    candidates = read_candidates(MUTATIONS)

    top = list(range(1, 31))
    # 30 evenly spaced over the rest, deterministic - no RNG, so this sample is
    # reproducible from the ranking file alone.
    rest_lo, rest_hi, n = 31, len(ranked), 30
    spread = sorted({rest_lo + round(i * (rest_hi - rest_lo) / (n - 1)) for i in range(n)})
    picks = top + spread
    rows = [ranked[r - 1] for r in picks]
    print(f"sample: {len(rows)} candidates - ranks 1-30 plus {len(spread)} spread over {rest_lo}-{rest_hi}")
    print(f"spread ranks: {spread}")

    scratch = Path(tempfile.mkdtemp(dir=OUT.parent, prefix="sdwork_"))
    jobs = [(r["variant_id"], candidates[r["variant_id"]], scratch) for r in rows]
    by_id = {r["variant_id"]: r for r in rows}

    fields = ["rank", "name", "source", "variant_id", "n_subs", "ddg_mean", "ddg_sd",
              "ddg_min", "ddg_max", "ddg_range", "n_runs", "ddg_runs", "ddg_in_foldx_csv"]
    results, failed, done = [], [], 0
    with OUT.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields)
        writer.writeheader()
        with ThreadPoolExecutor(max_workers=WORKERS) as pool:
            futures = {pool.submit(shifts_for, job): job[0] for job in jobs}
            for finished in as_completed(futures):
                got = finished.result()
                meta = by_id[got["variant_id"]]
                if "shifts" in got and len(got["shifts"]) >= 2:
                    s = got["shifts"]
                    row = {
                        "rank": meta["rank"], "name": meta["name"], "source": meta["source"],
                        "variant_id": got["variant_id"], "n_subs": meta["n_subs"],
                        "ddg_mean": round(statistics.mean(s), 6),
                        "ddg_sd": round(statistics.stdev(s), 6),
                        "ddg_min": round(min(s), 6), "ddg_max": round(max(s), 6),
                        "ddg_range": round(max(s) - min(s), 6),
                        "n_runs": len(s),
                        "ddg_runs": ";".join(f"{v:.6f}" for v in s),
                        "ddg_in_foldx_csv": meta["foldx_ddg"],
                    }
                    writer.writerow(row)
                    handle.flush()
                    results.append(row)
                else:
                    failed.append(got)
                done += 1
                print(f"sd: {done}/{len(jobs)} ({len(results)} measured, {len(failed)} failed)", flush=True)
    shutil.rmtree(scratch, ignore_errors=True)

    print(f"\nwrote {OUT}: {len(results)} rows, {len(failed)} failed")
    for f in failed:
        print(f"  FAILED {f['variant_id']}: {f.get('error')}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
