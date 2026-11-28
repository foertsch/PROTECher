"""The top 50 of the 432 ranking at FOLDX_RUNS=20, keeping every per-run ddG.

Same calls as production - the adapter's own run_foldx and parse_interaction_energy - so
these means are comparable with foldx.csv's. Writes its own CSV and touches neither
foldx.csv nor the n=3 sample.

At 20 runs the SE of the mean falls by sqrt(20/3) = 2.58 against the 3-run figure, which is
what decides whether the top of the ranking has any resolvable structure in it.
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
OUT = REPO / "runs/track_a_foldx/foldx_runs20_top50.csv"
TOP = int(os.environ.get("TOP_N", "50"))

RUNS = int(os.environ.get("FOLDX_RUNS", "20"))
WORKERS = int(os.environ.get("FOLDX_WORKERS", "12"))
GROUPS = os.environ.get("FOLDX_GROUPS", "AB,C")
BINARY = Path(os.environ["FOLDX_BINARY"]).resolve()
PDB = Path(os.environ["DESIGN_PIPELINE_FOLDX_PDB"]).resolve()


def shifts_for(job):
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
        ranked = list(csv.DictReader(handle))[:TOP]
    candidates = read_candidates(MUTATIONS)
    print(f"top {len(ranked)} at FOLDX_RUNS={RUNS}, {WORKERS} workers", flush=True)

    scratch = Path(tempfile.mkdtemp(dir=OUT.parent, prefix="r20work_"))
    jobs = [(r["variant_id"], candidates[r["variant_id"]], scratch) for r in ranked]
    by_id = {r["variant_id"]: r for r in ranked}

    fields = ["rank_n3", "name", "source", "variant_id", "n_subs", "ddg_mean", "ddg_sd",
              "ddg_se", "ddg_min", "ddg_max", "n_runs", "ddg_mean_n3", "ddg_runs"]
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
                    sd = statistics.stdev(s)
                    writer.writerow({
                        "rank_n3": meta["rank"], "name": meta["name"], "source": meta["source"],
                        "variant_id": got["variant_id"], "n_subs": meta["n_subs"],
                        "ddg_mean": round(statistics.mean(s), 6),
                        "ddg_sd": round(sd, 6),
                        "ddg_se": round(sd / len(s) ** 0.5, 6),
                        "ddg_min": round(min(s), 6), "ddg_max": round(max(s), 6),
                        "n_runs": len(s),
                        "ddg_mean_n3": meta["foldx_ddg"],
                        "ddg_runs": ";".join(f"{v:.6f}" for v in s),
                    })
                    handle.flush()
                    results.append(got["variant_id"])
                else:
                    failed.append(got)
                done += 1
                print(f"r20: {done}/{len(jobs)} ({len(results)} measured, {len(failed)} failed)", flush=True)
    shutil.rmtree(scratch, ignore_errors=True)
    print(f"\nwrote {OUT}: {len(results)} rows, {len(failed)} failed", flush=True)
    for f in failed:
        print(f"  FAILED {f['variant_id']}: {f.get('error')}", flush=True)
    return 0


if __name__ == "__main__":
    sys.exit(main())
