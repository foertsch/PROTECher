#!/usr/bin/env bash
# GPU job: how much faster is it to fold all variants in one Boltz process (msa_cache.fold_batch)
# than one process per design (the adapter's way)? Both use the parent's cached alignments, so no
# alignment server is involved. Needs the run test_01 (slurm/example_gpu_job.sh) for the five sequences
# and the parent's alignments.
#
# Runs, in this order, on the same node and GPU:
#   per-design  predict_structure under use_cached_msas        (sibling folder <run>_perdesign)
#   batch A     fold_batch, seed 0                             (this run's predictions, reported)
#   batch B     fold_batch again, seed 0                       (sibling folder <run>_repeat)
# Printed: wall time of each, then ipTM / pLDDT per design and the RMSD between structures.
# Folding in one process draws the diffusion noise in a different order than five processes do, so
# the scores are expected to differ by seed noise, not to be identical.
#
#   slurm/submit_with_report.sh msa_batch_01 slurm/msa_batch_job.sh -- --parent wt_model --label batched-fold
set -euo pipefail
. "${PROTECHER_CONFIG:-$HOME/.config/protecher/cluster.env}"
: "${PROTECHER_SCRATCH:?PROTECHER_SCRATCH is not set, copy slurm/cluster.env.example to ~/.config/protecher/cluster.env}"
. "$PROTECHER_SCRATCH/protecher_tools/env.sh"
cd "${CODE:-$PROTECHER_SCRATCH/protecher}"
: "${RUN_DIR:?RUN_DIR is not set, start this through slurm/submit_with_report.sh}"

.venv/bin/python - <<'PY'
import glob
import os
import shutil
import time
from pathlib import Path

import pandas as pd

from design_pipeline import Design
from design_pipeline.denovo import predict_structure
from design_pipeline.adapters._folding import fasta_chains
from design_pipeline.scoring import sc_rmsd
from msa_cache import MsaCache, fold_batch, use_cached_msas

run_dir = Path(os.environ["RUN_DIR"]).resolve()
src = Path("runs/test_01").resolve()
ids = ["wt", "v1", "v2", "v3", "v4"]
seqs = {i: src / "sequences" / f"{i}.fa" for i in ids}
designs = [Design(id=i, path=seqs[i], parent_id=None, producer="manual") for i in ids]

cache = MsaCache(run_dir / "msa_parent")
(msa_dir,) = glob.glob(str(src / "predictions" / "*" / "wt" / "boltz_results_wt" / "msa"))
cache.import_boltz_msas(msa_dir, fasta_chains(seqs["wt"])[1])
print("cache:", cache.info(), flush=True)

timings = {}

t0 = time.perf_counter()
with use_cached_msas(cache):
    per_design = predict_structure(designs, output_dir=run_dir.parent / (run_dir.name + "_perdesign") / "predictions",
                                   model_name="boltz", seed=0, use_msa_server=False)
timings["per-design"] = time.perf_counter() - t0

t0 = time.perf_counter()
batch_a = fold_batch(cache, designs, run_dir / "predictions" / "batch_a", seed=0)
timings["batch A"] = time.perf_counter() - t0

t0 = time.perf_counter()
batch_b = fold_batch(cache, designs, run_dir.parent / (run_dir.name + "_repeat") / "predictions" / "batch_b", seed=0)
timings["batch B"] = time.perf_counter() - t0

shutil.copy(src / "designs_meta.csv", run_dir / "designs_meta.csv")

n = len(designs)
print("\nwall time for", n, "designs (includes starting the processes and parsing the output):")
for name, seconds in timings.items():
    print(f"  {name:10s} {seconds:7.1f} s   {seconds / n:6.1f} s per design")
print(f"  speedup, batch A vs per-design: {timings['per-design'] / timings['batch A']:.2f}x, "
      f"batch B vs per-design: {timings['per-design'] / timings['batch B']:.2f}x")

rows = []
by = {name: {d.parent_id: d for d in folded} for name, folded in
      (("per-design", per_design), ("batch A", batch_a), ("batch B", batch_b))}
for i in ids:
    ref = by["per-design"][i]
    for name, folded in by.items():
        d = folded[i]
        rows.append({"design": i, "run": name, "ipTM": d.scores["iptm"], "pLDDT": d.scores["plddt"], "pTM": d.scores["ptm"],
                     "rmsd_to_per_design": 0.0 if name == "per-design" else sc_rmsd(d, ref)})
df = pd.DataFrame(rows)
df.to_csv(run_dir / "msa_batch_comparison.csv", index=False)
pd.set_option("display.width", 200)
print("\n" + df.to_string(index=False, float_format=lambda v: f"{v:.3f}"))
piv = df.pivot(index="design", columns="run", values="ipTM")
print("\nipTM, batch A minus per-design:", (piv["batch A"] - piv["per-design"]).round(4).to_dict())
print("ipTM, batch B minus batch A   :", (piv["batch B"] - piv["batch A"]).round(4).to_dict())
PY
