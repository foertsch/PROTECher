#!/usr/bin/env bash
# GPU job: does folding in one process (msa_cache.fold_batch) give the same scores as one process per
# design, once the random seed is averaged out? Follow-up to slurm/msa_batch_job.sh (run msa_batch_01),
# where seed 0 gave the same ipTM for the first design of a batch and different values for the rest:
# in one process the random numbers are drawn in a different order, so seed 0 is not the same draw.
#
# Same five sequences and the parent's cached alignments. Seeds 1, 2 and 3 are folded here both ways;
# seed 0 is read from msa_batch_01. The batch runs are this run's reported designs (ids end in the seed),
# the one-process-per-design runs go to the sibling folder <run>_solo.
#
#   slurm/submit_with_report.sh msa_batch_seeds_01 slurm/msa_batch_seeds_job.sh -- --parent wt_model --label batch-vs-solo-seeds
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

from design_pipeline import Design, predict_structure
from design_pipeline.adapters._folding import fasta_chains
from design_pipeline.cache import read_manifest
from msa_cache import MsaCache, fold_batch, use_cached_msas

run_dir = Path(os.environ["RUN_DIR"]).resolve()
solo_dir = run_dir.parent / (run_dir.name + "_solo")
src = Path("runs/test_01").resolve()
ids = ["wt", "v1", "v2", "v3", "v4"]
seqs = {i: src / "sequences" / f"{i}.fa" for i in ids}
designs = [Design(id=i, path=seqs[i], parent_id=None, producer="manual") for i in ids]

cache = MsaCache(run_dir / "msa_parent")
(msa_dir,) = glob.glob(str(src / "predictions" / "*" / "wt" / "boltz_results_wt" / "msa"))
cache.import_boltz_msas(msa_dir, fasta_chains(seqs["wt"])[1])

rows, seconds = [], {"solo": 0.0, "batch": 0.0}


def add(mode, seed, folded):
    for d in folded:
        rows.append({"design": d.parent_id, "mode": mode, "seed": seed, "iptm": d.scores["iptm"], "plddt": d.scores["plddt"]})


# seed 0 from the earlier run
(solo0,) = [Path(p).parent for p in glob.glob("runs/msa_batch_01_perdesign/predictions/*/DONE")]
add("solo", 0, read_manifest(solo0))
add("batch", 0, read_manifest(Path("runs/msa_batch_01/predictions/batch_a")))

for seed in (1, 2, 3):
    t0 = time.perf_counter()
    with use_cached_msas(cache):
        add("solo", seed, predict_structure(designs, output_dir=solo_dir / f"seed{seed}" / "predictions",
                                            model_name="boltz", seed=seed, use_msa_server=False))
    seconds["solo"] += time.perf_counter() - t0
    t0 = time.perf_counter()
    add("batch", seed, fold_batch(cache, designs, run_dir / "predictions" / f"batch_s{seed}", seed=seed))
    seconds["batch"] += time.perf_counter() - t0
    print(f"seed {seed} done", flush=True)

shutil.copy(src / "designs_meta.csv", run_dir / "designs_meta.csv")
df = pd.DataFrame(rows)
df.to_csv(run_dir / "msa_batch_seeds.csv", index=False)
pd.set_option("display.width", 200)

print("\nipTM per design and seed (solo = one process per design, batch = one process for all five)")
print(df.pivot_table(index=["design", "mode"], columns="seed", values="iptm").to_string(float_format=lambda v: f"{v:.3f}"))

stats = df.groupby(["design", "mode"]).iptm.agg(["mean", "std", "min", "max"]).unstack("mode")
print("\nipTM over seeds 0-3 (n=4 per cell)")
print(stats.to_string(float_format=lambda v: f"{v:.3f}"))

wide = df.pivot_table(index=["design", "seed"], columns="mode", values="iptm")
diff = wide["batch"] - wide["solo"]
within = df.groupby(["design", "mode"]).iptm.std().groupby("mode").mean()
print(f"\nbatch minus solo, same design and seed (n={len(diff)}): mean {diff.mean():+.4f}, sd {diff.std():.4f}, "
      f"mean |diff| {diff.abs().mean():.4f}, worst {diff.min():+.3f} / {diff.max():+.3f}")
print(f"average seed-to-seed sd within a design: solo {within['solo']:.4f}, batch {within['batch']:.4f}")
print(f"seed 1-3 wall time for 5 designs x 3 seeds: solo {seconds['solo']:.0f} s, batch {seconds['batch']:.0f} s "
      f"({seconds['solo'] / seconds['batch']:.2f}x)")
PY
