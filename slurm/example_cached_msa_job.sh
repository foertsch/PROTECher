#!/usr/bin/env bash
# Example GPU job: fold the parent and its variants with the parent's cached alignments, so the
# public alignment server is never called per variant (see msa_cache/cache.py).
#
# Step 1 keeps the alignments an earlier Boltz run saved for the parent (here run test_01).
#   For a new antibody use build_parent instead (step 3), which calls the server once.
# Step 2 folds all designs with the cache and use_msa_server=False.
# Step 3 shows build_parent on its own, with a second cache, as a check of that path.
#
#   slurm/submit_with_report.sh msa_helper_01 slurm/example_cached_msa_job.sh -- --parent wt_model --label msa-helper
set -euo pipefail
. /scratch/foertsch/protecher_tools/env.sh
cd "${CODE:-/scratch/foertsch/protecher}"
: "${RUN_DIR:?RUN_DIR is not set, start this through slurm/submit_with_report.sh}"

.venv/bin/python - <<'PY'
import glob
import os
import shutil
from pathlib import Path

from design_pipeline import Design, predict_structure
from design_pipeline.adapters._folding import fasta_chains
from msa_cache import MsaCache, build_parent, use_cached_msas

run_dir = Path(os.environ["RUN_DIR"]).resolve()
src = Path("runs/test_01").resolve()
ids = ["wt", "v1", "v2", "v3", "v4"]
seqs = {i: src / "sequences" / f"{i}.fa" for i in ids}
parent_chains = fasta_chains(seqs["wt"])[1]

# 1. keep the parent's alignments from the earlier run
cache = MsaCache(run_dir / "msa_parent")
(msa_dir,) = glob.glob(str(src / "predictions" / "*" / "wt" / "boltz_results_wt" / "msa"))
cache.import_boltz_msas(msa_dir, parent_chains)
print("cache:", cache.info())

# 2. fold every design with the cache, no server calls
designs = [Design(id=i, path=seqs[i], parent_id=None, producer="manual") for i in ids]
with use_cached_msas(cache):
    folded = predict_structure(designs, output_dir=run_dir / "predictions", model_name="boltz", seed=0, use_msa_server=False)
for d in sorted(folded, key=lambda d: d.parent_id):
    print(f"{d.parent_id:3s} ipTM {d.scores['iptm']:.3f} pLDDT {d.scores['plddt']:.2f} pTM {d.scores['ptm']:.3f}")
shutil.copy(src / "designs_meta.csv", run_dir / "designs_meta.csv")

# 3. build_parent on its own: one server call, then the same cache API
built = MsaCache(run_dir.parent / (run_dir.name + "_built"))
build_parent(built, parent_chains, work_dir=run_dir.parent / (run_dir.name + "_built_fold"))
print("built cache:", built.info())
PY
