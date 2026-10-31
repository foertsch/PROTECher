#!/usr/bin/env bash
# GPU job: does reusing the parent's server alignments (MSAs) change Boltz-2's predictions?
# Needs the run test_01 (slurm/example_gpu_job.sh), whose per-design server alignments Boltz saved.
#
# Conditions, all on the same five sequences (wt = parent, v1 to v4 = test variants):
#   server  the server was called for every sequence                    (run test_01, already done)
#   cached  the parent's alignments reused, with the query row swapped for each variant's sequence
#   none    no alignment (single-sequence mode), for wt, v1 and v4       (lower bound)
# Seed noise: the parent is folded with seeds 0, 1 and 2 using the cached alignments.
#
#   slurm/submit_with_report.sh msa_cache_01 slurm/msa_comparison_job.sh -- --parent wt_model --label cached-msa
set -euo pipefail
. "${PROTECHER_CONFIG:-$HOME/.config/protecher/cluster.env}"
: "${PROTECHER_SCRATCH:?PROTECHER_SCRATCH is not set, copy slurm/cluster.env.example to ~/.config/protecher/cluster.env}"
. "$PROTECHER_SCRATCH/protecher_tools/env.sh"
cd "${CODE:-$PROTECHER_SCRATCH/protecher}"
: "${RUN_DIR:?RUN_DIR is not set, start this through slurm/submit_with_report.sh}"

.venv/bin/python - <<'PY'
import csv
import glob
import os
import shutil
from pathlib import Path

import pandas as pd

from design_pipeline import Design
from design_pipeline.denovo import predict_structure
from design_pipeline.adapters import _boltz
from design_pipeline.adapters._folding import fasta_chains
from design_pipeline.cache import read_manifest
from design_pipeline.scoring import sc_rmsd

run_dir = Path(os.environ["RUN_DIR"]).resolve()
none_dir = run_dir.parent / (run_dir.name.replace("cache", "none") if "cache" in run_dir.name else run_dir.name + "_none")
src = Path("runs/test_01").resolve()
IDS = ["wt", "v1", "v2", "v3", "v4"]

# --- the sequences and the server alignments of the parent
seqs = {i: fasta_chains(src / "sequences" / f"{i}.fa")[1] for i in IDS}  # id -> [scFv, antigen]
msa_dirs = glob.glob(str(src / "predictions" / "*" / "wt" / "boltz_results_wt" / "msa"))
assert len(msa_dirs) == 1, msa_dirs
wt_scfv_csv, wt_antigen_csv = Path(msa_dirs[0]) / "wt_0.csv", Path(msa_dirs[0]) / "wt_1.csv"

# --- cached alignments: the parent's rows, query row swapped for the variant's own sequence
cache = run_dir / "msa_cache"
cache.mkdir(parents=True, exist_ok=True)
with open(wt_scfv_csv, newline="") as fh:
    scfv_rows = list(csv.reader(fh))
header, query = scfv_rows[0], scfv_rows[1]
col = header.index("sequence")
assert query[col] == seqs["wt"][0], "first alignment row is not the parent scFv"
cached = {}
for i in IDS:
    out0, out1 = cache / f"{i}_0.csv", cache / f"{i}_1.csv"
    rows = [list(r) for r in scfv_rows]
    rows[1][col] = seqs[i][0]
    assert len(rows[1][col]) == len(query[col])
    with open(out0, "w", newline="") as fh:
        csv.writer(fh).writerows(rows)
    shutil.copy(wt_antigen_csv, out1)
    cached[i] = [str(out0), str(out1)]
print("cached alignments written for", IDS, "| rows per scFv alignment:", len(scfv_rows) - 1)

# --- let the Boltz adapter point at alignment files (patched here, Kira's package is not edited)
MSAS: dict[str, list[str]] = {}
def write_with_msa(path, sequence):
    _, chains = fasta_chains(sequence.path)
    msas = MSAS[sequence.id]
    path.write_text("".join(f">{chr(ord('A') + k)}|protein|{msas[k]}\n{c}\n" for k, c in enumerate(chains)))
    return path
_boltz.write_boltz_fasta = write_with_msa

designs = {i: Design(id=i, path=src / "sequences" / f"{i}.fa", parent_id=None, producer="manual") for i in IDS}

def fold(ids, seed, out, msa_map):
    MSAS.clear()
    MSAS.update({i: msa_map[i] for i in ids})
    return predict_structure([designs[i] for i in ids], output_dir=out / "predictions",
                             model_name="boltz", seed=seed, use_msa_server=False)

results = []  # (condition, seed, Design)
for seed, ids in ((0, IDS), (1, ["wt"]), (2, ["wt"])):
    results += [("cached", seed, d) for d in fold(ids, seed, run_dir, cached)]
none_map = {i: ["empty", "empty"] for i in IDS}
results += [("none", 0, d) for d in fold(["wt", "v1", "v4"], 0, none_dir, none_map)]

# --- the server run, for comparison
server = {d.parent_id: d for d in read_manifest(next(Path(p).parent for p in glob.glob(str(src / "predictions" / "*" / "DONE"))))}
results = [("server", 0, d) for d in server.values()] + results

rows = []
for cond, seed, d in results:
    ref = server[d.parent_id]
    rows.append({
        "design": d.parent_id, "condition": cond, "seed": seed,
        "ipTM": d.scores["iptm"], "pLDDT": d.scores["plddt"], "pTM": d.scores["ptm"],
        "rmsd_to_server": 0.0 if cond == "server" else sc_rmsd(d, ref),
    })
df = pd.DataFrame(rows).sort_values(["design", "condition", "seed"])
df.to_csv(run_dir / "msa_comparison.csv", index=False)
pd.set_option("display.width", 200)
print(df.to_string(index=False, float_format=lambda v: f"{v:.3f}"))

wt = df[df.design == "wt"]
noise = wt[wt.condition == "cached"]
print("\nseed spread of the parent with cached alignments (seeds 0,1,2): "
      f"ipTM sd {noise.ipTM.std():.4f}, pLDDT sd {noise.pLDDT.std():.3f}, "
      f"max RMSD between seeds vs server {noise.rmsd_to_server.max():.2f} A")
piv = df[df.seed == 0].pivot(index="design", columns="condition", values="ipTM")
print("ipTM, cached minus server per design:", (piv["cached"] - piv["server"]).round(4).to_dict())
print("ipTM, none minus server per design:", (piv["none"] - piv["server"]).dropna().round(4).to_dict())
PY
cp "${CODE:-$PROTECHER_SCRATCH/protecher}/runs/test_01/designs_meta.csv" "$RUN_DIR/designs_meta.csv"
