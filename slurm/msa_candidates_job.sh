#!/usr/bin/env bash
# CPU job: fetch and store the Boltz-2 alignments for the antibody/antigen pairs the Track A challenge
# could announce (msa_cache/candidates.py), so tomorrow's first folds need no call to the public
# alignment server. No GPU and no fold: this is the step Boltz runs before folding
# (msa_cache.fetch_alignments). slurm/msa_candidates_check_job.sh folds each known complex from the
# stored alignments afterwards, which needs a GPU.
#
# Output, per candidate:  $CACHE_ROOT/<name>/{parent.json,parent_0.csv,parent_1.csv,candidate.json}
#   CACHE_ROOT  node-local scratch, where the jobs read it   (default $PROTECHER_SCRATCH/protecher_cache/msa_candidates)
#   BACKUP_DIR  home on the shared filesystem, a backup in case scratch is cleaned (default ~/protecher_msa_candidates)
#   ONLY        comma-separated candidate names to fetch, default all
# A candidate whose cache already holds the same chains is skipped, so the job can be rerun.
#
# Not a design run, so it is not submitted with submit_with_report.sh. Pinned to the node that holds the
# install, without a GPU, so it does not wait behind GPU jobs. From the laptop:
#   ssh login '. ~/.config/protecher/cluster.env; sbatch --parsable -p $PROTECHER_GPU_PARTITION -w $PROTECHER_GPU_NODE --cpus-per-task=4 --mem=16G \
#     -t 01:00:00 -J msa_candidates -o $HOME/protecher_logs/msa_candidates_%j.out' < slurm/msa_candidates_job.sh
set -euo pipefail
. "${PROTECHER_CONFIG:-$HOME/.config/protecher/cluster.env}"
: "${PROTECHER_SCRATCH:?PROTECHER_SCRATCH is not set, copy slurm/cluster.env.example to ~/.config/protecher/cluster.env}"
. "$PROTECHER_SCRATCH/protecher_tools/env.sh"
cd "${CODE:-$PROTECHER_SCRATCH/protecher}"

.venv/bin/python - <<'PY'
import json
import os
import shutil
import sys
import time
from pathlib import Path

from msa_cache import CANDIDATES, MsaCache, fetch_alignments, resolve

root = Path(os.environ.get("CACHE_ROOT", os.environ["PROTECHER_SCRATCH"] + "/protecher_cache/msa_candidates"))
backup = Path(os.environ.get("BACKUP_DIR", str(Path.home() / "protecher_msa_candidates")))
only = [n for n in os.environ.get("ONLY", "").split(",") if n]
root.mkdir(parents=True, exist_ok=True)
backup.mkdir(parents=True, exist_ok=True)

rows, failed = [], []
for candidate in CANDIDATES:
    if only and candidate.name not in only:
        continue
    started = time.perf_counter()
    try:
        info = resolve(candidate)
        chains = [info["scfv"], info["antigen"]]
        print(f"\n=== {candidate.name}: scFv {len(chains[0])} (VH {info['vh_length']} + VL {info['vl_length']}), "
              f"antigen {len(chains[1])} (UniProt {candidate.uniprot} {info['antigen_range']})", flush=True)
        if not info["pdb_antigen_inside_domain"]:
            print(f"  note: the {info['pdb_antigen_length']}-residue antigen in {candidate.pdb} is not a stretch of the UniProt domain")
        cache = MsaCache(root / candidate.name)
        if cache.has_parent() and cache.parent_chains == chains:
            print("  cache already holds these chains, skipped", flush=True)
        else:
            shutil.rmtree(root / candidate.name)
            cache = MsaCache(root / candidate.name)
            fetch_alignments(cache, chains, work_dir=root / "_work" / candidate.name)
        info.update({**cache.info(), "built_in_seconds": round(time.perf_counter() - started)})
        (root / candidate.name / "candidate.json").write_text(json.dumps(info, indent=1))
        shutil.copytree(root / candidate.name, backup / candidate.name, dirs_exist_ok=True)
        rows.append(info)
        print(f"  alignment rows {info['alignment_rows']}, {info['built_in_seconds']} s", flush=True)
    except Exception as error:  # one failed candidate must not stop the others
        failed.append(candidate.name)
        print(f"  FAILED {candidate.name}: {type(error).__name__}: {error}", file=sys.stderr, flush=True)

index = [{k: r[k] for k in ("name", "pdb", "uniprot", "chain_lengths", "alignment_rows")} for r in rows]
(root / "index.json").write_text(json.dumps(index, indent=1))
shutil.copy(root / "index.json", backup / "index.json")

print("\nname                 scFv  antigen  rows scFv / antigen")
for r in rows:
    print(f"{r['name']:20s} {r['chain_lengths'][0]:4d}  {r['chain_lengths'][1]:7d}  {r['alignment_rows'][0]:7d} / {r['alignment_rows'][1]:<7d}")
print(f"\nbackup in {backup}: {sum(p.stat().st_size for p in backup.rglob('*') if p.is_file()) / 1e6:.0f} MB")
if failed:
    print("FAILED:", ", ".join(failed))
    sys.exit(1)
PY
