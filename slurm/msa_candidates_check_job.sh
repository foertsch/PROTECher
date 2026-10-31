#!/usr/bin/env bash
# GPU job: check the alignments that slurm/msa_candidates_job.sh stored, using only the stored files.
#
# 1. Each candidate's known complex (scFv from the PDB Fab plus the whole UniProt extracellular domain)
#    is folded with its own cached alignments and no server call. Boltz-2 should give a confident
#    interface for a known drug and its target, as it did for trastuzumab with a HER2 fragment (ipTM 0.94).
# 2. Fragment slicing, on trastuzumab: the five sequences of run test_01 carry a 128-residue HER2 fragment.
#    They are folded with the alignment cut out of the cached 630-residue domain alignment, and compared
#    with the server run (test_01, server called per design) and with the run that reused test_01's own
#    fragment alignment (msa_helper_01). Same seed 0 throughout, so structures can be compared directly.
#
#   CACHE_ROOT  default $PROTECHER_SCRATCH/protecher_cache/msa_candidates
#   ONLY        comma-separated candidate names for step 1, default all
#
# Needs a free GPU on the node that holds the install. From the laptop:
#   ssh login '. ~/.config/protecher/cluster.env; sbatch --parsable -p $PROTECHER_GPU_PARTITION -w $PROTECHER_GPU_NODE --gres=gpu:1 --cpus-per-task=8 --mem=32G \
#     -t 01:00:00 -J msa_candidates_check -o $HOME/protecher_logs/msa_candidates_check_%j.out' < slurm/msa_candidates_check_job.sh
set -euo pipefail
. "${PROTECHER_CONFIG:-$HOME/.config/protecher/cluster.env}"
: "${PROTECHER_SCRATCH:?PROTECHER_SCRATCH is not set, copy slurm/cluster.env.example to ~/.config/protecher/cluster.env}"
. "$PROTECHER_SCRATCH/protecher_tools/env.sh"
cd "${CODE:-$PROTECHER_SCRATCH/protecher}"

.venv/bin/python - <<'PY'
import glob
import json
import os
import sys
from pathlib import Path

import pandas as pd

from design_pipeline import Design
from design_pipeline.denovo import predict_structure
from design_pipeline.adapters._folding import fasta_chains
from design_pipeline.cache import read_manifest
from design_pipeline.scoring import sc_rmsd
from msa_cache import MsaCache, use_cached_msas

root = Path(os.environ.get("CACHE_ROOT", os.environ["PROTECHER_SCRATCH"] + "/protecher_cache/msa_candidates"))
only = [n for n in os.environ.get("ONLY", "").split(",") if n]
names = [p.parent.name for p in sorted(root.glob("*/parent.json")) if not only or p.parent.name in only]
failed = []

print("1. known complexes folded from their own stored alignments\n")
rows = []
for name in names:
    try:
        cache = MsaCache(root / name)
        fasta = root / "_check" / name / "complex.fa"
        fasta.parent.mkdir(parents=True, exist_ok=True)
        fasta.write_text(f">{name}\n" + ":".join(cache.parent_chains) + "\n")
        design = Design(id=name, path=fasta, parent_id=None, producer="manual")
        with use_cached_msas(cache):
            (folded,) = predict_structure([design], output_dir=root / "_check" / name / "predictions",
                                          model_name="boltz", seed=0, use_msa_server=False)
        scores = {k: round(folded.scores[k], 4) for k in ("iptm", "plddt", "ptm") if k in folded.scores}
        info_path = root / name / "candidate.json"
        info = json.loads(info_path.read_text())
        info["validation_fold"] = scores
        info_path.write_text(json.dumps(info, indent=1))
        rows.append({"candidate": name, "scFv": len(cache.parent_chains[0]), "antigen": len(cache.parent_chains[1]), **scores})
        print(f"  {name}: {scores}", flush=True)
    except Exception as error:
        failed.append(name)
        print(f"  FAILED {name}: {type(error).__name__}: {error}", file=sys.stderr, flush=True)
if rows:
    print("\n" + pd.DataFrame(rows).to_string(index=False, float_format=lambda v: f"{v:.3f}"))

print("\n2. fragment slicing on trastuzumab (128-residue HER2 fragment cut from the 630-residue alignment)\n")
try:
    cache = MsaCache(root / "trastuzumab_her2")
    src = Path("runs/test_01").resolve()
    ids = ["wt", "v1", "v2", "v3", "v4"]
    designs = [Design(id=i, path=src / "sequences" / f"{i}.fa", parent_id=None, producer="manual") for i in ids]
    fragment = fasta_chains(designs[0].path)[1][1]
    offset = cache.parent_chains[1].find(fragment)
    print(f"  fragment of {len(fragment)} residues sits at position {offset + 1} of the cached {len(cache.parent_chains[1])}-residue domain")
    with use_cached_msas(cache):
        sliced = predict_structure(designs, output_dir=root / "_check" / "trastuzumab_fragment" / "predictions",
                                   model_name="boltz", seed=0, use_msa_server=False)
    server = {d.parent_id: d for d in read_manifest(Path(glob.glob(str(src / "predictions" / "*" / "DONE"))[0]).parent)}
    own = {d.parent_id: d for d in read_manifest(Path(glob.glob("runs/msa_helper_01/predictions/*/DONE")[0]).parent)}
    sliced = {d.parent_id: d for d in sliced}
    table = pd.DataFrame([{
        "design": i,
        "ipTM server": server[i].scores["iptm"], "ipTM own fragment msa": own[i].scores["iptm"], "ipTM sliced": sliced[i].scores["iptm"],
        "pLDDT sliced": sliced[i].scores["plddt"], "rmsd sliced vs server": sc_rmsd(sliced[i], server[i]),
    } for i in ids])
    table.to_csv(root / "_check" / "trastuzumab_fragment" / "slice_vs_server.csv", index=False)
    print(table.to_string(index=False, float_format=lambda v: f"{v:.3f}"))
    d = table["ipTM sliced"] - table["ipTM server"]
    print(f"\n  ipTM sliced minus server: mean {d.mean():+.4f}, worst {d.min():+.3f} / {d.max():+.3f}")
except Exception as error:
    failed.append("fragment slicing")
    print(f"  FAILED fragment slicing: {type(error).__name__}: {error}", file=sys.stderr, flush=True)

if failed:
    print("\nFAILED:", ", ".join(failed))
    sys.exit(1)
PY
