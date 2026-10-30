#!/usr/bin/env bash
# Example GPU job for slurm/submit_with_report.sh: the parent trastuzumab scFv plus four plumbing-test
# variants, folded with Boltz-2 together with a HER2 fragment. The variants are conservative swaps at
# framework positions, chosen to exercise the report (parent deltas, groups, mutation strings). They are
# NOT designs. Numbering is sequential within each domain, not Kabat.
# All sequences come from the public PDB entry 1N8Z.
#
#   slurm/submit_with_report.sh test_01 slurm/example_gpu_job.sh -- --parent wt_model --label plumbing-test
set -euo pipefail
. "${PROTECHER_CONFIG:-$HOME/.config/protecher/cluster.env}"
: "${PROTECHER_SCRATCH:?PROTECHER_SCRATCH is not set, copy slurm/cluster.env.example to ~/.config/protecher/cluster.env}"
. "$PROTECHER_SCRATCH/protecher_tools/env.sh"
cd "${CODE:-$PROTECHER_SCRATCH/protecher}"
: "${RUN_DIR:?RUN_DIR is not set, start this through slurm/submit_with_report.sh}"

.venv/bin/python - <<'PY'
import csv
import os
import random
import urllib.request
from pathlib import Path

from design_pipeline import Design, predict_structure

run_dir = Path(os.environ["RUN_DIR"])
(run_dir / "sequences").mkdir(parents=True, exist_ok=True)

text = urllib.request.urlopen("https://www.rcsb.org/fasta/entry/1N8Z/display", timeout=30).read().decode()
chains = {}
for block in text.split(">")[1:]:
    head, *seq = block.strip().split("\n")
    chains[head.split("|")[1].split()[1].rstrip(",")] = "".join(seq)
vl, vh = chains["A"][:107], chains["B"][:120]
antigen = chains["C"][-128:]
LINKER = "GGGGS" * 3
assert len(vh + LINKER + vl) == 242

SWAP = {"S": "T", "T": "S", "K": "R", "R": "K", "N": "D", "D": "E", "E": "D", "Q": "E", "A": "S", "V": "I", "I": "V"}
# framework stretches (sequential numbering, 1-based, inclusive), away from the CDRs
FRAMEWORK = {"VH": (vh, [(1, 25), (66, 94)]), "VL": (vl, [(1, 23), (57, 88)])}


def variant(seed: int) -> tuple[str, str]:
    rng = random.Random(seed)
    candidates = [
        (name, pos)
        for name, (seq, spans) in FRAMEWORK.items()
        for lo, hi in spans
        for pos in range(lo, hi + 1)
        if seq[pos - 1] in SWAP
    ]
    chosen = sorted(rng.sample(candidates, 5))
    new = {"VH": list(vh), "VL": list(vl)}
    labels = []
    for name, pos in chosen:
        old = new[name][pos - 1]
        new[name][pos - 1] = SWAP[old]
        labels.append(f"{name}:{old}{pos}{SWAP[old]}")
    return "".join(new["VH"]) + LINKER + "".join(new["VL"]), "+".join(labels)


entries = [("wt", vh + LINKER + vl, "", "parent")]
for i in range(1, 5):
    seq, muts = variant(seed=i)
    entries.append((f"v{i}", seq, muts, "test_a" if i <= 2 else "test_b"))

designs = []
for design_id, scfv, _, _ in entries:
    path = run_dir / "sequences" / f"{design_id}.fa"
    path.write_text(f">{design_id}\n{scfv}:{antigen}\n")
    designs.append(Design(id=design_id, path=path, parent_id=None, producer="manual"))

with open(run_dir / "designs_meta.csv", "w", newline="") as fh:
    writer = csv.writer(fh)
    writer.writerow(["design", "mutations", "group", "n_mutations"])
    for design_id, _, muts, group in entries:
        writer.writerow([design_id, muts, group, len(muts.split("+")) if muts else 0])

predict_structure(designs, output_dir=run_dir / "predictions", model_name="boltz", seed=0)
print("folded", len(designs), "sequences into", run_dir)
PY
