"""Several sequence-design variants in one round.

This needs no API support and no new argument: `design_sequences` returns a
plain list, so the arms concatenate. All arms then flow into one scoring pass,
and because `producer` and `meta` record which arm made each design, `select`
picks across arms on merit rather than per arm.

    python examples/fanout.py --dummy
"""

from __future__ import annotations

import sys

from collections import Counter

from design_pipeline import (
    design_sequences,
    make_protein_backbones,
    predict_structure,
    score_designs,
    select,
)

dummy = "--dummy" in sys.argv
target_pdb = "examples/target_stub.pdb" if dummy else "inputs/target.pdb"
bb_model = "dummy_backbones" if dummy else "rfdiffusion"
fold_model = "dummy_structures" if dummy else "colabfold"

# Offline every arm is the same fixture model, so the arms differ only by
# temperature; with the real models they differ by model_type too.
def mpnn(name: str) -> str:
    return "dummy_sequences" if dummy else name


backbones = make_protein_backbones(
    target_pdb,
    contigs="[10-40/A163-181/10-40]",
    num_designs=8,
    seed=0,
    model_name=bb_model,
)

seqs = [
    *design_sequences(backbones, model_name=mpnn("protein_mpnn"), temperature=0.1, seed=0),
    *design_sequences(backbones, model_name=mpnn("protein_mpnn"), temperature=0.3, seed=0),
    *design_sequences(backbones, model_name=mpnn("soluble_mpnn"), temperature=0.2, seed=0),
]
scored = score_designs(predict_structure(seqs, model_name=fold_model), reference=backbones)

best = select(scored, top_k=10, by="pae_interaction")

# The arm is recorded on the sequence design that made it, so a prediction
# reaches it through parent_id. That is what parent_id is for; `lineage` walks
# further when you need the whole chain.
arm = {s.id: (s.producer, s.meta["temperature"]) for s in seqs}

print("designs per arm:")
for (model, temperature), count in sorted(Counter(arm.values()).items()):
    print(f"  {model:<16} T={temperature}  {count}")

print(f"\ntop 10 of {len(scored)} by interface PAE:")
for d in best:
    model, temperature = arm[d.parent_id]
    print(f"  {d.scores['pae_interaction']:6.2f}  {model:<16} T={temperature}  {d.id}")
