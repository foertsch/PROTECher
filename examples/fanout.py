"""Several sequence-design variants in one round, folded by several models.

Two kinds of fan-out, and they look different on purpose:

* Sequence arms need no API support and no new argument: `design_sequences`
  returns a plain list, so the arms concatenate.
* Folding takes the list itself: `predict_structure(model_name=[...])` runs
  each named model over the same sequences in turn, so N sequences x M models
  gives N*M predictions.

Either way it is one list of designs and one scoring pass, and `producer` and
`meta` record which arm and which model made each design. A model list means a
mixed pool, and two predictors do not report the same metrics, so filter to one
producer before `select`.

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

# The fold step takes the list directly: every sequence is folded by every
# model named here, in order. Two predictors from different architecture
# families is the cross-model agreement check in docs/filtering_strategy.md —
# same sequences, two opinions. Offline the fixture counterparts stand in.
fold_models = ["dummy_chai", "dummy_boltz"] if dummy else ["chai", "boltz"]
preds = predict_structure(seqs, model_name=fold_models)
scored = score_designs(preds, reference=backbones)

# One pool spanning every predictor. Ranking across them would compare a metric
# one model reports against one it does not, so pick a producer first - a list
# comprehension, no helper - and rerun this part per model.
ranked = [d for d in scored if d.producer == fold_models[0]]
best = select(ranked, top_k=10, by="pae_interaction")

# The arm is recorded on the sequence design that made it, so a prediction
# reaches it through parent_id. That is what parent_id is for; `lineage` walks
# further when you need the whole chain.
arm = {s.id: (s.producer, s.meta["temperature"]) for s in seqs}

print("designs per arm:")
for (model, temperature), count in sorted(Counter(arm.values()).items()):
    print(f"  {model:<16} T={temperature}  {count}")

print(f"\npredictions: {len(preds)} = {len(seqs)} sequences x {len(fold_models)} model(s)")
print(f"\ntop 10 of {len(ranked)} {fold_models[0]} predictions by interface PAE:")
for d in best:
    model, temperature = arm[d.parent_id]
    print(f"  {d.scores['pae_interaction']:6.2f}  {model:<16} T={temperature}  {d.id}")
