"""Iterative refinement: generate, design, fold, score, keep the best, re-diffuse.

The loop is a `for` statement. There is no campaign object, no DAG and no
scheduler, and adding one would not buy anything the cache does not already
give: rerun this script and every completed round is skipped.

To run it offline with no GPU, pass the dummy models:

    python examples/refine.py --dummy

This is the optional de novo branch, not the default path. The public API is
affinity maturation (see examples/maturation.py); these steps live in
`design_pipeline.denovo` and are reached from `diffuse_loops`.
"""

from __future__ import annotations

import sys

from design_pipeline.denovo import (
    design_sequences,
    make_protein_backbones,
    predict_structure,
)
from design_pipeline.designs import select, to_dataframe
from design_pipeline.scoring import score_designs

contigs = "[10-40/A163-181/10-40]"
n_rounds = 3

# Swapping every step for a fixture model is a model_name, nothing else.
dummy = "--dummy" in sys.argv
target_pdb = "examples/target_stub.pdb" if dummy else "inputs/target.pdb"
models = (
    {"backbones": "dummy_backbones", "seqs": "dummy_sequences", "fold": "dummy_structures"}
    if dummy
    else {"backbones": "rfdiffusion", "seqs": "protein_mpnn", "fold": "colabfold"}
)

backbones = make_protein_backbones(
    target_pdb,
    contigs=contigs,
    num_designs=100,
    seed=0,
    model_name=models["backbones"],
)

for rnd in range(n_rounds):
    seqs = design_sequences(backbones, num_seqs=8, model_name=models["seqs"])
    preds = predict_structure(seqs, model_name=models["fold"])
    scored = score_designs(preds, reference=backbones)
    survivors = select(scored, top_k=8, by="pae_interaction")
    # `backbones` and `survivors` are both PDB files. Keeping the names
    # distinct is the point: reusing one would silently design the next round's
    # sequences against the folded prediction instead of the diffused backbone.
    backbones = make_protein_backbones(
        survivors, partial_T=20, num_designs=12, seed=rnd, model_name=models["backbones"]
    )

print(to_dataframe(scored).sort_values("pae_interaction").head(10))
