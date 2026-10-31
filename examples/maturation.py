"""Affinity maturation end to end: structure, scan, score, gate, pair, score.

The whole pipeline is six calls. Steps 3 and 6 are the same function: a pair is
an ordinary candidate, so nothing here adds two single scores together to guess
at one.

To run it offline with no GPU, no Docker and no network:

    python examples/maturation.py --dummy

With the real models, `--dummy` off, this wants the AntiFold and ESM-IF images
and a GPU; see INSTALL.md.
"""

from __future__ import annotations

import sys

from design_pipeline import (
    enumerate_variants,
    gate,
    get_structure,
    score_variants,
    select_pairs,
)

dummy = "--dummy" in sys.argv
out = "runs/maturation_dummy" if dummy else "runs/her2_v1"

# Step 1. The experimental structure, because coordinates beat a prediction.
# Offline, the fixture complex stands in for the download.
complex_ = (
    get_structure(out, model_name="dummy_complex")
    if dummy
    else get_structure(out, pdb_id="1N8Z", chains=["A", "B", "C"])
)
chain = "H" if dummy else "B"  # 1N8Z labels trastuzumab's heavy chain B
print(f"structure: {complex_.path.name}, {len(complex_)} model(s) in the ensemble")

# Step 2. Every single mutation of the chain: one row each, 19 per position.
variants = enumerate_variants(complex_, chain=chain)
print(f"variants:  {variants.candidates} candidates in {variants.path.name}")

# Step 3. Fan out over the cheap scorers. One column each, joined on
# variant_id; nothing is averaged into a combined score, because which column
# to trust comes out of the harness rather than out of a guess here.
cheap = ["dummy_scorer"] if dummy else ["antifold", "esm_if"]
scores = score_variants(complex_, variants, models=cheap)
print(f"scored:    {', '.join(scores.score_columns)}")

# Step 4. Drop the confident developability failures. `liabilities` is pure
# sequence regex and costs nothing; the container checks are opt-in.
survivors = gate(scores, checks=["liabilities"])
print(f"survivors: {survivors.candidates} of {scores.candidates}")

# Rank on the one column you decided to trust. This is the caller's job on
# purpose - the library ships no weights.
ranked = survivors.frame.sort_values(f"score_{cheap[0]}", ascending=False)
print(f"\ntop 10 singles by score_{cheap[0]}:")
for row in ranked.head(10).itertuples():
    print(f"  {getattr(row, 'score_' + cheap[0]):7.3f}  {row.variant_id}")

# Step 5. Combine the best singles into pairs at least 10 A apart, CA to CA.
pairs = select_pairs(
    survivors, max_pairs=20, min_distance=10.0, by=f"score_{cheap[0]}"
)
print(f"\npairs:     {pairs.candidates}")

# Step 6. Score the pairs. Same function, same table shape. `antifold` is a
# lookup and has no answer for two mutations, so the pairs go to `antifold_seq`,
# which scores each as a whole mutant and needs the chain ids (1N8Z: heavy B,
# light A). The per-variant scorer is capped, because it runs one pass per
# candidate. GearBind is left out: it has no image here, and a flag one model
# rejects fails the whole call, so it would need a call of its own.
expensive = ["dummy_scorer", "dummy_per_variant"] if dummy else ["antifold_seq"]
chains = {} if dummy else {"heavy": "B", "light": "A"}
pair_scores = score_variants(complex_, pairs, models=expensive, max_variants=20, **chains)

frame = pair_scores.frame.drop_duplicates("variant_id")
print(f"\npairs by score_{expensive[0]}:")
for row in frame.sort_values(f"score_{expensive[0]}", ascending=False).itertuples():
    print(f"  {getattr(row, 'score_' + expensive[0]):7.3f}  {row.variant_id}")

print(f"\neverything is on disk under {out}/")
