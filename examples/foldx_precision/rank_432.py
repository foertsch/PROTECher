"""qualified_v1 = z(-foldx_ddg) over all 432 track A designs, plus the figures asked for.

Ordering is on raw ddg (monotone in it); the z only labels the rows for presentation.
The original 115 are defined by the foldx.csv backup taken before this session added
anything, so "the 317" means exactly the designs that cut did not cover.
"""
import csv
import sys
from pathlib import Path

import numpy as np
from scipy.stats import spearmanr

REPO = Path(__file__).resolve().parents[2]
FINAL = REPO / "data/track_a_final/candidates.csv"
FOLDX = REPO / "runs/track_a_foldx/foldx.csv"
BASE115 = Path(sys.argv[1])                      # pre-session foldx.csv backup
OLD_RANK = REPO / "runs/rerank/ranking_qualified_v1.csv"
OUT = REPO / "runs/rerank/ranking_qualified_v1_432.csv"


def read(path):
    with open(path, newline="") as handle:
        return list(csv.DictReader(handle))


meta = {r["variant_id"]: r for r in read(FINAL)}
scores = {r["variant_id"]: r for r in read(FOLDX) if r.get("variant_id")}
old115 = {r["variant_id"] for r in read(BASE115) if r.get("variant_id")}
the317 = set(meta) - old115

print(f"candidates: {len(meta)}   scored: {len(scores)}   original 115: {len(old115)}   the 317: {len(the317)}")

# ---- 6. failures first: anything in the table with no row in foldx.csv
missing = [v for v in meta if v not in scores]
print(f"\n[6] failed to score: {len(missing)}")
for v in missing:
    print(f"    {meta[v]['name']:>5}  {v}")

ranked = sorted((v for v in meta if v in scores), key=lambda v: float(scores[v]["ddg"]))
ddg = np.array([float(scores[v]["ddg"]) for v in ranked])
neg = -ddg
# ddof=0: the package's own `_z` in examples/track_a_chem_rank.py uses
# statistics.pstdev, so the sample sd would put this column on a different
# scale from every ranking file already on disk. Ordering is unaffected either
# way - it is monotone in ddg - but the number has to mean the same thing.
z = (neg - neg.mean()) / neg.std(ddof=0)

with OUT.open("w", newline="") as handle:
    writer = csv.writer(handle)
    writer.writerow(["rank", "name", "source", "variant_id", "profile", "n_subs",
                     "priority_or_additive", "foldx_ddg", "ddg_fold", "score"])
    for i, v in enumerate(ranked):
        m = meta[v]
        writer.writerow([i + 1, m["name"], m["source"], v, "qualified_v1", m["n_subs"],
                         m["priority_or_additive"], round(float(scores[v]["ddg"]), 4),
                         round(float(scores[v]["ddg_fold"]), 4), round(float(z[i]), 4)])
print(f"\nwrote {OUT}: {len(ranked)} ranked rows")

# ---- old top 15 over the 115, from the shipped ranking file
old_rows = read(OLD_RANK)
old_top15 = [r["variant_id"] for r in old_rows[:15]]
recomputed = sorted(old115 & set(scores), key=lambda v: float(scores[v]["ddg"]))[:15]
print(f"\nold top 15 reproduced from ddg over the 115: {old_top15 == recomputed}")

# ---- 1. how many of the 317 beat the current top 15 on measured ddG
cutoff = float(scores[old_top15[14]]["ddg"])
best_old = float(scores[old_top15[0]]["ddg"])
beat_cutoff = [v for v in the317 if v in scores and float(scores[v]["ddg"]) < cutoff]
beat_best = [v for v in the317 if v in scores and float(scores[v]["ddg"]) < best_old]
print(f"\n[1] old top-15 cutoff (15th best of the 115) ddg = {cutoff:+.4f}")
print(f"    of the 317, better than that cutoff : {len(beat_cutoff)}")
print(f"    of the 317, better than the old #1 ({best_old:+.4f}): {len(beat_best)}")

# ---- 2. overlap of the new top 15 over 432 with the old top 15 over 115
new_top15 = ranked[:15]
overlap = [v for v in new_top15 if v in old_top15]
print(f"\n[2] new top 15 over 432 vs old top 15 over 115: {len(overlap)} of 15 shared")
print(f"    kept: {', '.join(meta[v]['name'] for v in overlap) or 'none'}")
print(f"    dropped out: {', '.join(meta[v]['name'] for v in old_top15 if v not in new_top15) or 'none'}")

# ---- 3. how many of the 32 ours make the new top 15
ours = {v for v, m in meta.items() if m["source"] == "ours"}
ours_top15 = [v for v in new_top15 if v in ours]
print(f"\n[3] of the 32 'ours', in the new top 15: {len(ours_top15)}")
print(f"    {', '.join(meta[v]['name'] for v in ours_top15) or 'none'}")
ours_ranks = [(i + 1, meta[v]["name"]) for i, v in enumerate(ranked) if v in ours]
print(f"    best 'ours' rank: {ours_ranks[0][0]} ({ours_ranks[0][1]}) of {len(ranked)}")

# ---- 4. Spearman of the additive proxy against measured ddg across the 115
pairs = [(float(meta[v]["priority_or_additive"]), float(scores[v]["ddg"]))
         for v in old115 if v in scores]
rho, p = spearmanr([a for a, _ in pairs], [b for _, b in pairs])
print(f"\n[4] Spearman(priority_or_additive, measured ddg) over the 115: rho = {rho:+.4f}  p = {p:.3g}  n = {len(pairs)}")
sources = {meta[v]["source"] for v in old115}
print(f"    sources inside the 115: {sorted(sources)}")

# ---- 5. ddg > 0 over all 432
pos = [v for v in ranked if float(scores[v]["ddg"]) > 0]
print(f"\n[5] ddg > 0 (weaker binding than wild type) over all {len(ranked)} scored: {len(pos)}")
for label, subset in (("the 115", old115), ("the 317", the317), ("ours", ours)):
    inset = [v for v in subset if v in scores]
    n = sum(1 for v in inset if float(scores[v]["ddg"]) > 0)
    print(f"    {label:>8}: {n} of {len(inset)}")

print(f"\nnew top 15 over 432")
print(f"{'#':>3} {'name':<6}{'source':<7}{'ddg':>9}{'z':>8}  variant")
for i, v in enumerate(new_top15):
    m = meta[v]
    print(f"{i+1:>3} {m['name']:<6}{m['source']:<7}{float(scores[v]['ddg']):9.4f}{z[i]:8.3f}  {v}")
