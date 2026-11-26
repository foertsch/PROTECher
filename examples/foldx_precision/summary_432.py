import csv
import statistics
from pathlib import Path

R = Path(__file__).resolve().parents[2]
meta = {r["variant_id"]: r for r in csv.DictReader(open(R / "data/track_a_final/candidates.csv", newline=""))}
sc = {r["variant_id"]: float(r["ddg"]) for r in csv.DictReader(open(R / "runs/track_a_foldx/foldx.csv", newline=""))}
# The original 115, i.e. the cut the 317 were measured against. Read from the
# shipped 115-row ranking rather than from a dated backup under $HOME, which
# existed only on the box this was first run on.
old = {r["variant_id"] for r in csv.DictReader(open(R / "runs/rerank/ranking_qualified_v1.csv", newline=""))}
the317 = set(meta) - old
ours = {v for v, m in meta.items() if m["source"] == "ours"}
ranked = sorted(sc, key=lambda v: sc[v])
pos = {v: i + 1 for i, v in enumerate(ranked)}


def summary(label, s):
    d = [sc[v] for v in s]
    print(f"{label:>10}  n={len(d):>3}  best={min(d):+.3f}  median={statistics.median(d):+.3f}  "
          f"mean={statistics.mean(d):+.3f}  worst={max(d):+.3f}")


for lab, s in (("all 432", set(sc)), ("the 115", old), ("the 317", the317),
               ("ours 32", ours), ("kira 400", set(meta) - ours)):
    summary(lab, s)

best317 = min(the317, key=lambda v: sc[v])
print(f"\nbest of the 317: {meta[best317]['name']} rank {pos[best317]} ddg {sc[best317]:+.4f}")
for k in (15, 25, 50, 100, 150):
    print(f"  of the 317, inside the overall top {k:>3}: {sum(1 for v in the317 if pos[v] <= k)}")

ours_ranks = sorted(pos[v] for v in ours)
kira_ranks = sorted(pos[v] for v in set(meta) - ours)
print(f"\nours ranks: {ours_ranks}")
print(f"  ours inside top 100: {sum(1 for r in ours_ranks if r <= 100)} of 32")
print(f"  median rank  ours {statistics.median(ours_ranks):.0f}   kira {statistics.median(kira_ranks):.0f}")
print(f"  ours with ddg < 0: {sum(1 for v in ours if sc[v] < 0)} of 32")
