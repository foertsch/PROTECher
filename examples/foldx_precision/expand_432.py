"""Expand data/track_a_final/candidates.csv (one row per variant) into the
per-mutation table the FoldX wrapper reads, via the package's own parse_id."""
import csv
from pathlib import Path

from design_pipeline.variants import parse_id

SRC = Path("data/track_a_final/candidates.csv")
OUT = Path("runs/track_a_foldx/candidates_all432.csv")
SCORED = Path("runs/track_a_foldx/foldx.csv")

with SRC.open(newline="") as handle:
    rows = list(csv.DictReader(handle))

with SCORED.open(newline="") as handle:
    scored = {r["variant_id"] for r in csv.DictReader(handle) if r.get("variant_id")}

with OUT.open("w", newline="") as handle:
    writer = csv.writer(handle)
    writer.writerow(["variant_id", "chain", "position", "wt", "mut"])
    n = 0
    for row in rows:
        for chain, position, wt, mut in parse_id(row["variant_id"]):
            writer.writerow([row["variant_id"], chain, position, wt, mut])
            n += 1

ids = [r["variant_id"] for r in rows]
assert len(set(ids)) == len(ids) == 432, (len(ids), len(set(ids)))
print(f"wrote {OUT}: {len(ids)} variants, {n} mutation rows")

ours = [r for r in rows if r["source"] == "ours"]
ours_pending = [r for r in ours if r["variant_id"] not in scored]
print(f"ours: {len(ours)} ({ours[0]['name']}..{ours[-1]['name']})")
print(f"  already scored: {len(ours) - len(ours_pending)}")
print(f"  pending this run: {len(ours_pending)}")
print(f"  names pending: {','.join(r['name'] for r in ours_pending)}")
print(f"total pending: {len(set(ids) - scored)} of 432 ({len(scored)} already scored)")
