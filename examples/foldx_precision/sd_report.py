"""The four figures asked for, from the 60-candidate per-run sample."""
import csv
import statistics
from pathlib import Path

R = Path(__file__).resolve().parents[2]
sd_rows = list(csv.DictReader(open(R / "runs/track_a_foldx/foldx_sd_sample.csv", newline="")))
rank_rows = list(csv.DictReader(open(R / "runs/rerank/ranking_qualified_v1_432.csv", newline="")))

ddg_all = [float(r["foldx_ddg"]) for r in rank_rows]
CUT = ddg_all[14]
sds = sorted(float(r["ddg_sd"]) for r in sd_rows)
n = len(sds)


def pct(values, q):
    return values[min(len(values) - 1, int(q * len(values)))]


print(f"sample n = {n}, 0 failed")
print(f"  median SD = {statistics.median(sds):.4f}")
print(f"  p90    SD = {pct(sds, 0.90):.4f}")
print(f"  mean   SD = {statistics.mean(sds):.4f}")
print(f"  min/max   = {sds[0]:.4f} / {sds[-1]:.4f}")
print(f"  SD == 0   = {sum(1 for s in sds if s == 0)} of {n}")
print(f"  SD > 0.1  = {sum(1 for s in sds if s > 0.1)} of {n}")

top30 = [r for r in sd_rows if int(r["rank"]) <= 30]
rest = [r for r in sd_rows if int(r["rank"]) > 30]
for label, group in (("ranks 1-30", top30), ("ranks 31-432", rest)):
    g = sorted(float(r["ddg_sd"]) for r in group)
    print(f"  {label:>13}: n={len(g)} median={statistics.median(g):.4f} p90={pct(g, 0.90):.4f} max={g[-1]:.4f}")

# does the mean of three runs reproduce what foldx.csv recorded?
drift = [abs(float(r["ddg_mean"]) - float(r["ddg_in_foldx_csv"])) for r in sd_rows]
print(f"\nremeasured mean vs the stored ddg: max |diff| = {max(drift):.4f}, median {statistics.median(drift):.4f}")
print(f"  within 0.01 of the stored value: {sum(1 for d in drift if d <= 0.01)} of {n}")

MED = statistics.median(sds)
P90 = pct(sds, 0.90)
print("\ndesigns within k*SD of the rank-15 cutoff "
      f"({rank_rows[14]['name']} at {CUT:+.4f})")
for label, s in (("median SD", MED), ("p90 SD", P90), ("mean SD", statistics.mean(sds))):
    n1 = sum(1 for d in ddg_all if abs(d - CUT) <= s)
    n2 = sum(1 for d in ddg_all if abs(d - CUT) <= 2 * s)
    print(f"  {label:>10} = {s:.4f} -> within 1 SD: {n1:>3}   within 2 SD: {n2:>3}")

# the tied band: ranks whose mean is not separable from rank 15's at 3 runs each
for label, s in (("median SD", MED), ("p90 SD", P90)):
    need = 1.96 * s * (2 / 3) ** 0.5
    tied = [i + 1 for i, d in enumerate(ddg_all) if abs(d - CUT) < need]
    print(f"\ntied band around rank 15 at {label} ({s:.4f}); resolvable gap = {need:.4f} kcal/mol")
    if tied:
        print(f"  ranks {min(tied)}-{max(tied)}  ({len(tied)} designs)")
        print(f"  ddg span {ddg_all[min(tied)-1]:+.4f} .. {ddg_all[max(tied)-1]:+.4f}"
              f"  = {ddg_all[max(tied)-1] - ddg_all[min(tied)-1]:.4f} kcal/mol wide")
    else:
        print("  rank 15 is separable from every other design")

print(f"\nadjacent-rank separation needs SD < {min(ddg_all[i+1]-ddg_all[i] for i in range(29)) and 0:.0f}", end="")
gaps = [ddg_all[i + 1] - ddg_all[i] for i in range(29)]
print(f"\n  median adjacent gap in the top 30 = {statistics.median(gaps):.4f}")
print(f"  that gap is resolvable only if SD < {statistics.median(gaps)/(1.96*(2/3)**0.5):.4f}")
print(f"  measured median SD is {MED:.4f} -> adjacent ranks "
      f"{'ARE' if MED < statistics.median(gaps)/(1.96*(2/3)**0.5) else 'are NOT'} separable")
