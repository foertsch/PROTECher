"""Is the merged liability rule right, and does it change what the pool contains?

`design_pipeline/adapters/_gates.py` holds one implementation of the sequence liability
rules. It used to have a narrower twin in `examples/track_a_run.py` - no `D[GSD]`, no `DP`,
`NG` where this one uses `N[GST]`, and a cysteine COUNT comparison instead of positions -
and the union lives in the library now, with the script calling it. This measures the
difference the merge made, on the populations that matter:

  A  the 4,313 single substitutions of the two variable domains: how many each rule flags,
     and where they disagree. A correctness measurement, not a filter.
  B  the BENCHMARK 29-variant direction panel: the merged rule must not fail a measured
     improvement, because a rule that vetoes molecules somebody already measured as better
     is wrong no matter how sound its motifs look.
  C  the 432 design sets: the union must not re-filter the pool.

Reads the chain sequences and the disulfide map off the structure exactly as
`_gates.LiabilityCheck.run` does, so this measures the shipping rule and not a copy.

    python examples/liability_union_check.py [structure.pdb]
"""

from __future__ import annotations

import csv
import re
import sys
from collections import Counter
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT))

from design_pipeline.adapters._gates import LIABILITIES, liability_reasons  # noqa: E402
from design_pipeline.structures import chain_sequence, disulfide_positions  # noqa: E402
from design_pipeline.variants import parse_id  # noqa: E402

STRUCTURE = ROOT / "data/structures/1N8Z_fv_trimmed_Repair.pdb"
SINGLES = ROOT / "runs/negative_controls/wt_argmax_region/scores/68d2f5b741fed874/candidates.csv"
PANEL = ROOT / "runs/negative_controls/direction_set.csv"
SETS = ROOT / "data/track_a_final/candidates.csv"
SHIPPED = ROOT / "runs/rerank/ranking_qualified_v1.csv"

# The narrower copy, as _gates.py's docstring records it. Kept here and nowhere else: it
# exists to be measured against, not to be used.
OLD_LIABILITIES = (
    ("n_glycosylation", r"N[^P][ST]"),
    ("deamidation", r"NG"),
)


def rows(path: Path) -> list[dict]:
    with path.open(newline="") as handle:
        return list(csv.DictReader(handle))


def old_reasons(chains, mutations) -> list[str]:
    """The narrower copy: two motifs, and a cysteine count comparison."""
    found: list[str] = []
    for chain in sorted({c for c, _, _, _ in mutations}):
        residues = list(chains[chain])
        for c, position, _wt, mut in mutations:
            if c == chain:
                residues[position - 1] = mut
        mutant = "".join(residues)
        for label, pattern in OLD_LIABILITIES:
            look = f"(?={pattern})"
            before = {m.start() for m in re.finditer(look, chains[chain])}
            after = {m.start() for m in re.finditer(look, mutant)}
            found += [f"{chain}:{label}@{s + 1}" for s in sorted(after - before)]
        if mutant.count("C") != chains[chain].count("C"):
            found.append(f"{chain}:cysteine_count")
    return sorted(found)


def main(argv: list[str] | None = None) -> int:
    structure = Path((argv or sys.argv[1:])[0]) if (argv or sys.argv[1:]) else STRUCTURE
    singles = rows(SINGLES)
    chains = sorted({r["chain"] for r in singles})
    wild = {c: chain_sequence(structure, c) for c in chains}
    paired = disulfide_positions(structure)

    print(f"structure {structure.name}")
    for c in chains:
        cys = [i + 1 for i, a in enumerate(wild[c]) if a == "C"]
        print(f"  chain {c}: {len(wild[c])} residues; cysteines {cys}; paired {paired.get(c)}")
    print(f"  merged motifs: {[l for l, _ in LIABILITIES]}")
    print(f"  narrower copy: {[l for l, _ in OLD_LIABILITIES]} + cysteine count")

    off = [r for r in singles if wild[r["chain"]][int(r["position"]) - 1] != r["wt"]]
    print(f"\nwt column against the structure: {len(singles) - len(off)}/{len(singles)} agree")
    if off:
        print(f"  MISMATCH at {off[0]} - positions are not 1-based into the chain sequence")
        return 1

    # ------------------------------------------------------------------ A
    print(f"\n=== A. the {len(singles)} single substitutions ===")
    merged, old = set(), set()
    merged_why, merged_only_why = Counter(), Counter()
    copy_only = []
    for r in singles:
        muts = [(r["chain"], int(r["position"]), r["wt"], r["mut"])]
        new_r, old_r = liability_reasons(wild, muts, paired), old_reasons(wild, muts)
        if new_r:
            merged.add(r["variant_id"])
            merged_why.update(x.split("@")[0] for x in new_r)
            if not old_r:
                merged_only_why.update(x.split("@")[0].split(":")[1] for x in new_r)
        if old_r:
            old.add(r["variant_id"])
            if not new_r:
                copy_only.append((r, old_r))

    print(f"merged rule flags : {len(merged)} of {len(singles)}  ({100*len(merged)/len(singles):.1f}%)")
    print(f"narrower copy     : {len(old)} of {len(singles)}")
    print(f"disagreements     : {len(merged ^ old)}  "
          f"(merged only {len(merged - old)}, copy only {len(old - merged)})")
    print("merged rule by reason:")
    for label, n in merged_why.most_common():
        print(f"  {label:<22} {n}")
    print(f"motifs the copy lacked: {dict(merged_only_why)}")

    # The copy catching anything the union misses is the interesting direction.
    print(f"\nflagged by the copy and NOT by the merged rule: {len(copy_only)}")
    for r, why in copy_only:
        at = int(why[0].split("@")[1])
        chain = r["chain"]
        print(f"  {r['variant_id']:<12} {chain}{at}:{wild[chain][at-1]}{wild[chain][at]}"
              f" -> {wild[chain][at-1]}{r['mut']}   copy: {why}")
    if copy_only:
        print("  Cause: `new_motifs` is occurrence-wise by START position, and the merged")
        print("  deamidation class N[GST] already matched before the substitution. NT -> NG")
        print("  at the same asparagine is a real jump in deamidation risk and the merged")
        print("  rule reports nothing, because the match position did not move.")

    # ------------------------------------------------------------------ B
    panel = rows(PANEL)
    up = [r for r in panel if r["direction"] == "1"]
    print(f"\n=== B. the {len(panel)}-variant direction panel "
          f"({len(up)} improvements + {len(panel)-len(up)} losses) ===")
    failed = [(r, liability_reasons(wild, list(parse_id(r["variant_id"])), paired)) for r in panel]
    failed = [(r, why) for r, why in failed if why]
    print(f"failed the merged rule: {len(failed)} of {len(panel)}")
    for r, why in failed:
        kind = "IMPROVEMENT" if r["direction"] == "1" else "loss"
        print(f"  {r['variant_id']:<40} [{kind}] label {r['label']}  {why}")
    bad = [r for r, _ in failed if r["direction"] == "1"]
    print(f"measured improvements vetoed: {len(bad)} of {len(up)}"
          + (f"  {[r['variant_id'] for r in bad]}" if bad else "  (none - the rule is safe here)"))

    # ------------------------------------------------------------------ C
    sets = rows(SETS)
    hits = [r for r in sets if liability_reasons(wild, list(parse_id(r["variant_id"])), paired)]
    shipped = {r["variant_id"] for r in rows(SHIPPED)}
    bad_shipped = [v for v in shipped if liability_reasons(wild, list(parse_id(v)), paired)]
    print(f"\n=== C. does the union re-filter the pool? ===")
    print(f"the {len(sets)} design sets: {len(hits)} fail")
    print(f"the shipped {len(shipped)}: {len(bad_shipped)} fail")
    if not hits:
        print("  -> no. The merge is a correctness fix, not a change to what was selected.")
    return 0


if __name__ == "__main__":
    sys.exit(main())
