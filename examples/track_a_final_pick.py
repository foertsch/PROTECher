"""Pick the final designs from the scored Boltz-2 screen: the best by ipSAE plus a little chemistry, then the best by chemistry.

Input is what `examples/track_a_chem_rank.py score` wrote (rows CSVs, one row per fold, the parent folded in every batch) and
the candidates.csv the screen was made from. A candidate fold that appears in several rows files (the same candidate and seed
folded twice) counts once, the first one read; every fold of the parent counts.

    python examples/track_a_final_pick.py --rows rows_*.csv --candidates candidates.csv --out final15.csv

The rule, in order:
  1. The first `--top` (10) candidates by score = z(ipSAE) + 0.15 * chemistry (`track_a_chem_rank.rank_weighted`). Chemistry
     is z(contact) + z(coulomb) over the seeds on the epitope; a candidate with no such seed gets 0.
  2. The next `--chem-picks` (5): the highest chemistry among the rest whose mean ipSAE is at least the parent's mean minus
     two SD, and which differ from every earlier pick by at least two substitutions. Chemistry here is the lead criterion, so
     the band keeps out candidates whose complex has collapsed.
Writes one row per pick with its group, and the parent's mean and SD that set the band.
"""

from __future__ import annotations

import argparse
import csv
import glob
import json
import math
import sys
from pathlib import Path
from typing import Sequence

sys.path.insert(0, str(Path(__file__).resolve().parent))

from track_a_chem_rank import rank_weighted  # noqa: E402

GROUP_SCORE = "ipSAE + chemistry"
GROUP_CHEM = "chemistry"
FIELDS = ["design", "name", "group", "variant_id", "seeds", "ipsae_max", "ipsae_min", "iptm", "chem_z", "contact", "coulomb", "score"]


def read_rows(patterns: Sequence[str]) -> list[dict]:
    """Rows of every file the patterns match, each candidate fold (name, seed) once. The parent is folded in every batch and
    each of those folds sets the band, so every one is kept."""
    rows, seen = [], set()
    for pattern in patterns:
        for path in sorted(glob.glob(pattern)):
            with open(path, newline="") as handle:
                for row in csv.DictReader(handle):
                    key = (row["name"], row["seed"])
                    if row["name"] == "parent" or key not in seen:
                        seen.add(key)
                        rows.append(row)
    return rows


def substitutions(variant_id: str) -> set[str]:
    return {s for s in variant_id.split("+") if s}


def pick(rows: list[dict], variants: dict[str, str], chem_weight: float = 0.15, top: int = 10, chem_picks: int = 5,
         band_sd: float = 2.0, min_different: int = 2) -> tuple[list[dict], dict]:
    """The picks in order (score group first), and the parent's band: mean, SD and the floor the chemistry picks must clear."""
    ranked, band = rank_weighted(rows, chem_weight, band_sd)
    floor = band["parent_mean"] - band_sd * band["parent_sd"]
    for o in ranked:
        o["variant_id"] = variants.get(o["name"], "")
    picks = [dict(o, group=GROUP_SCORE) for o in ranked[:top]]
    taken = [substitutions(p["variant_id"]) for p in picks]
    rest = sorted((o for o in ranked[top:] if o["ipsae_max"] >= floor), key=lambda o: -o["chem_z"])
    for o in rest:
        if sum(p["group"] == GROUP_CHEM for p in picks) == chem_picks:
            break
        subs = substitutions(o["variant_id"])
        if all(len(subs - t) >= min_different for t in taken):
            picks.append(dict(o, group=GROUP_CHEM))
            taken.append(subs)
    return picks, {**band, "floor": floor}


def main(argv: Sequence[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    parser.add_argument("--rows", nargs="+", required=True, help="rows CSVs from track_a_chem_rank.py score (globs allowed)")
    parser.add_argument("--candidates", required=True, help="candidates.csv with name and variant_id")
    parser.add_argument("--out", required=True, type=Path)
    parser.add_argument("--chem-weight", type=float, default=0.15)
    parser.add_argument("--top", type=int, default=10)
    parser.add_argument("--chem-picks", type=int, default=5)
    args = parser.parse_args(argv)

    with open(args.candidates, newline="") as handle:
        variants = {r["name"]: r["variant_id"] for r in csv.DictReader(handle)}
    picks, band = pick(read_rows(args.rows), variants, args.chem_weight, args.top, args.chem_picks)
    if len(picks) < args.top + args.chem_picks:
        print(f"only {len(picks)} picks: fewer candidates clear the band ({band['floor']:.3f}) with two new substitutions than asked for", file=sys.stderr)
    with open(args.out, "w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=FIELDS)
        writer.writeheader()
        for i, p in enumerate(picks, 1):
            writer.writerow({"design": f"design_{i:02d}", **{k: (round(v, 4) if isinstance(v, float) else v) for k, v in p.items() if k in FIELDS}})
    args.out.with_suffix(".parent.json").write_text(json.dumps({k: band[k] for k in ("parent_mean", "parent_sd", "parent_folds", "floor")}, indent=1))
    print(f"parent ipSAE (max) {band['parent_mean']:.3f}, SD {band['parent_sd']:.3f} over {band['parent_folds']} folds; chemistry picks need >= {band['floor']:.3f}")
    for i, p in enumerate(picks, 1):
        chem = "" if math.isnan(p["chem_z"]) else f"{p['chem_z']:6.2f}"
        print(f"{i:>3} {p['name']:<6}{p['group']:<20}{p['ipsae_max']:6.3f}{chem:>7}  {p['variant_id']}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
