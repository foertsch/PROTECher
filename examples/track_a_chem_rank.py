"""Score a Boltz-2 candidate screen with ipSAE and the chemistry terms of PR 24, then rank it.

Input is what `examples/track_a_screen.py fold` (cluster) or `cloud/modal_fold.py::fold_dir` (Modal) writes: batches of
folds in which the parent was folded as the record `parent`, so every candidate is compared with the parent of the same batch
and seed. A batch is one `boltz_results_inputs` folder; its seed is the `seed<k>` folder above it (0 when there is none).

    python examples/track_a_chem_rank.py score --root screen/modal_s1 --out rows_modal_s1.csv --workers 8
    python examples/track_a_chem_rank.py rank --rows rows_*.csv --candidates candidates.csv --out ranking.csv

`score` is the CPU-heavy part and runs one process per fold. `rank` follows docs/chem_score_spec.md: candidates whose mean
ipSAE (max direction) is at least the parent's mean minus two SD are ranked by chem = z(contact) + z(coulomb), over the seeds
whose pose passes the epitope filter; the rest follow, ordered by ipSAE. A candidate with no on-site seed has no chem score
and is ranked by ipSAE. Nothing is vetoed.
"""

from __future__ import annotations

import argparse
import csv
import glob
import json
import math
import re
import statistics
import sys
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path
from typing import Sequence

sys.path.insert(0, str(Path(__file__).resolve().parent))

ROW_FIELDS = ["name", "seed", "batch", "eligible", "contact", "coulomb", "epitope_overlap", "n_added_cation", "paratope_net_charge",
              "iptm", "plddt", "ipsae_max", "ipsae_min"]
PARENT = "parent"


def _fold_files(pae_path: Path) -> tuple[str, Path, Path]:
    name = pae_path.name[len("pae_"):-len("_model_0.npz")]
    return name, pae_path.parent / f"{name}_model_0.cif", pae_path.parent / f"confidence_{name}_model_0.json"


def _score_one(job: tuple) -> dict | None:
    pae_path, root, parent_cif = job
    import numpy as np

    from design_pipeline.scoring import charge_contact_report, ipsae
    from track_a_screen import cif_chains

    name, cif, confidence = _fold_files(Path(pae_path))
    if not cif.is_file() or not confidence.is_file():
        return None
    rel = str(Path(pae_path).relative_to(root))
    seed = re.search(r"seed(\d+)", rel)
    batch = str(Path(pae_path).parents[2].relative_to(root))
    summary = json.loads(confidence.read_text())
    scores = ipsae(np.load(pae_path)["pae"], cif_chains(cif))
    row = {"name": name, "seed": int(seed.group(1)) if seed else 0, "batch": batch, "iptm": round(summary["iptm"], 4),
           "plddt": round(100.0 * summary["complex_plddt"], 2), "ipsae_max": round(scores["max"], 4), "ipsae_min": round(scores["min"], 4)}
    if name != PARENT and parent_cif is not None:
        chem = charge_contact_report(cif, Path(parent_cif))
        if chem is not None:
            row.update({k: round(float(chem[k]), 4) for k in ("eligible", "contact", "coulomb", "epitope_overlap", "n_added_cation", "paratope_net_charge")})
    return row


def score(args: argparse.Namespace) -> int:
    root = Path(args.root).resolve()
    jobs = []
    for pae in sorted(root.rglob("pae_*_model_0.npz")):
        batch_predictions = pae.parents[1]
        parent_cif = batch_predictions / PARENT / f"{PARENT}_model_0.cif"
        jobs.append((str(pae), str(root), str(parent_cif) if parent_cif.is_file() else None))
    print(f"{len(jobs)} folds under {root}", flush=True)
    with ProcessPoolExecutor(max_workers=args.workers) as pool:
        rows = [r for r in pool.map(_score_one, jobs, chunksize=2) if r is not None]
    with open(args.out, "w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=ROW_FIELDS)
        writer.writeheader()
        writer.writerows(rows)
    print(f"{len(rows)} rows -> {args.out}")
    return 0


def _mean(values: Sequence[float]) -> float:
    values = [v for v in values if v is not None and not math.isnan(v)]
    return statistics.mean(values) if values else math.nan


def _z(values: dict[str, float]) -> dict[str, float]:
    mean, sd = statistics.mean(values.values()), statistics.pstdev(values.values())
    return {k: (v - mean) / (sd if sd else 1.0) for k, v in values.items()}


def rank_rows(rows: list[dict], band_sd: float = 2.0) -> tuple[list[dict], dict]:
    """One dict per candidate, best first, and the parent's mean and SD of ipSAE (max) that set the band."""
    def num(row, key):
        try:
            return float(row[key])
        except (KeyError, TypeError, ValueError):
            return math.nan

    parent = [num(r, "ipsae_max") for r in rows if r["name"] == PARENT]
    if len(parent) < 2:
        raise SystemExit("need the parent folded at two seeds or more to set the ipSAE band")
    parent_mean, parent_sd = statistics.mean(parent), statistics.stdev(parent)
    floor = parent_mean - band_sd * parent_sd
    by: dict[str, list[dict]] = {}
    for r in rows:
        if r["name"] != PARENT:
            by.setdefault(r["name"], []).append(r)
    out = []
    for name, rs in by.items():
        elig = [r for r in rs if num(r, "eligible") >= 1.0]
        out.append({"name": name, "seeds": len(rs), "ipsae_max": _mean([num(r, "ipsae_max") for r in rs]),
                    "ipsae_min": _mean([num(r, "ipsae_min") for r in rs]), "iptm": _mean([num(r, "iptm") for r in rs]),
                    "n_on_site": len(elig), "contact": _mean([num(r, "contact") for r in elig]),
                    "coulomb": _mean([num(r, "coulomb") for r in elig]), "para_net": _mean([num(r, "paratope_net_charge") for r in rs])})
    for o in out:
        o["in_band"] = o["ipsae_max"] >= floor
    scorable = {o["name"]: o for o in out if o["in_band"] and o["n_on_site"] > 0 and not math.isnan(o["contact"])}
    if len(scorable) >= 2:
        zc, zq = _z({n: o["contact"] for n, o in scorable.items()}), _z({n: o["coulomb"] for n, o in scorable.items()})
        for n, o in scorable.items():
            o["chem"] = zc[n] + zq[n]
    for o in out:
        o.setdefault("chem", math.nan)
    out.sort(key=lambda o: (0 if not math.isnan(o["chem"]) else 1 if o["in_band"] else 2,
                            -o["chem"] if not math.isnan(o["chem"]) else -o["ipsae_max"]))
    for i, o in enumerate(out, 1):
        o["rank"] = i
    return out, {"parent_mean": parent_mean, "parent_sd": parent_sd, "floor": floor, "parent_folds": len(parent)}


def rank_weighted(rows: list[dict], chem_weight: float, band_sd: float = 1.0) -> tuple[list[dict], dict]:
    """Rank by score = z(ipSAE) + chem_weight * chem, over all candidates, with no band cut.

    ipSAE is the mean over seeds (max direction). chem is z(contact) + z(coulomb) over the candidates that have an on-site
    seed; one without any gets 0 (neutral, never a penalty). A small weight keeps chemistry in the ranking without
    letting it override ipSAE. `in_band` is still reported against the parent's mean minus `band_sd` SD."""
    ranked, band = rank_rows(rows, band_sd)
    ip = _z({o["name"]: o["ipsae_max"] for o in ranked})
    have = {o["name"]: o for o in ranked if o["n_on_site"] > 0 and not math.isnan(o["contact"])}
    chem = {}
    if len(have) >= 2:
        zc, zq = _z({n: o["contact"] for n, o in have.items()}), _z({n: o["coulomb"] for n, o in have.items()})
        chem = {n: zc[n] + zq[n] for n in have}
    for o in ranked:
        o["chem_z"] = chem.get(o["name"], 0.0)
        o["score"] = ip[o["name"]] + chem_weight * o["chem_z"]
    ranked.sort(key=lambda o: -o["score"])
    for i, o in enumerate(ranked, 1):
        o["rank"] = i
    return ranked, band


def rank(args: argparse.Namespace) -> int:
    rows = []
    for pattern in args.rows:
        for path in glob.glob(pattern):
            with open(path, newline="") as handle:
                rows += list(csv.DictReader(handle))
    meta = {}
    if args.candidates:
        with open(args.candidates, newline="") as handle:
            meta = {r["name"]: r for r in csv.DictReader(handle)}
    ranked, band = rank_weighted(rows, args.chem_weight) if args.chem_weight is not None else rank_rows(rows)
    print(f"parent ipSAE (max) {band['parent_mean']:.3f}, SD {band['parent_sd']:.3f} over {band['parent_folds']} folds; band floor {band['floor']:.3f}")
    print(f"{len(ranked)} candidates; in band {sum(o['in_band'] for o in ranked)}; chem-scored {sum(not math.isnan(o['chem']) for o in ranked)}")
    fields = ["rank", "name", "source", "variant_id", "seeds", "ipsae_max", "ipsae_min", "iptm", "in_band", "n_on_site", "contact", "coulomb", "para_net", "chem", "chem_z", "score"]
    with open(args.out, "w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields)
        writer.writeheader()
        for o in ranked:
            m = meta.get(o["name"], {})
            writer.writerow({**{k: (round(v, 4) if isinstance(v, float) else v) for k, v in o.items() if k in fields},
                             "source": m.get("source", ""), "variant_id": m.get("variant_id", "")})
    print(f"{'#':>3} {'name':<6}{'ipSAE':>7}{'seeds':>6}{'on-site':>8}{'contact':>9}{'coulomb':>9}{'chem':>7}  variant")
    for o in ranked[: args.top]:
        m = meta.get(o["name"], {})
        print(f"{o['rank']:>3} {o['name']:<6}{o['ipsae_max']:7.3f}{o['seeds']:6d}{o['n_on_site']:8d}{o['contact']:9.3f}{o['coulomb']:9.3f}{o['chem']:7.2f}  {m.get('variant_id', '')[:70]}")
    return 0


def main(argv: Sequence[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    sub = parser.add_subparsers(dest="command", required=True)
    s = sub.add_parser("score", help="ipSAE and chemistry terms for every fold under --root")
    s.add_argument("--root", required=True)
    s.add_argument("--out", required=True)
    s.add_argument("--workers", type=int, default=4)
    s.set_defaults(run=score)
    r = sub.add_parser("rank", help="aggregate rows files over seeds and rank the candidates")
    r.add_argument("--rows", nargs="+", required=True, help="rows CSVs from `score` (globs allowed)")
    r.add_argument("--candidates", help="candidates.csv, to show the variant ids")
    r.add_argument("--out", required=True)
    r.add_argument("--top", type=int, default=25)
    r.add_argument("--chem-weight", type=float, help="rank by z(ipSAE) + W * chemistry over all candidates instead of the band rule (0.15 keeps chemistry in at a small weight)")
    r.set_defaults(run=rank)
    args = parser.parse_args(argv)
    return args.run(args)


if __name__ == "__main__":
    sys.exit(main())
