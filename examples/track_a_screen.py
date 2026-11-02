"""Screen scFv candidates with Boltz-2, one process per batch, and score every fold with ipTM, pLDDT and ipSAE.

Meant for a time-boxed screen of many candidates. `examples/track_a_run.py fold` starts one `boltz predict` per design, and
each one reloads the weights (about 55 s a fold on an L40S). This folds a whole slice of candidates in one process with
`msa_cache.fold_batch`, the way cloud/modal_fold.py does (5 folds in 95 s on Modal). The alignments are the parent's,
from the cache, so no server is called. One seed per candidate catches collapses and cannot rank candidates that are
within about 0.1 ipSAE of each other.

    python examples/track_a_screen.py fold  --fasta candidates.fasta --subset 0-89 --seeds 1 --out screen/fold_a
    python examples/track_a_screen.py score --out screen/fold_a        # also works on a run that was stopped part way

candidates.fasta holds scFv sequences only (`>name` then the 242 residues); the antigen fragment is added here. `fold`
scores at the end, and `score` reads whatever Boltz has written under --out, so a stopped batch still gives results.
"""

from __future__ import annotations

import argparse
import csv
import json
import re
import statistics
import sys
import time
from pathlib import Path
from typing import Sequence

FIELDS = ["name", "seed", "iptm", "plddt", "ipsae_max", "ipsae_min", "ipsae_ab", "ipsae_ba"]


def read_fasta(path: Path) -> list[tuple[str, str]]:
    records, name, parts = [], None, []
    for line in path.read_text().splitlines():
        if line.startswith(">"):
            if name is not None:
                records.append((name, "".join(parts)))
            name, parts = line[1:].strip(), []
        elif line.strip():
            parts.append(line.strip())
    if name is not None:
        records.append((name, "".join(parts)))
    return records


def cif_chains(path: Path):
    """The chain id of each CA atom of a Boltz mmCIF, in file order: one entry per residue, as the PAE rows are."""
    from biotite.structure.io import pdbx

    atoms = pdbx.get_structure(pdbx.CIFFile.read(str(path)), model=1)
    return atoms[atoms.atom_name == "CA"].chain_id


def score_folds(out: Path, pae_cutoff: float = 15.0) -> list[dict]:
    """One row per fold found under `out`: ipTM and pLDDT from Boltz, ipSAE (d0 per residue) from the PAE matrix."""
    import numpy as np

    from design_pipeline.scoring import ipsae

    rows = []
    for pae_path in sorted(out.rglob("pae_*_model_0.npz")):
        name = pae_path.name[len("pae_"):-len("_model_0.npz")]
        folder = pae_path.parent
        cif, confidence = folder / f"{name}_model_0.cif", folder / f"confidence_{name}_model_0.json"
        if not cif.is_file() or not confidence.is_file():
            continue
        seed = re.search(r"seed(\d+)", str(pae_path.relative_to(out)))
        pae = np.load(pae_path)["pae"]
        chains = cif_chains(cif)
        summary = json.loads(confidence.read_text())
        scores = ipsae(pae, chains, pae_cutoff=pae_cutoff)
        rows.append({
            "name": name, "seed": int(seed.group(1)) if seed else 0,
            "iptm": round(summary["iptm"], 4), "plddt": round(100.0 * summary["complex_plddt"], 2),
            "ipsae_max": round(scores["max"], 4), "ipsae_min": round(scores["min"], 4),
            "ipsae_ab": round(scores["ab"], 4), "ipsae_ba": round(scores["ba"], 4),
        })
    return rows


def write_and_report(rows: list[dict], out: Path) -> None:
    out.mkdir(parents=True, exist_ok=True)
    with (out / "screen_scores.csv").open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=FIELDS)
        writer.writeheader()
        writer.writerows(rows)
    by: dict[str, list[dict]] = {}
    for row in rows:
        by.setdefault(row["name"], []).append(row)
    ranked = sorted(by, key=lambda n: -statistics.mean(r["ipsae_max"] for r in by[n]))
    print(f"{len(rows)} folds of {len(by)} candidates -> {out / 'screen_scores.csv'}")
    print(f"{'name':<8}{'ipSAE max':>10}{'ipSAE min':>10}{'ipTM':>7}{'pLDDT':>7}  seeds")
    shown = ranked if len(ranked) <= 24 else ranked[:12] + ["..."] + ranked[-12:]
    for name in shown:
        if name == "...":
            print("...")
            continue
        r = by[name]
        print(f"{name:<8}{statistics.mean(x['ipsae_max'] for x in r):10.3f}{statistics.mean(x['ipsae_min'] for x in r):10.3f}"
              f"{statistics.mean(x['iptm'] for x in r):7.3f}{statistics.mean(x['plddt'] for x in r):7.1f}  {len(r)}")


def fold(args: argparse.Namespace) -> int:
    from design_pipeline import Design
    from msa_cache import MsaCache
    from msa_cache.batch import fold_batch

    out = Path(args.out)
    candidate = json.loads(Path(args.candidate).read_text())
    start, end = (int(x) for x in args.antigen_region.split("-"))
    antigen = candidate["antigen"][start - 1:end]
    records = read_fasta(Path(args.fasta))
    first, last = (int(x) for x in args.subset.split("-")) if args.subset else (0, len(records) - 1)
    records = records[first:last + 1]
    inputs = out / "designs"
    inputs.mkdir(parents=True, exist_ok=True)
    designs = []
    for name, scfv in records:
        path = inputs / f"{name}.fa"
        path.write_text(f">{name}\n{scfv}:{antigen}\n")
        designs.append(Design(id=name, path=path, parent_id=None, producer="track_a_screen"))
    print(f"folding {len(designs)} candidates ({designs[0].id} to {designs[-1].id}), antigen {start}-{end}, "
          f"seeds {args.seed_start} to {args.seed_start + args.seeds - 1}, one process per seed", flush=True)
    cache = MsaCache(args.msa_cache)
    seconds = {}
    for seed in range(args.seed_start, args.seed_start + args.seeds):
        began = time.perf_counter()
        fold_batch(cache, designs, out / f"seed{seed}", seed=seed)
        seconds[f"seed{seed}"] = round(time.perf_counter() - began, 1)
        print(f"seed {seed}: {seconds[f'seed{seed}']} s for {len(designs)} folds "
              f"({seconds[f'seed{seed}'] / len(designs):.1f} s each)", flush=True)
    (out / "fold_timing.json").write_text(json.dumps({"seconds": seconds, "candidates": len(designs)}, indent=1))
    write_and_report(score_folds(out), out)
    return 0


def score(args: argparse.Namespace) -> int:
    write_and_report(score_folds(Path(args.out)), Path(args.out))
    return 0


def main(argv: Sequence[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    sub = parser.add_subparsers(dest="command", required=True)
    f = sub.add_parser("fold", help="fold a slice of candidates, one boltz process per seed, then score")
    f.add_argument("--fasta", required=True)
    f.add_argument("--out", required=True)
    f.add_argument("--subset", help="FIRST-LAST record numbers, 0-based and inclusive (default all)")
    f.add_argument("--seeds", type=int, default=1)
    f.add_argument("--seed-start", type=int, default=0)
    f.add_argument("--antigen-region", default="480-607", help="START-END of the 630-residue HER2 ectodomain, 1-based")
    f.add_argument("--candidate", default="/scratch/foertsch/protecher_cache/msa_candidates/trastuzumab_her2/candidate.json")
    f.add_argument("--msa-cache", default="/scratch/foertsch/protecher_cache/msa_candidates/trastuzumab_her2")
    f.set_defaults(run=fold)
    s = sub.add_parser("score", help="ipTM, pLDDT and ipSAE of every fold found under --out")
    s.add_argument("--out", required=True)
    s.set_defaults(run=score)
    args = parser.parse_args(argv)
    return args.run(args)


if __name__ == "__main__":
    sys.exit(main())
