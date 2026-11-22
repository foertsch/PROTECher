"""ipSAE through the negative-control panel, as the seven scorers were.

`03_negative_controls.ipynb` judges a scorer on labelled negatives. ipSAE never
went through it: it is not a `score_variants` adapter, because it needs a fold
rather than a sequence or a side chain. It is, however, the column that carries
weight 1.0 in `examples/track_a_chem_rank.py`'s composite, and the seed-variance
measurement on `data/track_a_final/rows_*.csv` says most of that column is noise
(ICC 0.30 over 3 seeds, so one seed against another ranks at rho ~0.31). A
column that decides the shortlist has to clear the same bar as the scorers it
outranks.

So: fold every design in a panel, score each fold's ipSAE, and write one
`variant_id,score_ipsae` table the notebook reads as another scorer column -
which is what keeps the AUROC, its bootstrap interval and the permutation p
identical to the other seven rather than a second implementation of them.

**The aggregate is the mean across seeds, not the max.** `ipsae_max`/`ipsae_min`
in `ranking.csv` are the two chain directions (`ab` aligns on the scFv and
scores the antigen, `ba` the reverse), averaged over seeds; this script keeps
those as reported columns and scores on the per-seed mean of the larger
direction. The max over k folds is the largest of k draws and so is biased
upward by roughly 0.85 SD at k = 3, and that bias scales with the per-design
noise - which is exactly what differs between a design folded well and one
folded badly, so ranking on it rewards noise.

**The parent sequences come from the alignment cache, not from the crystal.**
`MsaCache.parent_chains` holds the two chains the alignments were actually built
for, and those are the sequences every fold has to use or the ipSAE is not
comparable to `ranking.csv`'s column. The crystal is the fallback and the check:
VH is chain B's variable domain, VL is chain A's, joined by the (G4S)3 linker,
and `scfv.verify_parent` confirms the scFv against 1N8Z either way. The crystal
cannot supply the antigen on its own - chain C is missing 10 of the 128 residues
of the 480-607 window, so a fragment cut from the coordinates is 118 residues
and a fold of it would answer a slightly different question.

    # on the box, with the parent alignments in place
    python examples/ipsae_panel.py fold  --panel direction --msa-cache data/msa_candidates/trastuzumab_her2
    python examples/ipsae_panel.py score --panel direction

`fold` needs a GPU and the Boltz environment; `score` reads whatever Boltz has
written, so a run that was stopped part way still gives a table.
"""

from __future__ import annotations

import argparse
import csv
import json
import re
import statistics
import time
from pathlib import Path

ANTIGEN_FIRST, ANTIGEN_LAST = 480, 607
"""The HER2 fragment, by chain C's own labels in 1N8Z: the part near the
paratope, the same window `examples/trim_fv_for_foldx.py` cuts."""

VH_LENGTH, VL_LENGTH = 120, 107
HEAVY, LIGHT, ANTIGEN = "B", "A", "C"

FOLD_FIELDS = ["variant_id", "name", "seed", "iptm", "plddt",
               "ipsae_ab", "ipsae_ba", "ipsae_max", "ipsae_min"]


THREE_TO_ONE = {
    "ALA": "A", "ARG": "R", "ASN": "N", "ASP": "D", "CYS": "C", "GLN": "Q",
    "GLU": "E", "GLY": "G", "HIS": "H", "ILE": "I", "LEU": "L", "LYS": "K",
    "MET": "M", "PHE": "F", "PRO": "P", "SER": "S", "THR": "T", "TRP": "W",
    "TYR": "Y", "VAL": "V",
}


def seqres(pdb: Path, chain: str) -> str:
    """A chain's full sequence from SEQRES, including residues with no coordinates.

    What the crystal was *made of*, rather than what diffracted. Chain C of 1N8Z
    is 607 residues of SEQRES and 581 of them are observed, so the antigen
    window has to come from here: a fragment cut from the ATOM records drops the
    26 unresolved residues and comes out short, and folding a sequence with
    internal deletions answers a question nobody asked. On 1N8Z the labels run
    1..607 on this chain, so SEQRES position i is label i; `parent_from_crystal`
    checks that against the coordinates rather than assuming it.
    """
    residues: list[str] = []
    for line in pdb.read_text(encoding="utf-8", errors="replace").splitlines():
        if line.startswith("SEQRES") and len(line) > 19 and line[11] == chain:
            residues.extend(line[19:].split())
    return "".join(THREE_TO_ONE.get(residue, "X") for residue in residues)


def parent_from_crystal(pdb: Path):
    """(scfv, antigen) rebuilt from the structure file. The fallback and the check.

    The scFv is the crystal's two variable domains joined by the linker. The
    antigen is the 480-607 window of chain C's SEQRES, and the window is only
    meaningful if that chain's labels are 1-based and contiguous over it, which
    is asserted here against the observed residues.
    """
    from design_pipeline import structures
    from design_pipeline.scfv import LINKER

    heavy = structures.chain_sequence(pdb, HEAVY)[:VH_LENGTH]
    light = structures.chain_sequence(pdb, LIGHT)[:VL_LENGTH]

    full = seqres(pdb, ANTIGEN)
    antigen = full[ANTIGEN_FIRST - 1 : ANTIGEN_LAST]
    if len(antigen) != ANTIGEN_LAST - ANTIGEN_FIRST + 1:
        raise SystemExit(
            f"chain {ANTIGEN} of {pdb} has {len(full)} SEQRES residues, too few "
            f"for the {ANTIGEN_FIRST}-{ANTIGEN_LAST} window"
        )
    # SEQRES carries no numbering, so this is where the assumption that position
    # i is label i gets tested: every observed residue of the window has to sit
    # at the same place in SEQRES as its own label says.
    labels = structures.chain_residue_ids(pdb, ANTIGEN)
    observed = structures.chain_sequence(pdb, ANTIGEN)
    for residue, label in zip(observed, labels):
        if not label.isdigit():
            continue
        index = int(label)
        if ANTIGEN_FIRST <= index <= ANTIGEN_LAST and full[index - 1] != residue:
            raise SystemExit(
                f"chain {ANTIGEN} label {label} holds {residue} but SEQRES "
                f"position {index} holds {full[index - 1]}: the numbering is "
                f"not 1-based over this chain, so the window is not what it "
                f"says it is"
            )
    return heavy + LINKER + light, antigen


def parent_scfv(pdb: Path, cache=None):
    """(layout, antigen, source) for the panel's folds.

    The alignment cache is authoritative when it has a parent: those are the
    sequences its alignments were built against, and `MsaCache` would otherwise
    be matching a chain it has no columns for. The crystal is used only when the
    cache has no parent, and `verify_parent` runs against 1N8Z in both cases, so
    a cache holding some other antibody's parent is an error here rather than a
    silent fold of the wrong molecule.
    """
    from design_pipeline.scfv import layout_from_parent, verify_parent

    crystal_scfv, crystal_antigen = parent_from_crystal(pdb)
    scfv, antigen, source = crystal_scfv, crystal_antigen, "the crystal"
    if cache is not None and cache.has_parent():
        chains = list(cache.parent_chains)
        if len(chains) != 2:
            raise SystemExit(
                f"the alignment cache holds {len(chains)} parent chain(s); this "
                f"panel folds an scFv and one antigen"
            )
        scfv, antigen = chains
        source = "the alignment cache"

    layout = layout_from_parent(scfv, first=HEAVY, second=LIGHT)
    verify_parent(layout, pdb)
    if not antigen:
        raise SystemExit(
            f"no antigen sequence: the cache has no parent and chain {ANTIGEN} "
            f"of {pdb} has no residues labelled {ANTIGEN_FIRST}-{ANTIGEN_LAST}"
        )
    if cache is not None and cache.has_parent() and antigen != crystal_antigen:
        print(
            f"note: the cache's antigen is {len(antigen)} residues and the "
            f"crystal's {ANTIGEN_FIRST}-{ANTIGEN_LAST} window resolves "
            f"{len(crystal_antigen)}; folding the cache's, which is what its "
            f"alignment covers"
        )
    return layout, antigen, source


def parent_fasta(args: argparse.Namespace) -> int:
    """Write the parent FASTA the alignments are built from, and nothing else.

    Separate from `fold` because building the alignments is the one step that
    talks to a server: `msa_cache.build_parent` folds this once with
    `use_msa_server=True` and keeps what comes back, and every variant alignment
    after that is derived from it locally. So this prints exactly what would be
    sent, for a decision to be made on the sequences rather than on trust.
    """
    from design_pipeline.scfv import layout_from_parent, verify_parent

    pdb = Path(args.pdb)
    scfv, antigen = parent_from_crystal(pdb)
    verify_parent(layout_from_parent(scfv, first=HEAVY, second=LIGHT), pdb)
    path = Path(args.out) / "parent.fa"
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(f">parent\n{scfv}:{antigen}\n", encoding="utf-8")
    print(f"parent scFv: {len(scfv)} residues (VH {VH_LENGTH} + linker + VL {VL_LENGTH}),"
          f" verified against {pdb}")
    print(f"antigen:     {len(antigen)} residues, chain {ANTIGEN} SEQRES "
          f"{ANTIGEN_FIRST}-{ANTIGEN_LAST}")
    print(f"\n-> {path}")
    print("\nBoth sequences are public: the scFv is the two variable domains of this")
    print("crystal and the antigen is its own SEQRES. No design is in this file.")
    print("To build the alignments from it (one call to the public server, then")
    print("nothing further leaves the machine):")
    print("\n    python -c \"from msa_cache import MsaCache; from msa_cache.cache import "
          "build_parent; import pathlib;"
          f" build_parent(MsaCache('{args.msa_cache}'), "
          f"pathlib.Path('{path}').read_text().split(chr(10))[1].split(':'), "
          f"work_dir='{Path(args.out) / 'parent_fold'}')\"")
    return 0


def read_panel(panel: str, out_dir: Path) -> list[tuple[str, int]]:
    """[(variant_id, label)] for a panel, from the table the notebook wrote.

    `direction` is Part 1's labelled set (1 = measured as improved). `chinery`
    is Part 3's clean split (1 = a fitted KD, 0 = `I.C.`), which Part 3 shows is
    a clean-kinetics label and not an affinity one - so the two are never pooled
    and the notebook reports them in separate columns.
    """
    # The clean split only. `chinery_set.csv` carries all 687 wells, and the
    # 434 `negative_caveated` are the `N/A` ones the notebook declines to treat
    # as demonstrated non-binders - folding them would buy 434 x 3 folds of
    # label noise pointing in an unknown direction.
    keep = {"positive", "negative_clean"}
    if panel == "direction":
        path = out_dir / "direction_set.csv"
        key, label, where = "variant_id", "direction", None
    elif panel == "chinery":
        path = out_dir / "chinery" / "chinery_set.csv"
        key, label, where = "variant_id", "binary", ("label_class", keep)
    else:
        raise SystemExit(f"unknown panel {panel!r}")
    if not path.is_file():
        raise SystemExit(
            f"{path} does not exist. Run the notebook first: it writes the "
            f"panel's labels, and taking them from the benchmark CSVs here "
            f"instead would be a second definition of the same set."
        )
    rows = list(csv.DictReader(path.open(encoding="utf-8")))
    if where is not None:
        column, allowed = where
        rows = [row for row in rows if row[column] in allowed]
    seen: dict[str, int] = {}
    for row in rows:
        variant = row[key]
        value = int(float(row[label]))
        if variant in seen and seen[variant] != value:
            raise SystemExit(
                f"{variant} appears with both labels in {path}; one design "
                f"gets one fold and so one label"
            )
        seen[variant] = value
    return sorted(seen.items())


def fold(args: argparse.Namespace) -> int:
    from design_pipeline import Design
    from design_pipeline.scfv import apply_id
    from msa_cache import MsaCache
    from msa_cache.batch import fold_batch

    out = Path(args.out)
    cache = MsaCache(args.msa_cache)
    if not cache.has_parent():
        raise SystemExit(
            f"{args.msa_cache} has no parent alignment, so every fold would be "
            f"single-sequence and its ipSAE would not be comparable to the "
            f"column in ranking.csv. Copy the parent alignments in first."
        )
    layout, antigen, source = parent_scfv(Path(args.pdb), cache)
    panel = read_panel(args.panel, Path(args.out_root))
    print(f"parent scFv {len(layout.scfv)} residues and antigen {len(antigen)}, "
          f"both from {source}; verified against {args.pdb}")
    print(f"{len(panel)} designs in the {args.panel} panel")

    designs, names = [], {}
    inputs = out / "designs"
    inputs.mkdir(parents=True, exist_ok=True)
    for index, (variant_id, _label) in enumerate(panel):
        # Boltz names its outputs after the record, and a variant id holds ':'
        # and '+'. So the record is a plain index and the mapping is written
        # beside it rather than parsed back out of a filename.
        name = f"d{index:04d}"
        names[name] = variant_id
        sequence = layout.scfv if variant_id in ("", "parent") else apply_id(layout, variant_id)
        path = inputs / f"{name}.fa"
        path.write_text(f">{name}\n{sequence}:{antigen}\n", encoding="utf-8")
        designs.append(Design(id=name, path=path, parent_id=None, producer="ipsae_panel"))
    (out / "names.json").write_text(json.dumps(names, indent=1), encoding="utf-8")

    executable = args.boltz or "boltz"
    print(f"folding {len(designs)} designs x {args.seeds} seed(s) with {executable}; "
          f"alignments from {args.msa_cache}", flush=True)
    for seed in range(args.seed_start, args.seed_start + args.seeds):
        began = time.perf_counter()
        fold_batch(cache, designs, out / f"seed{seed}", seed=seed, executable=executable)
        spent = time.perf_counter() - began
        print(f"seed {seed}: {spent:.0f} s for {len(designs)} folds "
              f"({spent / len(designs):.1f} s each)", flush=True)
    return score(args)


def score(args: argparse.Namespace) -> int:
    import numpy as np

    from design_pipeline.scoring import ipsae

    out = Path(args.out)
    names = json.loads((out / "names.json").read_text(encoding="utf-8"))
    panel = dict(read_panel(args.panel, Path(args.out_root)))

    rows = []
    for pae_path in sorted(out.rglob("pae_*_model_0.npz")):
        name = pae_path.name[len("pae_"):-len("_model_0.npz")]
        folder = pae_path.parent
        cif = folder / f"{name}_model_0.cif"
        confidence = folder / f"confidence_{name}_model_0.json"
        if not cif.is_file() or not confidence.is_file():
            continue
        seed = re.search(r"seed(\d+)", str(pae_path.relative_to(out)))
        chains = _cif_chains(cif)
        summary = json.loads(confidence.read_text(encoding="utf-8"))
        scores = ipsae(np.load(pae_path)["pae"], chains, pae_cutoff=args.pae_cutoff)
        rows.append({
            "variant_id": names.get(name, name),
            "name": name,
            "seed": int(seed.group(1)) if seed else 0,
            "iptm": round(summary["iptm"], 4),
            "plddt": round(100.0 * summary["complex_plddt"], 2),
            "ipsae_ab": round(scores["ab"], 4),
            "ipsae_ba": round(scores["ba"], 4),
            "ipsae_max": round(scores["max"], 4),
            "ipsae_min": round(scores["min"], 4),
        })
    if not rows:
        raise SystemExit(f"no folds found under {out}")

    out.mkdir(parents=True, exist_ok=True)
    with (out / "folds.csv").open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=FOLD_FIELDS)
        writer.writeheader()
        writer.writerows(rows)

    by: dict[str, list[dict]] = {}
    for row in rows:
        by.setdefault(row["variant_id"], []).append(row)

    # The scored column: the MEAN over seeds of the larger chain direction. Not
    # the max over seeds - see the module docstring - and the per-seed spread is
    # carried beside it so a reader can see how much of a gap is noise.
    table = []
    for variant_id, folds in sorted(by.items()):
        values = [f["ipsae_max"] for f in folds]
        table.append({
            "variant_id": variant_id,
            "score_ipsae": round(statistics.mean(values), 4),
            "seeds": len(values),
            "seed_sd": round(statistics.stdev(values), 4) if len(values) > 1 else "",
            "seed_spread": round(max(values) - min(values), 4),
            "ipsae_ab_mean": round(statistics.mean(f["ipsae_ab"] for f in folds), 4),
            "ipsae_ba_mean": round(statistics.mean(f["ipsae_ba"] for f in folds), 4),
            "iptm_mean": round(statistics.mean(f["iptm"] for f in folds), 4),
            "label": panel.get(variant_id, ""),
        })
    path = out / f"{args.panel}_ipsae.csv"
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(table[0]))
        writer.writeheader()
        writer.writerows(table)

    missing = sorted(set(panel) - set(by))
    spreads = [r["seed_spread"] for r in table if r["seeds"] > 1]
    print(f"{len(rows)} folds of {len(by)} designs -> {path}")
    if missing:
        print(f"{len(missing)} design(s) of the panel have no fold yet: "
              f"{', '.join(missing[:6])}{' ...' if len(missing) > 6 else ''}")
    if spreads:
        print(f"per-design seed spread of ipsae_max: median "
              f"{statistics.median(spreads):.3f}, max {max(spreads):.3f}")
        print("Read that against the spread between designs before believing any")
        print("ranking off this column: if they are comparable, the column is noise.")
    print(f"\n{'variant':<22}{'score_ipsae':>12}{'seeds':>7}{'spread':>8}{'label':>7}")
    shown = table if len(table) <= 30 else table[:15] + [None] + table[-15:]
    for row in shown:
        if row is None:
            print("...")
            continue
        print(f"{row['variant_id'][:22]:<22}{row['score_ipsae']:>12.3f}"
              f"{row['seeds']:>7}{row['seed_spread']:>8.3f}{str(row['label']):>7}")
    print(f"\nThe notebook reads {path.name} as another scorer column, so its AUROC,")
    print("bootstrap interval and permutation p come from the same code as the other")
    print("seven rather than from a second implementation here.")
    return 0


def _cif_chains(path: Path):
    """The chain id of each CA atom, in file order: one entry per PAE row."""
    from biotite.structure.io import pdbx

    atoms = pdbx.get_structure(pdbx.CIFFile.read(str(path)), model=1)
    return atoms[atoms.atom_name == "CA"].chain_id


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("command", choices=("parent-fasta", "fold", "score"))
    parser.add_argument("--panel", default="direction", choices=("direction", "chinery"))
    parser.add_argument("--pdb", default="data/structures/1N8Z.pdb")
    parser.add_argument("--out-root", default="runs/negative_controls",
                        help="where the notebook wrote the panel labels")
    parser.add_argument("--out", default=None,
                        help="default: <out-root>/ipsae/<panel>")
    parser.add_argument("--msa-cache", default="data/msa_candidates/trastuzumab_her2")
    parser.add_argument("--boltz", default=None,
                        help="the boltz executable; default: $DESIGN_PIPELINE_BOLTZ or PATH")
    parser.add_argument("--seeds", type=int, default=3)
    parser.add_argument("--seed-start", type=int, default=0)
    parser.add_argument("--pae-cutoff", type=float, default=15.0)
    args = parser.parse_args()
    if args.out is None:
        args.out = str(Path(args.out_root) / "ipsae" / args.panel)
    if args.boltz is None:
        import os

        args.boltz = os.environ.get("DESIGN_PIPELINE_BOLTZ") or None
    if args.command == "parent-fasta":
        return parent_fasta(args)
    return fold(args) if args.command == "fold" else score(args)


if __name__ == "__main__":
    raise SystemExit(main())
