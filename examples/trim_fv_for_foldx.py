"""Build the structure FoldX scores against: the Fv plus the antigen domain, trimmed.

`DESIGN_PIPELINE_FOLDX_PDB` has no default, and what you point it at decides both the
cost and the comparability of every ddG. This script writes the same structure the
alanine scan behind `design_pipeline.hotspots` used, so the scan's hotspot values and a
campaign's ddGs are on the same footing:

    chain A  1-107      VL   (not CL, 108-214)
    chain B  1-120      VH   (not CH1, 121-220)
    chain C  480-607    the HER2 domain near the paratope

227 antibody residues plus 118 antigen, 345 in total against the Fab's 552. That matters
because FoldX repacks the whole structure on every BuildModel, and a campaign pays that
seven times per candidate (one BuildModel at --numberOfRuns=3, then two AnalyseComplex
per model). Measured on the Nebius box: ~1 candidate/min on the 552-residue Fab against
~9-10 on this one, for candidates of 5-10 substitutions.

The boundaries are not guesses - they are what `trastuzumab_foldx_alanine_scan.csv`
records in its `domain` column, and `docs/track_a_filtering.md` describes as "an alanine
scan of the 1N8Z scFv with HER2 residues 480 to 607".

    python3 examples/trim_fv_for_foldx.py --pdb data/structures/1N8Z.pdb \
        --out data/structures/1N8Z_fv_trimmed.pdb

Then repair it once, which is FoldX's own prerequisite, and export the result:

    cd tools/foldx && ./foldx_<expiry> --command=RepairPDB \
        --pdb=1N8Z_fv_trimmed.pdb --pdb-dir=<abs data/structures> --output-dir=<same>
    export DESIGN_PIPELINE_FOLDX_PDB=<abs>/1N8Z_fv_trimmed_Repair.pdb

Needs biotite, which `core` already has.
"""

from __future__ import annotations

import argparse
from pathlib import Path

import biotite.structure as struc
import biotite.structure.io.pdb as pdb

# chain -> inclusive residue range to keep. None means "no upper bound".
FV_PLUS_ANTIGEN: dict[str, tuple[int, int]] = {
    "A": (1, 107),      # VL
    "B": (1, 120),      # VH
    "C": (480, 607),    # HER2, the part near the paratope
}


def trim(path: Path, keep: dict[str, tuple[int, int]]) -> struc.AtomArray:
    """The amino-acid atoms of `keep`, and nothing else - waters and ligands included."""
    array = pdb.PDBFile.read(str(path)).get_structure(model=1)
    array = array[struc.filter_amino_acids(array)]
    mask = None
    for chain, (low, high) in keep.items():
        chain_mask = (array.chain_id == chain) & (array.res_id >= low) & (array.res_id <= high)
        mask = chain_mask if mask is None else (mask | chain_mask)
    if mask is None or not mask.any():
        raise SystemExit(f"nothing to keep in {path}: check the chain ids against {sorted(set(array.chain_id))}")
    return array[mask]


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--pdb", type=Path, required=True, help="the full structure, e.g. data/structures/1N8Z.pdb")
    parser.add_argument("--out", type=Path, required=True, help="where to write the trimmed structure")
    args = parser.parse_args(argv)

    trimmed = trim(args.pdb, FV_PLUS_ANTIGEN)
    out = pdb.PDBFile()
    out.set_structure(trimmed)
    args.out.parent.mkdir(parents=True, exist_ok=True)
    out.write(str(args.out))

    for chain in sorted(set(trimmed.chain_id)):
        sub = trimmed[trimmed.chain_id == chain]
        print(f"chain {chain}: {len(set(sub.res_id))} residues ({min(sub.res_id)}-{max(sub.res_id)})")
    print(f"{len(set(zip(trimmed.chain_id, trimmed.res_id)))} residues, {trimmed.array_length()} atoms -> {args.out}")
    print("now run foldx --command=RepairPDB on it; see this file's docstring.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
