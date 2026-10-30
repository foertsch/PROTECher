"""Same-site guard for co-folded antibody-antigen complexes.

A co-folding model can dock the antibody on the wrong surface of the antigen and still report a plausible ipTM. A fold is
accepted only if (1) its ipTM is at least IPTM_MIN and (2) its pose agrees with a reference complex of the same antibody
and antigen: DockQ at least DOCKQ_MIN and fnat at least FNAT_MIN, from the DockQ package. A fold that fails is not used
for any structure-based score.

The reference is the experimental complex if there is one, otherwise the best wild-type fold (highest ipTM, at least 0.8,
on the smallest antigen construct that holds the epitope), chosen without looking at any crystal.

    python -m site_guard REFERENCE.cif MODEL.cif [MODEL.cif ...]

Chains are A = antibody (the scFv as one chain, or heavy and light merged) and B = antigen, as Boltz writes them; ipTM comes
from the confidence_<model>.json beside each model. Exit status 0 if every model passes, 1 if any fails, 2 on bad usage.
DockQ is the executable in the DOCKQ environment variable, else $PROTECHER_SCRATCH/protecher_tools/dockq/bin/DockQ, where
PROTECHER_SCRATCH is the node-local scratch directory (source the cluster config first). It lives in its own environment because
the package pins numpy below 2.

The thresholds come from one complex (trastuzumab with HER2, 24 folds): unrestrained wrong poses had ipTM up to 0.59 and right
poses 0.82 or more, but a fold restrained to a wrong epitope reached ipTM 0.80 with DockQ 0.01 against the reference, which
only the pose check catches. Recheck them on a new antibody before relying on them.
"""

from __future__ import annotations

import json
import os
import subprocess
import sys
import tempfile
from pathlib import Path

import biotite.structure as struc
import biotite.structure.io as strucio

DOCKQ_DEFAULT = os.path.join(os.environ.get("PROTECHER_SCRATCH", ""), "protecher_tools/dockq/bin/DockQ")
DOCKQ = os.environ.get("DOCKQ", DOCKQ_DEFAULT)
IPTM_MIN = 0.70
DOCKQ_MIN = 0.49
"""Lower edge of DockQ's own "medium quality" class."""
FNAT_MIN = 0.50
"""Half of the reference contacts recovered."""
MISMATCHES = 30
"""A design differs from the reference by a few mutations; DockQ must still map the chains."""


def to_pdb(path: Path, out: Path) -> Path:
    """The protein heavy atoms of a structure file as a PDB, which is what DockQ reads reliably from Boltz output."""
    atoms = strucio.load_structure(str(path))
    atoms = atoms[struc.filter_amino_acids(atoms)]
    strucio.save_structure(str(out), atoms[atoms.element != "H"])
    return out


def iptm_of(model: Path) -> float:
    name = model.name.replace(".cif", "").replace(".pdb", "")
    return float(json.loads((model.parent / f"confidence_{name}.json").read_text())["iptm"])


def check(model: Path, reference: Path, iptm: float | None = None, workdir: Path | None = None) -> dict:
    """Pass or fail for one fold, with the numbers and the reasons. Any failure to compare counts as a fail."""
    model, reference = Path(model), Path(reference)
    iptm = iptm_of(model) if iptm is None else iptm
    with tempfile.TemporaryDirectory(dir=workdir) as tmp:
        tmp = Path(tmp)
        m = to_pdb(model, tmp / "model.pdb") if model.suffix == ".cif" else model
        r = to_pdb(reference, tmp / "ref.pdb") if reference.suffix == ".cif" else reference
        out = tmp / "dockq.json"
        try:
            run = subprocess.run(
                [DOCKQ, str(m), str(r), "--mapping", "AB:AB", "--allowed_mismatches", str(MISMATCHES), "--json", str(out), "--n_cpu", "1"],
                capture_output=True,
                text=True,
            )
        except FileNotFoundError:
            return dict(iptm=iptm, passed=False, reasons=[f"DockQ not found at {DOCKQ}, set the DOCKQ environment variable or PROTECHER_SCRATCH"])
        if run.returncode != 0 or not out.exists():
            return dict(iptm=iptm, passed=False, reasons=["DockQ failed: " + (run.stderr or run.stdout)[-160:].replace("\n", " ")])
        res = list(json.loads(out.read_text())["best_result"].values())[0]
    reasons = []
    if iptm < IPTM_MIN:
        reasons.append(f"ipTM {iptm:.2f} < {IPTM_MIN}")
    if res["DockQ"] < DOCKQ_MIN:
        reasons.append(f"DockQ vs reference {res['DockQ']:.2f} < {DOCKQ_MIN}")
    if res["fnat"] < FNAT_MIN:
        reasons.append(f"fnat vs reference {res['fnat']:.2f} < {FNAT_MIN}")
    return dict(
        iptm=iptm,
        dockq_ref=res["DockQ"],
        fnat_ref=res["fnat"],
        irmsd_ref=res["iRMSD"],
        lrmsd_ref=res["LRMSD"],
        passed=not reasons,
        reasons=reasons,
    )


def main(argv: list[str]) -> int:
    if len(argv) < 3:
        print(__doc__)
        return 2
    reference, models = Path(argv[1]), [Path(p) for p in argv[2:]]
    failed = 0
    for model in models:
        result = check(model, reference)
        failed += not result["passed"]
        numbers = f"DockQ vs ref {result['dockq_ref']:.2f}  fnat {result['fnat_ref']:.2f}  " if "dockq_ref" in result else ""
        print(f"{'PASS' if result['passed'] else 'FAIL'}  {model.name}  ipTM {result['iptm']:.2f}  {numbers}" + "; ".join(result["reasons"]))
    return 1 if failed else 0


if __name__ == "__main__":
    sys.exit(main(sys.argv))
