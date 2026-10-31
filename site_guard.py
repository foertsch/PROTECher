"""Same-site guard for co-folded antibody-antigen complexes.

A co-folding model can dock the antibody on the wrong surface of the antigen and still report a plausible ipTM. A fold is
accepted only if (1) its ipTM is at least IPTM_MIN and (2) its pose agrees with a reference complex of the same antibody
and antigen: DockQ at least DOCKQ_MIN and fnat at least FNAT_MIN, from the DockQ package. A fold that fails is not used
for any structure-based score.

The reference is the experimental complex if there is one. Otherwise it is the best wild-type fold (highest ipTM, on the
smallest antigen construct that holds the epitope), chosen without looking at any crystal, and it needs a second opinion
before anything is gated against it: an ipTM of 0.8 or more does not make it right. `confirm_reference` asks whether the
best fold of a second co-folding model puts the antibody on the same site; if not, nothing is gated.

    python -m site_guard [--iptm-min X] [--confirm-with OTHER.cif] REFERENCE.cif MODEL.cif [MODEL.cif ...]

Chains are A = antibody (the scFv as one chain, or heavy and light merged) and B = antigen, as Boltz writes them; ipTM comes
from the confidence_<model>.json beside each model. With --confirm-with the reference is checked against OTHER.cif first (pose
only, no ipTM) and the run stops with status 1 if they disagree. --iptm-min replaces IPTM_MIN, which is Boltz's scale. Exit
status 0 if every model passes, 1 if any fails or the reference is not confirmed, 2 on bad usage.
DockQ is the executable in the DOCKQ environment variable, else $PROTECHER_SCRATCH/protecher_tools/dockq/bin/DockQ, where
PROTECHER_SCRATCH is the node-local scratch directory (source the cluster config first). It lives in its own environment because
the package pins numpy below 2.

The thresholds were set on one complex (trastuzumab with HER2, 24 folds): unrestrained wrong poses had ipTM up to 0.59 and right
poses 0.82 or more, but a fold restrained to a wrong epitope reached ipTM 0.80 with DockQ 0.01 against the reference, which
only the pose check catches. They were rechecked on five other complexes (pertuzumab, cetuximab, pembrolizumab, atezolizumab,
sacituzumab; Boltz-2 and Chai-1, 84 folds, two complexes right at the reference and three wrong):
- ipTM does not separate the sites. 6 of 32 wrong-site Boltz folds reached 0.81 to 0.94 (PD-L1, Trop2); the pose check does the work.
- A blind reference can be wrong with a high ipTM (PD-L1 0.94, Trop2 0.81). With the best full-antigen Boltz fold as the reference
  and the best full-antigen Chai fold as the yardstick, the reference matched it for both right complexes and for none of the
  three wrong ones, five out of five. The direction matters: DockQ is anchored on its second argument, and with the roles swapped
  one right pair scored 0.20.
- The guard compares with the reference, so a fold that repeats a wrong reference's pose passes the pose check. One Boltz fold did.
- IPTM_MIN is Boltz's scale. Right-site Chai folds scored 0.60 to 0.81 (3 of 12 below 0.70), wrong ones 0.52 at most.
Five complexes are five observations, not 84 folds: recheck again on a new antibody.
"""

from __future__ import annotations

import argparse
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


def check(
    model: Path, reference: Path, iptm: float | None = None, workdir: Path | None = None, iptm_min: float | None = IPTM_MIN
) -> dict:
    """Pass or fail for one fold, with the numbers and the reasons. Any failure to compare counts as a fail.

    `iptm_min` is the ipTM floor for this model; IPTM_MIN is Boltz's scale. None skips the ipTM test, and then no
    confidence file is read.
    """
    model, reference = Path(model), Path(reference)
    if iptm is None and iptm_min is not None:
        iptm = iptm_of(model)
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
    if iptm_min is not None and iptm < iptm_min:
        reasons.append(f"ipTM {iptm:.2f} < {iptm_min}")
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


def confirm_reference(reference: Path, other: Path, workdir: Path | None = None) -> dict:
    """Does a second model's best fold back a blind reference? Pose only: DockQ and fnat of `reference` against `other`.

    `other` is the yardstick, the best-ipTM fold of another co-folding model on the same antigen, and the reference is
    scored as the model against it. That is the direction the five-complex test used, and it is not interchangeable:
    DockQ is anchored on its second argument. No ipTM test and no confidence file: the other model's ipTM is on another scale.
    """
    result = check(reference, other, iptm_min=None, workdir=workdir)
    return dict(
        confirmed=result["passed"],
        dockq_other=result.get("dockq_ref"),
        fnat_other=result.get("fnat_ref"),
        reasons=[reason.replace("vs reference", "vs the other model's fold") for reason in result["reasons"]],
    )


def main(argv: list[str]) -> int:
    parser = argparse.ArgumentParser(prog="site_guard", add_help=False)
    parser.add_argument("paths", nargs="*")
    parser.add_argument("--iptm-min", type=float, default=IPTM_MIN)
    parser.add_argument("--confirm-with", type=Path)
    try:
        args = parser.parse_args(argv[1:])
    except SystemExit:
        print(__doc__)
        return 2
    if len(args.paths) < 2:
        print(__doc__)
        return 2
    reference, models = Path(args.paths[0]), [Path(p) for p in args.paths[1:]]
    if args.confirm_with is not None:
        confirmation = confirm_reference(reference, args.confirm_with)
        numbers = f"DockQ {confirmation['dockq_other']:.2f}  fnat {confirmation['fnat_other']:.2f}  " if confirmation["dockq_other"] is not None else ""
        print(f"{'CONFIRMED' if confirmation['confirmed'] else 'NOT CONFIRMED'}  {reference.name} vs {args.confirm_with.name}  {numbers}" + "; ".join(confirmation["reasons"]))
        if not confirmation["confirmed"]:
            return 1
    failed = 0
    for model in models:
        result = check(model, reference, iptm_min=args.iptm_min)
        failed += not result["passed"]
        numbers = f"DockQ vs ref {result['dockq_ref']:.2f}  fnat {result['fnat_ref']:.2f}  " if "dockq_ref" in result else ""
        print(f"{'PASS' if result['passed'] else 'FAIL'}  {model.name}  ipTM {result['iptm']:.2f}  {numbers}" + "; ".join(result["reasons"]))
    return 1 if failed else 0


if __name__ == "__main__":
    sys.exit(main(sys.argv))
