"""
ipsae_score.py -- per-fold ipSAE scoring for PROTECher.

Wraps Dunbrack's reference ipsae.py (pinned commit) to score ONE Boltz-2 fold
(coordinates + full PAE) for the scFv<->HER2 interface.

This is a STRUCTURE scorer: it consumes Boltz output, not a sequence. It is a
sibling of site_guard, NOT of the sequence scorers (antifold_seq, ablang2, esm).
Do not register it in the sequence-scorer registry.

Reproducibility note: the ipSAE number depends on WHICH commit of ipsae.py you
run (v3 -> v4 changed the d0 calculation). Vendor ONE ipsae.py, pin it, and
point $IPSAE_SRC at it. The commit hash is recorded in every row.

Usage (library):
    from ipsae_score import score_fold
    row = score_fold(cif_path, pae_npz_path, design_id="d017", seed=3)

Usage (CLI, for one fold):
    python ipsae_score.py model_0.cif pae_model_0.npz --design d017 --seed 3
"""
from __future__ import annotations

import argparse
import glob
import json
import os
import shutil
import subprocess
import tempfile
from dataclasses import dataclass, asdict
from pathlib import Path
from typing import Optional

import pandas as pd

# Point this at your vendored, pinned copy. Overridable via env.
IPSAE_SRC = os.environ.get("IPSAE_SRC", str(Path(__file__).resolve().parent / "ipsae.py"))


class NoInterface(Exception):
    """No scFv<->HER2 contact rows in the ipsae output (valid reject)."""


@dataclass
class IpsaeRow:
    design_id: str
    seed: int
    construct: str            # "full_ecd" | "domain_iv" | ...
    pae_cut: float
    dist_cut: float
    scfv_chain: str
    her2_chain: str
    ipsae_min: float          # min over directions -- the organisers' "ipSAE_min"
    ipsae_max: float          # the max row of ipsae.py
    ipsae_scfv_to_her2: float # scFv residues scored against HER2 (d0res)
    ipsae_her2_to_scfv: float
    iptm_af: float
    n_contacts: int           # dist1(AB) + dist1(BA): interface-size proxy
    pose_pass: Optional[bool] # filled by site_guard / rank layer; None here

    def as_dict(self) -> dict:
        return asdict(self)


def _pad(cut: float) -> str:
    # ipsae.py zero-pads cutoffs below 10 in the output filename (pae=8 -> "08").
    s = str(int(cut))
    return "0" + s if int(cut) < 10 else s


def _run_ipsae(pae_path, struct_path, pae_cut, dist_cut, ipsae_src, workdir) -> pd.DataFrame:
    # Copy inputs into an isolated workdir so ipsae.py writes its outputs there,
    # never next to the originals.
    pae_local = Path(workdir) / Path(pae_path).name
    struct_local = Path(workdir) / Path(struct_path).name
    shutil.copy(pae_path, pae_local)
    shutil.copy(struct_path, struct_local)
    subprocess.run(
        ["python", str(ipsae_src), str(pae_local), str(struct_local),
         str(int(pae_cut)), str(int(dist_cut))],
        cwd=workdir, check=True, capture_output=True, text=True,
    )
    # ipsae.py names the main table from the STRUCTURE stem with zero-padded
    # cutoffs: "<struct-without-ext>_<PP>_<DD>.txt" (and "..._byres.txt").
    name = struct_local.name
    base = name[:-4] if name.endswith((".cif", ".pdb")) else struct_local.stem
    expected = Path(workdir) / f"{base}_{_pad(pae_cut)}_{_pad(dist_cut)}.txt"
    if not expected.exists():
        cands = [f for f in glob.glob(str(Path(workdir) / "*.txt"))
                 if not f.endswith("_byres.txt")]
        if not cands:
            raise FileNotFoundError("ipsae.py produced no chain-pair table")
        expected = Path(cands[0])
    return pd.read_csv(expected, sep=r"\s+")


def _extract_pair(df: pd.DataFrame,
                  scfv_chain: Optional[str],
                  her2_chain: Optional[str]) -> dict:
    df.columns = [c.strip() for c in df.columns]
    chains = sorted(set(df["Chn1"]) | set(df["Chn2"]))

    # scFv is one chain, HER2 (typically) one chain -> exactly 2 chains,
    # so the pair is unambiguous. If HER2 folds as >2 chains, pass ids explicitly.
    if scfv_chain is None or her2_chain is None:
        if len(chains) == 2:
            scfv_chain = scfv_chain or chains[0]
            her2_chain = her2_chain or chains[1]
        else:
            raise ValueError(
                f"{len(chains)} chains {chains}; pass scfv_chain/her2_chain "
                f"explicitly (confirm the scFv chain id on one real fold first)."
            )

    asym = df[df["Type"] == "asym"]
    ab = asym[(asym["Chn1"] == scfv_chain) & (asym["Chn2"] == her2_chain)]
    ba = asym[(asym["Chn1"] == her2_chain) & (asym["Chn2"] == scfv_chain)]
    if ab.empty or ba.empty:
        raise NoInterface

    s2h = float(ab["ipSAE"].iloc[0])
    h2s = float(ba["ipSAE"].iloc[0])
    iptm = float(ab["ipTM_af"].iloc[0]) if "ipTM_af" in ab else 0.0
    ncontacts = (int(ab["dist1"].iloc[0]) if "dist1" in ab else 0) + \
                (int(ba["dist1"].iloc[0]) if "dist1" in ba else 0)
    return dict(
        ipsae_min=min(s2h, h2s), ipsae_max=max(s2h, h2s),
        ipsae_scfv_to_her2=s2h, ipsae_her2_to_scfv=h2s,
        iptm_af=iptm, n_contacts=ncontacts,
        scfv_chain=scfv_chain, her2_chain=her2_chain,
    )


def score_fold(
    struct_path: str,
    pae_path: str,
    design_id: str,
    seed: int,
    construct: str = "full_ecd",
    pae_cut: float = 10.0,
    dist_cut: float = 10.0,
    scfv_chain: Optional[str] = None,
    her2_chain: Optional[str] = None,
    ipsae_src: str = IPSAE_SRC,
) -> IpsaeRow:
    """Score one Boltz fold. pae_path is the .npz, struct_path the .cif.

    Deterministic given a fold: all seed noise comes from Boltz upstream, so
    results are safe to cache on (fold hash, pae_cut, dist_cut, ipsae_commit).
    """
    
    try:
        with tempfile.TemporaryDirectory() as wd:
            df = _run_ipsae(pae_path, struct_path, pae_cut, dist_cut, ipsae_src, wd)
        payload = _extract_pair(df, scfv_chain, her2_chain)
    except (NoInterface, subprocess.CalledProcessError, FileNotFoundError):
        # Wrong-site / no-contact fold: record a zero, don't crash the batch.
        payload = dict(
            ipsae_min=0.0, ipsae_max=0.0, ipsae_scfv_to_her2=0.0,
            ipsae_her2_to_scfv=0.0, iptm_af=0.0, n_contacts=0,
            scfv_chain=scfv_chain or "?", her2_chain=her2_chain or "?",
        )
    return IpsaeRow(
        design_id=design_id, seed=seed, construct=construct,
        pae_cut=pae_cut, dist_cut=dist_cut, pose_pass=None,
        **payload,
    )


if __name__ == "__main__":
    p = argparse.ArgumentParser(description="ipSAE for one Boltz fold")
    p.add_argument("struct"); p.add_argument("pae")
    p.add_argument("--design", required=True); p.add_argument("--seed", type=int, required=True)
    p.add_argument("--construct", default="full_ecd")
    p.add_argument("--pae-cut", type=float, default=10.0)
    p.add_argument("--dist-cut", type=float, default=10.0)
    p.add_argument("--scfv-chain", default=None); p.add_argument("--her2-chain", default=None)
    a = p.parse_args()
    row = score_fold(a.struct, a.pae, a.design, a.seed, a.construct,
                     a.pae_cut, a.dist_cut, a.scfv_chain, a.her2_chain)
    print(json.dumps(row.as_dict(), indent=2))
