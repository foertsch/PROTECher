"""charge_contact: the salt-bridge ramp and screened-Coulomb term on structures whose geometry is set by hand.

The pure-numpy tests (the ramp, the Coulomb sign and screening) need no biotite. The structure tests build tiny two-chain
complexes and need biotite, so they skip where it is missing (it pins against AntiFold, so it is not in every environment).
"""

from __future__ import annotations

import numpy as np
import pytest

from design_pipeline.scoring import (
    CONTACT_FULL,
    CONTACT_ZERO,
    _contact_score,
    _screened_coulomb,
    charge_contact_report,
)


# --------------------------------------------------------------------------- pure functions, no biotite
def test_contact_ramp_is_one_up_to_the_full_distance():
    assert _contact_score(2.0) == 1.0
    assert _contact_score(CONTACT_FULL) == 1.0


def test_contact_ramp_falls_linearly_to_zero():
    assert _contact_score(3.75) == pytest.approx(0.5)      # halfway between 3.5 and 4.0
    assert _contact_score(3.8) == pytest.approx(0.4)


def test_contact_ramp_is_zero_at_and_beyond_the_cutoff():
    assert _contact_score(CONTACT_ZERO) == 0.0
    assert _contact_score(4.5) == 0.0


def test_coulomb_opposite_charges_attract_is_negative():
    e = _screened_coulomb([[0, 0, 0]], [1.0], [[0, 0, 3.0]], [-1.0])
    assert e < 0


def test_coulomb_like_charges_repel_is_positive():
    e = _screened_coulomb([[0, 0, 0]], [1.0], [[0, 0, 3.0]], [1.0])
    assert e > 0


def test_coulomb_screening_weakens_with_distance():
    near = _screened_coulomb([[0, 0, 0]], [1.0], [[0, 0, 3.0]], [-1.0])
    far = _screened_coulomb([[0, 0, 0]], [1.0], [[0, 0, 9.0]], [-1.0])
    assert abs(far) < abs(near)


def test_coulomb_is_zero_with_no_charges():
    assert _screened_coulomb(np.empty((0, 3)), [], [[0, 0, 0]], [-1.0]) == 0.0


# --------------------------------------------------------------------------- structure tests, need biotite
pytest.importorskip("biotite")
from biotite.structure import Atom, array  # noqa: E402
import biotite.structure.io.pdbx as pdbx   # noqa: E402

# A side chain is only ever read through the atoms the scorer names, so a residue carries just those plus a CA.
_SIDE = {
    "ASP": {"CA": (0, 0, 0), "CB": (0.5, 0.5, 0), "CG": (1.0, 0, 0), "OD1": (1.5, 0.7, 0), "OD2": (1.5, -0.7, 0)},
    "LYS": {"CA": (0, 0, 0), "CB": (0.5, 0.5, 0), "CG": (1.0, 0, 0), "CD": (1.5, 0, 0), "CE": (2.0, 0, 0), "NZ": (2.5, 0, 0)},
    "ALA": {"CA": (0, 0, 0), "CB": (0.5, 0.5, 0)},
}


def _residue(chain, res_id, res_name, origin, nz_at=None):
    atoms = []
    for name, (x, y, z) in _SIDE[res_name].items():
        coord = np.array(origin, float) + np.array([x, y, z], float)
        if nz_at is not None and name == "NZ":
            coord = np.array(nz_at, float)
        atoms.append(Atom(coord, chain_id=chain, res_id=res_id, res_name=res_name, atom_name=name, element=name[0]))
    return atoms


def _write(atoms, path):
    f = pdbx.CIFFile()
    pdbx.set_structure(f, array(atoms))
    f.write(str(path))


def _parent(tmp_path):
    """scFv chain A: ASP1 near the antigen's ASP100, a far ALA2 spacer. Antigen chain B: ASP100 at the origin cluster."""
    atoms = []
    atoms += _residue("B", 100, "ASP", (0, 0, 0))                 # OD1 at ~(1.5, 0.7, 0)
    atoms += _residue("A", 1, "ASP", (0, 0, 4.0))                 # within 5 A of the antigen -> in the epitope
    atoms += _residue("A", 2, "ALA", (0, 0, 40.0))               # far away, never in the epitope
    path = tmp_path / "parent.cif"
    _write(atoms, path)
    return path


def test_parent_against_itself_scores_zero(tmp_path):
    parent = _parent(tmp_path)
    r = charge_contact_report(parent, parent)
    assert r["contact"] == 0.0
    assert r["coulomb"] == pytest.approx(0.0, abs=1e-9)
    assert r["n_added_cation"] == 0.0
    assert r["eligible"] == 1.0                                   # same structure, epitope overlap 1


def test_aspartate_to_lysine_near_a_carboxylate_rewards_contact_and_coulomb(tmp_path):
    parent = _parent(tmp_path)
    # design: ASP1 -> LYS1, its NZ placed 3.0 A from the antigen OD1 at (1.5, 0.7, 0)
    atoms = []
    atoms += _residue("B", 100, "ASP", (0, 0, 0))
    atoms += _residue("A", 1, "LYS", (0, 0, 4.0), nz_at=(1.5, 0.7, 3.0))
    atoms += _residue("A", 2, "ALA", (0, 0, 40.0))
    design = tmp_path / "design.cif"
    _write(atoms, design)

    r = charge_contact_report(design, parent)
    assert r["n_added_cation"] == 1.0
    assert r["contact"] == pytest.approx(1.0)                     # NZ...OD1 is 3.0 A, below the full-contact distance
    assert r["coulomb"] > 0                                       # -1/-1 repulsion replaced by +1/-1 attraction
    assert r["eligible"] == 1.0
    assert r["paratope_net_charge"] == pytest.approx(1.0)         # the one epitope scFv residue is now a Lys


def test_a_cation_beyond_the_cutoff_scores_no_contact_but_still_counts(tmp_path):
    parent = _parent(tmp_path)
    atoms = []
    atoms += _residue("B", 100, "ASP", (0, 0, 0))
    atoms += _residue("A", 1, "LYS", (0, 0, 4.0), nz_at=(1.5, 0.7, 5.0))   # NZ...OD1 is 5.0 A
    atoms += _residue("A", 2, "ALA", (0, 0, 40.0))
    design = tmp_path / "design_far.cif"
    _write(atoms, design)

    r = charge_contact_report(design, parent)
    assert r["n_added_cation"] == 1.0
    assert r["contact"] == 0.0


def test_an_off_site_pose_is_marked_ineligible(tmp_path):
    parent = _parent(tmp_path)
    # the whole scFv is shoved 50 A away: it no longer touches the antigen, so the epitope is empty
    atoms = []
    atoms += _residue("B", 100, "ASP", (0, 0, 0))
    atoms += _residue("A", 1, "LYS", (0, 0, 54.0), nz_at=(1.5, 0.7, 53.0))
    atoms += _residue("A", 2, "ALA", (0, 0, 90.0))
    design = tmp_path / "design_offsite.cif"
    _write(atoms, design)

    r = charge_contact_report(design, parent)
    assert r["epitope_overlap"] == 0.0
    assert r["eligible"] == 0.0
    assert r["contact"] == 0.0                                    # NZ is nowhere near a carboxylate either


def test_none_when_the_reference_is_not_a_two_chain_complex(tmp_path):
    parent = _parent(tmp_path)
    monomer = array(_residue("A", 1, "ASP", (0, 0, 0)))
    f = tmp_path / "monomer.cif"
    _write(monomer, f)
    assert charge_contact_report(parent, f) is None
