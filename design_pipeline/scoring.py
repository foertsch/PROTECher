"""Scoring: the only place in this package where research judgment lives.

Pure numpy + biotite. No model, no GPU, no network.

Adding a metric means writing one function of (prediction, reference) and
putting it in METRICS. A metric returns None when it does not apply to a
design — a monomer has no interface PAE — and that simply means no score.
"""

from __future__ import annotations

import json
from dataclasses import replace
from pathlib import Path
from typing import Callable, Mapping, Sequence

import numpy as np

from .designs import Design, Scores

Metric = Callable[[Design, "Design | None"], float | None]


def _ca_atoms(path: Path):
    """CA atoms of the first model, with B-factors."""
    from biotite.structure.io.pdb import PDBFile

    structure = PDBFile.read(str(path)).get_structure(
        model=1, extra_fields=["b_factor"]
    )
    return structure[structure.atom_name == "CA"]


def sc_rmsd(prediction: Design, reference: Design | None) -> float | None:
    """Self-consistency RMSD: predicted structure against the backbone it came
    from, CA atoms, after superposition.

    Returns None when the two have different CA counts. Mapping unequal lengths
    needs RFdiffusion's contig map to say which residues correspond, and
    guessing (truncate to the shorter) would quietly produce a number that
    means nothing. Skipped until we need it.
    """
    if reference is None:
        return None
    import biotite.structure as struc

    pred = _ca_atoms(prediction.path)
    ref = _ca_atoms(reference.path)
    if pred.array_length() != ref.array_length():
        return None
    # biotite >= 2.0 returns the transformation rather than the fitted copy.
    fitted = struc.superimpose(ref, pred).apply(pred)
    return float(struc.rmsd(ref, fitted))


def mean_plddt(prediction: Design, reference: Design | None) -> float | None:
    """Mean pLDDT. Prefers the model's own number, falls back to the B-factor
    column, which is where ColabFold writes per-residue pLDDT."""
    if "plddt" in prediction.scores:
        return float(prediction.scores["plddt"])
    atoms = _ca_atoms(prediction.path)
    if not hasattr(atoms, "b_factor"):
        return None
    return float(np.mean(atoms.b_factor))


def pae_interaction(prediction: Design, reference: Design | None) -> float | None:
    """Mean PAE across the interface, averaged over both off-diagonal blocks.

    Needs the PAE matrix, which lives in the predictor's scores JSON rather
    than in the PDB; the adapter records its path in meta.
    """
    pae = _load_pae(prediction)
    if pae is None:
        return None
    chains = _ca_atoms(prediction.path).chain_id
    if len(set(chains)) < 2 or len(chains) != pae.shape[0]:
        return None
    first = chains[0]
    mask = chains == first
    # Two chains only: everything that is not the first chain is the partner.
    # A real multi-chain complex would need a chain pair argument; skipped.
    inter = np.concatenate([pae[mask][:, ~mask].ravel(), pae[~mask][:, mask].ravel()])
    return float(inter.mean())


def _load_pae(prediction: Design) -> np.ndarray | None:
    candidate = prediction.meta.get("scores_json")
    path = Path(candidate) if candidate else None
    if path is None or not path.is_file():
        siblings = sorted(prediction.path.parent.glob("*scores*.json"))
        path = siblings[0] if siblings else None
    if path is None or not path.is_file():
        return None
    try:
        payload = json.loads(path.read_text())
    except (OSError, json.JSONDecodeError):
        return None
    pae = payload.get("pae") or payload.get("predicted_aligned_error")
    return None if pae is None else np.asarray(pae, dtype=float)


# ------------------------------------------------------------------ chemistry score (charge_contact)
# A soft, rank-only score for the acidic HER2 epitope (E558 D560 D570 D596 E598 in the Boltz poses, docs/chem_score_spec.md).
# It rewards substitutions that add a positive side chain near an epitope carboxylate (short-range salt bridge) and a more
# favourable long-range electrostatic configuration (steering). It never vetoes. All geometry is from biotite; the Coulomb
# term is our own screened sum (a few lines of numpy), because APBS is heavier than the challenge day allows.

CONTACT_FULL = 3.5          # a N...O distance at or below this scores a full contact
CONTACT_ZERO = 4.0          # the whiteboard cutoff: at or beyond this a contact scores nothing
HIS_WEIGHT = 0.5            # histidine is a weak, pH-dependent cation (general knowledge, not checked here)
DEBYE_LAMBDA = 8.0         # screening length of the Coulomb term, in A
EPITOPE_CUTOFF = 5.0       # an antigen residue is in the epitope if a heavy atom is within this of the scFv
POSE_OVERLAP_MIN = 0.5     # a seed counts only if its epitope overlaps the parent's by at least this (Jaccard)

_CARBOXYL_O = {"ASP": ("OD1", "OD2"), "GLU": ("OE1", "OE2")}
_CATION_N = {"LYS": ("NZ",), "ARG": ("NE", "NH1", "NH2"), "HIS": ("ND1", "NE2")}
_CHARGE_SITE = {"ASP": ("CG", -1.0), "GLU": ("CD", -1.0), "LYS": ("NZ", 1.0), "ARG": ("CZ", 1.0)}
_FORMAL_CHARGE = {"ASP": -1, "GLU": -1, "LYS": 1, "ARG": 1}


def _contact_score(distance: float) -> float:
    """1 below CONTACT_FULL, falling linearly to 0 at CONTACT_ZERO, 0 beyond. Pure; no biotite."""
    if distance <= CONTACT_FULL:
        return 1.0
    if distance >= CONTACT_ZERO:
        return 0.0
    return (CONTACT_ZERO - distance) / (CONTACT_ZERO - CONTACT_FULL)


def _screened_coulomb(coords_a, charges_a, coords_b, charges_b, lam: float = DEBYE_LAMBDA) -> float:
    """Sum of q_i q_j exp(-r/lam) / r over all pairs (a, b). Pure numpy; no biotite, no units.

    Attraction between opposite charges is negative, so a more favourable configuration is more negative.
    """
    a = np.asarray(coords_a, dtype=float)
    b = np.asarray(coords_b, dtype=float)
    if a.size == 0 or b.size == 0:
        return 0.0
    r = np.linalg.norm(a[:, None, :] - b[None, :, :], axis=-1)
    r = np.where(r < 1e-6, 1e-6, r)
    qq = np.outer(np.asarray(charges_a, float), np.asarray(charges_b, float))
    return float((qq * np.exp(-r / lam) / r).sum())


def _load_complex(path: Path):
    """All non-hydrogen, non-hetero atoms of model 1, as an AtomArray, from a .cif or .pdb. biotite is imported here so
    the module still loads without it."""
    path = Path(path)
    if path.suffix.lower() in (".cif", ".bcif"):
        import biotite.structure.io.pdbx as pdbx

        atoms = pdbx.get_structure(pdbx.CIFFile.read(str(path)), model=1)
    else:
        from biotite.structure.io.pdb import PDBFile

        atoms = PDBFile.read(str(path)).get_structure(model=1)
    atoms = atoms[~atoms.hetero]
    return atoms[atoms.element != "H"]


def _residue_names(chain) -> dict[int, str]:
    ca = chain[chain.atom_name == "CA"]
    return dict(zip(ca.res_id.tolist(), ca.res_name.tolist()))


def _min_distance(coords_a, coords_b) -> float:
    a, b = np.asarray(coords_a, float), np.asarray(coords_b, float)
    if a.size == 0 or b.size == 0:
        return float("inf")
    return float(np.linalg.norm(a[:, None, :] - b[None, :, :], axis=-1).min())


def _epitope(scfv, antigen, cutoff: float = EPITOPE_CUTOFF) -> set[int]:
    """Antigen res_ids with a heavy atom within `cutoff` of any scFv atom."""
    if len(scfv) == 0 or len(antigen) == 0:
        return set()
    d = np.linalg.norm(scfv.coord[:, None, :] - antigen.coord[None, :, :], axis=-1)
    return set(antigen.res_id[(d <= cutoff).any(axis=0)].tolist())


def _carboxyl_oxygens(antigen):
    names = np.array([a for pair in _CARBOXYL_O.values() for a in pair])
    return antigen[np.isin(antigen.res_name, list(_CARBOXYL_O)) & np.isin(antigen.atom_name, names)].coord


def _charge_sites(chain, res_ids):
    """(coords, charges) for the formal charge centres of the ASP/GLU/LYS/ARG residues among `res_ids`."""
    coords, charges = [], []
    for rid in res_ids:
        res = chain[chain.res_id == rid]
        if len(res) == 0:
            continue
        name = res.res_name[0]
        site = _CHARGE_SITE.get(name)
        if site is None:
            continue
        atom = res[res.atom_name == site[0]]
        if len(atom):
            coords.append(atom.coord[0])
            charges.append(site[1])
    return (np.array(coords) if coords else np.empty((0, 3))), np.array(charges)


def _added_cation_contact(design_chain, mutated, antigen) -> tuple[float, int]:
    """Contact term and count, over positions that became LYS/ARG/HIS in the design."""
    oxy = _carboxyl_oxygens(antigen)
    total, n = 0.0, 0
    for rid in mutated:
        res = design_chain[design_chain.res_id == rid]
        name = res.res_name[0]
        atoms = _CATION_N.get(name)
        if atoms is None:
            continue
        n += 1
        n_coords = res[np.isin(res.atom_name, atoms)].coord
        total += (HIS_WEIGHT if name == "HIS" else 1.0) * _contact_score(_min_distance(n_coords, oxy))
    return total, n


def _split_two_chains(prediction, parent):
    """(scfv_chain_id, antigen_chain_id). The scFv is the chain whose sequence differs between the two structures,
    or, failing that (parent against parent), the longer chain."""
    ids = [c for c in dict.fromkeys(prediction.chain_id.tolist()) if c in set(parent.chain_id.tolist())]
    if len(ids) != 2:
        return None
    differing = [c for c in ids if _residue_names(prediction[prediction.chain_id == c]) != _residue_names(parent[parent.chain_id == c])]
    if len(differing) == 1:
        scfv = differing[0]
    else:
        scfv = max(ids, key=lambda c: (prediction.chain_id == c).sum())
    antigen = next(c for c in ids if c != scfv)
    return scfv, antigen


def charge_contact_report(prediction_path: Path, parent_path: Path, pose_overlap_min: float = POSE_OVERLAP_MIN) -> dict | None:
    """The chemistry terms for one predicted complex against the parent's prediction of the same seed.

    Returns None unless both structures have the same two chains. Otherwise a dict with `contact` (sum of the salt-bridge
    ramp over added cationic substitutions), `coulomb` (parent minus design screened electrostatic energy of the mutated
    charge sites against the antigen, so larger is more favourable), `epitope_overlap` with the parent, `eligible`
    (1.0 if the overlap is at least `pose_overlap_min`), `n_added_cation`, `n_mutations`, and `paratope_net_charge`.
    The terms are reported for every pose; `eligible` marks the ones a caller should average. Nothing is vetoed here.
    """
    pred, par = _load_complex(prediction_path), _load_complex(parent_path)
    chains = _split_two_chains(pred, par)
    if chains is None:
        return None
    scfv_id, antigen_id = chains
    d_scfv, d_antigen = pred[pred.chain_id == scfv_id], pred[pred.chain_id == antigen_id]
    p_scfv, p_antigen = par[par.chain_id == scfv_id], par[par.chain_id == antigen_id]

    d_names, p_names = _residue_names(d_scfv), _residue_names(p_scfv)
    mutated = [rid for rid in d_names if rid in p_names and d_names[rid] != p_names[rid]]

    contact, n_cation = _added_cation_contact(d_scfv, mutated, d_antigen)

    d_sites = _charge_sites(d_scfv, mutated)
    p_sites = _charge_sites(p_scfv, mutated)
    d_ant_sites = _antigen_charge_sites(d_antigen)
    p_ant_sites = _antigen_charge_sites(p_antigen)
    e_design = _screened_coulomb(d_sites[0], d_sites[1], d_ant_sites[0], d_ant_sites[1])
    e_parent = _screened_coulomb(p_sites[0], p_sites[1], p_ant_sites[0], p_ant_sites[1])
    coulomb = e_parent - e_design

    epi_d = _epitope(d_scfv, d_antigen)
    epi_p = _epitope(p_scfv, p_antigen)
    overlap = len(epi_d & epi_p) / len(epi_d | epi_p) if (epi_d or epi_p) else 0.0

    para = [d_names[rid] for rid in _epitope(d_antigen, d_scfv)]
    net = sum(_FORMAL_CHARGE.get(n, 0) for n in para)

    return {
        "contact": contact,
        "coulomb": coulomb,
        "epitope_overlap": overlap,
        "eligible": 1.0 if overlap >= pose_overlap_min else 0.0,
        "n_added_cation": float(n_cation),
        "n_mutations": float(len(mutated)),
        "paratope_net_charge": float(net),
    }


def _antigen_charge_sites(antigen):
    """(coords, charges) for the ASP/GLU/LYS/ARG formal charge centres of a whole antigen chain."""
    res_ids = sorted(set(antigen.res_id[np.isin(antigen.res_name, list(_CHARGE_SITE))].tolist()))
    return _charge_sites(antigen, res_ids)


def _charge_contact_metric(key: str) -> Metric:
    def metric(prediction: Design, reference: Design | None) -> float | None:
        if reference is None:
            return None
        try:
            report = charge_contact_report(prediction.path, reference.path)
        except Exception:
            return None            # a metric yields no score where it does not apply; a bad file is not fatal
        return None if report is None else report[key]

    metric.__name__ = f"charge_contact_{key}"
    metric.__doc__ = f"charge_contact {key} against the reference complex; None without a two-chain reference (docs/chem_score_spec.md)."
    return metric


METRICS: dict[str, Metric] = {
    "sc_rmsd": sc_rmsd,
    "mean_plddt": mean_plddt,
    "pae_interaction": pae_interaction,
    **{f"charge_contact_{key}": _charge_contact_metric(key) for key in ("contact", "coulomb", "eligible", "epitope_overlap", "paratope_net_charge")},
}


def score_designs(
    predictions: Sequence[Design],
    reference: Sequence[Design],
    metrics: Mapping[str, Metric] | None = None,
) -> list[Design]:
    """Return the predictions with `scores` populated.

    Each prediction is paired with its reference backbone through the ancestry
    the steps record: `parent_id` points at the sequence it was folded from, so
    reaching the backbone two steps up uses the ancestor ids in meta.
    """
    index = {d.id: d for d in reference}
    chosen = METRICS if metrics is None else metrics
    scored: list[Design] = []
    for prediction in predictions:
        ref = _reference_for(prediction, index)
        values = {
            name: value
            for name, metric in chosen.items()
            if (value := metric(prediction, ref)) is not None
        }
        scored.append(
            replace(
                prediction,
                scores=Scores({**prediction.scores, **values}),
                meta={**prediction.meta, "reference_id": ref.id if ref else None},
            )
        )
    return scored


def _reference_for(prediction: Design, index: Mapping[str, Design]) -> Design | None:
    ancestors = prediction.meta.get("ancestors") or ()
    for candidate in (prediction.parent_id, *reversed(tuple(ancestors))):
        if candidate is not None and candidate in index:
            return index[candidate]
    return None
