"""CPU geometry screening of explicitly mapped antibody structures.

These are inspectable geometry flags, not force-field energies, a MolProbity
clashscore, or predictions of expression, affinity or thermodynamic stability.
No atom building, inference, downloads, or mutation of input files takes place.
"""
from __future__ import annotations

import hashlib
from functools import lru_cache
from pathlib import Path

import numpy as np
import biotite.structure as struc
from biotite.structure import info


PARAMETERS = {
    "clash_overlap_A": 0.6,
    "donor_acceptor_overlap_allowance_A": 0.4,
    "salt_bridge_max_A": 4.0,
    "polar_contact_max_A": 3.5,
    "hbond_H_acceptor_max_A": 2.5,
    "hbond_angle_min_deg": 120.0,
    "peptide_min_A": 1.1,
    "peptide_max_A": 1.6,
    "disulfide_min_A": 1.8,
    "disulfide_max_A": 2.3,
    "disulfide_abs_dihedral_min_deg": 45.0,
    "disulfide_abs_dihedral_max_deg": 135.0,
    "sasa_probe_A": 1.4,
    "sasa_points": 1000,
}
THREE_TO_ONE = dict(zip(
    "ALA ARG ASN ASP CYS GLN GLU GLY HIS ILE LEU LYS MET PHE PRO SER THR TRP TYR VAL".split(),
    "ARNDCQEGHILKMFPSTWYV",
))
CHARGED = {"ASP": (-1, {"OD1", "OD2"}), "GLU": (-1, {"OE1", "OE2"}),
           "LYS": (1, {"NZ"}), "ARG": (1, {"NE", "NH1", "NH2"})}
DONORS = {"ARG": {"NE", "NH1", "NH2"}, "LYS": {"NZ"}, "ASN": {"ND2"},
          "GLN": {"NE2"}, "TRP": {"NE1"}, "SER": {"OG"}, "THR": {"OG1"}, "TYR": {"OH"}}
ACCEPTORS = {"ASP": {"OD1", "OD2"}, "GLU": {"OE1", "OE2"}, "ASN": {"OD1"},
             "GLN": {"OE1"}, "SER": {"OG"}, "THR": {"OG1"}, "TYR": {"OH"}}
METRICS = (
    "clashes_binder", "clashes_interface", "peptide_outliers_binder",
    "disulfide_geometry_outliers", "disulfides_unassessed",
    "salt_bridges_binder", "salt_bridges_interface",
    "like_charge_contacts_binder", "like_charge_contacts_interface",
    "polar_contacts_binder", "polar_contacts_interface",
    "hydrogen_bonds_binder", "hydrogen_bonds_interface",
    "binder_sasa_A2", "binder_hydrophobic_sidechain_sasa_A2",
    "binder_buried_on_binding_A2", "missing_heavy_atoms_binder",
)
EVENT_FIELDS = ("kind", "scope", "atom1", "scfv_position1", "atom2", "scfv_position2",
                "distance_A", "overlap_A", "angle_deg", "note")


def sequence_hash(sequence: str) -> str:
    return hashlib.sha256(sequence.encode()).hexdigest()


def load_structure(path: Path, model: int = 1):
    """Select one model and the highest-occupancy alternate conformer."""
    if type(model) is not int or model < 1:
        raise ValueError("model must be a positive integer")
    if path.suffix.lower() == ".pdb":
        from biotite.structure.io.pdb import PDBFile
        return PDBFile.read(path).get_structure(model=model, altloc="occupancy")
    if path.suffix.lower() in {".cif", ".mmcif", ".bcif"}:
        from biotite.structure.io import pdbx
        cls = pdbx.BinaryCIFFile if path.suffix.lower() == ".bcif" else pdbx.CIFFile
        # The manifest always refers to author chain and residue identifiers.
        return pdbx.get_structure(cls.read(path), model=model, altloc="occupancy", use_author_fields=True)
    raise ValueError("Structure must be PDB, CIF/mmCIF or BCIF")


def _segment(segment, binder):
    expected = {"chain", "start", "end"} | ({"sequence_start"} if binder else set())
    if set(segment) != expected:
        raise ValueError(f"Segment needs exactly {sorted(expected)}")
    if not isinstance(segment["chain"], str) or not segment["chain"]:
        raise ValueError("Segment chain must be a nonempty string")
    for key in expected - {"chain"}:
        if type(segment[key]) is not int:
            raise ValueError(f"Segment {key} must be an integer")
    if segment["end"] < segment["start"]:
        raise ValueError("Segment end precedes start")
    if binder and segment["sequence_start"] < 1:
        raise ValueError("sequence_start must be positive")


def select_structure(atoms, record, sequence):
    """Validate every mapped binder residue, without assuming a chain order.

    Only substitution-only 242-aa scFvs are supported. Fv-only selections omit
    positions 121-135 explicitly. Missing binder residues are not silently scored.
    Target gaps are retained in the coverage warning and comparison signature.
    """
    scope = record["scope"]
    if scope not in {"fv_only", "full_scfv"} or len(sequence) != 242 or set(sequence) - set(THREE_TO_ONE.values()):
        raise ValueError("Expected scope fv_only/full_scfv and a 242-aa candidate sequence")
    wanted = set(range(1, 243))
    if scope == "fv_only":
        wanted -= set(range(121, 136))
    selected = np.zeros(len(atoms), bool)
    binder_mask = np.zeros(len(atoms), bool)
    positions = np.zeros(len(atoms), int)
    mapped, target_tokens, gaps = set(), [], []
    binder_segments, antigen_segments = record["binder_segments"], record.get("antigen_segments", [])
    if not isinstance(binder_segments, list) or not binder_segments or not isinstance(antigen_segments, list):
        raise ValueError("binder_segments must be a nonempty list; antigen_segments must be a list")
    for is_binder, segments in ((True, binder_segments), (False, antigen_segments)):
        for seg in segments:
            _segment(seg, is_binder)
            tokens = []
            for res_id in range(seg["start"], seg["end"] + 1):
                mask = (atoms.chain_id == seg["chain"]) & (atoms.res_id == res_id)
                if np.any(selected & mask):
                    raise ValueError("Overlapping structure segments")
                if np.any(atoms.ins_code[mask] != ""):
                    raise ValueError("Insertion codes require explicit remapping before screening")
                if not mask.any():
                    if is_binder:
                        raise ValueError(f"Missing binder residue {seg['chain']}:{res_id}")
                    gaps.append(f"{seg['chain']}:{res_id}")
                    tokens.append("?")
                    continue
                residues = set(atoms.res_name[mask])
                if len(residues) != 1 or next(iter(residues)) not in THREE_TO_ONE:
                    raise ValueError(f"Nonstandard/ambiguous protein residue {seg['chain']}:{res_id}")
                aa = THREE_TO_ONE[next(iter(residues))]
                if is_binder:
                    pos = seg["sequence_start"] + res_id - seg["start"]
                    if pos not in wanted or pos in mapped:
                        raise ValueError("Duplicate or out-of-scope scFv position")
                    if sequence[pos-1] != aa:
                        raise ValueError(f"Structure/sequence mismatch at scFv position {pos}")
                    mapped.add(pos)
                    positions[mask] = pos
                    binder_mask[mask] = True
                else:
                    tokens.append(aa)
                selected |= mask
            if not is_binder:
                target_tokens.append("".join(tokens))
    if mapped != wanted:
        raise ValueError("Binder segments do not cover exactly the declared scFv/Fv scope")
    if set(atoms.chain_id[binder_mask]) & set(atoms.chain_id[selected & ~binder_mask]):
        raise ValueError("Binder and antigen must have separate chains")
    if antigen_segments and not np.any(selected & ~binder_mask):
        raise ValueError("Antigen selection has no coordinates")
    selected_atoms = atoms[selected]
    if not np.isfinite(selected_atoms.coord).all():
        raise ValueError("Non-finite selected coordinates")
    keys = list(zip(selected_atoms.chain_id, selected_atoms.res_id, selected_atoms.atom_name))
    if len(set(keys)) != len(keys):
        raise ValueError("Duplicate selected atom identities")
    return selected_atoms, binder_mask[selected], positions[selected], {
        "scope": scope, "antigen_missing_residues": gaps,
        "antigen_signature": sequence_hash("|".join(target_tokens)) if target_tokens else None,
        "mapped_binder_residues": len(mapped),
    }


@lru_cache(None)
def _template(name):
    template = info.residue(name)
    return {str(a.atom_name): str(a.element) for a in template if a.element not in {"H", "D"}}


def _residues(atoms):
    groups = {}
    for i, (chain, number) in enumerate(zip(atoms.chain_id, atoms.res_id)):
        groups.setdefault((str(chain), int(number)), []).append(i)
    return list(groups.values())


def _types(atoms):
    donors, acceptors = set(), set()
    for i, atom in enumerate(atoms):
        if (atom.atom_name == "N" and atom.res_name != "PRO") or atom.atom_name in DONORS.get(atom.res_name, set()):
            donors.add(i)
        if atom.atom_name in {"O", "OXT"} or atom.atom_name in ACCEPTORS.get(atom.res_name, set()):
            acceptors.add(i)
    return donors, acceptors


def _neighbours(atoms, radius):
    cells = struc.CellList(atoms, cell_size=radius)
    neighbours = cells.get_atoms(atoms.coord, radius=radius)
    for i, row in enumerate(neighbours):
        for j in row:
            if j > i:
                yield i, int(j), float(np.linalg.norm(atoms.coord[i]-atoms.coord[j]))


def analyze_geometry(atoms, binder_mask, positions):
    """Analyze selected protein atoms; mask/positions come from select_structure()."""
    all_atoms, all_binder, all_positions = atoms, binder_mask, positions
    heavy = ~np.isin(atoms.element, ["H", "D"])
    atoms, binder, positions = atoms[heavy], binder_mask[heavy], positions[heavy]
    if not len(atoms) or not binder.any() or not np.isfinite(atoms.coord).all():
        raise ValueError("Finite heavy-atom binder coordinates are required")
    groups = _residues(atoms)
    events, flags, missing = [], [], []
    metrics = {name: 0 for name in METRICS}
    has_target = bool((~binder).any())
    for name in METRICS:
        if ("interface" in name or "on_binding" in name) and not has_target:
            metrics[name] = None

    def label(i):
        return f"{atoms.chain_id[i]}:{atoms.res_id[i]}:{atoms.res_name[i]}:{atoms.atom_name[i]}"

    def scope(i, j):
        return "binder" if binder[i] and binder[j] else "interface" if binder[i] != binder[j] else "antigen"

    def event(kind, i, j=None, distance=None, overlap=None, angle=None, note=""):
        events.append(dict(zip(EVENT_FIELDS, [kind, "binder" if j is None and binder[i] else
                      "antigen" if j is None else scope(i, j), label(i), int(positions[i]) or None,
                      label(j) if j is not None else None, int(positions[j]) or None if j is not None else None,
                      distance, overlap, angle, note])))

    missing_binder, missing_target = 0, 0
    for inds in groups:
        i = inds[0]
        template = _template(str(atoms.res_name[i]))
        names = set(atoms.atom_name[inds])
        for j in inds:
            if atoms.atom_name[j] not in template or atoms.element[j] != template[atoms.atom_name[j]]:
                raise ValueError(f"Unknown atom or element mismatch: {label(j)}")
        absent = sorted(set(template) - names - {"OXT"})
        if absent:
            if binder[i]:
                missing_binder += len(absent)
            else:
                missing_target += len(absent)
            event("missing_heavy_atoms", i, note=",".join(absent))
    metrics["missing_heavy_atoms_binder"] = missing_binder
    if missing_binder or missing_target:
        missing.append("incomplete_heavy_atoms")

    # Topological exclusions: CCD intra-residue bonds, explicit peptide bonds,
    # and geometrically plausible disulfides. Exclude up to three bond steps.
    bonds = struc.connect_via_residue_names(atoms, inter_residue=False).as_array()
    adjacency = [set() for _ in atoms]

    def add_bond(i, j):
        adjacency[i].add(j)
        adjacency[j].add(i)

    for i, j, _ in bonds:
        add_bond(int(i), int(j))
    lookup = {(str(atoms.chain_id[g[0]]), int(atoms.res_id[g[0]])):
              {str(atoms.atom_name[i]): i for i in g} for g in groups}
    for (chain, number), names in lookup.items():
        next_names = lookup.get((chain, number+1), {})
        if "C" in names and "N" in next_names:
            i, j = names["C"], next_names["N"]
            add_bond(i, j)
            distance = float(np.linalg.norm(atoms.coord[i]-atoms.coord[j]))
            if not PARAMETERS["peptide_min_A"] <= distance <= PARAMETERS["peptide_max_A"]:
                event("peptide_distance_outlier", i, j, distance)
                if binder[i] and binder[j]:
                    metrics["peptide_outliers_binder"] += 1
    sulfurs = np.flatnonzero((atoms.res_name == "CYS") & (atoms.atom_name == "SG"))
    for n, i in enumerate(sulfurs):
        for j in sulfurs[n+1:]:
            distance = float(np.linalg.norm(atoms.coord[i]-atoms.coord[j]))
            if PARAMETERS["disulfide_min_A"] <= distance <= PARAMETERS["disulfide_max_A"]:
                add_bond(int(i), int(j))
    exclusions = []
    for i in range(len(atoms)):
        seen, frontier = {i}, {i}
        for _ in range(3):
            frontier = {j for k in frontier for j in adjacency[k]} - seen
            seen |= frontier
        exclusions.append(seen)

    # Assess the two canonical scFv disulfides, including a broad torsion check.
    for first, second in ((22, 96), (158, 223)):
        a = np.flatnonzero(positions == first)
        b = np.flatnonzero(positions == second)
        names_a, names_b = {atoms.atom_name[i]: i for i in a}, {atoms.atom_name[i]: i for i in b}
        if not all(n in names_a and n in names_b for n in ("SG", "CB")):
            metrics["disulfides_unassessed"] += 1
            if len(a):
                event("disulfide_unassessed", int(a[0]), int(b[0]) if len(b) else None,
                      note="Requires CYS CB/SG at both canonical positions")
            continue
        i, j = names_a["SG"], names_b["SG"]
        distance = float(np.linalg.norm(atoms.coord[i]-atoms.coord[j]))
        angle = float(np.degrees(struc.dihedral(atoms.coord[names_a["CB"]], atoms.coord[i],
                                             atoms.coord[j], atoms.coord[names_b["CB"]])))
        good = (np.isfinite(angle) and PARAMETERS["disulfide_min_A"] <= distance <= PARAMETERS["disulfide_max_A"]
                and PARAMETERS["disulfide_abs_dihedral_min_deg"] <= abs(angle) <= PARAMETERS["disulfide_abs_dihedral_max_deg"])
        event("canonical_disulfide", i, j, distance, angle=angle if np.isfinite(angle) else None,
              note="within_screening_window" if good else "geometry_review_required")
        metrics["disulfide_geometry_outliers"] += int(not good)
    if metrics["disulfides_unassessed"]:
        missing.append("canonical_disulfide_geometry_incomplete")

    donors, acceptors = _types(atoms)
    radii = np.array([info.vdw_radius_single(str(e)) for e in atoms.element], dtype=float)
    if not np.isfinite(radii).all():
        raise ValueError("Unknown van der Waals radius")
    for i, j, distance in _neighbours(atoms, 4.0):
        where = scope(i, j)
        if where == "antigen" or j in exclusions[i]:
            continue
        polar = (i in donors and j in acceptors) or (j in donors and i in acceptors)
        overlap = float(radii[i]+radii[j]-distance)
        adjusted = overlap - (PARAMETERS["donor_acceptor_overlap_allowance_A"] if polar else 0)
        if adjusted >= PARAMETERS["clash_overlap_A"]:
            metrics[f"clashes_{where}"] += 1
            event("heavy_atom_clash", i, j, distance, overlap=overlap,
                  note=f"adjusted_overlap_A={adjusted:.4f}; not a MolProbity clashscore")
        if polar and distance <= PARAMETERS["polar_contact_max_A"]:
            metrics[f"polar_contacts_{where}"] += 1
            event("polar_distance_contact", i, j, distance, note="Distance only; not a confirmed hydrogen bond")

    # One contact per charged RESIDUE pair; do not count every O-N atom pair.
    charged = []
    for inds in groups:
        name = str(atoms.res_name[inds[0]])
        if name in CHARGED:
            charge, tip_names = CHARGED[name]
            tips = [i for i in inds if atoms.atom_name[i] in tip_names]
            if tips:
                charged.append((charge, tips))
    for n, (q1, tips1) in enumerate(charged):
        for q2, tips2 in charged[n+1:]:
            where = scope(tips1[0], tips2[0])
            if where == "antigen":
                continue
            distances = np.linalg.norm(atoms.coord[tips1, None, :] - atoms.coord[tips2][None, :, :], axis=-1)
            first, second = np.unravel_index(np.argmin(distances), distances.shape)
            distance = float(distances[first, second])
            if distance <= PARAMETERS["salt_bridge_max_A"]:
                opposite = q1 != q2
                metric = "salt_bridges" if opposite else "like_charge_contacts"
                metrics[f"{metric}_{where}"] += 1
                event("salt_bridge_candidate" if opposite else "like_charge_contact", tips1[first], tips2[second], distance,
                      note="Assumed Asp/Glu negative, Lys/Arg positive; His excluded; no energy inferred")

    # Hydrogen-aware geometry is only reported when hydrogens actually exist.
    original_indices = np.flatnonzero(heavy)
    excluded_pairs = {(int(original_indices[i]), int(original_indices[j]))
                      for i, excluded in enumerate(exclusions) for j in excluded}
    h_status, h_counts = _hydrogen_bonds(all_atoms, all_binder, all_positions, events, excluded_pairs)
    for where in ("binder", "interface"):
        metrics[f"hydrogen_bonds_{where}"] = h_counts[where] if (where == "binder" or has_target) else None
    if h_status != "available_explicit_hydrogens":
        missing.append(h_status)

    free_sasa = struc.sasa(atoms[binder], vdw_radii="Single", probe_radius=PARAMETERS["sasa_probe_A"],
                           point_number=PARAMETERS["sasa_points"])
    if not np.isfinite(free_sasa).all():
        raise ValueError("Non-finite SASA")
    metrics["binder_sasa_A2"] = float(free_sasa.sum())
    free = atoms[binder]
    hydrophobic = np.isin(free.res_name, ["ALA", "VAL", "ILE", "LEU", "MET", "PHE", "TRP", "TYR", "CYS"])
    hydrophobic &= ~np.isin(free.atom_name, ["N", "CA", "C", "O", "OXT"])
    metrics["binder_hydrophobic_sidechain_sasa_A2"] = float(free_sasa[hydrophobic].sum())
    if has_target:
        bound = struc.sasa(atoms, vdw_radii="Single", probe_radius=PARAMETERS["sasa_probe_A"],
                           point_number=PARAMETERS["sasa_points"])
        if not np.isfinite(bound).all():
            raise ValueError("Non-finite complex SASA")
        metrics["binder_buried_on_binding_A2"] = float(free_sasa.sum()-bound[binder].sum())
    for name in ("clashes_binder", "clashes_interface", "peptide_outliers_binder", "disulfide_geometry_outliers"):
        if metrics[name]:
            flags.append(name)
    return {"metrics": metrics, "flags": flags, "missing_checks": missing,
            "hydrogen_bond_status": h_status, "events": events,
            "heavy_atoms_binder": int(binder.sum()), "missing_heavy_atoms_antigen": missing_target,
            "histidine_residues_excluded_from_charge_analysis": sum(atoms.res_name[g[0]] == "HIS" for g in groups)}


def _hydrogen_bonds(atoms, binder, positions, events, excluded_pairs=frozenset()):
    donors, acceptors = _types(atoms)
    hydrogens = np.flatnonzero(np.isin(atoms.element, ["H", "D"]))
    if not len(hydrogens):
        return "unavailable_no_explicit_hydrogens", {"binder": None, "interface": None}
    counts = {"binder": 0, "interface": 0}
    covered = set()
    for d in donors:
        # Same residue and a covalent-sized distance; atom-name chemistry is
        # already restricted. Do not infer missing hydrogen coordinates.
        hs = [h for h in hydrogens if atoms.chain_id[h] == atoms.chain_id[d]
              and atoms.res_id[h] == atoms.res_id[d]
              and 0.7 <= np.linalg.norm(atoms.coord[h]-atoms.coord[d]) <= 1.25]
        if hs:
            covered.add(d)
        for a in acceptors:
            if d == a or (d, a) in excluded_pairs or (atoms.chain_id[d] == atoms.chain_id[a] and atoms.res_id[d] == atoms.res_id[a]):
                continue
            if not binder[d] and not binder[a]:
                continue
            for h in hs:
                distance = float(np.linalg.norm(atoms.coord[h]-atoms.coord[a]))
                angle = float(np.degrees(struc.angle(atoms.coord[d], atoms.coord[h], atoms.coord[a])))
                if distance <= PARAMETERS["hbond_H_acceptor_max_A"] and angle >= PARAMETERS["hbond_angle_min_deg"]:
                    where = "binder" if binder[d] and binder[a] else "interface"
                    counts[where] += 1
                    labels = [f"{atoms.chain_id[i]}:{atoms.res_id[i]}:{atoms.res_name[i]}:{atoms.atom_name[i]}" for i in (d, a)]
                    events.append(dict(zip(EVENT_FIELDS, ["hydrogen_bond_geometry", where, labels[0], int(positions[d]) or None,
                                  labels[1], int(positions[a]) or None, distance, None, angle,
                                  f"H={atoms.atom_name[h]}; restricted donor/acceptor chemistry"])))
    status = "available_explicit_hydrogens" if covered == donors else "partial_explicit_hydrogens"
    return status, counts
