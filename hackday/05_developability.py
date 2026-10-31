"""Local, parent-relative risk screening for trastuzumab scFv substitutions.

Sequence metrics, heuristic patches on the *parent* Fv, and cached AntiFold
additive log-probability differences. No model inference, downloads or GPU.
These are descriptive proxies, not expression yields, melting temperatures,
validated TAP metrics or a calibrated pass/fail prediction of developability.
See hackday/DEVELOPABILITY.md for inputs, assumptions and output columns.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
from collections import defaultdict
from datetime import datetime, timezone
from importlib.metadata import PackageNotFoundError, version
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parent.parent
LINKER_LEN, VH_LEN, VL_LEN = 15, 120, 107
# Exact sequence supplied by the organisers; never infer the reference from a candidate.
PARENT = (
    "EVQLVESGGGLVQPGGSLRLSCAASGFNIKDTYIHWVRQAPGKGLEWVARIYPTNGYTRYADSVKGRFTISADTSKNT"
    "AYLQMNSLRAEDTAVYYCSRWGGDGFYAMDYWGQGTLVTVSSGGGGSGGGGSGGGGSDIQMTQSPSSLSASVGDR"
    "VTITCRASQDVNTAVAWYQQKPGKAPKLLIYSASFLYSGVPSRFSGSRSGTDFTLTISSLQPEDFATYYCQQHYTTP"
    "PTFGQGTKVEIK"
)
PARENT_NAMES = {"PARENT", "parent", "trastuzumab_scFv_parent", "trastuzumab_scfv_parent"}
CONSERVED_CYS = (22, 96, 158, 223)  # 1-based full-scFv numbering
AAS = "ACDEFGHIKLMNPQRSTVWY"
DEFAULT_STRUCTURE = ROOT / "data/structures/1N8Z.pdb"
DEFAULT_ANTIFOLD = ROOT / "data/scores/antifold/1N8Z_BAC.csv"

THREE_TO_ONE = {
    "ALA": "A", "ARG": "R", "ASN": "N", "ASP": "D", "CYS": "C", "GLN": "Q", "GLU": "E",
    "GLY": "G", "HIS": "H", "ILE": "I", "LEU": "L", "LYS": "K", "MET": "M", "PHE": "F",
    "PRO": "P", "SER": "S", "THR": "T", "TRP": "W", "TYR": "Y", "VAL": "V",
}

# Kyte-Doolittle, for GRAVY and for the hydrophobic patch term.
KD = {"A": 1.8, "R": -4.5, "N": -3.5, "D": -3.5, "C": 2.5, "Q": -3.5, "E": -3.5,
      "G": -0.4, "H": -3.2, "I": 4.5, "L": 3.8, "K": -3.9, "M": 1.9, "F": 2.8,
      "P": -1.6, "S": -0.8, "T": -0.7, "W": -0.9, "Y": -1.3, "V": 4.2}
CHARGE = {"D": -1.0, "E": -1.0, "K": 1.0, "R": 1.0, "H": 0.1}
PKA = {"D": 3.9, "E": 4.3, "C": 8.3, "Y": 10.1, "H": 6.0, "K": 10.5, "R": 12.5}
MAX_SASA = {"A": 121, "R": 265, "N": 187, "D": 187, "C": 148, "Q": 214, "E": 214,
            "G": 97, "H": 216, "I": 195, "L": 191, "K": 230, "M": 203, "F": 228,
            "P": 154, "S": 143, "T": 163, "W": 264, "Y": 255, "V": 165}

CDR_RANGES = {"VH": [(31, 35), (50, 65), (99, 109)], "VL": [(24, 34), (50, 56), (89, 97)]}


# ---------------------------------------------------------------- sequence block
def pi_of(seq: str) -> float:
    def charge(ph: float) -> float:
        q = 1 / (1 + 10 ** (ph - 9.7)) - 1 / (1 + 10 ** (2.35 - ph))
        for aa in "KRH":
            q += seq.count(aa) / (1 + 10 ** (ph - PKA[aa]))
        for aa in "DECY":
            q -= seq.count(aa) / (1 + 10 ** (PKA[aa] - ph))
        return q
    lo, hi = 0.0, 14.0
    for _ in range(100):
        mid = (lo + hi) / 2
        lo, hi = (mid, hi) if charge(mid) > 0 else (lo, mid)
    return (lo + hi) / 2


def net_charge(seq: str, ph: float = 7.4) -> float:
    q = 1 / (1 + 10 ** (ph - 9.7)) - 1 / (1 + 10 ** (2.35 - ph))
    for aa in "KRH":
        q += seq.count(aa) / (1 + 10 ** (ph - PKA[aa]))
    for aa in "DECY":
        q -= seq.count(aa) / (1 + 10 ** (PKA[aa] - ph))
    return q


def instability_index(seq: str) -> float:
    # Biopython uses all 400 dipeptide weights; do not fill missing pairs with 1.
    from Bio.SeqUtils.ProtParam import ProteinAnalysis
    return ProteinAnalysis(seq).instability_index()


def longest_run(seq: str, members: str) -> int:
    best = cur = 0
    for aa in seq:
        cur = cur + 1 if aa in members else 0
        best = max(best, cur)
    return best


def motifs(seq: str) -> set[str]:
    out = set()
    for i in range(len(seq) - 1):
        pair = seq[i:i + 2]
        if pair in ("NG", "NS", "NT"):
            out.add(f"deamid:{pair}@{i + 1}")
        if pair in ("DG", "DS", "DD"):
            out.add(f"isomer:{pair}@{i + 1}")
        if pair == "DP":
            out.add(f"frag:DP@{i + 1}")
    for i in range(len(seq) - 2):
        if seq[i] == "N" and seq[i + 1] != "P" and seq[i + 2] in "ST":
            out.add(f"glyc:N{i + 1}")
    return out


# ---------------------------------------------------------------- structure block
def load_fv(path: Path = DEFAULT_STRUCTURE):
    """Parent Fv from 1N8Z: VH (chain B 1-120) then VL (chain A 1-107), heavy atoms."""
    from biotite.structure.io.pdb import PDBFile
    s = PDBFile.read(path).get_structure(model=1)
    s = s[s.element != "H"]
    vh = s[(s.chain_id == "B") & (s.res_id >= 1) & (s.res_id <= VH_LEN)]
    vl = s[(s.chain_id == "A") & (s.res_id >= 1) & (s.res_id <= VL_LEN)]
    fv = vh + vl
    if np.any(fv.ins_code != "") or np.any(fv.hetero):
        raise ValueError("Expected standard residues without insertion codes in B1:120/A1:107")
    # scFv index per atom: VH p -> p, VL p -> p + VH_LEN + LINKER_LEN
    idx = np.where(fv.chain_id == "B", fv.res_id, fv.res_id + VH_LEN + LINKER_LEN)
    return fv, idx


def surface_table(fv, idx):
    """Per-residue: scFv index, side-chain centre, relative exposure in the free Fv."""
    from biotite.structure import sasa
    atom_sasa = sasa(fv, vdw_radii="Single")
    rows = []
    for i in np.unique(idx):
        m = idx == i
        res = fv[m]
        aa = THREE_TO_ONE.get(res.res_name[0])
        if aa is None or len(set(res.res_name)) != 1:
            raise ValueError(f"Ambiguous residue at scFv position {i}")
        if len(res[res.atom_name == "CA"]) != 1:
            raise ValueError(f"Expected one CA at scFv position {i}")
        side = res[~np.isin(res.atom_name, ["N", "CA", "C", "O"])]
        centre = side.coord.mean(axis=0) if len(side) else res[res.atom_name == "CA"].coord[0]
        rel = atom_sasa[m][~np.isin(res.atom_name, ["N", "C", "O"])].sum() / MAX_SASA[aa]
        rows.append({"scfv": int(i), "wt": aa, "centre": centre, "rel_sasa": float(rel)})
    expected = set(range(1, 121)) | set(range(136, 243))
    if {r["scfv"] for r in rows} != expected:
        raise ValueError("Reference structure must cover all 227 VH/VL residues")
    for row in rows:
        if row["wt"] != PARENT[row["scfv"] - 1]:
            raise ValueError(f"Structure/parent mismatch at scFv position {row['scfv']}")
        if not np.isfinite(row["rel_sasa"]) or not np.isfinite(row["centre"]).all():
            raise ValueError("Non-finite structure coordinates or exposure")
    return rows


def cdr_vicinity(surf) -> set[int]:
    """CDR residues plus every residue with a side-chain centre within 7.5 A of one."""
    cdr = set()
    for i in range(1, VH_LEN + 1):
        if any(lo <= i <= hi for lo, hi in CDR_RANGES["VH"]):
            cdr.add(i)
    for i in range(1, VL_LEN + 1):
        if any(lo <= i <= hi for lo, hi in CDR_RANGES["VL"]):
            cdr.add(i + VH_LEN + LINKER_LEN)
    centres = {r["scfv"]: r["centre"] for r in surf}
    near = set(cdr)
    for r in surf:
        for c in cdr:
            if c in centres and np.linalg.norm(r["centre"] - centres[c]) <= 7.5:
                near.add(r["scfv"])
                break
    return near


def patch_metrics(seq: str, surf, vicinity: set[int], radius: float = 7.5):
    """Largest hydrophobic, positive and negative patch in the CDR vicinity.

    A residue's patch score is the sum, over exposed neighbours within `radius`, of the
    neighbour's hydrophobicity (or charge), each weighted by its own exposure. The metric
    is the largest such score over the vicinity. These are custom heuristics,
    NOT the published TAP PSH/PPC/PNC definitions or calibrated clinical cutoffs.
    """
    pts = [r for r in surf if r["scfv"] in vicinity and r["rel_sasa"] >= 0.075]
    if not pts:
        return 0.0, 0.0, 0.0
    coords = np.array([p["centre"] for p in pts])
    aas = [seq[p["scfv"] - 1] for p in pts]
    expo = np.array([min(p["rel_sasa"], 1.0) for p in pts])
    d = np.linalg.norm(coords[:, None, :] - coords[None, :, :], axis=-1)
    near = (d <= radius)

    h = np.array([max(KD[a], 0.0) for a in aas]) * expo
    pos = np.array([max(CHARGE.get(a, 0.0), 0.0) for a in aas]) * expo
    neg = np.array([max(-CHARGE.get(a, 0.0), 0.0) for a in aas]) * expo
    return float((near @ h).max()), float((near @ pos).max()), float((near @ neg).max())


def fv_charge_symmetry(seq: str) -> float:
    vh, vl = seq[:VH_LEN], seq[VH_LEN + LINKER_LEN:]
    return net_charge(vh) * net_charge(vl)


def read_sequences(path: Path, *, csv_no_header: bool = False) -> dict[str, str]:
    """Accept FASTA or a CSV with name,sequence; never silently overwrite IDs."""
    if csv_no_header and path.suffix.lower() != ".csv":
        raise ValueError("--csv-no-header requires a .csv input")
    records: list[tuple[str, str]] = []
    with path.open(encoding="utf-8-sig", newline="") as fh:
        if csv_no_header:
            for row in csv.reader(fh):
                if not row or not any(cell.strip() for cell in row):
                    continue
                if len(row) != 2:
                    raise ValueError("Headerless sequence CSV requires exactly two columns: name,sequence")
                if [cell.strip().lower() for cell in row] == ["name", "sequence"]:
                    raise ValueError("CSV has a header; omit --csv-no-header")
                records.append((row[0].strip(), row[1]))
        elif path.suffix.lower() == ".csv":
            reader = csv.DictReader(fh)
            if not {"name", "sequence"} <= set(reader.fieldnames or []):
                raise ValueError("Sequence CSV requires name,sequence columns")
            for row in reader:
                if row.get("name") is None or row.get("sequence") is None or None in row:
                    raise ValueError("Malformed sequence CSV row")
                records.append((row["name"].strip(), row["sequence"]))
        else:
            name, parts = None, []
            for raw in fh:
                line = raw.strip()
                if not line:
                    continue
                if line.startswith(">"):
                    if name is not None:
                        records.append((name, "".join(parts)))
                    name, parts = line[1:].strip(), []
                elif name is None:
                    raise ValueError("FASTA sequence before the first header")
                else:
                    parts.append(line)
            if name is not None:
                records.append((name, "".join(parts)))
    if not records:
        raise ValueError("No sequences found")
    sequences: dict[str, str] = {}
    for name, seq in records:
        if not name or name in sequences:
            raise ValueError(f"Empty or duplicate sequence identifier: {name!r}")
        sequences[name] = "".join(seq.split()).upper()
    return sequences


def prepare_sequences(sequences: dict[str, str], parent_name: str) -> dict[str, str]:
    """Always include one verified reference row, even for candidate-only input."""
    if not parent_name.strip():
        raise ValueError("Parent identifier cannot be empty")
    candidates = {}
    for name, seq in sequences.items():
        if name in PARENT_NAMES | {parent_name}:
            if seq != PARENT:
                raise ValueError(f"Reference {name!r} does not match the official parent")
        else:
            candidates[name] = seq
    return {parent_name: PARENT, **candidates}


def read_antifold(path: Path) -> dict[int, dict[str, float]]:
    """Read once and verify reference identities and unique residue positions."""
    lp = {}
    with path.open(newline="") as fh:
        reader = csv.DictReader(fh)
        if not {"pdb_chain", "pdb_pos", "pdb_res", *AAS} <= set(reader.fieldnames or []):
            raise ValueError("AntiFold CSV lacks position, reference residue or AA columns")
        for row in reader:
            chain = row["pdb_chain"]
            if chain not in ("B", "A"):
                continue
            pos = int(row["pdb_pos"])
            if not 1 <= pos <= (VH_LEN if chain == "B" else VL_LEN):
                continue
            if row.get("pdb_posins", str(pos)) != str(pos):
                raise ValueError("AntiFold insertion-code numbering is unsupported")
            key = pos if chain == "B" else pos + VH_LEN + LINKER_LEN
            if key in lp or row["pdb_res"] != PARENT[key - 1]:
                raise ValueError(f"Duplicate or mismatched AntiFold residue at {key}")
            values = {aa: float(row[aa]) for aa in AAS}
            if not all(math.isfinite(v) for v in values.values()):
                raise ValueError(f"Non-finite AntiFold scores at {key}")
            lp[key] = values
    if not lp:
        raise ValueError("AntiFold CSV contains no mapped VH/VL residues")
    return lp


def antifold_delta(seq: str, parent: str, lp: dict | None) -> float | None:
    """Additive proxy on a fixed backbone, not a joint multi-mutant prediction."""
    if lp is None or len(seq) != len(parent):
        return None
    changes = [(i, a, b) for i, (a, b) in enumerate(zip(parent, seq), 1) if a != b]
    if any(i not in lp or a not in lp[i] or b not in lp[i] for i, a, b in changes):
        return None  # A partial sum would look like a complete score.
    return sum(lp[i][b] - lp[i][a] for i, a, b in changes)


def sequence_hash(seq: str) -> str:
    return hashlib.sha256(seq.encode()).hexdigest()


def read_ipsae(path: Path, sequences: dict[str, str]) -> dict[str, float]:
    """Explicit sequence match prevents joining scores from an older design version."""
    scores = {}
    with path.open(encoding="utf-8-sig", newline="") as fh:
        reader = csv.DictReader(fh)
        if not {"design", "sequence", "ipsae"} <= set(reader.fieldnames or []):
            raise ValueError("ipSAE CSV requires design,sequence,ipsae columns")
        for row in reader:
            name = (row.get("design") or "").strip()
            if name not in sequences or name in scores:
                raise ValueError(f"Unknown or duplicate ipSAE design: {name!r}")
            seq = "".join((row.get("sequence") or "").split()).upper()
            if seq != sequences[name]:
                raise ValueError(f"ipSAE sequence mismatch for {name!r}")
            value = float(row["ipsae"])
            if not math.isfinite(value) or not 0 <= value <= 1:
                raise ValueError(f"ipSAE must be finite and between 0 and 1: {name}")
            scores[name] = value
    return scores


METRIC_COLUMNS = (
    "pI", "d_pI", "charge_7.4", "d_charge", "gravy", "d_gravy", "instability",
    "d_instability", "hydrophobic_run", "d_hydrophobic_run", "PSH", "d_PSH",
    "PPC", "d_PPC", "PNC", "d_PNC", "SFvCSP", "d_SFvCSP", "antifold_d",
)


def liability_details(seq: str, new_motifs: list[str]) -> str:
    """Explain which changed positions intersect each newly detected motif."""
    details = []
    for motif in new_motifs:
        glyco = motif.startswith("glyc:")
        start = int(motif.split(":N")[1] if glyco else motif.split("@")[1])
        positions = range(start, start + (3 if glyco else 2))
        changes = [f"{PARENT[i-1]}{i}{seq[i-1]}" for i in positions if PARENT[i-1] != seq[i-1]]
        details.append(f"{motif} overlaps substitutions {', '.join(changes)} (full-scFv numbering)")
    return " | ".join(details)


def review_notes(row: dict) -> str:
    """Readable explanations; keep biological flags separate from missing evidence."""
    if row["input_status"] == "invalid":
        return f"Sequence not evaluated: {row['input_errors']}."
    if row["is_reference"]:
        return "Official parent reference; the five-substitution rule does not apply."
    notes = []
    rules = {
        "fewer_than_5_substitutions": "Fewer than five substitutions relative to the official parent.",
        "linker_changed": "The protected linker sequence has changed.",
        "conserved_cysteine_changed": "At least one protected cysteine has changed.",
        "length_change_requires_alignment_and_rule_review":
            "Length differs from the parent; positional checks require alignment and rule review.",
    }
    notes.extend(rules[issue] for issue in row["rule_issues"].split(";") if issue)
    flags = set(row["flags"].split(";"))
    if "new_liability" in flags:
        notes.append("New sequence liability motif(s): " + row["liability_details"] + ".")
    if "new_glycosylation_motif" in flags:
        notes.append("Potential glycosylation motif; interpretation depends on the expression system.")
    if "cys_changed" in flags:
        notes.append("Cysteine positions changed; check disulfide connectivity.")
    if "introduced_cysteine" in flags:
        notes.append("Additional cysteine at full-scFv position(s): " + row["introduced_cysteines"] + ".")
    if "pI_shift" in flags:
        notes.append(f"Estimated pI shifts by {row['d_pI']:+.2f}; this is not an expression failure prediction.")
    if "hydrophobic_run" in flags:
        notes.append("The longest hydrophobic stretch is longer than in the parent.")
    for flag, label in (("hydrophobic_patch", "hydrophobic"), ("positive_patch", "positive-charge"),
                        ("negative_patch", "negative-charge")):
        if flag in flags:
            notes.append(f"The {label} patch heuristic exceeds the parent-relative tolerance.")
    if "approx_buried" in flags:
        notes.append("A mutation maps to a buried position; fixed-parent structural estimates need caution.")
    if "duplicate_sequence" in flags:
        notes.append("Identical sequence also appears as: " + row["duplicate_sequence_ids"] +
                     ". Coordinate planned runs; duplicates may be intentional controls or seed replicates.")
    if row["data_warnings"]:
        notes.append("Some evidence was not evaluated; see data_warnings. Sequence checks remain usable before GPU runs.")
    if not notes:
        notes.append("No warnings from the implemented checks; not experimental validation.")
    return " ".join(notes)


def screen(sequences: dict[str, str], *, parent_name: str = "PARENT", surf=None,
           lp=None, tol: float = 1.0, ipsae=None, resource_warnings=()) -> list[dict]:
    """Keep every input row. Biological warnings do not silently discard candidates."""
    if not math.isfinite(tol) or tol < 0:
        raise ValueError("Patch tolerance must be finite and nonnegative")
    sequences = prepare_sequences(sequences, parent_name)
    sequence_ids = defaultdict(list)
    for name, seq in sequences.items():
        if name != parent_name:
            sequence_ids[seq].append(name)
    p_pi, p_q = pi_of(PARENT), net_charge(PARENT)
    p_gravy = sum(KD[a] for a in PARENT) / len(PARENT)
    p_run = longest_run(PARENT, "AVILMFWC")
    try:
        p_instability = instability_index(PARENT)
    except ImportError:
        p_instability = None
    vicinity = cdr_vicinity(surf) if surf is not None else None
    patches = patch_metrics(PARENT, surf, vicinity) if surf is not None else None
    buried = {r["scfv"] for r in surf if r["rel_sasa"] < 0.25} if surf else set()
    rows = []
    for name, seq in sequences.items():
        reference = name == parent_name
        row = dict(design=name, sequence=seq, sequence_sha256=sequence_hash(seq),
                   is_reference=reference, len=len(seq), input_status="valid",
                   input_errors="", rules_status="reference" if reference else "pass",
                   rule_issues="", n_mut=None, mutations="", linker_unchanged=None,
                   conserved_cysteines_unchanged=None, introduced_cysteines="")
        row.update(dict.fromkeys(METRIC_COLUMNS))
        row.update(instability_status="unavailable", structure_status="unavailable",
                   antifold_status="unavailable", ipsae=(ipsae or {}).get(name),
                   ipsae_status="available" if name in (ipsae or {}) else "not_provided",
                   new_motifs="", liability_details="", duplicate_sequence_ids="",
                   flags="", data_warnings="", assessment="", review_notes="")
        data_warnings = list(resource_warnings)
        invalid = sorted(set(seq) - set(AAS))
        if not seq or invalid:
            row.update(input_status="invalid", input_errors=("empty_sequence" if not seq
                       else "invalid_amino_acids:" + "".join(invalid)),
                       rules_status="not_evaluated", assessment="invalid_input")
            row["review_notes"] = review_notes(row)
            rows.append(row)
            continue

        flags, issues = [], []
        duplicates = [other for other in sequence_ids[seq] if other != name] if not reference else []
        if duplicates:
            flags.append("duplicate_sequence")
            row["duplicate_sequence_ids"] = ";".join(duplicates)
        same_length = len(seq) == len(PARENT)
        muts = []
        if same_length:
            muts = [i for i, (a, b) in enumerate(zip(PARENT, seq), 1) if a != b]
            row.update(n_mut=len(muts),
                       mutations=";".join(f"{PARENT[i-1]}{i}{seq[i-1]}" for i in muts),
                       linker_unchanged=seq[120:135] == PARENT[120:135],
                       conserved_cysteines_unchanged=all(seq[i-1] == "C" for i in CONSERVED_CYS))
            if not reference and len(muts) < 5:
                issues.append("fewer_than_5_substitutions")
            if not row["linker_unchanged"]:
                issues.append("linker_changed")
            if not row["conserved_cysteines_unchanged"]:
                issues.append("conserved_cysteine_changed")
            added_cys = [i for i in muts if seq[i-1] == "C"]
            row["introduced_cysteines"] = ";".join(map(str, added_cys))
            if added_cys:
                flags.append("introduced_cysteine")
            if any(PARENT[i-1] == "C" or seq[i-1] == "C" for i in muts):
                flags.append("cys_changed")
            new = sorted(motifs(seq) - motifs(PARENT))
            row["new_motifs"] = ";".join(new)
            row["liability_details"] = liability_details(seq, new)
            if new:
                flags.append("new_liability")
            if any(m.startswith("glyc:") for m in new):
                flags.append("new_glycosylation_motif")
        else:
            # Indel rules have not been agreed; do not invent a hard length limit.
            row["rules_status"] = "review_required"
            issues.append("length_change_requires_alignment_and_rule_review")
            data_warnings.append("positional_comparisons_not_evaluated")
        if issues:
            if same_length:
                row["rules_status"] = "fail"
            row["rule_issues"] = ";".join(issues)

        pi, q = pi_of(seq), net_charge(seq)
        gravy = sum(KD[a] for a in seq) / len(seq)
        run = longest_run(seq, "AVILMFWC")
        row.update(pI=round(pi, 2), d_pI=round(pi-p_pi, 2))
        row.update({"charge_7.4": round(q, 2), "d_charge": round(q-p_q, 2),
                    "gravy": round(gravy, 3), "d_gravy": round(gravy-p_gravy, 3),
                    "hydrophobic_run": run, "d_hydrophobic_run": run-p_run})
        if abs(pi-p_pi) >= 1.0:
            flags.append("pI_shift")
        if run > p_run:
            flags.append("hydrophobic_run")
        if p_instability is not None:
            instability = instability_index(seq)
            row.update(instability=round(instability, 2),
                       d_instability=round(instability-p_instability, 2),
                       instability_status="available_biopython")
        else:
            data_warnings.append("instability_unavailable_install_biopython")

        if same_length and surf is not None:
            values = patch_metrics(seq, surf, vicinity)
            for metric, value, parent_value, flag in zip(
                ("PSH", "PPC", "PNC"), values, patches,
                ("hydrophobic_patch", "positive_patch", "negative_patch")
            ):
                row[metric], row["d_" + metric] = round(value, 2), round(value-parent_value, 2)
                if value > parent_value + tol:
                    flags.append(flag)
            symmetry = fv_charge_symmetry(seq)
            row.update(SFvCSP=round(symmetry, 2),
                       d_SFvCSP=round(symmetry-fv_charge_symmetry(PARENT), 2),
                       structure_status="parent_backbone_heuristic")
            if set(muts) & buried:
                flags.append("approx_buried")
        else:
            data_warnings.append("structure_not_evaluated")
        af = antifold_delta(seq, PARENT, lp)
        if af is not None:
            row.update(antifold_d=round(af, 2), antifold_status="available_additive_proxy")
        else:
            data_warnings.append("antifold_unavailable_or_unmapped_mutations")
        row["flags"] = ";".join(flags)
        row["data_warnings"] = ";".join(dict.fromkeys(data_warnings))
        row["assessment"] = (
            "reference" if reference else "rule_violation" if row["rules_status"] == "fail"
            else "needs_review" if flags or data_warnings or issues
            else "no_warnings_detected"
        )
        row["review_notes"] = review_notes(row)
        rows.append(row)
    return rows


def file_info(path: Path | None) -> dict | None:
    if path is None:
        return None
    return {"path": str(path.resolve()), "sha256": hashlib.sha256(path.read_bytes()).hexdigest()
            if path.is_file() else None}


def attach_physics(rows: list[dict], path: Path) -> None:
    """Attach the companion CPU report using exact IDs AND sequence identity.

    Parent observations never populate a candidate's row. Missing structures and
    incomplete checks are explicit; neither changes a sequence-rule verdict.
    """
    from design_pipeline.physics import METRICS
    by_name = {row["design"]: row for row in rows}
    with path.open(newline="") as handle:
        reader = csv.DictReader(handle)
        required = {"design", "sequence", "sequence_sha256", "physics_status", "physics_scope",
                    "physics_flags", "physics_missing_checks", "physics_error", *METRICS}
        if not required <= set(reader.fieldnames or []):
            raise ValueError("Physics CSV must be produced by 06_physics_screen.py")
        imported = {}
        allowed_statuses = {"not_evaluated_no_structure", "invalid_structure", "review_required",
                            "partial", "no_geometry_warnings"}
        for record in reader:
            name = record["design"]
            if name in PARENT_NAMES and name not in by_name:
                name = next(row["design"] for row in rows if row["is_reference"])
            if name not in by_name or name in imported:
                raise ValueError(f"Unknown or duplicate physics ID: {name}")
            if record["sequence"] != by_name[name]["sequence"] or record["sequence_sha256"] != sequence_hash(record["sequence"]):
                raise ValueError(f"Physics sequence/hash mismatch for {name}")
            if record["physics_status"] not in allowed_statuses:
                raise ValueError(f"Unknown physics status for {name}")
            for metric in (*METRICS, *("d_"+m for m in METRICS)):
                value = record.get(metric)
                if value and not math.isfinite(float(value)):
                    raise ValueError(f"Non-finite physics metric {metric} for {name}")
            imported[name] = record
    for row in rows:
        record = imported.get(row["design"], {})
        row.update(physics_status=record.get("physics_status", "not_evaluated_no_structure"),
                   physics_scope=record.get("physics_scope", ""), physics_flags=record.get("physics_flags", ""),
                   physics_missing_checks=record.get("physics_missing_checks", "candidate_structure_not_provided"),
                   physics_error=record.get("physics_error", ""))
        for metric in (*METRICS, *("d_"+m for m in METRICS)):
            value = record.get(metric)
            row["physics_"+metric] = float(value) if value else None
        if not row["is_reference"]:
            warning = row["physics_status"] != "no_geometry_warnings"
            if warning:
                row["data_warnings"] = ";".join(filter(None, (row["data_warnings"], "physics_"+row["physics_status"])))
                if row["assessment"] == "no_warnings_detected":
                    row["assessment"] = "needs_review"
            row["review_notes"] += f" Physics geometry: {row['physics_status']}."
            if row["physics_flags"]:
                row["review_notes"] += " Geometry flags: " + row["physics_flags"] + "."
            if row["physics_missing_checks"]:
                row["review_notes"] += " Missing physics checks: " + row["physics_missing_checks"] + "."


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--input", "--fasta", dest="input", default="data/designs/variants.fasta")
    ap.add_argument("--csv-no-header", action="store_true",
                    help="input is a two-column name,sequence CSV without a header (exporter default)")
    ap.add_argument("--parent", default="PARENT", help="reference row label; sequence is fixed")
    ap.add_argument("--structure", default=str(DEFAULT_STRUCTURE), help="1N8Z numbering: B=VH, A=VL")
    ap.add_argument("--antifold", default=str(DEFAULT_ANTIFOLD), help="cached residue score CSV")
    ap.add_argument("--sequence-only", action="store_true", help="skip structure and AntiFold files")
    ap.add_argument("--ipsae", help="optional CSV: design,sequence,ipsae")
    ap.add_argument("--physics", help="optional physics_review.csv from 06_physics_screen.py")
    ap.add_argument("--tol", type=float, default=1.0, help="uncalibrated absolute patch tolerance")
    ap.add_argument("--out", default="data/designs/developability_review.csv")
    args = ap.parse_args(argv)
    if not math.isfinite(args.tol) or args.tol < 0:
        ap.error("--tol must be finite and nonnegative")
    input_path, output = ROOT / args.input, ROOT / args.out
    structure_path = None if args.sequence_only else ROOT / args.structure
    af_path = None if args.sequence_only else ROOT / args.antifold
    ipsae_path = ROOT / args.ipsae if args.ipsae else None
    physics_path = ROOT / args.physics if args.physics else None
    metadata_path = output.with_suffix(".metadata.json")
    # Avoid overwriting any input through a mistaken --out argument.
    sources = [input_path, structure_path, af_path, ipsae_path, physics_path, Path(__file__)]
    if output.resolve() == metadata_path.resolve() or any(
        target.resolve() == source.resolve()
        for target in (output, metadata_path) for source in sources if source is not None
    ):
        ap.error("Output and metadata paths must differ from each other and all inputs")
    try:
        sequences = prepare_sequences(read_sequences(input_path, csv_no_header=args.csv_no_header), args.parent)
        ipsae = read_ipsae(ipsae_path, sequences) if ipsae_path else {}
    except (OSError, ValueError, TypeError, csv.Error) as exc:
        ap.error(str(exc))
    surf, lp, warnings = None, None, []
    if not args.sequence_only:
        try:
            surf = surface_table(*load_fv(structure_path))
        except (OSError, ValueError, ImportError) as exc:
            warnings.append(f"structure_unavailable:{exc}")
        try:
            lp = read_antifold(af_path)
        except (OSError, ValueError, TypeError, csv.Error) as exc:
            warnings.append(f"antifold_unavailable:{exc}")
    rows = screen(sequences, parent_name=args.parent, surf=surf, lp=lp,
                  tol=args.tol, ipsae=ipsae, resource_warnings=warnings)
    if physics_path:
        try:
            attach_physics(rows, physics_path)
        except (OSError, ValueError, TypeError, csv.Error) as exc:
            ap.error(str(exc))
    versions = {}
    for package in ("numpy", "biotite", "biopython"):
        try:
            versions[package] = version(package)
        except PackageNotFoundError:
            versions[package] = None
    metadata = {
        "schema_version": 2, "created_utc": datetime.now(timezone.utc).isoformat(),
        "parent_sequence": PARENT, "parent_sha256": sequence_hash(PARENT),
        "input": file_info(input_path), "script": file_info(Path(__file__)),
        "structure": file_info(structure_path), "antifold": file_info(af_path),
        "ipsae": file_info(ipsae_path), "physics": file_info(physics_path), "package_versions": versions,
        "parameters": vars(args), "resource_warnings": warnings,
        "limitations": [
            "Risk screening only; no expression yield, melting temperature or folding ddG prediction.",
            "PSH/PPC/PNC/SFvCSP are custom parent-backbone heuristics, not validated TAP metrics.",
            "Patch tolerance is an absolute heuristic, not standard deviations or a clinical cutoff.",
            "No mutant side-chain rebuilding or linker structural modelling.",
            "AntiFold is an additive fixed-backbone proxy; multi-mutation interactions are not modelled.",
            "pI and charge are simple sequence estimates with fixed pKa values, not measured values.",
            "Instability index is descriptive; no universal threshold is used to reject scFvs.",
            "Missing ipSAE is separate from developability assessment; no ipSAE ranking is performed.",
            "Attached physics geometry is structure-specific, may be incomplete, and is not an energy or expression prediction.",
        ],
    }
    output.parent.mkdir(parents=True, exist_ok=True)
    with output.open("w", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)
    metadata_path.write_text(json.dumps(metadata, indent=2) + "\n")
    print(f"{'design':<28} {'rules':<16} {'assessment':<22} flags")
    for row in rows:
        print(f"{row['design']:<28} {row['rules_status']:<16} {row['assessment']:<22} {row['flags'] or '-'}")
    print(f"\nWrote {output}\nMetadata: {metadata_path}")
    print("Empty numeric cells mean unavailable. No warnings is not experimental validation.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
