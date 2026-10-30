"""msa_cache.candidates: cutting variable domains, telling the entities of a PDB entry apart, resolving a candidate."""

import json

import pytest

from msa_cache import candidates
from msa_cache.candidates import (
    CANDIDATES,
    LINKER,
    VH_END,
    VL_END,
    Candidate,
    CandidateError,
    build_scfv,
    classify,
    resolve,
    uniprot_extracellular,
    variable_domain,
)

VH = "Q" * 100 + "WGQGTLVTVSS"
VL = "D" * 95 + "FGQGTKVEIK"
HEAVY = VH + "ASTKGPSVFPLAPSS"
LIGHT = VL + "RTVAAPSVFIFPP"
ANTIGEN = "ACDEFGHIKL"


def test_variable_domain_stops_after_the_j_region():
    assert variable_domain(HEAVY, VH_END) == VH and len(variable_domain(HEAVY, VH_END)) == 111
    assert variable_domain(LIGHT, VL_END) == VL and len(variable_domain(LIGHT, VL_END)) == 105


def test_cetuximab_style_j_regions_are_recognised():
    assert variable_domain("Q" * 50 + "WGAGTTVTVSA" + "AKTTPP", VH_END).endswith("WGAGTTVTVSA")
    assert variable_domain("D" * 50 + "FGAGTKLELK" + "RADAAP", VL_END).endswith("FGAGTKLELK")


def test_a_chain_without_a_j_region_is_an_error():
    with pytest.raises(CandidateError, match="heavy chain"):
        variable_domain("Q" * 80, VH_END, "heavy chain")


def test_scfv_is_vh_then_a_15_residue_linker_then_vl():
    scfv = build_scfv(HEAVY, LIGHT)
    assert scfv == VH + LINKER + VL and len(LINKER) == 15 and len(scfv) == 111 + 15 + 105


def test_classify_needs_exactly_one_of_each():
    entities = [("Receptor", ANTIGEN), ("Fab heavy chain", HEAVY), ("Fab light chain", LIGHT)]
    assert classify(entities) == (HEAVY, LIGHT, ANTIGEN)
    with pytest.raises(CandidateError, match="9XXX"):
        classify(entities[:2], "9XXX")
    with pytest.raises(CandidateError):
        classify(entities + [("another heavy chain", HEAVY)])


def uniprot_json(sequence="MABCDEFGHIJKLMNOP", extracellular=(3, 9)):
    feature = lambda kind, desc, a, b: {  # noqa: E731
        "type": kind, "description": desc, "location": {"start": {"value": a}, "end": {"value": b}}
    }
    return json.dumps(
        {
            "sequence": {"value": sequence},
            "features": [feature("Signal", "", 1, 2), feature("Topological domain", "Extracellular", *extracellular),
                         feature("Transmembrane", "Helical", 10, 12), feature("Topological domain", "Cytoplasmic", 13, 17)],
        }
    )


def test_uniprot_extracellular_picks_the_extracellular_segment(monkeypatch):
    monkeypatch.setattr(candidates, "_get", lambda url: uniprot_json())
    assert uniprot_extracellular("P00000") == ("BCDEFGH", (3, 9))  # residues 3 to 9 of MABCDEFGHIJKLMNOP


def test_uniprot_record_without_an_extracellular_domain_is_an_error(monkeypatch):
    monkeypatch.setattr(candidates, "_get", lambda url: json.dumps({"sequence": {"value": "MAAA"}, "features": []}))
    with pytest.raises(CandidateError, match="extracellular"):
        uniprot_extracellular("P00000")


def fake_get(uniprot: str):
    """RCSB gives a Fab and the 10-residue ANTIGEN; UniProt gives whatever `uniprot` holds."""

    def get(url: str) -> str:
        if "rcsb.org" in url:
            return (f">1XYZ_1|Chain C|Receptor protein|Homo sapiens (9606)\n{ANTIGEN}\n"
                    f">1XYZ_2|Chain B|Test Fab heavy chain|Homo sapiens (9606)\n{HEAVY}\n"
                    f">1XYZ_3|Chain A|Test Fab light chain|Homo sapiens (9606)\n{LIGHT}\n")
        return uniprot

    return get


def test_resolve_reports_whether_the_pdb_antigen_lies_inside_the_domain(monkeypatch):
    candidate = Candidate("test_antigen", "1XYZ", "P00000")
    # a 7-residue extracellular domain cannot hold the 10-residue antigen
    monkeypatch.setattr(candidates, "_get", fake_get(uniprot_json("M" + "XX" + ANTIGEN + "YYYY", (3, 9))))
    info = resolve(candidate)
    assert info["scfv"] == VH + LINKER + VL and info["vh_length"] == 111 and info["vl_length"] == 105
    assert info["antigen_range"] == [3, 9] and info["pdb_antigen_length"] == len(ANTIGEN)
    assert info["pdb_antigen_inside_domain"] is False

    monkeypatch.setattr(candidates, "_get", fake_get(uniprot_json("M" + ANTIGEN + "YY", (2, 11))))
    info = resolve(candidate)
    assert info["antigen"] == ANTIGEN and info["pdb_antigen_inside_domain"] is True


def test_candidate_names_and_entries_are_unique_and_ranked_trastuzumab_first():
    assert CANDIDATES[0].name == "trastuzumab_her2"
    assert len({c.name for c in CANDIDATES}) == len(CANDIDATES) == len({c.pdb for c in CANDIDATES})
