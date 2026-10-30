"""Sequences of the antibody and antigen pairs that the Track A challenge could announce.

The research ranks them (docs/research/06): trastuzumab with HER2 first, because VH 120 + a 15-residue
(G4S)3 linker + VL 107 is exactly the announced 242 residues. Each candidate's scFv is built from the
public Fab chains of a PDB entry, VH first as in Adaptyv's layout. The antigen is the whole extracellular
domain from UniProt, so it covers any shorter construct. Nothing here is announced, it is a hedge.

    for candidate in CANDIDATES:
        info = resolve(candidate)      # scfv, antigen, lengths, whether the PDB antigen lies inside the domain
"""

from __future__ import annotations

import json
import re
import urllib.request
from dataclasses import dataclass

LINKER = "GGGGS" * 3
VH_END = re.compile(r"WG.G[TS][A-Z]{1,2}VTVS[SAP]")
"""The J region that closes a heavy variable domain, WGQGTLVTVSS in trastuzumab."""
VL_END = re.compile(r"FG.GT[KR][LVA][EDT][IL][KR]")
"""The J region that closes a kappa light variable domain, FGQGTKVEIK in pertuzumab."""


class CandidateError(RuntimeError):
    """A PDB entry or UniProt record that does not look the way this module assumes."""


@dataclass(frozen=True)
class Candidate:
    name: str
    pdb: str
    """Entry holding the Fab and the antigen; its heavy chain, light chain and antigen are separate entities."""
    uniprot: str
    """Accession of the antigen, whose extracellular domain is used."""


CANDIDATES = (
    Candidate("trastuzumab_her2", "1N8Z", "P04626"),
    Candidate("pertuzumab_her2", "1S78", "P04626"),
    Candidate("sacituzumab_trop2", "9PI9", "P09758"),
    Candidate("cetuximab_egfr", "1YY9", "P00533"),
    Candidate("pembrolizumab_pd1", "5GGS", "Q15116"),
    Candidate("atezolizumab_pdl1", "5X8L", "Q9NZQ7"),
)


def _get(url: str) -> str:
    with urllib.request.urlopen(url, timeout=30) as response:
        return response.read().decode()


def pdb_entities(pdb: str) -> list[tuple[str, str]]:
    """(title, sequence) for each entity of a PDB entry, from RCSB's FASTA."""
    entities = []
    for block in _get(f"https://www.rcsb.org/fasta/entry/{pdb}/display").split(">")[1:]:
        head, *lines = block.strip().split("\n")
        parts = head.split("|")
        entities.append((parts[2] if len(parts) > 2 else head, "".join(lines)))
    return entities


def classify(entities: list[tuple[str, str]], pdb: str = "") -> tuple[str, str, str]:
    """(heavy, light, antigen) from entity titles. Anything but one of each is an error, not a guess."""
    heavy = [s for t, s in entities if re.search(r"heavy", t, re.I)]
    light = [s for t, s in entities if re.search(r"light", t, re.I)]
    other = [s for t, s in entities if not re.search(r"heavy|light", t, re.I)]
    if (len(heavy), len(light), len(other)) != (1, 1, 1):
        titles = [t for t, _ in entities]
        raise CandidateError(f"{pdb}: expected one heavy chain, one light chain and one antigen, got titles {titles}")
    return heavy[0], light[0], other[0]


def variable_domain(chain: str, end: re.Pattern, what: str = "chain") -> str:
    """The chain up to and including the J region, which drops the constant domain."""
    found = end.search(chain)
    if found is None:
        raise CandidateError(f"no J region found in the {what}, cannot cut the variable domain")
    return chain[: found.end()]


def build_scfv(heavy: str, light: str) -> str:
    return variable_domain(heavy, VH_END, "heavy chain") + LINKER + variable_domain(light, VL_END, "light chain")


def uniprot_extracellular(accession: str) -> tuple[str, tuple[int, int]]:
    """The extracellular domain of a UniProt entry and its 1-based inclusive range."""
    record = json.loads(_get(f"https://rest.uniprot.org/uniprotkb/{accession}.json"))
    for feature in record.get("features", []):
        if feature.get("type") == "Topological domain" and feature.get("description", "").startswith("Extracellular"):
            start, end = feature["location"]["start"]["value"], feature["location"]["end"]["value"]
            return record["sequence"]["value"][start - 1 : end], (start, end)
    raise CandidateError(f"{accession}: no extracellular topological domain in the UniProt record")


def resolve(candidate: Candidate) -> dict:
    """Everything the cache build needs for one candidate, plus the checks that say if it can be trusted."""
    heavy, light, pdb_antigen = classify(pdb_entities(candidate.pdb), candidate.pdb)
    vh = variable_domain(heavy, VH_END, "heavy chain")
    vl = variable_domain(light, VL_END, "light chain")
    ecd, (start, end) = uniprot_extracellular(candidate.uniprot)
    return {
        "name": candidate.name,
        "pdb": candidate.pdb,
        "uniprot": candidate.uniprot,
        "scfv": vh + LINKER + vl,
        "vh_length": len(vh),
        "vl_length": len(vl),
        "antigen": ecd,
        "antigen_range": [start, end],
        "pdb_antigen_length": len(pdb_antigen),
        "pdb_antigen_inside_domain": pdb_antigen in ecd,
    }
