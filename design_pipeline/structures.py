"""The structure a campaign starts from, and the few things we read off it.

Affinity maturation starts from a molecule that already exists, so step 1
*obtains* a structure rather than generating one: `Complex` is what
`get_structure` hands back, and every later step reads what it needs out of the
files it points at. Nothing is held in memory between steps.

Reading is deliberately thin: a chain's one-letter sequence, its residue ids,
and CA coordinates. Everything else a scorer wants, the scorer reads itself.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from pathlib import Path
from typing import Any, Mapping

import numpy as np

# Three-letter to one-letter, the 20 plus the modified residues that turn up in
# antibody crystal structures. Anything else reads as X: an unknown residue is
# a position we cannot mutate, not a reason to refuse the whole file.
AA3TO1: Mapping[str, str] = {
    "ALA": "A", "ARG": "R", "ASN": "N", "ASP": "D", "CYS": "C",
    "GLN": "Q", "GLU": "E", "GLY": "G", "HIS": "H", "ILE": "I",
    "LEU": "L", "LYS": "K", "MET": "M", "PHE": "F", "PRO": "P",
    "SER": "S", "THR": "T", "TRP": "W", "TYR": "Y", "VAL": "V",
    "MSE": "M", "HSD": "H", "HSE": "H", "HSP": "H",
    "CSO": "C", "CME": "C", "MLY": "K", "SEP": "S", "TPO": "T", "PTR": "Y",
}

AMINO_ACIDS = "ACDEFGHIKLMNPQRSTVWY"
"""The 20, one letter each, alphabetical. The default mutation set."""

CIF_SUFFIXES = (".cif", ".mmcif")
"""Read as text mmCIF. What the co-folding models write."""

PDB_SUFFIXES = (".pdb", ".ent")
"""Read as PDB. What `pdb_download` writes."""

STRUCTURE_SUFFIXES = PDB_SUFFIXES + CIF_SUFFIXES
"""Every suffix `load` can dispatch on. Anything that writes a structure file
names it with one of these, so that the name never disagrees with the bytes."""


@dataclass(frozen=True)
class Complex:
    """What step 1 returns: structure files on disk plus the resolved config.

    `paths` is the whole ensemble. A downloaded crystal structure is one file;
    a predictor asked for 20 seeds returns 20, because ranking of
    antibody-antigen complexes is unreliable and downstream code has to be able
    to see the spread rather than trust the top pick. `path` is the
    representative - the one scorers run on unless told otherwise - and is
    simply the first, which for `pdb_download` is the only one and for a
    predictor is the tool's own ranking.

    `output_dir` is the campaign directory the caller passed. Later steps
    default their own output under it, so a campaign stays in one tree.
    """

    paths: tuple[Path, ...]
    output_dir: Path
    config: Mapping[str, Any] = field(default_factory=dict)
    """Resolved: model_name, the kwargs as passed on, the cache key."""
    confidence: Mapping[str, Any] | None = None
    """How much to trust these coordinates, when they were predicted: the
    model-reported scores, the seed agreement, and the decoy control. Also
    written beside the structure as confidence.json. None when the caller asked
    for it to be skipped. See design_pipeline.confidence - and note the
    asymmetry recorded there: low ipTM reliably means do not trust this, while
    high ipTM does not reliably mean this is correct."""

    def __post_init__(self) -> None:
        if not self.paths:
            raise ValueError("a Complex needs at least one structure file")
        object.__setattr__(self, "paths", tuple(Path(p) for p in self.paths))
        object.__setattr__(self, "output_dir", Path(self.output_dir))

    @property
    def path(self) -> Path:
        """The representative structure. Pass `structure=` to a later step to
        score a different seed of the ensemble."""
        return self.paths[0]

    def __len__(self) -> int:
        return len(self.paths)


def load(structure: Path):
    """CA-and-all atoms of the first model, from PDB or mmCIF.

    The co-folding models write mmCIF and `pdb_download` writes PDB, so both
    are read here rather than converted at every call site.
    """
    structure = Path(structure)
    suffix = structure.suffix.lower()
    if suffix in CIF_SUFFIXES:
        from biotite.structure.io import pdbx

        return pdbx.get_structure(pdbx.CIFFile.read(str(structure)), model=1)
    if suffix in PDB_SUFFIXES:
        from biotite.structure.io.pdb import PDBFile

        return PDBFile.read(str(structure)).get_structure(model=1)
    # Rather than assume PDB and let biotite fail on the first line of an
    # mmCIF file, which surfaces as a parse error inside someone else's
    # library and sends you looking at the coordinates instead of at the name.
    raise ValueError(
        f"cannot tell the format of {structure.name} from its suffix; "
        f"expected one of {', '.join(STRUCTURE_SUFFIXES)}"
    )


def residue_starts(atoms, chain: str):
    """Index of the first atom of each residue of `chain`, in file order.

    Via biotite's own residue segmentation, which breaks on insertion codes as
    well as residue numbers - antibodies are full of 100A, 100B - so a CDR
    insertion counts as its own position instead of collapsing onto its
    neighbour.
    """
    import biotite.structure as struc

    sub = atoms[(atoms.chain_id == chain) & struc.filter_amino_acids(atoms)]
    if sub.array_length() == 0:
        present = ", ".join(sorted(set(atoms.chain_id))) or "none"
        raise KeyError(
            f"chain {chain!r} holds no amino acids in this structure; "
            f"chains present: {present}"
        )
    return sub, struc.get_residue_starts(sub)


def chain_ids(structure: Path) -> tuple[str, ...]:
    """Chain ids holding amino acids, in file order."""
    import biotite.structure as struc

    atoms = load(structure)
    atoms = atoms[struc.filter_amino_acids(atoms)]
    seen: list[str] = []
    for chain in atoms.chain_id:
        # str(), not the numpy scalar biotite hands back: these ids end up in
        # JSON, in error messages and in dict keys.
        label = str(chain)
        if label not in seen:
            seen.append(label)
    return tuple(seen)


def chain_sequence(structure: Path, chain: str) -> str:
    """One-letter sequence of `chain`, in file order.

    This is what `position` in the variant table indexes, 1-based: the chain as
    given, not the PDB residue numbers, which skip and repeat. IMGT labels are
    the cross-antibody key; this is the within-file one.
    """
    sub, starts = residue_starts(load(structure), chain)
    return "".join(AA3TO1.get(sub.res_name[i].upper(), "X") for i in starts)


def chain_residue_ids(structure: Path, chain: str) -> tuple[str, ...]:
    """PDB residue label per position, e.g. ("1", "2", "100A", ...).

    Needed whenever we have to speak to a tool in the file's own numbering -
    RFdiffusion contigs, a Rosetta resfile - rather than in chain positions.
    """
    sub, starts = residue_starts(load(structure), chain)
    return tuple(
        f"{sub.res_id[i]}{sub.ins_code[i].strip()}".strip() for i in starts
    )


def ca_coordinates(structure: Path, chain: str) -> np.ndarray:
    """(n_positions, 3) CA coordinates for `chain`, aligned with `chain_sequence`.

    A residue with no CA atom - a gap modelled as side-chain only, which does
    happen in crystal structures - gets NaN, and distance filters treat NaN as
    "cannot tell" rather than as zero.
    """
    sub, starts = residue_starts(load(structure), chain)
    bounds = list(starts) + [sub.array_length()]
    coords = np.full((len(starts), 3), np.nan)
    for i, (begin, end) in enumerate(zip(bounds, bounds[1:])):
        residue = sub[begin:end]
        ca = residue[residue.atom_name == "CA"]
        if ca.array_length():
            coords[i] = ca.coord[0]
    return coords


def read_fasta(path: Path) -> list[str]:
    """Every chain in a FASTA someone else wrote, in file order.

    Both conventions, and a mix of them: one record per chain, which is what
    you get out of a database or an alignment tool, and the one record with
    chains joined by ':' that `write_fasta` produces. A record is split on
    ':' and the records are concatenated, so

        >heavy          >complex
        EVQL...         EVQL...:DIQM...:TQVC...
        >light
        DIQM...

    both read as three chains. Order is what carries meaning downstream -
    `chain_ids` labels them positionally - so the record names are read only
    to name a record in an error.

    Validated rather than trusted, because the cost of a bad character here is
    a fold that runs for an hour and comes back wrong: a gapped alignment
    ('-'), a translated nucleotide file ('*') or a PDB handed over by mistake
    all fail at this line instead of somewhere inside the model.
    """
    path = Path(path)
    if not path.is_file():
        raise FileNotFoundError(f"no such FASTA: {path}")
    text = path.read_text()
    # [1:] drops whatever precedes the first header, which in a valid FASTA
    # is nothing and in a PDB handed over by mistake is the entire file - so
    # that one lands on the error below rather than parsing as one nameless
    # record.
    records = [r for r in text.split(">")[1:] if r.strip()]
    if not records:
        raise ValueError(
            f"{path.name} holds no FASTA records (nothing after a '>' header). "
            f"A structure file goes to source=, not to fasta=."
        )

    chains: list[str] = []
    for index, record in enumerate(records, start=1):
        header, _, body = record.partition("\n")
        name = header.split()[0] if header.split() else f"record {index}"
        # Joining on nothing rather than stripping line by line: a FASTA body
        # is wrapped at 60 or 80 columns, and the wrap is not a chain break.
        parts = [part.upper() for part in "".join(body.split()).split(":")]
        parts = [part for part in parts if part]
        if not parts:
            raise ValueError(f"{path.name}: record {name!r} has no sequence")
        for part in parts:
            unexpected = sorted({c for c in part if not c.isalpha()})
            if unexpected:
                raise ValueError(
                    f"{path.name}: record {name!r} holds "
                    f"{', '.join(repr(c) for c in unexpected)}, which is not a "
                    f"residue. Un-gap an alignment and drop any stop codon "
                    f"before folding it."
                )
        chains.extend(parts)
    return chains


def write_fasta(path: Path, name: str, chains: list[str]) -> Path:
    """One record, chains joined with ':'.

    LigandMPNN's convention, which ColabFold and the co-folding adapters in
    this package already read, so a sequence-only campaign can be handed to any
    of them without a second input format.
    """
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(f">{name}\n{':'.join(chains)}\n")
    return path


def relabel_chains(
    structure: Path,
    mapping: Mapping[str, str],
    out: Path,
    *,
    keep_only: bool = True,
) -> Path:
    """Write a copy of `structure` with its chains renamed.

    A deposited antibody complex calls its chains whatever the depositors
    chose - 1N8Z's heavy chain is B - while a prediction from this package
    calls them H, L and A. Every comparison in `metrics` matches chains by id
    across the two files, so one of them has to be relabelled, and relabelling
    the reference is the safe direction: the prediction is what gets scored.

    `keep_only` drops the chains nobody named, which is usually what you want -
    a crystal structure often holds a second copy of the complex, and leaving
    it in gives `interface_residues` a partner nobody asked about.
    """
    import biotite.structure as struc
    from biotite.structure.io.pdb import PDBFile

    atoms = load(structure)
    atoms = atoms[struc.filter_amino_acids(atoms)]
    if keep_only:
        atoms = atoms[np.isin(atoms.chain_id, list(mapping))]
        if atoms.array_length() == 0:
            raise KeyError(
                f"{structure.name} has none of chain(s) {', '.join(mapping)}"
            )
    renamed = atoms.chain_id.copy()
    for source, target in mapping.items():
        renamed[atoms.chain_id == source] = target
    atoms.chain_id = renamed

    out.parent.mkdir(parents=True, exist_ok=True)
    out_file = PDBFile()
    out_file.set_structure(atoms)
    out_file.write(str(out))
    return out
