"""Crystal chain positions and the single-chain Fv they belong to.

The pipeline scores, gates and combines substitutions in the numbering of a crystal structure. A Fab has a heavy and a
light chain, a variant id such as `B:S97A` or `A:Q27E` counts 1-based along the chain as given, and the constant
domains (CH1, CL) are part of those chains. The sequence that is folded and submitted is an scFv: one chain, the
variable heavy domain, a linker, the variable light domain. The two numberings agree on the first domain and differ by a
fixed offset on the second. For the trastuzumab scFv (VH 120, linker 15, VL 107) a light-chain position p is scFv
position p + 135, and the heavy chain keeps its numbers.

    layout = layout_from_parent(PARENT_SCFV)                  # VH first; heavy chain "B", light chain "A" as in 1N8Z
    verify_parent(layout, "1N8Z.cif")                         # the crystal's V domains are the parent's
    variant = apply_id(layout, "B:S97A+B:G103Y+A:Q27E+...")   # the scFv with those substitutions
    crystal_mutations(layout, variant)                        # the other way: an scFv's differences as chain positions

Nothing here runs a model. Every conversion checks the wild-type residue against the parent, so a shifted number is an
error and not a different residue. A position in a constant domain or in the linker has no counterpart and raises: a design
can only change the variable domains, and `enumerate_variants` scans the whole crystal chain, constant domain included.
The linker is left as it is, so a substitution inside it is refused in both directions.
"""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
from typing import Iterable, Mapping

from .structures import AMINO_ACIDS, chain_sequence
from .variants import Mutation, parse_id

LINKER = "GGGGSGGGGSGGGGS"
"""The 15-residue (G4S)3 linker of the trastuzumab scFv."""


@dataclass(frozen=True)
class ScfvLayout:
    """Where each crystal chain's variable domain sits in the scFv parent.

    `offsets[chain]` is how many scFv residues come before the chain's first one, so crystal position p of that chain is
    scFv position `offsets[chain] + p`. `lengths[chain]` is the length of its variable domain, which bounds p.
    """

    scfv: str
    offsets: Mapping[str, int]
    lengths: Mapping[str, int]
    linker: tuple[int, int]
    """(0-based start, length) of the linker inside `scfv`."""

    def domain(self, chain: str) -> str:
        start = self.offsets[chain]
        return self.scfv[start : start + self.lengths[chain]]


def layout_from_parent(scfv: str, first: str = "B", second: str = "A", linker: str = LINKER) -> ScfvLayout:
    """The layout of an scFv parent: `first` is the chain of the domain at the start, `second` the one after the linker.

    The defaults are the VH-first orientation of the trastuzumab scFv with the heavy chain as chain B and the light chain
    as chain A, as in 1N8Z. For a VL-first scFv pass `first="A", second="B"`. The linker has to occur exactly once. A
    linker that is longer than the one given leaves its extra residues in the second domain, which `verify_parent` catches.
    """
    unknown = sorted(set(scfv) - set(AMINO_ACIDS))
    if unknown:
        raise ValueError(f"the scFv has characters that are not amino acids: {unknown}")
    if first == second:
        raise ValueError(f"the two domains need two different chains, got {first!r} twice")
    found = scfv.count(linker)
    if found != 1:
        raise ValueError(f"expected the linker {linker} exactly once in the scFv, found it {found} times")
    start = scfv.index(linker)
    head, tail = scfv[:start], scfv[start + len(linker) :]
    if not head or not tail:
        raise ValueError("the linker has to sit between two domains, but one side of it is empty")
    return ScfvLayout(
        scfv=scfv,
        offsets={first: 0, second: len(head) + len(linker)},
        lengths={first: len(head), second: len(tail)},
        linker=(start, len(linker)),
    )


def scfv_position(layout: ScfvLayout, chain: str, position: int) -> int:
    """The 1-based scFv position of crystal `chain` position `position`. Raises outside the variable domain."""
    if chain not in layout.offsets:
        raise ValueError(f"chain {chain!r} is not in the layout, which has {sorted(layout.offsets)}")
    length = layout.lengths[chain]
    if not 1 <= position <= length:
        raise ValueError(
            f"{chain}:{position} is outside the variable domain of chain {chain} (positions 1 to {length}). That is a "
            f"constant domain or a position that does not exist; the scFv has no counterpart for it."
        )
    return layout.offsets[chain] + position


def chain_position(layout: ScfvLayout, index: int) -> tuple[str, int]:
    """(crystal chain, position) of the 1-based scFv position `index`. Raises in the linker."""
    if not 1 <= index <= len(layout.scfv):
        raise ValueError(f"scFv position {index} is outside the sequence, which has {len(layout.scfv)} residues")
    for chain, offset in layout.offsets.items():
        if offset < index <= offset + layout.lengths[chain]:
            return chain, index - offset
    start, length = layout.linker
    raise ValueError(f"scFv position {index} is in the linker (positions {start + 1} to {start + length}), which is not substituted")


def apply_mutations(layout: ScfvLayout, mutations: Iterable[Mutation]) -> str:
    """The scFv parent with every `(chain, position, wt, mut)` applied, positions in the crystal's chain numbering.

    Each wild-type residue is checked against the parent at the mapped position. Two mutations at one position, a
    mutation to the residue that is already there, and a position outside the variable domains are errors.
    """
    residues = list(layout.scfv)
    seen: set[int] = set()
    for chain, position, wt, mut in mutations:
        index = scfv_position(layout, chain, position)
        if index in seen:
            raise ValueError(f"two mutations at {chain}:{position}")
        seen.add(index)
        found = layout.scfv[index - 1]
        if found != wt:
            raise ValueError(
                f"{chain}:{wt}{position} does not match the parent, which has {found} at scFv position {index}. Either "
                f"the numbering differs from the crystal's or this is the wrong chain."
            )
        if len(mut) != 1 or mut not in AMINO_ACIDS or mut == wt:
            raise ValueError(f"{chain}:{wt}{position}{mut} is not a substitution to another amino acid")
        residues[index - 1] = mut
    return "".join(residues)


def apply_id(layout: ScfvLayout, variant_id: str) -> str:
    """`apply_mutations` for a variant id such as `B:S97A+A:Q27E`."""
    return apply_mutations(layout, parse_id(variant_id))


def substitutions(parent: str, variant: str) -> list[tuple[int, str, str]]:
    """(1-based position, parent residue, variant residue) wherever two scFv sequences differ. Substitutions only.

    The count of the list is the number of substitutions a design makes. Sequences of different length raise: an
    insertion or a deletion is not a substitution.
    """
    if len(parent) != len(variant):
        raise ValueError(f"the variant has {len(variant)} residues and the parent {len(parent)}; only substitutions are allowed")
    return [(i + 1, a, b) for i, (a, b) in enumerate(zip(parent, variant)) if a != b]


def crystal_mutations(layout: ScfvLayout, variant: str) -> list[Mutation]:
    """The differences between an scFv and the parent as `(chain, position, wt, mut)` in the crystal's numbering.

    The way back: for scoring a submitted or announced variant with the crystal-based scorers. A difference in the linker
    raises.
    """
    return [(*chain_position(layout, index), wt, mut) for index, wt, mut in substitutions(layout.scfv, variant)]


def verify_parent(layout: ScfvLayout, structure: Path | str) -> None:
    """Raise unless each chain of `structure` starts with the scFv's domain for that chain.

    This is what makes the offset trustworthy: the scFv's domains are the crystal's variable domains, residue for residue.
    Where they are not, for example a reformatted announced scFv with a changed framework residue, the crystal numbers
    still map but the crystal is no longer the parent, and this says where.
    """
    for chain in layout.offsets:
        crystal = chain_sequence(Path(structure), chain)
        domain = layout.domain(chain)
        window = crystal[: len(domain)]
        if window == domain:
            continue
        first = next((i for i, (a, b) in enumerate(zip(window, domain)) if a != b), len(window))
        here = crystal[first] if first < len(crystal) else "(end of chain)"
        raise ValueError(
            f"chain {chain} of {structure} does not start with the scFv's domain: first difference at position {first + 1} "
            f"(crystal {here}, scFv {domain[first]}). The crystal chain has {len(crystal)} residues, the domain {len(domain)}."
        )
