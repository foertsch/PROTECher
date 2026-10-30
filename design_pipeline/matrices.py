"""One (L x 20) log-probability matrix, and the one way we score off it.

Several scorers in this package are the same thing wearing different clothes.
AntiFold, AbLang2, AntiBERTy, ESM-1v, ESM-2 and ESM-IF all answer the same
question - how likely is each of the 20 residues at each position of this chain
- and the score is always the same subtraction:

    score = M[position, mut] - M[position, wt]

A log-likelihood ratio against the wild type, which is already
higher-is-better, so `score` and `raw` agree for every model here. One
implementation, because six implementations of one subtraction are six chances
to get the indexing wrong.

The indexing is the dangerous part, and it is the reason this is a module
rather than a function in each adapter. A row is a position, and every tool
numbers positions its own way: AntiFold reports the PDB's own numbering, which
on an IMGT-renumbered antibody is neither 1-based nor contiguous; a masked
language model numbers the string it was handed, which for a paired model is
the heavy and light chains joined. A row that is off by one still yields a
plausible number - a real log-odds for the wrong position - and nothing
downstream can tell. So the matrix carries the wild-type sequence its rows are
aligned to, and every lookup is checked against it.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Mapping, Sequence

import numpy as np

from .structures import AMINO_ACIDS

Mutation = tuple[str, int, str, str]
"""(chain, position, wt, mut), as in design_pipeline.variants."""


class IndexMismatchError(ValueError):
    """A matrix row does not hold the residue the variant table says is there.

    Raised rather than scored. A shifted matrix produces numbers that look
    exactly like the real ones, so this is the one failure in the package that
    has to be loud.
    """


@dataclass(frozen=True)
class LogProbMatrix:
    """Per-position log probabilities for one chain, and what sits at each one.

    `values` is (L, 20): row i is the chain's position `offset + i + 1`, column
    j the residue `alphabet[j]`. `sequence` is the wild-type sequence those
    rows describe, in the same order - the thing that makes a misalignment
    detectable, since it is checked against the variant table's own `wt`.

    `offset` is how far into the chain the first row sits, so a tool that
    scores only the variable domain still lines up with positions counted from
    the start of the chain.
    """

    chain: str
    values: np.ndarray
    sequence: str
    offset: int = 0
    alphabet: str = AMINO_ACIDS

    def __post_init__(self) -> None:
        values = np.asarray(self.values, dtype=float)
        object.__setattr__(self, "values", values)
        if values.ndim != 2 or values.shape[1] != len(self.alphabet):
            raise ValueError(
                f"{self.chain}: expected an (L, {len(self.alphabet)}) matrix, "
                f"got {values.shape}"
            )
        if values.shape[0] != len(self.sequence):
            raise ValueError(
                f"{self.chain}: {values.shape[0]} matrix rows for "
                f"{len(self.sequence)} residues of sequence; the two have to "
                f"describe the same positions"
            )
        if self.offset < 0:
            raise ValueError(f"{self.chain}: offset {self.offset} is negative")

    @property
    def positions(self) -> range:
        """The chain positions this matrix covers, 1-based and inclusive."""
        return range(self.offset + 1, self.offset + len(self.sequence) + 1)

    def wild_type_at(self, position: int) -> str | None:
        """The residue this matrix says is at `position`, or None if outside."""
        index = position - 1 - self.offset
        return self.sequence[index] if 0 <= index < len(self.sequence) else None

    def log_prob(self, position: int, residue: str) -> float | None:
        index = position - 1 - self.offset
        column = self.alphabet.find(residue)
        if not (0 <= index < self.values.shape[0]) or column < 0:
            return None
        value = float(self.values[index, column])
        return None if not np.isfinite(value) else value

    def score(self, position: int, wt: str, mut: str) -> float | None:
        """M[position, mut] - M[position, wt], or None where it is undefined.

        `wt` is checked against the sequence first. A position this matrix does
        not cover is None - a legitimate "no answer" - but a position it does
        cover holding a different residue is a misalignment and raises.
        """
        here = self.wild_type_at(position)
        if here is None:
            return None
        if here != wt:
            raise IndexMismatchError(
                f"chain {self.chain} position {position}: the variant table "
                f"says the wild type is {wt!r}, this matrix has {here!r}. The "
                f"rows are misaligned with the chain, so every score off it "
                f"would be a real number for the wrong position."
            )
        mutant = self.log_prob(position, mut)
        wild = self.log_prob(position, wt)
        if mutant is None or wild is None:
            return None
        return mutant - wild

    def assert_matches(self, sequence: str) -> None:
        """Every row holds the residue `sequence` has at that position.

        The round-trip check as an assertion rather than a test helper, so an
        adapter can call it the moment it has built a matrix, before any score
        comes off it.
        """
        window = sequence[self.offset : self.offset + len(self.sequence)]
        if window == self.sequence:
            return
        first = next(
            (i for i, (a, b) in enumerate(zip(window, self.sequence)) if a != b),
            min(len(window), len(self.sequence)),
        )
        raise IndexMismatchError(
            f"chain {self.chain}: matrix rows do not match the chain sequence. "
            f"First difference at position {self.offset + first + 1}: the chain "
            f"has {window[first : first + 1]!r}, the matrix has "
            f"{self.sequence[first : first + 1]!r}. The rows cover "
            f"{len(self.sequence)} residues from offset {self.offset}; the "
            f"chain has {len(sequence)}."
        )


def from_position_mapping(
    chain: str,
    log_probs: Mapping[int, Mapping[str, float]],
    sequence: str,
    *,
    alphabet: str = AMINO_ACIDS,
) -> LogProbMatrix:
    """A matrix from {position: {residue: log probability}}.

    What a masked language model produces: one distribution per position it was
    asked about, nothing for the rest. The gaps stay NaN, which reads back as
    "no score" rather than as zero - a zero here is a log-probability of 1.0,
    the most confident claim the scale can make.

    Positions are 1-based into `sequence`, the chain as the variant table
    counts it.
    """
    values = np.full((len(sequence), len(alphabet)), np.nan)
    for position, row in log_probs.items():
        index = position - 1
        if not (0 <= index < len(sequence)):
            raise IndexMismatchError(
                f"chain {chain}: position {position} is outside a chain of "
                f"{len(sequence)} residues"
            )
        for residue, value in row.items():
            column = alphabet.find(residue)
            if column >= 0:
                values[index, column] = float(value)
    return LogProbMatrix(
        chain=chain, values=values, sequence=sequence, alphabet=alphabet
    )


def score_candidates(
    matrices: Mapping[str, LogProbMatrix],
    candidates: Mapping[str, Sequence[Mutation]],
) -> list[tuple[str, float | None, float | None]]:
    """The rows `write_scores` wants, for every candidate, from the matrices.

    Shared by every matrix scorer, which is what makes two rules hold
    everywhere at once rather than once per adapter:

    * A candidate with several mutations gets no score. Adding two marginals is
      the additivity shortcut this package refuses; scoring the real mutant
      means a pass per candidate, which is what a per-variant scorer is for.
    * A position no matrix covers gets no score, not a zero.

    The score is already higher-is-better, so `raw` is the same number.
    """
    rows: list[tuple[str, float | None, float | None]] = []
    for variant_id, mutations in candidates.items():
        if len(mutations) != 1:
            rows.append((variant_id, None, None))
            continue
        chain, position, wt, mut = mutations[0]
        matrix = matrices.get(chain)
        if matrix is None:
            rows.append((variant_id, None, None))
            continue
        score = matrix.score(position, wt, mut)
        rows.append((variant_id, score, score))
    return rows
