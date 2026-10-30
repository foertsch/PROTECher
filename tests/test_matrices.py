"""The shared matrix reader: known matrix in, known scores out.

No models and no GPU. The matrix is written by hand so the expected scores can
be worked out on paper, which is the only way to tell a correct implementation
from one that is off by a row.
"""

from __future__ import annotations

import numpy as np
import pytest

from design_pipeline.matrices import (
    IndexMismatchError,
    LogProbMatrix,
    from_position_mapping,
    score_candidates,
)
from design_pipeline.structures import AMINO_ACIDS

SEQUENCE = "ACDEF"
"""Five positions, one per row, each a different residue - so a matrix shifted
by one row cannot accidentally still match."""


def matrix(sequence: str = SEQUENCE, offset: int = 0) -> LogProbMatrix:
    """A matrix whose every cell is worked out from its own indices.

    `values[i, j] = -(i + 1) * 100 - j`, so a cell names the row and column it
    came from: any score read off it says, in its own value, which position and
    which residue were actually looked up. A shifted row shows up as a score
    off by exactly 100.
    """
    values = np.array(
        [
            [-(row + 1) * 100.0 - column for column in range(len(AMINO_ACIDS))]
            for row in range(len(sequence))
        ]
    )
    return LogProbMatrix(
        chain="H", values=values, sequence=sequence, offset=offset
    )


def test_a_score_is_the_mutant_minus_the_wild_type() -> None:
    """The one subtraction, on numbers chosen so the answer is checkable.

    Position 3 is row 2, so every cell there is -300 - column. D is column 2 of
    the alphabet and W is column 18, so the score is
    (-300 - 18) - (-300 - 2) = -16.
    """
    assert AMINO_ACIDS.index("D") == 2
    assert AMINO_ACIDS.index("W") == 18

    assert matrix().score(3, "D", "W") == pytest.approx(-16.0)
    # A substitution the matrix likes more than the wild type scores positive,
    # one it likes less scores negative. A flipped sign would swap these two
    # and nothing about the magnitudes would look wrong.
    assert matrix().score(3, "D", "A") == pytest.approx(2.0)
    # The reverse direction is not askable: position 3 holds D, so wt="W"
    # there is the misalignment this module refuses to score through.
    with pytest.raises(IndexMismatchError):
        matrix().score(3, "W", "D")


def test_the_wild_type_the_matrix_implies_matches_the_sequence_everywhere() -> None:
    """The real check. An index offset gives wrong scores rather than an error,
    so the alignment is asserted position by position rather than inferred."""
    subject = matrix()

    assert list(subject.positions) == [1, 2, 3, 4, 5]
    for position, expected in zip(subject.positions, SEQUENCE):
        assert subject.wild_type_at(position) == expected

    subject.assert_matches(SEQUENCE)


def test_a_matrix_shifted_by_one_row_is_an_error_not_a_score() -> None:
    """The failure this module exists for: a plausible number for the wrong
    position. The variant table says position 3 is D; a matrix built from the
    sequence one residue along has C there."""
    shifted = LogProbMatrix(
        chain="H", values=matrix().values[:4], sequence=SEQUENCE[1:]
    )

    with pytest.raises(IndexMismatchError, match="misaligned"):
        shifted.score(3, "D", "W")

    with pytest.raises(IndexMismatchError, match="First difference at position 1"):
        shifted.assert_matches(SEQUENCE)


def test_an_offset_matrix_scores_the_positions_it_actually_covers() -> None:
    """A tool that scores only part of the chain still lines up: the offset is
    carried, not guessed, so position 3 is this matrix's first row."""
    window = matrix(sequence=SEQUENCE[2:], offset=2)

    assert list(window.positions) == [3, 4, 5]
    assert window.wild_type_at(3) == "D"
    assert window.wild_type_at(2) is None  # outside: no answer, not an error
    window.assert_matches(SEQUENCE)

    # Row 0 of this matrix is -100 - column, and it is position 3.
    assert window.score(3, "D", "W") == pytest.approx((-100 - 18) - (-100 - 2))
    assert window.score(2, "C", "W") is None


def test_the_matrix_has_to_describe_the_sequence_it_was_given() -> None:
    with pytest.raises(ValueError, match="matrix rows for"):
        LogProbMatrix(chain="H", values=matrix().values[:3], sequence=SEQUENCE)

    with pytest.raises(ValueError, match=r"expected an \(L, 20\) matrix"):
        LogProbMatrix(chain="H", values=np.zeros((5, 4)), sequence=SEQUENCE)


def test_a_position_the_model_did_not_answer_is_empty_not_zero() -> None:
    """A zero log-probability is a probability of 1.0 - the most confident
    claim the scale can make - so an unanswered position has to stay empty."""
    sparse = from_position_mapping(
        "H", {2: {"C": -1.0, "W": -4.0}}, SEQUENCE
    )

    assert sparse.score(2, "C", "W") == pytest.approx(-3.0)
    assert sparse.score(1, "A", "W") is None  # never asked about
    assert np.isnan(sparse.values[0]).all()


def test_a_mapping_outside_the_chain_is_caught_at_build_time() -> None:
    with pytest.raises(IndexMismatchError, match="outside a chain"):
        from_position_mapping("H", {99: {"A": -1.0}}, SEQUENCE)


def test_score_candidates_answers_one_row_per_candidate() -> None:
    candidates = {
        "H:D3W": (("H", 3, "D", "W"),),
        "H:A1C": (("H", 1, "A", "C"),),
        "H:D3W+H:A1C": (("H", 3, "D", "W"), ("H", 1, "A", "C")),
        "L:D3W": (("L", 3, "D", "W"),),
    }

    rows = dict((variant_id, score) for variant_id, score, _raw in
                score_candidates({"H": matrix()}, candidates))

    assert rows["H:D3W"] == pytest.approx(-16.0)
    assert rows["H:A1C"] == pytest.approx((-100 - 1) - (-100 - 0))
    # A pair gets no score: adding two marginals is the additivity shortcut
    # this package refuses.
    assert rows["H:D3W+H:A1C"] is None
    # A chain with no matrix is no answer rather than an error.
    assert rows["L:D3W"] is None
    assert len(rows) == len(candidates)


def test_raw_and_score_agree_because_the_ratio_is_already_oriented() -> None:
    """A log-likelihood ratio against the wild type is higher-is-better as it
    stands, so there is no sign to flip here - and nothing should flip it."""
    for _variant_id, score, raw in score_candidates(
        {"H": matrix()}, {"H:D3W": (("H", 3, "D", "W"),)}
    ):
        assert score == raw
