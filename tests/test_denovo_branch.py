"""The optional de novo branch: `diffuse_loops` and the pieces under it.

Kept, demoted. What is tested here is that it produces an ordinary variant
table and refuses to produce a misleading one; the models themselves need a
container runtime, so the fixture adapters stand in and the length check they
trip is the interesting half.
"""

from __future__ import annotations

from pathlib import Path

import pytest

from design_pipeline import get_structure, numbering
from design_pipeline.denovo import _cdr_span, _contigs, _diff, diffuse_loops


@pytest.fixture
def complex_(tmp_path: Path):
    return get_structure(tmp_path / "runs", model_name="dummy_complex")


def test_the_branch_is_not_on_the_default_path() -> None:
    """`diffuse_loops` is the only way in, and the five steps do not touch
    rfdiffusion or mpnn."""
    import design_pipeline as dp

    assert "diffuse_loops" in dp.__all__
    assert "make_protein_backbones" not in dp.__all__
    assert "design_sequences" not in dp.__all__


def test_a_diff_is_every_position_that_changed() -> None:
    mutations = _diff("ACDEF", "AWDEY", "H", Path("x.fa"))

    assert mutations == [("H", 2, "C", "W"), ("H", 5, "F", "Y")]


def test_a_length_change_is_refused_rather_than_aligned() -> None:
    """A variant table is a position-by-position diff; an indel leaves nothing
    to diff, and guessing an alignment would put mutations at wrong
    positions."""
    with pytest.raises(ValueError, match="original length"):
        _diff("ACDEF", "ACDE", "H", Path("x.fa"))


def test_the_contig_string_pins_the_loop_to_its_own_length(complex_) -> None:
    contigs = _contigs(complex_.path, "H", 5, 12)

    assert contigs.startswith("[") and contigs.endswith("]")
    assert "8-8" in contigs  # the rebuilt span, fixed length
    assert "H1-4" in contigs and "H13-" in contigs  # the flanks stay
    assert "/0 " in contigs  # the other chains are context


def test_a_cdr_name_needs_numbering_and_says_so(complex_) -> None:
    if numbering.available():
        pytest.skip("anarci is installed, so H3 resolves")
    with pytest.raises(ValueError, match="span="):
        _cdr_span(complex_.path, "H3")


def test_a_bad_cdr_name_is_rejected(complex_) -> None:
    with pytest.raises(ValueError, match="not a CDR"):
        _cdr_span(complex_.path, "H9")


def test_the_fixture_models_trip_the_length_check(complex_) -> None:
    """The fixture backbone and sequence adapters write a fixed-size design,
    so this is the error path rather than a table. The real pair rebuilds the
    loop at its original length; see the module docstring on how far outside
    the scorers' validated range those variants sit."""
    with pytest.raises(ValueError, match="length"):
        diffuse_loops(
            complex_,
            span=(5, 12),
            num_designs=1,
            backbone_model="dummy_backbones",
            sequence_model="dummy_sequences",
        )
