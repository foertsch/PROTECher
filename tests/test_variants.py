"""The variant table: ids, the scan, and pairing."""

from __future__ import annotations

from pathlib import Path

import pytest

from design_pipeline import enumerate_variants, get_structure, select_pairs
from design_pipeline.adapters._pdb import FIXTURE_CHAINS
from design_pipeline.structures import AMINO_ACIDS
from design_pipeline.variants import SCHEMA, make_id, mutant_sequence, parse_id


@pytest.fixture
def complex_(tmp_path: Path):
    return get_structure(tmp_path / "runs", model_name="dummy_complex")


def test_an_id_is_canonical_whichever_order_the_pair_arrives_in() -> None:
    """One candidate, one cache entry: a pair discovered either way round is
    the same id."""
    forward = make_id([("H", 30, "T", "Y"), ("L", 91, "S", "A")])
    backward = make_id([("L", 91, "S", "A"), ("H", 30, "T", "Y")])

    assert forward == backward == "H:T30Y+L:S91A"
    assert parse_id(forward) == (("H", 30, "T", "Y"), ("L", 91, "S", "A"))


@pytest.mark.parametrize("bad", ["", "H30Y", "H:30Y", "H:T30", "H:T30Y+"])
def test_a_malformed_id_is_rejected(bad: str) -> None:
    with pytest.raises(ValueError):
        parse_id(bad)


def test_the_default_is_the_full_scan(complex_) -> None:
    """Every position of the chain against the other 19 amino acids."""
    variants = enumerate_variants(complex_, chain="H")
    length = len(FIXTURE_CHAINS["H"])

    assert variants.candidates == length * (len(AMINO_ACIDS) - 1)
    assert variants.rows == variants.candidates  # singles: one row each
    assert list(variants.frame.columns) == list(SCHEMA)
    assert variants.path.suffix == ".parquet"


def test_the_scan_never_includes_the_wild_type_residue(complex_) -> None:
    frame = enumerate_variants(complex_, chain="H").frame

    assert not (frame["wt"] == frame["mut"]).any()


def test_positions_and_mutations_narrow_the_scan(complex_) -> None:
    variants = enumerate_variants(
        complex_, chain="H", positions=[3, 10], mutations="AW"
    )
    frame = variants.frame

    assert variants.candidates == 4
    assert sorted(set(frame["position"])) == [3, 10]
    assert sorted(set(frame["mut"])) == ["A", "W"]


def test_a_position_outside_the_chain_is_an_error(complex_) -> None:
    with pytest.raises(IndexError, match="outside chain"):
        enumerate_variants(complex_, chain="H", positions=[9999])


def test_only_restricts_the_table_and_accepts_a_pair(complex_) -> None:
    """How the harness scores a measured set, double mutants included."""
    wanted = ["H:E1A", f"H:{FIXTURE_CHAINS['H'][2]}3W+L:{FIXTURE_CHAINS['L'][1]}2G"]
    variants = enumerate_variants(complex_, chain="H", only=wanted)

    assert variants.candidates == 2
    assert variants.rows == 3  # the pair contributes two rows
    assert set(variants.frame["chain"]) == {"H", "L"}


def test_only_checks_the_wild_type_against_the_structure(complex_) -> None:
    """The one error a ground-truth set cannot be trusted to avoid: numbering
    that is IMGT or Kabat rather than chain positions."""
    with pytest.raises(ValueError, match="does not match"):
        enumerate_variants(complex_, chain="H", only=["H:W1A"])


def test_imgt_is_empty_without_anarci(complex_) -> None:
    """Degrading to raw sequence indexing, with a warning, is the documented
    behaviour; joins are on chain positions regardless."""
    from design_pipeline import numbering

    frame = enumerate_variants(complex_, chain="H", positions=[1]).frame
    if numbering.available():
        pytest.skip("anarci is installed, so the labels are real")
    assert (frame["imgt"] == "").all()


def test_mutant_sequence_applies_every_mutation_of_the_candidate(complex_) -> None:
    """A pair is scored as one mutant, which is what makes summing two single
    scores unnecessary as well as forbidden."""
    wild, mutant = mutant_sequence(
        complex_.path, [("H", 1, "E", "A"), ("H", 3, "Q", "W")], "H"
    )

    assert wild == FIXTURE_CHAINS["H"]
    assert mutant[0] == "A" and mutant[2] == "W"
    assert sum(a != b for a, b in zip(wild, mutant)) == 2


def test_pairs_are_far_enough_apart_and_capped(complex_) -> None:
    singles = enumerate_variants(complex_, chain="H", mutations="A")
    pairs = select_pairs(singles, max_pairs=6, min_distance=10.0)
    frame = pairs.frame

    assert pairs.candidates == 6
    assert pairs.rows == 12
    assert (frame["distance"] >= 10.0).all()
    assert pairs.config["distance_units"] == "angstrom"
    # Every candidate is two rows under one id, in the same schema as a single.
    assert set(frame.groupby("variant_id").size()) == {2}


def test_pairing_spreads_across_the_top_rather_than_fixing_the_best_single(
    complex_,
) -> None:
    """Triangular order, so 20 pairs are not 20 variations on candidate one."""
    singles = enumerate_variants(complex_, chain="H", mutations="A")
    pairs = select_pairs(singles, max_pairs=6, min_distance=1.0)
    first_positions = {
        tuple(sorted(group["position"]))[0]
        for _, group in pairs.frame.groupby("variant_id")
    }

    assert len(first_positions) > 1


def test_score_columns_do_not_survive_pairing(tmp_path, complex_) -> None:
    """The additivity guard: a pair must not inherit a number computed for one
    of its halves under a `score_` name."""
    from design_pipeline import score_variants

    singles = enumerate_variants(complex_, chain="H", mutations="A")
    scored = score_variants(complex_, singles, models=["dummy_scorer"])
    pairs = select_pairs(scored, max_pairs=3, min_distance=10.0)

    assert "score_dummy_scorer" in scored.columns
    assert pairs.score_columns == ()
    assert not any(c.startswith("score_") for c in pairs.columns)


def test_pairing_a_table_that_already_holds_pairs_is_refused(complex_) -> None:
    singles = enumerate_variants(complex_, chain="H", mutations="A")
    pairs = select_pairs(singles, max_pairs=3, min_distance=10.0)

    with pytest.raises(ValueError, match="single mutations"):
        select_pairs(pairs, max_pairs=3)


def test_an_impossible_distance_is_an_error_rather_than_an_empty_table(
    complex_,
) -> None:
    singles = enumerate_variants(complex_, chain="H", positions=[1, 2], mutations="A")

    with pytest.raises(ValueError, match="min_distance"):
        select_pairs(singles, max_pairs=5, min_distance=10_000.0)


def test_by_sorts_on_one_column_and_only_one(tmp_path, complex_) -> None:
    """Ranking stays the caller's business: `by` picks a column, it never
    combines them."""
    from design_pipeline import score_variants

    singles = enumerate_variants(complex_, chain="H", mutations="A")
    scored = score_variants(complex_, singles, models=["dummy_scorer"])

    pairs = select_pairs(scored, max_pairs=2, min_distance=1.0, by="score_dummy_scorer")
    best = scored.frame.sort_values("score_dummy_scorer", ascending=False)
    top_position = int(best.iloc[0]["position"])
    paired_positions = set(pairs.frame["position"])

    assert top_position in paired_positions
    with pytest.raises(KeyError, match="no column"):
        select_pairs(scored, by="score_nonexistent")
