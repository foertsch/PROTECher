"""Step 4: the gate. Failures are dropped, never traded off against score."""

from __future__ import annotations

from pathlib import Path

import pytest

from design_pipeline import enumerate_variants, gate, get_structure, score_variants
from design_pipeline.adapters._gates import new_motifs


@pytest.fixture
def complex_(tmp_path: Path):
    return get_structure(tmp_path / "runs", model_name="dummy_complex")


@pytest.mark.parametrize(
    "wild, mutant, expected",
    [
        # N-x-S/T, x not proline.
        ("AAAST", "ANAST", "n_glycosylation@2"),
        ("ANPSA", "ANPSA", None),  # x is proline: not a sequon
        ("AAAAA", "ANGAA", "deamidation@2"),
        ("AAAAA", "ADGAA", "isomerization@2"),
        ("AAAAA", "ADPAA", "fragmentation@2"),
    ],
)
def test_the_regexes_find_the_motifs_a_mutation_creates(
    wild: str, mutant: str, expected: str | None
) -> None:
    found = new_motifs(wild, mutant)

    if expected is None:
        assert found == []
    else:
        assert expected in found


def test_a_motif_that_was_already_there_is_not_this_variants_fault() -> None:
    assert new_motifs("ANGS", "ANGA") == []


def test_liabilities_drops_the_candidates_that_create_a_motif(complex_) -> None:
    """Pure sequence regex, milliseconds, and it catches the most common
    own-goal: a new glycosylation sequon bought with an affinity gain."""
    variants = enumerate_variants(complex_, chain="H")
    survivors = gate(variants, checks=["liabilities"])

    assert survivors.candidates < variants.candidates
    assert survivors.frame["pass_liabilities"].all()


def test_an_introduced_cysteine_is_a_failure(complex_) -> None:
    variants = enumerate_variants(complex_, chain="H", positions=[5], mutations="C")
    with pytest.warns(RuntimeWarning, match="dropped every"):
        survivors = gate(variants, checks=["liabilities"])

    assert survivors.candidates == 0
    annotated = Path(survivors.config["annotated"])
    import pandas as pd

    assert "free_cysteine" in pd.read_parquet(annotated)["reason_liabilities"].iloc[0]


def test_the_annotated_table_keeps_the_failures_and_their_reasons(complex_) -> None:
    """So a gate that dropped too much can be re-tuned without re-running the
    tools."""
    import pandas as pd

    variants = enumerate_variants(complex_, chain="H")
    survivors = gate(variants, checks=["liabilities"])
    annotated = pd.read_parquet(survivors.config["annotated"])

    assert len(annotated) == variants.rows
    failures = annotated[~annotated["pass_liabilities"]]
    assert len(failures) > 0
    assert (failures["reason_liabilities"] != "").all()


def test_several_checks_all_have_to_pass(complex_) -> None:
    variants = enumerate_variants(complex_, chain="H", mutations="AW")
    one = gate(variants, checks=["liabilities"])
    both = gate(variants, checks=["liabilities", "dummy_check"])

    assert both.candidates < one.candidates
    assert {"pass_liabilities", "pass_dummy_check"} <= set(both.columns)


def test_scores_survive_the_gate_untouched(complex_) -> None:
    """The gate filters rows; it does not touch the columns the caller ranks
    on, and it never weighs a failure against a score."""
    variants = enumerate_variants(complex_, chain="H", mutations="AW")
    scores = score_variants(complex_, variants, models=["dummy_scorer"])
    survivors = gate(scores, checks=["liabilities"])

    before = dict(zip(scores.frame["variant_id"], scores.frame["score_dummy_scorer"]))
    after = dict(
        zip(survivors.frame["variant_id"], survivors.frame["score_dummy_scorer"])
    )
    assert after and all(before[k] == v for k, v in after.items())


def test_the_gate_finds_the_structure_on_the_table(complex_) -> None:
    """`gate(scores, checks=[...])` is the documented call, so the table has to
    carry the complex the checks read."""
    variants = enumerate_variants(complex_, chain="H", positions=[1], mutations="A")

    assert variants.structure == complex_.path
    assert gate(variants, checks=["liabilities"]).structure == complex_.path


def test_a_table_with_no_structure_asks_for_one(tmp_path, complex_) -> None:
    from dataclasses import replace

    variants = enumerate_variants(complex_, chain="H", positions=[1], mutations="A")
    detached = replace(variants, structure=None)

    with pytest.raises(ValueError, match="needs the complex"):
        gate(detached, checks=["liabilities"])
    assert gate(detached, checks=["liabilities"], complex_=complex_).candidates == 1


def test_gating_everything_away_warns_rather_than_failing_silently(complex_) -> None:
    variants = enumerate_variants(complex_, chain="H", positions=[5], mutations="C")

    with pytest.warns(RuntimeWarning, match="dropped every"):
        survivors = gate(variants, checks=["liabilities"])
    assert survivors.rows == 0
