"""The validation harness, on the shipped fixture systems: offline, no GPU."""

from __future__ import annotations

from pathlib import Path

import numpy as np
import pytest

from design_pipeline.harness import (
    FIXTURE_SYSTEMS,
    SYSTEMS,
    System,
    magnitude_plot,
    report,
    run_systems,
)
from design_pipeline.metrics import precision_at_k, rankdata, spearman


def test_the_fixtures_run_end_to_end_with_no_network(tmp_path: Path) -> None:
    results = run_systems(models=["dummy_scorer"], output_dir=tmp_path / "harness")

    assert results.systems == FIXTURE_SYSTEMS
    assert results.rows > 0
    frame = results.frame
    assert set(frame["system"]) == set(FIXTURE_SYSTEMS)
    assert {"measured_affinity", "oriented", "score_dummy_scorer"} <= set(frame.columns)
    # One row per candidate, so a double mutant is not counted twice.
    assert not frame.duplicated(["system", "variant_id"]).any()


def test_the_report_breaks_results_out_per_scorer(tmp_path: Path) -> None:
    """Per scorer, never for the ensemble only: choosing which column to trust
    is the entire point, and an ensemble-level number hides it."""
    results = run_systems(
        models=["dummy_scorer", "dummy_per_variant"],
        max_variants=40,
        output_dir=tmp_path / "harness",
    )
    text = report(results, k=5)

    for scorer in ("dummy_scorer", "dummy_per_variant"):
        for system in FIXTURE_SYSTEMS:
            assert any(
                scorer in line and system in line for line in text.splitlines()
            ), f"{scorer} / {system} missing from the report"
    assert "spearman" in text and "p@5" in text
    assert "(mean)" in text  # per scorer, across systems


def test_the_report_also_asks_whether_magnitude_predicts_effect_size(
    tmp_path: Path,
) -> None:
    """The second table, and the question only the harness can answer: if a
    mutation scoring twice as high does not bind twice as tightly, only the
    ordering is usable and candidates have to be picked by rank."""
    results = run_systems(models=["dummy_scorer"], output_dir=tmp_path / "harness")
    text = report(results)

    assert "score magnitude against measured effect size" in text
    assert "r2" in text and "slope" in text
    # Both tables, in one report: rank order first, then magnitude.
    assert text.index("spearman") < text.index("score magnitude")


def test_the_magnitude_plot_is_written_and_needs_no_display(tmp_path: Path) -> None:
    pytest.importorskip("matplotlib", reason="the notebooks extra")
    results = run_systems(models=["dummy_scorer"], output_dir=tmp_path / "harness")
    out = tmp_path / "plots" / "magnitude.png"

    figure = magnitude_plot(results, out)

    assert out.is_file() and out.stat().st_size > 0
    # One panel per score column.
    assert len(figure.axes) >= 1


def test_the_harness_statistics_are_the_ones_in_metrics(tmp_path: Path) -> None:
    """The harness contains none of its own: a Spearman here and a Spearman in
    a notebook plot have to be the same Spearman. Checked by reading the
    source, because the failure mode is somebody inlining a mean."""
    import design_pipeline.harness as harness

    source = Path(harness.__file__).read_text()

    assert "def _nanmean" not in source
    for inlined in ("np.corrcoef", "np.polyfit", "scipy"):
        assert inlined not in source, f"{inlined} belongs in metrics, not here"


def test_a_scorer_with_no_skill_reports_a_correlation_near_zero(tmp_path) -> None:
    """What the fixtures can actually check. The numbers are invented, so a
    hash-keyed scorer should land near zero rather than crash."""
    results = run_systems(models=["dummy_scorer"], output_dir=tmp_path / "harness")
    frame = results.frame
    rho = spearman(frame["score_dummy_scorer"], frame["oriented"])

    assert abs(rho) < 0.5


def test_the_harness_is_the_ordinary_pipeline(tmp_path: Path) -> None:
    """Not a second code path: the run directories it leaves behind are the
    ordinary content-addressed ones, and a rerun resumes from them."""
    first = run_systems(models=["dummy_scorer"], output_dir=tmp_path / "harness")
    second = run_systems(models=["dummy_scorer"], output_dir=tmp_path / "harness")

    assert first.rows == second.rows
    assert (tmp_path / "harness" / "fixture_a" / "scores").is_dir()
    assert (tmp_path / "harness" / "fixture_a" / "variants").is_dir()


def test_an_unknown_system_names_the_ones_that_exist(tmp_path: Path) -> None:
    with pytest.raises(KeyError, match="fixture_a"):
        run_systems(["no_such_system"], output_dir=tmp_path / "harness")


def test_a_real_system_says_which_download_it_needs(tmp_path, monkeypatch) -> None:
    """The measurements are not redistributable, so the error has to name the
    source rather than just the missing file."""
    from design_pipeline.errors import MissingDependencyError
    from design_pipeline.harness.systems import GROUND_TRUTH_ENV

    monkeypatch.setenv(GROUND_TRUTH_ENV, str(tmp_path / "nothing"))
    with pytest.raises(MissingDependencyError, match="FLAb"):
        SYSTEMS["1n8z_her2"].truth_path()


def test_a_system_can_be_added_from_outside_the_package(tmp_path) -> None:
    SYSTEMS["mine"] = System(
        name="mine",
        chain="H",
        truth="fixture_a.csv",
        source="fixture_a.pdb",
        chains=("H",),
        higher_is_better=True,
        fixture=True,
    )
    try:
        results = run_systems(
            ["mine"], models=["dummy_scorer"], output_dir=tmp_path / "harness"
        )
    finally:
        del SYSTEMS["mine"]

    assert results.rows > 0


def test_ranks_average_ties() -> None:
    assert list(rankdata([10.0, 20.0, 20.0, 30.0])) == [1.0, 2.5, 2.5, 4.0]


def test_spearman_is_one_on_a_monotone_relationship() -> None:
    assert spearman([1.0, 2.0, 3.0, 4.0], [10.0, 20.0, 30.0, 40.0]) == pytest.approx(1)
    assert spearman([1.0, 2.0, 3.0, 4.0], [40.0, 30.0, 20.0, 10.0]) == pytest.approx(-1)


def test_an_uncomputable_correlation_is_nan_rather_than_zero() -> None:
    """Zero would read as 'no skill', which is a claim; NaN is the absence of
    one."""
    assert np.isnan(spearman([1.0, 2.0], [1.0, 2.0]))  # too few points
    assert np.isnan(spearman([1.0, 1.0, 1.0], [1.0, 2.0, 3.0]))  # constant column
    assert np.isnan(spearman([1.0, 2.0, np.nan], [1.0, np.nan, 3.0]))


def test_precision_at_k_is_the_overlap_of_the_two_top_tens() -> None:
    predicted = [5.0, 4.0, 3.0, 2.0, 1.0]
    assert precision_at_k(predicted, predicted, k=2) == 1.0
    assert precision_at_k(predicted, [1.0, 2.0, 3.0, 4.0, 5.0], k=2) == 0.0
    assert np.isnan(precision_at_k([1.0, 2.0], [1.0, 2.0], k=10))
