"""The Category 1 metrics - the ones that need a known answer.

Against hand-computed fixtures. The Category 2 diagnostics are in
test_diagnostics.py, which is the same separation the module makes.
"""

from __future__ import annotations

import json
import math
from pathlib import Path

import numpy as np
import pytest

from design_pipeline import metrics
from design_pipeline.adapters._pdb import (
    ANTIBODY_SPACING,
    DOCKED,
    FIXTURE_CHAINS,
    MISPLACED,
    write_ca_pdb,
)

DOCKED_ORIGINS = {
    "H": (0.0, 0.0, 0.0),
    "L": (ANTIBODY_SPACING, 0.0, 0.0),
    "A": DOCKED,
}
WRONG_ORIGINS = {**DOCKED_ORIGINS, "A": MISPLACED}


@pytest.fixture
def reference(tmp_path: Path) -> Path:
    return write_ca_pdb(tmp_path / "ref.pdb", FIXTURE_CHAINS, origins=DOCKED_ORIGINS)


@pytest.fixture
def wrong_pose(tmp_path: Path) -> Path:
    """The same complex with the antigen on the opposite face: a different
    epitope, not a jittered copy of the right answer."""
    return write_ca_pdb(tmp_path / "wrong.pdb", FIXTURE_CHAINS, origins=WRONG_ORIGINS)


# --------------------------------------------------------------------------
# Residue correspondence
# --------------------------------------------------------------------------


def test_identical_sequences_map_one_to_one() -> None:
    assert metrics.align_positions("ACDEF", "ACDEF") == {1: 1, 2: 2, 3: 3, 4: 4, 5: 5}


def test_a_deletion_shifts_the_mapping_and_skips_the_gap() -> None:
    """The real case: a prediction numbered 1..N against a deposited chain with
    a disordered loop missing from the model."""
    mapping = metrics.align_positions("ACDEFGH", "ACDGH")

    assert mapping[1] == 1 and mapping[2] == 2 and mapping[3] == 3
    assert mapping[6] == 4 and mapping[7] == 5
    # E and F are in the reference and not in the prediction, so they are
    # absent rather than mapped onto a neighbour.
    assert 4 not in mapping and 5 not in mapping


def test_a_leading_insertion_is_handled() -> None:
    mapping = metrics.align_positions("CDEF", "ACDEF")

    assert mapping == {1: 2, 2: 3, 3: 4, 4: 5}


def test_identity_is_over_the_shorter_sequence() -> None:
    assert metrics.sequence_identity("ACDEF", "ACDEF") == 1.0
    assert metrics.sequence_identity("ACDEF", "ACDEY") == pytest.approx(0.8)
    assert metrics.sequence_identity("ACDE", "ACDEFGHI") == 1.0
    assert metrics.sequence_identity("", "ACDE") == 0.0


# --------------------------------------------------------------------------
# Structure comparison
# --------------------------------------------------------------------------


def test_the_interface_is_the_residues_in_contact(reference) -> None:
    epitope = metrics.interface_residues(reference, ["A"], ["H", "L"])
    paratope = metrics.interface_residues(reference, ["H", "L"], ["A"])

    assert epitope and paratope
    assert {chain for chain, _ in epitope} == {"A"}
    assert {chain for chain, _ in paratope} <= {"H", "L"}
    # 1-based into the chain as given, the same convention as the variant table.
    assert all(1 <= position <= len(FIXTURE_CHAINS[chain]) for chain, position in epitope)


def test_a_wider_cutoff_finds_at_least_as_much(reference) -> None:
    near = metrics.interface_residues(reference, ["A"], ["H", "L"], cutoff=4.0)
    far = metrics.interface_residues(reference, ["A"], ["H", "L"], cutoff=8.0)

    assert near <= far


def test_the_default_split_reads_h_and_l_as_the_antibody(reference) -> None:
    assert metrics.interface_residues(reference) == metrics.interface_residues(
        reference, ["H", "L"], ["A"]
    )


def test_a_structure_with_no_partner_chain_is_an_error(tmp_path) -> None:
    alone = write_ca_pdb(tmp_path / "one.pdb", {"H": FIXTURE_CHAINS["H"]})

    with pytest.raises(ValueError, match="no interface"):
        metrics.interface_residues(alone)


def test_a_structure_compared_to_itself_scores_perfectly(reference) -> None:
    assert metrics.interface_rmsd(reference, reference) == pytest.approx(0.0, abs=1e-4)
    assert metrics.epitope_recall(reference, reference, ["A"], ["H", "L"]) == 1.0


def test_a_wrong_pose_loses_the_epitope(reference, wrong_pose) -> None:
    """The metric that matters: the antigen is still in contact, so a
    coordinate RMSD can look survivable while the paratope sits on the wrong
    patch - which is what silently poisons every scorer downstream."""
    recall = metrics.epitope_recall(wrong_pose, reference, ["A"], ["H", "L"])

    assert recall < 0.5
    assert metrics.interface_rmsd(wrong_pose, reference) > metrics.interface_rmsd(
        reference, reference
    )


def test_epitope_recall_is_nan_when_there_is_no_epitope(tmp_path, reference) -> None:
    apart = write_ca_pdb(
        tmp_path / "apart.pdb", FIXTURE_CHAINS, origins={**DOCKED_ORIGINS, "A": (3.0, 60.0, 0.0)}
    )
    assert math.isnan(metrics.epitope_recall(apart, apart, ["A"], ["H", "L"]))
    # Against a reference that does have one, nothing is recalled.
    assert metrics.epitope_recall(apart, reference, ["A"], ["H", "L"]) == 0.0


def test_interface_rmsd_is_nan_when_too_little_corresponds(tmp_path) -> None:
    """Under three corresponding residues the superposition is
    underdetermined, so the answer is NaN rather than a number."""
    pair = {"H": FIXTURE_CHAINS["H"], "A": FIXTURE_CHAINS["A"]}
    origins = {"H": (0.0, 0.0, 0.0), "A": DOCKED}
    reference = write_ca_pdb(tmp_path / "two.pdb", pair, origins=origins)
    # One residue per chain: two corresponding positions in total.
    stub = write_ca_pdb(tmp_path / "stub.pdb", {"H": "A", "A": "M"}, origins=origins)

    assert metrics.interface_residues(reference, ["A"], ["H"])
    assert math.isnan(metrics.interface_rmsd(stub, reference, ["H"], ["A"]))


# --------------------------------------------------------------------------
# DockQ
# --------------------------------------------------------------------------


def test_the_quality_bands_are_the_published_ones() -> None:
    assert metrics.quality(0.85) == "high"
    assert metrics.quality(metrics.HIGH) == "high"
    assert metrics.quality(0.60) == "medium"
    assert metrics.quality(metrics.MEDIUM) == "medium"
    assert metrics.quality(0.30) == "acceptable"
    assert metrics.quality(metrics.ACCEPTABLE) == "acceptable"
    assert metrics.quality(0.10) == "incorrect"
    assert metrics.quality(float("nan")) == "unknown"


def test_dockq_output_is_reduced_to_the_best_interface() -> None:
    """A Fab against an antigen has an H-antigen and an L-antigen interface;
    the question is whether the complex is right, not which chain holds more
    of the contact."""
    summary = metrics.summarize_dockq(
        {
            "HA": {"DockQ": 0.31, "fnat": 0.4, "iRMSD": 2.2, "LRMSD": 7.0},
            "LA": {"DockQ": 0.67, "fnat": 0.7, "iRMSD": 1.1, "LRMSD": 3.0},
        }
    )

    assert summary["dockq"] == 0.67
    assert summary["quality"] == "medium"
    assert summary["irmsd"] == 1.1
    assert summary["interface"] == "LA"
    assert summary["interfaces"] == 2


def test_dockq_component_names_that_moved_are_still_read() -> None:
    """They have changed between DockQ versions, which is the reason this
    parse is a separate, tested function rather than inline."""
    summary = metrics.summarize_dockq({"x": {"dockq": 0.5, "iRMS": 2.0, "LRMS": 6.0}})

    assert summary["dockq"] == 0.5
    assert summary["irmsd"] == 2.0
    assert summary["lrmsd"] == 6.0


def test_no_interfaces_is_not_a_crash() -> None:
    summary = metrics.summarize_dockq({})

    assert math.isnan(summary["dockq"])
    assert summary["quality"] == "unknown"


def test_dockq_says_which_extra_to_install(reference) -> None:
    from design_pipeline.errors import MissingDependencyError

    pytest.importorskip  # noqa: B018 - documents why this is conditional
    try:
        import DockQ  # noqa: F401
    except ModuleNotFoundError:
        with pytest.raises(MissingDependencyError, match="metrics"):
            metrics.dockq(reference, reference)
    else:  # pragma: no cover - only where the extra is installed
        assert "dockq" in metrics.dockq(reference, reference)


# --------------------------------------------------------------------------
# Ranking
# --------------------------------------------------------------------------


def test_ranks_average_ties() -> None:
    assert list(metrics.rankdata([10.0, 20.0, 20.0, 30.0])) == [1.0, 2.5, 2.5, 4.0]


def test_spearman_is_one_on_a_monotone_relationship() -> None:
    ascending = [1.0, 2.0, 3.0, 4.0]
    assert metrics.spearman(ascending, [10.0, 20.0, 30.0, 40.0]) == pytest.approx(1)
    assert metrics.spearman(ascending, [40.0, 30.0, 20.0, 10.0]) == pytest.approx(-1)


def test_an_uncomputable_correlation_is_nan_rather_than_zero() -> None:
    assert math.isnan(metrics.spearman([1.0, 2.0], [1.0, 2.0]))
    assert math.isnan(metrics.spearman([1.0, 1.0, 1.0], [1.0, 2.0, 3.0]))
    assert math.isnan(metrics.spearman([1.0, 2.0, np.nan], [1.0, np.nan, 3.0]))


def test_precision_at_k_is_the_overlap_of_the_two_top_ks() -> None:
    descending = [5.0, 4.0, 3.0, 2.0, 1.0]
    assert metrics.precision_at_k(descending, descending, k=2) == 1.0
    assert metrics.precision_at_k(descending, [1.0, 2.0, 3.0, 4.0, 5.0], k=2) == 0.0
    assert math.isnan(metrics.precision_at_k([1.0, 2.0], [1.0, 2.0], k=10))


def test_top_k_recall_is_a_shortlist_against_what_worked() -> None:
    """The realistic question: you order 4 variants hoping to catch the 2 good
    ones. With k == n_true it is precision@k; larger k is where they differ."""
    predicted = [5.0, 4.0, 3.0, 2.0, 1.0]
    measured = [1.0, 2.0, 5.0, 4.0, 3.0]

    assert metrics.top_k_recall(predicted, measured, k=2, n_true=2) == (
        metrics.precision_at_k(predicted, measured, k=2)
    )
    # The two best measured are at indices 2 and 3; a shortlist of 4 catches both.
    assert metrics.top_k_recall(predicted, measured, k=4, n_true=2) == 1.0
    assert math.isnan(metrics.top_k_recall(predicted, measured, k=99, n_true=2))


def test_the_mean_skips_nans_and_is_nan_when_there_are_only_nans() -> None:
    assert metrics.nanmean([1.0, float("nan"), 3.0]) == pytest.approx(2.0)
    assert math.isnan(metrics.nanmean([float("nan")]))
    assert math.isnan(metrics.nanmean([]))


def test_a_linear_scorer_agrees_in_magnitude_as_well_as_rank() -> None:
    predicted = [1.0, 2.0, 3.0, 4.0, 5.0]
    measured = [2.0, 4.0, 6.0, 8.0, 10.0]

    fit = metrics.magnitude_agreement(predicted, measured)

    assert fit["slope"] == pytest.approx(2.0)
    assert fit["intercept"] == pytest.approx(0.0)
    assert fit["r2"] == pytest.approx(1.0)
    assert fit["rank_only"] is False


def test_a_monotone_but_curved_scorer_is_flagged_rank_only() -> None:
    """The question the harness exists to answer at step 3: this column orders
    every candidate correctly and its magnitudes are meaningless, so the
    shortlist has to be cut by rank and not by a threshold."""
    predicted = [1.0, 2.0, 3.0, 4.0, 5.0, 6.0]
    measured = [1.0, 1.1, 1.2, 1.3, 1.4, 60.0]

    fit = metrics.magnitude_agreement(predicted, measured)

    assert fit["spearman"] == pytest.approx(1.0)
    assert fit["pearson"] < 0.8
    assert fit["rank_only"] is True


def test_magnitude_agreement_is_nan_on_a_constant_column() -> None:
    fit = metrics.magnitude_agreement([1.0, 1.0, 1.0], [1.0, 2.0, 3.0])

    assert math.isnan(fit["slope"]) and math.isnan(fit["r2"])
    assert fit["rank_only"] is False


def test_score_columns_skip_the_raw_ones() -> None:
    """A `_raw` column is the model's own sign convention, so correlating it
    against anything would silently flip half an ensemble."""
    import pandas as pd

    frame = pd.DataFrame(
        columns=["variant_id", "score_esm_if", "score_esm_if_raw", "measured_affinity"]
    )

    assert metrics.score_columns(frame) == ["score_esm_if"]
    assert metrics.scorer_name("score_esm_if") == "esm_if"
    assert metrics.scorer_name("oriented") == "oriented"


# --------------------------------------------------------------------------
# The schema
# --------------------------------------------------------------------------


def test_metrics_round_trip(tmp_path) -> None:
    path = metrics.save(tmp_path / "run", {"dockq": 0.51, "quality": "medium", "n": 3})

    assert path.name == metrics.FILENAME
    assert metrics.load(tmp_path / "run") == {"dockq": 0.51, "quality": "medium", "n": 3}


def test_a_nested_payload_is_refused(tmp_path) -> None:
    """Flatness is the whole point: a notebook and the harness both read any
    step's numbers without knowing which step wrote them."""
    with pytest.raises(TypeError, match="flat map"):
        metrics.save(tmp_path / "run", {"confidence": {"iptm": 0.8}})


def test_a_missing_file_says_to_use_load_all(tmp_path) -> None:
    with pytest.raises(FileNotFoundError, match="load_all"):
        metrics.load(tmp_path)


def test_load_all_sweeps_a_tree_and_skips_what_has_none(tmp_path) -> None:
    metrics.save(tmp_path / "a", {"x": 1})
    metrics.save(tmp_path / "b" / "c", {"x": 2})
    (tmp_path / "empty").mkdir()
    (tmp_path / "broken").mkdir()
    (tmp_path / "broken" / metrics.FILENAME).write_text("{not json")

    found = metrics.load_all(tmp_path)

    assert set(found) == {tmp_path / "a", tmp_path / "b" / "c"}
    assert found[tmp_path / "a"] == {"x": 1}


def test_a_scorer_run_leaves_metrics_behind(tmp_path) -> None:
    """Every step that produces numbers writes one, which is what makes
    load_all() enough to collect a campaign."""
    from design_pipeline import enumerate_variants, get_structure, score_variants

    complex_ = get_structure(tmp_path / "runs", model_name="dummy_complex")
    variants = enumerate_variants(complex_, chain="H", positions=[1, 2], mutations="A")
    score_variants(complex_, variants, models=["dummy_scorer"])

    collected = metrics.load_all(tmp_path / "runs")
    scorer_runs = [m for m in collected.values() if m.get("model_name") == "dummy_scorer"]

    assert scorer_runs and scorer_runs[0]["candidates"] == 2
    assert scorer_runs[0]["per_variant"] is False
