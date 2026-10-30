"""Structure confidence: what gets recorded, what warns, what it recommends."""

from __future__ import annotations

import json
import warnings
from contextlib import contextmanager
from pathlib import Path

import pytest

from design_pipeline import confidence, get_structure, metrics
from design_pipeline.adapters._pdb import FIXTURE_CHAINS

ANTIBODY = f"{FIXTURE_CHAINS['H']}:{FIXTURE_CHAINS['L']}"
ANTIGEN = FIXTURE_CHAINS["A"]
OTHER_ANTIGEN = "MKTIIALSYIFCLVFAQKLPGNDNSTATLCLGHHAVPNGTIVKTITNDQ"


@contextmanager
def quiet():
    """The fixture predictor trips the low-confidence warnings by design; the
    test that cares about them asserts them directly."""
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", RuntimeWarning)
        yield


def predict(tmp_path: Path, **kwargs):
    with quiet():
        return get_structure(
            tmp_path / "runs",
            model_name="dummy_predictor",
            antibody_seq=ANTIBODY,
            antigen_seq=ANTIGEN,
            chain_ids=("H", "L", "A"),
            **kwargs,
        )


def test_a_prediction_records_what_the_model_said_and_what_the_seeds_did(tmp_path):
    complex_ = predict(tmp_path, num_seeds=8)
    recorded = complex_.confidence

    assert recorded is not None
    assert 0.0 <= recorded["iptm"] <= 1.0
    assert recorded["n_seeds"] == 8
    assert 0.0 < recorded["dominant_cluster_fraction"] <= 1.0
    assert recorded["n_clusters"] >= 1
    assert recorded["mean_interchain_pae"] is not None
    assert recorded["experimental"] is False


def test_plddt_is_recorded_and_marked_uninformative(tmp_path) -> None:
    """Both chains fold well individually whether or not the docking is right,
    so pLDDT on a complex is close to a constant. Kept for the record and
    flagged, because somebody will otherwise read it as a docking score."""
    recorded = predict(tmp_path, num_seeds=2).confidence

    assert recorded["plddt"] is not None
    assert recorded["plddt_uninformative"] is True
    assert "plddt" in confidence.UNINFORMATIVE


def test_the_ensemble_is_clustered_not_averaged(tmp_path) -> None:
    """The documented bottleneck is identifying the right model among those
    generated, so the spread is the evidence, not the mean."""
    complex_ = predict(tmp_path, num_seeds=6)

    assert len(complex_.paths) == 6
    fraction, clusters = metrics.seed_agreement(
        list(complex_.paths), ["H", "L"], ["A"]
    )
    assert complex_.confidence["dominant_cluster_fraction"] == pytest.approx(fraction)
    assert complex_.confidence["n_clusters"] == clusters


def test_the_decoy_control_is_run_on_identical_settings(tmp_path) -> None:
    complex_ = predict(tmp_path, num_seeds=4, decoy_antigen_seq=OTHER_ANTIGEN)
    recorded = complex_.confidence

    assert recorded["decoy_iptm"] is not None
    assert recorded["decoy_gap"] == pytest.approx(
        recorded["iptm"] - recorded["decoy_iptm"]
    )
    # Its own cache entry, beside the real one rather than inside it.
    assert (tmp_path / "runs" / "decoy").is_dir()


def test_a_decoy_without_a_known_antibody_split_is_refused(tmp_path) -> None:
    with pytest.raises(ValueError, match="decoy_antigen_seq needs antibody_seq"):
        get_structure(
            tmp_path / "runs",
            model_name="dummy_predictor",
            fasta=None,
            decoy_antigen_seq=OTHER_ANTIGEN,
        )


def test_confidence_is_written_beside_the_structure_and_into_metrics(tmp_path):
    complex_ = predict(tmp_path, num_seeds=3)
    run_dir = complex_.path.parent.parent

    stored = json.loads((run_dir / confidence.FILENAME).read_text())
    assert stored["iptm"] == complex_.confidence["iptm"]
    # Mirrored into the one schema every step writes.
    assert metrics.load(run_dir)["confidence_iptm"] == complex_.confidence["iptm"]


def test_a_cache_hit_reuses_the_recorded_confidence(tmp_path) -> None:
    """The seed clustering is pairwise over the ensemble, so recomputing it on
    every hit would make a resumed campaign slower than a fresh one."""
    first = predict(tmp_path, num_seeds=6)
    second = predict(tmp_path, num_seeds=6)

    assert first.confidence == second.confidence


def test_it_can_be_skipped(tmp_path) -> None:
    complex_ = predict(tmp_path, num_seeds=2, confidence=False)

    assert complex_.confidence is None
    assert confidence.recommend_weights(None) # a dict, not a crash


def test_an_experimental_structure_has_no_model_scores(tmp_path) -> None:
    from design_pipeline.adapters._pdb import write_ca_pdb

    source = write_ca_pdb(tmp_path / "x.pdb", FIXTURE_CHAINS)
    complex_ = get_structure(tmp_path / "runs", source=source)

    assert complex_.confidence["experimental"] is True
    assert complex_.confidence["iptm"] is None
    assert complex_.confidence["mean_interchain_pae"] is None


# --------------------------------------------------------------------------
# The rules
# --------------------------------------------------------------------------


def test_the_warnings_only_ever_point_down() -> None:
    """Low ipTM reliably means do not trust this; high ipTM does not reliably
    mean this is correct. So there is no branch that says a structure is good."""
    with pytest.warns(RuntimeWarning, match="ipTM is below"):
        assert confidence.warn({"iptm": 0.2}) != []
    with pytest.warns(RuntimeWarning, match="seeds agree"):
        assert confidence.warn({"dominant_cluster_fraction": 0.2}) != []
    with pytest.warns(RuntimeWarning, match="decoy"):
        assert confidence.warn({"decoy_gap": 0.01}) != []
    # Nothing fires upward: there is no warning, and no claim, for a good one.
    with warnings.catch_warnings():
        warnings.simplefilter("error")
        assert confidence.warn({"iptm": 0.95, "dominant_cluster_fraction": 1.0}) == []


def test_a_missing_signal_is_not_a_warning() -> None:
    """Most of these fields can legitimately be absent - no seeds, no decoy, a
    downloaded structure - and treating absence as failure would fire on every
    honest run."""
    assert confidence.warn({}) == []
    assert confidence.warn({"iptm": None, "decoy_gap": None}) == []


@pytest.mark.parametrize(
    "recorded, expected",
    [
        ({"iptm": 0.9, "dominant_cluster_fraction": 0.9}, 1.0),
        ({"iptm": 0.3}, 0.4),
        ({"dominant_cluster_fraction": 0.3}, 0.3),
        ({"decoy_gap": 0.02}, 0.2),
        # The lowest triggered weight wins.
        ({"iptm": 0.3, "decoy_gap": 0.02}, 0.2),
    ],
)
def test_the_rule_table_downweights_structure_based_scorers(recorded, expected):
    weights = confidence.recommend_weights(recorded, ["antifold", "gearbind"])

    assert weights == {"antifold": expected, "gearbind": expected}


def test_sequence_scorers_come_up_when_the_structure_is_doubted() -> None:
    """If the complex is a fiction, they are the only columns left that mean
    anything; when it is trustworthy they sit below the structure-based ones,
    because they know nothing about the interface."""
    trusted = confidence.recommend_weights({"iptm": 0.9}, ["antifold", "esm1v"])
    doubted = confidence.recommend_weights({"iptm": 0.2}, ["antifold", "esm1v"])

    assert trusted["antifold"] > trusted["esm1v"]
    assert doubted["esm1v"] > doubted["antifold"]
    assert doubted["esm1v"] == 1.0


def test_it_recommends_and_does_not_apply(tmp_path) -> None:
    """Deciding what to rank on stays the caller's. Nothing in the pipeline
    reads these weights."""
    complex_ = predict(tmp_path, num_seeds=4)
    weights = confidence.recommend_weights(complex_.confidence)

    assert isinstance(weights, dict)
    assert all(isinstance(value, float) for value in weights.values())
    # The variant table has no weight column and no combined score.
    from design_pipeline import enumerate_variants, score_variants

    variants = enumerate_variants(complex_, chain="H", positions=[1], mutations="A")
    scored = score_variants(complex_, variants, models=["dummy_scorer"])
    assert not any("weight" in column for column in scored.columns)


def test_with_no_confidence_it_declines_to_guess() -> None:
    weights = confidence.recommend_weights(None, ["antifold", "esm1v"])

    assert weights == {"antifold": 1.0, "esm1v": 1.0}
    assert confidence.weight_reasons(None) == ["no confidence recorded for this structure"]


def test_an_unknown_scorer_is_assumed_to_need_the_structure() -> None:
    """The safe default: one that turns out to need coordinates would
    otherwise keep full weight on a structure nobody trusts."""
    weights = confidence.recommend_weights({"iptm": 0.2}, ["not_registered"])

    assert weights["not_registered"] == 0.4


def test_the_reasons_are_readable(tmp_path) -> None:
    reasons = confidence.weight_reasons({"decoy_gap": 0.01, "iptm": 0.1})

    assert len(reasons) == 2
    assert any("decoy" in reason for reason in reasons)
    assert any("ipTM" in reason for reason in reasons)


def test_read_returns_none_when_nothing_was_written(tmp_path) -> None:
    assert confidence.read(tmp_path) is None
    (tmp_path / confidence.FILENAME).write_text("{not json")
    assert confidence.read(tmp_path) is None
