"""The Category 2 diagnostics - the ones that need no answer key.

Against hand-written fixtures, by pipeline step. Every case here is a bug the
diagnostic is supposed to catch, written out as the smallest table that has it:
a flipped sign, two columns that are one opinion, a shortlist that is fifty
bets on two positions, a gate that kills exactly the variants the scorer liked.

The Category 1 metrics are in test_metrics.py.
"""

from __future__ import annotations

import math
from pathlib import Path

import pandas as pd
import pytest

from design_pipeline.adapters._pdb import (
    ANTIBODY_SPACING,
    DOCKED,
    FIXTURE_CHAINS,
    MISPLACED,
    write_ca_pdb,
)
from design_pipeline.metrics import diagnostics

DOCKED_ORIGINS = {
    "H": (0.0, 0.0, 0.0),
    "L": (ANTIBODY_SPACING, 0.0, 0.0),
    "A": DOCKED,
}
WRONG_ORIGINS = {**DOCKED_ORIGINS, "A": MISPLACED}


@pytest.fixture
def complex_path(tmp_path: Path) -> Path:
    return write_ca_pdb(tmp_path / "ref.pdb", FIXTURE_CHAINS, origins=DOCKED_ORIGINS)


@pytest.fixture
def wrong_pose(tmp_path: Path) -> Path:
    """The same complex with the antigen on the opposite face."""
    return write_ca_pdb(tmp_path / "wrong.pdb", FIXTURE_CHAINS, origins=WRONG_ORIGINS)


def singles(**columns: list[float]) -> pd.DataFrame:
    """A minimal variant table: four candidates on chain H, one row each."""
    ids = ["H:E1A", "H:V2A", "H:Q3A", "H:L4A"]
    frame = pd.DataFrame(
        {
            "variant_id": ids,
            "chain": ["H"] * 4,
            "position": [1, 2, 3, 4],
            "wt": ["E", "V", "Q", "L"],
            "mut": ["A"] * 4,
        }
    )
    for name, values in columns.items():
        frame[name] = values
    return frame


# --------------------------------------------------------------------------
# The caveat itself
# --------------------------------------------------------------------------


def test_the_caveat_says_these_do_not_measure_correctness() -> None:
    """It is in the module because it is the thing most likely to be forgotten:
    a table of green diagnostics read as evidence the predictions are good."""
    assert "cannot say whether any prediction is correct" in diagnostics.CAVEAT
    assert "ground-truth" in diagnostics.CAVEAT


def test_there_is_no_function_blending_the_two_categories() -> None:
    """A single score mixing a diagnostic with a ground-truth metric would be
    read as the second, which is why the module deliberately has none."""
    assert not [
        name
        for name in dir(diagnostics)
        if name in ("overall_score", "quality_score", "combined_score", "summary_score")
    ]


# --------------------------------------------------------------------------
# Step 1: the structure
# --------------------------------------------------------------------------


def test_one_prediction_agrees_with_itself_and_says_nothing(complex_path) -> None:
    assert diagnostics.seed_agreement([complex_path]) == (1.0, 1)


def test_identical_seeds_form_one_cluster(complex_path, tmp_path) -> None:
    copy = write_ca_pdb(tmp_path / "copy.pdb", FIXTURE_CHAINS, origins=DOCKED_ORIGINS)

    assert diagnostics.seed_agreement([complex_path, copy]) == (1.0, 1)


def test_a_split_ensemble_reports_the_dominant_fraction(
    complex_path, wrong_pose, tmp_path
) -> None:
    """Two seeds on one pose and one on another: the fraction is the evidence,
    and it is the one signal here that is not the model grading itself."""
    second = write_ca_pdb(
        tmp_path / "second.pdb", FIXTURE_CHAINS, origins=DOCKED_ORIGINS
    )
    fraction, clusters = diagnostics.seed_agreement(
        [complex_path, second, wrong_pose]
    )

    assert clusters == 2
    assert fraction == pytest.approx(2 / 3)


def test_no_predictions_is_an_error() -> None:
    with pytest.raises(ValueError, match="at least one"):
        diagnostics.seed_agreement([])


def test_the_decoy_gap_is_the_difference_and_nan_when_absent() -> None:
    assert diagnostics.decoy_gap(0.8, 0.3) == pytest.approx(0.5)
    assert diagnostics.decoy_gap(0.3, 0.8) == pytest.approx(-0.5)
    assert math.isnan(diagnostics.decoy_gap(0.8, None))
    assert math.isnan(diagnostics.decoy_gap(None, 0.3))


# --------------------------------------------------------------------------
# Step 2: enumerate
# --------------------------------------------------------------------------


def test_a_full_scan_is_nineteen_candidates_per_position() -> None:
    frame = pd.DataFrame(
        {
            "variant_id": [f"H:E1{m}" for m in "ACDEFGHIKLMNPQRST"],
            "chain": ["H"] * 17,
            "position": [1] * 17,
        }
    )
    counted = diagnostics.variant_count(frame, sequence_length=1)

    assert counted["candidates"] == 17
    assert counted["expected"] == 19
    assert counted["complete"] is False


def test_variant_count_says_nothing_about_completeness_with_nothing_to_compare() -> None:
    counted = diagnostics.variant_count(singles())

    assert counted["candidates"] == 4
    assert counted["positions"] == 4
    assert counted["expected"] is None and counted["complete"] is None


def test_a_short_scan_can_be_made_to_raise() -> None:
    """For a caller that wants the run to stop: a shortlist chosen from a set
    nobody intended is invisible downstream, because every later step is
    perfectly happy with a short table."""
    with pytest.raises(ValueError, match="full scan"):
        diagnostics.variant_count(singles(), sequence_length=4, strict=True)

    # And a table that does add up passes strict unchanged.
    assert diagnostics.variant_count(
        singles(), sequence_length=4, alphabet=1, strict=True
    )["complete"] is True


def test_the_position_buckets_partition_the_chain(complex_path) -> None:
    buckets = diagnostics.position_buckets(complex_path, "H")
    classes = (
        "interface",
        "second_shell",
        "cdr_non_contact",
        "vh_vl_interface",
        "buried_core",
        "surface_framework",
    )

    assert sum(buckets[name] for name in classes) == buckets["positions"]
    assert buckets["positions"] == len(FIXTURE_CHAINS["H"])
    # The fixture antigen is docked against the antibody, so there is one.
    assert buckets["interface"] > 0
    # And the light chain sits next to the heavy one.
    assert buckets["vh_vl_interface"] > 0


def test_cdr_positions_are_declared_rather_than_guessed(complex_path) -> None:
    """Deriving them needs ANARCI, and nothing in this module reaches for a
    tool - so without them the bucket is empty and says so."""
    without = diagnostics.position_buckets(complex_path, "H")
    assert without["cdr_known"] is False
    assert without["cdr_non_contact"] == 0

    with_cdr = diagnostics.position_buckets(
        complex_path, "H", cdr_positions=range(1, len(FIXTURE_CHAINS["H"]) + 1)
    )
    assert with_cdr["cdr_known"] is True
    # Every non-contact position is now a CDR position, so nothing is left
    # anywhere below CDR in the precedence order.
    assert with_cdr["buried_core"] == 0
    assert with_cdr["surface_framework"] == 0
    assert with_cdr["interface"] == without["interface"]


def test_a_chain_that_is_not_the_antibody_is_an_error(complex_path) -> None:
    with pytest.raises(KeyError, match="not on the antibody side"):
        diagnostics.position_buckets(complex_path, "A")


def test_a_docked_complex_has_an_interface_and_a_separated_one_does_not(
    complex_path, tmp_path
) -> None:
    """The cheapest place in the pipeline to catch a wrong chain assignment:
    the interface is computed before anything expensive runs on it."""
    docked = diagnostics.interface_size(complex_path)

    assert docked["n_residues_side"] > 0 and docked["n_residues_partner"] > 0
    assert docked["buried_area"] > 0
    assert docked["buried_area_per_side"] == pytest.approx(docked["buried_area"] / 2)
    assert docked["side"] == "H,L" and docked["partner"] == "A"


def test_a_four_residue_interface_is_reported_as_implausible(complex_path) -> None:
    """Not an unusual epitope - a wrong pose or a wrong chain. The fixture's
    contact is small, so the plausibility floor is what is being tested."""
    tiny = diagnostics.interface_size(complex_path, min_residues=1000)
    assert tiny["plausible"] is False

    generous = diagnostics.interface_size(complex_path, min_residues=1)
    assert generous["plausible"] is True


# --------------------------------------------------------------------------
# Step 3: score
# --------------------------------------------------------------------------


def test_two_columns_that_are_one_opinion_correlate_near_one() -> None:
    """The single most informative diagnostic: an ensemble of two columns at
    0.8 is an ensemble of one, and the width of the table hides it."""
    frame = singles(
        score_a=[1.0, 2.0, 3.0, 4.0],
        score_b=[1.1, 2.2, 3.3, 4.4],
        score_c=[4.0, 3.0, 2.0, 1.0],
    )
    pairwise = diagnostics.inter_model_spearman(
        frame, ["score_a", "score_b", "score_c"]
    )

    assert pairwise[("a", "b")] == pytest.approx(1.0)
    assert pairwise[("a", "c")] == pytest.approx(-1.0)
    # Keyed by sorted scorer names, so column order cannot change the keys.
    assert set(pairwise) == {("a", "b"), ("a", "c"), ("b", "c")}


def test_the_score_columns_are_found_without_being_named() -> None:
    frame = singles(score_a=[1.0, 2.0, 3.0, 4.0], score_b=[1.0, 2.0, 3.0, 4.0])

    assert set(diagnostics.inter_model_spearman(frame)) == {("a", "b")}


def test_a_column_name_that_does_not_exist_is_an_error() -> None:
    with pytest.raises(KeyError, match="score_typo"):
        diagnostics.inter_model_spearman(singles(), ["score_typo"])


def test_a_mis_signed_column_calls_most_mutations_beneficial() -> None:
    """Most mutations to a matured antibody are deleterious, so a column
    calling three quarters of them improving has its sign the wrong way round -
    caught here rather than after it inverts the shortlist."""
    frame = singles(
        score_sane=[-1.0, -2.0, -0.5, 1.0],
        score_flipped=[1.0, 2.0, 0.5, -1.0],
    )
    shares = diagnostics.beneficial_fraction(frame)

    assert shares["sane"] == pytest.approx(0.25)
    assert shares["flipped"] == pytest.approx(0.75)
    assert shares["flipped"] > diagnostics.PLAUSIBLE_BENEFICIAL


def test_a_column_that_declined_has_no_beneficial_fraction() -> None:
    frame = singles(score_absent=[float("nan")] * 4)

    assert math.isnan(diagnostics.beneficial_fraction(frame)["absent"])


def test_a_matrix_scorer_should_prefer_the_wild_type() -> None:
    """On a molecule maturation already optimised. A low rate invalidates the
    run rather than weakening it: the model is not reading this antibody."""
    frame = pd.DataFrame(
        {
            "variant_id": ["H:E1E", "H:E1A", "H:V2V", "H:V2A"],
            "chain": ["H"] * 4,
            "position": [1, 1, 2, 2],
            "wt": ["E", "E", "V", "V"],
            "mut": ["E", "A", "V", "A"],
            "score_good": [2.0, 1.0, 2.0, 1.0],
            "score_bad": [1.0, 2.0, 1.0, 2.0],
        }
    )

    assert diagnostics.wt_argmax_rate(frame, "score_good")["rate"] == 1.0
    assert diagnostics.wt_argmax_rate(frame, "score_bad")["rate"] == 0.0


def test_a_table_without_wild_type_rows_cannot_answer_the_argmax_question() -> None:
    """`enumerate_variants` skips wt-to-wt because it is not a variant, so this
    needs the scorer's own matrix and says so rather than guessing."""
    rate = diagnostics.wt_argmax_rate(singles(score_a=[1.0, 2.0, 3.0, 4.0]), "score_a")

    assert rate["skipped"] == 4 and rate["scored"] == 0
    assert math.isnan(rate["rate"])


def test_a_hundred_mutations_on_two_positions_are_two_bets() -> None:
    """The number the diversity constraint is actually about: `effective_
    positions`, not the length of the list."""
    clustered = pd.DataFrame(
        {
            "variant_id": [f"H:E1{m}" for m in "ACDEF"] + [f"H:V2{m}" for m in "ACDEF"],
            "chain": ["H"] * 10,
            "position": [1] * 5 + [2] * 5,
            "score_a": [1.0] * 10,
        }
    )
    spread = pd.DataFrame(
        {
            "variant_id": [f"H:E{i}A" for i in range(1, 11)],
            "chain": ["H"] * 10,
            "position": list(range(1, 11)),
            "score_a": [1.0] * 10,
        }
    )

    tight = diagnostics.positional_entropy(clustered, "score_a", top=10)
    wide = diagnostics.positional_entropy(spread, "score_a", top=10)

    assert tight["positions"] == 2
    assert tight["effective_positions"] == pytest.approx(2.0)
    assert tight["entropy"] == pytest.approx(1.0)
    assert wide["effective_positions"] == pytest.approx(10.0)
    assert wide["normalized"] == pytest.approx(1.0)


def test_positional_entropy_takes_only_the_top_of_the_list() -> None:
    frame = singles(score_a=[4.0, 3.0, 2.0, 1.0])
    top = diagnostics.positional_entropy(frame, "score_a", top=2)

    assert top["n"] == 2 and top["positions"] == 2
    assert top["top_position_share"] == pytest.approx(0.5)


def test_scores_that_move_with_the_seed_are_measuring_the_structure() -> None:
    """The link back to step 1, and the question step 1 cannot answer: whether
    the structure was good enough to score on at all."""
    stable = {
        "seed0": singles(score_a=[1.0, 2.0, 3.0, 4.0]),
        "seed1": singles(score_a=[1.1, 2.1, 3.1, 4.1]),
    }
    unstable = {
        "seed0": singles(score_a=[1.0, 2.0, 3.0, 4.0]),
        "seed1": singles(score_a=[4.0, 3.0, 2.0, 1.0]),
    }

    assert diagnostics.seed_rank_stability(stable, "score_a")["mean_spearman"] == (
        pytest.approx(1.0)
    )
    shaky = diagnostics.seed_rank_stability(unstable, "score_a")
    assert shaky["mean_spearman"] == pytest.approx(-1.0)
    assert shaky["shared_candidates"] == 4
    assert shaky["seeds"] == 2


def test_one_seed_cannot_disagree_with_itself() -> None:
    with pytest.raises(ValueError, match="at least two"):
        diagnostics.seed_rank_stability([singles(score_a=[1.0, 2.0, 3.0, 4.0])], "score_a")


# --------------------------------------------------------------------------
# Step 4: the gate
# --------------------------------------------------------------------------


def test_the_pass_rate_is_per_check_and_never_pooled() -> None:
    """A liabilities check dropping most of a scan means the regex is wrong,
    and a pooled survival rate would hide which check did it."""
    frame = singles()
    frame["pass_liabilities"] = [True, True, True, False]
    frame["pass_charge"] = [True, False, False, False]

    rates = diagnostics.pass_rate(frame)

    assert rates == {"liabilities": pytest.approx(0.75), "charge": pytest.approx(0.25)}
    assert diagnostics.check_names(frame) == ["charge", "liabilities"]


def test_a_check_nobody_ran_is_an_error() -> None:
    with pytest.raises(KeyError, match="pass_absent"):
        diagnostics.pass_rate(singles(), ["absent"])


def test_a_gate_that_kills_the_best_variants_correlates_negatively() -> None:
    """A real tension between binding and developability, surfaced rather than
    silently absorbed: what comes out is then the second-best-binding set.

    The check kills exactly the top two, which is as adverse as a gate can be,
    and the correlation is -2/sqrt(5) rather than -1: a verdict column is two
    blocks of tied ranks, so a perfect inversion cannot reach the full -1. The
    sign is what is being read here, not the magnitude.
    """
    frame = singles(score_a=[4.0, 3.0, 2.0, 1.0])
    frame["pass_liabilities"] = [False, False, True, True]

    correlations = diagnostics.gate_score_correlation(frame)

    assert correlations[("liabilities", "a")] == pytest.approx(-2 / 5**0.5)


def test_a_check_that_passed_everything_has_no_correlation() -> None:
    frame = singles(score_a=[4.0, 3.0, 2.0, 1.0])
    frame["pass_liabilities"] = [True] * 4

    assert math.isnan(
        diagnostics.gate_score_correlation(frame)[("liabilities", "a")]
    )


# --------------------------------------------------------------------------
# Step 5: pairs
# --------------------------------------------------------------------------


def pair_table(scores: list[float], distances: list[float]) -> pd.DataFrame:
    ids = ["H:E1A+H:Q3A", "H:V2A+H:L4A"]
    return pd.DataFrame(
        {
            "variant_id": ids,
            "chain": ["H"] * 2,
            "position": [1, 2],
            "score_a": scores,
            "distance": distances,
        }
    )


def test_an_additive_scorer_shows_no_epistasis() -> None:
    """If this is near zero across every pair, the scorer is additive by
    construction and the second, expensive scoring run bought nothing."""
    parents = singles(score_a=[1.0, 2.0, 3.0, 4.0])
    # 1 + 3 and 2 + 4, exactly.
    pairs = pair_table([4.0, 6.0], [12.0, 15.0])

    magnitude = diagnostics.epistasis_magnitude(pairs, parents, "score_a")

    assert magnitude["pairs"] == 2
    assert magnitude["mean"] == pytest.approx(0.0)
    assert magnitude["max_abs"] == pytest.approx(0.0)


def test_a_non_additive_pair_is_measured_directly() -> None:
    parents = singles(score_a=[1.0, 2.0, 3.0, 4.0])
    pairs = pair_table([10.0, 6.0], [12.0, 15.0])

    table = diagnostics.epistasis(pairs, parents, "score_a")

    assert list(table["additive"]) == [4.0, 6.0]
    assert list(table["epistasis"]) == [6.0, 0.0]
    assert diagnostics.epistasis_magnitude(pairs, parents, "score_a")["max_abs"] == 6.0


def test_a_pair_whose_singles_are_missing_is_dropped_not_halved() -> None:
    """A sum over one mutation would read as an epistasis equal to the missing
    single's score, which is worse than no answer."""
    parents = singles(score_a=[1.0, float("nan"), 3.0, 4.0])
    pairs = pair_table([4.0, 6.0], [12.0, 15.0])

    table = diagnostics.epistasis(pairs, parents, "score_a")

    assert list(table["variant_id"]) == ["H:E1A+H:Q3A"]


def test_distance_predicts_additivity_when_the_pairing_rule_works() -> None:
    """The 10 A rule is a premise about our own data, and this is the test of
    it: bigger distance, smaller departure from additivity."""
    parents = singles(score_a=[1.0, 2.0, 3.0, 4.0])
    pairs = pd.DataFrame(
        {
            "variant_id": [
                "H:E1A+H:Q3A",
                "H:V2A+H:L4A",
                "H:E1A+H:V2A",
                "H:Q3A+H:L4A",
            ],
            "score_a": [4.0, 6.0, 4.0, 10.0],
            "distance": [30.0, 20.0, 12.0, 10.0],
        }
    )

    tested = diagnostics.distance_vs_epistasis(pairs, parents, "score_a")

    assert tested["pairs"] == 4
    assert tested["spearman_abs"] < 0


def test_a_pairs_table_with_no_distances_cannot_answer_it() -> None:
    parents = singles(score_a=[1.0, 2.0, 3.0, 4.0])
    pairs = pair_table([4.0, 6.0], [float("nan")] * 2)

    tested = diagnostics.distance_vs_epistasis(pairs, parents, "score_a")

    assert tested["pairs"] == 0 and math.isnan(tested["spearman_abs"])


# --------------------------------------------------------------------------
# Step 6: the shortlist
# --------------------------------------------------------------------------


def test_fifty_candidates_on_two_positions_are_two_bets() -> None:
    frame = pd.DataFrame(
        {
            "variant_id": [f"H:E1{m}" for m in "ACDEF"] + ["H:V2A"],
            "chain": ["H"] * 6,
            "position": [1] * 5 + [2],
        }
    )
    diversity = diagnostics.position_diversity(frame)

    assert diversity["candidates"] == 6
    assert diversity["positions"] == 2
    assert diversity["max_per_position"] == 5
    assert diversity["largest_position_share"] == pytest.approx(5 / 6)


def test_the_chain_balance_is_the_share_of_mutations_on_each_chain() -> None:
    frame = pd.DataFrame(
        {
            "variant_id": ["H:E1A", "H:V2A", "H:Q3A", "L:D1A"],
            "chain": ["H", "H", "H", "L"],
            "position": [1, 2, 3, 1],
        }
    )

    assert diagnostics.chain_balance(frame) == {
        "H": pytest.approx(0.75),
        "L": pytest.approx(0.25),
    }


def test_the_edit_distance_histogram_separates_singles_from_pairs() -> None:
    """It catches the table that was meant to be the paired set and is still
    all singles."""
    frame = pd.DataFrame(
        {
            "variant_id": ["H:E1A", "H:V2A+H:L4A", "H:V2A+H:L4A", "H:Q3A"],
            "chain": ["H"] * 4,
            "position": [1, 2, 4, 3],
        }
    )

    assert diagnostics.edit_distance_histogram(frame) == {1: 2, 2: 1}


def test_the_scorers_can_agree_about_the_shortlist_and_not_the_scan() -> None:
    """A different question from `inter_model_spearman` over the whole scan,
    and the one that decides whether the submission's composition is an
    artifact of which column it was sorted on."""
    shortlist = singles(
        score_a=[1.0, 2.0, 3.0, 4.0],
        score_b=[1.0, 2.0, 3.0, 4.0],
        score_c=[4.0, 3.0, 2.0, 1.0],
    )
    agreement = diagnostics.final_rank_agreement(shortlist)

    assert agreement["scorers"] == 3 and agreement["pairs"] == 3
    assert agreement["max_spearman"] == pytest.approx(1.0)
    assert agreement["min_spearman"] == pytest.approx(-1.0)
    assert agreement["mean_spearman"] == pytest.approx(-1 / 3)


# --------------------------------------------------------------------------
# Into the schema
# --------------------------------------------------------------------------


def test_a_diagnostic_flattens_into_the_metrics_schema(tmp_path) -> None:
    """The shapes these return suit the analysis; `metrics.json` is flat. The
    translation happens once, here, not in every caller."""
    from design_pipeline import metrics

    frame = singles(score_a=[1.0, 2.0, 3.0, 4.0], score_b=[1.0, 2.0, 3.0, 4.0])
    flat = diagnostics.flatten(diagnostics.inter_model_spearman(frame), "redundancy")

    assert flat == {"redundancy_a_b": pytest.approx(1.0)}
    # And it is accepted by the schema, which refuses anything nested.
    metrics.save(tmp_path / "run", flat)
    assert metrics.load(tmp_path / "run")["redundancy_a_b"] == pytest.approx(1.0)


def test_flatten_walks_nested_maps_and_tuple_keys() -> None:
    flat = diagnostics.flatten(
        {"pass": {"liabilities": 0.75}, ("gate", "esm_if"): -1.0, "n": 4}
    )

    assert flat == {"pass_liabilities": 0.75, "gate_esm_if": -1.0, "n": 4}
