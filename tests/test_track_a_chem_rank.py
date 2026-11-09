"""examples/track_a_chem_rank.py: the ranking rules on synthetic rows. Scoring real folds needs Boltz output and is not run here."""

from __future__ import annotations

import importlib.util
import io
import math
import sys
from pathlib import Path

import pytest

SCRIPT = Path(__file__).resolve().parent.parent / "examples" / "track_a_chem_rank.py"
_spec = importlib.util.spec_from_file_location("track_a_chem_rank", SCRIPT)
chem = importlib.util.module_from_spec(_spec)
sys.modules["track_a_chem_rank"] = chem
_spec.loader.exec_module(chem)


def row(name, seed, ipsae, eligible=1.0, contact=0.0, coulomb=0.0):
    return {"name": name, "seed": seed, "ipsae_max": ipsae, "ipsae_min": ipsae - 0.1, "iptm": 0.8, "eligible": eligible,
            "contact": contact, "coulomb": coulomb, "paratope_net_charge": 0.0}


PARENT_ROWS = [row("parent", s, v) for s, v in enumerate((0.80, 0.70, 0.75))]     # mean 0.75, SD 0.05, floor 0.65


def test_candidates_inside_the_band_are_ranked_by_chemistry_and_the_rest_follow_by_ipsae():
    rows = PARENT_ROWS + [
        row("good", 0, 0.74, contact=1.0, coulomb=2.0), row("good", 1, 0.74, contact=1.0, coulomb=2.0),
        row("plain", 0, 0.76, contact=0.0, coulomb=0.0), row("plain", 1, 0.76, contact=0.0, coulomb=0.0),
        row("low", 0, 0.30, contact=5.0, coulomb=9.0), row("low", 1, 0.30, contact=5.0, coulomb=9.0),
    ]
    ranked, band = chem.rank_rows(rows)
    assert band["floor"] == pytest.approx(0.65)
    assert [o["name"] for o in ranked] == ["good", "plain", "low"]       # chemistry decides inside the band, ipSAE outside
    assert not ranked[2]["in_band"] and math.isnan(ranked[2]["chem"])


def test_only_on_site_seeds_are_averaged_and_a_candidate_with_none_is_left_to_ipsae():
    rows = PARENT_ROWS + [
        row("a", 0, 0.74, eligible=1.0, contact=1.0, coulomb=1.0), row("a", 1, 0.74, eligible=0.0, contact=9.0, coulomb=9.0),
        row("b", 0, 0.72, eligible=0.0, contact=3.0, coulomb=3.0), row("b", 1, 0.72, eligible=0.0, contact=3.0, coulomb=3.0),
        row("c", 0, 0.74, eligible=1.0, contact=0.0, coulomb=0.0),
    ]
    ranked, _ = chem.rank_rows(rows)
    by = {o["name"]: o for o in ranked}
    assert by["a"]["contact"] == 1.0 and by["a"]["n_on_site"] == 1           # the off-site seed's 9.0 is not averaged in
    assert by["b"]["n_on_site"] == 0 and math.isnan(by["b"]["chem"])
    assert ranked[-1]["name"] == "b"


def test_the_band_needs_the_parent_at_two_seeds():
    with pytest.raises(SystemExit):
        chem.rank_rows([row("parent", 0, 0.7), row("x", 0, 0.7)])


def test_a_small_chemistry_weight_reorders_near_ties_and_cannot_beat_a_large_ipsae_gap():
    rows = PARENT_ROWS + [
        row("top", 0, 0.90, contact=0.0, coulomb=0.0), row("top", 1, 0.90, contact=0.0, coulomb=0.0),
        row("near", 0, 0.60, contact=0.0, coulomb=0.0), row("near", 1, 0.60, contact=0.0, coulomb=0.0),
        row("near_chem", 0, 0.60, contact=1.0, coulomb=1.0), row("near_chem", 1, 0.60, contact=1.0, coulomb=1.0),
        row("low_chem", 0, 0.20, contact=9.0, coulomb=9.0), row("low_chem", 1, 0.20, contact=9.0, coulomb=9.0),
    ]
    ranked, _ = chem.rank_weighted(rows, 0.15)
    assert [o["name"] for o in ranked][:3] == ["top", "near_chem", "near"]       # chemistry breaks the tie, not the gap
    assert ranked[-1]["name"] == "low_chem"
    unweighted, _ = chem.rank_weighted(rows, 0.0)
    assert unweighted[0]["name"] == "top"


# ------------------------------------------------------------------ the score's inputs and the correlation audit


def test_ipsae_min_is_reported_but_never_enters_the_score():
    """The guard for the exclusion: perturbing ipsae_min must move the column and nothing else."""
    rows = PARENT_ROWS + [
        row("a", 0, 0.80, contact=1.0, coulomb=1.0), row("a", 1, 0.80, contact=1.0, coulomb=1.0),
        row("b", 0, 0.70, contact=2.0, coulomb=2.0), row("b", 1, 0.70, contact=2.0, coulomb=2.0),
        row("c", 0, 0.60, contact=0.0, coulomb=0.0), row("c", 1, 0.60, contact=0.0, coulomb=0.0),
    ]
    base, _ = chem.rank_weighted([dict(r) for r in rows], 0.15)
    shuffled, _ = chem.rank_weighted([dict(r, ipsae_min=0.9 - 0.3 * i) for i, r in enumerate(rows)], 0.15)
    assert [o["name"] for o in base] == [o["name"] for o in shuffled]
    assert [o["score"] for o in base] == [pytest.approx(o["score"]) for o in shuffled]
    assert {o["name"]: o["ipsae_min"] for o in base} != {o["name"]: o["ipsae_min"] for o in shuffled}
    assert "ipsae_min" not in chem.SCORE_INPUTS and set(chem.SCORE_INPUTS) <= set(base[0])


def test_the_score_inputs_are_exactly_the_columns_the_score_reads():
    """Freeze the formula's inputs, so a column added to `score` without being added to SCORE_INPUTS fails here."""
    assert chem.SCORE_INPUTS == ("ipsae_max", "contact", "coulomb")
    assert set(chem.score_weights(0.15)) == set(chem.SCORE_INPUTS)
    assert chem.score_weights(0.15) == {"ipsae_max": 1.0, "contact": 0.15, "coulomb": 0.15}


def records(**columns):
    return [dict(zip(columns, values)) for values in zip(*columns.values())]


def test_a_duplicated_score_input_is_warned_on_and_an_independent_one_is_not():
    duplicated = records(ipsae_max=[0.1, 0.3, 0.5, 0.7, 0.9], contact=[0.2, 0.6, 1.0, 1.4, 1.8],
                         coulomb=[1.0, -2.0, 0.5, 3.0, -1.0])
    lines = chem.check_score_correlations(duplicated, stream=io.StringIO())
    assert len(lines) == 1 and "ipsae_max and contact" in lines[0]
    assert chem.check_score_correlations(records(ipsae_max=[0.1, 0.9, 0.4, 0.5, 0.2],
                                                contact=[1.0, 0.2, 0.9, 0.1, 0.5],
                                                coulomb=[-1.0, 0.3, 2.0, -0.4, 0.8]), stream=io.StringIO()) == []


def test_the_warning_goes_to_the_stream_and_names_the_threshold():
    stream = io.StringIO()
    chem.check_score_correlations(records(ipsae_max=[1.0, 2.0, 3.0, 4.0], contact=[2.0, 4.0, 6.0, 8.0],
                                         coulomb=[0.0, 1.0, 0.0, 1.0]), 0.9, stream=stream)
    text = stream.getvalue()
    assert "0.90" in text and "ipsae_max and contact" in text and "r +1.000" in text


def test_spearman_catches_a_monotone_pair_that_pearson_softens():
    curved = records(ipsae_max=[1.0, 2.0, 3.0, 4.0, 5.0], contact=[1.0, 2.0, 4.0, 8.0, 160.0],
                     coulomb=[0.0, 1.0, 0.0, 1.0, 0.0])
    pair = ("ipsae_max", "contact")
    assert chem.correlation_matrix(curved, chem.SCORE_INPUTS, "pearson")[pair] < 0.9
    assert chem.correlation_matrix(curved, chem.SCORE_INPUTS, "spearman")[pair] == pytest.approx(1.0)
    assert chem.check_score_correlations(curved, stream=io.StringIO())      # rho alone is enough to flag it


def test_the_matrix_reads_csv_strings_bools_and_blanks_and_is_pairwise_complete():
    rows = [{"ipsae_max": "0.80", "in_band": "True", "contact": "1.0"},
            {"ipsae_max": "0.70", "in_band": "False", "contact": ""},
            {"ipsae_max": "0.60", "in_band": "False", "contact": "nan"},
            {"ipsae_max": "0.50", "in_band": False, "contact": "0.0"}]
    matrix = chem.correlation_matrix(rows, ("ipsae_max", "in_band", "contact"))
    assert matrix[("ipsae_max", "in_band")] == pytest.approx(0.775, abs=1e-3)       # True/False read as 1 and 0
    assert math.isnan(matrix[("ipsae_max", "contact")])                             # two usable rows: not comparable
    assert len(chem._both(rows, "ipsae_max", "in_band")) == 4


def test_a_constant_or_absent_column_is_not_comparable_rather_than_an_error():
    flat = records(ipsae_max=[0.1, 0.2, 0.3, 0.4], contact=[1.0, 1.0, 1.0, 1.0])
    assert math.isnan(chem.correlation_matrix(flat, ("ipsae_max", "contact"))[("ipsae_max", "contact")])
    assert chem.correlation_matrix(flat, ("ipsae_max", "para_net")) == {}            # absent column is skipped
    assert chem.check_score_correlations(flat, stream=io.StringIO()) == []


def test_effective_weights_sum_to_one_and_shift_when_a_column_is_dropped():
    rows = records(ipsae_max=[0.1, 0.3, 0.5, 0.7, 0.9], contact=[0.2, 0.6, 1.0, 1.4, 1.8],
                   coulomb=[1.0, -2.0, 0.5, 3.0, -1.0])
    before = chem.effective_weights(rows, 0.15)
    assert sum(before.values()) == pytest.approx(1.0)
    assert before["ipsae_max"] > before["contact"] > 0.0
    after = chem.effective_weights(rows, 0.15, drop=("contact",))
    assert "contact" not in after and sum(after.values()) == pytest.approx(1.0)
    assert after["ipsae_max"] > before["ipsae_max"]           # the duplicated chemistry term is gone


def test_the_formatted_matrix_marks_the_pairs_over_the_limit():
    text = chem.format_correlation_matrix(records(ipsae_max=[1.0, 2.0, 3.0, 4.0], contact=[2.0, 4.0, 6.0, 8.0],
                                                  coulomb=[0.0, 1.0, 0.0, 1.0]), chem.SCORE_INPUTS)
    assert "+1.000*" in text and text.splitlines()[0].split()[:1] == ["ipsae_ma"]


def test_rank_warns_on_its_own_output_without_being_asked():
    rows = PARENT_ROWS + [
        row(n, s, 0.70 + 0.01 * i, contact=float(i), coulomb=2.0 * i)       # coulomb is exactly 2 * contact
        for i, n in enumerate(("a", "b", "c", "d", "e")) for s in (0, 1)
    ]
    ranked, _ = chem.rank_weighted(rows, 0.15)
    assert any("contact and coulomb" in line for line in chem.check_score_correlations(ranked, stream=io.StringIO()))
