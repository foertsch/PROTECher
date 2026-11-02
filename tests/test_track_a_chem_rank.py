"""examples/track_a_chem_rank.py: the ranking rules on synthetic rows. Scoring real folds needs Boltz output and is not run here."""

from __future__ import annotations

import importlib.util
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
