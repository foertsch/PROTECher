"""ipSAE on small PAE matrices whose answers can be worked out by hand."""

from __future__ import annotations

import numpy as np
import pytest

from design_pipeline.scoring import METRICS, _d0, ipsae


def two_chains(n_a: int, n_b: int, inter: float) -> tuple[np.ndarray, np.ndarray]:
    """A PAE matrix with a small error inside each chain and `inter` between them."""
    chains = np.array(["A"] * n_a + ["B"] * n_b)
    pae = np.where(chains[:, None] == chains[None, :], 1.0, inter)
    return pae, chains


def test_d0_floor_and_formula():
    assert _d0(2) == pytest.approx(1.0)                     # n is taken as at least 26, which gives 0.958, floored to 1
    assert _d0(26) == pytest.approx(1.0)
    assert _d0(100) == pytest.approx(1.24 * 85 ** (1 / 3) - 1.8)
    assert _d0(370) == pytest.approx(1.24 * 355 ** (1 / 3) - 1.8)


def test_zero_error_between_chains_scores_one():
    pae, chains = two_chains(40, 30, 0.0)
    scores = ipsae(pae, chains)
    assert scores == {"ab": 1.0, "ba": 1.0, "max": 1.0, "min": 1.0}


def test_error_above_the_cutoff_scores_zero():
    pae, chains = two_chains(40, 30, 30.0)
    assert ipsae(pae, chains) == {"ab": 0.0, "ba": 0.0, "max": 0.0, "min": 0.0}


def test_a_pae_equal_to_the_cutoff_is_not_counted():
    pae, chains = two_chains(5, 5, 15.0)
    assert ipsae(pae, chains, pae_cutoff=15.0)["max"] == 0.0
    assert ipsae(pae, chains, pae_cutoff=15.5)["max"] > 0.0


def test_hand_worked_two_by_two_complex():
    # chain A = residues 0, 1; chain B = residues 2, 3. Rows align, columns are scored.
    pae = np.full((4, 4), 20.0)
    pae[0, 2], pae[0, 3] = 1.0, 2.0                          # A0 -> B: both pass, n = 2 so d0 = 1
    pae[2, 0], pae[2, 1] = 3.0, 3.0                          # B0 -> A: both pass
    scores = ipsae(pae, ["A", "A", "B", "B"])
    assert scores["ab"] == pytest.approx((1 / (1 + 1) + 1 / (1 + 4)) / 2)       # 0.35
    assert scores["ba"] == pytest.approx(1 / (1 + 9))                           # 0.1
    assert scores["max"] == pytest.approx(0.35) and scores["min"] == pytest.approx(0.1)


def test_the_best_aligned_residue_decides_and_a_bad_one_does_not_drag_it_down():
    pae = np.full((4, 4), 20.0)
    pae[0, 2] = pae[0, 3] = 0.0                              # A0 sees B perfectly
    pae[1, 2] = pae[1, 3] = 14.0                             # A1 sees it badly but under the cutoff
    assert ipsae(pae, ["A", "A", "B", "B"])["ab"] == pytest.approx(1.0)


def test_d0_grows_with_the_number_of_passing_residues():
    few = np.full((40, 40), 20.0)
    many = np.full((40, 40), 20.0)
    few[0, 20:23] = 5.0                                      # 3 passing partners: d0 = 1
    many[0, 20:40] = 5.0                                     # 20 passing partners: n is below 26 too, d0 = 1
    chains = ["A"] * 20 + ["B"] * 20
    assert ipsae(many, chains)["ab"] == pytest.approx(ipsae(few, chains)["ab"])
    wide = np.full((80, 80), 20.0)
    wide[0, 40:80] = 5.0                                     # 40 passing partners: d0 > 1, the same PAE scores higher
    assert ipsae(wide, ["A"] * 40 + ["B"] * 40)["ab"] > ipsae(few, chains)["ab"]
    n = 40
    expected = 1 / (1 + (5.0 / (1.24 * (n - 15) ** (1 / 3) - 1.8)) ** 2)
    assert ipsae(wide, ["A"] * 40 + ["B"] * 40)["ab"] == pytest.approx(expected)


def test_the_two_directions_are_not_the_same_for_an_asymmetric_pae():
    pae, chains = two_chains(30, 30, 20.0)
    pae[0, 30:] = 2.0                                        # only A rows see B well
    scores = ipsae(pae, chains)
    assert scores["ab"] > 0.0 and scores["ba"] == 0.0
    assert scores["max"] == scores["ab"] and scores["min"] == 0.0


def test_chn_and_dom_variants_use_their_own_d0():
    pae, chains = two_chains(60, 60, 8.0)
    res, chn, dom = (ipsae(pae, chains, d0=v)["ab"] for v in ("res", "chn", "dom"))
    d0_res = 1.24 * (60 - 15) ** (1 / 3) - 1.8               # 60 passing partners per residue
    d0_chn = 1.24 * (120 - 15) ** (1 / 3) - 1.8              # both chains
    assert res == pytest.approx(1 / (1 + (8.0 / d0_res) ** 2))
    assert chn == pytest.approx(1 / (1 + (8.0 / d0_chn) ** 2))
    assert dom == pytest.approx(chn)                         # every residue passes, so the domain is both chains


def test_bad_input_is_refused():
    pae, chains = two_chains(5, 5, 1.0)
    with pytest.raises(ValueError):
        ipsae(pae, ["A"] * 10)                               # one chain only
    with pytest.raises(ValueError):
        ipsae(pae, ["A"] * 3 + ["B"] * 3 + ["C"] * 4)        # three chains
    with pytest.raises(ValueError):
        ipsae(pae[:9, :9], chains)                           # shape and residue count differ
    with pytest.raises(ValueError):
        ipsae(pae, chains, d0="nope")


def test_the_four_columns_are_registered_and_none_without_a_pae(tmp_path):
    from design_pipeline import Design

    assert {"ipsae_ab", "ipsae_ba", "ipsae_max", "ipsae_min"} <= set(METRICS)
    design = Design(id="x", path=tmp_path / "x.pdb", parent_id=None, producer="test")
    assert METRICS["ipsae_max"](design, None) is None
