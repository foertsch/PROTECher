"""hotspots: the FoldX alanine-scan positions of the trastuzumab scFv, checked against the parent and against scfv.py.

Every wild-type letter and every scFv number is looked up in the parent sequence, and the scFv numbers are recomputed with
the layout from `scfv`, which does not use the table.
"""

import pytest

from design_pipeline import hotspots, scfv

PARENT = (
    "EVQLVESGGGLVQPGGSLRLSCAASGFNIKDTYIHWVRQAPGKGLEWVARIYPTNGYTRYADSVKGRFTISADTSKNTAYLQMNSLRAEDTAVYYCSRWGGDGFYAMDYWGQGTLVTVSS"
    "GGGGSGGGGSGGGGS"
    "DIQMTQSPSSLSASVGDRVTITCRASQDVNTAVAWYQQKPGKAPKLLIYSASFLYSGVPSRFSGSRSGTDFTLTISSLQPEDFATYYCQQHYTTPPTFGQGTKVEIK"
)
VH, VL = PARENT.split(scfv.LINKER)
DOMAIN = {"VH": VH, "VL": VL}


def of(status):
    return [h for h in hotspots.HOTSPOTS if h.status == status]


def test_ten_hotspots_five_in_each_domain_all_losing_at_least_one_kcal():
    ten = of("hotspot")
    assert len(ten) == 10
    assert sorted(h.domain for h in ten) == ["VH"] * 5 + ["VL"] * 5
    assert all(h.ddg_bind >= 1.0 for h in ten)


def test_the_near_misses_and_the_excluded_position_stay_under_one_kcal():
    assert {(h.domain, h.position) for h in of("near_miss")} == {("VH", 59), ("VL", 92)}
    assert {(h.domain, h.position) for h in of("excluded")} == {("VL", 66)}
    assert all(h.ddg_bind < 1.0 for h in of("near_miss") + of("excluded"))


def test_the_alanine_improves_positions_carry_no_energy_value():
    assert {(h.domain, h.position) for h in of("alanine_improves")} == {("VL", 28), ("VL", 50), ("VH", 103)}
    assert all(h.ddg_bind is None for h in of("alanine_improves"))


def test_every_wild_type_letter_is_the_parents():
    for h in hotspots.HOTSPOTS:
        assert DOMAIN[h.domain][h.position - 1] == h.wt, f"{h.domain} {h.position}"
        assert PARENT[h.scfv_position - 1] == h.wt, f"scFv {h.scfv_position}"


def test_the_scfv_numbers_match_the_layout_and_the_heavy_chain_is_b_the_light_chain_a():
    layout = scfv.layout_from_parent(PARENT)
    for h in hotspots.HOTSPOTS:
        assert h.chain == {"VH": "B", "VL": "A"}[h.domain]
        assert scfv.scfv_position(layout, h.chain, h.position) == h.scfv_position


def test_no_position_is_listed_twice():
    keys = [(h.domain, h.position) for h in hotspots.HOTSPOTS]
    assert len(keys) == len(set(keys)) == 16


def test_positions_returns_the_scfv_numbers_of_one_status():
    assert hotspots.positions() == frozenset({33, 50, 57, 99, 105, 165, 166, 188, 226, 229})
    assert hotspots.positions("near_miss") == frozenset({59, 227})
    with pytest.raises(ValueError, match="status"):
        hotspots.positions("hotspots")


def test_a_crystal_style_id_prefix_names_the_chain_and_wild_type():
    by_position = {(h.domain, h.position): h for h in hotspots.HOTSPOTS}
    assert by_position[("VH", 33)].crystal_id == "B:Y33"
    assert by_position[("VL", 30)].crystal_id == "A:N30"
