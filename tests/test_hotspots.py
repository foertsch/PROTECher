"""hotspots: the FoldX alanine-scan positions of the trastuzumab scFv, checked against the parent and against scfv.py.

Every wild-type letter and every scFv number is looked up in the parent sequence, and the scFv numbers are recomputed with
the layout from `scfv`, which does not use the table.
"""

import csv
from pathlib import Path

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


# The per-position table: the raw 3-run scan plus the values quoted in docs/track_a_filtering.md.
TABLE = Path(__file__).resolve().parents[1] / "trastuzumab_foldx_alanine_scan.csv"
DOMAIN_OF_CHAIN = {"B": "VH", "A": "VL"}


def table_rows():
    with TABLE.open(newline="", encoding="utf-8") as handle:
        return list(csv.DictReader(handle))


def test_the_table_has_every_position_of_the_two_domains():
    rows = table_rows()
    assert len(rows) == 227
    assert sum(r["chain"] == "B" for r in rows) == 120 and sum(r["chain"] == "A" for r in rows) == 107
    assert sum(r["status"] == "ok" for r in rows) == 210
    assert sum(r["status"] == "wt_is_ala" for r in rows) == 17


def test_every_row_matches_the_parent_and_the_scfv_numbering():
    layout = scfv.layout_from_parent(PARENT)
    for r in table_rows():
        domain = DOMAIN_OF_CHAIN[r["chain"]]
        assert r["domain"] == domain
        assert DOMAIN[domain][int(r["resseq"]) - 1] == r["wt"], f"{domain} {r['resseq']}"
        assert int(r["scfv_position"]) == scfv.scfv_position(layout, r["chain"], int(r["resseq"]))


def test_the_raw_columns_are_the_untouched_three_run_scan():
    """VH R50 is +0.34 here and +1.98 at 10 runs. Both are kept, in different columns, so nobody 'fixes' the raw one."""
    by = {(r["chain"], int(r["resseq"])): r for r in table_rows()}
    raw = by[("B", 50)]
    assert float(raw["ddg_bind"]) == pytest.approx(0.34, abs=0.01) and raw["n_runs"] == "3"
    assert float(raw["ddg_bind_doc"]) == pytest.approx(1.98) and raw["runs_doc"] == "10"
    assert {r["n_runs"] for r in table_rows() if r["status"] == "ok"} == {"3"}


def test_the_doc_columns_agree_with_the_hotspots_module():
    by = {(r["domain"], int(r["resseq"])): r for r in table_rows()}
    with_value = {(h.domain, h.position): h for h in hotspots.HOTSPOTS if h.ddg_bind is not None}
    filled = {key for key, r in by.items() if r["ddg_bind_doc"] != ""}
    assert filled == set(with_value)
    for key, h in with_value.items():
        assert float(by[key]["ddg_bind_doc"]) == pytest.approx(h.ddg_bind)
        assert by[key]["runs_doc"] == ("" if h.runs is None else str(h.runs))


def test_the_doc_values_at_or_above_one_are_exactly_the_ten_hotspots():
    above = {int(r["scfv_position"]) for r in table_rows() if r["ddg_bind_doc"] != "" and float(r["ddg_bind_doc"]) >= 1.0}
    assert above == hotspots.positions()
