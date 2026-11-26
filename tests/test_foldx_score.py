"""The FoldX scoring script: its notation and parsing, and the adapter wiring. FoldX itself is not needed."""

from __future__ import annotations

import csv
import importlib.util
import sys
from pathlib import Path

import pytest

SCRIPT = Path(__file__).resolve().parent.parent / "design_pipeline" / "adapters" / "_foldx_score.py"
_spec = importlib.util.spec_from_file_location("_foldx_score", SCRIPT)
foldx_score = importlib.util.module_from_spec(_spec)
sys.modules["_foldx_score"] = foldx_score
_spec.loader.exec_module(foldx_score)


def test_one_line_per_candidate_with_every_mutation_in_it():
    assert foldx_score.mutation_line([("A", 50, "S", "F")]) == "SA50F;"
    assert foldx_score.mutation_line([("A", 28, "D", "K"), ("B", 31, "D", "K")]) == "DA28K,DB31K;"


def test_candidates_group_by_variant_id(tmp_path):
    table = tmp_path / "candidates.csv"
    table.write_text("variant_id,chain,position,wt,mut\nA:S50F,A,50,S,F\nA:D28K+B:D31K,A,28,D,K\nA:D28K+B:D31K,B,31,D,K\n")
    found = foldx_score.read_candidates(table)
    assert list(found) == ["A:S50F", "A:D28K+B:D31K"]
    assert found["A:D28K+B:D31K"] == [("A", 28, "D", "K"), ("B", 31, "D", "K")]


def test_interaction_energy_is_read_from_the_last_data_line(tmp_path):
    path = tmp_path / "ac.fxout"
    path.write_text("junk\nPdb\tGroup1\tGroup2\tIntraclashes Group1\tInteraction Energy\nx.pdb\tAB\tC\t0\t-12.5\nx.pdb\tAB\tC\t0\t-13.25\n")
    assert foldx_score.parse_interaction_energy(path) == pytest.approx(-13.25)


def test_total_energies_one_per_model(tmp_path):
    path = tmp_path / "Dif.fxout"
    path.write_text("FoldX\n\nPdb\ttotal energy\tBackbone\nm_1_0.pdb\t1.5\t0\nm_1_1.pdb\t2.5\t0\n")
    assert foldx_score.parse_total_energies(path) == [1.5, 2.5]


def test_the_table_carries_the_spread_and_names_the_folding_term_for_what_it_is():
    assert foldx_score.FIELDS == ["variant_id", "ddg", "ddg_sd", "ddg_folding", "n_runs"]
    assert "ddg_fold" not in foldx_score.FIELDS   # the old name read as an aggregate of ddg; it is a separate quantity


def test_a_table_under_the_old_header_is_rewritten_and_keeps_every_row(tmp_path):
    """An append of today's five columns under a four-column header would silently misalign the file."""
    out_csv = tmp_path / "foldx.csv"
    out_csv.write_text("variant_id,ddg,ddg_fold,n_runs\nA:S50F,-1.5,2.25,3\nA:D28K,0.5,-1.0,3\n")

    scored, rewritten = foldx_score.read_scored(out_csv)

    assert rewritten and list(scored) == ["A:S50F", "A:D28K"]
    rows = list(csv.DictReader(out_csv.open(newline="")))
    assert [r["variant_id"] for r in rows] == ["A:S50F", "A:D28K"]
    assert rows[0]["ddg"] == "-1.5" and rows[0]["ddg_folding"] == "2.25" and rows[0]["n_runs"] == "3"
    # Not 0: an older run's per-run spread is gone, and 0 would claim its three runs agreed exactly.
    assert rows[0]["ddg_sd"] == "" and rows[1]["ddg_sd"] == ""


def test_a_table_already_on_the_current_header_is_left_alone(tmp_path):
    out_csv = tmp_path / "foldx.csv"
    text = "variant_id,ddg,ddg_sd,ddg_folding,n_runs\nA:S50F,-1.5,0.125,2.25,3\n"
    out_csv.write_text(text)

    scored, rewritten = foldx_score.read_scored(out_csv)

    assert not rewritten and list(scored) == ["A:S50F"]
    assert out_csv.read_text() == text


def test_no_table_yet_is_not_a_rewrite(tmp_path):
    assert foldx_score.read_scored(tmp_path / "foldx.csv") == ({}, False)


def test_the_foldx_adapter_is_wired_and_scores_per_variant():
    from design_pipeline.adapters import _scorers

    assert _scorers.foldx.tool is _scorers.foldx_process
    assert _scorers.foldx.per_variant and _scorers.foldx.sign == -1.0 and _scorers.foldx.value_column == "ddg"
    assert _scorers.foldx_process.output_globs == ("foldx.csv",)
