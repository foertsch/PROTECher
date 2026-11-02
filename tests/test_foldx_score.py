"""The FoldX scoring script: its notation and parsing, and the adapter wiring. FoldX itself is not needed."""

from __future__ import annotations

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


def test_the_foldx_adapter_is_wired_and_scores_per_variant():
    from design_pipeline.adapters import _scorers

    assert _scorers.foldx.tool is _scorers.foldx_process
    assert _scorers.foldx.per_variant and _scorers.foldx.sign == -1.0 and _scorers.foldx.value_column == "ddg"
    assert _scorers.foldx_process.output_globs == ("foldx.csv",)
