"""examples/track_a_screen.py: reading a Boltz batch's files and scoring them. Folding itself needs a GPU and is not run."""

from __future__ import annotations

import csv
import importlib.util
import json
import sys
from pathlib import Path

import numpy as np
import pytest

SCRIPT = Path(__file__).resolve().parent.parent / "examples" / "track_a_screen.py"
_spec = importlib.util.spec_from_file_location("track_a_screen", SCRIPT)
screen = importlib.util.module_from_spec(_spec)
sys.modules["track_a_screen"] = screen
_spec.loader.exec_module(screen)


def write_cif(path: Path, n_a: int, n_b: int) -> None:
    import biotite.structure as struc
    from biotite.structure.io import pdbx

    atoms = struc.AtomArray(n_a + n_b)
    atoms.coord = np.arange(3 * (n_a + n_b), dtype=float).reshape(-1, 3)
    atoms.chain_id[:] = ["A"] * n_a + ["B"] * n_b
    atoms.res_id[:] = list(range(1, n_a + 1)) + list(range(1, n_b + 1))
    atoms.res_name[:] = "ALA"
    atoms.atom_name[:] = "CA"
    atoms.element[:] = "C"
    cif = pdbx.CIFFile()
    pdbx.set_structure(cif, atoms)
    cif.write(str(path))


def write_fold(root: Path, seed: int, name: str, inter_pae: float, iptm: float, plddt: float, n_a=30, n_b=20) -> None:
    folder = root / f"seed{seed}" / "boltz_results_inputs" / "predictions" / name
    folder.mkdir(parents=True)
    chains = np.array([0] * n_a + [1] * n_b)
    pae = np.where(chains[:, None] == chains[None, :], 1.0, inter_pae)
    np.savez(folder / f"pae_{name}_model_0.npz", pae=pae)
    write_cif(folder / f"{name}_model_0.cif", n_a, n_b)
    (folder / f"confidence_{name}_model_0.json").write_text(json.dumps({"iptm": iptm, "complex_plddt": plddt}))


def test_fasta_records_come_back_in_file_order(tmp_path):
    path = tmp_path / "c.fasta"
    path.write_text(">o01\nAAAA\n>k001\nCC\nCC\n")
    assert screen.read_fasta(path) == [("o01", "AAAA"), ("k001", "CCCC")]


def test_folds_are_scored_from_the_files_boltz_wrote(tmp_path):
    write_fold(tmp_path, 0, "good", inter_pae=0.0, iptm=0.9, plddt=0.88)
    write_fold(tmp_path, 1, "bad", inter_pae=30.0, iptm=0.3, plddt=0.7)
    rows = {r["name"]: r for r in screen.score_folds(tmp_path)}
    assert rows["good"]["ipsae_max"] == pytest.approx(1.0) and rows["good"]["iptm"] == 0.9
    assert rows["good"]["plddt"] == pytest.approx(88.0) and rows["good"]["seed"] == 0
    assert rows["bad"]["ipsae_max"] == 0.0 and rows["bad"]["seed"] == 1


def test_a_fold_with_a_missing_file_is_skipped_so_a_stopped_batch_still_scores(tmp_path):
    write_fold(tmp_path, 0, "done", inter_pae=0.0, iptm=0.9, plddt=0.9)
    write_fold(tmp_path, 0, "half", inter_pae=0.0, iptm=0.9, plddt=0.9)
    (tmp_path / "seed0" / "boltz_results_inputs" / "predictions" / "half" / "confidence_half_model_0.json").unlink()
    assert [r["name"] for r in screen.score_folds(tmp_path)] == ["done"]


def test_the_report_writes_one_row_per_fold(tmp_path):
    write_fold(tmp_path, 0, "a", inter_pae=0.0, iptm=0.9, plddt=0.9)
    screen.write_and_report(screen.score_folds(tmp_path), tmp_path)
    with (tmp_path / "screen_scores.csv").open() as handle:
        assert [r["name"] for r in csv.DictReader(handle)] == ["a"]
