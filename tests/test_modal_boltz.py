"""examples/modal_boltz.py: the helpers and the dry run, offline.

Modal is not a dependency of this repo, so the tests install a stand-in `modal` module whose
decorators return the function unchanged. What is tested is the logic around the cloud call, not
the cloud call itself.
"""

from __future__ import annotations

import importlib.util
import json
import sys
import types
from pathlib import Path

import pytest

SCRIPT = Path(__file__).resolve().parent.parent / "examples" / "modal_boltz.py"


class _Chain:
    def __getattr__(self, name):
        return lambda *args, **kwargs: _Chain()


def _fake_modal() -> types.ModuleType:
    module = types.ModuleType("modal")

    class App:
        def __init__(self, *args, **kwargs):
            pass

        def function(self, *args, **kwargs):
            return lambda fn: fn

        def local_entrypoint(self, *args, **kwargs):
            return lambda fn: fn

    module.App = App
    module.Volume = types.SimpleNamespace(from_name=lambda *a, **k: _Chain())
    module.Image = types.SimpleNamespace(debian_slim=lambda *a, **k: _Chain())
    return module


@pytest.fixture
def mb(monkeypatch):
    monkeypatch.setitem(sys.modules, "modal", _fake_modal())
    spec = importlib.util.spec_from_file_location("modal_boltz_under_test", SCRIPT)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


SCFV = "EVQLVESGGG" * 3
ANTIGEN = "ACDEFGHIKLMNPQRSTVWY" * 5


def test_boltz_fasta_labels_chains_a_then_b(mb):
    assert mb.boltz_fasta(["AAA", "GGG"]) == ">A|protein\nAAA\n>B|protein\nGGG\n"


def test_confidence_row_reports_plddt_on_the_0_to_100_scale(mb):
    text = json.dumps({"complex_plddt": 0.921, "ptm": 0.94, "iptm": 0.93})
    row = mb.confidence_row("parent", 2, text)
    assert row == {"complex": "parent", "seed": 2, "iptm": 0.93, "ptm": 0.94, "plddt": pytest.approx(92.1)}


def test_confidence_row_tolerates_missing_fields(mb):
    row = mb.confidence_row("v1", 0, json.dumps({"iptm": 0.5}))
    assert row["plddt"] is None and row["ptm"] is None and row["iptm"] == 0.5


def test_dry_run_counts_folds_and_starts_nothing(mb, tmp_path: Path, capsys):
    scfv = tmp_path / "scfv.fa"
    scfv.write_text(f">parent\n{SCFV}\n>v1\n{SCFV[:-1]}K\n")
    mb.fold = types.SimpleNamespace(starmap=lambda jobs: pytest.fail("a dry run must not fold"))
    mb.main(scfv=str(scfv), antigen=ANTIGEN, antigen_region="21-60", seeds=3,
            out=str(tmp_path / "out"), dry_run=True)
    printed = capsys.readouterr().out
    assert "2 complex(es) x 3 seed(s) = 6 fold(s)" in printed
    assert "antigen residues 21-60 of 100 (40 kept)" in printed
    assert "dry run" in printed
    assert not (tmp_path / "out").exists()


def test_main_writes_files_and_a_summary_from_the_fold_results(mb, tmp_path: Path, capsys):
    scfv = tmp_path / "scfv.fa"
    scfv.write_text(f">parent\n{SCFV}\n")

    def starmap(jobs):
        for name, chains, seed in jobs:
            assert chains == [SCFV, ANTIGEN] and name == "parent"
            yield {
                "name": name, "seed": seed, "cif": "data_x\n",
                "confidence": json.dumps({"complex_plddt": 0.9, "ptm": 0.9, "iptm": 0.9 + 0.01 * seed}),
                "pae": b"npz",
            }

    mb.fold = types.SimpleNamespace(starmap=starmap)
    mb.main(scfv=str(scfv), antigen=ANTIGEN, seeds=2, out=str(tmp_path / "out"))
    assert (tmp_path / "out" / "predictions" / "parent_s1" / "model_0.cif").read_text() == "data_x\n"
    assert (tmp_path / "out" / "predictions" / "parent_s0" / "pae.npz").read_bytes() == b"npz"
    lines = (tmp_path / "out" / "summary.csv").read_text().splitlines()
    assert lines[0] == "complex,seed,iptm,ptm,plddt,cif" and len(lines) == 3
    assert "parent: ipTM mean 0.905" in capsys.readouterr().out


def test_bad_seed_count_exits(mb, tmp_path: Path):
    with pytest.raises(SystemExit, match="at least 1"):
        mb.main(scfv=SCFV, antigen=ANTIGEN, seeds=0, out=str(tmp_path), dry_run=True)
