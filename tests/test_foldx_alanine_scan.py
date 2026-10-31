"""examples/foldx_alanine_scan.py: the parsers, and the orchestration against a fake FoldX.

FoldX itself is licensed per user and not part of the test suite. The fake writes the files the real
one writes (shapes taken from FoldX 5.1 output), with energies chosen so the expected answer is known.
"""

from __future__ import annotations

import csv
import importlib.util
import stat
import sys
from pathlib import Path

import pytest

SCRIPT = Path(__file__).resolve().parent.parent / "examples" / "foldx_alanine_scan.py"
_spec = importlib.util.spec_from_file_location("foldx_alanine_scan", SCRIPT)
scan = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(scan)

PDB = (
    "ATOM      1  N   TYR A   5      11.0  12.0  13.0  1.00 20.00           N\n"
    "ATOM      2  CA  TYR A   5      12.0  12.0  13.0  1.00 20.00           C\n"
    "ATOM      3  CA  ALA A   6      13.0  12.0  13.0  1.00 20.00           C\n"
    "ATOM      4  CA  GLY B   9      14.0  12.0  13.0  1.00 20.00           C\n"
    "ATOM      5  CA  LYS C  30      15.0  12.0  13.0  1.00 20.00           C\n"
    "END\n"
)

FAKE = '''#!{python}
import pathlib, shutil, sys
args = dict(a[2:].split("=", 1) for a in sys.argv[1:] if a.startswith("--") and "=" in a)
if args["command"] == "BuildModel":
    name = args["pdb"]; stem = name[:-4]; out = pathlib.Path(args["output-dir"])
    runs = int(args["numberOfRuns"])
    assert pathlib.Path(args["mutant-file"]).read_text().strip().endswith("A;")
    for i in range(runs):
        shutil.copy(pathlib.Path(args["pdb-dir"]) / name, out / f"{{stem}}_1_{{i}}.pdb")
        shutil.copy(pathlib.Path(args["pdb-dir"]) / name, out / f"WT_{{stem}}_1_{{i}}.pdb")
    rows = "".join(f"{{stem}}_1_{{i}}.pdb\\t{{1.0 + i}}\\t0\\n" for i in range(runs))
    (out / f"Dif_{{stem}}.fxout").write_text(
        "FoldX\\n\\nOutput type: BuildModel\\nPdb\\ttotal energy\\tBackbone Hbond\\n" + rows)
else:
    name = args["pdb"]; stem = name[:-4]; out = pathlib.Path(args["output-dir"])
    energy = -19.0 if name.startswith("WT_") else -17.5
    (out / f"Interaction_{{stem}}_AC.fxout").write_text(
        "FoldX\\n\\nOutput type: AnalyseComplex\\n"
        "Pdb\\tGroup1\\tGroup2\\tIntraclashesGroup1\\tIntraclashesGroup2\\tInteraction Energy\\n"
        f"{{name}}\\tAB\\tC\\t1\\t2\\t{{energy}}\\n")
'''


@pytest.fixture
def fake_foldx(tmp_path: Path) -> Path:
    path = tmp_path / "foldx_fake"
    path.write_text(FAKE.format(python=sys.executable))
    path.chmod(path.stat().st_mode | stat.S_IXUSR)
    return path


@pytest.fixture
def repaired(tmp_path: Path) -> Path:
    path = tmp_path / "struct_Repair.pdb"
    path.write_text(PDB)
    return path


def test_read_positions_takes_ca_atoms_of_the_chosen_chains(repaired: Path):
    assert scan.read_positions(repaired, {"A", "B"}) == [("A", 5, "Y"), ("A", 6, "A"), ("B", 9, "G")]
    assert scan.read_positions(repaired, {"C"}) == [("C", 30, "K")]


def test_mutation_line_uses_foldx_notation():
    assert scan.mutation_line("Y", "B", 57) == "YB57A;"


def test_interaction_energy_is_read_by_column_name(tmp_path: Path):
    f = tmp_path / "i.fxout"
    f.write_text(
        "Output type: AnalyseComplex\n"
        "Pdb\tGroup1\tGroup2\tIntraclashesGroup1\tIntraclashesGroup2\tInteraction Energy\tBackbone Hbond\n"
        "x.pdb\tAB\tC\t25.4\t12.3\t-17.9384\t-1.4\n"
    )
    assert scan.parse_interaction_energy(f) == pytest.approx(-17.9384)


def test_total_energies_are_read_per_model_row(tmp_path: Path):
    f = tmp_path / "d.fxout"
    f.write_text("Output type: BuildModel\nPdb\ttotal energy\tBackbone Hbond\na_1.pdb\t1.5\t0\na_2.pdb\t2.5\t0\n")
    assert scan.parse_total_energies(f) == [1.5, 2.5]


def test_scan_writes_one_row_per_position_with_the_expected_shift(
    tmp_path: Path, fake_foldx: Path, repaired: Path
):
    out = tmp_path / "out"
    code = scan.main(
        ["--foldx", str(fake_foldx), "--pdb", str(repaired), "--scan-chains", "A,B", "--groups", "AB,C",
         "--runs", "3", "--workers", "2", "--out", str(out)]
    )
    assert code == 0
    with (out / "ddg_alanine_scan.csv").open() as handle:
        rows = {(r["chain"], r["resseq"]): r for r in csv.DictReader(handle)}
    assert set(rows) == {("A", "5"), ("A", "6"), ("B", "9")}
    tyr = rows[("A", "5")]
    assert tyr["status"] == "ok" and tyr["n_runs"] == "3"
    assert float(tyr["ddg_bind"]) == pytest.approx(1.5)  # -17.5 mutant against -19.0 wild type
    assert float(tyr["ddg_fold"]) == pytest.approx(2.0)  # mean of 1.0, 2.0, 3.0
    assert rows[("A", "6")]["status"] == "wt_is_ala"  # nothing to mutate
    assert not list(out.glob("work_*"))  # scratch is cleaned


def test_a_failing_foldx_marks_the_position_and_exits_nonzero(tmp_path: Path, repaired: Path):
    broken = tmp_path / "foldx_broken"
    broken.write_text(f"#!{sys.executable}\nprint('no output written')\n")
    broken.chmod(broken.stat().st_mode | stat.S_IXUSR)
    out = tmp_path / "out"
    code = scan.main(
        ["--foldx", str(broken), "--pdb", str(repaired), "--scan-chains", "B", "--groups", "AB,C",
         "--runs", "1", "--out", str(out)]
    )
    assert code == 1
    with (out / "ddg_alanine_scan.csv").open() as handle:
        (row,) = list(csv.DictReader(handle))
    assert row["status"].startswith("failed: foldx BuildModel did not write")


def test_limit_scans_only_the_first_positions(tmp_path: Path, fake_foldx: Path, repaired: Path):
    out = tmp_path / "out"
    scan.main(
        ["--foldx", str(fake_foldx), "--pdb", str(repaired), "--scan-chains", "A,B", "--groups", "AB,C",
         "--runs", "1", "--limit", "1", "--out", str(out)]
    )
    with (out / "ddg_alanine_scan.csv").open() as handle:
        assert len(list(csv.DictReader(handle))) == 1
