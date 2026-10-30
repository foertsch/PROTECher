"""site_guard: the thresholds, the failure modes and the command line, with a stub in place of the DockQ executable."""

import importlib
import json
import sys

import biotite.structure.io as strucio
import pytest

import site_guard


@pytest.fixture
def dockq(tmp_path, monkeypatch):
    """A stand-in DockQ that writes whatever the test puts in the STUB_RESULT environment variable."""
    stub = tmp_path / "DockQ"
    stub.write_text(
        f"#!{sys.executable}\n"
        "import json, os, sys\n"
        "if os.environ.get('STUB_FAIL'):\n"
        "    sys.exit(3)\n"
        "out = sys.argv[sys.argv.index('--json') + 1]\n"
        "json.dump({'best_result': {'AB': json.loads(os.environ['STUB_RESULT'])}}, open(out, 'w'))\n"
    )
    stub.chmod(0o755)
    monkeypatch.setattr(site_guard, "DOCKQ", str(stub))

    def returns(dockq_value, fnat, irmsd=1.0, lrmsd=2.0):
        monkeypatch.setenv("STUB_RESULT", json.dumps({"DockQ": dockq_value, "fnat": fnat, "iRMSD": irmsd, "LRMSD": lrmsd}))

    return returns


@pytest.fixture
def models(tmp_path):
    """Two model files and a reference as PDB paths, which check() passes to DockQ without converting."""
    for name in ("ref.pdb", "m_model_0.pdb"):
        (tmp_path / name).write_text("END\n")
    return tmp_path / "m_model_0.pdb", tmp_path / "ref.pdb"


def test_a_confident_fold_on_the_reference_pose_passes(dockq, models):
    dockq(0.89, 0.93)
    result = site_guard.check(*models, iptm=0.91)
    assert result["passed"] and result["reasons"] == []
    assert result["dockq_ref"] == 0.89 and result["fnat_ref"] == 0.93


def test_a_confident_fold_on_another_surface_fails_on_the_pose_alone(dockq, models):
    """The case that motivates the pose check: ipTM 0.80 on the wrong epitope, DockQ 0.01 against the reference."""
    dockq(0.01, 0.0, irmsd=24.0, lrmsd=65.0)
    result = site_guard.check(*models, iptm=0.80)
    assert not result["passed"]
    assert any("DockQ vs reference" in r for r in result["reasons"]) and any("fnat" in r for r in result["reasons"])
    assert not any("ipTM" in r for r in result["reasons"])


def test_the_right_pose_with_a_low_iptm_fails_on_iptm_alone(dockq, models):
    dockq(0.80, 0.90)
    result = site_guard.check(*models, iptm=0.59)
    assert not result["passed"] and len(result["reasons"]) == 1 and "ipTM 0.59" in result["reasons"][0]


def test_every_threshold_is_inclusive(dockq, models):
    dockq(site_guard.DOCKQ_MIN, site_guard.FNAT_MIN)
    assert site_guard.check(*models, iptm=site_guard.IPTM_MIN)["passed"]
    dockq(site_guard.DOCKQ_MIN - 0.01, site_guard.FNAT_MIN)
    assert not site_guard.check(*models, iptm=site_guard.IPTM_MIN)["passed"]


def test_a_comparison_that_cannot_be_made_counts_as_a_fail(dockq, models, monkeypatch):
    monkeypatch.setenv("STUB_FAIL", "1")
    result = site_guard.check(*models, iptm=0.95)
    assert not result["passed"] and "DockQ failed" in result["reasons"][0]
    monkeypatch.setattr(site_guard, "DOCKQ", "/nonexistent/DockQ")
    result = site_guard.check(*models, iptm=0.95)
    assert not result["passed"] and "DockQ not found" in result["reasons"][0]


def test_iptm_is_read_from_the_confidence_file_beside_the_model(tmp_path):
    (tmp_path / "confidence_x_model_0.json").write_text(json.dumps({"iptm": 0.42, "ptm": 0.9}))
    assert site_guard.iptm_of(tmp_path / "x_model_0.cif") == 0.42
    assert site_guard.iptm_of(tmp_path / "x_model_0.pdb") == 0.42


def atom(serial, name, res, chain, number, xyz, element, record="ATOM"):
    return (f"{record:<6}{serial:>5} {' ' + name.ljust(3) if len(element) == 1 else name:<4} {res:>3} {chain}{number:>4}    "
            f"{xyz[0]:>8.3f}{xyz[1]:>8.3f}{xyz[2]:>8.3f}{1.0:>6.2f}{0.0:>6.2f}          {element:>2}\n")


def test_conversion_keeps_protein_heavy_atoms_and_chains_and_drops_hydrogens_and_water(tmp_path):
    lines, serial = [], 1
    for number, chain in ((1, "A"), (2, "A"), (1, "B")):
        for k, (name, element) in enumerate([("N", "N"), ("CA", "C"), ("C", "C"), ("O", "O"), ("CB", "C"), ("H", "H")]):
            lines.append(atom(serial, name, "ALA", chain, number, (number + k * 1.4, 0.0, 0.0 if chain == "A" else 20.0), element))
            serial += 1
    lines.append(atom(serial, "O", "HOH", "A", 3, (5.0, 5.0, 5.0), "O", record="HETATM"))
    (tmp_path / "in.pdb").write_text("".join(lines) + "END\n")
    out = site_guard.to_pdb(tmp_path / "in.pdb", tmp_path / "out.pdb")
    atoms = strucio.load_structure(str(out))
    assert set(atoms.chain_id) == {"A", "B"}
    assert "H" not in set(atoms.element) and "HOH" not in set(atoms.res_name)
    assert len(atoms) == 3 * 5


def test_the_command_line_exits_1_if_any_model_fails_and_prints_one_line_per_model(dockq, models, tmp_path, capsys):
    model, reference = models
    (tmp_path / "confidence_m_model_0.json").write_text(json.dumps({"iptm": 0.9}))
    dockq(0.85, 0.9)
    assert site_guard.main(["site_guard", str(reference), str(model)]) == 0
    assert capsys.readouterr().out.startswith("PASS  m_model_0.pdb  ipTM 0.90")
    dockq(0.02, 0.0)
    assert site_guard.main(["site_guard", str(reference), str(model), str(model)]) == 1
    lines = capsys.readouterr().out.strip().splitlines()
    assert len(lines) == 2 and all(line.startswith("FAIL") for line in lines)
    assert site_guard.main(["site_guard"]) == 2


def test_the_default_dockq_path_comes_from_protecher_scratch_unless_dockq_is_set(monkeypatch):
    with monkeypatch.context() as m:
        m.delenv("DOCKQ", raising=False)
        m.setenv("PROTECHER_SCRATCH", "/x/scratch")
        importlib.reload(site_guard)
        assert site_guard.DOCKQ == "/x/scratch/protecher_tools/dockq/bin/DockQ"
        m.setenv("DOCKQ", "/y/DockQ")
        importlib.reload(site_guard)
        assert site_guard.DOCKQ == "/y/DockQ"
    importlib.reload(site_guard)
