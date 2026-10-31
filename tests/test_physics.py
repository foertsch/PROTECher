"""CPU screening must distinguish missing evidence from successful checks."""
import copy
import csv
import importlib.util
import json
from pathlib import Path

import numpy as np
import pytest
import biotite.structure as struc
from biotite.structure.io.pdb import PDBFile

from design_pipeline import physics as p

ROOT = Path(__file__).resolve().parents[1]
spec = importlib.util.spec_from_file_location("physics_cli", ROOT / "hackday/06_physics_screen.py")
cli = importlib.util.module_from_spec(spec)
spec.loader.exec_module(cli)


def atoms(rows):
    """Rows: chain, residue number/name, atom name, element, xyz."""
    a = struc.AtomArray(len(rows))
    for field, values in zip(("chain_id", "res_id", "res_name", "atom_name", "element", "coord"), zip(*rows)):
        setattr(a, field, np.array(values))
    return a


def analyze(a, binder=None, positions=None):
    return p.analyze_geometry(a, np.ones(len(a), bool) if binder is None else np.array(binder),
                              np.zeros(len(a), int) if positions is None else np.array(positions))


def record():
    return copy.deepcopy(json.loads((ROOT / "hackday/physics_parent.json").read_text())["structures"][0])


@pytest.fixture
def fv_atoms():
    # Intentionally incomplete CA-only coordinates: validates mapping, never
    # serves as a claim of realistic geometry or a biological control.
    inverse = {aa: name for name, aa in p.THREE_TO_ONE.items()}
    rows = []
    for pos, aa in enumerate(cli.dev.PARENT, 1):
        if 121 <= pos <= 135:
            continue
        chain, number = ("B", pos) if pos <= 120 else ("A", pos-135)
        rows.append((chain, number, inverse[aa], "CA", "C", [pos*4., 0., 0.]))
    return atoms(rows)


@pytest.fixture
def mapping(fv_atoms):
    r = record()
    r["antigen_segments"] = []
    return r


def test_exact_binder_mapping_and_explicit_fv_scope(fv_atoms, mapping):
    a, binder, pos, coverage = p.select_structure(fv_atoms, mapping, cli.dev.PARENT)
    assert len(a) == 227 and binder.all()
    assert set(pos) == set(range(1, 121)) | set(range(136, 243))
    assert coverage["scope"] == "fv_only" and coverage["antigen_signature"] is None


@pytest.mark.parametrize("change", ["missing", "wrong_sequence", "duplicate", "nan", "insertion"])
def test_invalid_structure_mapping_is_rejected(fv_atoms, mapping, change):
    a = fv_atoms.copy()
    if change == "missing": a = a[1:]
    if change == "wrong_sequence": a.res_name[0] = "ALA"
    if change == "duplicate": a = a + a[:1]
    if change == "nan": a.coord[0, 0] = np.nan
    if change == "insertion": a.ins_code[0] = "A"
    with pytest.raises(ValueError):
        p.select_structure(a, mapping, cli.dev.PARENT)


def test_full_scfv_scope_does_not_accept_missing_linker(fv_atoms, mapping):
    mapping["scope"] = "full_scfv"
    with pytest.raises(ValueError, match="cover exactly"):
        p.select_structure(fv_atoms, mapping, cli.dev.PARENT)


def test_full_scfv_explicit_mapping_includes_linker():
    inverse = {aa: name for name, aa in p.THREE_TO_ONE.items()}
    a = atoms([("X", pos, inverse[aa], "CA", "C", [pos*4., 0, 0])
               for pos, aa in enumerate(cli.dev.PARENT, 1)])
    r = {"scope": "full_scfv", "binder_segments": [{"chain": "X", "start": 1, "end": 242, "sequence_start": 1}]}
    _, binder, pos, coverage = p.select_structure(a, r, cli.dev.PARENT)
    assert binder.all() and set(pos) == set(range(1, 243))
    assert coverage["mapped_binder_residues"] == 242


@pytest.mark.parametrize("suffix", [".pdb", ".cif", ".bcif"])
def test_structure_formats_keep_author_chain_ids(tmp_path, fv_atoms, suffix):
    from biotite.structure.io import pdbx
    path = tmp_path / ("structure"+suffix)
    if suffix == ".pdb":
        f = PDBFile(); f.set_structure(fv_atoms)
    else:
        f = pdbx.CIFFile() if suffix == ".cif" else pdbx.BinaryCIFFile()
        pdbx.set_structure(f, fv_atoms)
    f.write(path)
    loaded = p.load_structure(path)
    assert np.array_equal(loaded.chain_id, fv_atoms.chain_id)
    assert np.array_equal(loaded.res_id, fv_atoms.res_id)


def test_antigen_gap_changes_signature(fv_atoms, mapping):
    target = atoms([("C", 480, "ALA", "CA", "C", [0, 0, 0])])
    mapping["antigen_segments"] = [{"chain": "C", "start": 480, "end": 481}]
    *_, coverage = p.select_structure(fv_atoms+target, mapping, cli.dev.PARENT)
    assert coverage["antigen_missing_residues"] == ["C:481"]
    assert coverage["antigen_signature"] == p.sequence_hash("A?")


def test_covalent_neighbours_not_clashes_but_distinct_chains_are():
    a = atoms([("A", 1, "ALA", "CA", "C", [0, 0, 0]),
               ("A", 1, "ALA", "CB", "C", [1.5, 0, 0])])
    assert analyze(a)["metrics"]["clashes_binder"] == 0
    a.chain_id[1] = "B"
    result = analyze(a, [True, False])
    assert result["metrics"]["clashes_interface"] == 1


def test_covalent_peptide_distance_is_checked_separately():
    a = atoms([("A", 1, "ALA", "C", "C", [0, 0, 0]),
               ("A", 2, "ALA", "N", "N", [2.5, 0, 0])])
    result = analyze(a)
    assert result["metrics"]["peptide_outliers_binder"] == 1
    assert result["metrics"]["clashes_binder"] == 0


def test_salt_contacts_are_per_residue_pair_and_charge_specific():
    a = atoms([("A", 1, "ASP", "OD1", "O", [0, 0, 0]),
               ("A", 1, "ASP", "OD2", "O", [0, .3, 0]),
               ("B", 1, "LYS", "NZ", "N", [3, 0, 0]),
               ("B", 2, "LYS", "NZ", "N", [6, 0, 0]),
               ("B", 3, "HIS", "NE2", "N", [0, 1, 0])])
    result = analyze(a)
    assert result["metrics"]["salt_bridges_binder"] == 1
    assert result["metrics"]["like_charge_contacts_binder"] == 1
    assert result["histidine_residues_excluded_from_charge_analysis"] == 1


def test_hydrogen_bonds_missing_not_zero_and_polar_contacts_not_bonds():
    a = atoms([("A", 1, "ASN", "ND2", "N", [0, 0, 0]),
               ("B", 1, "ASP", "OD1", "O", [2.8, 0, 0])])
    result = analyze(a)
    assert result["metrics"]["hydrogen_bonds_binder"] is None
    assert result["metrics"]["polar_contacts_binder"] == 1


@pytest.mark.parametrize("hydrogen,expected", [([1, 0, 0], 1), ([0, 1, 0], 0)])
def test_hydrogen_bonds_require_angle_and_correct_chemistry(hydrogen, expected):
    a = atoms([("A", 1, "ASN", "ND2", "N", [0, 0, 0]),
               ("A", 1, "ASN", "HD21", "H", hydrogen),
               ("B", 1, "ASP", "OD1", "O", [2.8, 0, 0])])
    result = analyze(a)
    assert result["metrics"]["hydrogen_bonds_binder"] == expected
    a.res_name[2], a.atom_name[2], a.element[2] = "LYS", "NZ", "N"
    assert analyze(a)["metrics"]["hydrogen_bonds_binder"] == 0


def test_partial_hydrogens_are_reported_as_partial():
    a = atoms([("A", 1, "ASN", "ND2", "N", [0, 0, 0]),
               ("A", 1, "ASN", "HD21", "H", [1, 0, 0]),
               ("B", 1, "LYS", "NZ", "N", [10, 0, 0])])
    assert analyze(a)["hydrogen_bond_status"] == "partial_explicit_hydrogens"


def test_expected_disulfide_geometry_and_stretched_bond():
    a = atoms([("A", 1, "CYS", "CB", "C", [0, 1.8, 0]),
               ("A", 1, "CYS", "SG", "S", [0, 0, 0]),
               ("A", 8, "CYS", "SG", "S", [2.05, 0, 0]),
               ("A", 8, "CYS", "CB", "C", [2.05, 0, 1.8])])
    result = analyze(a, positions=[22, 22, 96, 96])
    assert result["metrics"]["disulfide_geometry_outliers"] == 0
    assert result["metrics"]["clashes_binder"] == 0
    assert result["metrics"]["disulfides_unassessed"] == 1  # VL absent
    a.coord[2:] += [1, 0, 0]
    assert analyze(a, positions=[22, 22, 96, 96])["metrics"]["disulfide_geometry_outliers"] == 1


def test_surface_area_and_no_target_are_explicit():
    a = atoms([("A", 1, "ALA", "CB", "C", [0, 0, 0])])
    result = analyze(a)
    assert result["metrics"]["binder_sasa_A2"] > 0
    assert result["metrics"]["binder_hydrophobic_sidechain_sasa_A2"] > 0
    assert result["metrics"]["binder_buried_on_binding_A2"] is None
    assert result["metrics"]["clashes_interface"] is None
    assert "incomplete_heavy_atoms" in result["missing_checks"]


@pytest.fixture
def campaign(tmp_path, fv_atoms, mapping):
    pdb = tmp_path / "parent.pdb"
    f = PDBFile(); f.set_structure(fv_atoms); f.write(pdb)
    mapping["path"] = "parent.pdb"
    manifest = tmp_path / "manifest.json"
    manifest.write_text(json.dumps({"schema_version": 1, "structures": [mapping]}))
    seq = "A"+cli.dev.PARENT[1:]
    source = tmp_path / "candidates.csv"
    source.write_text(f"name,sequence\nvariant,{seq}\n")
    return source, manifest, tmp_path / "out", mapping


def test_end_to_end_missing_candidates_never_inherit_parent_scores(campaign):
    source, manifest, out, _ = campaign
    rows, meta = cli.run(source, manifest, out)
    assert rows[0]["physics_status"] == "partial"
    candidate = rows[1]
    assert candidate["physics_status"] == "not_evaluated_no_structure"
    assert all(candidate[m] is None for m in p.METRICS)
    assert meta["inputs"][0]["sha256"]
    assert "unavailable" in (out / "physics_report.md").read_text()
    dev_rows = cli.dev.screen({"variant": candidate["sequence"]})
    cli.dev.attach_physics(dev_rows, out / "physics_review.csv")
    assert dev_rows[1]["physics_clashes_binder"] is None
    assert dev_rows[1]["physics_status"] == "not_evaluated_no_structure"


def test_reused_parent_coordinates_for_changed_sequence_fail(campaign):
    source, manifest, out, mapping = campaign
    candidate = copy.deepcopy(mapping)
    candidate.update(design="variant", sequence_sha256=p.sequence_hash("A"+cli.dev.PARENT[1:]))
    manifest.write_text(json.dumps({"schema_version": 1, "structures": [mapping, candidate]}))
    rows, _ = cli.run(source, manifest, out)
    assert rows[1]["physics_status"] == "invalid_structure"
    assert "mismatch" in rows[1]["physics_error"]


def test_manifest_rejects_stale_sequence_hash(campaign):
    source, manifest, out, mapping = campaign
    mapping["sequence_sha256"] = "stale"
    manifest.write_text(json.dumps({"schema_version": 1, "structures": [mapping]}))
    with pytest.raises(ValueError, match="hash mismatch"):
        cli.run(source, manifest, out)


def test_import_rejects_stale_sequences(campaign):
    source, manifest, out, _ = campaign
    cli.run(source, manifest, out)
    rows = cli.dev.screen({"variant": "C"+cli.dev.PARENT[1:]})
    with pytest.raises(ValueError, match="mismatch"):
        cli.dev.attach_physics(rows, out / "physics_review.csv")


@pytest.mark.parametrize("change", ["unknown_id", "duplicate", "inf"])
def test_physics_import_rejects_invalid_records(campaign, change):
    source, manifest, out, _ = campaign
    cli.run(source, manifest, out)
    path = out / "physics_review.csv"
    with path.open() as handle:
        reader = csv.DictReader(handle)
        rows = list(reader); fields = reader.fieldnames
    if change == "unknown_id": rows[1]["design"] = "unrelated"
    if change == "duplicate": rows.append(rows[1].copy())
    if change == "inf": rows[0]["clashes_binder"] = "inf"
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields)
        writer.writeheader(); writer.writerows(rows)
    with pytest.raises(ValueError):
        cli.dev.attach_physics(cli.dev.screen({"variant": "A"+cli.dev.PARENT[1:]}), path)


def test_cli_candidate_structure_errors_have_nonzero_exit(campaign):
    source, manifest, out, mapping = campaign
    mapping["path"] = "absent.pdb"
    manifest.write_text(json.dumps({"schema_version": 1, "structures": [mapping]}))
    assert cli.main(["--input", str(source), "--manifest", str(manifest), "--out-dir", str(out)]) == 1
    row = next(csv.DictReader((out / "physics_review.csv").open()))
    assert row["physics_status"] == "invalid_structure" and not row["clashes_binder"]


def test_cannot_overwrite_manifest(campaign):
    source, manifest, out, mapping = campaign
    path = out / "physics_metadata.json"
    out.mkdir(); path.write_text(manifest.read_text())
    # Paths in the new manifest remain valid if absolute.
    mapping["path"] = str(manifest.parent / "parent.pdb")
    path.write_text(json.dumps({"schema_version": 1, "structures": [mapping]}))
    with pytest.raises(ValueError, match="overwrite"):
        cli.run(source, path, out)


def test_local_parent_smoke(tmp_path):
    path = ROOT / "data/structures/1N8Z.pdb"
    if not path.exists():
        pytest.skip("Optional local reference structure; no download during tests")
    rows, meta = cli.run(None, ROOT / "hackday/physics_parent.json", tmp_path)
    row = rows[0]
    assert row["physics_status"] != "invalid_structure"
    assert row["disulfides_unassessed"] == 0
    assert row["hydrogen_bonds_binder"] is None
    assert meta["structure_details"]["PARENT"]["mapped_binder_residues"] == 227
