"""Regression checks for local candidate screening; no models, GPU or network."""

import csv
import importlib.util
import json
from pathlib import Path

import pytest

SCRIPT = Path(__file__).resolve().parents[1] / "hackday/05_developability.py"
spec = importlib.util.spec_from_file_location("developability", SCRIPT)
dev = importlib.util.module_from_spec(spec)
spec.loader.exec_module(dev)


def mutate(changes):
    seq = list(dev.PARENT)
    for pos, aa in changes.items():
        seq[pos - 1] = aa
    return "".join(seq)


def candidate(seq, **kwargs):
    return dev.screen({"variant": seq}, **kwargs)[1]


def test_official_reference_and_zero_deltas():
    assert len(dev.PARENT) == 242
    assert dev.PARENT[120:135] == "GGGGSGGGGSGGGGS"
    assert [i for i, aa in enumerate(dev.PARENT, 1) if aa == "C"] == [22, 96, 158, 223]
    row = dev.screen({}, lp={})[0]
    assert row["assessment"] == "reference"
    assert row["n_mut"] == 0
    for key in ("d_pI", "d_charge", "d_gravy", "d_hydrophobic_run", "antifold_d"):
        assert row[key] == 0
    if row["instability"] is not None:
        assert row["d_instability"] == 0


def test_instability_matches_published_biopython_example():
    pytest.importorskip("Bio")
    # Expected 41.98 in Biopython's official ProtParam documentation.
    seq = ("MAEGEITTFTALTEKFNLPPGNYKKPKLLYCSNGGHFLRILPDGTVDGT"
           "RDRSDQHIQLQLSAESVGEVYIKSTETGQYLAMDTSGLLYGSQTPSEEC"
           "LFLERLEENHYNTYTSKKHAEKNWFVGLKKNGSCKRGPRTHYGQKAILF"
           "LPLPV")
    assert dev.instability_index(seq) == pytest.approx(41.98, abs=0.005)


@pytest.mark.parametrize("count,status", [(0, "fail"), (4, "fail"), (5, "pass")])
def test_minimum_five_substitutions(count, status):
    row = candidate(mutate({i: "A" for i in range(1, count + 1)}))
    assert row["n_mut"] == count
    assert row["rules_status"] == status
    assert row["ipsae"] is None


def test_cysteine_swap_with_unchanged_count_is_detected():
    seq = mutate({22: "A", 21: "C"})
    assert seq.count("C") == dev.PARENT.count("C")
    row = candidate(seq)
    assert "conserved_cysteine_changed" in row["rule_issues"]
    assert "introduced_cysteine" in row["flags"]
    assert row["introduced_cysteines"] == "21"


def test_linker_mutation_is_a_rule_violation_and_not_a_partial_antifold_score():
    row = candidate(mutate({121: "A"}), lp={1: {"E": -1, "A": -2}})
    assert "linker_changed" in row["rule_issues"]
    assert row["antifold_d"] is None


def test_new_glycosylation_motif_is_reported_as_a_warning():
    row = candidate(mutate({30: "S"}))
    assert "glyc:N28" in row["new_motifs"]
    assert "new_liability" in row["flags"]
    assert "new_glycosylation_motif" in row["flags"]
    assert not dev.screen({})[0]["new_motifs"]


@pytest.mark.parametrize("sequence", ["", "ACDX", "ACD*", "ACD-"])
def test_invalid_sequences_are_retained_without_scores(sequence):
    row = candidate(sequence)
    assert row["assessment"] == "invalid_input"
    assert row["rules_status"] == "not_evaluated"
    assert row["pI"] is None
    assert row["antifold_d"] is None


def test_length_change_does_not_get_truncated_mutation_count():
    row = candidate(dev.PARENT[:-1], lp={})
    assert row["rules_status"] == "review_required"
    assert row["n_mut"] is None
    assert row["linker_unchanged"] is None
    assert row["antifold_d"] is None
    assert row["PSH"] is None
    assert row["pI"] is not None


def test_missing_or_partial_antifold_is_not_zero():
    seq = mutate({1: "A", 2: "A"})
    assert dev.antifold_delta(seq, dev.PARENT, None) is None
    lp = {1: {"E": -1, "A": -3}}
    assert dev.antifold_delta(seq, dev.PARENT, lp) is None
    lp[2] = {"V": -2, "A": -5}
    assert dev.antifold_delta(seq, dev.PARENT, lp) == -5
    lp[2]["A"] = 0
    assert dev.antifold_delta(seq, dev.PARENT, lp) == 0  # real zero is retained


def test_missing_biopython_leaves_explicit_missing_value(monkeypatch):
    def missing(_seq):
        raise ImportError("not installed")
    monkeypatch.setattr(dev, "instability_index", missing)
    row = candidate(dev.PARENT)
    assert row["instability"] is None
    assert "instability_unavailable" in row["data_warnings"]


@pytest.mark.parametrize("suffix,contents", [
    (".fasta", ">a\nACD\n>a\nEFG\n"),
    (".csv", "name,sequence\na,ACD\na,EFG\n"),
    (".fasta", "ACD\n>a\nEFG\n"),
    (".fasta", ">\nACD\n"),
    (".csv", "name,sequence\na\n"),
    (".csv", "design,sequence\na,ACD\n"),
])
def test_ambiguous_input_is_rejected(tmp_path, suffix, contents):
    path = tmp_path / ("input" + suffix)
    path.write_text(contents)
    with pytest.raises(ValueError):
        dev.read_sequences(path)


def test_candidate_only_csv_and_wrapped_fasta(tmp_path):
    path = tmp_path / "input.csv"
    path.write_text("name,sequence\na,acd efg\n")
    assert dev.read_sequences(path) == {"a": "ACDEFG"}
    path = tmp_path / "input.fasta"
    path.write_text(">a\nacd\nefg\n")
    assert dev.read_sequences(path) == {"a": "ACDEFG"}
    assert list(dev.prepare_sequences(dev.read_sequences(path), "PARENT")) == ["PARENT", "a"]


def test_incorrect_parent_cannot_silently_replace_official_reference():
    with pytest.raises(ValueError, match="official parent"):
        dev.prepare_sequences({"PARENT": mutate({1: "A"})}, "PARENT")
    assert dev.prepare_sequences({"trastuzumab_scFv_parent": dev.PARENT}, "PARENT") == {
        "PARENT": dev.PARENT
    }


def test_antifold_reference_and_duplicate_validation(tmp_path):
    path = tmp_path / "af.csv"
    row = {"pdb_chain": "B", "pdb_pos": 1, "pdb_res": "E", **dict.fromkeys(dev.AAS, -1)}
    def write(rows):
        with path.open("w", newline="") as fh:
            writer = csv.DictWriter(fh, fieldnames=list(row))
            writer.writeheader()
            writer.writerows(rows)
    write([row])
    assert dev.read_antifold(path)[1]["E"] == -1
    for rows in ([row, row], [{**row, "pdb_res": "A"}], [{**row, "A": "nan"}]):
        write(rows)
        with pytest.raises(ValueError):
            dev.read_antifold(path)


@pytest.mark.parametrize("problem", [None, "mismatch", "missing", "nonfinite"])
def test_structure_mapping_is_checked_before_scoring(monkeypatch, problem):
    import numpy as np
    import biotite.structure as struc
    # Synthetic CA-only geometry isolates coverage/identity validation from SASA.
    positions = list(range(1, 121)) + list(range(136, 243))
    if problem == "missing":
        positions.pop()
    atoms = struc.AtomArray(len(positions))
    atoms.coord[:] = 0
    atoms.atom_name[:] = "CA"
    reverse = {one: three for three, one in dev.THREE_TO_ONE.items()}
    atoms.res_name = np.array([reverse[dev.PARENT[i-1]] for i in positions])
    if problem == "mismatch":
        atoms.res_name[0] = "ALA"
    if problem == "nonfinite":
        atoms.coord[0] = np.nan
    monkeypatch.setattr(struc, "sasa", lambda *_args, **_kwargs: np.ones(len(positions)))
    if problem:
        with pytest.raises(ValueError):
            dev.surface_table(atoms, np.array(positions))
    else:
        surf = dev.surface_table(atoms, np.array(positions))
        row = dev.screen({}, surf=surf, lp={})[0]
        for metric in ("d_PSH", "d_PPC", "d_PNC", "d_SFvCSP"):
            assert row[metric] == 0


@pytest.mark.parametrize("design,sequence,value", [
    ("a", "ACD", "0.7"), ("unknown", "ACE", "0.7"),
    ("a", "ACE", "nan"), ("a", "ACE", "1.1"),
])
def test_ipsae_rejects_stale_unknown_and_invalid_scores(tmp_path, design, sequence, value):
    path = tmp_path / "scores.csv"
    path.write_text(f"design,sequence,ipsae\n{design},{sequence},{value}\n")
    with pytest.raises(ValueError):
        dev.read_ipsae(path, {"a": "ACE"})


def test_sequence_only_cli_and_ipsae_join(tmp_path, monkeypatch):
    sequence = mutate({i: "A" for i in range(1, 6)})
    source, scores, out = (tmp_path / name for name in ("seq.csv", "ip.csv", "report.csv"))
    source.write_text(f"name,sequence\nv1,{sequence}\n")
    scores.write_text(f"design,sequence,ipsae\nv1,{sequence},0.75\n")
    def forbidden(*_args):
        pytest.fail("Sequence-only mode must not load structures or AntiFold")
    monkeypatch.setattr(dev, "load_fv", forbidden)
    monkeypatch.setattr(dev, "read_antifold", forbidden)
    assert dev.main(["--input", str(source), "--ipsae", str(scores),
                     "--sequence-only", "--out", str(out)]) == 0
    with out.open() as fh:
        rows = list(csv.DictReader(fh))
    assert len(rows) == 2
    assert rows[1]["rules_status"] == "pass"
    assert rows[1]["ipsae"] == "0.75"
    assert rows[1]["antifold_d"] == ""
    assert rows[1]["PSH"] == ""
    meta = json.loads(out.with_suffix(".metadata.json").read_text())
    assert meta["input"]["sha256"]
    assert meta["structure"] is None


def test_missing_resources_do_not_prevent_sequence_report(tmp_path):
    source, out = tmp_path / "input.fasta", tmp_path / "report.csv"
    source.write_text(f">v1\n{dev.PARENT}\n")
    dev.main(["--input", str(source), "--structure", str(tmp_path / "missing.pdb"),
              "--antifold", str(tmp_path / "missing.csv"), "--out", str(out)])
    with out.open() as fh:
        row = list(csv.DictReader(fh))[1]
    assert row["pI"]
    assert row["antifold_d"] == ""
    assert "structure_unavailable" in row["data_warnings"]


def test_output_cannot_overwrite_input(tmp_path):
    source = tmp_path / "input.csv"
    text = f"name,sequence\nv1,{dev.PARENT}\n"
    source.write_text(text)
    with pytest.raises(SystemExit):
        dev.main(["--input", str(source), "--out", str(source)])
    assert source.read_text() == text


@pytest.mark.parametrize("tol", [-1, float("nan"), float("inf")])
def test_invalid_patch_tolerance_is_rejected(tol):
    with pytest.raises(ValueError):
        dev.screen({}, tol=tol)


@pytest.mark.parametrize("contents", ["name,sequence\nv1,ACD\n", "v1,ACD,extra\n", "v1\n"])
def test_headerless_csv_refuses_headers_and_wrong_column_counts(tmp_path, contents):
    path = tmp_path / "input.csv"
    path.write_text(contents)
    with pytest.raises(ValueError):
        dev.read_sequences(path, csv_no_header=True)


def test_headerless_mode_cannot_be_used_for_fasta(tmp_path):
    path = tmp_path / "input.fasta"
    path.write_text(">v1\nACD\n")
    with pytest.raises(ValueError, match=".csv"):
        dev.read_sequences(path, csv_no_header=True)


def test_default_csv_mode_does_not_silently_drop_first_candidate(tmp_path):
    path = tmp_path / "input.csv"
    path.write_text("v1,ACD\nv2,ACE\n")
    with pytest.raises(ValueError):
        dev.read_sequences(path)
    assert dev.read_sequences(path, csv_no_header=True) == {"v1": "ACD", "v2": "ACE"}


def test_duplicate_sequences_keep_both_ids_and_explain_redundant_runs():
    seq = mutate({i: "A" for i in range(1, 6)})
    rows = dev.screen({"first": seq, "second": seq, "unique": mutate({30: "S"})})
    assert len(rows) == 4
    for row, other in zip(rows[1:3], ("second", "first")):
        assert row["rules_status"] == "pass"
        assert "duplicate_sequence" in row["flags"]
        assert row["duplicate_sequence_ids"] == other
        assert "intentional controls or seed replicates" in row["review_notes"]
    assert not rows[0]["duplicate_sequence_ids"]
    assert "duplicate_sequence" not in rows[3]["flags"]


@pytest.mark.parametrize("mutation,motif", [("S", "glyc:N28"), ("D", "isomer:DD@30")])
def test_notes_identify_mutation_creating_new_motif(mutation, motif):
    row = candidate(mutate({30: mutation}))
    assert motif in row["liability_details"]
    assert f"K30{mutation}" in row["liability_details"]
    assert row["liability_details"] in row["review_notes"]


def test_motif_details_include_each_changed_residue_in_the_motif():
    seq = mutate({2: "N", 3: "A", 4: "S"})
    row = candidate(seq)
    assert "glyc:N2 overlaps substitutions V2N, Q3A, L4S" in row["liability_details"]


@pytest.mark.parametrize("format", ["fasta", "csv_header", "csv_no_header"])
def test_team_exporter_outputs_feed_screening_without_manual_edits(tmp_path, format):
    from design_pipeline import scfv
    spec = importlib.util.spec_from_file_location("track_a_export", SCRIPT.parents[1] / "examples/track_a_export.py")
    export = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(export)
    # Independent full-scFv positions, crossing both domains.
    seq = mutate({30: "S", 66: "R", 97: "A", 159: "K", 202: "R"})
    designs = export.build(scfv.layout_from_parent(dev.PARENT),
                           ["B:K30S+B:G66R+B:S97A+A:R24K+A:S67R"])
    assert designs[0].sequence == seq
    path, out = tmp_path / ("input.fasta" if format == "fasta" else "input.csv"), tmp_path / "report.csv"
    args = ["--input", str(path), "--out", str(out), "--sequence-only"]
    if format == "fasta":
        export.write_fasta(path, ("parent", dev.PARENT), designs)
    else:
        export.write_submission(path, designs, header=format == "csv_header")
        if format == "csv_no_header":
            args.append("--csv-no-header")
    assert dev.main(args) == 0
    with out.open() as fh:
        rows = list(csv.DictReader(fh))
    assert len(rows) == 2  # no duplicate reference from exporter's lowercase parent
    assert rows[1]["design"] == designs[0].name
    assert rows[1]["sequence"] == seq
    assert rows[1]["rules_status"] == "pass"
    assert rows[1]["n_mut"] == "5"
    assert "R159K" in rows[1]["mutations"]
    assert "new_glycosylation_motif" in rows[1]["flags"]
