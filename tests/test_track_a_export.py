"""examples/track_a_export.py: designs in crystal numbering to scFv sequences and a submission file, on the trastuzumab scFv.

Expected sequences are built without the mapping (edit VH and VL as strings, rejoin), so the script is not tested
against its own offset.
"""

from __future__ import annotations

import csv
import importlib.util
from pathlib import Path

import pytest

from design_pipeline import scfv

SCRIPT = Path(__file__).resolve().parent.parent / "examples" / "track_a_export.py"
_spec = importlib.util.spec_from_file_location("track_a_export", SCRIPT)
export = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(export)

VH = "EVQLVESGGGLVQPGGSLRLSCAASGFNIKDTYIHWVRQAPGKGLEWVARIYPTNGYTRYADSVKGRFTISADTSKNTAYLQMNSLRAEDTAVYYCSRWGGDGFYAMDYWGQGTLVTVSS"
VL = "DIQMTQSPSSLSASVGDRVTITCRASQDVNTAVAWYQQKPGKAPKLLIYSASFLYSGVPSRFSGSRSGTDFTLTISSLQPEDFATYYCQQHYTTPPTFGQGTKVEIK"
PARENT = VH + scfv.LINKER + VL


def mut(chain: str, position: int, new: str | None = None) -> str:
    """A variant-id token whose wild type is read from the domain strings, not from the layout."""
    wt = (VH if chain == "B" else VL)[position - 1]
    new = new or ("A" if wt != "A" else "G")
    return f"{chain}:{wt}{position}{new}"


def design(*tokens: str) -> str:
    return "+".join(tokens)


FIVE = design(mut("B", 97), mut("B", 103, "Y"), mut("B", 104, "Y"), mut("A", 27, "E"), mut("A", 91, "R"))
OTHER_FIVE = design(mut("B", 31), mut("B", 50), mut("B", 55), mut("A", 30), mut("A", 53))


@pytest.fixture
def layout():
    return scfv.layout_from_parent(PARENT)


def by_hand(*edits: tuple[str, int, str]) -> str:
    """The expected scFv: change VH and VL on their own, then join them with the linker."""
    vh, vl = list(VH), list(VL)
    for chain, position, new in edits:
        (vh if chain == "B" else vl)[position - 1] = new
    return "".join(vh) + scfv.LINKER + "".join(vl)


def test_read_parent_from_fasta_a_name_sequence_line_a_header_and_a_bare_sequence(tmp_path):
    fasta = tmp_path / "p.fa"
    fasta.write_text(f">trastuzumab_scFv_parent some text\n{PARENT[:100]}\n{PARENT[100:]}\n")
    line = tmp_path / "p.csv"
    line.write_text(f"trastuzumab_scFv_parent,{PARENT}\n")
    headed = tmp_path / "h.csv"
    headed.write_text(f"name,sequence\ntrastuzumab_scFv_parent,{PARENT}\n")
    assert export.read_parent(str(fasta)) == ("trastuzumab_scFv_parent", PARENT)
    assert export.read_parent(str(line)) == ("trastuzumab_scFv_parent", PARENT)
    assert export.read_parent(str(headed)) == ("trastuzumab_scFv_parent", PARENT)
    assert export.read_parent(f"trastuzumab_scFv_parent,{PARENT}") == ("trastuzumab_scFv_parent", PARENT)
    assert export.read_parent(" " + PARENT.lower() + " ") == ("parent", PARENT)


def test_read_parent_names_a_residue_that_is_not_an_amino_acid():
    with pytest.raises(ValueError, match="not a standard amino acid: X at 3"):
        export.read_parent("EVXQ")
    with pytest.raises(ValueError, match="empty"):
        export.read_parent("name,   ")


def test_read_designs_orders_by_the_design_column_and_keeps_each_id_once(tmp_path):
    table = tmp_path / "designs.csv"
    rows = [(2, FIVE, "B"), (2, FIVE, "A"), (1, OTHER_FIVE, "B"), (1, OTHER_FIVE, "A")]
    table.write_text("design,variant_id,chain\n" + "".join(f"{d},{v},{c}\n" for d, v, c in rows))
    assert export.read_designs(str(table)) == [OTHER_FIVE, FIVE]


def test_read_designs_reads_a_parquet_and_refuses_a_table_without_variant_id(tmp_path):
    pd = pytest.importorskip("pandas")
    pytest.importorskip("pyarrow")
    path = tmp_path / "designs.parquet"
    pd.DataFrame({"variant_id": [FIVE, FIVE, OTHER_FIVE], "design": [1, 1, 2]}).to_parquet(path)
    assert export.read_designs(str(path)) == [FIVE, OTHER_FIVE]
    pd.DataFrame({"id": [FIVE]}).to_csv(tmp_path / "bad.csv", index=False)
    with pytest.raises(ValueError, match="no variant_id column"):
        export.read_designs(str(tmp_path / "bad.csv"))


def test_a_design_becomes_the_scfv_with_its_substitutions_at_the_right_places(layout):
    designs = export.build(layout, [FIVE, OTHER_FIVE])
    assert [d.name for d in designs] == ["design_01", "design_02"]
    assert designs[0].sequence == by_hand(("B", 97, "A"), ("B", 103, "Y"), ("B", 104, "Y"), ("A", 27, "E"), ("A", 91, "R"))
    assert [(p, w, m) for p, w, m in designs[0].substitutions] == [(97, "S", "A"), (103, "G", "Y"), (104, "F", "Y"), (162, "Q", "E"), (226, "H", "R")]
    assert all(len(d.sequence) == 242 and len(d.substitutions) == 5 for d in designs)


def test_fewer_than_five_substitutions_is_named_and_the_threshold_can_be_changed(layout):
    four = design(mut("B", 97), mut("B", 103, "Y"), mut("A", 27, "E"), mut("A", 91, "R"))
    with pytest.raises(ValueError, match=r"design_01 .*: 4 substitutions, the rule is at least 5"):
        export.build(layout, [four])
    assert len(export.build(layout, [four], min_substitutions=4)) == 1
    with pytest.raises(ValueError, match="more than the 5 allowed"):
        export.build(layout, [design(FIVE, mut("B", 31))], max_substitutions=5)


def test_a_constant_domain_position_wrong_wild_type_and_repeats_are_all_reported_at_once(layout):
    in_ch1 = design(mut("B", 97), mut("B", 103, "Y"), mut("B", 104, "Y"), mut("A", 27, "E"), "B:K213A")
    wrong_wt = design(mut("B", 97), mut("B", 103, "Y"), mut("B", 104, "Y"), mut("A", 27, "E"), "A:W91R")
    with pytest.raises(ValueError) as caught:
        export.build(layout, [FIVE, in_ch1, wrong_wt, FIVE])
    message = str(caught.value)
    assert message.startswith("3 problem(s)")
    assert "design_02" in message and "outside the variable domain" in message
    assert "design_03" in message and "does not match the parent" in message
    assert "design_04" in message and "the same sequence as design_01" in message


def test_a_sequence_carrying_the_strep_tag_is_refused(layout, monkeypatch):
    monkeypatch.setattr(export, "STREP_MOTIF", "QQHYTTPP")  # a stretch of this parent's light chain
    with pytest.raises(ValueError, match="the parent carries QQHYTTPP"):
        export.build(layout, [FIVE])


def run(tmp_path, *args, parent=PARENT):
    parent_file = tmp_path / "parent.csv"
    parent_file.write_text(f"trastuzumab_scFv_parent,{parent}\n")
    return export.main(["--parent", str(parent_file), "--out", str(tmp_path / "out"), *args])


def test_the_command_line_writes_the_three_files(tmp_path, capsys):
    assert run(tmp_path, "--ids", FIVE, OTHER_FIVE) == 0
    out = tmp_path / "out"
    rows = list(csv.reader((out / "submission.csv").open()))
    assert rows[0] == ["design_01", by_hand(("B", 97, "A"), ("B", 103, "Y"), ("B", 104, "Y"), ("A", 27, "E"), ("A", 91, "R"))]
    assert len(rows) == 2 and rows[1][0] == "design_02"
    fasta = (out / "scfv_for_folding.fasta").read_text().splitlines()
    assert fasta[0] == ">trastuzumab_scFv_parent" and fasta[1] == PARENT and fasta[2] == ">design_01"
    report = list(csv.DictReader((out / "designs_report.csv").open()))
    assert report[0]["n_substitutions"] == "5" and report[0]["n_in_chain_A"] == "2" and report[0]["n_in_chain_B"] == "3"
    assert report[0]["scfv_substitutions"] == "S97A;G103Y;F104Y;Q162E;H226R"
    text = capsys.readouterr().out
    assert "NOT checked against a crystal" in text and "2 designs, 5 to 5 substitutions each" in text


def test_header_flag_and_a_design_table_as_input(tmp_path):
    table = tmp_path / "designs.csv"
    table.write_text(f"design,variant_id\n1,{FIVE}\n2,{OTHER_FIVE}\n")
    assert run(tmp_path, "--designs", str(table), "--header", "--name-prefix", "t") == 0
    rows = list(csv.reader((tmp_path / "out" / "submission.csv").open()))
    assert rows[0] == ["name", "sequence"] and rows[1][0] == "t_01" and rows[2][0] == "t_02"


def test_nothing_is_written_when_a_design_is_refused(tmp_path, capsys):
    assert run(tmp_path, "--ids", FIVE, "B:K213A") == 1
    assert not (tmp_path / "out").exists()
    assert "not written: 1 problem(s)" in capsys.readouterr().err


def test_the_structure_check_runs_when_a_crystal_is_given_and_says_where_it_differs(tmp_path, monkeypatch, capsys):
    crystal = {"B": VH + "ASTKGP", "A": VL + "RTVAAP"}
    monkeypatch.setattr(scfv, "chain_sequence", lambda structure, chain: crystal[chain])
    assert run(tmp_path, "--ids", FIVE, "--structure", "1N8Z.cif") == 0
    assert "matches the variable domains of 1N8Z.cif" in capsys.readouterr().out
    crystal["A"] = VL[:49] + ("G" if VL[49] != "G" else "A") + VL[50:] + "RTVAAP"
    assert run(tmp_path, "--ids", FIVE, "--structure", "1N8Z.cif") == 1
    assert "first difference at position 50" in capsys.readouterr().err


def test_a_light_chain_first_scfv(tmp_path):
    flipped = VL + scfv.LINKER + VH
    assert run(tmp_path, "--ids", design(*[mut("B", p) for p in (31, 50, 55)], mut("A", 30), mut("A", 53)), "--light-first", parent=flipped) == 0
    sequence = list(csv.reader((tmp_path / "out" / "submission.csv").open()))[0][1]
    vl, vh = list(VL), list(VH)
    for p in (31, 50, 55):
        vh[p - 1] = "A" if vh[p - 1] != "A" else "G"
    for p in (30, 53):
        vl[p - 1] = "A" if vl[p - 1] != "A" else "G"
    assert sequence == "".join(vl) + scfv.LINKER + "".join(vh)


def test_ids_and_a_design_table_together_are_refused(tmp_path):
    with pytest.raises(SystemExit):
        run(tmp_path, "--ids", FIVE, "--designs", "x.csv")
