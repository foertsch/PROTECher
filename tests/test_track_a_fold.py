"""examples/track_a_fold.py: input handling offline, and the dummy fold end to end.

The pure functions need nothing installed. The last two tests go through
`get_adapter` with the `dummy_boltz` fixture, so they need the core
dependencies like the rest of the suite.
"""

from __future__ import annotations

import csv
import importlib.util
from pathlib import Path

import pytest

SCRIPT = Path(__file__).resolve().parent.parent / "examples" / "track_a_fold.py"
_spec = importlib.util.spec_from_file_location("track_a_fold", SCRIPT)
fold = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(fold)

SCFV = "EVQLVESGGG" * 3
ANTIGEN = "ACDEFGHIKLMNPQRSTVWY" * 5  # 100 residues


def write_fasta(path: Path, *records: tuple[str, str]) -> Path:
    path.write_text("".join(f">{name} some description\n{seq}\n" for name, seq in records))
    return path


def test_read_sequences_from_fasta_and_from_a_bare_string(tmp_path: Path):
    fasta = write_fasta(tmp_path / "s.fa", ("parent", "EVQL"), ("v1", "evql\nkkk"))
    assert fold.read_sequences(str(fasta)) == [("parent", "EVQL"), ("v1", "EVQLKKK")]
    assert fold.read_sequences(" EVQL KKK ") == [("seq", "EVQLKKK")]


def test_a_sequence_longer_than_a_file_name_is_still_a_sequence():
    assert fold.read_sequences("A" * 5000) == [("seq", "A" * 5000)]


def test_non_standard_residue_is_named_with_its_position():
    with pytest.raises(ValueError, match=r"v1: not a standard amino acid: X at 3"):
        fold.clean("v1", "EVXQL")


def test_empty_sequence_and_headerless_fasta_are_rejected(tmp_path: Path):
    with pytest.raises(ValueError, match="empty sequence"):
        fold.clean("v1", "  ")
    bare = tmp_path / "bare.fa"
    bare.write_text("EVQL\n")
    with pytest.raises(ValueError, match="before the first '>' header"):
        fold.read_sequences(str(bare))


def test_regions_and_residue_lists_parse():
    assert fold.parse_region("480-607") == (480, 607)
    assert fold.parse_residues("557, 558,560-562") == [557, 558, 560, 561, 562]
    with pytest.raises(ValueError):
        fold.parse_region("607-480")
    with pytest.raises(ValueError):
        fold.parse_region("480")
    with pytest.raises(ValueError):
        fold.parse_residues("0,5")


def test_trim_by_region_is_one_based_and_inclusive():
    kept, start, end = fold.trim_antigen(ANTIGEN, region=(1, 20))
    assert (kept, start, end) == (ANTIGEN[:20], 1, 20)
    kept, start, end = fold.trim_antigen(ANTIGEN, region=(41, 60))
    assert kept == ANTIGEN[40:60] and (start, end) == (41, 60)


def test_trim_by_epitope_pads_and_clips_at_the_ends():
    kept, start, end = fold.trim_antigen(ANTIGEN, epitope=[50, 55], pad=10)
    assert (start, end) == (40, 65) and len(kept) == 26
    _, start, end = fold.trim_antigen(ANTIGEN, epitope=[3, 98], pad=10)
    assert (start, end) == (1, 100)


def test_no_trim_keeps_the_whole_antigen():
    assert fold.trim_antigen(ANTIGEN) == (ANTIGEN, 1, 100)


def test_trim_rejects_bad_requests():
    with pytest.raises(ValueError, match="not both"):
        fold.trim_antigen(ANTIGEN, region=(1, 10), epitope=[5])
    with pytest.raises(ValueError, match="runs past"):
        fold.trim_antigen(ANTIGEN, region=(90, 110))
    with pytest.raises(ValueError, match="past the antigen"):
        fold.trim_antigen(ANTIGEN, epitope=[120])


def test_complex_fasta_joins_the_chains_with_a_colon():
    assert fold.complex_fasta("v1", "AAA", "GGG") == ">v1\nAAA:GGG\n"


def test_warnings_for_a_strep_tag_and_a_long_antigen():
    tagged = [("v1", SCFV + fold.STREP_MOTIF)]
    assert any("Twin-Strep" in w for w in fold.warnings_for(tagged, 100))
    assert any("trim it" in w for w in fold.warnings_for([("v1", SCFV)], fold.LONG_ANTIGEN + 1))
    assert fold.warnings_for([("v1", SCFV)], 100) == []


def test_write_inputs_makes_one_two_chain_fasta_per_scfv(tmp_path: Path):
    written = fold.write_inputs(tmp_path, [("parent", SCFV), ("v 1", SCFV)], ANTIGEN, (1, 100))
    assert [name for name, _ in written] == ["parent", "v_1"]
    assert written[0][1].read_text() == f">parent\n{SCFV}:{ANTIGEN}\n"
    assert (tmp_path / "antigen_trimmed.fa").read_text() == f">antigen_1-100\n{ANTIGEN}\n"
    assert sorted(p.name for p in (tmp_path / "inputs").iterdir()) == ["parent.fa", "v_1.fa"]  # only complexes in inputs/


def test_the_inputs_are_the_files_the_batch_fold_on_modal_reads(tmp_path: Path):
    """`cloud/modal_fold.py::fold_dir` stages every `*.fa` in a folder through msa_cache.remote.read_chains."""
    from msa_cache.remote import read_chains

    written = fold.write_inputs(tmp_path, [("parent", SCFV), ("v1", SCFV)], ANTIGEN, (1, 100))
    assert [read_chains(path) for _, path in written] == [[SCFV, ANTIGEN], [SCFV, ANTIGEN]]


def test_two_scfv_records_with_one_name_are_rejected(tmp_path: Path):
    with pytest.raises(ValueError, match="unique name"):
        fold.write_inputs(tmp_path, [("a b", SCFV), ("a_b", SCFV)], ANTIGEN, (1, 100))


def test_dry_run_writes_inputs_and_folds_nothing(tmp_path: Path, monkeypatch, capsys):
    def boom(*args, **kwargs):
        raise AssertionError("a dry run must not fold")

    monkeypatch.setattr(fold, "run_folds", boom)
    scfv = write_fasta(tmp_path / "scfv.fa", ("parent", SCFV), ("v1", SCFV[:-1] + "K"))
    antigen = write_fasta(tmp_path / "ag.fa", ("antigen", ANTIGEN))
    out = tmp_path / "run"
    code = fold.main(
        ["--scfv", str(scfv), "--antigen", str(antigen), "--antigen-region", "21-60",
         "--out", str(out), "--dry-run"]
    )
    assert code == 0
    assert (out / "inputs" / "parent.fa").read_text() == f">parent\n{SCFV}:{ANTIGEN[20:60]}\n"
    assert not (out / "summary.csv").exists()
    assert "dry run" in capsys.readouterr().out


def test_bad_input_exits_with_a_message(tmp_path: Path, capsys):
    scfv = write_fasta(tmp_path / "scfv.fa", ("parent", "EVQLX"))
    with pytest.raises(SystemExit) as raised:
        fold.main(["--scfv", str(scfv), "--antigen", ANTIGEN, "--dry-run", "--out", str(tmp_path)])
    assert raised.value.code == 2
    assert "not a standard amino acid" in capsys.readouterr().err


def test_dummy_fold_writes_one_row_per_scfv_and_seed(tmp_path: Path):
    scfv = write_fasta(tmp_path / "scfv.fa", ("parent", SCFV), ("v1", SCFV[:-1] + "K"))
    out = tmp_path / "run"
    code = fold.main(
        ["--scfv", str(scfv), "--antigen", ANTIGEN, "--epitope", "50,55", "--pad", "10",
         "--seeds", "2", "--out", str(out), "--dummy"]
    )
    assert code == 0
    with (out / "summary.csv").open() as handle:
        rows = list(csv.DictReader(handle))
    assert len(rows) == 4
    assert {(r["complex"], r["seed"]) for r in rows} == {
        ("parent", "0"), ("parent", "1"), ("v1", "0"), ("v1", "1")
    }
    assert all(r["iptm"] and Path(r["pdb"]).is_file() for r in rows)


def test_seed_spread_needs_two_seeds_for_a_standard_deviation():
    rows = [
        {"complex": "parent", "iptm": 0.90},
        {"complex": "parent", "iptm": 0.94},
        {"complex": "v1", "iptm": 0.80},
    ]
    spread = {name: (n, sd) for name, n, _, sd in fold.seed_spread(rows)}
    assert spread["parent"][0] == 2 and spread["parent"][1] == pytest.approx(0.02828, abs=1e-4)
    assert spread["v1"] == (1, None)
