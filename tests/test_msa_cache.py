"""msa_cache: importing the parent's alignments, swapping the query row, and patching the Boltz adapter."""

import csv
import sys
from pathlib import Path

import pytest

from design_pipeline.adapters import _boltz
from design_pipeline.designs import Design
from msa_cache import CacheError, MsaCache, fetch_alignments, use_cached_msas
from msa_cache.__main__ import main
from msa_cache.cache import _boltz_python

SCFV = "EVQLVESGGGLVQPGGSLRLSCAASGFNIKDTYIHWVRQAPGKGLEWVARIYPTNGYTRYADSV"
ANTIGEN = "DECVGEGLACHQLCARGHCWGPGPTQCVNCSQFLRGQECVEECRVLQ"


def write_msa(directory: Path, name: str, index: int, query: str, others: list[str]) -> Path:
    directory.mkdir(parents=True, exist_ok=True)
    path = directory / f"{name}_{index}.csv"
    rows = [["key", "sequence"], ["0", query]] + [[str(i + 1), s] for i, s in enumerate(others)]
    with open(path, "w", newline="") as fh:
        csv.writer(fh).writerows(rows)
    return path


@pytest.fixture
def boltz_msas(tmp_path):
    msa_dir = tmp_path / "boltz_results_wt" / "msa"
    write_msa(msa_dir, "wt", 0, SCFV, ["-" * 10 + SCFV[10:], SCFV[:20] + "-" * (len(SCFV) - 20)])
    write_msa(msa_dir, "wt", 1, ANTIGEN, ["N" * len(ANTIGEN)])
    return msa_dir


@pytest.fixture
def cache(tmp_path, boltz_msas):
    c = MsaCache(tmp_path / "cache")
    c.import_boltz_msas(boltz_msas, [SCFV, ANTIGEN])
    return c


def mutate(seq: str, *positions: int) -> str:
    out = list(seq)
    for p in positions:
        out[p] = "A" if out[p] != "A" else "G"
    return "".join(out)


def read_rows(path: Path) -> list[list[str]]:
    with open(path, newline="") as fh:
        return list(csv.reader(fh))


def test_empty_cache_says_what_is_missing(tmp_path):
    cache = MsaCache(tmp_path / "c")
    assert not cache.has_parent()
    with pytest.raises(CacheError, match="no parent alignments"):
        cache.parent_chains


def test_import_keeps_alignments_and_chains(cache):
    assert cache.has_parent() and cache.parent_chains == [SCFV, ANTIGEN]
    info = cache.info()
    assert info["chain_lengths"] == [len(SCFV), len(ANTIGEN)] and info["alignment_rows"] == [3, 2]


def test_import_rejects_alignment_of_another_sequence(tmp_path, boltz_msas):
    with pytest.raises(CacheError, match="first alignment row"):
        MsaCache(tmp_path / "c").import_boltz_msas(boltz_msas, [mutate(SCFV, 3), ANTIGEN])


def test_import_rejects_missing_chain_file(tmp_path, boltz_msas):
    (boltz_msas / "wt_1.csv").unlink()
    with pytest.raises(CacheError, match="expected one"):
        MsaCache(tmp_path / "c").import_boltz_msas(boltz_msas, [SCFV, ANTIGEN])


def test_identical_chain_uses_the_parent_file(cache):
    assert cache.alignment_for(1, ANTIGEN) == cache.dir / "parent_1.csv"
    assert cache.alignment_for(0, SCFV) == cache.dir / "parent_0.csv"


def test_variant_gets_its_own_query_row_and_the_parents_other_rows(cache):
    variant = mutate(SCFV, 5, 20, 40, 41, 50)
    path = cache.alignment_for(0, variant)
    rows = read_rows(path)
    assert rows[1] == ["0", variant]  # query row swapped
    assert rows[2:] == read_rows(cache.dir / "parent_0.csv")[2:]  # the rest untouched
    assert rows[0] == ["key", "sequence"]
    assert cache.alignment_for(0, variant) == path  # second call reuses the file


def test_chain_too_far_from_the_parent_is_refused(cache):
    far = mutate(SCFV, *range(0, 30))  # about half the positions changed
    assert cache.alignment_for(0, far) is None
    assert cache.alignment_for(0, SCFV + "A" * 10) is None  # 10 extra residues are more than 10% of the chain
    assert cache.alignment_for(5, SCFV) is None  # no such chain


def test_adapter_patch_writes_alignment_paths_and_restores(cache, tmp_path):
    variant = mutate(SCFV, 5, 20, 40, 41, 50)
    fasta = tmp_path / "v1.fa"
    fasta.write_text(f">v1\n{variant}:{ANTIGEN}\n")
    design = Design(id="v1", path=fasta, parent_id=None, producer="manual")
    original = _boltz.write_boltz_fasta
    with use_cached_msas(cache):
        out = _boltz.write_boltz_fasta(tmp_path / "v1.fasta", design)
    assert _boltz.write_boltz_fasta is original
    lines = out.read_text().splitlines()
    assert lines[0].startswith(">A|protein|") and lines[2].startswith(">B|protein|")
    assert lines[1] == variant and lines[3] == ANTIGEN
    first_msa = Path(lines[0].split("|")[2])
    assert first_msa.is_file() and read_rows(first_msa)[1][1] == variant
    assert lines[2].split("|")[2] == str((cache.dir / "parent_1.csv").resolve())


def test_adapter_patch_restores_even_after_an_error(cache, tmp_path):
    fasta = tmp_path / "far.fa"
    fasta.write_text(f">far\n{mutate(SCFV, *range(0, 30))}:{ANTIGEN}\n")
    design = Design(id="far", path=fasta, parent_id=None, producer="manual")
    original = _boltz.write_boltz_fasta
    with pytest.raises(CacheError, match="design far: chain 1"):
        with use_cached_msas(cache):
            _boltz.write_boltz_fasta(tmp_path / "far.fasta", design)
    assert _boltz.write_boltz_fasta is original


def test_server_fallback_leaves_the_far_chain_without_a_path(cache, tmp_path):
    fasta = tmp_path / "far.fa"
    far = mutate(SCFV, *range(0, 30))
    fasta.write_text(f">far\n{far}:{ANTIGEN}\n")
    design = Design(id="far", path=fasta, parent_id=None, producer="manual")
    with use_cached_msas(cache, fallback="server"):
        lines = _boltz.write_boltz_fasta(tmp_path / "far.fasta", design).read_text().splitlines()
    assert lines[0] == ">A|protein" and lines[2].startswith(">B|protein|")


def test_unknown_fallback_is_rejected(cache):
    with pytest.raises(ValueError):
        with use_cached_msas(cache, fallback="empty"):
            pass


def test_command_line_import_and_info(tmp_path, boltz_msas, capsys):
    parent = tmp_path / "parent.fa"
    parent.write_text(f">wt\n{SCFV}:{ANTIGEN}\n")
    cache_dir = tmp_path / "c"
    assert main(["import", str(boltz_msas), str(parent), "--cache", str(cache_dir)]) == 0
    assert main(["info", "--cache", str(cache_dir)]) == 0
    assert '"chain_lengths"' in capsys.readouterr().out


@pytest.fixture
def window_cache(tmp_path):
    """A parent whose antigen alignment has gaps and a lowercase insertion, as a server alignment does."""
    msa_dir = tmp_path / "boltz_results_p" / "msa"
    write_msa(msa_dir, "p", 0, SCFV, [SCFV])
    write_msa(
        msa_dir,
        "p",
        1,
        ANTIGEN,
        [
            ANTIGEN[:20] + "xyz" + ANTIGEN[20:],  # an insertion after column 20
            "-" * 5 + ANTIGEN[5:],  # starts with gaps
            "-" * 30 + ANTIGEN[30:],  # aligns to nothing in the first 30 columns
        ],
    )
    c = MsaCache(tmp_path / "window_cache")
    c.import_boltz_msas(msa_dir, [SCFV, ANTIGEN])
    return c


def test_fragment_of_the_parent_gets_the_matching_columns_and_keeps_insertions(window_cache):
    fragment = ANTIGEN[10:30]
    rows = read_rows(window_cache.alignment_for(1, fragment))
    assert rows[0] == ["key", "sequence"] and rows[1] == ["0", fragment]
    # the insertion stays with its column, and the row that aligns to nothing in the window is dropped
    assert rows[2:] == [["1", ANTIGEN[10:20] + "xyz" + ANTIGEN[20:30]], ["2", fragment]]


def test_a_few_extra_residues_at_an_end_become_gap_columns(window_cache):
    tagged = ANTIGEN + "HHHH"  # 4 of 51 residues lie beyond the parent's chain
    rows = read_rows(window_cache.alignment_for(1, tagged))
    assert rows[1][1] == tagged
    assert rows[2][1] == ANTIGEN[:20] + "xyz" + ANTIGEN[20:] + "----"


def test_a_substitution_inside_a_fragment_is_tolerated(window_cache):
    fragment = mutate(ANTIGEN[10:30], 3)
    rows = read_rows(window_cache.alignment_for(1, fragment))
    assert rows[1][1] == fragment and len(rows) == 4


def test_what_does_not_sit_on_the_parent_is_refused(window_cache):
    assert window_cache.alignment_for(1, "W" * 20) is None
    assert window_cache.alignment_for(1, ANTIGEN + "H" * 10) is None  # more than 10% beyond the parent


def test_a_row_that_does_not_match_the_parents_length_is_an_error(tmp_path):
    msa_dir = tmp_path / "m"
    write_msa(msa_dir, "p", 0, SCFV, [SCFV])
    write_msa(msa_dir, "p", 1, ANTIGEN, ["A" * 10])
    c = MsaCache(tmp_path / "c")
    c.import_boltz_msas(msa_dir, [SCFV, ANTIGEN])
    with pytest.raises(CacheError, match="columns"):
        c.alignment_for(1, ANTIGEN[10:30])


def test_adapter_patch_points_a_fragment_antigen_at_the_sliced_alignment(window_cache, tmp_path):
    fragment = ANTIGEN[10:30]
    fasta = tmp_path / "v1.fa"
    fasta.write_text(f">v1\n{SCFV}:{fragment}\n")
    design = Design(id="v1", path=fasta, parent_id=None, producer="manual")
    with use_cached_msas(window_cache):
        lines = _boltz.write_boltz_fasta(tmp_path / "v1.fasta", design).read_text().splitlines()
    assert lines[3] == fragment
    assert read_rows(Path(lines[2].split("|")[2]))[1][1] == fragment


@pytest.fixture
def fake_boltz_module(tmp_path, monkeypatch):
    """A stand-in `boltz.main.compute_msa` that writes two-row alignments, found through PYTHONPATH."""
    package = tmp_path / "fake_site" / "boltz"
    package.mkdir(parents=True)
    (package / "__init__.py").write_text("")
    (package / "main.py").write_text(
        "import os\n"
        "def compute_msa(data, target_id, msa_dir, msa_server_url, msa_pairing_strategy):\n"
        "    if os.environ.get('FAKE_BOLTZ_FAIL'):\n"
        "        raise SystemExit(3)\n"
        "    for name, seq in data.items():\n"
        "        (msa_dir / f'{name}.csv').write_text('key,sequence\\n0,' + seq + '\\n-1,' + seq.lower() + '\\n')\n"
    )
    monkeypatch.setenv("PYTHONPATH", str(tmp_path / "fake_site"))


def test_fetch_alignments_keeps_what_boltz_saved_for_each_chain(tmp_path, fake_boltz_module):
    cache = MsaCache(tmp_path / "cache")
    fetch_alignments(cache, [SCFV, ANTIGEN], tmp_path / "work", python=sys.executable)
    assert cache.has_parent() and cache.parent_chains == [SCFV, ANTIGEN]
    assert cache.info()["alignment_rows"] == [2, 2]
    assert read_rows(cache.dir / "parent_1.csv")[1] == ["0", ANTIGEN]


def test_fetch_alignments_reports_a_failed_fetch(tmp_path, fake_boltz_module, monkeypatch):
    monkeypatch.setenv("FAKE_BOLTZ_FAIL", "1")
    with pytest.raises(CacheError, match="exit code 3"):
        fetch_alignments(MsaCache(tmp_path / "cache"), [SCFV, ANTIGEN], tmp_path / "work", python=sys.executable)


def test_boltz_python_is_read_from_the_shebang(tmp_path, monkeypatch):
    script = tmp_path / "bin" / "boltz"
    script.parent.mkdir()
    script.write_text("#!/opt/venv/bin/python3\nprint('x')\n")
    script.chmod(0o755)
    monkeypatch.setenv("PATH", str(script.parent))
    assert _boltz_python() == "/opt/venv/bin/python3"
    script.write_text("#!/usr/bin/env python3\n")
    with pytest.raises(CacheError, match="pass python="):
        _boltz_python()
    monkeypatch.setenv("PATH", str(tmp_path))
    with pytest.raises(CacheError, match="not on PATH"):
        _boltz_python()
