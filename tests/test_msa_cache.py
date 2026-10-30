"""msa_cache: importing the parent's alignments, swapping the query row, and patching the Boltz adapter."""

import csv
from pathlib import Path

import pytest

from design_pipeline.adapters import _boltz
from design_pipeline.designs import Design
from msa_cache import CacheError, MsaCache, use_cached_msas
from msa_cache.__main__ import main

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
    assert cache.alignment_for(0, SCFV + "A") is None  # different length
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
