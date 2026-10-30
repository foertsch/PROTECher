"""msa_cache.batch: one Boltz process for many designs, with a stand-in `boltz` executable."""

import inspect
import os
import stat
from pathlib import Path

import pytest

from design_pipeline.adapters._boltz import BoltzAdapter
from design_pipeline.cache import read_manifest
from design_pipeline.designs import Design
from design_pipeline.errors import MissingDependencyError, RunFailedError
from msa_cache import CacheError, batch, fold_batch, write_batch_inputs
from test_msa_cache import ANTIGEN, SCFV, boltz_msas, cache, mutate  # noqa: F401  (fixtures)

FAKE_BOLTZ = """#!/usr/bin/env bash
echo "$@" >> "$BOLTZ_CALLS"
inputs=$2
while [ $# -gt 0 ]; do [ "$1" = "--out_dir" ] && out=$2; shift; done
for f in "$inputs"/*.fasta; do
  id=$(basename "$f" .fasta)
  [ "$id" = "${BOLTZ_SKIP:-}" ] && continue
  mkdir -p "$out/boltz_results_inputs/predictions/$id"
done
"""


def sequence_design(directory: Path, design_id: str, scfv: str, antigen: str = ANTIGEN) -> Design:
    directory.mkdir(parents=True, exist_ok=True)
    path = directory / f"{design_id}.fa"
    path.write_text(f">{design_id}\n{scfv}:{antigen}\n")
    return Design(id=design_id, path=path, parent_id=None, producer="manual")


@pytest.fixture
def designs(tmp_path):
    seqs = tmp_path / "seqs"
    return [
        sequence_design(seqs, "wt", SCFV),
        sequence_design(seqs, "v1", mutate(SCFV, 5, 20, 40, 41, 50)),
        sequence_design(seqs, "v2", mutate(SCFV, 7, 21, 42, 43, 52)),
    ]


@pytest.fixture
def fake_boltz(tmp_path, monkeypatch):
    """A `boltz` on PATH that records its arguments and writes an empty prediction folder per input."""
    bin_dir = tmp_path / "bin"
    bin_dir.mkdir()
    exe = bin_dir / "boltz"
    exe.write_text(FAKE_BOLTZ)
    exe.chmod(exe.stat().st_mode | stat.S_IEXEC)
    calls = tmp_path / "calls.txt"
    monkeypatch.setenv("PATH", f"{bin_dir}{os.pathsep}{os.environ['PATH']}")
    monkeypatch.setenv("BOLTZ_CALLS", str(calls))
    monkeypatch.setenv("BOLTZ_CACHE", str(tmp_path / "boltz_cache"))

    def collect(prediction_dir, output_dir, design):
        pdb = prediction_dir / f"{design.id}_model_0.pdb"
        pdb.write_text("ATOM\n")
        return Design(
            id=f"{pdb.stem}-{output_dir.name[:8]}",
            path=pdb,
            parent_id=design.id,
            producer="boltz",
            scores={"iptm": 0.9, "plddt": 90.0},
            meta={},
        )

    monkeypatch.setattr(batch, "_collect", collect)
    return calls


def test_inputs_point_every_chain_at_the_cache(cache, designs, tmp_path):
    files = write_batch_inputs(cache, designs, tmp_path / "inputs")
    assert [f.name for f in files] == ["wt.fasta", "v1.fasta", "v2.fasta"]
    wt = files[0].read_text().splitlines()
    assert wt == [
        f">A|protein|{(cache.dir / 'parent_0.csv').resolve()}",
        SCFV,
        f">B|protein|{(cache.dir / 'parent_1.csv').resolve()}",
        ANTIGEN,
    ]
    v1 = files[1].read_text().splitlines()
    assert v1[0].startswith(">A|protein|") and "variants" in v1[0] and v1[1] == mutate(SCFV, 5, 20, 40, 41, 50)
    assert v1[2] == wt[2]  # the antigen is unchanged, so it shares the parent's alignment


def test_a_design_too_far_from_the_parent_fails_before_anything_is_written(cache, designs, tmp_path):
    far = sequence_design(tmp_path / "seqs", "far", mutate(SCFV, *range(0, 30)))
    inputs = tmp_path / "inputs"
    with pytest.raises(CacheError, match="design far: chain 1"):
        write_batch_inputs(cache, [*designs, far], inputs)
    assert not inputs.exists()


def test_duplicate_and_unusable_ids_are_refused(cache, designs, tmp_path):
    with pytest.raises(CacheError, match="unique"):
        write_batch_inputs(cache, [designs[0], designs[0]], tmp_path / "inputs")
    odd = Design(id="a/b", path=designs[0].path, parent_id=None, producer="manual")
    with pytest.raises(CacheError, match="file name"):
        write_batch_inputs(cache, [odd], tmp_path / "inputs")


def test_one_boltz_process_folds_every_design(cache, designs, fake_boltz, tmp_path):
    out = tmp_path / "run" / "batch"
    folded = fold_batch(cache, designs, out, seed=7)
    calls = fake_boltz.read_text().splitlines()
    assert len(calls) == 1  # one process for three designs
    command = calls[0]
    assert f"predict {out / 'inputs'} " in command
    for flag in ("--seed 7", "--write_full_pae", "--output_format mmcif", f"--cache {tmp_path / 'boltz_cache'}"):
        assert flag in command
    assert "--use_msa_server" not in command
    assert [d.parent_id for d in folded] == ["wt", "v1", "v2"]
    assert all(d.producer == "boltz" and d.scores["iptm"] == 0.9 for d in folded)
    assert [d.id for d in read_manifest(out)] == [d.id for d in folded]  # run_report reads this file


def test_a_missing_prediction_fails_the_whole_call(cache, designs, fake_boltz, tmp_path, monkeypatch):
    monkeypatch.setenv("BOLTZ_SKIP", "v2")
    with pytest.raises(RunFailedError, match="'v2'"):
        fold_batch(cache, designs, tmp_path / "out")
    assert not (tmp_path / "out" / "DONE").exists()


def test_without_boltz_on_path_the_error_says_so(cache, designs, tmp_path, monkeypatch):
    monkeypatch.setenv("PATH", str(tmp_path))
    with pytest.raises(MissingDependencyError, match="not on PATH"):
        fold_batch(cache, designs, tmp_path / "out")


def test_nothing_to_fold_is_an_error(cache, tmp_path):
    with pytest.raises(CacheError, match="no designs"):
        fold_batch(cache, [], tmp_path / "out")


def test_a_second_batch_into_the_same_folder_replaces_the_first(cache, designs, fake_boltz, tmp_path):
    out = tmp_path / "out"
    fold_batch(cache, designs, out)
    fold_batch(cache, designs[:1], out)
    assert [p.name for p in (out / "boltz_results_inputs" / "predictions").iterdir()] == ["wt"]
    assert [p.name for p in (out / "inputs").iterdir()] == ["wt.fasta"]


def test_adapter_collect_keeps_the_signature_this_module_calls():
    """If this fails, the adapter changed: update batch._collect rather than the test."""
    assert list(inspect.signature(BoltzAdapter._collect).parameters) == ["self", "work", "output_dir", "sequence"]
