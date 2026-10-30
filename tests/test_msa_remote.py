"""msa_cache.remote: staging a batch, folding it from the staged folder, and fixing paths after a download."""

import json

import pytest

from design_pipeline.cache import read_manifest
from msa_cache import CacheError
from msa_cache.remote import ORDER, localize_manifest, read_chains, run_staged, stage
from test_msa_batch import designs, fake_boltz, sequence_design  # noqa: F401  (fixtures)
from test_msa_cache import ANTIGEN, SCFV, boltz_msas, cache, mutate  # noqa: F401  (fixtures)


@pytest.fixture
def design_dir(designs):
    return designs[0].path.parent


def test_read_chains_splits_on_colon_and_refuses_several_records(tmp_path):
    (tmp_path / "a.fa").write_text(f">a\n{SCFV}:{ANTIGEN}\n")
    assert read_chains(tmp_path / "a.fa") == [SCFV, ANTIGEN]
    (tmp_path / "b.fa").write_text(">a\nAAA\n>b\nCCC\n")
    with pytest.raises(CacheError, match="2 FASTA records"):
        read_chains(tmp_path / "b.fa")


def test_stage_copies_designs_and_the_parents_alignments_only(cache, design_dir, tmp_path):
    cache.alignment_for(0, mutate(SCFV, 5, 20, 40, 41, 50))  # makes a variants/ file in the source cache
    ids = stage(design_dir, cache.dir, tmp_path / "stage")
    assert ids == ["v1", "v2", "wt"]
    staged = tmp_path / "stage"
    assert sorted(p.name for p in (staged / "designs").iterdir()) == [ORDER, "v1.fa", "v2.fa", "wt.fa"]
    assert sorted(p.name for p in (staged / "cache").iterdir()) == ["parent.json", "parent_0.csv", "parent_1.csv"]
    assert (staged / "designs" / ORDER).read_text().split() == ids


def test_stage_keeps_the_requested_folding_order(cache, design_dir, tmp_path):
    assert stage(design_dir, cache.dir, tmp_path / "stage", order=["wt", "v2", "v1"]) == ["wt", "v2", "v1"]
    assert (tmp_path / "stage" / "designs" / ORDER).read_text().split() == ["wt", "v2", "v1"]


def test_stage_refuses_before_anything_is_copied(cache, design_dir, tmp_path):
    with pytest.raises(CacheError, match="does not match the designs"):
        stage(design_dir, cache.dir, tmp_path / "stage", order=["wt", "v1"])
    with pytest.raises(CacheError, match="no parent alignments"):
        stage(design_dir, tmp_path / "empty", tmp_path / "stage")
    with pytest.raises(CacheError, match="no \\*.fa files"):
        stage(tmp_path / "nothing_here", cache.dir, tmp_path / "stage")
    sequence_design(design_dir, "far", mutate(SCFV, *range(0, 30)))
    with pytest.raises(CacheError, match="design far: chain 1"):
        stage(design_dir, cache.dir, tmp_path / "stage")
    assert not (tmp_path / "stage").exists()


def test_run_staged_folds_in_the_staged_order_with_one_process(cache, design_dir, fake_boltz, tmp_path):
    stage(design_dir, cache.dir, tmp_path / "run", order=["v2", "wt", "v1"])
    folded = run_staged(tmp_path / "run", seed=3)
    assert [d.parent_id for d in folded] == ["v2", "wt", "v1"]
    assert len(fake_boltz.read_text().splitlines()) == 1 and "--seed 3" in fake_boltz.read_text()
    assert [d.parent_id for d in read_manifest(tmp_path / "run" / "out")] == ["v2", "wt", "v1"]


def test_localize_manifest_rewrites_only_paths_under_the_remote_folder(tmp_path):
    out = tmp_path / "out"
    out.mkdir()
    (out / "DONE").write_text(json.dumps([{
        "id": "a", "path": "x.pdb", "parent_id": "wt", "producer": "boltz", "scores": {"iptm": 0.9},
        "meta": {"cif": "/work/run1/out/p/a.cif", "scores_json": "/work/run1/out", "other": "/work/run10/out/b", "seed": 0, "note": "x"},
    }]))
    assert localize_manifest(out, "/work/run1/out/") == 2
    meta = json.loads((out / "DONE").read_text())[0]["meta"]
    assert meta["cif"] == f"{out.resolve()}/p/a.cif" and meta["scores_json"] == str(out.resolve())
    assert meta["other"] == "/work/run10/out/b" and meta["seed"] == 0 and meta["note"] == "x"
