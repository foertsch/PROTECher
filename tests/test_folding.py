"""The co-folding adapters: Chai-1, Boltz-2, Protenix.

Offline. The real tools are stubbed at their one seam each — the imported
`run_inference` for Chai-1, `stream` for the two that shell out — and the stub
writes the files that tool writes, so the conversion to PDB and the score
normalization run for real.
"""

from __future__ import annotations

import json
from dataclasses import dataclass
from pathlib import Path
from typing import Any

import numpy as np
import pytest

from design_pipeline import (
    REGISTRY,
    Design,
    MissingDependencyError,
    MissingGPUError,
    RunFailedError,
    get_adapter,
    predict_structure,
    score_designs,
)
from design_pipeline.adapters import _boltz, _chai, _dummy, _folding, _protenix
from design_pipeline.adapters._boltz import BoltzAdapter
from design_pipeline.adapters._chai import ChaiAdapter
from design_pipeline.adapters._protenix import ProtenixAdapter

COLABFOLD_KEYS = {"plddt", "ptm", "iptm"}
"""What the ColabFold adapter puts in `scores`, and so what every other
predictor has to put there too: `select` sorts one name across a mixed pool."""

TOKENS = _dummy.CHAIN_A + _dummy.CHAIN_B


@pytest.fixture
def sequence(tmp_path: Path) -> Design:
    """One sequence design, two chains, as the MPNN step writes it."""
    fasta = tmp_path / "seq.fa"
    fasta.write_text(f">bb0_sample1\n{'A' * _dummy.CHAIN_A}:{'G' * _dummy.CHAIN_B}\n")
    return Design(id="bb0_sample1", path=fasta, parent_id="bb0", producer="test")


@pytest.fixture
def on_gpu(monkeypatch):
    """A Linux host with an Ampere card, so the guard is not what is under test."""
    monkeypatch.setattr(_folding.sys, "platform", "linux")
    monkeypatch.setattr(_folding, "gpus", lambda: [("NVIDIA A100-SXM4-40GB", 8.0)])


def write_prediction_cif(path: Path, token: str) -> Path:
    """What all three tools write: an mmCIF with pLDDT in the B-factor column."""
    path.parent.mkdir(parents=True, exist_ok=True)
    staging = path.with_suffix(".staging.pdb")
    _dummy.write_pdb(staging, token=token, spread=1.2)
    _dummy.write_cif(staging, path)
    staging.unlink()
    return path


def pae_matrix(token: str) -> list[list[float]]:
    return [
        [3.0 + _dummy._jitter(f"{token}:{i}:{j}") for j in range(TOKENS)]
        for i in range(TOKENS)
    ]


# --- registry and defaults ---------------------------------------------------


def test_the_three_names_resolve_to_their_adapters() -> None:
    assert isinstance(get_adapter("chai"), ChaiAdapter)
    assert isinstance(get_adapter("boltz"), BoltzAdapter)
    assert isinstance(get_adapter("protenix"), ProtenixAdapter)


def test_each_has_a_dummy_counterpart() -> None:
    for name in ("dummy_chai", "dummy_boltz", "dummy_protenix"):
        adapter = get_adapter(name)
        assert isinstance(adapter, _dummy.DummyFold)
        assert adapter.producer == name


def test_colabfold_is_still_registered_and_unchanged() -> None:
    from design_pipeline.adapters._colabfold import ColabFoldAdapter

    assert isinstance(get_adapter("colabfold"), ColabFoldAdapter)
    assert "colabfold" in REGISTRY


def test_chai_is_the_default_folding_model() -> None:
    import inspect

    default = inspect.signature(predict_structure).parameters["model_name"].default
    assert default == "chai"


def test_importing_the_package_imports_no_model() -> None:
    import sys

    import design_pipeline  # noqa: F401

    assert not {
        m
        for m in sys.modules
        if m.split(".")[0] in {"torch", "jax", "chai_lab", "boltz", "protenix"}
    }


# --- the platform and GPU guard ----------------------------------------------


def test_the_guard_names_linux_off_linux(monkeypatch) -> None:
    monkeypatch.setattr(_folding.sys, "platform", "win32")
    with pytest.raises(MissingGPUError) as exc:
        _folding.require_gpu("chai", "chai")
    message = str(exc.value)
    assert "Linux" in message and "win32" in message
    # The common reason to meet this is wanting to run on a laptop.
    assert "dummy_chai" in message


def test_the_guard_reports_no_gpu(monkeypatch) -> None:
    monkeypatch.setattr(_folding.sys, "platform", "linux")
    monkeypatch.setattr(_folding, "gpus", lambda: [])
    with pytest.raises(MissingGPUError, match="nvidia-smi"):
        _folding.require_gpu("boltz", "boltz")


def test_the_guard_rejects_a_pre_ampere_card_for_bfloat16(monkeypatch) -> None:
    monkeypatch.setattr(_folding.sys, "platform", "linux")
    monkeypatch.setattr(_folding, "gpus", lambda: [("Tesla T4", 7.5)])
    with pytest.raises(MissingGPUError) as exc:
        _folding.require_gpu("chai", "chai", bfloat16=True)
    assert "bfloat16" in str(exc.value) and "Tesla T4" in str(exc.value)
    # Only Chai-1 needs it; Boltz runs on the same card.
    _folding.require_gpu("boltz", "boltz")


def test_an_unknown_capability_does_not_block(monkeypatch) -> None:
    """A driver too old for the compute_cap query still lists its cards, and
    refusing to run on a GPU we cannot interrogate would be worse."""
    monkeypatch.setattr(_folding.sys, "platform", "linux")
    monkeypatch.setattr(_folding, "gpus", lambda: [("NVIDIA H100", None)])
    _folding.require_gpu("chai", "chai", bfloat16=True)


def test_each_adapter_checks_the_gpu_before_the_tool(monkeypatch, tmp_path, sequence):
    """The guard comes first: no weights download, no missing-executable error."""
    monkeypatch.setattr(_folding.sys, "platform", "darwin")
    for adapter in (ChaiAdapter(), BoltzAdapter(), ProtenixAdapter()):
        with pytest.raises(MissingGPUError):
            adapter.run([sequence], tmp_path / "out")


# --- the shared helpers ------------------------------------------------------


def test_cif_to_pdb_keeps_the_b_factor_column(tmp_path) -> None:
    """Scoring falls back to the B-factor column for pLDDT, so the conversion
    has to carry it across."""
    cif = write_prediction_cif(tmp_path / "pred.cif", token="t")
    pdb = _folding.cif_to_pdb(cif)

    assert pdb.suffix == ".pdb" and pdb.is_file()
    from biotite.structure.io.pdb import PDBFile

    atoms = PDBFile.read(str(pdb)).get_structure(model=1, extra_fields=["b_factor"])
    assert atoms.array_length() == TOKENS
    assert sorted(set(atoms.chain_id)) == ["A", "B"]
    assert float(atoms.b_factor.mean()) > 0


def test_fasta_chains_splits_on_the_chain_separator(tmp_path, sequence) -> None:
    name, chains = _folding.fasta_chains(sequence.path)
    assert name == "bb0_sample1"
    assert [len(c) for c in chains] == [_dummy.CHAIN_A, _dummy.CHAIN_B]


def test_fasta_chains_refuses_a_multi_record_file(tmp_path) -> None:
    fasta = tmp_path / "two.fa"
    fasta.write_text(">a\nAAAA\n>b\nGGGG\n")
    with pytest.raises(ValueError, match="one design per file"):
        _folding.fasta_chains(fasta)


def test_work_dir_clears_a_leftover_attempt(tmp_path, sequence) -> None:
    first = _folding.work_dir(tmp_path, sequence)
    (first / "half_written.cif").write_text("junk")
    again = _folding.work_dir(tmp_path, sequence)
    assert again == first
    assert list(again.iterdir()) == []


# --- chai --------------------------------------------------------------------


@dataclass
class FakeCandidates:
    """What `run_inference` returns, with the fields the adapter reads."""

    cif_paths: list[Path]
    ranking_data: list[Any]
    plddt: Any
    pae: Any

    def sorted(self) -> "FakeCandidates":
        return self


@dataclass
class FakeRanking:
    aggregate_score: float
    ptm_scores: Any


@dataclass
class FakePTM:
    complex_ptm: float
    interface_ptm: float


@pytest.fixture
def chai_ran(monkeypatch):
    """Stub `run_inference`: write the mmCIF and the npz, return the scores."""
    calls: list[dict[str, Any]] = []

    def run_inference(**kwargs):
        calls.append(kwargs)
        work = Path(kwargs["output_dir"])
        # chai_lab's own assertion, kept here: it refuses to write into a
        # directory that is not empty, so the converted input must not be in it.
        assert not work.exists() or not any(work.iterdir())
        cif = write_prediction_cif(work / "pred.model_idx_0.cif", token="chai")
        np.savez(work / "scores.model_idx_0.npz", ptm=np.array(0.81))
        return FakeCandidates(
            cif_paths=[cif],
            ranking_data=[
                FakeRanking(
                    aggregate_score=0.77,
                    ptm_scores=FakePTM(complex_ptm=0.81, interface_ptm=0.64),
                )
            ],
            # Chai-1 reports pLDDT per token as a fraction.
            plddt=np.full((1, TOKENS), 0.9),
            pae=np.array([pae_matrix("chai")]),
        )

    monkeypatch.setattr(_chai, "load_run_inference", lambda: run_inference)
    return calls


def test_chai_normalizes_scores_and_hands_back_a_pdb(
    tmp_path, sequence, on_gpu, chai_ran
) -> None:
    out = tmp_path / "fold"
    out.mkdir()
    designs = ChaiAdapter().run([sequence], out)

    assert len(designs) == 1
    design = designs[0]
    assert design.producer == "chai"
    assert design.parent_id == sequence.id
    assert design.path.suffix == ".pdb" and design.path.is_file()
    assert Path(design.meta["cif"]).is_file()
    assert Path(design.meta["raw_scores"]).suffix == ".npz"

    assert COLABFOLD_KEYS <= set(design.scores)
    # 0.9 per token on Chai-1's scale is 90 on ColabFold's.
    assert design.scores["plddt"] == pytest.approx(90.0)
    assert design.scores["ptm"] == pytest.approx(0.81)
    assert design.scores["iptm"] == pytest.approx(0.64)
    # Its own ranking score comes through prefixed, not as a ColabFold key.
    assert design.scores["chai_aggregate_score"] == pytest.approx(0.77)


def test_chai_writes_the_pae_matrix_where_scoring_reads_it(
    tmp_path, sequence, on_gpu, chai_ran
) -> None:
    """The reason this adapter runs in process: its command line drops the PAE."""
    out = tmp_path / "fold"
    out.mkdir()
    design = ChaiAdapter().run([sequence], out)[0]

    payload = json.loads(Path(design.meta["scores_json"]).read_text())
    assert len(payload["pae"]) == TOKENS
    assert payload["ptm"] == pytest.approx(0.81)
    assert len(payload["plddt"]) == TOKENS and payload["plddt"][0] == pytest.approx(90.0)


def test_chai_is_msa_free_unless_asked(tmp_path, sequence, on_gpu, chai_ran) -> None:
    out = tmp_path / "fold"
    out.mkdir()
    ChaiAdapter().run([sequence], out)
    assert chai_ran[0]["use_msa_server"] is False
    assert chai_ran[0]["seed"] == 0
    # The tool's own default sample count stands unless the caller sets ours.
    assert "num_diffn_samples" not in chai_ran[0]

    ChaiAdapter().run([sequence], out, use_msa_server=True, num_diffusion_samples=2)
    assert chai_ran[1]["use_msa_server"] is True
    assert chai_ran[1]["num_diffn_samples"] == 2


def test_chai_converts_the_fasta_to_one_record_per_chain(
    tmp_path, sequence, on_gpu, chai_ran
) -> None:
    out = tmp_path / "fold"
    out.mkdir()
    ChaiAdapter().run([sequence], out)
    written = Path(chai_ran[0]["fasta_file"]).read_text()
    assert written.count(">") == 2
    assert ">protein|name=bb0_sample1-A" in written
    assert ":" not in written


def test_chai_says_what_to_do_with_extra_args(
    tmp_path, sequence, on_gpu, chai_ran
) -> None:
    with pytest.raises(ValueError, match="keyword arguments"):
        ChaiAdapter().run([sequence], tmp_path / "fold", extra_args=["--fast"])


def test_chai_missing_names_the_extra(monkeypatch) -> None:
    def no_module():
        raise MissingDependencyError(
            "chai_lab is not installed. Install the optional extra with "
            "`uv pip install -e '.[chai]'`."
        )

    monkeypatch.setattr(_chai, "load_run_inference", no_module)
    with pytest.raises(MissingDependencyError, match=r"\[chai\]"):
        _chai.load_run_inference()


# --- boltz -------------------------------------------------------------------


@pytest.fixture
def boltz_ran(monkeypatch, tmp_path):
    """Record the command; write what `boltz predict` writes."""
    calls: list[list[str]] = []
    monkeypatch.setattr(_boltz.shutil, "which", lambda name: f"/usr/bin/{name}")
    monkeypatch.setenv("BOLTZ_CACHE", str(tmp_path / "boltz_cache"))

    def run(cmd):
        cmd = [str(c) for c in cmd]
        calls.append(cmd)
        out = Path(cmd[cmd.index("--out_dir") + 1])
        record = Path(cmd[2]).stem
        predictions = out / f"boltz_results_{record}" / "predictions" / record
        write_prediction_cif(predictions / f"{record}_model_0.cif", token="boltz")
        (predictions / f"confidence_{record}_model_0.json").write_text(
            json.dumps(
                {
                    "confidence_score": 0.71,
                    "ptm": 0.66,
                    "iptm": 0.52,
                    # Boltz reports pLDDT as a fraction.
                    "complex_plddt": 0.83,
                    "complex_iplddt": 0.79,
                }
            )
        )
        np.savez_compressed(
            predictions / f"pae_{record}_model_0.npz",
            pae=np.array(pae_matrix("boltz")),
        )
        np.savez_compressed(
            predictions / f"plddt_{record}_model_0.npz",
            plddt=np.full(TOKENS, 0.83),
        )

    monkeypatch.setattr(_boltz, "stream", run)
    return calls


def test_boltz_normalizes_scores_and_hands_back_a_pdb(
    tmp_path, sequence, on_gpu, boltz_ran
) -> None:
    out = tmp_path / "fold"
    out.mkdir()
    design = BoltzAdapter().run([sequence], out)[0]

    assert design.producer == "boltz"
    assert design.parent_id == sequence.id
    assert design.path.suffix == ".pdb" and design.path.is_file()
    assert Path(design.meta["cif"]).suffix == ".cif"
    assert Path(design.meta["cif"]).is_file()
    assert Path(design.meta["raw_scores"]).name.startswith("confidence_")

    assert COLABFOLD_KEYS <= set(design.scores)
    assert design.scores["plddt"] == pytest.approx(83.0)  # 0.83 -> 0-100
    assert design.scores["ptm"] == pytest.approx(0.66)
    assert design.scores["iptm"] == pytest.approx(0.52)
    assert design.scores["boltz_confidence_score"] == pytest.approx(0.71)

    payload = json.loads(Path(design.meta["scores_json"]).read_text())
    assert len(payload["pae"]) == TOKENS
    assert payload["plddt"][0] == pytest.approx(83.0)


def test_boltz_command(tmp_path, sequence, on_gpu, boltz_ran) -> None:
    out = tmp_path / "fold"
    out.mkdir()
    BoltzAdapter().run([sequence], out, recycling_steps=5)
    cmd = boltz_ran[0]

    assert cmd[0].endswith("boltz") and cmd[1] == "predict"
    assert cmd[2].endswith(".fasta")
    assert cmd[cmd.index("--output_format") + 1] == "mmcif"
    # Boltz expects MSAs, so the server is on unless turned off.
    assert "--use_msa_server" in cmd
    # pae_interaction needs the matrix, which the summary does not carry.
    assert "--write_full_pae" in cmd
    assert cmd[cmd.index("--seed") + 1] == "0"
    assert cmd[cmd.index("--recycling_steps") + 1] == "5"
    assert cmd[cmd.index("--cache") + 1].endswith("boltz_cache")

    BoltzAdapter().run([sequence], out, use_msa_server=False)
    assert "--use_msa_server" not in boltz_ran[1]


def test_boltz_converts_the_fasta_to_boltz_headers(
    tmp_path, sequence, on_gpu, boltz_ran
) -> None:
    out = tmp_path / "fold"
    out.mkdir()
    BoltzAdapter().run([sequence], out)
    written = Path(boltz_ran[0][2]).read_text()
    assert ">A|protein" in written and ">B|protein" in written
    assert ":" not in written


def test_boltz_creates_its_weights_cache(tmp_path, sequence, on_gpu, boltz_ran) -> None:
    """The first run downloads several GB; the directory exists before it starts."""
    out = tmp_path / "fold"
    out.mkdir()
    BoltzAdapter().run([sequence], out)
    assert (tmp_path / "boltz_cache").is_dir()


def test_boltz_rejects_a_relative_cache(monkeypatch, tmp_path, sequence, on_gpu):
    monkeypatch.setattr(_boltz.shutil, "which", lambda name: f"/usr/bin/{name}")
    monkeypatch.setenv("BOLTZ_CACHE", "relative/cache")
    with pytest.raises(MissingDependencyError, match="absolute"):
        BoltzAdapter().run([sequence], tmp_path / "fold")


def test_boltz_missing_names_the_extra(monkeypatch, tmp_path, sequence, on_gpu):
    monkeypatch.setattr(_boltz.shutil, "which", lambda name: None)
    with pytest.raises(MissingDependencyError) as exc:
        BoltzAdapter().run([sequence], tmp_path / "fold")
    assert "[boltz]" in str(exc.value) and "dummy_boltz" in str(exc.value)


def test_boltz_reports_an_empty_run(monkeypatch, tmp_path, sequence, on_gpu):
    monkeypatch.setattr(_boltz.shutil, "which", lambda name: f"/usr/bin/{name}")
    monkeypatch.setattr(_boltz, "stream", lambda cmd: None)
    with pytest.raises(RunFailedError, match="no mmCIF"):
        BoltzAdapter().run([sequence], tmp_path / "fold")


# --- protenix ----------------------------------------------------------------


@pytest.fixture
def protenix_ran(monkeypatch):
    """Record the command; write what `protenix predict` writes."""
    calls: list[list[str]] = []
    monkeypatch.setattr(_protenix.shutil, "which", lambda name: f"/usr/bin/{name}")

    def run(cmd):
        cmd = [str(c) for c in cmd]
        calls.append(cmd)
        out = Path(cmd[cmd.index("--out_dir") + 1])
        seed = cmd[cmd.index("--seeds") + 1]
        job = json.loads(Path(cmd[cmd.index("--input") + 1]).read_text())[0]
        name = job["name"]
        predictions = out / "default" / name / f"seed_{seed}" / "predictions"
        write_prediction_cif(predictions / f"{name}_sample_0.cif", token="protenix")
        (predictions / f"{name}_summary_confidence_sample_0.json").write_text(
            # Protenix reports pLDDT on the 0-100 scale already.
            json.dumps(
                {"plddt": 88.5, "ptm": 0.74, "iptm": 0.61, "ranking_score": 0.69}
            )
        )
        (predictions / f"{name}_full_data_sample_0.json").write_text(
            json.dumps({"token_pair_pae": pae_matrix("protenix")})
        )

    monkeypatch.setattr(_protenix, "stream", run)
    return calls


def test_protenix_normalizes_scores_and_hands_back_a_pdb(
    tmp_path, sequence, on_gpu, protenix_ran
) -> None:
    out = tmp_path / "fold"
    out.mkdir()
    design = ProtenixAdapter().run([sequence], out)[0]

    assert design.producer == "protenix"
    assert design.parent_id == sequence.id
    assert design.path.suffix == ".pdb" and design.path.is_file()
    assert Path(design.meta["cif"]).is_file()
    assert "summary_confidence" in Path(design.meta["raw_scores"]).name
    assert design.meta["checkpoint"] == _protenix.BASE_CHECKPOINT

    assert COLABFOLD_KEYS <= set(design.scores)
    assert design.scores["plddt"] == pytest.approx(88.5)
    assert design.scores["protenix_ranking_score"] == pytest.approx(0.69)

    payload = json.loads(Path(design.meta["scores_json"]).read_text())
    assert len(payload["pae"]) == TOKENS


def test_protenix_converts_the_fasta_to_its_json_input(
    tmp_path, sequence, on_gpu, protenix_ran
) -> None:
    out = tmp_path / "fold"
    out.mkdir()
    ProtenixAdapter().run([sequence], out)
    cmd = protenix_ran[0]
    job = json.loads(Path(cmd[cmd.index("--input") + 1]).read_text())

    assert len(job) == 1 and job[0]["name"] == sequence.id
    chains = [entry["proteinChain"] for entry in job[0]["sequences"]]
    assert [len(c["sequence"]) for c in chains] == [_dummy.CHAIN_A, _dummy.CHAIN_B]
    assert {c["count"] for c in chains} == {1}


def test_protenix_takes_a_checkpoint(tmp_path, sequence, on_gpu, protenix_ran) -> None:
    out = tmp_path / "fold"
    out.mkdir()
    ProtenixAdapter().run([sequence], out)
    cmd = protenix_ran[0]
    assert cmd[cmd.index("--model_name") + 1] == _protenix.BASE_CHECKPOINT
    # Its bool options take a value, so True is not a bare switch here.
    assert cmd[cmd.index("--use_msa") + 1] == "True"
    assert cmd[cmd.index("--need_atom_confidence") + 1] == "True"

    design = ProtenixAdapter().run(
        [sequence], out, checkpoint=_protenix.MINI_CHECKPOINT, use_msa=False
    )[0]
    cmd = protenix_ran[1]
    assert cmd[cmd.index("--model_name") + 1] == _protenix.MINI_CHECKPOINT
    assert cmd[cmd.index("--use_msa") + 1] == "False"
    assert design.meta["checkpoint"] == _protenix.MINI_CHECKPOINT


def test_protenix_missing_names_the_extra(monkeypatch, tmp_path, sequence, on_gpu):
    monkeypatch.setattr(_protenix.shutil, "which", lambda name: None)
    with pytest.raises(MissingDependencyError) as exc:
        ProtenixAdapter().run([sequence], tmp_path / "fold")
    assert "[protenix]" in str(exc.value) and "dummy_protenix" in str(exc.value)


# --- the fixtures, and the list form across two models -----------------------


def test_the_dummies_report_the_same_keys_as_the_real_adapters(tmp_path, sequence):
    """What the offline suite is worth depends on the fixtures matching."""
    out = tmp_path / "fold"
    out.mkdir()
    design = _dummy.DummyFold(producer="dummy_chai").run([sequence], out)[0]

    assert COLABFOLD_KEYS <= set(design.scores)
    assert 0 <= design.scores["plddt"] <= 100
    assert design.path.suffix == ".pdb"
    assert Path(design.meta["cif"]).suffix == ".cif"
    assert json.loads(Path(design.meta["scores_json"]).read_text())["pae"]


def test_a_model_list_folds_every_sequence_with_both_models(tmp_path) -> None:
    """Two predictors over one pool of sequences, which is what the list form
    is for: both producers come back, each prediction on its own sequence."""
    from design_pipeline import design_sequences, make_protein_backbones

    target = tmp_path / "target.pdb"
    target.write_text("ATOM      1  CA  GLY A   1       0.000   0.000   0.000\n")
    backbones = make_protein_backbones(
        target, tmp_path / "bb", model_name="dummy_backbones", num_designs=2
    )
    seqs = design_sequences(
        backbones, tmp_path / "seqs", model_name="dummy_sequences", num_seqs=1
    )
    preds = predict_structure(
        seqs, tmp_path / "fold", model_name=["dummy_chai", "dummy_protenix"]
    )

    assert len(preds) == 2 * len(seqs)
    assert [d.producer for d in preds] == ["dummy_chai"] * 2 + ["dummy_protenix"] * 2
    assert [d.parent_id for d in preds] == [s.id for s in seqs] * 2
    assert len({d.id for d in preds}) == len(preds)
    assert all(d.path.suffix == ".pdb" and d.path.is_file() for d in preds)

    # Both producers reach the metrics, which is the point of normalizing keys.
    scored = score_designs(preds, reference=backbones)
    assert all("pae_interaction" in d.scores for d in scored)
    assert all("sc_rmsd" in d.scores for d in scored)
