"""The MPNN sequence-design step: one pip-installed tool, three models."""

from __future__ import annotations

import inspect
from pathlib import Path

import pytest

from design_pipeline import (
    MissingDependencyError,
    MissingLigandError,
    Design,
    design_sequences,
    get_adapter,
)
from design_pipeline.adapters import _mpnn
from design_pipeline.adapters._mpnn import MODEL_TYPES, MPNNAdapter, header_scores

FASTA = (
    ">bb0, T=0.1, seed=0, num_res=22, num_ligand_res=0, use_ligand_context=True\n"
    "MKTAYIAKQRQISFVKSHFSRQ\n"
    ">bb0, id=1, T=0.1, seed=0, overall_confidence=0.4301,"
    " ligand_confidence=0.4301, seq_rec=0.3656\n"
    "MKTAYIAKQRQISFVKSHFSRX\n"
    ">bb0, id=2, T=0.1, seed=0, overall_confidence=0.5102,"
    " ligand_confidence=0.5102, seq_rec=0.4011\n"
    "MKTAYIAKQRQISFVKSHFSRY\n"
)


def backbone(tmp_path: Path, het: str = "", id: str = "bb0") -> Design:
    pdb = tmp_path / f"{id}.pdb"
    pdb.write_text(
        "ATOM      1  CA  GLY A   1       0.000   0.000   0.000  1.00 50.00           C\n"
        + het
        + "END\n"
    )
    return Design(id=id, path=pdb, parent_id=None, producer="test")


@pytest.fixture
def ran(monkeypatch, tmp_path):
    """Record commands; pretend the tool wrote its seqs/ FASTA."""
    calls: list[list[str]] = []
    monkeypatch.setattr(_mpnn.shutil, "which", lambda name: f"/usr/bin/{name}")

    def run(cmd):
        calls.append([str(c) for c in cmd])
        out = Path(cmd[cmd.index("--out_folder") + 1])
        stem = Path(cmd[cmd.index("--pdb_path") + 1]).stem
        (out / "seqs").mkdir(parents=True, exist_ok=True)
        (out / "seqs" / f"{stem}.fa").write_text(FASTA)

    monkeypatch.setattr(_mpnn, "stream", run)
    return calls


# --- registry and defaults ---------------------------------------------------


def test_three_names_resolve_to_one_adapter() -> None:
    adapters = {name: get_adapter(name) for name in MODEL_TYPES}
    assert {type(a) for a in adapters.values()} == {MPNNAdapter}
    assert [a.model_type for a in adapters.values()] == list(MODEL_TYPES)
    # Same executable, same code path: only --model_type differs.
    assert {a.executable for a in adapters.values()} == {"ligandmpnn"}


def test_default_model_is_protein_mpnn() -> None:
    default = inspect.signature(design_sequences).parameters["model_name"].default
    assert default == "protein_mpnn"


def test_the_container_adapter_is_gone() -> None:
    from design_pipeline import REGISTRY

    assert "proteinmpnn" not in REGISTRY
    assert "protein_mpnn" in REGISTRY


# --- command assembly --------------------------------------------------------


def test_command_carries_model_type_seed_and_paths(ran, tmp_path) -> None:
    bb = backbone(tmp_path)
    out = tmp_path / "out"
    out.mkdir()

    get_adapter("soluble_mpnn").run([bb], out, num_seqs=4, temperature=0.3, seed=7)

    cmd = ran[0]
    assert cmd[0].endswith("ligandmpnn")
    assert cmd[cmd.index("--model_type") + 1] == "soluble_mpnn"
    assert cmd[cmd.index("--seed") + 1] == "7"
    assert cmd[cmd.index("--pdb_path") + 1] == str(bb.path)
    assert cmd[cmd.index("--out_folder") + 1] == str(out)
    assert cmd[cmd.index("--temperature") + 1] == "0.3"
    # num_seqs is batch_size x number_of_batches.
    assert cmd[cmd.index("--batch_size") + 1] == "4"
    assert cmd[cmd.index("--number_of_batches") + 1] == "1"


def test_explicit_batching_overrides_num_seqs(ran, tmp_path) -> None:
    out = tmp_path / "out"
    out.mkdir()
    get_adapter("protein_mpnn").run(
        [backbone(tmp_path)], out, num_seqs=15, batch_size=3, number_of_batches=5
    )
    cmd = ran[0]
    assert cmd[cmd.index("--batch_size") + 1] == "3"
    assert cmd[cmd.index("--number_of_batches") + 1] == "5"
    assert "--num_seqs" not in cmd


def test_unanticipated_kwargs_and_extra_args_pass_through(ran, tmp_path) -> None:
    out = tmp_path / "out"
    out.mkdir()
    get_adapter("protein_mpnn").run(
        [backbone(tmp_path)],
        out,
        checkpoint_protein_mpnn="proteinmpnn_v_48_020.pt",
        bias_AA="W:1.0",
        extra_args=["--zero_indexed", "1"],
    )
    cmd = ran[0]
    assert cmd[cmd.index("--checkpoint_protein_mpnn") + 1] == "proteinmpnn_v_48_020.pt"
    assert cmd[cmd.index("--bias_AA") + 1] == "W:1.0"
    assert cmd[-2:] == ["--zero_indexed", "1"]


def test_one_invocation_per_backbone(ran, tmp_path) -> None:
    out = tmp_path / "out"
    out.mkdir()
    backbones = [backbone(tmp_path, id=f"bb{i}") for i in range(3)]
    designs = get_adapter("protein_mpnn").run(backbones, out, num_seqs=2)
    assert len(ran) == 3
    assert [d.parent_id for d in designs] == ["bb0", "bb0", "bb1", "bb1", "bb2", "bb2"]


# --- outputs -----------------------------------------------------------------


def test_fasta_is_split_per_sequence_with_scores(ran, tmp_path) -> None:
    out = tmp_path / "out"
    out.mkdir()
    designs = get_adapter("protein_mpnn").run([backbone(tmp_path)], out, num_seqs=2)

    # The input sequence record is not a design.
    assert len(designs) == 2
    assert all(d.path.suffix == ".fa" for d in designs)
    assert all(d.path.is_file() for d in designs)
    assert designs[0].path.read_text().endswith("MKTAYIAKQRQISFVKSHFSRX\n")
    assert designs[0].scores["overall_confidence"] == 0.4301
    assert designs[0].scores["seq_rec"] == 0.3656
    assert designs[1].scores["overall_confidence"] == 0.5102
    assert designs[0].producer == "protein_mpnn"
    assert designs[0].meta["model_type"] == "protein_mpnn"
    assert designs[0].meta["temperature"] == 0.1
    # path stays on the FASTA so downstream folding is unaffected.
    assert designs[0].meta["fasta"].endswith("bb0.fa")
    assert "packed_pdb" not in designs[0].meta


def test_packed_structure_is_recorded_when_present(monkeypatch, tmp_path) -> None:
    out = tmp_path / "out"
    (out / "packed").mkdir(parents=True)
    monkeypatch.setattr(_mpnn.shutil, "which", lambda name: f"/usr/bin/{name}")

    def run(cmd):
        (out / "seqs").mkdir(parents=True, exist_ok=True)
        (out / "seqs" / "bb0.fa").write_text(FASTA)
        (out / "packed" / "bb0_packed_1.pdb").write_text("ATOM\n")

    monkeypatch.setattr(_mpnn, "stream", run)
    designs = get_adapter("protein_mpnn").run(
        [backbone(tmp_path)], out, num_seqs=2, pack_side_chains=1
    )
    assert designs[0].meta["packed_pdb"].endswith("bb0_packed_1.pdb")


def test_header_parsing_skips_non_numeric_fields() -> None:
    scores = header_scores(
        ">bb0, id=1, T=0.1, overall_confidence=0.43, use_ligand_context=True,"
        " model_path=/x/y.pt"
    )
    assert scores == {"id": 1.0, "T": 0.1, "overall_confidence": 0.43}


def test_silent_failure_is_an_error(monkeypatch, tmp_path) -> None:
    monkeypatch.setattr(_mpnn.shutil, "which", lambda name: f"/usr/bin/{name}")
    monkeypatch.setattr(_mpnn, "stream", lambda cmd: None)
    out = tmp_path / "out"
    out.mkdir()
    with pytest.raises(_mpnn.RunFailedError):
        get_adapter("protein_mpnn").run([backbone(tmp_path)], out)


# --- ligand guard ------------------------------------------------------------


def test_ligand_mpnn_refuses_a_backbone_with_no_ligand(ran, tmp_path) -> None:
    """RFdiffusion output carries no ligand, so without this the model would
    silently degrade to unconditioned design."""
    out = tmp_path / "out"
    out.mkdir()
    with pytest.raises(MissingLigandError) as exc:
        get_adapter("ligand_mpnn").run([backbone(tmp_path)], out, num_seqs=2)
    message = str(exc.value)
    assert "ligand_mpnn" in message
    assert "HETATM" in message
    assert "protein_mpnn" in message  # tells you what to use instead
    assert not ran  # refused before launching the model


def test_water_alone_is_not_ligand_context(ran, tmp_path) -> None:
    water = (
        "HETATM   10  O   HOH A 100      1.000   1.000   1.000  1.00  0.00           O\n"
    )
    out = tmp_path / "out"
    out.mkdir()
    with pytest.raises(MissingLigandError):
        get_adapter("ligand_mpnn").run([backbone(tmp_path, het=water)], out)


def test_ligand_mpnn_runs_with_a_ligand(ran, tmp_path) -> None:
    ligand = (
        "HETATM   10  C1  ATP A 200      1.000   1.000   1.000  1.00  0.00           C\n"
    )
    out = tmp_path / "out"
    out.mkdir()
    designs = get_adapter("ligand_mpnn").run(
        [backbone(tmp_path, het=ligand)], out, num_seqs=2
    )
    assert ran[0][ran[0].index("--model_type") + 1] == "ligand_mpnn"
    assert [d.producer for d in designs] == ["ligand_mpnn", "ligand_mpnn"]


def test_other_models_do_not_require_a_ligand(ran, tmp_path) -> None:
    out = tmp_path / "out"
    out.mkdir()
    for name in ("protein_mpnn", "soluble_mpnn"):
        assert get_adapter(name).run([backbone(tmp_path)], out / name, num_seqs=2)


# --- missing dependency ------------------------------------------------------


def test_missing_tool_names_the_extra_and_the_pin_conflict(monkeypatch, tmp_path) -> None:
    monkeypatch.setattr(_mpnn.shutil, "which", lambda name: None)
    with pytest.raises(MissingDependencyError) as exc:
        get_adapter("protein_mpnn").run([backbone(tmp_path)], tmp_path / "out")
    message = str(exc.value)
    assert "[mpnn]" in message
    assert "numpy==1.23.5" in message  # the conflict people will hit
    assert "dummy_sequences" in message


def test_executable_can_point_at_another_environment(monkeypatch, tmp_path) -> None:
    """The documented way around the dependency pins."""
    import dataclasses

    other = tmp_path / "env" / "ligandmpnn"
    other.parent.mkdir()
    other.write_text("#!/bin/sh\n")
    monkeypatch.setattr(_mpnn.shutil, "which", lambda name: None)

    calls: list[list[str]] = []

    def run(cmd):
        calls.append([str(c) for c in cmd])
        out = Path(cmd[cmd.index("--out_folder") + 1])
        (out / "seqs").mkdir(parents=True, exist_ok=True)
        (out / "seqs" / "bb0.fa").write_text(FASTA)

    monkeypatch.setattr(_mpnn, "stream", run)
    adapter = dataclasses.replace(get_adapter("protein_mpnn"), executable=str(other))
    out = tmp_path / "out"
    out.mkdir()
    assert adapter.run([backbone(tmp_path)], out, num_seqs=2)
    assert calls[0][0] == str(other)
