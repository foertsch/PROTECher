"""Container command assembly, with the runtime and subprocess stubbed out."""

from __future__ import annotations

from pathlib import Path

import pytest

from design_pipeline import (  # type: ignore[attr-defined]
    Design,
    MissingDependencyError,
    MissingRuntimeError,
    adapters,
)
from design_pipeline import containers, routes
from design_pipeline.adapters import _colabfold
from design_pipeline.containers import detect_runtime, format_kwargs


def design(path: Path, id: str = "parent") -> Design:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text("ATOM\n")
    return Design(id=id, path=path, parent_id=None, producer="test")


@pytest.fixture(autouse=True)
def no_runtime_override(monkeypatch):
    """Clear DESIGN_PIPELINE_RUNTIME for every test in this module.

    Autouse because it is now the first thing `detect_runtime` reads, so a
    developer with it exported would otherwise see the detection tests below
    fail for a reason that has nothing to do with the code. The tests that
    exercise the variable set it themselves.
    """
    monkeypatch.delenv(containers.RUNTIME_VAR, raising=False)


@pytest.fixture
def docker_only(monkeypatch):
    monkeypatch.setattr(
        containers.shutil,
        "which",
        lambda name: "/usr/bin/docker" if name == "docker" else None,
    )


@pytest.fixture
def both_runtimes(monkeypatch):
    monkeypatch.setattr(containers.shutil, "which", lambda name: f"/usr/bin/{name}")


def fake_stream(monkeypatch, writes):
    """Record commands and let the test pretend the model wrote something."""
    calls: list[list[str]] = []

    def run(cmd):
        calls.append([str(c) for c in cmd])
        writes(len(calls))

    monkeypatch.setattr(containers, "stream", run)
    return calls


def test_apptainer_is_preferred(both_runtimes) -> None:
    assert detect_runtime() == "apptainer"


def test_docker_is_used_when_it_is_the_only_one(docker_only) -> None:
    assert detect_runtime() == "docker"


def test_missing_runtime_says_what_to_install(monkeypatch) -> None:
    monkeypatch.setattr(containers.shutil, "which", lambda name: None)
    with pytest.raises(MissingRuntimeError) as exc:
        detect_runtime()
    message = str(exc.value)
    assert "apptainer" in message and "docker" in message and "dummy" in message


def test_the_runtime_variable_overrides_the_preference(both_runtimes, monkeypatch) -> None:
    """The point of the variable: apptainer used to win silently on a box with
    both, and docker is usually the one that works there."""
    monkeypatch.setenv(containers.RUNTIME_VAR, "docker")

    assert detect_runtime() == "docker"


def test_the_runtime_variable_is_case_insensitive(both_runtimes, monkeypatch) -> None:
    monkeypatch.setenv(containers.RUNTIME_VAR, "  Docker  ")

    assert detect_runtime() == "docker"


def test_an_unknown_runtime_names_the_two_that_work(both_runtimes, monkeypatch) -> None:
    monkeypatch.setenv(containers.RUNTIME_VAR, "podman")

    with pytest.raises(MissingRuntimeError) as exc:
        detect_runtime()

    message = str(exc.value)
    assert "podman" in message and "apptainer" in message and "docker" in message


def test_a_named_runtime_that_is_absent_does_not_fall_back(docker_only, monkeypatch) -> None:
    """Asking for apptainer and silently getting docker is not a favour: it is a
    different command line and different mount semantics."""
    monkeypatch.setenv(containers.RUNTIME_VAR, "apptainer")

    with pytest.raises(MissingRuntimeError) as exc:
        detect_runtime()

    assert "not on PATH" in str(exc.value)


def test_rfdiffusion_fresh_generation(tmp_path, monkeypatch, docker_only) -> None:
    target = tmp_path / "inputs" / "target.pdb"
    target.parent.mkdir()
    target.write_text("ATOM\n")
    out = tmp_path / "out"
    out.mkdir()
    calls = fake_stream(monkeypatch, lambda n: (out / f"design_{n}.pdb").write_text("x"))

    designs = adapters.rfdiffusion.run(
        [],
        out,
        input_pdb=target,
        num_designs=10,
        contigs="[10-40/A163-181/10-40]",
        seed=0,
        unanticipated_flag="7",
    )

    assert len(calls) == 1
    cmd = calls[0]
    assert cmd[:2] == ["docker", "run"]
    assert "--rm" in cmd and cmd[cmd.index("--gpus") + 1] == "all"
    assert routes.RFDIFFUSION_IMAGE in cmd
    assert f"{(tmp_path / 'inputs').as_posix()}:/mnt/in0" in cmd
    assert f"{out.resolve().as_posix()}:/mnt/out" in cmd
    # The target is mounted and rewritten to its container path.
    assert "inference.input_pdb=/mnt/in0/target.pdb" in cmd
    assert "inference.output_prefix=/mnt/out/design" in cmd
    assert "inference.model_directory_path=/app/RFdiffusion/models" in cmd
    assert "inference.num_designs=10" in cmd
    assert "contigmap.contigs=[10-40/A163-181/10-40]" in cmd
    assert "inference.seed=0" in cmd
    # An unknown kwarg reaches the model without a change to the package.
    assert "unanticipated_flag=7" in cmd
    assert [d.parent_id for d in designs] == [None]


def test_rfdiffusion_partial_diffusion_runs_per_design(
    tmp_path, monkeypatch, docker_only
) -> None:
    parents = [design(tmp_path / "in" / f"d{i}.pdb", id=f"d{i}") for i in range(2)]
    out = tmp_path / "out"
    out.mkdir()
    calls = fake_stream(monkeypatch, lambda n: (out / f"design_{n}.pdb").write_text("x"))

    designs = adapters.rfdiffusion.run(parents, out, partial_T=20, num_designs=1)

    assert len(calls) == 2
    # partial_T is passed raw, not converted to a noise fraction.
    assert "diffuser.partial_T=20" in calls[0]
    assert "inference.input_pdb=/mnt/in0/d0.pdb" in calls[0]
    assert "inference.input_pdb=/mnt/in0/d1.pdb" in calls[1]
    assert [d.parent_id for d in designs] == ["d0", "d1"]


def test_apptainer_flags(tmp_path, monkeypatch, both_runtimes) -> None:
    parent = design(tmp_path / "in" / "d.pdb")
    out = tmp_path / "out"
    out.mkdir()
    calls = fake_stream(monkeypatch, lambda n: (out / f"design_{n}.pdb").write_text("x"))

    adapters.rfdiffusion.run([parent], out, num_designs=1)

    cmd = calls[0]
    assert cmd[:2] == ["apptainer", "run"]
    assert "--nv" in cmd
    assert "--bind" in cmd
    assert f"docker://{routes.RFDIFFUSION_IMAGE}" in cmd


def test_hydra_and_flag_styles() -> None:
    assert format_kwargs({"a": 1, "b": None}, "hydra", {}) == ["a=1"]
    assert format_kwargs({"a": 1}, "flag", {"a": "alpha"}) == ["--alpha", "1"]
    assert format_kwargs({"flagged": True}, "flag", {}) == ["--flagged"]
    assert format_kwargs({"host_url": "u"}, "flag", {}, dash_flags=True) == [
        "--host-url",
        "u",
    ]
    assert format_kwargs({"hotspot_res": ["A59", "A83"]}, "hydra", {}) == [
        "hotspot_res=[A59,A83]"
    ]


def test_silent_failure_is_an_error(tmp_path, monkeypatch, both_runtimes) -> None:
    parent = design(tmp_path / "in" / "d.pdb")
    out = tmp_path / "out"
    out.mkdir()
    fake_stream(monkeypatch, lambda n: None)

    with pytest.raises(containers.RunFailedError):
        adapters.rfdiffusion.run([parent], out, num_designs=1)


def test_colabfold_missing_names_the_extra(tmp_path, monkeypatch) -> None:
    monkeypatch.setattr(_colabfold.shutil, "which", lambda name: None)
    fasta = design(tmp_path / "s.fa")
    with pytest.raises(MissingDependencyError) as exc:
        _colabfold.ColabFoldAdapter().run([fasta], tmp_path / "out")
    assert "design-pipeline[colabfold]" in str(exc.value)
    assert "dummy_structures" in str(exc.value)


def test_colabfold_command(tmp_path, monkeypatch) -> None:
    monkeypatch.setattr(
        _colabfold.shutil, "which", lambda name: "/usr/bin/colabfold_batch"
    )
    calls: list[list[str]] = []
    out = tmp_path / "out"
    out.mkdir()

    def run(cmd):
        calls.append([str(c) for c in cmd])
        (out / "s_unrelaxed_rank_001_model_1.pdb").write_text("ATOM\n")
        (out / "s_scores_rank_001_model_1.json").write_text('{"plddt": [80, 90]}')

    monkeypatch.setattr(_colabfold, "stream", run)
    fasta = design(tmp_path / "s.fa", id="seq1")

    designs = _colabfold.ColabFoldAdapter().run([fasta], out, recycles=3, seed=7)

    cmd = calls[0]
    assert cmd[0].endswith("colabfold_batch")
    # Public MSA server by default; host_url lets a local one be pointed at.
    assert cmd[cmd.index("--host-url") + 1] == "https://api.colabfold.com"
    assert cmd[cmd.index("--num-recycle") + 1] == "3"
    assert cmd[cmd.index("--random-seed") + 1] == "7"
    assert len(designs) == 1
    assert designs[0].parent_id == "seq1"
    assert designs[0].scores["plddt"] == 85.0
    assert designs[0].meta["scores_json"].endswith("s_scores_rank_001_model_1.json")


def test_colabfold_accepts_a_local_msa_server(tmp_path, monkeypatch) -> None:
    monkeypatch.setattr(
        _colabfold.shutil, "which", lambda name: "/usr/bin/colabfold_batch"
    )
    calls: list[list[str]] = []
    out = tmp_path / "out"
    out.mkdir()

    def run(cmd):
        calls.append([str(c) for c in cmd])
        (out / "s_unrelaxed_rank_001.pdb").write_text("ATOM\n")

    monkeypatch.setattr(_colabfold, "stream", run)
    _colabfold.ColabFoldAdapter().run(
        [design(tmp_path / "s.fa")], out, host_url="http://localhost:8080"
    )
    assert calls[0][calls[0].index("--host-url") + 1] == "http://localhost:8080"
