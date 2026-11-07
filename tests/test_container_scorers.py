"""The containerised scorers and gate checks, with the runtime faked.

What matters here is the boundary: the candidates CSV reaches the container as
a mounted path, the tool's own column names are read back, and the sign is
normalized so `score` always means better binding.
"""

from __future__ import annotations

import subprocess
import sys
from dataclasses import replace
from pathlib import Path

import pytest

from design_pipeline import Design, containers
from design_pipeline.adapters import _gates, _scorers
from design_pipeline.adapters._tabular import (
    SCORES,
    write_candidates,
)


def scorer(adapter, image="test/scorer", *, gpu=True):
    """A shipped scorer with an image attached, for a boundary test.

    The registry ships no image names any more - every containerised scorer has
    `tool=None` until you build one and set it - so a test of the container
    boundary has to supply the image the way a user does.
    """
    return replace(adapter, tool=_scorers.matrix_container(image, gpu=gpu))


def check(adapter, image="test/check", *, gpu=True):
    """The same, for a gate check."""
    return replace(adapter, tool=_gates.check_container(image, gpu=gpu))


@pytest.fixture
def docker_only(monkeypatch):
    monkeypatch.setattr(
        containers.shutil,
        "which",
        lambda name: "/usr/bin/docker" if name == "docker" else None,
    )


@pytest.fixture
def candidates(tmp_path: Path) -> Path:
    return write_candidates(
        tmp_path / "in" / "candidates.csv",
        [
            {"variant_id": "H:T30Y", "chain": "H", "position": 30, "wt": "T", "mut": "Y"},
            {"variant_id": "H:S35A", "chain": "H", "position": 35, "wt": "S", "mut": "A"},
        ],
    )


@pytest.fixture
def structure(tmp_path: Path) -> Design:
    path = tmp_path / "in" / "complex.pdb"
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text("ATOM\n")
    return Design(id="complex", path=path, parent_id=None, producer="test")


def fake_stream(monkeypatch, writes):
    calls: list[list[str]] = []

    def run(cmd):
        calls.append([str(c) for c in cmd])
        writes(len(calls))

    monkeypatch.setattr(containers, "stream", run)
    return calls


def test_a_ddg_scorer_is_flipped_so_higher_is_better(
    tmp_path, monkeypatch, docker_only, candidates, structure
) -> None:
    """GearBind predicts a ddG of binding, where the improvement is negative.
    The table never shows that: `score` is normalized and `raw` keeps it."""
    out = tmp_path / "out"
    out.mkdir()
    calls = fake_stream(
        monkeypatch,
        lambda n: (out / "gearbind.csv").write_text(
            "variant_id,ddg\nH:T30Y,-1.5\nH:S35A,0.75\n"
        ),
    )

    produced = scorer(_scorers.gearbind).run([structure], out, candidates=candidates)
    rows = dict(
        line.split(",")[:2]
        for line in produced[0].path.read_text().splitlines()[1:]
    )
    raws = {
        line.split(",")[0]: line.split(",")[2]
        for line in produced[0].path.read_text().splitlines()[1:]
    }

    assert produced[0].path.name == SCORES
    assert float(rows["H:T30Y"]) == 1.5  # negative ddG, better binding
    assert float(rows["H:S35A"]) == -0.75
    assert float(raws["H:T30Y"]) == -1.5  # the tool's own number, kept
    assert len(calls) == 1


def test_the_candidates_csv_reaches_the_container_as_a_mounted_path(
    tmp_path, monkeypatch, docker_only, candidates, structure
) -> None:
    out = tmp_path / "out"
    out.mkdir()
    calls = fake_stream(
        monkeypatch,
        lambda n: (out / "antifold.csv").write_text("variant_id,score\nH:T30Y,0.2\n"),
    )

    scorer(_scorers.antifold).run([structure], out, candidates=candidates)
    command = calls[0]

    assert "docker" in command[0]
    assert "--gpus" in command  # a matrix scorer on the GPU
    # The flag the adapter maps `candidates` onto, pointing inside the container.
    index = command.index("--mutations_csv")
    assert command[index + 1].startswith("/mnt/")
    assert command[index + 1].endswith("candidates.csv")
    assert str(tmp_path) not in " ".join(command[index : index + 2])


def test_a_candidate_the_tool_skipped_is_no_score_rather_than_a_zero(
    tmp_path, monkeypatch, docker_only, candidates, structure
) -> None:
    out = tmp_path / "out"
    out.mkdir()
    fake_stream(
        monkeypatch,
        lambda n: (out / "antifold.csv").write_text("variant_id,score\nH:T30Y,0.2\n"),
    )

    produced = scorer(_scorers.antifold).run([structure], out, candidates=candidates)
    lines = produced[0].path.read_text().splitlines()

    assert lines[0] == "variant_id,score,raw"
    assert "H:S35A,," in produced[0].path.read_text()
    assert len(lines) == 3  # both candidates get a row


def test_the_tools_column_names_are_overridable_per_call(
    tmp_path, monkeypatch, docker_only, candidates, structure
) -> None:
    """They move between versions, so a rename must not need a patch to this
    package."""
    out = tmp_path / "out"
    out.mkdir()
    fake_stream(
        monkeypatch,
        lambda n: (out / "antifold.csv").write_text("id,llr\nH:T30Y,0.2\n"),
    )

    produced = scorer(_scorers.antifold).run(
        [structure], out, candidates=candidates, id_column="id", value_column="llr"
    )

    assert "H:T30Y,0.2" in produced[0].path.read_text()


def test_an_unreadable_output_says_which_columns_were_expected(
    tmp_path, monkeypatch, docker_only, candidates, structure
) -> None:
    out = tmp_path / "out"
    out.mkdir()
    fake_stream(
        monkeypatch,
        lambda n: (out / "antifold.csv").write_text("something,else\n1,2\n"),
    )

    with pytest.raises(ValueError, match="id_column"):
        scorer(_scorers.antifold).run([structure], out, candidates=candidates)


def test_a_gate_container_fails_only_on_the_wrong_side_of_its_threshold(
    tmp_path, monkeypatch, docker_only, candidates, structure
) -> None:
    """A destabilizing prediction past the threshold is a confident failure;
    a marginal one is not, which is why the shipped number is loose."""
    out = tmp_path / "out"
    out.mkdir()
    fake_stream(
        monkeypatch,
        lambda n: (out / "aggrescan3d.csv").write_text(
            "variant_id,a3d_score\nH:T30Y,2.5\nH:S35A,0.2\n"
        ),
    )

    produced = check(_gates.aggrescan3d).run([structure], out, candidates=candidates)
    rows = {
        line.split(",")[0]: line.split(",")[1]
        for line in produced[0].path.read_text().splitlines()[1:]
    }

    assert rows["H:T30Y"] == "false"
    assert rows["H:S35A"] == "true"


def test_a_solubility_check_fails_below_its_threshold(
    tmp_path, monkeypatch, docker_only, candidates, structure
) -> None:
    out = tmp_path / "out"
    out.mkdir()
    fake_stream(
        monkeypatch,
        lambda n: (out / "netsolp.csv").write_text(
            "variant_id,solubility\nH:T30Y,0.1\nH:S35A,0.9\n"
        ),
    )

    produced = check(_gates.netsolp).run([structure], out, candidates=candidates)
    text = produced[0].path.read_text()

    assert "H:T30Y,false" in text
    assert "H:S35A,true" in text


def test_a_check_with_no_answer_does_not_drop_the_candidate(
    tmp_path, monkeypatch, docker_only, candidates, structure
) -> None:
    """Gate only on confident failures: no value is not a failure."""
    out = tmp_path / "out"
    out.mkdir()
    fake_stream(
        monkeypatch,
        lambda n: (out / "netsolp.csv").write_text(
            "variant_id,solubility\nH:T30Y,0.9\n"
        ),
    )

    produced = check(_gates.netsolp).run([structure], out, candidates=candidates)
    text = produced[0].path.read_text()

    assert "H:S35A,true" in text
    assert "no value" in text


def test_the_threshold_is_overridable_per_call(
    tmp_path, monkeypatch, docker_only, candidates, structure
) -> None:
    """The shipped numbers are first guesses; the harness is what settles
    them."""
    out = tmp_path / "out"
    out.mkdir()
    fake_stream(
        monkeypatch,
        lambda n: (out / "aggrescan3d.csv").write_text(
            "variant_id,a3d_score\nH:T30Y,2.0\nH:S35A,0.5\n"
        ),
    )

    strict = check(_gates.aggrescan3d).run(
        [structure], out, candidates=candidates, threshold=0.1
    )

    assert strict[0].path.read_text().count(",false,") == 2


def test_every_scorer_and_check_declares_its_cost_class() -> None:
    """`score_variants` reads this to decide what to cap, so a new adapter that
    forgets it would be treated as cheap."""
    from design_pipeline import REGISTRY
    from design_pipeline.adapters._scorers import MATRIX_SCORERS, PER_VARIANT_SCORERS

    for name in MATRIX_SCORERS:
        assert getattr(REGISTRY[name], "per_variant") is False, name
    for name in PER_VARIANT_SCORERS:
        assert getattr(REGISTRY[name], "per_variant") is True, name


def test_antifold_is_run_through_antifold_score_and_never_gets_a_bare_candidates_flag(
    tmp_path, candidates, structure
) -> None:
    """AntiFold's own command line has no `--candidates`; the adapter used to
    forward one and the run died with an argparse error."""
    command = _scorers.antifold_process._command(
        "af-python", [structure], tmp_path / "out", (), {"candidates": str(candidates)}
    )

    assert command[:2] == ["af-python", "-c"]
    assert "antifold_score" in command[2]
    assert "--candidates" not in command
    assert command[command.index("--mutations_csv") + 1] == str(candidates)
    assert command[command.index("--route") + 1] == "table"
    assert command[command.index("--pdb") + 1] == str(structure.path)


def test_antifold_scores_come_back_through_the_scorer_contract(
    tmp_path, monkeypatch, candidates, structure
) -> None:
    out = tmp_path / "out"
    out.mkdir()
    calls = fake_stream(
        monkeypatch,
        lambda n: (out / SCORES).write_text(
            "variant_id,score,raw,note\nH:T30Y,0.2,0.1,\nH:S35A,,,2 mutations\n"
        ),
    )

    produced = _scorers.antifold.run(
        [structure], out, candidates=candidates, executable="af-python"
    )
    rows = produced[0].path.read_text().splitlines()

    assert len(calls) == 1 and calls[0][0] == "af-python"
    assert rows[1].split(",")[:2] == ["H:T30Y", "0.2"]
    assert not any(line.startswith("H:S35A") and line.split(",")[1] for line in rows[1:])


def test_the_antifold_score_command_is_found_from_any_working_directory(tmp_path) -> None:
    """No PYTHONPATH and a working directory that is not the repo: the `-c`
    bootstrap has to put the repo root on the path itself."""
    import os

    env = {k: v for k, v in os.environ.items() if k != "PYTHONPATH"}
    done = subprocess.run(
        [sys.executable, "-c", _scorers._ANTIFOLD_SCORE, "--help"],
        cwd=tmp_path,
        env=env,
        capture_output=True,
        text=True,
    )

    assert done.returncode == 0, done.stderr
    assert "--mutations_csv" in done.stdout
