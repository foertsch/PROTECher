"""ThermoMPNN's gate check, with the subprocess faked.

The thing worth testing here is the join. ThermoMPNN answers in its own
sequence index rather than in chain positions or PDB residue numbers, so the
chain below starts at author residue 3 and has a hole where residue 5 would
be: a check that simply used `position - 1` would read the wrong rows and
still look plausible.
"""

from __future__ import annotations

from pathlib import Path

import pytest

from design_pipeline import Design, containers
from design_pipeline.adapters import _gates
from design_pipeline.adapters._tabular import write_candidates

def atom(serial: int, residue: str, chain: str, resseq: int) -> str:
    """One CA line, in the fixed columns biotite insists on."""
    return (
        f"ATOM  {serial:>5d}  CA  {residue:>3s} {chain}{resseq:>4d}    "
        f"{10.0 + serial:8.3f}{12.0:8.3f}{13.0:8.3f}  1.00 20.00"
        f"{'':10s} C\n"
    )


# Chain A holds author residues 3, 4 and 6. ThermoMPNN keys each residue by
# `author - 1`, walks that span and emits a '-' for the missing number, so its
# indices are 0, 1 and 3 for chain positions 1, 2 and 3.
PDB = atom(1, "TYR", "A", 3) + atom(2, "SER", "A", 4) + atom(3, "THR", "A", 6)

HEADER = ",ddG_pred,position,wildtype,mutation,pdb,chain,Model,Dataset\n"


@pytest.fixture
def structure(tmp_path: Path) -> Design:
    path = tmp_path / "in" / "complex.pdb"
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(PDB)
    return Design(id="complex", path=path, parent_id=None, producer="test")


@pytest.fixture
def candidates(tmp_path: Path) -> Path:
    return write_candidates(
        tmp_path / "in" / "candidates.csv",
        [
            {"variant_id": "A:Y1F", "chain": "A", "position": 1, "wt": "Y", "mut": "F"},
            {"variant_id": "A:T3A", "chain": "A", "position": 3, "wt": "T", "mut": "A"},
        ],
    )


@pytest.fixture
def installed(tmp_path: Path, monkeypatch) -> Path:
    """Both variables set and a script on disk, as a real install looks."""
    src = tmp_path / "ThermoMPNN"
    script = src / "analysis" / "custom_inference.py"
    script.parent.mkdir(parents=True, exist_ok=True)
    script.write_text("# stand-in\n")
    monkeypatch.setenv("DESIGN_PIPELINE_THERMOMPNN", "/opt/envs/thermompnn/bin/python")
    monkeypatch.setenv("THERMOMPNN_SRC", str(src))
    return script


def fake_run(monkeypatch, rows: str):
    """Write `rows` into whatever --out_dir the command names."""
    calls: list[list[str]] = []

    def run(cmd):
        command = [str(c) for c in cmd]
        calls.append(command)
        out = Path(command[command.index("--out_dir") + 1])
        out.mkdir(parents=True, exist_ok=True)
        (out / "ThermoMPNN_inference_complex.csv").write_text(HEADER + rows)

    monkeypatch.setattr(containers, "stream", run)
    return calls


def gated(produced) -> dict[str, str]:
    lines = produced[0].path.read_text().splitlines()[1:]
    return {line.split(",")[0]: line.split(",")[1] for line in lines}


def test_the_index_map_follows_gaps_and_insertion_codes() -> None:
    """Reconstructed rather than offset: an antibody chain has both, and one
    offset is wrong after the first of either."""
    assert _gates.thermompnn_index(("3", "4", "6")) == {1: 0, 2: 1, 3: 3}
    assert _gates.thermompnn_index(("100", "100A", "100B", "101")) == {
        1: 0,
        2: 1,
        3: 2,
        4: 3,
    }
    assert _gates.thermompnn_index(("1", "2", "3")) == {1: 0, 2: 1, 3: 2}


def test_a_candidate_is_joined_on_thermompnns_own_index(
    tmp_path, monkeypatch, installed, candidates, structure
) -> None:
    """Position 3 of the chain is ThermoMPNN's index 3, not its index 2. The
    destabilizing row sits at index 3, so a check that used `position - 1`
    would find the padding gap, record no value, and pass the candidate."""
    out = tmp_path / "out"
    fake_run(monkeypatch, "0,0.2,0,Y,F,complex,A,m,d\n1,2.9,3,T,A,complex,A,m,d\n")

    produced = _gates.thermompnn.run([structure], out, candidates=candidates)
    rows = gated(produced)

    assert rows["A:T3A"] == "false"
    assert rows["A:Y1F"] == "true"
    assert "ddG_pred=2.900 > 1.5" in produced[0].path.read_text()


def test_the_script_is_the_command_and_the_chain_is_a_flag(
    tmp_path, monkeypatch, installed, candidates, structure
) -> None:
    """No console script exists, so the interpreter runs the file, and the
    chain is passed because the script defaults to 'A' otherwise."""
    out = tmp_path / "out"
    calls = fake_run(
        monkeypatch, "0,0.2,0,Y,F,complex,A,m,d\n1,0.1,3,T,A,complex,A,m,d\n"
    )

    _gates.thermompnn.run([structure], out, candidates=candidates)

    assert len(calls) == 1, "one chain among the candidates, one call"
    command = calls[0]
    assert Path(command[0]) == Path("/opt/envs/thermompnn/bin/python")
    assert command[1] == str(installed)
    assert command[command.index("--chain") + 1] == "A"
    assert "--mutations_csv" not in command, "the script takes no mutation list"


def test_one_call_per_chain_the_candidates_touch(
    tmp_path, monkeypatch, installed, structure
) -> None:
    """Site saturation is per chain, and the script writes the same filename
    every time, so each chain gets its own directory."""
    both = write_candidates(
        tmp_path / "in" / "both.csv",
        [
            {"variant_id": "A:Y1F", "chain": "A", "position": 1, "wt": "Y", "mut": "F"},
            {"variant_id": "B:K1A", "chain": "B", "position": 1, "wt": "K", "mut": "A"},
        ],
    )
    structure.path.write_text(PDB + atom(4, "LYS", "B", 1))
    out = tmp_path / "out"
    calls = fake_run(monkeypatch, "0,0.2,0,Y,F,complex,A,m,d\n1,2.9,0,K,A,complex,B,m,d\n")

    produced = _gates.thermompnn.run([structure], out, candidates=both)

    assert [c[c.index("--chain") + 1] for c in calls] == ["A", "B"]
    assert (out / "chain_A").is_dir() and (out / "chain_B").is_dir()
    assert gated(produced)["B:K1A"] == "false"


def test_a_multi_mutation_candidate_gets_no_value_and_survives(
    tmp_path, monkeypatch, installed, structure
) -> None:
    """The script cannot be asked about a combination, and adding up two of its
    single-mutation numbers would be inventing an answer."""
    pair = write_candidates(
        tmp_path / "in" / "pair.csv",
        [
            {"variant_id": "A:pair", "chain": "A", "position": 1, "wt": "Y", "mut": "F"},
            {"variant_id": "A:pair", "chain": "A", "position": 3, "wt": "T", "mut": "A"},
        ],
    )
    out = tmp_path / "out"
    calls = fake_run(monkeypatch, "0,2.9,0,Y,F,complex,A,m,d\n")

    produced = _gates.thermompnn.run([structure], out, candidates=pair)
    text = produced[0].path.read_text()

    assert calls == [], "nothing to ask it about, so nothing is run"
    assert "A:pair,true" in text
    assert "no value" in text


def test_a_wild_type_disagreement_raises_rather_than_gating(
    tmp_path, monkeypatch, installed, candidates, structure
) -> None:
    """A wrong numbering would otherwise pass and fail candidates on values
    belonging to other positions, which is worse than not running."""
    out = tmp_path / "out"
    fake_run(monkeypatch, "0,2.9,0,W,F,complex,A,m,d\n")

    with pytest.raises(ValueError, match="disagrees about position 1 of chain A"):
        _gates.thermompnn.run([structure], out, candidates=candidates)


def test_the_source_folder_is_named_when_it_is_not_set(
    tmp_path, monkeypatch, candidates, structure
) -> None:
    monkeypatch.setenv("DESIGN_PIPELINE_THERMOMPNN", "/opt/envs/thermompnn/bin/python")
    monkeypatch.delenv("THERMOMPNN_SRC", raising=False)

    from design_pipeline.errors import MissingRuntimeError

    with pytest.raises(MissingRuntimeError, match="THERMOMPNN_SRC"):
        _gates.thermompnn.run([structure], tmp_path / "out", candidates=candidates)


def test_a_clone_without_the_script_is_named(
    tmp_path, monkeypatch, candidates, structure
) -> None:
    monkeypatch.setenv("DESIGN_PIPELINE_THERMOMPNN", "/opt/envs/thermompnn/bin/python")
    monkeypatch.setenv("THERMOMPNN_SRC", str(tmp_path / "not-a-clone"))

    from design_pipeline.errors import MissingRuntimeError

    with pytest.raises(MissingRuntimeError, match="custom_inference.py"):
        _gates.thermompnn.run([structure], tmp_path / "out", candidates=candidates)
