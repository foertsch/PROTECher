"""antifold_seq: AntiFold's whole-sequence likelihood as a scorer, one score per variant_id.

AntiFold itself is not imported here (it needs torch, Python 3.13 and its own environment). `FakeComplex` stands in
for `score_sequence_in_complex` with the same contract, a chain sequence in and an average log-likelihood out, and is
deliberately not additive within a chain, because that is the property this scorer exists for. The script is also run
as a subprocess against a fake `antifold` package, which is what proves it works by path with nothing of ours imported.
The real function was run on the cluster, on 1N8Z and on two public sets of measured variants.
"""

from __future__ import annotations

import ast
import os
import subprocess
import sys
import textwrap
import types
from pathlib import Path

import numpy as np
import pytest

from design_pipeline import Design, containers, routes
from design_pipeline.adapters import _antifold_seq as seq
from design_pipeline.adapters import _scorers
from design_pipeline.adapters._tabular import SCORES, read_tool_table, write_candidates
from design_pipeline.errors import MissingRuntimeError
from design_pipeline.registry import REGISTRY

HEAVY, LIGHT = "EVQLK", "DIQM"


class FakeComplex:
    """Stands in for ComplexScorer. A change at position i costs 0.5 + 0.1 i nats, and every pair of changes in one
    chain costs 0.3 more, so several mutations are worse than the sum of the singles."""

    def __init__(self):
        self.sequences = {"B": HEAVY, "A": LIGHT}
        self.calls: list[tuple[str, str]] = []

    def log_likelihood(self, chain: str, sequence: str) -> float:
        self.calls.append((chain, sequence))
        changed = [i for i, (a, b) in enumerate(zip(self.sequences[chain], sequence)) if a != b]
        total = -sum(0.5 + 0.1 * i for i in changed) - 0.3 * len(changed) * (len(changed) - 1) / 2
        return total / len(sequence)


def rows_of(candidates):
    fake = FakeComplex()
    return fake, {row[0]: row for row in seq.score_with_complex(fake, candidates)}


# ---- scoring -----------------------------------------------------------------------------------------------------


def test_a_single_mutation_is_the_change_in_total_log_likelihood_of_its_chain():
    _, rows = rows_of({"v": [("B", 3, "Q", "E")]})
    _, score, raw, note = rows["v"]
    assert note == "" and score == pytest.approx(-(0.5 + 0.1 * 2))
    assert raw == pytest.approx(score)  # the fake's wild type has log-likelihood 0, so the mutant's total is the change


def test_several_mutations_in_one_chain_are_scored_as_a_whole_not_as_a_sum_of_singles():
    _, rows = rows_of({"a": [("B", 1, "E", "D")], "b": [("B", 3, "Q", "E")], "both": [("B", 1, "E", "D"), ("B", 3, "Q", "E")]})
    assert rows["both"][1] == pytest.approx(-(0.5 + 0.7) - 0.3)
    assert rows["both"][1] < rows["a"][1] + rows["b"][1] - 0.2


def test_mutations_in_both_chains_add_the_two_chain_scores():
    _, rows = rows_of({"h": [("B", 2, "V", "L")], "l": [("A", 2, "I", "V")], "hl": [("B", 2, "V", "L"), ("A", 2, "I", "V")]})
    assert rows["hl"][1] == pytest.approx(rows["h"][1] + rows["l"][1])


def test_the_order_of_the_mutations_in_a_candidate_does_not_matter():
    _, rows = rows_of({"x": [("B", 1, "E", "D"), ("B", 4, "L", "I")], "y": [("B", 4, "L", "I"), ("B", 1, "E", "D")]})
    assert rows["x"][1] == pytest.approx(rows["y"][1])


def test_one_pass_per_mutated_chain_and_the_same_mutant_chain_is_scored_once():
    fake, rows = rows_of({
        "x": [("B", 1, "E", "D")],
        "same": [("B", 1, "E", "D")],  # the identical mutant chain again
        "hl": [("B", 2, "V", "L"), ("A", 2, "I", "V")],
    })
    assert fake.calls[:2] == [("B", HEAVY), ("A", LIGHT)]  # wild types first, once per chain
    assert len(fake.calls) == 2 + 1 + 2  # two wild types, one pass shared by x and same, one per chain for hl
    assert rows["x"][1] == rows["same"][1]


@pytest.mark.parametrize(
    "mutations, expected",
    [
        ([("B", 2, "A", "G")], "wild type A at B2, the structure has V"),
        ([("B", 9, "E", "G")], "position 9 is outside chain B, which has 5 residues"),
        ([("B", 0, "E", "G")], "position 0 is outside chain B"),
        ([("H", 1, "E", "G")], "chain 'H' is not in the structure"),
        ([("B", 1, "E", "X")], "'X' is not one of the 20 amino acids"),
        ([("B", 1, "E", "D"), ("B", 1, "E", "G")], "B1 is mutated twice in one candidate"),
    ],
)
def test_a_candidate_that_does_not_fit_gets_no_score_and_a_reason_and_the_others_are_scored(mutations, expected):
    _, rows = rows_of({"v": mutations, "ok": [("B", 1, "E", "D")]})
    assert rows["v"][1] is None and rows["v"][2] is None and expected in rows["v"][3]
    assert rows["ok"][1] is not None


def test_check_mutations_returns_an_empty_string_for_a_fitting_candidate():
    assert seq.check_mutations({"B": HEAVY}, [("B", 1, "E", "D"), ("B", 2, "V", "L")]) == ""


# ---- the wrapper around AntiFold, with stand-in modules --------------------------------------------------------


@pytest.fixture
def fake_antifold(monkeypatch):
    seen = {}

    def load_model(checkpoint_path=None, device="cuda"):
        seen["load_model"] = (checkpoint_path, device)
        return "model"

    def load_complex_coords(path, chains):
        seen["load_complex_coords"] = (path, chains)
        return {np.str_(c): np.zeros((3, 3, 3)) for c in chains}, {np.str_("B"): np.str_(HEAVY), np.str_("A"): np.str_(LIGHT)}, {}, {}

    def score_sequence_in_complex(model, alphabet, coords, chain, sequence):
        seen["score"] = (model, alphabet, sorted(coords), chain, sequence)
        return -1.5, -9.9  # (over the full chain, without residues lacking coordinates)

    modules = {
        "antifold": types.ModuleType("antifold"),
        "antifold.esm": types.ModuleType("antifold.esm"),
        "antifold.antiscripts": types.SimpleNamespace(load_model=load_model),
        "antifold.esm_multichain_util_custom": types.SimpleNamespace(
            load_complex_coords=load_complex_coords, score_sequence_in_complex=score_sequence_in_complex
        ),
    }
    modules["antifold.esm"].data = types.SimpleNamespace(Alphabet=types.SimpleNamespace(from_architecture=lambda name: f"alphabet:{name}"))
    modules["antifold"].esm = modules["antifold.esm"]  # `import antifold.esm` only finds this if the parent carries it
    for name, module in modules.items():
        monkeypatch.setitem(sys.modules, name, module)
    return seen


def test_complex_scorer_passes_weights_device_and_chain_order_to_antifold(fake_antifold, tmp_path):
    scorer = seq.ComplexScorer(tmp_path / "1N8Z.pdb", ["B", "A"], device="cpu", weights=tmp_path / "model.pt")
    assert fake_antifold["load_model"] == (str(tmp_path / "model.pt"), "cpu")
    assert fake_antifold["load_complex_coords"] == (str(tmp_path / "1N8Z.pdb"), ["B", "A"])
    assert scorer.sequences == {"B": HEAVY, "A": LIGHT} and all(type(c) is str for c in scorer.sequences)
    assert scorer.log_likelihood("B", HEAVY) == -1.5  # the full-chain average, not the one without missing coordinates
    assert fake_antifold["score"] == ("model", "alphabet:invariant_gvp", ["A", "B"], "B", HEAVY)


def test_complex_scorer_uses_antifolds_default_weights_when_none_are_given(fake_antifold, tmp_path):
    seq.ComplexScorer(tmp_path / "x.pdb", ["B", "A"])
    assert fake_antifold["load_model"] == (None, "cpu")


def test_complex_scorer_says_what_to_do_when_antifold_is_not_importable(monkeypatch, tmp_path):
    monkeypatch.setitem(sys.modules, "antifold", None)  # an import of it now raises ModuleNotFoundError
    with pytest.raises(seq.ScoreError, match="0.4.0 or later.*DESIGN_PIPELINE_ANTIFOLD"):
        seq.ComplexScorer(tmp_path / "x.pdb", ["B", "A"])


# ---- the command line --------------------------------------------------------------------------------------------


def test_main_writes_the_scores_scores_pairs_and_keeps_the_gpus_away(tmp_path, monkeypatch, capsys):
    made = {}

    class Stub(FakeComplex):
        def __init__(self, pdb, chains, device, weights):
            super().__init__()
            made.update(pdb=pdb, chains=chains, device=device, weights=weights, cuda=os.environ.get("CUDA_VISIBLE_DEVICES"))

    monkeypatch.setattr(seq, "ComplexScorer", Stub)
    monkeypatch.setenv("CUDA_VISIBLE_DEVICES", "0")
    write_candidates(tmp_path / "c.csv", [
        {"variant_id": "pair", "chain": "B", "position": 1, "wt": "E", "mut": "D"},
        {"variant_id": "pair", "chain": "B", "position": 3, "wt": "Q", "mut": "E"},
        {"variant_id": "bad", "chain": "B", "position": 2, "wt": "A", "mut": "G"},
    ])
    args = ["--pdb", "x.pdb", "--heavy", "B", "--light", "A", "--out_dir", str(tmp_path / "out"), "--mutations_csv", str(tmp_path / "c.csv"), "--weights", "w.pt"]
    assert seq.main(args) == 0
    assert made == {"pdb": Path("x.pdb"), "chains": ["B", "A"], "device": "cpu", "weights": Path("w.pt"), "cuda": ""}
    out = tmp_path / "out" / seq.OUTPUT
    assert read_tool_table([out], "variant_id", "score") == {"pair": pytest.approx(-(0.5 + 0.7) - 0.3)}  # the pair is scored, the bad row is not
    assert out.read_text().splitlines()[0] == "variant_id,score,raw,note"
    printed = capsys.readouterr().out
    assert "1 of 2 candidates scored" in printed and "no score for bad: wild type A at B2" in printed


@pytest.mark.parametrize("extra", [[], ["--heavy", "B"], ["--light", "A"]])
def test_main_needs_both_chain_ids(tmp_path, extra):
    with pytest.raises(SystemExit):
        seq.main(["--pdb", "x.pdb", "--out_dir", str(tmp_path), "--mutations_csv", str(tmp_path / "c.csv"), *extra])


FAKE_ANTIFOLD = {
    "antifold/__init__.py": "",
    "antifold/esm/__init__.py": textwrap.dedent("""
        import types
        data = types.SimpleNamespace(Alphabet=types.SimpleNamespace(from_architecture=lambda name: name))
    """),
    "antifold/antiscripts.py": textwrap.dedent("""
        import os
        def load_model(checkpoint_path=None, device="cuda"):
            open(os.environ["FAKE_SEEN"], "a").write(f"load_model {checkpoint_path} {device} cuda={os.environ.get('CUDA_VISIBLE_DEVICES')!r}\\n")
            return "model"
    """),
    "antifold/esm_multichain_util_custom.py": textwrap.dedent(f"""
        import numpy as np
        BASE = {{"B": {HEAVY!r}, "A": {LIGHT!r}}}
        def load_complex_coords(path, chains):
            return {{c: np.zeros((len(BASE[c]), 3, 3)) for c in chains}}, dict(BASE), {{}}, {{}}
        def score_sequence_in_complex(model, alphabet, coords, chain, sequence):
            changed = [i for i, (a, b) in enumerate(zip(BASE[chain], sequence)) if a != b]
            total = -sum(0.5 + 0.1 * i for i in changed) - 0.3 * len(changed) * (len(changed) - 1) / 2
            return total / len(sequence), 0.0
    """),
}


def test_the_script_runs_by_path_with_only_antifold_importable(tmp_path):
    """The adapter starts it with AntiFold's own interpreter, where this package is not installed."""
    fake = tmp_path / "fake_env"
    for name, text in FAKE_ANTIFOLD.items():
        (fake / name).parent.mkdir(parents=True, exist_ok=True)
        (fake / name).write_text(text)
    write_candidates(tmp_path / "c.csv", [
        {"variant_id": "single", "chain": "B", "position": 3, "wt": "Q", "mut": "E"},
        {"variant_id": "pair", "chain": "B", "position": 1, "wt": "E", "mut": "D"},
        {"variant_id": "pair", "chain": "A", "position": 2, "wt": "I", "mut": "V"},
    ])
    seen = tmp_path / "seen.txt"
    done = subprocess.run(
        [sys.executable, str(Path(seq.__file__)), "--pdb", str(tmp_path / "x.pdb"), "--heavy", "B", "--light", "A",
         "--out_dir", str(tmp_path / "out"), "--mutations_csv", str(tmp_path / "c.csv"), "--weights", "w.pt"],
        env={**os.environ, "PYTHONPATH": str(fake), "FAKE_SEEN": str(seen)}, capture_output=True, text=True, cwd=tmp_path,
    )
    assert done.returncode == 0, done.stderr
    scores = read_tool_table([tmp_path / "out" / seq.OUTPUT], "variant_id", "score")
    assert scores["single"] == pytest.approx(-0.7) and scores["pair"] == pytest.approx(-0.5 + -0.6)
    assert seen.read_text().strip() == "load_model w.pt cpu cuda=''"


def test_the_script_says_so_and_exits_when_antifold_is_missing(tmp_path):
    write_candidates(tmp_path / "c.csv", [{"variant_id": "v", "chain": "B", "position": 1, "wt": "E", "mut": "D"}])
    done = subprocess.run(
        [sys.executable, "-S", str(Path(seq.__file__)), "--pdb", "x.pdb", "--heavy", "B", "--light", "A",
         "--out_dir", str(tmp_path / "out"), "--mutations_csv", str(tmp_path / "c.csv")],
        env={"PATH": os.environ.get("PATH", "")}, capture_output=True, text=True, cwd=tmp_path,  # -S: no site-packages
    )
    assert done.returncode != 0 and "antifold_seq: AntiFold 0.4.0 or later is not importable" in done.stderr


def test_the_script_imports_nothing_from_this_package():
    tree = ast.parse(Path(seq.__file__).read_text(encoding="utf-8"))
    imported = {n.module.split(".")[0] for n in ast.walk(tree) if isinstance(n, ast.ImportFrom) and n.module and n.level == 0}
    imported |= {a.name.split(".")[0] for n in ast.walk(tree) if isinstance(n, ast.Import) for a in n.names}
    assert not any(name.startswith("design_pipeline") for name in imported)
    assert not any(isinstance(n, ast.ImportFrom) and n.level > 0 for n in ast.walk(tree))  # no relative imports either


# ---- the adapter -------------------------------------------------------------------------------------------------


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


def fake_stream(monkeypatch, out: Path):
    calls: list[list[str]] = []

    def run(cmd):
        calls.append([str(c) for c in cmd])
        (out / seq.OUTPUT).write_text("variant_id,score,raw,note\nH:T30Y,-1.5,-10.0,\nH:S35A,,,wild type S at H35, the structure has T\n")

    monkeypatch.setattr(containers, "stream", run)
    return calls


def test_the_adapter_starts_the_script_by_path_with_the_antifold_interpreter(tmp_path, monkeypatch, candidates, structure):
    out = tmp_path / "out"
    out.mkdir()
    calls = fake_stream(monkeypatch, out)
    monkeypatch.setenv("DESIGN_PIPELINE_ANTIFOLD", "/envs/antifold/bin/python")

    _scorers.antifold_seq.run([structure], out, candidates=candidates, heavy="H", light="L", weights="/w/model.pt")
    command = calls[0]

    assert command[0] == "/envs/antifold/bin/python"  # reached by path, like the matrix scorer
    assert command[1] == str(seq.__file__) and Path(command[1]).is_file()
    for flag, value in [("--pdb", str(structure.path)), ("--out_dir", str(out)), ("--mutations_csv", str(candidates)), ("--heavy", "H"), ("--light", "L"), ("--weights", "/w/model.pt")]:
        assert command[command.index(flag) + 1] == value


def test_an_explicit_executable_wins_over_the_environment_variable(tmp_path, monkeypatch, candidates, structure):
    out = tmp_path / "out"
    out.mkdir()
    calls = fake_stream(monkeypatch, out)
    monkeypatch.setenv("DESIGN_PIPELINE_ANTIFOLD", "/envs/antifold/bin/python")

    _scorers.antifold_seq.run([structure], out, candidates=candidates, heavy="H", light="L", executable="/other/bin/python")

    assert calls[0][0] == "/other/bin/python"


def test_the_tools_file_and_the_normalised_scores_do_not_collide(tmp_path, monkeypatch, candidates, structure):
    out = tmp_path / "out"
    out.mkdir()
    fake_stream(monkeypatch, out)
    monkeypatch.setenv("DESIGN_PIPELINE_ANTIFOLD", "/envs/antifold/bin/python")

    produced = _scorers.antifold_seq.run([structure], out, candidates=candidates, heavy="H", light="L")
    text = produced[0].path.read_text()

    assert produced[0].path.name == SCORES and (out / seq.OUTPUT).is_file() and seq.OUTPUT != SCORES
    assert text.splitlines()[0] == "variant_id,score,raw"
    assert "H:T30Y,-1.5,-1.5" in text  # sign +1: score and raw agree
    assert "H:S35A,," in text  # a candidate the tool could not score is no score, not a zero


def test_without_a_configured_environment_the_error_carries_the_route(tmp_path, monkeypatch, candidates, structure):
    monkeypatch.delenv("DESIGN_PIPELINE_ANTIFOLD", raising=False)
    monkeypatch.setitem(sys.modules, "antifold", None)  # not importable here either

    with pytest.raises(MissingRuntimeError, match="antifold_seq is not installed.*0.4.0 or later"):
        _scorers.antifold_seq.run([structure], tmp_path / "out", candidates=candidates, heavy="H", light="L")


# ---- registration ------------------------------------------------------------------------------------------------


def test_it_is_registered_next_to_antifold_but_is_not_a_default_scorer():
    assert REGISTRY["antifold_seq"] is _scorers.antifold_seq and "antifold" in REGISTRY
    assert "antifold_seq" not in _scorers.DEFAULT_SCORERS  # the default campaign is the matrix scorer's


def test_it_is_per_variant_so_the_uncapped_table_guard_applies():
    assert _scorers.antifold_seq.per_variant is True and _scorers.antifold.per_variant is False


def test_its_route_reads_the_same_variable_as_the_matrix_scorer():
    route = routes.ROUTES["antifold_seq"]
    assert route.env == routes.ENV_B and route.verified and route.module == "antifold"
    assert route.env_var == routes.ROUTES["antifold"].env_var == _scorers.antifold_seq_process.env_var
    assert "0.4.0 or later" in route.how
