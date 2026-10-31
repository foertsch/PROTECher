"""antifold_score: candidates scored with AntiFold's own functions, one score per variant_id.

AntiFold itself is not imported here (it needs torch and its own Python). Two stand-ins carry the same contracts
as the real functions, so these tests pin what the glue does with them:

* `fake_scorer` for `get_sequence_sampled_global_score`: a sequence over the table's rows in, the mean negative
  log-likelihood over the rows out (twice, as the real one returns two numbers).
* `FakeComplex` for `score_sequence_in_complex`: a chain sequence in, an average log-likelihood out, deliberately
  not additive within a chain, because that is the property the complex route exists for.

The real functions are exercised on the cluster, in the AntiFold environment.
"""

import os
import sys
import textwrap
import types
from pathlib import Path

import numpy as np
import pandas as pd
import pytest

from antifold_score import ComplexScorer, ScoreError
from antifold_score import scorer as scorer_module
from antifold_score.scorer import (
    AMINO, SCORES, check_mutations, main, read_candidates, run_antifold, score_with_complex, score_with_table, write_scores,
)
from design_pipeline.adapters import _tabular

HEAVY, LIGHT = "EVQLK", "DIQM"


@pytest.fixture
def table() -> pd.DataFrame:
    rng = np.random.default_rng(0)
    logits = rng.normal(size=(len(HEAVY) + len(LIGHT), len(AMINO)))
    log_probs = logits - np.log(np.exp(logits).sum(axis=1, keepdims=True))
    frame = pd.DataFrame(log_probs, columns=list(AMINO))
    frame.insert(0, "pdb_res", list(HEAVY + LIGHT))
    frame.insert(0, "pdb_chain", ["B"] * len(HEAVY) + ["A"] * len(LIGHT))
    return frame


def fake_scorer(sequence: str, table: pd.DataFrame) -> tuple[float, float]:
    values = table[list(AMINO)].values
    mean = -float(np.mean([values[i, AMINO.index(aa)] for i, aa in enumerate(sequence)]))
    return mean, mean


def lookup(table, row: int, mut: str) -> float:
    """log P(mut) - log P(wt) at one row, straight from the table."""
    values = table[list(AMINO)].values
    return values[row, AMINO.index(mut)] - values[row, AMINO.index(table["pdb_res"][row])]


def one(table, candidates):
    (row,) = score_with_table(table, candidates, scorer=fake_scorer)
    return row


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


# ---- table route -------------------------------------------------------------------------------------------------


def test_a_single_mutation_scores_as_the_log_likelihood_ratio_in_nats(table):
    variant_id, score, raw, note = one(table, {"v": [("B", 3, "Q", "E")]})
    assert note == ""
    assert score == pytest.approx(lookup(table, 2, "E"))
    assert raw == pytest.approx(fake_scorer(HEAVY[:2] + "E" + HEAVY[3:] + LIGHT, table)[0])


def test_positions_count_within_each_chain_from_one(table):
    assert one(table, {"v": [("A", 1, "D", "E")]})[1] == pytest.approx(lookup(table, len(HEAVY), "E"))
    assert one(table, {"v": [("B", 5, "K", "R")]})[1] == pytest.approx(lookup(table, 4, "R"))


def test_higher_is_better_means_the_table_likes_the_mutant_more(table):
    others = [aa for aa in AMINO if aa != "E"]
    best = max(others, key=lambda aa: lookup(table, 0, aa))
    worst = min(others, key=lambda aa: lookup(table, 0, aa))
    assert one(table, {"v": [("B", 1, "E", best)]})[1] > one(table, {"v": [("B", 1, "E", worst)]})[1]


@pytest.mark.parametrize(
    "mutation, expected",
    [
        (("B", 2, "A", "G"), "wild type A at B2, the structure has V"),
        (("B", 9, "E", "G"), "position 9 is outside chain B, which has 5 residues"),
        (("B", 0, "E", "G"), "position 0 is outside chain B"),
        (("H", 1, "E", "G"), "chain 'H' is not in the structure"),
        (("B", 1, "E", "X"), "'X' is not one of the 20 amino acids"),
    ],
)
def test_a_candidate_that_does_not_fit_gets_no_score_and_a_reason(table, mutation, expected):
    _, score, raw, note = one(table, {"v": [mutation]})
    assert score is None and raw is None
    assert expected in note


def test_several_mutations_get_no_score_on_the_table_route_because_a_lookup_would_add_singles(table):
    _, score, _, note = one(table, {"pair": [("B", 1, "E", "D"), ("A", 2, "I", "V")]})
    assert score is None and "2 mutations" in note


def test_one_bad_candidate_does_not_stop_the_others(table):
    rows = score_with_table(table, {"ok": [("B", 1, "E", "D")], "bad": [("B", 2, "A", "G")], "ok2": [("A", 2, "I", "V")]}, scorer=fake_scorer)
    assert [r[0] for r in rows] == ["ok", "bad", "ok2"]
    assert [r[1] is None for r in rows] == [False, True, False]


def test_a_table_AntiFold_cannot_score_is_refused(table):
    with pytest.raises(ScoreError, match="pdb_chain"):
        score_with_table(table.drop(columns=["pdb_chain"]), {}, scorer=fake_scorer)
    odd = table.copy()
    odd.loc[0, "pdb_res"] = "X"
    with pytest.raises(ScoreError, match="cannot take"):
        score_with_table(odd, {}, scorer=fake_scorer)


# ---- complex route -----------------------------------------------------------------------------------------------


def complex_rows(candidates):
    fake = FakeComplex()
    return fake, {row[0]: row for row in score_with_complex(fake, candidates)}


def test_a_single_mutation_is_the_change_in_total_log_likelihood_of_its_chain():
    _, rows = complex_rows({"v": [("B", 3, "Q", "E")]})
    _, score, raw, note = rows["v"]
    assert note == "" and score == pytest.approx(-(0.5 + 0.1 * 2))
    assert raw == pytest.approx(score)  # the fake's wild type has log-likelihood 0, so the mutant's total is the change


def test_several_mutations_in_one_chain_are_scored_as_a_whole_not_as_a_sum_of_singles():
    _, rows = complex_rows({"a": [("B", 1, "E", "D")], "b": [("B", 3, "Q", "E")], "both": [("B", 1, "E", "D"), ("B", 3, "Q", "E")]})
    singles = rows["a"][1] + rows["b"][1]
    assert rows["both"][1] == pytest.approx(-(0.5 + 0.7) - 0.3)
    assert rows["both"][1] < singles - 0.2


def test_mutations_in_both_chains_add_the_two_chain_scores():
    _, rows = complex_rows({"h": [("B", 2, "V", "L")], "l": [("A", 2, "I", "V")], "hl": [("B", 2, "V", "L"), ("A", 2, "I", "V")]})
    assert rows["hl"][1] == pytest.approx(rows["h"][1] + rows["l"][1])


def test_the_order_of_the_mutations_in_a_candidate_does_not_matter():
    _, rows = complex_rows({"x": [("B", 1, "E", "D"), ("B", 4, "L", "I")], "y": [("B", 4, "L", "I"), ("B", 1, "E", "D")]})
    assert rows["x"][1] == pytest.approx(rows["y"][1])


def test_one_pass_per_chain_and_the_same_mutant_chain_is_scored_once():
    fake, rows = complex_rows({
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
        ([("H", 1, "E", "G")], "chain 'H' is not in the structure"),
        ([("B", 1, "E", "X")], "'X' is not one of the 20 amino acids"),
        ([("B", 1, "E", "D"), ("B", 1, "E", "G")], "B1 is mutated twice in one candidate"),
    ],
)
def test_the_complex_route_refuses_a_candidate_that_does_not_fit_and_says_why(mutations, expected):
    _, rows = complex_rows({"v": mutations, "ok": [("B", 1, "E", "D")]})
    assert rows["v"][1] is None and rows["v"][2] is None and expected in rows["v"][3]
    assert rows["ok"][1] is not None  # the others are scored


def test_check_mutations_returns_an_empty_string_for_a_fitting_candidate():
    assert check_mutations({"B": HEAVY}, [("B", 1, "E", "D"), ("B", 2, "V", "L")]) == ""


# ---- AntiFold wrapper, with stand-in modules ---------------------------------------------------------------------


@pytest.fixture
def fake_antifold(monkeypatch):
    seen = {}

    def load_model(checkpoint_path=None, device="cuda"):
        seen["load_model"] = (checkpoint_path, device)
        return "model"

    def load_complex_coords(path, chains):
        seen["load_complex_coords"] = (path, chains)
        coords = {np.str_(c): np.zeros((3, 3, 3)) for c in chains}
        return coords, {np.str_("B"): np.str_(HEAVY), np.str_("A"): np.str_(LIGHT)}, {}, {}

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
    scorer = ComplexScorer(tmp_path / "1N8Z.pdb", ["B", "A"], device="cpu", weights=tmp_path / "model.pt")
    assert fake_antifold["load_model"] == (str(tmp_path / "model.pt"), "cpu")
    assert fake_antifold["load_complex_coords"] == (str(tmp_path / "1N8Z.pdb"), ["B", "A"])
    assert scorer.sequences == {"B": HEAVY, "A": LIGHT} and all(type(c) is str for c in scorer.sequences)
    assert scorer.log_likelihood("B", HEAVY) == -1.5  # the full-chain average, not the one without missing coordinates
    assert fake_antifold["score"] == ("model", "alphabet:invariant_gvp", ["A", "B"], "B", HEAVY)


def test_complex_scorer_uses_antifolds_default_weights_when_none_are_given(fake_antifold, tmp_path):
    ComplexScorer(tmp_path / "x.pdb", ["B", "A"])
    assert fake_antifold["load_model"] == (None, "cpu")


def test_complex_scorer_says_what_to_do_when_antifold_is_not_installed(monkeypatch, tmp_path):
    monkeypatch.setitem(sys.modules, "antifold", None)  # an import of it now raises ModuleNotFoundError
    with pytest.raises(ScoreError, match="AntiFold environment"):
        ComplexScorer(tmp_path / "x.pdb", ["B", "A"])


# ---- files and command line --------------------------------------------------------------------------------------


def test_output_is_read_back_by_design_pipelines_own_readers(table, tmp_path):
    candidates = tmp_path / "candidates.csv"
    _tabular.write_candidates(candidates, [
        {"variant_id": "a", "chain": "B", "position": 3, "wt": "Q", "mut": "E"},
        {"variant_id": "p", "chain": "B", "position": 1, "wt": "E", "mut": "D"},
        {"variant_id": "p", "chain": "A", "position": 2, "wt": "I", "mut": "V"},
    ])
    assert read_candidates(candidates) == {k: list(v) for k, v in _tabular.read_candidates(candidates).items()}
    path = write_scores(tmp_path / SCORES, score_with_table(table, read_candidates(candidates), scorer=fake_scorer))
    assert _tabular.read_tool_table([path], "variant_id", "score") == {"a": pytest.approx(lookup(table, 2, "E"))}
    assert path.read_text().splitlines()[0] == "variant_id,score,raw,note"
    both = write_scores(tmp_path / "both.csv", score_with_complex(FakeComplex(), read_candidates(candidates)))
    assert set(_tabular.read_tool_table([both], "variant_id", "score")) == {"a", "p"}  # the pair is scored on this route


def test_main_table_route_with_an_existing_table_writes_scores_csv(table, tmp_path, monkeypatch, capsys):
    monkeypatch.setattr(scorer_module, "default_scorer", lambda: fake_scorer)
    table.to_csv(tmp_path / "table.csv", index=False)
    _tabular.write_candidates(tmp_path / "c.csv", [{"variant_id": "a", "chain": "A", "position": 2, "wt": "I", "mut": "V"}])
    args = ["--pdb", "x.pdb", "--out_dir", str(tmp_path / "out"), "--mutations_csv", str(tmp_path / "c.csv"), "--route", "table", "--table", str(tmp_path / "table.csv")]
    assert main(args) == 0
    assert _tabular.read_tool_table([tmp_path / "out" / SCORES], "variant_id", "score")["a"] == pytest.approx(lookup(table, len(HEAVY) + 1, "V"))
    assert "1 of 1 candidates scored (table route)" in capsys.readouterr().out


def test_main_complex_route_is_the_default_scores_pairs_and_keeps_the_gpus_away(tmp_path, monkeypatch, capsys):
    made = {}

    class Stub(FakeComplex):
        def __init__(self, pdb, chains, device, weights):
            super().__init__()
            made.update(pdb=pdb, chains=chains, device=device, weights=weights, cuda=os.environ.get("CUDA_VISIBLE_DEVICES"))

    monkeypatch.setattr(scorer_module, "ComplexScorer", Stub)
    monkeypatch.setenv("CUDA_VISIBLE_DEVICES", "0")
    src = tmp_path / "src"
    (src / "models").mkdir(parents=True)
    (src / "models" / "model.pt").write_text("weights")
    _tabular.write_candidates(tmp_path / "c.csv", [
        {"variant_id": "pair", "chain": "B", "position": 1, "wt": "E", "mut": "D"},
        {"variant_id": "pair", "chain": "B", "position": 3, "wt": "Q", "mut": "E"},
    ])
    args = ["--pdb", "x.pdb", "--heavy", "B", "--light", "A", "--out_dir", str(tmp_path / "out"), "--mutations_csv", str(tmp_path / "c.csv"), "--antifold_src", str(src)]
    assert main(args) == 0
    assert made == {"pdb": Path("x.pdb"), "chains": ["B", "A"], "device": "cpu", "weights": src / "models" / "model.pt", "cuda": ""}
    assert _tabular.read_tool_table([tmp_path / "out" / SCORES], "variant_id", "score")["pair"] == pytest.approx(-(0.5 + 0.7) - 0.3)
    assert "1 of 1 candidates scored (complex route)" in capsys.readouterr().out


@pytest.mark.parametrize(
    "extra",
    [[], ["--heavy", "B"], ["--table", "t.csv", "--heavy", "B", "--light", "A"]],
)
def test_main_complex_route_needs_both_chains_and_refuses_a_table(tmp_path, extra):
    with pytest.raises(SystemExit):
        main(["--pdb", "x.pdb", "--out_dir", str(tmp_path), "--mutations_csv", str(tmp_path / "c.csv"), *extra])


def test_main_table_route_needs_the_antifold_folder_unless_a_table_is_given(tmp_path, monkeypatch):
    monkeypatch.delenv("ANTIFOLD_SRC", raising=False)
    (tmp_path / "c.csv").write_text("variant_id,chain,position,wt,mut\n")
    with pytest.raises(SystemExit):
        main(["--pdb", "x.pdb", "--out_dir", str(tmp_path), "--mutations_csv", str(tmp_path / "c.csv"), "--route", "table"])


def test_main_table_route_takes_both_chain_ids_or_neither(tmp_path):
    (tmp_path / "c.csv").write_text("variant_id,chain,position,wt,mut\n")
    with pytest.raises(SystemExit):
        main(["--pdb", "x.pdb", "--out_dir", str(tmp_path), "--mutations_csv", str(tmp_path / "c.csv"), "--route", "table",
              "--antifold_src", str(tmp_path), "--heavy", "B"])


def test_run_antifold_starts_it_in_its_source_folder_with_cuda_hidden(table, tmp_path, monkeypatch):
    src = tmp_path / "src"
    (src / "antifold").mkdir(parents=True)
    table.to_csv(tmp_path / "fixture.csv", index=False)
    (src / "antifold" / "main.py").write_text(textwrap.dedent(f"""
        import os, shutil, sys
        a = sys.argv
        out = a[a.index("--out_dir") + 1]
        os.makedirs(out, exist_ok=True)
        shutil.copy({str(tmp_path / "fixture.csv")!r}, os.path.join(out, "1N8Z_BA.csv"))
        open(os.path.join(out, "seen.txt"), "w").write("|".join([os.getcwd(), os.environ.get("CUDA_VISIBLE_DEVICES", "unset"), a[a.index("--pdb_file") + 1], a[a.index("--device") + 1]]))
    """))
    (tmp_path / "1N8Z.pdb").write_text("")
    monkeypatch.chdir(tmp_path)  # a relative --pdb must still reach AntiFold, which runs somewhere else
    result = run_antifold("1N8Z.pdb", "B", "A", tmp_path / "out", src, executable=sys.executable)
    assert result == (tmp_path / "out").resolve() / "1N8Z_BA.csv"
    cwd, cuda, pdb, device = (tmp_path / "out" / "seen.txt").read_text().split("|")
    assert (cwd, cuda, device) == (str(src.resolve()), "", "cpu")
    assert pdb == str((tmp_path / "1N8Z.pdb").resolve())


def test_a_failing_antifold_run_reports_its_last_output(tmp_path):
    src = tmp_path / "src"
    (src / "antifold").mkdir(parents=True)
    (src / "antifold" / "main.py").write_text("import sys; print('boom: no weights'); sys.exit(3)\n")
    with pytest.raises(ScoreError, match=r"(?s)exit 3.*boom: no weights"):
        run_antifold("x.pdb", "B", "A", tmp_path / "out", src, executable=sys.executable)


def fake_antifold_main(src: Path, tables: tuple[str, ...]) -> None:
    """An AntiFold main.py that writes the named csv files and the command line it was given."""
    (src / "antifold").mkdir(parents=True)
    (src / "antifold" / "main.py").write_text(textwrap.dedent(f"""
        import os, sys
        out = sys.argv[sys.argv.index("--out_dir") + 1]
        os.makedirs(out, exist_ok=True)
        for name in {tables!r}:
            open(os.path.join(out, name), "w").write("pdb_chain,pdb_res\\n")
        open(os.path.join(out, "argv.txt"), "w").write(" ".join(sys.argv[1:]))
    """))


def test_run_antifold_without_chain_ids_leaves_them_to_antifold_and_takes_the_one_table(tmp_path):
    src = tmp_path / "src"
    fake_antifold_main(src, ("1N8Z_BA.csv",))
    result = run_antifold("1N8Z.pdb", None, None, tmp_path / "out", src, executable=sys.executable)
    assert result == (tmp_path / "out").resolve() / "1N8Z_BA.csv"
    argv = (tmp_path / "out" / "argv.txt").read_text()
    assert "--heavy_chain" not in argv and "--light_chain" not in argv
    assert "--pdb_file" in argv and "--device cpu" in argv


def test_run_antifold_without_chain_ids_refuses_two_tables_rather_than_pick_one(tmp_path):
    src = tmp_path / "src"
    fake_antifold_main(src, ("1N8Z_BA.csv", "1N8Z_HL.csv"))
    with pytest.raises(ScoreError, match="tables found"):
        run_antifold("1N8Z.pdb", None, None, tmp_path / "out", src, executable=sys.executable)


def test_run_antifold_takes_both_chain_ids_or_neither(tmp_path):
    with pytest.raises(ScoreError, match="both chain ids or neither"):
        run_antifold("x.pdb", "B", None, tmp_path / "out", tmp_path, executable=sys.executable)
