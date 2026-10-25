"""The round loop from examples/refine.py, three rounds, offline.

The loop is a `for` statement, so what is tested here is that the five
functions compose: ancestry survives two steps, scoring finds the right
reference each round, selection is deterministic, and a rerun resumes.
"""

from __future__ import annotations

from pathlib import Path

import pytest

from design_pipeline import (
    design_sequences,
    make_protein_backbones,
    predict_structure,
    score_designs,
    select,
)

MODELS = {
    "bb": "dummy_backbones",
    "seq": "dummy_sequences",
    "fold": "dummy_structures",
}


def run_loop(root: Path, target: Path, n_rounds: int = 3) -> list[list]:
    backbones = make_protein_backbones(
        target, root / "bb0", model_name=MODELS["bb"], num_designs=4, seed=0
    )
    rounds = []
    for rnd in range(n_rounds):
        seqs = design_sequences(
            backbones, root / f"seqs{rnd}", model_name=MODELS["seq"], num_seqs=2
        )
        preds = predict_structure(
            seqs, root / f"fold{rnd}", model_name=MODELS["fold"]
        )
        scored = score_designs(preds, reference=backbones)
        survivors = select(scored, top_k=2, by="pae_interaction")
        rounds.append(survivors)
        backbones = make_protein_backbones(
            survivors,
            root / f"bb{rnd + 1}",
            model_name=MODELS["bb"],
            partial_T=20,
            num_designs=3,
            seed=rnd,
        )
    return rounds


@pytest.fixture
def target(tmp_path: Path) -> Path:
    pdb = tmp_path / "target.pdb"
    pdb.write_text("ATOM      1  CA  GLY A   1       0.000   0.000   0.000\n")
    return pdb


def test_three_rounds_run_and_narrow(tmp_path, target) -> None:
    rounds = run_loop(tmp_path / "run", target)

    assert len(rounds) == 3
    for survivors in rounds:
        assert len(survivors) == 2
        for design in survivors:
            assert design.path.is_file()
            assert design.scores["pae_interaction"] > 0
            assert design.scores["sc_rmsd"] >= 0
            assert design.meta["reference_id"] is not None


def test_each_round_scores_against_that_round_backbones(tmp_path, target) -> None:
    rounds = run_loop(tmp_path / "run", target)
    references = [{d.meta["reference_id"] for d in r} for r in rounds]
    # A later round must not be scoring against round 0's backbones.
    assert references[0] & references[1] == set()
    assert references[1] & references[2] == set()


def test_survivors_are_the_lowest_pae(tmp_path, target) -> None:
    root = tmp_path / "run"
    backbones = make_protein_backbones(
        target, root / "bb0", model_name=MODELS["bb"], num_designs=4
    )
    seqs = design_sequences(backbones, root / "s", model_name=MODELS["seq"], num_seqs=2)
    scored = score_designs(
        predict_structure(seqs, root / "f", model_name=MODELS["fold"]),
        reference=backbones,
    )
    survivors = select(scored, top_k=3, by="pae_interaction")
    best = sorted(d.scores["pae_interaction"] for d in scored)[:3]
    assert [d.scores["pae_interaction"] for d in survivors] == best


def test_rerunning_the_loop_resumes_from_cache(tmp_path, target) -> None:
    root = tmp_path / "run"
    first = run_loop(root, target)
    second = run_loop(root, target)

    assert [[d.id for d in r] for r in second] == [[d.id for d in r] for r in first]
    assert [[d.scores["pae_interaction"] for d in r] for r in second] == [
        [d.scores["pae_interaction"] for d in r] for r in first
    ]
    # Every step of the second pass was a hit, including the later rounds,
    # which only line up if select is deterministic.
    assert all(d.meta["cached"] for r in second for d in r)


def test_the_loop_is_deterministic_in_a_fresh_tree(tmp_path, target) -> None:
    first = run_loop(tmp_path / "a", target)
    second = run_loop(tmp_path / "b", target)
    assert [[d.id for d in r] for r in second] == [[d.id for d in r] for r in first]
