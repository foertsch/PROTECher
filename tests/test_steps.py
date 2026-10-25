"""The three model steps, on the fixture adapters: offline, no GPU."""

from __future__ import annotations

from dataclasses import dataclass, replace
from pathlib import Path
from typing import Any, Sequence

import pytest

from design_pipeline import (
    REGISTRY,
    Design,
    UnknownModelError,
    design_sequences,
    make_protein_backbones,
    predict_structure,
)
from design_pipeline.adapters import dummy_structures

BB = "dummy_backbones"
SEQ = "dummy_sequences"
FOLD = "dummy_structures"
FOLD_B = "fixture_fold_b"


@pytest.fixture
def target(tmp_path: Path) -> Path:
    pdb = tmp_path / "target.pdb"
    pdb.write_text("ATOM      1  CA  GLY A   1       0.000   0.000   0.000\n")
    return pdb


def test_fresh_generation_has_no_parent(tmp_path, target) -> None:
    backbones = make_protein_backbones(
        target, tmp_path / "bb", model_name=BB, num_designs=3
    )
    assert len(backbones) == 3
    assert all(d.parent_id is None for d in backbones)
    assert all(d.path.is_file() for d in backbones)
    assert all(d.producer == BB for d in backbones)


def test_each_step_sets_parent_id_and_ancestors(tmp_path, target) -> None:
    backbones = make_protein_backbones(
        target, tmp_path / "bb", model_name=BB, num_designs=2
    )
    seqs = design_sequences(backbones, tmp_path / "seqs", model_name=SEQ, num_seqs=4)
    preds = predict_structure(seqs, tmp_path / "fold", model_name=FOLD)

    assert len(seqs) == 8  # 2 backbones x 4 sequences
    assert len(preds) == 8
    by_id = {d.id: d for d in backbones + seqs}
    for seq in seqs:
        assert by_id[seq.parent_id] in backbones
    for pred in preds:
        assert by_id[pred.parent_id] in seqs
        # parent_id alone only reaches the sequence; the backbone is one more up.
        assert pred.meta["ancestors"][0] in {d.id for d in backbones}


def test_partial_diffusion_is_the_same_function(tmp_path, target) -> None:
    backbones = make_protein_backbones(
        target, tmp_path / "bb", model_name=BB, num_designs=2
    )
    refined = make_protein_backbones(
        backbones, tmp_path / "bb2", model_name=BB, partial_T=20, num_designs=3
    )
    assert len(refined) == 6  # 2 parents x 3 designs
    assert {d.parent_id for d in refined} == {d.id for d in backbones}


def test_ids_are_unique_across_rounds(tmp_path, target) -> None:
    first = make_protein_backbones(target, tmp_path / "bb", model_name=BB, num_designs=2)
    second = make_protein_backbones(
        first, tmp_path / "bb", model_name=BB, partial_T=20, num_designs=2
    )
    assert not {d.id for d in first} & {d.id for d in second}


def test_second_run_is_a_cache_hit_with_scores_intact(tmp_path, target) -> None:
    backbones = make_protein_backbones(
        target, tmp_path / "bb", model_name=BB, num_designs=2
    )
    first = design_sequences(backbones, tmp_path / "seqs", model_name=SEQ, num_seqs=2)
    second = design_sequences(backbones, tmp_path / "seqs", model_name=SEQ, num_seqs=2)

    assert [d.meta["cached"] for d in first] == [False, False, False, False]
    assert all(d.meta["cached"] for d in second)
    assert [d.id for d in second] == [d.id for d in first]
    assert [d.path for d in second] == [d.path for d in first]
    # Scores are parsed from model output a resumed run never sees.
    assert [dict(d.scores) for d in second] == [dict(d.scores) for d in first]
    assert [d.parent_id for d in second] == [d.parent_id for d in first]


def test_identical_inputs_from_different_rounds_do_not_share_a_cache_entry(
    tmp_path, target
) -> None:
    """Two designs can hold byte-identical files and different ancestry.

    A content-only key would hand the second round the first round's records,
    so its predictions would be scored against the wrong backbone.
    """
    one = make_protein_backbones(
        target, tmp_path / "a", model_name=BB, num_designs=1, seed=0
    )
    two = make_protein_backbones(
        target, tmp_path / "b", model_name=BB, num_designs=1, seed=1
    )
    # Same file names, same bytes, different provenance.
    assert one[0].path.name == two[0].path.name
    assert one[0].id != two[0].id

    first = design_sequences(one, tmp_path / "seqs", model_name=SEQ, num_seqs=1)
    second = design_sequences(two, tmp_path / "seqs", model_name=SEQ, num_seqs=1)
    assert second[0].meta["cached"] is False
    assert first[0].parent_id == one[0].id
    assert second[0].parent_id == two[0].id


def test_two_sequence_arms_share_one_fold_run(tmp_path, target) -> None:
    """The fanout pattern: two design arms over the same backbones, folded in
    one pass. Both arms name their sequences after the same backbones, so the
    FASTA record names have to be unique or the fold outputs overwrite.
    """
    backbones = make_protein_backbones(
        target, tmp_path / "bb", model_name=BB, num_designs=2
    )
    cold = design_sequences(
        backbones, tmp_path / "seqs", model_name=SEQ, num_seqs=1, temperature=0.1
    )
    hot = design_sequences(
        backbones, tmp_path / "seqs", model_name=SEQ, num_seqs=1, temperature=0.3
    )
    assert len({d.id for d in cold + hot}) == 4

    preds = predict_structure(cold + hot, tmp_path / "fold", model_name=FOLD)
    assert len(preds) == 4
    assert len({d.id for d in preds}) == 4
    assert len({d.path for d in preds}) == 4
    # Each prediction still points back at the arm that produced it.
    assert {d.parent_id for d in preds} == {d.id for d in cold + hot}


@dataclass(frozen=True)
class RelabelledFold:
    """A second structure predictor, so the list form has two models to run.

    The fixture adapter under another producer name, registered the way a user
    registers their own model rather than shipped in the package.
    """

    producer: str

    def run(
        self, inputs: Sequence[Design], output_dir: Path, **kwargs: Any
    ) -> list[Design]:
        produced = dummy_structures.run(inputs, output_dir, **kwargs)
        return [replace(d, producer=self.producer) for d in produced]


@pytest.fixture
def second_fold_model():
    REGISTRY[FOLD_B] = RelabelledFold(FOLD_B)
    yield FOLD_B
    del REGISTRY[FOLD_B]


@pytest.fixture
def sequences(tmp_path, target) -> list[Design]:
    backbones = make_protein_backbones(
        target, tmp_path / "bb", model_name=BB, num_designs=2
    )
    return design_sequences(backbones, tmp_path / "seqs", model_name=SEQ, num_seqs=1)


def test_one_model_as_a_string_or_a_single_element_list_match(
    tmp_path, sequences
) -> None:
    plain = predict_structure(sequences, tmp_path / "a", model_name=FOLD)
    listed = predict_structure(sequences, tmp_path / "b", model_name=[FOLD])

    assert [d.id for d in listed] == [d.id for d in plain]
    assert [d.parent_id for d in listed] == [d.parent_id for d in plain]
    assert [d.producer for d in listed] == [d.producer for d in plain]
    assert [dict(d.scores) for d in listed] == [dict(d.scores) for d in plain]
    assert [d.path.name for d in listed] == [d.path.name for d in plain]
    # The run directory is the cache key, so equal keys mean the list form
    # resolved to the same config and resumes the same run.
    assert [d.meta["cache_key"] for d in listed] == [d.meta["cache_key"] for d in plain]


def test_a_list_folds_every_sequence_with_every_model(
    tmp_path, sequences, second_fold_model
) -> None:
    preds = predict_structure(
        sequences, tmp_path / "fold", model_name=[FOLD, second_fold_model]
    )

    assert len(preds) == 2 * len(sequences)  # N sequences x M models
    # Order follows the model list: the first model's predictions, then the second's.
    assert [d.producer for d in preds] == (
        [FOLD] * len(sequences) + [second_fold_model] * len(sequences)
    )
    # Each prediction keeps the sequence it was folded from, and one sequence is
    # legitimately the parent of several predictions.
    assert [d.parent_id for d in preds] == [d.id for d in sequences] * 2
    assert len({d.id for d in preds}) == len(preds)
    assert len({d.path for d in preds}) == len(preds)
    assert all(d.path.is_file() for d in preds)


def test_the_cache_is_per_model(tmp_path, sequences, second_fold_model) -> None:
    first = predict_structure(sequences, tmp_path / "fold", model_name=FOLD)
    assert all(d.meta["cached"] is False for d in first)

    mixed = predict_structure(
        sequences, tmp_path / "fold", model_name=[FOLD, second_fold_model]
    )
    # The model that already ran resumes; the one that has not still runs.
    assert all(d.meta["cached"] for d in mixed if d.producer == FOLD)
    assert not any(d.meta["cached"] for d in mixed if d.producer == second_fold_model)
    assert [d.id for d in mixed if d.producer == FOLD] == [d.id for d in first]

    again = predict_structure(
        sequences, tmp_path / "fold", model_name=[FOLD, second_fold_model]
    )
    assert all(d.meta["cached"] for d in again)
    assert [d.id for d in again] == [d.id for d in mixed]


def test_an_empty_model_list_is_rejected(tmp_path, sequences) -> None:
    with pytest.raises(ValueError, match="empty list"):
        predict_structure(sequences, tmp_path / "fold", model_name=[])


def test_a_repeated_model_is_rejected(tmp_path, sequences, second_fold_model) -> None:
    with pytest.raises(ValueError) as exc:
        predict_structure(
            sequences, tmp_path / "fold", model_name=[FOLD, second_fold_model, FOLD]
        )
    assert FOLD in str(exc.value)


def test_a_failing_model_in_a_list_raises_rather_than_returning_a_partial_result(
    tmp_path, sequences
) -> None:
    with pytest.raises(UnknownModelError):
        predict_structure(
            sequences, tmp_path / "fold", model_name=[FOLD, "alphafold9000"]
        )


def test_a_different_seed_is_a_different_run(tmp_path, target) -> None:
    """Without the seed in the key, "another 100 designs" returns the cached 100."""
    a = make_protein_backbones(target, tmp_path / "bb", model_name=BB, num_designs=2, seed=0)
    b = make_protein_backbones(target, tmp_path / "bb", model_name=BB, num_designs=2, seed=1)
    assert a[0].meta["cache_key"] != b[0].meta["cache_key"]
    assert all(d.meta["cached"] is False for d in b)
    assert {d.path for d in a} != {d.path for d in b}


def test_changed_input_content_changes_the_key(tmp_path, target) -> None:
    first = make_protein_backbones(target, tmp_path / "bb", model_name=BB, num_designs=1)
    target.write_text("ATOM      1  CA  ALA A   1       1.000   0.000   0.000\n")
    second = make_protein_backbones(target, tmp_path / "bb", model_name=BB, num_designs=1)
    assert first[0].meta["cache_key"] != second[0].meta["cache_key"]


def test_cache_false_reruns(tmp_path, target) -> None:
    make_protein_backbones(target, tmp_path / "bb", model_name=BB, num_designs=1)
    again = make_protein_backbones(
        target, tmp_path / "bb", model_name=BB, num_designs=1, cache=False
    )
    assert again[0].meta["cached"] is False


def test_extra_args_reach_the_key(tmp_path, target) -> None:
    plain = make_protein_backbones(target, tmp_path / "bb", model_name=BB, num_designs=1)
    tweaked = make_protein_backbones(
        target,
        tmp_path / "bb",
        model_name=BB,
        num_designs=1,
        extra_args=["denoiser.noise_scale_ca=0"],
    )
    assert plain[0].meta["cache_key"] != tweaked[0].meta["cache_key"]


def test_unanticipated_kwarg_changes_the_key(tmp_path, target) -> None:
    plain = make_protein_backbones(target, tmp_path / "bb", model_name=BB, num_designs=1)
    exotic = make_protein_backbones(
        target, tmp_path / "bb", model_name=BB, num_designs=1, some_new_flag=3
    )
    assert plain[0].meta["cache_key"] != exotic[0].meta["cache_key"]


def test_deleted_outputs_are_a_miss(tmp_path, target) -> None:
    first = make_protein_backbones(target, tmp_path / "bb", model_name=BB, num_designs=1)
    first[0].path.unlink()
    again = make_protein_backbones(target, tmp_path / "bb", model_name=BB, num_designs=1)
    assert again[0].meta["cached"] is False


def test_missing_input_file_is_reported(tmp_path) -> None:
    with pytest.raises(FileNotFoundError):
        make_protein_backbones(tmp_path / "nope.pdb", tmp_path / "bb", model_name=BB)


def test_unknown_model_lists_the_known_ones(tmp_path, target) -> None:
    with pytest.raises(UnknownModelError) as exc:
        make_protein_backbones(target, tmp_path / "bb", model_name="alphafold9000")
    message = str(exc.value)
    assert "rfdiffusion" in message and "dummy_backbones" in message


def test_importing_design_pipeline_pulls_in_no_model() -> None:
    import sys

    assert not {m for m in sys.modules if m.split(".")[0] in {"torch", "jax"}}
