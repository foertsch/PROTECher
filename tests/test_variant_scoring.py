"""Step 3: the fan-out, the join, the per-scorer cache and the cost classes."""

from __future__ import annotations

import math
from dataclasses import dataclass
from pathlib import Path
from typing import Any, Sequence

import pytest

from design_pipeline import (
    REGISTRY,
    Design,
    MissingDependencyError,
    RunFailedError,
    enumerate_variants,
    get_structure,
    score_variants,
    select_pairs,
)
from design_pipeline import api
from design_pipeline.adapters._scorers import DummyScorer
from design_pipeline.adapters._tabular import SCORES, read_candidates, write_scores


@pytest.fixture
def complex_(tmp_path: Path):
    return get_structure(tmp_path / "runs", model_name="dummy_complex")


@pytest.fixture
def singles(complex_):
    return enumerate_variants(complex_, chain="H", mutations="AW")


@dataclass
class Counting:
    """A scorer that records how many times it actually ran."""

    name: str = "counting"
    per_variant: bool = False
    calls: int = 0
    seen: int = 0

    def run(
        self,
        inputs: Sequence[Design],
        output_dir: Path,
        *,
        candidates: Path | str | None = None,
        extra_args: Sequence[str] | None = None,
        **kwargs: Any,
    ) -> list[Design]:
        self.calls += 1
        wanted = read_candidates(Path(str(candidates)))
        self.seen = len(wanted)
        path = write_scores(
            output_dir / SCORES, ((v, float(i), float(i)) for i, v in enumerate(wanted))
        )
        return [Design(id="scores", path=path, parent_id=None, producer=self.name)]


@pytest.fixture
def counting():
    adapter = Counting()
    REGISTRY["counting"] = adapter
    try:
        yield adapter
    finally:
        del REGISTRY["counting"]


def test_each_model_appends_one_column_and_its_raw_twin(complex_, singles) -> None:
    scores = score_variants(
        complex_, singles, models=["dummy_scorer", "dummy_per_variant"], max_variants=5
    )

    assert "score_dummy_scorer" in scores.columns
    assert "score_dummy_scorer_raw" in scores.columns
    assert scores.score_columns == ("score_dummy_scorer", "score_dummy_per_variant")
    # The join never multiplies rows.
    assert scores.rows == singles.rows


def test_nothing_is_averaged_into_a_combined_column(complex_, singles) -> None:
    """No weighted ensemble, no learned ranker, no default weights: ranking is
    the caller's business and the weights come out of the harness."""
    scores = score_variants(complex_, singles, models=["dummy_scorer"])

    assert not any(
        c in scores.columns for c in ("score", "combined", "consensus_score", "rank")
    )


def test_a_scorer_runs_once_per_candidate_set_and_is_cached(
    complex_, singles, counting
) -> None:
    first = score_variants(complex_, singles, models=["counting"])
    second = score_variants(complex_, singles, models=["counting"])

    assert counting.calls == 1
    assert first.frame["score_counting"].equals(second.frame["score_counting"])


def test_adding_a_model_does_not_invalidate_the_one_that_already_ran(
    complex_, singles, counting
) -> None:
    """The reason the key is per scorer: the harness re-scores held-out systems
    while the scorer set keeps changing, and a per-step key would throw away
    every finished column each time somebody added one."""
    score_variants(complex_, singles, models=["counting"])
    scores = score_variants(complex_, singles, models=["counting", "dummy_scorer"])

    assert counting.calls == 1
    assert "score_dummy_scorer" in scores.columns
    assert scores.frame["score_counting"].notna().all()


def test_a_matrix_scorer_sees_the_whole_table(complex_, singles, counting) -> None:
    """One forward pass answers every single mutation, so capping a cheap
    scorer would throw away columns that cost nothing."""
    score_variants(complex_, singles, models=["counting"], max_variants=3)

    assert counting.seen == singles.candidates


def test_a_per_variant_scorer_is_capped(complex_, singles) -> None:
    capped = score_variants(
        complex_, singles, models=["dummy_per_variant"], max_variants=4
    )
    answered = capped.frame["score_dummy_per_variant"].notna().sum()

    assert answered == 4
    # The rest are NaN rather than invented.
    assert answered < capped.rows


def test_an_uncapped_per_variant_scorer_refuses_a_large_table(
    monkeypatch, complex_, singles
) -> None:
    """Seconds to minutes per candidate. Refusing beats silently queueing a
    week of compute."""
    monkeypatch.setattr(api, "PER_VARIANT_LIMIT", 5)

    with pytest.raises(ValueError, match="max_variants"):
        score_variants(complex_, singles, models=["dummy_per_variant"])


def test_a_pair_is_scored_as_a_candidate_not_as_a_sum(complex_) -> None:
    """The rule the whole design hangs on: there is no code path that
    estimates a multi-mutation score from single-mutation scores."""
    singles = enumerate_variants(complex_, chain="H", mutations="A")
    scored_singles = score_variants(complex_, singles, models=["dummy_scorer"])
    pairs = select_pairs(scored_singles, max_pairs=3, min_distance=10.0)
    scored_pairs = score_variants(complex_, pairs, models=["dummy_scorer"])

    lookup = dict(
        zip(
            scored_singles.frame["variant_id"],
            scored_singles.frame["score_dummy_scorer"],
        )
    )
    frame = scored_pairs.frame
    for variant_id, group in frame.groupby("variant_id"):
        halves = [
            f"{row.chain}:{row.wt}{row.position}{row.mut}" for row in group.itertuples()
        ]
        total = sum(lookup[half] for half in halves)
        pair_score = float(group["score_dummy_scorer"].iloc[0])
        assert not math.isclose(pair_score, total, abs_tol=1e-9)
        # One score for the candidate, repeated on both of its rows.
        assert group["score_dummy_scorer"].nunique() == 1


def test_a_scorer_that_answers_an_unasked_id_is_an_error(complex_, singles) -> None:
    """A drifted join key is worth an error: silently dropping the column would
    look like a model with no opinion."""

    @dataclass
    class Wrong:
        name: str = "wrong"
        per_variant: bool = False

        def run(self, inputs, output_dir, *, candidates=None, extra_args=None, **kw):
            path = write_scores(output_dir / SCORES, [("H:Z999Z", 1.0, 1.0)])
            return [Design(id="s", path=path, parent_id=None, producer=self.name)]

    REGISTRY["wrong"] = Wrong()
    try:
        scores = score_variants(complex_, singles, models=["wrong"])
    finally:
        del REGISTRY["wrong"]
    # The row is not in the table, so the left join drops it and the column is
    # entirely NaN: the scorer answered about nothing we asked.
    assert scores.frame["score_wrong"].isna().all()


def test_a_duplicated_row_is_an_error_rather_than_a_row_explosion(
    complex_, singles
) -> None:
    @dataclass
    class Doubled:
        name: str = "doubled"
        per_variant: bool = False

        def run(self, inputs, output_dir, *, candidates=None, extra_args=None, **kw):
            first = next(iter(read_candidates(Path(str(candidates)))))
            path = write_scores(
                output_dir / SCORES, [(first, 1.0, 1.0), (first, 2.0, 2.0)]
            )
            return [Design(id="s", path=path, parent_id=None, producer=self.name)]

    REGISTRY["doubled"] = Doubled()
    try:
        with pytest.raises(RunFailedError, match="duplicate"):
            score_variants(complex_, singles, models=["doubled"])
    finally:
        del REGISTRY["doubled"]


def test_a_repeated_model_is_rejected(complex_, singles) -> None:
    with pytest.raises(ValueError, match="more than once"):
        score_variants(complex_, singles, models=["dummy_scorer", "dummy_scorer"])


def test_the_sign_convention_is_normalized_at_the_adapter_boundary(tmp_path) -> None:
    """`score` is higher-is-better whatever the tool reports; `raw` keeps the
    tool's own number so a flipped convention is recoverable without a re-run."""
    from design_pipeline.adapters._scorers import _signed

    assert _signed(2.5, -1.0) == -2.5
    assert _signed(2.5, 1.0) == 2.5
    assert _signed(None, -1.0) is None


def test_consensus_needs_a_table_and_says_where_to_get_one(complex_, singles) -> None:
    with pytest.raises(MissingDependencyError, match="OAS"):
        score_variants(complex_, singles, models=["consensus"])


def test_consensus_scores_a_single_and_declines_a_pair(tmp_path, complex_) -> None:
    """Combining two positional marginals is multiplying marginals, which is
    the additivity shortcut by another name. An empty cell is the honest
    answer."""
    table = tmp_path / "oas.csv"
    table.write_text(
        "position,aa,frequency\n1,E,0.5\n1,A,0.1\n3,Q,0.4\n3,W,0.2\n"
    )
    singles = enumerate_variants(
        complex_, chain="H", positions=[1, 3], mutations="AW"
    )
    pairs = select_pairs(singles, max_pairs=1, min_distance=1.0)

    scored = score_variants(
        complex_, singles, models=["consensus"], frequencies=str(table)
    )
    frame = scored.frame.set_index("variant_id")
    assert frame.loc["H:E1A", "score_consensus"] == pytest.approx(
        math.log((0.1 + 1e-4) / (0.5 + 1e-4))
    )
    # A position the table does not cover is no score, not a zero.
    assert math.isnan(frame.loc["H:E1W", "score_consensus"])

    scored_pairs = score_variants(
        complex_, pairs, models=["consensus"], frequencies=str(table)
    )
    assert scored_pairs.frame["score_consensus"].isna().all()


def test_a_scorer_may_be_registered_from_outside_the_package(complex_, singles):
    """The registry is a dict. Adding a model is one assignment, no entry
    points, no discovery."""
    REGISTRY["mine"] = DummyScorer(name="mine")
    try:
        scores = score_variants(complex_, singles, models=["mine"])
    finally:
        del REGISTRY["mine"]

    assert scores.frame["score_mine"].notna().all()
