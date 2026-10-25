"""The record, and the three things you do to a list of them."""

from __future__ import annotations

from pathlib import Path

import pytest

from design_pipeline import Design, lineage, select, to_dataframe


def make(
    id: str,
    parent: str | None = None,
    producer: str = "test",
    **scores: float,
) -> Design:
    return Design(
        id=id,
        path=Path(f"{id}.pdb"),
        parent_id=parent,
        producer=producer,
        scores=scores,
    )


def test_score_miss_lists_what_is_present() -> None:
    design = make("a", plddt=80.0, ptm=0.7)
    with pytest.raises(KeyError) as exc:
        design.scores["pae_interaction"]
    message = str(exc.value)
    assert "pae_interaction" in message
    assert "plddt" in message and "ptm" in message


def test_score_miss_on_an_unscored_design_says_none() -> None:
    with pytest.raises(KeyError, match="none"):
        make("a").scores["plddt"]


def test_design_is_frozen() -> None:
    design = make("a")
    with pytest.raises(Exception):
        design.id = "b"  # type: ignore[misc]


def test_lineage_walks_across_steps() -> None:
    pool = [
        make("bb"),
        make("seq", parent="bb"),
        make("pred", parent="seq"),
        make("unrelated"),
    ]
    chain = lineage(pool, "pred")
    assert [d.id for d in chain] == ["bb", "seq", "pred"]


def test_lineage_stops_where_the_chain_leaves_the_pool() -> None:
    pool = [make("pred", parent="seq_not_here")]
    assert [d.id for d in lineage(pool, "pred")] == ["pred"]


def test_lineage_survives_a_cycle() -> None:
    # Not expected, but an infinite loop here would hang a campaign.
    pool = [make("a", parent="b"), make("b", parent="a")]
    assert len(lineage(pool, "a")) == 2


def test_lineage_unknown_id() -> None:
    with pytest.raises(KeyError):
        lineage([make("a")], "nope")


def test_select_sorts_truncates_and_filters() -> None:
    pool = [make("a", pae=9.0), make("b", pae=3.0), make("c", pae=5.0)]
    assert [d.id for d in select(pool, 2, by="pae")] == ["b", "c"]
    assert [d.id for d in select(pool, 2, by="pae", ascending=False)] == ["a", "c"]
    assert [d.id for d in select(pool, 3, by="pae", where=lambda d: d.id != "b")] == [
        "c",
        "a",
    ]


def test_select_breaks_ties_on_id_in_both_directions() -> None:
    # A nondeterministic selection would shift every downstream cache key.
    pool = [make(i, pae=1.0) for i in ("d", "b", "c", "a")]
    assert [d.id for d in select(pool, 4, by="pae")] == ["a", "b", "c", "d"]
    assert [d.id for d in select(pool, 4, by="pae", ascending=False)] == [
        "a",
        "b",
        "c",
        "d",
    ]
    shuffled = [pool[2], pool[0], pool[3], pool[1]]
    assert [d.id for d in select(shuffled, 2, by="pae")] == ["a", "b"]


def test_select_missing_score_names_the_design_and_the_keys() -> None:
    pool = [make("a", pae=1.0), make("b", plddt=80.0)]
    with pytest.raises(KeyError) as exc:
        select(pool, 1, by="pae")
    message = str(exc.value)
    assert "'b'" in message
    assert "plddt" in message


def test_select_on_a_mixed_pool_points_at_filtering_by_producer() -> None:
    """A `predict_structure` call with a list of models returns a mixed pool,
    and the models do not report the same metrics."""
    pool = [
        make("a", producer="colabfold", pae_interaction=7.0),
        make("b", producer="boltz", affinity=0.4),
    ]
    with pytest.raises(KeyError) as exc:
        select(pool, 1, by="pae_interaction")
    message = str(exc.value)
    assert "producer" in message
    assert "boltz" in message and "colabfold" in message


def test_select_on_a_single_producer_pool_does_not_mention_filtering() -> None:
    pool = [make("a", pae=1.0), make("b", plddt=80.0)]
    with pytest.raises(KeyError) as exc:
        select(pool, 1, by="pae")
    assert "producer" not in str(exc.value)


def test_to_dataframe_flattens_scores() -> None:
    pytest.importorskip("pandas")
    pool = [make("a", plddt=80.0), make("b", plddt=90.0, ptm=0.8)]
    frame = to_dataframe(pool)
    assert list(frame["id"]) == ["a", "b"]
    assert list(frame["plddt"]) == [80.0, 90.0]
    assert frame["ptm"].isna()[0]
