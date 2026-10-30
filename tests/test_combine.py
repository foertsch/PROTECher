"""combine: sets of 5 to 10 substitutions from the gated singles, and the designs chosen from the scored sets.

The rules of `pick_designs` are tested on tables written by hand, so each test names the one situation it is about and
nothing else is in play. The end-to-end test uses the offline fixtures: the dummy complex and the dummy scorer.
"""

from __future__ import annotations

import sys
from pathlib import Path
from typing import Mapping

import numpy as np
import pandas as pd
import pytest

from design_pipeline import enumerate_variants, get_structure, score_variants, select_pairs
from design_pipeline.combine import _biopython_pi, _hamming, frontier, pick_designs, propose_sets
from design_pipeline.errors import MissingDependencyError
from design_pipeline.structures import ca_coordinates
from design_pipeline.variants import SCHEMA, VariantTable, make_id, parse_id, write_table

# A warning this module raises that a test did not ask for is a failure; the missing-ANARCI one is the suite's own.
pytestmark = pytest.mark.filterwarnings("error::RuntimeWarning", "ignore:anarci is not installed:RuntimeWarning")

H_SEQ = "ACDEFGHIKLMNPQRSTVWY"  # 20 positions, one of each letter, so a position's residue is easy to read


@pytest.fixture
def complex_(tmp_path: Path):
    return get_structure(tmp_path / "runs", model_name="dummy_complex")


def mutation(position: int, mut: str | None = None, chain: str = "H") -> tuple[str, int, str, str]:
    """A substitution at a position of H_SEQ: to `mut`, or to A (to W where the residue is already A)."""
    wt = H_SEQ[position - 1]
    return (chain, position, wt, mut or ("A" if wt != "A" else "W"))


def candidate(*positions: int) -> str:
    return make_id([mutation(p) for p in positions])


def table(tmp_path: Path, sets: Mapping[str, float | None], column: str = "score_x", structure=None) -> VariantTable:
    """A variant table written by hand: each id with its score, one row per mutation."""
    rows = []
    for variant_id, score in sets.items():
        for chain, position, wt, mut in parse_id(variant_id):
            rows.append({"variant_id": variant_id, "chain": chain, "position": position, "imgt": "-", "wt": wt, "mut": mut, column: score})
    return write_table(pd.DataFrame(rows), tmp_path / "t.parquet", output_dir=tmp_path, structure=structure, config={})


def ids_of(designs: VariantTable) -> list[str]:
    frame = designs.frame.drop_duplicates("variant_id").sort_values("design")
    return list(frame["variant_id"])


# ---- propose_sets ------------------------------------------------------------------------------------------------


@pytest.fixture
def singles(complex_):
    base = enumerate_variants(complex_, chain="H", mutations="AW")
    return score_variants(complex_, base, models=["dummy_scorer"])


def test_every_set_has_the_size_it_was_drawn_for_and_ids_are_unique(singles) -> None:
    sets = propose_sets(singles, sizes=(5, 6, 7), per_size=15, by="score_dummy_scorer", pool=30)
    frame = sets.frame

    per_id = frame.groupby("variant_id").agg(n=("position", "size"), declared=("n_mutations", "first"))
    assert (per_id["n"] == per_id["declared"]).all()
    assert sorted(per_id["n"].unique()) == [5, 6, 7]
    assert sets.candidates == len(per_id) == sum(sets.config["found"].values())
    assert all(count <= 15 for count in sets.config["found"].values())


def test_a_set_never_changes_a_position_twice_and_the_ids_are_canonical(singles) -> None:
    frame = propose_sets(singles, sizes=(5, 8), per_size=25, by="score_dummy_scorer", pool=30).frame

    for variant_id, group in frame.groupby("variant_id"):
        assert group["position"].is_unique
        assert variant_id == make_id(list(zip(group["chain"], group["position"], group["wt"], group["mut"])))


def test_the_sets_are_drawn_from_the_best_singles_only(singles) -> None:
    ranked = singles.frame.sort_values("score_dummy_scorer", ascending=False, kind="stable").head(12)
    shortlist = set(ranked["variant_id"])
    frame = propose_sets(singles, sizes=(5,), per_size=20, by="score_dummy_scorer", pool=12).frame

    used = {make_id([(r.chain, r.position, r.wt, r.mut)]) for r in frame.itertuples()}
    assert used <= shortlist


def test_the_same_seed_gives_the_same_sets_and_another_seed_gives_others(singles) -> None:
    def drawn(seed):
        return set(propose_sets(singles, sizes=(6,), per_size=20, by="score_dummy_scorer", pool=30, seed=seed).frame["variant_id"])

    assert drawn(1) == drawn(1)
    assert drawn(1) != drawn(2)


def test_no_two_substitutions_in_a_set_are_closer_than_min_distance(singles) -> None:
    sets = propose_sets(singles, sizes=(5,), per_size=20, by="score_dummy_scorer", pool=30, min_distance=7.5)
    coords = ca_coordinates(sets.structure, "H")

    assert sets.config["distance_units"] == "angstrom"
    for _id, group in sets.frame.groupby("variant_id"):
        points = coords[group["position"].to_numpy() - 1]
        for i in range(len(points)):
            for j in range(i):
                assert np.linalg.norm(points[i] - points[j]) >= 7.5


def test_a_distance_nothing_satisfies_is_an_error_that_says_what_to_relax(singles) -> None:
    with pytest.raises(ValueError, match="relax a limit, lower min_distance|Relax a limit, lower min_distance"):
        propose_sets(singles, sizes=(5,), per_size=5, by="score_dummy_scorer", pool=30, min_distance=1000.0)


def test_a_limit_caps_the_substitutions_a_set_takes_from_one_region(singles) -> None:
    regions = {("H", position): "H3" for position in range(1, 31)}
    frame = propose_sets(
        singles, sizes=(5,), per_size=20, by="score_dummy_scorer", pool=60, regions=regions, limits={"H3": 2}
    ).frame

    in_region = frame[frame["position"] <= 30].groupby("variant_id").size()
    assert (in_region <= 2).all() and len(in_region) > 0


def test_limits_without_regions_is_refused(singles) -> None:
    with pytest.raises(ValueError, match="limits needs regions"):
        propose_sets(singles, sizes=(5,), by="score_dummy_scorer", limits={"H3": 1})


def test_a_table_that_already_holds_pairs_is_refused(complex_, singles) -> None:
    pairs = select_pairs(singles, max_pairs=5, min_distance=0.0)

    with pytest.raises(ValueError, match="already have more than one"):
        propose_sets(pairs, by="score_dummy_scorer")


def test_scores_are_not_carried_over_and_the_columns_are_the_schema_plus_the_size(singles) -> None:
    sets = propose_sets(singles, sizes=(5,), per_size=5, by="score_dummy_scorer", pool=30)

    assert list(sets.frame.columns) == [*SCHEMA, "n_mutations"]
    assert not sets.score_columns


def test_singles_without_a_score_are_left_out_and_too_few_for_the_size_is_an_error(tmp_path) -> None:
    scored = table(tmp_path, {candidate(1): 1.0, candidate(3): float("nan"), candidate(5): 2.0, candidate(7): 3.0})

    with pytest.raises(ValueError, match="only 3 singles have a 'score_x' score"):
        propose_sets(scored, sizes=(4,), by="score_x")


def test_an_unknown_ranking_column_names_the_ones_there_are(singles) -> None:
    with pytest.raises(KeyError, match="no column 'nope'"):
        propose_sets(singles, by="nope")


def test_sizes_below_two_are_refused(singles) -> None:
    with pytest.raises(ValueError, match="at least 2"):
        propose_sets(singles, sizes=(1,), by="score_dummy_scorer")


# ---- pick_designs ------------------------------------------------------------------------------------------------


def test_designs_come_out_best_first_and_are_numbered(tmp_path) -> None:
    scored = table(tmp_path, {candidate(1, 2, 3, 4, 5): 1.0, candidate(6, 7, 8, 9, 10): 3.0, candidate(11, 12, 13, 14, 15): 2.0})

    designs = pick_designs(scored, by="score_x", n=3, min_uses=1, max_pi_shift=None)

    assert ids_of(designs) == [candidate(6, 7, 8, 9, 10), candidate(11, 12, 13, 14, 15), candidate(1, 2, 3, 4, 5)]
    assert list(designs.frame["design"].drop_duplicates()) == [1, 2, 3]
    assert set(designs.frame["n_mutations"]) == {5}
    assert "score_x" in designs.frame.columns  # each set's own score, kept as it is


def test_designs_must_differ_from_each_other_at_min_differences_positions(tmp_path) -> None:
    a, near, far = candidate(1, 2, 3, 4, 5), candidate(1, 2, 3, 4, 6), candidate(1, 2, 3, 6, 7)
    scored = table(tmp_path, {a: 10.0, near: 9.0, far: 8.0})

    designs = pick_designs(scored, by="score_x", n=2, min_uses=1, max_pi_shift=None)

    assert ids_of(designs) == [a, far]  # `near` differs from `a` at only 2 positions
    for left, right in [(a, far)]:
        assert _hamming(parse_id(left), parse_id(right)) >= 3


def test_hamming_counts_a_position_once_and_treats_no_change_as_the_parent() -> None:
    assert _hamming([("H", 1, "A", "G")], [("H", 1, "A", "W")]) == 1  # same position, different residues
    assert _hamming([("H", 1, "A", "G")], [("H", 2, "C", "G")]) == 2  # each leaves the other at the parent's residue
    assert _hamming([("H", 1, "A", "G"), ("H", 2, "C", "G")], [("H", 1, "A", "G"), ("H", 2, "C", "G")]) == 0


def test_no_substitution_is_used_in_more_than_max_uses_designs(tmp_path) -> None:
    shared = [candidate(1, *rest) for rest in [(2, 3, 4, 5), (6, 7, 8, 9), (10, 11, 12, 13)]]
    other = candidate(14, 15, 16, 17, 18)
    scored = table(tmp_path, {shared[0]: 10.0, shared[1]: 9.0, shared[2]: 8.0, other: 7.0})

    with pytest.warns(RuntimeWarning, match="chose 3 of 15"):
        designs = pick_designs(scored, by="score_x", n=15, max_uses=2, min_uses=1, max_pi_shift=None)

    assert ids_of(designs) == [shared[0], shared[1], other]  # the third design with H:A1G would be its third use
    assert designs.config["uses"][make_id([mutation(1)])] == 2


def test_a_use_limit_for_a_region_caps_its_substitutions_below_max_uses(tmp_path) -> None:
    regions = {("H", 1): "H3"}
    shared = [candidate(1, *rest) for rest in [(2, 3, 4, 5), (6, 7, 8, 9), (10, 11, 12, 13)]]
    scored = table(tmp_path, {shared[0]: 10.0, shared[1]: 9.0, shared[2]: 8.0})

    with pytest.warns(RuntimeWarning, match="chose 2 of 15"):
        designs = pick_designs(scored, by="score_x", n=15, regions=regions, use_limits={"H3": 2}, min_uses=1, max_pi_shift=None)

    assert ids_of(designs) == shared[:2]


def test_use_limits_without_regions_is_refused(tmp_path) -> None:
    with pytest.raises(ValueError, match="use_limits needs regions"):
        pick_designs(table(tmp_path, {candidate(1, 2, 3, 4, 5): 1.0}), by="score_x", use_limits={"H3": 1})


def test_equal_scores_go_to_the_candidate_with_fewer_cdr_changes(tmp_path) -> None:
    in_loop, outside = candidate(1, 2, 3, 4, 20), candidate(5, 6, 7, 8, 9)
    scored = table(tmp_path, {in_loop: 5.0, outside: 5.0})

    designs = pick_designs(scored, by="score_x", n=2, regions={("H", 20): "H3"}, cdr_labels={"H3"}, min_uses=1, max_pi_shift=None)

    assert ids_of(designs) == [outside, in_loop]


def test_without_a_quota_the_smallest_sets_win_and_with_one_each_size_gets_its_share(tmp_path) -> None:
    five_a, five_b = candidate(1, 2, 3, 4, 5), candidate(6, 7, 8, 9, 10)
    six_a, six_b = candidate(11, 12, 13, 14, 15, 16), candidate(1, 7, 13, 17, 18, 19)
    scored = table(tmp_path, {five_a: -10.0, five_b: -11.0, six_a: -20.0, six_b: -21.0})

    # Read each result before the next call: the table is the file on disk, and every call writes the same default path.
    plain = ids_of(pick_designs(scored, by="score_x", n=2, min_uses=1, max_pi_shift=None))
    shared = ids_of(pick_designs(scored, by="score_x", n=2, quota={5: 1, 6: 1}, min_uses=1, max_pi_shift=None))
    sixes_only = ids_of(pick_designs(scored, by="score_x", n=2, quota={6: 2}, min_uses=1, max_pi_shift=None))

    assert plain == [five_a, five_b]
    assert shared == [five_a, six_a]
    assert sixes_only == [six_a, six_b]  # a size the quota does not name gets nothing


def test_candidates_the_scorer_could_not_score_are_skipped_and_counted(tmp_path) -> None:
    scored = table(tmp_path, {candidate(1, 2, 3, 4, 5): 1.0, candidate(6, 7, 8, 9, 10): float("nan")})

    with pytest.warns(RuntimeWarning, match="chose 1 of 15"):
        designs = pick_designs(scored, by="score_x", min_uses=1, max_pi_shift=None)

    assert designs.config["unscored"] == 1 and ids_of(designs) == [candidate(1, 2, 3, 4, 5)]


def test_substitutions_in_fewer_than_min_uses_designs_are_reported_not_enforced(tmp_path) -> None:
    a, b = candidate(1, 2, 3, 4, 5), candidate(1, 2, 3, 6, 7)
    scored = table(tmp_path, {a: 2.0, b: 1.0})

    with pytest.warns(RuntimeWarning, match="fewer than 2 designs"):
        designs = pick_designs(scored, by="score_x", n=2, max_pi_shift=None)

    assert ids_of(designs) == [a, b]
    assert make_id([mutation(4)]) in designs.config["underused"]  # only in one design
    assert make_id([mutation(1)]) not in designs.config["underused"]  # in both


def test_nothing_to_choose_from_is_an_error(tmp_path) -> None:
    scored = table(tmp_path, {candidate(1, 2, 3, 4, 5): float("nan")})

    with pytest.warns(RuntimeWarning), pytest.raises(ValueError, match="no design could be chosen"):
        pick_designs(scored, by="score_x", max_pi_shift=None)


def test_an_unknown_ranking_column_is_an_error_here_too(tmp_path) -> None:
    with pytest.raises(KeyError, match="no column 'nope'"):
        pick_designs(table(tmp_path, {candidate(1, 2, 3, 4, 5): 1.0}), by="nope")


# ---- the isoelectric point rule ----------------------------------------------------------------------------------


def charge(sequence: str) -> float:
    """A stand-in for the isoelectric point: more K and R, higher; more D and E, lower."""
    return sum(c in "KR" for c in sequence) - sum(c in "DE" for c in sequence)


def test_a_design_that_moves_the_isoelectric_point_too_far_is_not_chosen(tmp_path) -> None:
    # Position 9 is K in H_SEQ: K to A loses one positive charge, K to E loses it and adds a negative one. The other
    # positions (A, C, F, G, H, I, L, M, N, P) carry no charge, so the shift is the one at position 9 alone.
    mild = make_id([mutation(p) for p in (1, 2, 5, 6, 7)] + [mutation(9, "A")])
    harsh = make_id([mutation(p) for p in (8, 10, 11, 12, 13)] + [mutation(9, "E")])
    scored = table(tmp_path, {harsh: 10.0, mild: 9.0})

    with pytest.warns(RuntimeWarning, match="chose 1 of 2"):
        designs = pick_designs(
            scored, by="score_x", n=2, min_uses=1, max_pi_shift=1.5, sequences={"H": H_SEQ}, pi_function=charge
        )

    assert ids_of(designs) == [mild]  # a shift of 1 passes and a shift of 2 does not, although it scored better


def test_the_rule_can_be_skipped(tmp_path) -> None:
    scored = table(tmp_path, {candidate(1, 2, 3, 4, 5): 1.0})

    designs = pick_designs(scored, by="score_x", n=1, min_uses=1, max_pi_shift=None, pi_function=lambda s: 1 / 0)

    assert len(ids_of(designs)) == 1  # the function was never called


def test_the_rule_needs_sequences_from_somewhere(tmp_path) -> None:
    with pytest.raises(ValueError, match="needs the chain sequences"):
        pick_designs(table(tmp_path, {candidate(1, 2, 3, 4, 5): 1.0}), by="score_x", pi_function=charge)


def test_a_substitution_that_does_not_match_the_sequence_is_an_error(tmp_path) -> None:
    scored = table(tmp_path, {candidate(1, 2, 3, 4, 5): 1.0})

    with pytest.raises(ValueError, match="does not match the sequence of chain H"):
        pick_designs(scored, by="score_x", sequences={"H": "G" * 20}, pi_function=charge)


def test_without_biopython_the_error_says_how_to_get_it_or_skip_the_rule(tmp_path, monkeypatch) -> None:
    monkeypatch.setitem(sys.modules, "Bio", None)  # an import of it now fails
    scored = table(tmp_path, {candidate(1, 2, 3, 4, 5): 1.0})

    with pytest.raises(MissingDependencyError, match="pip install biopython.*max_pi_shift=None"):
        pick_designs(scored, by="score_x", sequences={"H": H_SEQ})


def test_biopythons_isoelectric_point_moves_the_way_the_residues_say() -> None:
    pytest.importorskip("Bio")

    assert _biopython_pi("AAAKKKAAA") > _biopython_pi("AAAAAAAAA") > _biopython_pi("AAAEEEAAA")


# ---- frontier ----------------------------------------------------------------------------------------------------


def test_the_frontier_is_the_best_candidate_at_each_size_and_skips_unscored_ones(tmp_path) -> None:
    sets = {
        candidate(1, 2, 3, 4, 5): 3.0,
        candidate(6, 7, 8, 9, 10): 5.0,
        candidate(1, 2, 3, 4, 5, 6): 4.0,
        candidate(11, 12, 13, 14, 15, 16, 17): float("nan"),
    }
    scored = table(tmp_path, sets)

    best = frontier(scored, by="score_x")
    two = frontier(scored, by="score_x", keep=2)

    assert list(best["n_mutations"]) == [5, 6]
    assert list(best["variant_id"]) == [candidate(6, 7, 8, 9, 10), candidate(1, 2, 3, 4, 5, 6)]
    assert list(two["n_mutations"]) == [5, 5, 6]
    assert list(two["score_x"]) == [5.0, 3.0, 4.0]


# ---- the whole path, offline -------------------------------------------------------------------------------------


def test_singles_to_sets_to_scores_to_designs_with_the_offline_fixtures(complex_) -> None:
    base = enumerate_variants(complex_, chain="H", mutations="AW")
    scored_singles = score_variants(complex_, base, models=["dummy_scorer"])
    sets = propose_sets(scored_singles, sizes=(5, 6), per_size=20, by="score_dummy_scorer", pool=40, seed=3)
    scored = score_variants(complex_, sets, models=["dummy_scorer"])

    designs = pick_designs(
        scored, by="score_dummy_scorer", n=4, min_uses=1, max_uses=4, max_pi_shift=None, quota={5: 2, 6: 2}
    )
    curve = frontier(scored, by="score_dummy_scorer")

    chosen = ids_of(designs)
    assert len(chosen) == 4 and len(set(chosen)) == 4
    for i, left in enumerate(chosen):
        assert all(_hamming(parse_id(left), parse_id(right)) >= 3 for right in chosen[:i])
    assert set(designs.frame["n_mutations"]) == {5, 6}
    assert list(curve["n_mutations"]) == [5, 6]
    assert scored.candidates == sets.candidates  # the sets were scored whole, one row of score per candidate
