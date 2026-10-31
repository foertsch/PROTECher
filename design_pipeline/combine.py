"""Combine single substitutions into variants of 5 to 10, and choose designs among the scored ones.

`select_pairs` stops at two. Track A needs at least five substitutions per design, so this module takes the gated
singles to sets of five to ten, and then picks the designs that go forward. It is three plain functions and a
loop, like the rest of `variants`: nothing here runs a model, so nothing here is an adapter.

    sets = propose_sets(survivors, by="score_antifold", regions=regions)     # singles -> candidate sets
    scored = score_variants(complex_, sets, models=["antifold_seq"], ...)    # the caller scores them whole
    designs = pick_designs(scored, by="score_antifold_seq", regions=regions,
                           quota={5: 3, 6: 3, 7: 3, 8: 2, 9: 2, 10: 2})      # one entry per size the sets hold
    table = frontier(scored, by="score_antifold_seq")                        # best score at each mutation count

Scoring is not done here, on purpose. A set's score is not a function of its singles' scores, so nothing in this
module carries a single's score forward or adds two together; `propose_sets` drops every `score_` column, as
`select_pairs` does, and the sets are scored as whole candidates by `score_variants`. Of the scorers that can do
that, `antifold_seq` is the one that scores a candidate with several mutations at all. Six sizes at 300 sets each
come to 1,800 candidates, under the 2,000 that `score_variants` accepts for a per-variant scorer without a cap.

Two kinds of rule, in two places. A rule about one set (positions distinct, mutations not on top of each other,
no more than so many in one region) is enforced while the set is drawn. A rule about the designs taken together
(each mutation in so many designs, designs that differ enough, a quota per size) can only be enforced while they are
chosen, so it lives in `pick_designs`. Region labels come in as a plain mapping from (chain, position) to a label,
because the classes the plan uses (paratope, CDR-H3, second shell, framework) come from the position map and not
from anything in this package.
"""

from __future__ import annotations

import warnings
from collections import Counter
from pathlib import Path
from typing import TYPE_CHECKING, Callable, Iterable, Mapping, Sequence

import numpy as np

from .errors import MissingDependencyError
from .structures import chain_sequence
from .variants import (
    Mutation,
    VariantTable,
    _coordinates,
    _distance,
    _imgt_lookup,
    _rows_for,
    make_id,
    parse_id,
    write_table,
)

if TYPE_CHECKING:
    import pandas as pd

Regions = Mapping[tuple[str, int], str]
"""(chain, position) -> a label such as "H3" or "framework". A position with no entry has no label and is never limited."""

ATTEMPTS_PER_SET = 20
"""Random draws tried for each set wanted before a size is given up on. A draw that cannot reach the size, because
the pool is too crowded or a limit is hit, counts as an attempt."""


def propose_sets(
    variants: VariantTable,
    sizes: Iterable[int] = range(5, 11),
    per_size: int = 300,
    output_dir: Path | str | None = None,
    *,
    by: str,
    pool: int = 45,
    min_distance: float = 8.0,
    regions: Regions | None = None,
    max_per_set: Mapping[str, int] | None = None,
    skip: Iterable[str] = (),
    per_position: int | None = 1,
    bias: float = 0.0,
    seed: int = 0,
    structure: Path | str | None = None,
) -> VariantTable:
    """Draw candidate sets of several substitutions from the best singles. Returns them as a variant table.

    `by` names the score column that ranks the singles, so the ranking is the caller's and this function never
    combines scores. At most `per_position` singles are kept for each position, the best-scoring ones (`None` keeps
    them all), and the best `pool` of those are the shortlist. The pool sets how often a substitution can recur
    across the designs picked later: fifteen designs of about six substitutions are some 90 slots, and the plan wants
    each substitution in two to four designs, which asks for some 22 to 45 distinct substitutions. The default pool
    is 45, from a run on 1N8Z with real `antifold_seq` scores (the constant domain skipped, `per_position=1`, 8 Å,
    300 sets per size, no isoelectric point rule, 15 designs split equally over the sizes). A pool of 30 filled 13 of
    the 15 designs with one seed and 15 with another, and in that second run 7 of the 15 were size 5 instead of the
    even split; 2 and 4 substitutions were in fewer than two designs. A pool of 45 filled all 15 with the 3, 3, 3, 2,
    2, 2 split in both seeds, with 8 and 7 substitutions in fewer than two designs; a pool of 60 filled them too, with
    12. So 45 is the smallest pool tried that fills the designs reliably, and the price is coverage: the plan's
    minimum of two uses per substitution is reported, not enforced. That is two seeds, one structure and one scorer.
    Re-run it on your own scores.

    `per_position` matters on a real table. On 1N8Z chain B, with the table-route AntiFold score and the constant
    domain left out, the best 30 singles sit on 8 positions, one of them with ten substitutions in the pool. A set
    changes a position once, so no set of 9 or more can be drawn from that pool, and at the old 4 Å default nothing
    beyond size 5. With one substitution per position the same pool spans 30 positions, and each size from 5 to 10
    yields its full 100 sets at 4, 6 and 8 Å. A warning is raised when the pool has fewer positions than the
    largest size. Every set is a random greedy draw from the shortlist: the singles are put in a random order and
    taken in turn when they fit, until the set has its size. A draw that cannot reach the size is dropped, and a set
    seen twice is kept once, so the number returned per size can fall short of `per_size`; `config["found"]` says
    by how much.

    A set obeys these, checked as it grows:

    * one substitution per position;
    * no two substitutions closer than `min_distance` angstroms, CA to CA (8 by default), read off the structure
      the table carries (residues along the chain without one, as in `select_pairs`). Scoring a set whole fixes
      the arithmetic, not the physics: `antifold_seq` scores a mutant against the one fixed parent backbone, so it
      cannot see the rearrangement that two interacting changes cause, and keeping substitutions apart is the only
      guard against that here. The earlier default of 4 Å excluded only adjacent residues, because consecutive CAs
      are about 3.8 Å apart. On 1N8Z, 8 Å still gives all sizes from 5 to 10 their full 100 sets with
      `per_size=100`, `pool=30` and `per_position=1`; at 10 Å size 10 falls to 28. Raise it to be stricter, or pass 0
      to switch the rule off (a pair that cannot be measured is still dropped);
    * at most `max_per_set[label]` substitutions whose position has that label in `regions`, for example
      `max_per_set={"H3": 1}`. Positions with no label are not limited.

    `skip` names region labels whose singles are dropped before the pool is cut, for example `skip=("CH1",)` for
    the constant domain of a Fab structure, which an scFv does not have. The order is: drop the skipped, keep the
    best `pool` of the rest, then draw sets under `max_per_set`. Without it the pool is cut first, and on a Fab
    structure the best singles by a structure-based score can be mostly constant domain (93 of the top 100 on 1N8Z
    chain B with the table-route AntiFold score), which leaves a pool that cannot make a set of the size wanted.

    Inside the pool the draw is uniform by default: rank 1 and rank 30 are taken equally often, and the score only
    decides who is in the pool. `bias` changes that. A set is built by taking the pool in a random order, and with
    `bias > 0` the order is drawn with weights proportional to `rank ** -bias`, so higher-ranked singles come
    earlier and so recur more across the sets; `bias=0` is the uniform draw and gives the sets it always did. The
    greedy draw, which takes a single if it fits the ones already taken, also under-represents a single that
    conflicts with many others in the pool. `config["conflicts"]` counts, for each single in the pool, how many
    others it can never be in a set with (the same position, or closer than `min_distance`).

    The same seed gives the same sets. Score columns are dropped, because a set's score is not a function of its
    singles'; score the sets with `score_variants`.
    """
    import pandas as pd

    frame = variants.frame
    sizes = sorted({int(k) for k in sizes})
    if not sizes or min(sizes) < 2:
        raise ValueError(f"sizes must be whole numbers of at least 2, got {sizes}")
    per_candidate = frame.groupby("variant_id")["position"].size()
    multi = per_candidate[per_candidate > 1]
    if len(multi):
        raise ValueError(
            f"propose_sets takes single mutations and {len(multi)} of these candidates already have more than "
            f"one (e.g. {multi.index[0]!r}). Pass the singles table, or the gated singles."
        )
    if by not in frame.columns:
        raise KeyError(f"no column {by!r} to rank on; present: {', '.join(sorted(frame.columns))}")
    max_per_set = dict(max_per_set or {})
    skip = set(skip)
    if max_per_set and regions is None:
        raise ValueError("max_per_set needs regions: the labels the limits are counted over")
    if skip and regions is None:
        raise ValueError("skip needs regions: the labels it names")
    if bias < 0:
        raise ValueError(f"bias must not be negative, got {bias}")
    if per_position is not None and per_position < 1:
        raise ValueError(f"per_position must be at least 1 or None, got {per_position}")
    regions = dict(regions or {})
    with_score = frame[frame[by].notna()]
    skipped = 0
    if skip:
        keep = [regions.get((str(c), int(p))) not in skip for c, p in zip(with_score["chain"], with_score["position"])]
        skipped = len(keep) - sum(keep)
        with_score = with_score[keep]
    ordered = with_score.sort_values(by, ascending=False, kind="stable")
    if per_position is not None:
        ordered = ordered.groupby(["chain", "position"], sort=False).head(per_position)
    ranked = ordered.head(pool)
    if len(ranked) < max(sizes):
        raise ValueError(
            f"only {len(ranked)} singles have a {by!r} score"
            f"{f' and are outside skip={sorted(skip)}' if skip else ''}, fewer than the {max(sizes)} a set of that "
            f"size needs. Widen the singles, or lower the largest size."
        )
    positions_in_pool = len(set(zip(ranked["chain"], ranked["position"])))
    if positions_in_pool < max(sizes):
        warnings.warn(
            f"the pool of {len(ranked)} singles sits on {positions_in_pool} positions, fewer than the {max(sizes)} "
            f"a set of the largest size needs, and a set changes a position once, so those sizes cannot be drawn. "
            f"Keep fewer substitutions per position (per_position=1), or widen the pool.",
            RuntimeWarning,
            stacklevel=2,
        )

    target = Path(structure) if structure is not None else variants.structure
    singles = [
        (str(r.variant_id), str(r.chain), int(r.position), str(r.wt), str(r.mut)) for r in ranked.itertuples()
    ]
    labels = [regions.get((s[1], s[2])) for s in singles]
    coords = _coordinates(target, {s[1] for s in singles})
    n = len(singles)
    compatible = np.zeros((n, n), dtype=bool)
    for i in range(n):
        for j in range(i):
            if (singles[i][1], singles[i][2]) == (singles[j][1], singles[j][2]):
                continue  # two substitutions of one position are not a set
            distance = _distance(coords, singles[i], singles[j])
            ok = distance is not None and distance >= min_distance
            compatible[i, j] = compatible[j, i] = ok
    conflicts = {singles[i][0]: int(n - 1 - compatible[i].sum()) for i in range(n)}

    rng = np.random.default_rng(seed)
    rows: list[dict] = []
    found: dict[int, int] = {}
    label_lookup = _imgt_lookup(ranked)
    # `ranked` is best first, so index i has rank i + 1. A weighted order without replacement: each single gets the
    # key u ** (rank ** bias) for a uniform u, and the largest keys come first. With bias 0 the keys are all u, which
    # is a uniform order, so that case keeps the plain permutation and the sets it always drew.
    exponents = np.arange(1, n + 1, dtype=float) ** bias
    for size in sizes:
        drawn: dict[str, list[int]] = {}
        attempts = 0
        while len(drawn) < per_size and attempts < per_size * ATTEMPTS_PER_SET:
            attempts += 1
            chosen: list[int] = []
            in_label: Counter[str] = Counter()
            order = rng.permutation(n) if not bias else np.argsort(-(rng.random(n) ** exponents), kind="stable")
            for i in order:
                if not all(compatible[i, j] for j in chosen):
                    continue
                label = labels[i]
                if label in max_per_set and in_label[label] >= max_per_set[label]:
                    continue
                chosen.append(int(i))
                if label is not None:
                    in_label[label] += 1
                if len(chosen) == size:
                    break
            if len(chosen) == size:
                mutations = [(singles[i][1], singles[i][2], singles[i][3], singles[i][4]) for i in chosen]
                drawn.setdefault(make_id(mutations), chosen)
        found[size] = len(drawn)
        for chosen in drawn.values():
            mutations = [(singles[i][1], singles[i][2], singles[i][3], singles[i][4]) for i in chosen]
            for row in _rows_for(mutations, label_lookup):
                rows.append({**row, "n_mutations": size})

    if not rows:
        raise ValueError(
            f"no set of sizes {sizes} could be drawn from these {n} singles under min_distance={min_distance} "
            f"and max_per_set={max_per_set}. Relax a limit, lower min_distance, or widen the singles."
        )
    out = Path(output_dir) if output_dir is not None else variants.output_dir / "sets"
    return write_table(
        pd.DataFrame(rows),
        out / "sets.parquet",
        output_dir=variants.output_dir,
        structure=target,
        config={
            "sizes": sizes,
            "per_size": per_size,
            "found": found,
            "pool": n,
            "by": by,
            "min_distance": min_distance,
            "max_per_set": max_per_set,
            "skip": sorted(skip),
            "skipped": skipped,
            "per_position": per_position,
            "positions_in_pool": positions_in_pool,
            "bias": bias,
            "conflicts": conflicts,
            "seed": seed,
            "from": str(variants.path),
            "distance_units": "angstrom" if target is not None else "residues",
        },
    )


def pick_designs(
    scored: VariantTable,
    *,
    by: str,
    n: int = 15,
    min_differences: int = 3,
    max_uses: int = 4,
    min_uses: int = 2,
    regions: Regions | None = None,
    max_per_design: Mapping[str, int] | None = None,
    cdr_labels: Iterable[str] = (),
    quota: Mapping[int, int] | None = None,
    max_pi_shift: float | tuple[float, float] | None = 1.0,
    sequences: Mapping[str, str] | None = None,
    pi_function: Callable[[str], float] | None = None,
    output_dir: Path | str | None = None,
) -> VariantTable:
    """Choose `n` designs from scored sets, ranked within each size, under the rules that hold across designs.

    Scores are not comparable across sizes (every extra substitution costs score, and the set is scored whole), so
    the candidates are ranked within their size and the sizes are served in turn, not by one global order. Within
    a size, candidates go in order of `by`, highest first; ties go to fewer substitutions in `cdr_labels`, then to
    the id, so the same table always gives the same designs. The sizes then take turns, smallest first: each takes
    its best candidate that passes the rules below, then the next size does, round after round, until `n` designs
    are chosen or no size has a candidate left. All sizes draw on the same `max_per_design` and `max_uses` budget, so
    the smallest size has first call on it in each round; the designs are numbered in the order they were chosen.
    A candidate is kept when all of these hold:

    * it differs from every design already kept at `min_differences` positions or more, counting a position as
      different when the two designs put different residues there, the parent's included;
    * none of its substitutions is already in `max_uses` designs, and one whose label is in `max_per_design` is
      in no more designs than that, for example `max_per_design={"H3": 2}` for the plan's CDR-H3 rule;
    * its size still has room in `quota`, a mapping from number of substitutions to designs. Without a quota the
      sizes in the table share the designs equally, the smallest sizes taking the extra one when `n` does not
      divide, and a size with nothing left gives its turn to the others. With one, each size gets at most its
      entry, and a size named with 0 gets none. A size in the table that the quota leaves out gets none too, and a
      warning names it and the candidates that cannot be chosen;
    * its isoelectric point is within `max_pi_shift` of the parent's: a number for a symmetric bound, or
      `(down, up)` to allow a larger fall than rise, or the reverse. The point is computed on the chains joined
      in name order, from `sequences` (chain to sequence) or else from the structure the table carries. For a Fab
      structure that includes the constant domains, which an scFv does not have and which dilute the shift, so
      a warning is raised when the structure is used; pass the scFv's own chains in `sequences`. `pi_function`
      replaces the calculation; the default is Biopython's, imported when first needed. Pass
      `max_pi_shift=None` to skip the rule.

    The plan also wants each substitution in at least `min_uses` designs. A greedy pass cannot guarantee a lower
    bound, so it is reported rather than enforced: `config["underused"]` lists the substitutions in fewer designs,
    and a warning says so, as it does when fewer than `n` designs could be chosen. Scores of the chosen sets are
    kept as they are, because each is the set's own, scored whole.
    """
    import pandas as pd

    frame = scored.frame
    if by not in frame.columns:
        raise KeyError(f"no column {by!r} to rank on; present: {', '.join(sorted(frame.columns))}")
    regions = dict(regions or {})
    max_per_design = dict(max_per_design or {})
    cdr_labels = set(cdr_labels)
    if max_per_design and not regions:
        raise ValueError("max_per_design needs regions: the labels the limits are counted over")
    pi_bounds = _pi_bounds(max_pi_shift) if max_pi_shift is not None else None

    values = frame.groupby("variant_id", sort=False)[by].first()
    by_candidate_size: dict[int, list] = {}
    for variant_id, score in values.items():
        if pd.isna(score):
            continue  # a candidate the scorer could not score cannot be ranked
        mutations = parse_id(str(variant_id))
        cdr_changes = sum(regions.get((m[0], m[1])) in cdr_labels for m in mutations)
        by_candidate_size.setdefault(len(mutations), []).append((-float(score), cdr_changes, str(variant_id), mutations))
    for group in by_candidate_size.values():
        group.sort(key=lambda c: c[:3])
    sizes = sorted(by_candidate_size)
    scored_candidates = sum(len(group) for group in by_candidate_size.values())
    unscored = len(values) - scored_candidates

    if quota is not None:
        left_out = [size for size in sizes if size not in quota]
        if left_out:
            ignored = sum(len(by_candidate_size[size]) for size in left_out)
            warnings.warn(
                f"quota names no slots for size {', '.join(map(str, left_out))}, so the {ignored} scored candidates of "
                f"{'that size' if len(left_out) == 1 else 'those sizes'} cannot be chosen. Name a size with 0 to leave "
                f"it out on purpose.",
                RuntimeWarning,
                stacklevel=2,
            )

    shift_of: Callable[[Sequence[Mutation]], float] | None = None
    if pi_bounds is not None:
        shift_of = _pi_shift_function(scored, sequences, pi_function)

    chosen: list[tuple[str, Sequence[Mutation]]] = []
    uses: Counter[Mutation] = Counter()
    by_size: Counter[int] = Counter()
    cursor = {size: 0 for size in sizes}  # a rule only ever gets stricter, so a rejected candidate stays rejected
    while len(chosen) < n:
        progressed = False
        for size in sizes:
            if len(chosen) >= n:
                break
            if quota is not None and by_size[size] >= quota.get(size, 0):
                continue
            group = by_candidate_size[size]
            while cursor[size] < len(group):
                _negative, _cdr, variant_id, mutations = group[cursor[size]]
                cursor[size] += 1
                if any(_hamming(mutations, other) < min_differences for _, other in chosen):
                    continue
                if any(
                    uses[m] >= min(max_uses, max_per_design.get(regions.get((m[0], m[1])), max_uses)) for m in mutations
                ):
                    continue
                if shift_of is not None and not -pi_bounds[0] <= shift_of(mutations) <= pi_bounds[1]:
                    continue
                chosen.append((variant_id, mutations))
                by_size[size] += 1
                uses.update(mutations)
                progressed = True
                break
        if not progressed:
            break

    underused = sorted(m for m, count in uses.items() if count < min_uses)
    if len(chosen) < n:
        warnings.warn(
            f"chose {len(chosen)} of {n} designs: the rules left no more candidates among {scored_candidates} scored. "
            f"Propose more sets, relax min_differences or max_uses, or widen the quota.",
            RuntimeWarning,
            stacklevel=2,
        )
    if underused:
        warnings.warn(
            f"{len(underused)} substitutions are in fewer than {min_uses} designs (e.g. {make_id([underused[0]])}); "
            f"see config['underused']. The greedy draw under-represents a substitution that conflicts with many "
            f"others in the pool, so look at config['conflicts'] of the table propose_sets wrote before widening "
            f"the pool or lowering max_uses.",
            RuntimeWarning,
            stacklevel=2,
        )
    if not chosen:
        raise ValueError(f"no design could be chosen from {scored_candidates} scored candidates ({unscored} unscored)")

    rank = {variant_id: i + 1 for i, (variant_id, _) in enumerate(chosen)}
    picked = frame[frame["variant_id"].isin(rank)].copy()
    picked["design"] = picked["variant_id"].map(rank)
    picked["n_mutations"] = picked["variant_id"].map(picked.groupby("variant_id")["position"].size())
    picked = picked.sort_values(["design"], kind="stable")
    out = Path(output_dir) if output_dir is not None else scored.output_dir / "designs"
    return write_table(
        picked,
        out / "designs.parquet",
        output_dir=scored.output_dir,
        structure=scored.structure,
        config={
            "n": n,
            "chosen": len(chosen),
            "by": by,
            "min_differences": min_differences,
            "max_uses": max_uses,
            "min_uses": min_uses,
            "max_per_design": max_per_design,
            "quota": dict(quota) if quota is not None else None,
            "by_size": dict(sorted(by_size.items())),
            "max_pi_shift": list(pi_bounds) if pi_bounds is not None and isinstance(max_pi_shift, (tuple, list)) else max_pi_shift,
            "unscored": unscored,
            "underused": [make_id([m]) for m in underused],
            "uses": {make_id([m]): c for m, c in sorted(uses.items())},
            "from": str(scored.path),
        },
    )


def frontier(scored: VariantTable, by: str, keep: int = 1) -> "pd.DataFrame":
    """The best `keep` candidates at each number of substitutions, for the jury's score-against-count plot.

    One row per kept candidate, with `n_mutations`, `variant_id` and the `by` column, ordered by size and then by
    score. This is the plan's "best-scoring sets for each mutation count"; reading where the curve flattens is how
    the mutation-count quota is set. A candidate with no score is left out.
    """
    frame = scored.frame
    if by not in frame.columns:
        raise KeyError(f"no column {by!r} to rank on; present: {', '.join(sorted(frame.columns))}")
    per_candidate = frame.groupby("variant_id", sort=False).agg(n_mutations=("position", "size"), **{by: (by, "first")})
    per_candidate = per_candidate.dropna(subset=[by]).reset_index()
    ordered = per_candidate.sort_values(["n_mutations", by], ascending=[True, False], kind="stable")
    return ordered.groupby("n_mutations", sort=True).head(keep).reset_index(drop=True)[["n_mutations", "variant_id", by]]


def _hamming(left: Sequence[Mutation], right: Sequence[Mutation]) -> int:
    """Positions at which two candidates put different residues, with the parent's residue where one has no change."""
    a = {(c, p): m for c, p, _wt, m in left}
    b = {(c, p): m for c, p, _wt, m in right}
    return sum(a.get(position) != b.get(position) for position in a.keys() | b.keys())


def _pi_bounds(limit: float | tuple[float, float]) -> tuple[float, float]:
    """(largest fall, largest rise) the isoelectric point may make, from a single number or a `(down, up)` pair."""
    if isinstance(limit, (tuple, list)):
        if len(limit) != 2:
            raise ValueError(f"max_pi_shift takes a number or (down, up), got {limit!r}")
        down, up = float(limit[0]), float(limit[1])
    else:
        down = up = float(limit)
    if down < 0 or up < 0:
        raise ValueError(f"max_pi_shift must not be negative, got {limit!r}")
    return down, up


def _pi_shift_function(
    scored: VariantTable, sequences: Mapping[str, str] | None, pi_function: Callable[[str], float] | None
) -> Callable[[Sequence[Mutation]], float]:
    """A function from a candidate's mutations to the signed shift of the isoelectric point from the parent's."""
    calculate = pi_function or _biopython_pi
    chains = dict(sequences) if sequences is not None else {}
    if not chains:
        if scored.structure is None:
            raise ValueError(
                "the isoelectric point rule needs the chain sequences: pass sequences={chain: sequence}, "
                "or a table that carries its structure. max_pi_shift=None skips the rule."
            )
        wanted = sorted(set(scored.frame["chain"].astype(str)))
        chains = {chain: chain_sequence(scored.structure, chain) for chain in wanted}
        lengths = ", ".join(f"{chain}: {len(sequence)} residues" for chain, sequence in chains.items())
        warnings.warn(
            f"the isoelectric point is computed on the whole chains of the structure ({lengths}). For a Fab that "
            f"includes the constant domains, which an scFv does not have and which dilute the shift; pass the "
            f"scFv's own chains as sequences={{chain: sequence}}.",
            RuntimeWarning,
            stacklevel=3,  # the caller of pick_designs
        )
    order = sorted(chains)
    parent = calculate("".join(chains[c] for c in order))

    def shift(mutations: Sequence[Mutation]) -> float:
        residues = {c: list(chains[c]) for c in order}
        for chain, position, wt, mut in mutations:
            if chain not in residues or not 0 < position <= len(residues[chain]) or residues[chain][position - 1] != wt:
                raise ValueError(f"{make_id([(chain, position, wt, mut)])} does not match the sequence of chain {chain}")
            residues[chain][position - 1] = mut
        return calculate("".join("".join(residues[c]) for c in order)) - parent

    return shift


def _biopython_pi(sequence: str) -> float:
    try:
        from Bio.SeqUtils.ProtParam import ProteinAnalysis
    except ModuleNotFoundError as exc:
        raise MissingDependencyError(
            "the isoelectric point rule uses Biopython, which is not installed here: `pip install biopython`, "
            "or pass max_pi_shift=None to skip the rule, or pi_function= to bring your own."
        ) from exc
    return float(ProteinAnalysis(sequence).isoelectric_point())
