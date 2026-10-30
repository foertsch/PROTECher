"""The variant table: one row per mutation, one `variant_id` per candidate.

This is the second artifact type, alongside `Complex`. Parquet on disk, six
columns:

    variant_id  str   stable join key across every scorer
    chain       str   H or L
    position    int   1-based index into the chain as given
    imgt        str   IMGT label from ANARCI, for cross-antibody transfer
    wt          str   original residue, one letter
    mut         str   replacement residue, one letter

A pair is two rows sharing one `variant_id`. Same schema, no special case, no
second table type - which is what lets `score_variants` be called on singles
and on pairs with no branch in between.

`enumerate_variants` and `select_pairs` are plain functions: a loop and a
distance filter. Nothing here runs a model, so nothing here is an adapter.
"""

from __future__ import annotations

import re
from dataclasses import dataclass, field
from pathlib import Path
from typing import TYPE_CHECKING, Any, Iterable, Mapping, Sequence

import numpy as np

from . import numbering
from .structures import AMINO_ACIDS, Complex, ca_coordinates, chain_sequence

if TYPE_CHECKING:
    import pandas as pd

SCHEMA: tuple[str, ...] = ("variant_id", "chain", "position", "imgt", "wt", "mut")
"""The columns every variant table has. Scorers append `score_<model>`
alongside; nothing ever removes one of these."""

_TOKEN = re.compile(r"^(?P<chain>[^:+]+):(?P<wt>[A-Z])(?P<position>\d+)(?P<mut>[A-Z])$")

Mutation = tuple[str, int, str, str]
"""(chain, position, wt, mut). Position is 1-based into the chain as given."""


def make_id(mutations: Iterable[Mutation]) -> str:
    """The canonical `variant_id` for a candidate: H:T30Y, or H:T30Y+L:S91A.

    Sorted, so the same pair discovered in either order is one candidate and
    one cache entry. The id carries every mutation, which is what makes it a
    join key a scorer can be handed on its own.
    """
    tokens = sorted(
        (chain, int(position), wt, mut) for chain, position, wt, mut in mutations
    )
    if not tokens:
        raise ValueError("a candidate needs at least one mutation")
    return "+".join(f"{c}:{wt}{p}{mut}" for c, p, wt, mut in tokens)


def parse_id(variant_id: str) -> tuple[Mutation, ...]:
    """The mutations a `variant_id` names. Raises on anything malformed."""
    mutations: list[Mutation] = []
    for token in variant_id.split("+"):
        match = _TOKEN.match(token.strip())
        if match is None:
            raise ValueError(
                f"{variant_id!r} is not a variant id: expected tokens like "
                f"H:T30Y joined with +, got {token!r}"
            )
        mutations.append(
            (match["chain"], int(match["position"]), match["wt"], match["mut"])
        )
    return tuple(mutations)


@dataclass(frozen=True)
class VariantTable:
    """A parquet file of variant rows, plus where it came from.

    `rows` and `candidates` differ once pairs are in play: 40 rows, 20
    candidates. Both are recorded here so a repr and a sanity check cost no
    disk read; `frame` is the one thing that does.

    `structure` is the complex the positions index. It is carried along every
    derived table because `gate` and `select_pairs` need the chain sequence and
    the coordinates, and the API hands them only the table.
    """

    path: Path
    output_dir: Path
    rows: int
    candidates: int
    columns: tuple[str, ...]
    structure: Path | None = None
    config: Mapping[str, Any] = field(default_factory=dict)

    def __post_init__(self) -> None:
        object.__setattr__(self, "path", Path(self.path))
        object.__setattr__(self, "output_dir", Path(self.output_dir))
        if self.structure is not None:
            object.__setattr__(self, "structure", Path(self.structure))
        object.__setattr__(self, "columns", tuple(self.columns))

    @property
    def frame(self) -> "pd.DataFrame":
        """Read the table. Not cached: the file on disk is the state."""
        return read_table(self.path)

    @property
    def score_columns(self) -> tuple[str, ...]:
        """`score_<model>` columns, the raw ones excluded.

        What the caller ranks on. There is deliberately no combined column and
        no default weighting: which of these to trust comes out of the
        validation harness, not out of this package.
        """
        return tuple(
            c for c in self.columns if c.startswith("score_") and not c.endswith("_raw")
        )

    def __len__(self) -> int:
        return self.rows


def read_table(path: Path) -> "pd.DataFrame":
    import pandas as pd

    frame = pd.read_parquet(path)
    missing = [c for c in SCHEMA if c not in frame.columns]
    if missing:
        raise ValueError(f"{path} is not a variant table: missing {missing}")
    return frame


def write_table(
    frame: "pd.DataFrame",
    path: Path,
    *,
    output_dir: Path,
    structure: Path | None,
    config: Mapping[str, Any],
) -> VariantTable:
    """Write the parquet and return the handle to it.

    Column order is fixed - the schema first, then whatever was appended - so a
    diff between two runs of the harness is readable.
    """
    missing = [c for c in SCHEMA if c not in frame.columns]
    if missing:
        raise ValueError(f"cannot write a variant table without {missing}")
    ordered = list(SCHEMA) + [c for c in frame.columns if c not in SCHEMA]
    frame = frame[ordered].reset_index(drop=True)
    path.parent.mkdir(parents=True, exist_ok=True)
    frame.to_parquet(path, index=False)
    return VariantTable(
        path=path,
        output_dir=output_dir,
        rows=len(frame),
        candidates=int(frame["variant_id"].nunique()),
        columns=tuple(frame.columns),
        structure=structure,
        config=config,
    )


def _rows_for(
    mutations: Sequence[Mutation],
    labels: Mapping[str, Sequence[str]],
) -> list[dict[str, Any]]:
    """The table rows for one candidate, one per mutation, sharing an id."""
    variant_id = make_id(mutations)
    rows = []
    for chain, position, wt, mut in sorted(mutations):
        chain_labels = labels.get(chain) or ()
        rows.append(
            {
                "variant_id": variant_id,
                "chain": chain,
                "position": int(position),
                "imgt": (
                    chain_labels[position - 1]
                    if 0 < position <= len(chain_labels)
                    else numbering.MISSING
                ),
                "wt": wt,
                "mut": mut,
            }
        )
    return rows


def enumerate_variants(
    complex_: Complex,
    chain: str,
    output_dir: Path | str | None = None,
    *,
    positions: Iterable[int] | None = None,
    mutations: str | Iterable[str] | None = None,
    only: Iterable[str] | None = None,
    structure: Path | str | None = None,
    scheme: str = "imgt",
) -> VariantTable:
    """Every single mutation of `chain`, as a variant table.

    The default is the full scan: every position of the chain against the other
    19 amino acids, so a 242-residue heavy chain gives 4,598 candidates.
    Narrow it with `positions` (1-based, into the chain as given) or
    `mutations` (the replacement alphabet).

    `only` restricts the table to exactly the candidates whose ids are listed,
    which is how the validation harness scores a measured set - and it is the
    one place multi-mutation ids are accepted here, because published
    ground-truth sets contain double mutants. Every id is checked against the
    structure, so a wt residue that disagrees is an error rather than a row
    that silently scores the wrong thing.

    `structure` picks a seed out of a predicted ensemble; the default is the
    complex's representative.
    """
    target = Path(structure) if structure is not None else complex_.path
    out = (
        Path(output_dir)
        if output_dir is not None
        else complex_.output_dir / "variants"
    )

    sequences = {chain: chain_sequence(target, chain)}
    labels: dict[str, Sequence[str]] = {
        chain: numbering.imgt_labels(sequences[chain], scheme=scheme)
    }

    rows: list[dict[str, Any]] = []
    if only is not None:
        for variant_id in dict.fromkeys(only):  # de-duplicated, order kept
            parsed = parse_id(variant_id)
            for mut_chain, position, wt, _mut in parsed:
                if mut_chain not in sequences:
                    sequences[mut_chain] = chain_sequence(target, mut_chain)
                    labels[mut_chain] = numbering.imgt_labels(
                        sequences[mut_chain], scheme=scheme
                    )
                _check_wt(target, mut_chain, position, wt, sequences[mut_chain])
            rows += _rows_for(parsed, labels)
        name = "requested"
        wanted: list[int] = []
    else:
        sequence = sequences[chain]
        wanted = (
            [int(p) for p in positions]
            if positions is not None
            else list(range(1, len(sequence) + 1))
        )
        alphabet = tuple(mutations) if mutations is not None else tuple(AMINO_ACIDS)
        for position in wanted:
            if not 0 < position <= len(sequence):
                raise IndexError(
                    f"position {position} is outside chain {chain!r}, which has "
                    f"{len(sequence)} residues"
                )
            wt = sequence[position - 1]
            for mut in alphabet:
                # The wt residue itself is not a variant. An X - a residue this
                # package could not read - is not a position we are willing to
                # mutate, because we would not know what we were replacing.
                if mut == wt or wt == "X":
                    continue
                rows += _rows_for([(chain, position, wt, mut)], labels)
        name = f"{chain}_singles"

    import pandas as pd

    frame = pd.DataFrame(rows, columns=list(SCHEMA))
    if frame.empty:
        raise ValueError(
            "no variants to enumerate: check `positions`, `mutations` and `only`"
        )
    return write_table(
        frame,
        out / f"{name}.parquet",
        output_dir=complex_.output_dir,
        structure=target,
        config={
            "chain": chain,
            "scheme": scheme,
            "positions": sorted(wanted) if positions is not None else None,
            "mutations": None if mutations is None else list(mutations),
            "only": None if only is None else int(frame["variant_id"].nunique()),
        },
    )


def _check_wt(
    structure: Path, chain: str, position: int, wt: str, sequence: str
) -> None:
    if not 0 < position <= len(sequence):
        raise IndexError(
            f"position {position} is outside chain {chain!r} of {structure}, "
            f"which has {len(sequence)} residues. Ground-truth sets are often "
            f"numbered in IMGT or in the PDB's own numbering; variant ids here "
            f"index the chain as given, 1-based."
        )
    found = sequence[position - 1]
    if found != wt:
        raise ValueError(
            f"{chain}:{wt}{position} does not match {structure}, which has "
            f"{found} at position {position} of chain {chain}. Either the "
            f"numbering differs or this is the wrong chain."
        )


def mutant_sequence(
    structure: Path, mutations: Sequence[Mutation], chain: str
) -> tuple[str, str]:
    """(wild-type, mutant) sequence of `chain` with every mutation applied.

    Every mutation of the candidate, not one of them: a scorer that needs the
    mutant sequence needs all of it. This is also the only honest way to score
    a pair - one pass over the full mutant - which is why nothing in this
    package adds two single scores together.
    """
    wild = chain_sequence(structure, chain)
    residues = list(wild)
    for mut_chain, position, wt, mut in mutations:
        if mut_chain != chain:
            continue
        _check_wt(structure, chain, position, wt, wild)
        residues[position - 1] = mut
    return wild, "".join(residues)


def select_pairs(
    variants: VariantTable,
    max_pairs: int = 20,
    min_distance: float = 10.0,
    output_dir: Path | str | None = None,
    *,
    by: str | None = None,
    structure: Path | str | None = None,
) -> VariantTable:
    """Combine single mutations into pairs that are far enough apart.

    `min_distance` is in angstroms, CA to CA, read off the structure the table
    carries - across chains as well as within one. Without a structure it falls
    back to separation along the chain in residues, and a cross-chain pair then
    has no distance at all and is kept; the units change with it, which is why
    the structure is the normal path.

    Candidates are taken in the table's existing row order, so ranking stays
    the caller's business: sort the survivors on the column you trust, or pass
    `by` to sort on one score column here. Pairs are then generated in
    triangular order - (1,2), (1,3), (2,3), (1,4) ... - so `max_pairs` spread
    across the top of the list instead of all containing the single best one.

    Score columns are dropped. A pair's score is not a function of its
    singles', and carrying the singles' numbers forward under a `score_` name
    is exactly the additivity shortcut this package refuses to take. Score the
    pairs instead: `score_variants(complex_, pairs, models=[...])`.
    """
    import pandas as pd

    frame = variants.frame
    sizes = frame.groupby("variant_id")["position"].size()
    multi = sizes[sizes > 1]
    if len(multi):
        raise ValueError(
            f"select_pairs takes single mutations and {len(multi)} of these "
            f"candidates already have more than one (e.g. {multi.index[0]!r}). "
            "Pass the singles table, or the gated singles."
        )
    if by is not None:
        if by not in frame.columns:
            raise KeyError(
                f"no column {by!r} to sort on; present: "
                f"{', '.join(sorted(frame.columns))}"
            )
        frame = frame.sort_values(by, ascending=False, kind="stable")

    target = Path(structure) if structure is not None else variants.structure
    singles = [
        (str(row.variant_id), str(row.chain), int(row.position), str(row.wt), str(row.mut))
        for row in frame.itertuples()
    ]
    coords = _coordinates(target, {chain for _, chain, _, _, _ in singles})
    labels = _imgt_lookup(frame)

    # Triangular order: every pair whose larger index is j comes before any
    # pair whose larger index is j + 1.
    order = sorted(
        ((i, j) for j in range(len(singles)) for i in range(j)),
        key=lambda pair: (pair[1], pair[0]),
    )

    rows: list[dict[str, Any]] = []
    kept = 0
    for i, j in order:
        if kept >= max_pairs:
            break
        left, right = singles[i], singles[j]
        if (left[1], left[2]) == (right[1], right[2]):
            continue  # two substitutions of one position are not a pair
        distance = _distance(coords, left, right)
        if distance is None or distance < min_distance:
            continue
        mutations = [
            (left[1], left[2], left[3], left[4]),
            (right[1], right[2], right[3], right[4]),
        ]
        for row in _rows_for(mutations, labels):
            rows.append({**row, "distance": distance})
        kept += 1

    if not rows:
        raise ValueError(
            f"no pair of these {len(singles)} singles is at least "
            f"{min_distance} apart. Lower min_distance, or widen the set of "
            f"positions the singles come from."
        )
    out = Path(output_dir) if output_dir is not None else variants.output_dir / "pairs"
    return write_table(
        pd.DataFrame(rows),
        out / "pairs.parquet",
        output_dir=variants.output_dir,
        structure=target,
        config={
            "max_pairs": max_pairs,
            "min_distance": min_distance,
            "by": by,
            "from": str(variants.path),
            "distance_units": "angstrom" if target is not None else "residues",
        },
    )


def _imgt_lookup(frame: "pd.DataFrame") -> Mapping[str, Sequence[str]]:
    """chain -> labels indexed by position - 1, built from the rows we have.

    Sparse: positions nobody enumerated stay empty. Enough, because a pair is
    built only out of positions that are in the table.
    """
    lookup: dict[str, list[str]] = {}
    for row in frame.itertuples():
        labels = lookup.setdefault(str(row.chain), [])
        while len(labels) < int(row.position):
            labels.append(numbering.MISSING)
        labels[int(row.position) - 1] = str(row.imgt)
    return lookup


def _coordinates(
    structure: Path | None, chains: Iterable[str]
) -> Mapping[str, np.ndarray] | None:
    if structure is None:
        return None
    return {chain: ca_coordinates(structure, chain) for chain in chains}


def _distance(
    coords: Mapping[str, np.ndarray] | None,
    left: tuple[str, str, int, str, str],
    right: tuple[str, str, int, str, str],
) -> float | None:
    """Angstroms between the two mutated positions, or None if unmeasurable.

    None means "cannot tell" - an unmodelled CA, a chain the structure does not
    have - and a pair we cannot measure is dropped rather than assumed to be
    far apart.
    """
    if coords is None:
        if left[1] != right[1]:
            # No structure and two different chains: nothing comparable to
            # measure, so the pair is kept and the filter does not apply.
            return float("inf")
        return float(abs(left[2] - right[2]))
    try:
        a = coords[left[1]][left[2] - 1]
        b = coords[right[1]][right[2] - 1]
    except (KeyError, IndexError):
        return None
    if not np.isfinite(a).all() or not np.isfinite(b).all():
        return None
    return float(np.linalg.norm(a - b))
