"""The one record every step speaks in, and the three things you do to a list of them.

Collections are plain `list[Design]`. Nothing wraps them.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from pathlib import Path
from typing import TYPE_CHECKING, Any, Callable, Iterable, Sequence

if TYPE_CHECKING:  # pandas is an analysis-only extra, never imported at runtime here
    import pandas as pd


class Scores(dict[str, float]):
    """A dict whose miss tells you what the model actually produced.

    Metrics differ by model (ColabFold: plddt/ptm/pae, Boltz adds affinity,
    ProteinMPNN gives a log-probability), so this stays untyped — but a typo in
    a metric name should not look like a missing design.
    """

    def __missing__(self, key: str) -> float:
        present = ", ".join(sorted(self)) or "none"
        raise KeyError(f"no score {key!r} on this design; present: {present}")


@dataclass(frozen=True)
class Design:
    id: str
    path: Path
    parent_id: str | None
    """The design this came from, across steps. Load-bearing: scoring pairs a
    prediction with the backbone it descends from through this."""
    producer: str
    scores: dict[str, float] = field(default_factory=Scores)
    meta: dict[str, Any] = field(default_factory=dict)
    """Round number, seed, ancestor ids, free-form."""

    def __post_init__(self) -> None:
        # Coerce so an adapter can hand us a plain dict and lookups still
        # explain themselves. Frozen, hence object.__setattr__.
        if not isinstance(self.scores, Scores):
            object.__setattr__(self, "scores", Scores(self.scores))
        if not isinstance(self.path, Path):
            object.__setattr__(self, "path", Path(self.path))


def lineage(designs: Sequence[Design], design_id: str) -> list[Design]:
    """The parent chain, oldest first, ending at `design_id`.

    Only follows links present in `designs`; a chain that leaves the pool stops
    there rather than guessing.
    """
    index = {d.id: d for d in designs}
    if design_id not in index:
        raise KeyError(f"no design {design_id!r} in this pool of {len(index)}")
    chain: list[Design] = []
    seen: set[str] = set()
    current: str | None = design_id
    while current is not None and current in index and current not in seen:
        seen.add(current)
        design = index[current]
        chain.append(design)
        current = design.parent_id
    chain.reverse()
    return chain


def select(
    designs: Sequence[Design],
    top_k: int,
    by: str,
    ascending: bool = True,
    where: Callable[[Design], bool] | None = None,
) -> list[Design]:
    """Sort by a score and truncate.

    Deterministic: stable sort with an explicit tiebreak on `id`, because a
    nondeterministic selection shifts every downstream cache key.
    """
    pool = [d for d in designs if where is None or where(d)]
    missing = [d.id for d in pool if by not in d.scores]
    if missing:
        present = sorted({k for d in pool for k in d.scores})
        raise KeyError(
            f"{len(missing)} design(s) have no score {by!r} "
            f"(e.g. {missing[0]!r}); present across this pool: "
            f"{', '.join(present) or 'none'}"
        )
    # Two stable sorts rather than a tuple key: reverse=True on a tuple would
    # also reverse the id tiebreak, so ties would reorder with the direction.
    ordered = sorted(pool, key=lambda d: d.id)
    ordered.sort(key=lambda d: d.scores[by], reverse=not ascending)
    return ordered[:top_k]


def to_dataframe(designs: Iterable[Design]) -> "pd.DataFrame":
    """One row per design, `scores` flattened into columns.

    pandas is an extra (`design-pipeline[analysis]`): analysis wants it, the pipeline
    does not, so it is imported here rather than at module import.
    """
    try:
        import pandas as pd
    except ModuleNotFoundError as exc:  # pragma: no cover - exercised by hand
        raise ModuleNotFoundError(
            "to_dataframe needs pandas. Install it with "
            "`uv pip install 'design-pipeline[analysis]'`."
        ) from exc

    rows = []
    for d in designs:
        row: dict[str, Any] = {
            "id": d.id,
            "parent_id": d.parent_id,
            "producer": d.producer,
            "path": str(d.path),
        }
        row.update(d.scores)
        # meta is free-form, so it is prefixed rather than merged into scores.
        row.update({f"meta_{k}": v for k, v in d.meta.items()})
        rows.append(row)
    return pd.DataFrame(rows)
