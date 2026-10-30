"""Content-addressed run directories.

The key is the resolved config, the seed, and the content hashes of the input
files. Chaining those content hashes across rounds gives resume for free,
provided `select` stays deterministic.

The seed is in the key on purpose: these models are stochastic, and without it
"another 100 designs with the same config" silently returns the cached 100.

The DONE marker holds the design records, not just the paths, because `scores`
and `parent_id` are parsed from model output that a resumed run never sees.
"""

from __future__ import annotations

import hashlib
import json
from pathlib import Path
from typing import Any, Iterable, Mapping, Sequence

from .designs import Design

DONE = "DONE"


def _hash_file(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as fh:
        for chunk in iter(lambda: fh.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def cache_key(
    config: Mapping[str, Any],
    inputs: Sequence[Design],
    files: Iterable[Path],
    seed: int | None,
) -> str:
    payload = {
        "config": json.loads(json.dumps(config, default=str, sort_keys=True)),
        "seed": seed,
        # Design id as well as content: two designs from different rounds can
        # hold byte-identical files (the same sequence rediscovered, say) and
        # different ancestry, and a hit must not hand back the other one's
        # parent_id. Ids are derived from the producing run's key, so they are
        # deterministic and a rerun still resumes.
        "inputs": sorted((d.id, _hash_file(d.path)) for d in inputs),
        # Name included: many models key their output names off the input name.
        "files": sorted((p.name, _hash_file(p)) for p in files),
    }
    return hashlib.sha256(json.dumps(payload, sort_keys=True).encode()).hexdigest()[:16]


def read_manifest(run_dir: Path) -> list[Design] | None:
    """The designs of a finished run, or None if it is not finished."""
    marker = run_dir / DONE
    if not marker.is_file():
        return None
    try:
        records = json.loads(marker.read_text())
    except json.JSONDecodeError:
        return None
    designs = [
        Design(
            id=r["id"],
            path=run_dir / r["path"],
            parent_id=r["parent_id"],
            producer=r["producer"],
            scores=r["scores"],
            meta=r["meta"],
        )
        for r in records
    ]
    # A marker pointing at files somebody deleted is a miss, not a hit.
    if not designs or not all(d.path.exists() for d in designs):
        return None
    return designs


def write_manifest(run_dir: Path, designs: Sequence[Design]) -> None:
    records = [
        {
            "id": d.id,
            "path": d.path.relative_to(run_dir).as_posix(),
            "parent_id": d.parent_id,
            "producer": d.producer,
            "scores": dict(d.scores),
            "meta": d.meta,
        }
        for d in designs
    ]
    (run_dir / DONE).write_text(json.dumps(records, indent=1, default=str))


def tool_key(
    *,
    structure_hash: str,
    candidate_hash: str,
    model_name: str,
    resolved: Mapping[str, Any],
) -> str:
    """The key for one scorer or one gate check over one candidate set.

    Per scorer, not per pipeline step. Adding GearBind to a `models=[...]` list
    must not invalidate the AntiFold run that already finished - the harness
    re-scores dozens of held-out systems while the scorer set keeps changing,
    and a key that spanned the whole fan-out would throw all of it away every
    time somebody added a column.

    The candidate hash covers the ids actually handed to this adapter, so a
    `max_variants` cap is part of what identifies the run rather than a
    footgun that returns a short result under the full set's key.
    """
    payload = {
        "structure": structure_hash,
        "candidates": candidate_hash,
        "model_name": model_name,
        "resolved": json.loads(json.dumps(resolved, default=str, sort_keys=True)),
    }
    return hashlib.sha256(json.dumps(payload, sort_keys=True).encode()).hexdigest()[:16]


def hash_file(path: Path) -> str:
    """Content hash of one file. The structure half of `tool_key`."""
    return _hash_file(Path(path))


def hash_ids(ids: Iterable[str]) -> str:
    """Content hash of a candidate set. Order-independent, duplicates collapsed:
    the same set asked for twice is one run."""
    joined = "\n".join(sorted(set(ids)))
    return hashlib.sha256(joined.encode()).hexdigest()[:16]
