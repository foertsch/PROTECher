"""IMGT position labels, via ANARCI, degrading to nothing in particular.

IMGT labels are what lets a result transfer between antibodies: "H3 position
107" means the same thing in two different molecules, where "position 103 of
this chain" does not. Within one campaign they are decoration - the variant
table joins on `variant_id`, which is built from chain positions - so when
ANARCI is not installed the pipeline runs without them and says so once.

ANARCI needs HMMER, which makes it awkward to install on a laptop and awkward
in CI, so it is an extra (`design-pipeline[numbering]`) rather than a core
dependency.
"""

from __future__ import annotations

import warnings
from functools import lru_cache

MISSING = ""
"""The `imgt` value for a position ANARCI did not label. Empty rather than a
sentinel string: a gate or a group-by on it should look empty, not meaningful."""


@lru_cache(maxsize=1)
def available() -> bool:
    try:
        import anarci  # noqa: F401
    except Exception:  # pragma: no cover - depends on the environment
        return False
    return True


@lru_cache(maxsize=1)
def _warn_once() -> None:  # pragma: no cover - depends on the environment
    warnings.warn(
        "anarci is not installed, so the variant table's `imgt` column will be "
        "empty and results will not carry across antibodies. Install it with "
        "`uv pip install 'design-pipeline[numbering]'` (it needs HMMER on "
        "PATH). Everything else works: joins are on chain positions.",
        RuntimeWarning,
        stacklevel=2,
    )


def imgt_labels(sequence: str, scheme: str = "imgt") -> tuple[str, ...]:
    """One label per residue of `sequence`, `MISSING` where there is none.

    Residues outside the variable domain - a Fab's constant region, a linker,
    the antigen if it is handed in by mistake - have no IMGT number and get
    `MISSING`. So does every residue when ANARCI is absent.
    """
    labels = [MISSING] * len(sequence)
    if not available():
        _warn_once()
        return tuple(labels)

    from anarci import run_anarci

    # run_anarci rather than number(): it reports where in the input the
    # variable domain starts and ends, which is what aligns the numbering back
    # onto our positions. number() returns the domain alone.
    _, numbered, _, _ = run_anarci([("query", sequence)], scheme=scheme)
    for domain in numbered[0] or ():
        annotation, start, _end = domain
        index = start
        for (position, insertion), residue in annotation:
            # ANARCI pads to the scheme's full length with '-'; those are
            # positions the scheme has and this antibody does not.
            if residue == "-":
                continue
            if index < len(labels):
                labels[index] = f"{position}{insertion.strip()}"
            index += 1
    return tuple(labels)
