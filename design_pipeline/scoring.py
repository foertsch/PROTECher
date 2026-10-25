"""Scoring: the only place in this package where research judgment lives.

Pure numpy + biotite. No model, no GPU, no network.

Adding a metric means writing one function of (prediction, reference) and
putting it in METRICS. A metric returns None when it does not apply to a
design — a monomer has no interface PAE — and that simply means no score.
"""

from __future__ import annotations

import json
from dataclasses import replace
from pathlib import Path
from typing import Callable, Mapping, Sequence

import numpy as np

from .designs import Design, Scores

Metric = Callable[[Design, "Design | None"], float | None]


def _ca_atoms(path: Path):
    """CA atoms of the first model, with B-factors."""
    from biotite.structure.io.pdb import PDBFile

    structure = PDBFile.read(str(path)).get_structure(
        model=1, extra_fields=["b_factor"]
    )
    return structure[structure.atom_name == "CA"]


def sc_rmsd(prediction: Design, reference: Design | None) -> float | None:
    """Self-consistency RMSD: predicted structure against the backbone it came
    from, CA atoms, after superposition.

    Returns None when the two have different CA counts. Mapping unequal lengths
    needs RFdiffusion's contig map to say which residues correspond, and
    guessing (truncate to the shorter) would quietly produce a number that
    means nothing. Skipped until we need it.
    """
    if reference is None:
        return None
    import biotite.structure as struc

    pred = _ca_atoms(prediction.path)
    ref = _ca_atoms(reference.path)
    if pred.array_length() != ref.array_length():
        return None
    # biotite >= 2.0 returns the transformation rather than the fitted copy.
    fitted = struc.superimpose(ref, pred).apply(pred)
    return float(struc.rmsd(ref, fitted))


def mean_plddt(prediction: Design, reference: Design | None) -> float | None:
    """Mean pLDDT. Prefers the model's own number, falls back to the B-factor
    column, which is where ColabFold writes per-residue pLDDT."""
    if "plddt" in prediction.scores:
        return float(prediction.scores["plddt"])
    atoms = _ca_atoms(prediction.path)
    if not hasattr(atoms, "b_factor"):
        return None
    return float(np.mean(atoms.b_factor))


def pae_interaction(prediction: Design, reference: Design | None) -> float | None:
    """Mean PAE across the interface, averaged over both off-diagonal blocks.

    Needs the PAE matrix, which lives in the predictor's scores JSON rather
    than in the PDB; the adapter records its path in meta.
    """
    pae = _load_pae(prediction)
    if pae is None:
        return None
    chains = _ca_atoms(prediction.path).chain_id
    if len(set(chains)) < 2 or len(chains) != pae.shape[0]:
        return None
    first = chains[0]
    mask = chains == first
    # Two chains only: everything that is not the first chain is the partner.
    # A real multi-chain complex would need a chain pair argument; skipped.
    inter = np.concatenate([pae[mask][:, ~mask].ravel(), pae[~mask][:, mask].ravel()])
    return float(inter.mean())


def _load_pae(prediction: Design) -> np.ndarray | None:
    candidate = prediction.meta.get("scores_json")
    path = Path(candidate) if candidate else None
    if path is None or not path.is_file():
        siblings = sorted(prediction.path.parent.glob("*scores*.json"))
        path = siblings[0] if siblings else None
    if path is None or not path.is_file():
        return None
    try:
        payload = json.loads(path.read_text())
    except (OSError, json.JSONDecodeError):
        return None
    pae = payload.get("pae") or payload.get("predicted_aligned_error")
    return None if pae is None else np.asarray(pae, dtype=float)


METRICS: dict[str, Metric] = {
    "sc_rmsd": sc_rmsd,
    "mean_plddt": mean_plddt,
    "pae_interaction": pae_interaction,
}


def score_designs(
    predictions: Sequence[Design],
    reference: Sequence[Design],
    metrics: Mapping[str, Metric] | None = None,
) -> list[Design]:
    """Return the predictions with `scores` populated.

    Each prediction is paired with its reference backbone through the ancestry
    the steps record: `parent_id` points at the sequence it was folded from, so
    reaching the backbone two steps up uses the ancestor ids in meta.
    """
    index = {d.id: d for d in reference}
    chosen = METRICS if metrics is None else metrics
    scored: list[Design] = []
    for prediction in predictions:
        ref = _reference_for(prediction, index)
        values = {
            name: value
            for name, metric in chosen.items()
            if (value := metric(prediction, ref)) is not None
        }
        scored.append(
            replace(
                prediction,
                scores=Scores({**prediction.scores, **values}),
                meta={**prediction.meta, "reference_id": ref.id if ref else None},
            )
        )
    return scored


def _reference_for(prediction: Design, index: Mapping[str, Design]) -> Design | None:
    ancestors = prediction.meta.get("ancestors") or ()
    for candidate in (prediction.parent_id, *reversed(tuple(ancestors))):
        if candidate is not None and candidate in index:
            return index[candidate]
    return None
