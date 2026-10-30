"""The adapter protocol and the registry.

The registry is a plain dict. Adding a model is one line here, or one
assignment from outside the package — no entry points, no discovery.
"""

from __future__ import annotations

from pathlib import Path
from typing import Any, Protocol, Sequence, runtime_checkable

from . import adapters
from .designs import Design
from .errors import UnknownModelError


@runtime_checkable
class ModelAdapter(Protocol):
    def run(
        self,
        inputs: Sequence[Design],
        output_dir: Path,
        **kwargs: Any,
    ) -> list[Design]:
        """Write outputs into output_dir and return one Design per output.

        The adapter sets `parent_id` and parses its own output format into
        `scores`. Synchronous on purpose: an executor can be slid underneath
        this signature later.
        """
        ...


REGISTRY: dict[str, ModelAdapter] = {
    # Step 1, structure. The download is the default: experimental coordinates
    # beat a prediction, and the predictors here are the fallback for a target
    # with no solved complex.
    "pdb_download": adapters.pdb_download,
    "alphafold3": adapters.alphafold3,
    # boltz2 and chai1 are the open-weight alternatives. Registered under both
    # names because the in-process adapters already here are the same models;
    # the container builds the spec also mentions would add a second way to run
    # one model, so they are skipped.
    "boltz2": adapters.boltz,
    "chai1": adapters.chai,
    "dummy_complex": adapters.dummy_complex,
    # An ensemble with confidence, for the validation notebook: honours
    # chain_ids and num_seeds, puts most seeds on one pose and the rest on
    # another, and reports ipTM and a PAE matrix.
    "dummy_predictor": adapters.dummy_predictor,
    # Step 3, scoring. Matrix scorers first - one forward pass over the
    # complex, then a lookup per single mutation - then the per-variant ones,
    # which build and relax a structure per candidate.
    "antifold": adapters.antifold,
    "esm_if": adapters.esm_if,
    "abmpnn": adapters.abmpnn,
    "proteinmpnn_score": adapters.proteinmpnn_score,
    "consensus": adapters.consensus,
    # Sequence-only scorers: the same position-by-residue matrix, read off the
    # sequence rather than the structure, so they are the only columns still
    # meaningful when the predicted complex is wrong. requires_structure=False
    # is the tag that selects them.
    "ablang2": adapters.ablang2,
    "esm1v": adapters.esm1v,
    "esm2_650m": adapters.esm2_650m,
    "antiberty": adapters.antiberty,
    "gearbind": adapters.gearbind,
    "rosetta_ddg": adapters.rosetta_ddg,
    "foldx": adapters.foldx,
    "dummy_scorer": adapters.dummy_scorer,
    "dummy_per_variant": adapters.dummy_per_variant,
    "dummy_sequence_scorer": adapters.dummy_sequence_scorer,
    # Step 4, gate.
    "liabilities": adapters.liabilities,
    "thermompnn": adapters.thermompnn,
    "aggrescan3d": adapters.aggrescan3d,
    "netsolp": adapters.netsolp,
    "dummy_check": adapters.dummy_check,
    # The optional de novo branch, reachable through diffuse_loops only.
    "rfdiffusion": adapters.rfdiffusion,
    # One MPNN codebase, three models, selected by --model_type.
    "protein_mpnn": adapters.protein_mpnn,
    "ligand_mpnn": adapters.ligand_mpnn,
    "soluble_mpnn": adapters.soluble_mpnn,
    # Co-folding models. Chai-1 is the default; ColabFold stays, with the
    # non-commercial AlphaFold2 weights that are the reason it is not.
    "chai": adapters.chai,
    "boltz": adapters.boltz,
    "protenix": adapters.protenix,
    "colabfold": adapters.colabfold,
    "dummy_backbones": adapters.dummy_backbones,
    "dummy_sequences": adapters.dummy_sequences,
    "dummy_structures": adapters.dummy_structures,
    "dummy_chai": adapters.dummy_chai,
    "dummy_boltz": adapters.dummy_boltz,
    "dummy_protenix": adapters.dummy_protenix,
}


def get_adapter(model_name: str) -> ModelAdapter:
    try:
        return REGISTRY[model_name]
    except KeyError:
        known = ", ".join(sorted(REGISTRY))
        raise UnknownModelError(
            f"No adapter registered for {model_name!r}. Known models: {known}. "
            f"Register your own with design_pipeline.REGISTRY['{model_name}'] = adapter."
        ) from None
