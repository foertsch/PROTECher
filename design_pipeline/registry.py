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
