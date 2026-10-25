"""design_pipeline: iterative protein design as five functions over `list[Design]`.

    from design_pipeline import (
        make_protein_backbones, design_sequences,
        predict_structure, score_designs, select,
    )

The round loop is a `for` statement in your own script, not a library feature.
See examples/refine.py.

Importing this imports no model and no heavy dependency.
"""

from __future__ import annotations

from .api import design_sequences, make_protein_backbones, predict_structure
from .designs import Design, Scores, lineage, select, to_dataframe
from .errors import (
    MissingDependencyError,
    MissingLigandError,
    MissingRuntimeError,
    PipelineError,
    RunFailedError,
    UnknownModelError,
)
from .registry import REGISTRY, ModelAdapter, get_adapter
from .scoring import METRICS, score_designs

__all__ = [
    # the five steps
    "make_protein_backbones",
    "design_sequences",
    "predict_structure",
    "score_designs",
    "select",
    # the record
    "Design",
    "Scores",
    "lineage",
    "to_dataframe",
    # swapping models out
    "REGISTRY",
    "ModelAdapter",
    "get_adapter",
    "METRICS",
    # errors
    "PipelineError",
    "UnknownModelError",
    "MissingRuntimeError",
    "MissingDependencyError",
    "MissingLigandError",
    "RunFailedError",
]
