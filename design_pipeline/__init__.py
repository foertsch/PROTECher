"""design_pipeline: antibody affinity maturation as five functions over files.

    from design_pipeline import (
        get_structure, enumerate_variants, score_variants, gate, select_pairs,
    )

    complex_ = get_structure("runs/her2_v1", pdb_id="1N8Z")
    variants = enumerate_variants(complex_, chain="H")          # 4,598 rows
    scores = score_variants(complex_, variants, models=["antifold", "esm_if"])
    survivors = gate(scores, checks=["liabilities"])
    pairs = select_pairs(survivors, max_pairs=20, min_distance=10.0)
    pair_scores = score_variants(complex_, pairs, models=["antifold", "gearbind"])

That is the whole surface. `score_variants` is called twice, on singles and on
pairs, because a pair is an ordinary candidate - and nothing in this package
estimates a pair's score by adding its singles' scores together.

Everything is files on disk: the two artifact types, `Complex` and
`VariantTable`, hold paths plus the config that produced them, and nothing is
held in memory between steps. The round loop, if you want one, is a `for`
statement in your own script.

Importing this imports no model and no heavy dependency.
"""

from __future__ import annotations

from .api import gate, get_structure, score_variants
from .combine import frontier, pick_designs, propose_sets
from .denovo import diffuse_loops
from .errors import (
    MissingDependencyError,
    MissingGPUError,
    MissingLigandError,
    MissingRuntimeError,
    PipelineError,
    RunFailedError,
    UnknownModelError,
)
from .designs import Design
from .registry import REGISTRY, ModelAdapter, get_adapter
from .structures import Complex
from .variants import SCHEMA, VariantTable, enumerate_variants, select_pairs

__all__ = [
    # the five steps
    "get_structure",
    "enumerate_variants",
    "score_variants",
    "gate",
    "select_pairs",
    # past pairs: sets of five or more, and the designs chosen from the scored sets
    "propose_sets",
    "pick_designs",
    "frontier",
    # the optional de novo branch
    "diffuse_loops",
    # the two artifact types
    "Complex",
    "VariantTable",
    "SCHEMA",
    # swapping models out: an adapter speaks Design
    "REGISTRY",
    "Design",
    "ModelAdapter",
    "get_adapter",
    # errors
    "PipelineError",
    "UnknownModelError",
    "MissingRuntimeError",
    "MissingDependencyError",
    "MissingGPUError",
    "MissingLigandError",
    "RunFailedError",
]
