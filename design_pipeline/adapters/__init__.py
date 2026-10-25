"""Adapter instances.

Importing this imports no model: the container adapter is plain data, and the
pip-backed adapters (MPNN, ColabFold) check for their dependency inside run().
Module names are underscored because the public names here are the adapter
*instances*.
"""

from __future__ import annotations

from ._colabfold import ColabFoldAdapter
from ._dummy import DummyBackbones, DummySequences, DummyStructures
from ._mpnn import ligand_mpnn, protein_mpnn, soluble_mpnn
from ._rfdiffusion import rfdiffusion

colabfold = ColabFoldAdapter()
dummy_backbones = DummyBackbones()
dummy_sequences = DummySequences()
dummy_structures = DummyStructures()

__all__ = [
    "rfdiffusion",
    "protein_mpnn",
    "ligand_mpnn",
    "soluble_mpnn",
    "colabfold",
    "dummy_backbones",
    "dummy_sequences",
    "dummy_structures",
]
