"""Adapter instances.

Importing this imports no model: the container adapter is plain data, and the
pip-backed adapters (MPNN, ColabFold) check for their dependency inside run().
Module names are underscored because the public names here are the adapter
*instances*.
"""

from __future__ import annotations

from ._boltz import BoltzAdapter
from ._chai import ChaiAdapter
from ._colabfold import ColabFoldAdapter
from ._dummy import DummyBackbones, DummyFold, DummySequences, DummyStructures
from ._mpnn import ligand_mpnn, protein_mpnn, soluble_mpnn
from ._protenix import ProtenixAdapter
from ._rfdiffusion import rfdiffusion

colabfold = ColabFoldAdapter()
chai = ChaiAdapter()
boltz = BoltzAdapter()
protenix = ProtenixAdapter()
dummy_backbones = DummyBackbones()
dummy_sequences = DummySequences()
dummy_structures = DummyStructures()
# One fixture class, one instance per co-folding model, so an offline pool
# mixes producers the way a real model list does.
dummy_chai = DummyFold(producer="dummy_chai")
dummy_boltz = DummyFold(producer="dummy_boltz")
dummy_protenix = DummyFold(producer="dummy_protenix")

__all__ = [
    "rfdiffusion",
    "protein_mpnn",
    "ligand_mpnn",
    "soluble_mpnn",
    "chai",
    "boltz",
    "protenix",
    "colabfold",
    "dummy_backbones",
    "dummy_sequences",
    "dummy_structures",
    "dummy_chai",
    "dummy_boltz",
    "dummy_protenix",
]
