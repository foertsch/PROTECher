"""Adapter instances.

Importing this imports no model: the container adapter is plain data, and the
pip-backed adapters (MPNN, ColabFold) check for their dependency inside run().
Module names are underscored because the public names here are the adapter
*instances*.
"""

from __future__ import annotations

from ._af3 import alphafold3
from ._boltz import BoltzAdapter
from ._chai import ChaiAdapter
from ._colabfold import ColabFoldAdapter
from ._dummy import DummyBackbones, DummyFold, DummySequences, DummyStructures
from ._gates import (
    aggrescan3d,
    dummy_check,
    liabilities,
    netsolp,
    thermompnn,
)
from ._mpnn import ligand_mpnn, protein_mpnn, soluble_mpnn
from ._pdb import dummy_complex, dummy_predictor, pdb_download
from ._plm import ablang2, antiberty, esm1v, esm2_650m
from ._protenix import ProtenixAdapter
from ._rfdiffusion import rfdiffusion
from ._scorers import (
    abmpnn,
    antifold,
    consensus,
    dummy_per_variant,
    dummy_scorer,
    dummy_sequence_scorer,
    esm_if,
    foldx,
    gearbind,
    proteinmpnn_score,
    rosetta_ddg,
)

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
    # step 1, structure
    "pdb_download",
    "alphafold3",
    "chai",
    "boltz",
    "protenix",
    "colabfold",
    "dummy_complex",
    "dummy_predictor",
    # step 3, scoring
    "antifold",
    "esm_if",
    "abmpnn",
    "proteinmpnn_score",
    "consensus",
    # sequence-only: the same matrix, from the sequence alone
    "ablang2",
    "esm1v",
    "esm2_650m",
    "antiberty",
    "gearbind",
    "rosetta_ddg",
    "foldx",
    "dummy_scorer",
    "dummy_per_variant",
    "dummy_sequence_scorer",
    # step 4, gate
    "liabilities",
    "thermompnn",
    "aggrescan3d",
    "netsolp",
    "dummy_check",
    # the optional de novo branch
    "rfdiffusion",
    "protein_mpnn",
    "ligand_mpnn",
    "soluble_mpnn",
    "dummy_backbones",
    "dummy_sequences",
    "dummy_structures",
    "dummy_chai",
    "dummy_boltz",
    "dummy_protenix",
]
