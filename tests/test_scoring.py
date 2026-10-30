"""Scoring: real numbers off real (tiny) files, numpy + biotite only."""

from __future__ import annotations

from pathlib import Path

import pytest

from design_pipeline import Design
from design_pipeline.denovo import (
    design_sequences,
    make_protein_backbones,
    predict_structure,
)
from design_pipeline.scoring import score_designs
from design_pipeline.adapters._dummy import CHAIN_A, CHAIN_B, write_pdb
from design_pipeline.scoring import pae_interaction, sc_rmsd


@pytest.fixture
def scored(tmp_path: Path):
    target = tmp_path / "target.pdb"
    target.write_text("ATOM      1  CA  GLY A   1       0.000   0.000   0.000\n")
    backbones = make_protein_backbones(
        target, tmp_path / "bb", model_name="dummy_backbones", num_designs=2
    )
    seqs = design_sequences(
        backbones, tmp_path / "seqs", model_name="dummy_sequences", num_seqs=2
    )
    preds = predict_structure(seqs, tmp_path / "fold", model_name="dummy_structures")
    return backbones, score_designs(preds, reference=backbones)


def test_all_three_metrics_land(scored) -> None:
    _, designs = scored
    for design in designs:
        assert design.scores["sc_rmsd"] > 0
        assert 50 <= design.scores["mean_plddt"] <= 90
        assert design.scores["pae_interaction"] > 0


def test_predictions_are_paired_to_their_own_backbone(scored) -> None:
    backbones, designs = scored
    ids = {d.id for d in backbones}
    assert {d.meta["reference_id"] for d in designs} <= ids
    # Two backbones, both actually used as references.
    assert len({d.meta["reference_id"] for d in designs}) == 2


def test_scoring_does_not_touch_the_predictions(scored) -> None:
    _, designs = scored
    # score_designs returns new records; the originals keep their own scores.
    assert all("sc_rmsd" in d.scores for d in designs)
    assert all(d.producer == "dummy_structures" for d in designs)


def test_sc_rmsd_is_zero_against_itself(tmp_path) -> None:
    pdb = tmp_path / "a.pdb"
    write_pdb(pdb, token="x")
    design = Design(id="a", path=pdb, parent_id=None, producer="t")
    assert sc_rmsd(design, design) == pytest.approx(0.0, abs=1e-6)


def test_sc_rmsd_skips_a_length_mismatch(tmp_path) -> None:
    """Mapping unequal lengths needs the contig map; we return nothing instead
    of a number that means nothing."""
    full = tmp_path / "full.pdb"
    write_pdb(full, token="x")
    short = tmp_path / "short.pdb"
    write_pdb(short, token="x", residues=(3, 3))
    pred = Design(id="p", path=full, parent_id=None, producer="t")
    ref = Design(id="r", path=short, parent_id=None, producer="t")
    assert sc_rmsd(pred, ref) is None


def test_pae_interaction_uses_only_the_interface(tmp_path) -> None:
    """The fixture writes intra-chain PAE around 3-7 and inter-chain 12-16, so
    an interface mean must land in the upper band."""
    import json

    pdb = tmp_path / "p.pdb"
    write_pdb(pdb, token="x")
    scores_json = tmp_path / "p_scores.json"
    n = CHAIN_A + CHAIN_B
    pae = [
        [1.0 if (i < CHAIN_A) == (j < CHAIN_A) else 20.0 for j in range(n)]
        for i in range(n)
    ]
    scores_json.write_text(json.dumps({"pae": pae}))
    design = Design(
        id="p",
        path=pdb,
        parent_id=None,
        producer="t",
        meta={"scores_json": str(scores_json)},
    )
    assert pae_interaction(design, None) == pytest.approx(20.0)


def test_pae_interaction_is_skipped_without_a_matrix(tmp_path) -> None:
    pdb = tmp_path / "p.pdb"
    write_pdb(pdb, token="x")
    design = Design(id="p", path=pdb, parent_id=None, producer="t")
    assert pae_interaction(design, None) is None


def test_metrics_are_extensible(scored) -> None:
    _, designs = scored
    extra = score_designs(
        designs, reference=[], metrics={"length": lambda pred, ref: 42.0}
    )
    assert extra[0].scores["length"] == 42.0
