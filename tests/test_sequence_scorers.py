"""The sequence-only scorers: masked marginals, and the structure-free tag.

The three real backends need torch and model weights, so what is tested here is
the part that is this package's own: the masked-marginal arithmetic, one
forward pass per position rather than per variant, the refusal to score a
double mutant, and the tag that lets a caller pick the columns that survive a
wrong structure. The backends themselves are tested for the error they raise
when their dependency is absent.
"""

from __future__ import annotations

import contextlib
import math
from dataclasses import dataclass, field
from pathlib import Path
from types import SimpleNamespace
from typing import Mapping, Sequence

import numpy as np
import pytest

from design_pipeline import (
    REGISTRY,
    MissingDependencyError,
    enumerate_variants,
    get_structure,
    score_variants,
    select_pairs,
)
from design_pipeline.adapters import _plm
from design_pipeline.adapters._pdb import FIXTURE_CHAINS
from design_pipeline.adapters._plm import MaskedMarginalScorer


@dataclass
class StubLM:
    """A masked LM with a known answer, and a record of what it was asked.

    `log_probs` returns a distribution that depends only on the residue, so
    every expected score is arithmetic the test can do itself.
    """

    calls: list[tuple[str, tuple[int, ...]]] = field(default_factory=list)
    seen: list[Mapping[str, str]] = field(default_factory=list)

    def log_probs(
        self,
        sequences: Mapping[str, str],
        chain: str,
        positions: Sequence[int],
    ) -> Mapping[int, Mapping[str, float]]:
        self.calls.append((chain, tuple(positions)))
        self.seen.append(dict(sequences))
        # W is the model's favourite residue, A its least; everything else sits
        # in between, so log P(mut) - log P(wt) is predictable.
        table = {residue: -3.0 for residue in "ACDEFGHIKLMNPQRSTVWY"}
        table["W"] = -1.0
        table["A"] = -5.0
        return {position: dict(table) for position in positions}


@pytest.fixture
def complex_(tmp_path: Path):
    return get_structure(tmp_path / "runs", model_name="dummy_complex")


@pytest.fixture
def stub():
    backend = StubLM()
    REGISTRY["stub_plm"] = MaskedMarginalScorer(name="stub_plm", backend=backend)
    try:
        yield backend
    finally:
        del REGISTRY["stub_plm"]


def test_the_score_is_the_log_likelihood_ratio_against_the_wild_type(
    complex_, stub
) -> None:
    variants = enumerate_variants(complex_, chain="H", positions=[3], mutations="AW")
    scored = score_variants(complex_, variants, models=["stub_plm"])
    frame = scored.frame.set_index("variant_id")
    wt = FIXTURE_CHAINS["H"][2]

    # wt is Q, which the stub puts at -3.0.
    assert frame.loc[f"H:{wt}3W", "score_stub_plm"] == pytest.approx(-1.0 - -3.0)
    assert frame.loc[f"H:{wt}3A", "score_stub_plm"] == pytest.approx(-5.0 - -3.0)
    # Already higher-is-better, so raw and normalized agree.
    assert frame.loc[f"H:{wt}3W", "score_stub_plm_raw"] == pytest.approx(2.0)


def test_one_forward_pass_per_position_not_per_variant(complex_, stub) -> None:
    """Masked marginals, not pseudo-perplexity: the 19 substitutions at a
    position come out of the one pass that masked it."""
    variants = enumerate_variants(complex_, chain="H", positions=[1, 5, 9])
    score_variants(complex_, variants, models=["stub_plm"])

    assert len(stub.calls) == 1  # one call, carrying every position
    chain, positions = stub.calls[0]
    assert chain == "H"
    assert positions == (1, 5, 9)
    assert variants.candidates == 3 * 19


def test_a_double_mutant_gets_no_score_rather_than_a_sum(complex_, stub) -> None:
    """Two positional marginals multiplied is the additivity shortcut by
    another name, and the honest alternative is the pass-per-variant this class
    exists to avoid. An empty cell is the answer."""
    singles = enumerate_variants(complex_, chain="H", mutations="W")
    pairs = select_pairs(singles, max_pairs=3, min_distance=10.0)
    scored = score_variants(complex_, pairs, models=["stub_plm"])

    assert scored.frame["score_stub_plm"].isna().all()


def test_the_light_chain_comes_along_for_the_paired_models(complex_, stub) -> None:
    """A paired antibody model conditions the heavy chain on the light one, so
    a scan of H still has to reach the model with L attached."""
    variants = enumerate_variants(complex_, chain="H", positions=[1], mutations="W")
    score_variants(complex_, variants, models=["stub_plm"])

    assert stub.seen[0] == {"H": FIXTURE_CHAINS["H"], "L": FIXTURE_CHAINS["L"]}
    # The antigen is not an input to a sequence-only scorer.
    assert "A" not in stub.seen[0]


def test_antibody_chains_can_be_named_explicitly(complex_, stub) -> None:
    variants = enumerate_variants(complex_, chain="H", positions=[1], mutations="W")
    score_variants(
        complex_, variants, models=["stub_plm"], antibody_chains=["H"]
    )

    assert stub.seen[0] == {"H": FIXTURE_CHAINS["H"]}


def test_the_structure_is_genuinely_optional(tmp_path) -> None:
    """The fan-out hands every scorer the structure because that is how the
    candidates were built. This one reads only sequences off it, and
    `sequences=` skips the file entirely for identical scores."""
    from design_pipeline.adapters._tabular import write_candidates
    from design_pipeline.designs import Design

    candidates = write_candidates(
        tmp_path / "candidates.csv",
        [{"variant_id": "H:Q3W", "chain": "H", "position": 3, "wt": "Q", "mut": "W"}],
    )
    scorer = MaskedMarginalScorer(name="stub", backend=StubLM())
    produced = scorer.run(
        [], tmp_path / "out", candidates=candidates, sequences={"H": "ACQDEF"}
    )

    assert "H:Q3W,2.0" in produced[0].path.read_text()


def test_without_a_structure_or_sequences_it_says_so(tmp_path) -> None:
    from design_pipeline.adapters._tabular import write_candidates

    candidates = write_candidates(
        tmp_path / "candidates.csv",
        [{"variant_id": "H:Q3W", "chain": "H", "position": 3, "wt": "Q", "mut": "W"}],
    )
    scorer = MaskedMarginalScorer(name="stub", backend=StubLM())

    with pytest.raises(ValueError, match="sequences"):
        scorer.run([], tmp_path / "out", candidates=candidates)


def test_a_position_the_backend_declined_is_no_score(complex_) -> None:
    @dataclass
    class Silent:
        def log_probs(self, sequences, chain, positions):
            return {}

    REGISTRY["silent_plm"] = MaskedMarginalScorer(name="silent_plm", backend=Silent())
    try:
        variants = enumerate_variants(complex_, chain="H", positions=[1], mutations="W")
        scored = score_variants(complex_, variants, models=["silent_plm"])
    finally:
        del REGISTRY["silent_plm"]

    assert scored.frame["score_silent_plm"].isna().all()


@pytest.mark.parametrize(
    "name", ["ablang2", "esm1v", "esm2_650m", "antiberty"]
)
def test_every_sequence_scorer_is_registered_and_structure_free(name: str) -> None:
    adapter = REGISTRY[name]

    assert adapter.requires_structure is False
    # The matrix slot: one column, no cap, no change to the Protocol.
    assert adapter.per_variant is False
    assert name in _plm.SEQUENCE_SCORERS


def test_the_structure_free_subset_is_one_line() -> None:
    """What the tag is for: when the predicted complex may be wrong, these are
    the only columns still meaningful."""
    free = [
        name
        for name, adapter in REGISTRY.items()
        if not getattr(adapter, "requires_structure", True)
    ]

    assert set(_plm.SEQUENCE_SCORERS) <= set(free)
    # consensus is a positional prior and liabilities is a regex: neither opens
    # the coordinates either.
    assert "consensus" in free and "liabilities" in free
    assert "antifold" not in free and "gearbind" not in free


@pytest.mark.parametrize(
    "backend",
    [
        _plm.HuggingFaceMaskedLM("facebook/esm2_t33_650M_UR50D"),
        _plm.AbLang2MaskedLM(),
        _plm.AntiBertyMaskedLM(),
    ],
)
def test_a_backend_without_its_dependency_names_the_extra(backend) -> None:
    for module in ("torch", "ablang2", "antiberty"):
        pytest.importorskip  # noqa: B018
    try:
        import torch  # noqa: F401

        installed = True
    except ModuleNotFoundError:
        installed = False
    if installed:  # pragma: no cover - only where the extra is installed
        pytest.skip("the plm extra is installed, so there is no error to check")

    with pytest.raises(MissingDependencyError, match="plm"):
        backend.log_probs({"H": "ACDEF"}, "H", [1])


def test_logits_are_normalized_before_they_are_compared() -> None:
    """AbLang2 and AntiBERTy return logits, not log probabilities. The ratio
    against the wild type is unchanged by the normalizer, but a column of
    unnormalized logits is not comparable between positions - and `report`
    correlates the column as a whole."""
    total = _plm._logsumexp([0.0, 0.0])

    assert total == pytest.approx(math.log(2.0))
    assert _plm._logsumexp([]) == 0.0
    # Large values must not overflow.
    assert math.isfinite(_plm._logsumexp([1000.0, 1000.0]))


# ESM's vocabulary, abridged. <cls> has id 0, which Python treats as falsy.
ESM_IDS = {"<cls>": 0, "<pad>": 1, "<eos>": 2, "<mask>": 32} | {
    residue: i for i, residue in enumerate("LAGVSERTIDPKQNFYMHWC", start=4)
}


class _Tensor:
    """The few torch.Tensor operations HuggingFaceMaskedLM.log_probs uses, on numpy."""

    def __init__(self, array) -> None:
        self.array = np.asarray(array)

    @property
    def shape(self):
        return self.array.shape

    def __getitem__(self, key):
        value = self.array[key]
        return _Tensor(value) if isinstance(value, np.ndarray) else value

    def __setitem__(self, key, value) -> None:
        self.array[key] = value

    def repeat(self, rows: int, columns: int) -> "_Tensor":
        return _Tensor(np.tile(self.array, (rows, columns)))

    def to(self, _device) -> "_Tensor":
        return self

    def float(self) -> "_Tensor":
        return self


def _fake_torch() -> SimpleNamespace:
    def log_softmax(x: _Tensor, dim: int = -1) -> _Tensor:
        shifted = x.array - x.array.max(axis=dim, keepdims=True)
        return _Tensor(shifted - np.log(np.exp(shifted).sum(axis=dim, keepdims=True)))

    return SimpleNamespace(no_grad=contextlib.nullcontext, log_softmax=log_softmax)


def _fake_esm(sequence: str, *, with_cls: bool):
    """A tokenizer and a perfect masked LM over `sequence`.

    The model answers only at a masked index, with the residue that was really
    there; everywhere else it is uniform. So a mask that lands on the wrong
    index makes the model name the wrong residue.
    """
    original = [ESM_IDS[r] for r in sequence]
    if with_cls:
        original = [ESM_IDS["<cls>"], *original, ESM_IDS["<eos>"]]

    class Tokenizer:
        cls_token_id = ESM_IDS["<cls>"] if with_cls else None
        mask_token_id = ESM_IDS["<mask>"]

        def __call__(self, _sequence, return_tensors=None):
            return {"input_ids": _Tensor([original])}

        def convert_tokens_to_ids(self, token):
            return ESM_IDS[token]

    class Model:
        device = "cpu"

        def __call__(self, ids: _Tensor):
            logits = np.zeros((*ids.shape, max(ESM_IDS.values()) + 1))
            for row, index in zip(*np.nonzero(ids.array == ESM_IDS["<mask>"])):
                logits[row, index, original[index]] = 10.0
            return SimpleNamespace(logits=_Tensor(logits))

    return _fake_torch(), Tokenizer(), Model()


@pytest.mark.parametrize("with_cls", [True, False])
def test_the_mask_lands_on_the_residue_it_scores(monkeypatch, with_cls: bool) -> None:
    """ESM's <cls> token has id 0. `cls_token_id or -1` read that as "no cls
    token", so the offset was 0 and every mask landed one residue early: each
    row scored its neighbour, and the wild-type check could not see it because
    it compares against the sequence, not against what the model was shown."""
    sequence = "ACDEFGHIKLMNPQRSTVWY"  # no two neighbours alike
    loaded = _fake_esm(sequence, with_cls=with_cls)
    monkeypatch.setattr(_plm.HuggingFaceMaskedLM, "_load", lambda self: loaded)
    backend = _plm.HuggingFaceMaskedLM("fake/esm")

    found = backend.log_probs({"H": sequence}, "H", list(range(1, len(sequence) + 1)))

    for position, residue in enumerate(sequence, start=1):
        row = found[position]
        scored = max(row, key=row.get)
        assert scored == residue, f"position {position} scored {scored}, not {residue}"
