"""Sequence-only scorers: a protein language model, scored by masked marginals.

The third cost class. A matrix scorer reads a probability matrix off one pass
over the *structure*; these read the same shape off passes over the *sequence*
alone - no coordinates, no antigen. So they occupy the matrix slot exactly: one
column per model, `per_variant=False`, no change to the adapter Protocol and no
change to the public API. What distinguishes them is the `requires_structure`
tag, which is False here and True for everything that touches coordinates.

They exist primarily as insurance. When no solved complex exists and the
predicted one has the paratope on the wrong patch, every structure-based column
in the table is scoring a fiction, and these are the only columns still
meaningful. Weight them below the structure-based scorers when the structure is
trustworthy - `design_pipeline.confidence.recommend_weights` does exactly that,
from the seed agreement and the decoy control.

Scoring is masked marginals: mask position i, one forward pass, take
log P(mut) - log P(wt) from it. Not pseudo-perplexity - that needs a pass per
*variant* rather than per *position*, which is 19x the compute for no added
signal at single-mutation resolution.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from pathlib import Path
from typing import Any, Mapping, Protocol, Sequence

from ..containers import make_id
from ..matrices import check_vocabulary, from_position_mapping, score_candidates
from ..designs import Design
from ..errors import MissingDependencyError
from ..structures import AMINO_ACIDS, chain_sequence
from ._tabular import CANDIDATES, SCORES, read_candidates, write_scores

PairedChains = ("H", "L")
"""Chain ids treated as the antibody when nothing says otherwise. A heuristic
on naming, and the reason `antibody_chains=` is a kwarg: the paired models want
the whole Fv even when only the heavy chain is being mutated."""


class MaskedLM(Protocol):
    """How one adapter talks to three model libraries.

    Not a scorer Protocol - the scorer above it satisfies the one existing
    `ModelAdapter` - but the same role `Collect` plays for `ContainerAdapter`:
    the part that differs between ESM, AbLang2 and AntiBERTy is data and a
    tokenizer, not a subclass of the pipeline's interface.

    `sequences` is every antibody chain, keyed by chain id, because the paired
    antibody models condition the heavy chain on the light one. `chain` is the
    one being mutated and `positions` are 1-based into it.
    """

    def log_probs(
        self,
        sequences: Mapping[str, str],
        chain: str,
        positions: Sequence[int],
    ) -> Mapping[int, Mapping[str, float]]:
        """{position: {residue: log probability}} with that position masked."""
        ...


@dataclass(frozen=True)
class MaskedMarginalScorer:
    """Satisfies the ModelAdapter Protocol. One class, one instance per model."""

    name: str
    backend: MaskedLM
    per_variant: bool = False
    requires_structure: bool = False

    def run(
        self,
        inputs: Sequence[Design],
        output_dir: Path,
        *,
        candidates: Path | str | None = None,
        sequences: Mapping[str, str] | None = None,
        antibody_chains: Sequence[str] | None = None,
        extra_args: Sequence[str] | None = None,
        **kwargs: Any,
    ) -> list[Design]:
        """`sequences` makes the structure genuinely optional.

        The fan-out hands every scorer the structure, because that is how the
        candidate set was built; this one reads only the chain sequences off it
        and never the coordinates. Pass `sequences={"H": ..., "L": ...}` to skip
        the file entirely - the scores are identical either way, which is the
        whole point of the class.
        """
        if candidates is None:
            raise ValueError(f"{self.name} needs candidates=<path to {CANDIDATES}>")
        wanted = read_candidates(Path(candidates))
        chains = _chains(wanted, sequences, inputs, antibody_chains)

        # One pass per position, shared by the 19 substitutions at it. Collect
        # the positions first so a chain is visited once.
        needed: dict[str, list[int]] = {}
        for mutations in wanted.values():
            if len(mutations) != 1:
                continue
            chain, position, _wt, _mut = mutations[0]
            needed.setdefault(chain, []).append(position)

        # The backend answers {position: {residue: log probability}}; the
        # matrix is what carries the chain sequence alongside it, so a position
        # that does not hold the residue the candidate table claims raises
        # instead of scoring. Both rules a matrix scorer obeys - no additivity
        # for a multi-mutation candidate, no zero for a position the model did
        # not answer - live in `score_candidates`, shared with AntiFold.
        matrices = {
            chain: from_position_mapping(
                chain,
                self.backend.log_probs(chains, chain, sorted(set(positions))),
                chains[chain],
                # Every backend here checks its own vocabulary against the
                # alphabet first, so a missing residue at a position it did
                # answer for is a bug and not an abstention.
                require_full_rows=True,
            )
            for chain, positions in needed.items()
        }
        rows = score_candidates(matrices, wanted)

        path = write_scores(output_dir / SCORES, rows)
        return [
            Design(
                id=make_id(path, output_dir),
                path=path,
                parent_id=inputs[0].id if inputs else None,
                producer=self.name,
            )
        ]


def _chains(
    wanted: Mapping[str, Sequence[tuple[str, int, str, str]]],
    sequences: Mapping[str, str] | None,
    inputs: Sequence[Design],
    antibody_chains: Sequence[str] | None,
) -> dict[str, str]:
    """The antibody sequences to condition on, from kwargs or off the structure.

    Defaults to the chains being mutated plus any H/L the structure has, so a
    scan of the heavy chain still reaches a paired model with its light chain
    attached. Name them explicitly when the ids are not H and L.
    """
    mutated = {chain for mutations in wanted.values() for chain, _, _, _ in mutations}
    if sequences is not None:
        if antibody_chains:
            missing = [c for c in antibody_chains if c not in sequences]
            if missing:
                raise KeyError(f"sequences= has no chain(s) {', '.join(missing)}")
            return {c: sequences[c] for c in antibody_chains}
        return dict(sequences)

    if not inputs:
        raise ValueError(
            "a sequence scorer needs either the complex or sequences={chain: seq}"
        )
    structure = inputs[0].path
    from ..structures import chain_ids

    present = set(chain_ids(structure))
    if antibody_chains:
        chosen = list(antibody_chains)
    else:
        chosen = sorted(mutated | (set(PairedChains) & present))
        # The order of this dict is the order a paired model reads the chains
        # in, and sorting is not that order: on 1N8Z it puts the light chain
        # ('A') before the heavy one ('B'). Where the ids are H and L, sorting
        # happens to be heavy-first and PairedChains pins it anyway. Where they
        # are not, there is nothing here that knows which is which, so this
        # refuses rather than guessing - the bug AbLang2's `order` used to hide.
        if len(chosen) > 1 and not (set(PairedChains) & set(chosen)):
            raise ValueError(
                f"chains {', '.join(chosen)} are not named H and L, so nothing "
                f"here can tell which is the heavy chain, and a paired model "
                f"reads them in the order they arrive. Pass "
                f"antibody_chains=[heavy, light] explicitly."
            )
    return {chain: chain_sequence(structure, chain) for chain in chosen}


@dataclass(frozen=True)
class HuggingFaceMaskedLM:
    """ESM-1v and ESM-2, through `transformers`.

    One forward pass per masked position, batched. ESM-1v is the one of the two
    built for this task - it was trained and published as a zero-shot variant
    effect predictor, and masked marginals is the protocol its paper uses - so
    it is the default of the pair; ESM-2 is the general baseline, scored
    identically.

    These see the mutated chain alone. Conditioning a general protein model on
    a ':'-joined Fv is not something it was trained for, and the antibody
    models below do the paired case properly.
    """

    model_id: str
    batch_size: int = 8
    device: str | None = None
    """None picks cuda when torch reports it, else cpu. A full 242-position
    scan on cpu is minutes, not hours, which is tolerable for one column."""

    def log_probs(
        self,
        sequences: Mapping[str, str],
        chain: str,
        positions: Sequence[int],
    ) -> Mapping[int, Mapping[str, float]]:
        torch, tokenizer, model = self._load()
        sequence = sequences[chain]
        encoded = tokenizer(sequence, return_tensors="pt")
        ids = encoded["input_ids"]
        # ESM prepends a cls token; locating it beats hardcoding the offset,
        # because the tokenizers of these two models do not agree about it.
        # Compare to None, not `or -1`: ESM's cls id is 0, which `or` skips.
        cls_id = tokenizer.cls_token_id
        offset = 1 if cls_id is not None and int(ids[0, 0]) == int(cls_id) else 0
        # The column-axis guard. `convert_tokens_to_ids` answers with the
        # *unknown* id for a token it does not have, so a residue this model
        # lacks would come back as a real id pointing at the wrong
        # distribution; `reserved` is what turns that into an error.
        vocabulary = check_vocabulary(
            self.model_id,
            {
                residue: tokenizer.convert_tokens_to_ids(residue)
                for residue in AMINO_ACIDS
            },
            # getattr, because a tokenizer is only required to have the ids it
            # uses: a stub or a model without a pad token has neither.
            reserved={
                "unknown": getattr(tokenizer, "unk_token_id", None),
                "mask": getattr(tokenizer, "mask_token_id", None),
                "padding": getattr(tokenizer, "pad_token_id", None),
                "cls": cls_id,
            },
        )

        found: dict[int, Mapping[str, float]] = {}
        with torch.no_grad():
            for start in range(0, len(positions), self.batch_size):
                batch = list(positions[start : start + self.batch_size])
                stacked = ids.repeat(len(batch), 1)
                for row, position in enumerate(batch):
                    index = offset + position - 1
                    if index >= stacked.shape[1]:
                        raise IndexError(
                            f"position {position} is past the end of chain "
                            f"{chain!r} ({len(sequence)} residues) as tokenized"
                        )
                    stacked[row, index] = tokenizer.mask_token_id
                logits = model(stacked.to(model.device)).logits
                for row, position in enumerate(batch):
                    index = offset + position - 1
                    logp = torch.log_softmax(logits[row, index].float(), dim=-1)
                    found[position] = {
                        residue: float(logp[token])
                        for residue, token in vocabulary.items()
                    }
        return found

    def _load(self):
        try:
            import torch
            from transformers import AutoModelForMaskedLM, AutoTokenizer
        except ModuleNotFoundError as exc:
            raise MissingDependencyError(
                f"{self.model_id} needs torch and transformers. Install them "
                f"with `uv pip install 'design-pipeline[plm]'`. To run the "
                f"pipeline with no model at all, use model_name='dummy_scorer'."
            ) from exc
        tokenizer = AutoTokenizer.from_pretrained(self.model_id)
        model = AutoModelForMaskedLM.from_pretrained(self.model_id)
        device = self.device or ("cuda" if torch.cuda.is_available() else "cpu")
        return torch, tokenizer, model.to(device).eval()


@dataclass(frozen=True)
class AbLang2MaskedLM:
    """AbLang2, trained on paired OAS sequences.

    The best prior of the four for framework and germline positions, because it
    has seen what an antibody is allowed to look like there - which a general
    protein model has not, and which no structure-based scorer knows about at
    all.

    AbLang2 takes the pair as TWO strings - `model([[heavy, light]])` - and does
    its own formatting: `add_extra_tokens` builds `f"<{heavy}>|<{light}>"`, with
    the empty `<>` dropped when one chain is absent. Handing it the joined
    string instead raises `ValueError: not enough values to unpack (expected 2,
    got 1)` from `heavy, light = seq`, which is what this adapter used to do, and
    which crashed any campaign that reached `ablang2` - the default scorer.

    That formatting also fixes the row layout of `mode="likelihood"`, which
    returns per-position logits over the 26-token vocabulary. Measured on
    ablang2 0.2.1 with a 120-residue heavy and a 107-residue light chain: 232
    rows, i.e. `len(heavy) + len(light) + 5`, laid out as

        0            `<`
        1 .. nH      the heavy residues, so heavy position i is row i
        nH+1         `>`
        nH+2         `|`
        nH+3         `<`
        nH+4 ..      the light residues, so light position j is row nH+3+j
        last         `>`

    so the offset before a chain is one token for its own `<`, plus
    `len(chain) + 4` for every chain ahead of it.
    """

    device: str | None = None

    def log_probs(
        self,
        sequences: Mapping[str, str],
        chain: str,
        positions: Sequence[int],
    ) -> Mapping[int, Mapping[str, float]]:

        model = self._load()
        # The first chain of `order` is the one AbLang2 reads as the heavy one:
        # its tokenizer does `heavy, light = seq` and builds
        # f"<{heavy}>|<{light}>". So the fallback is the caller's own order -
        # `_chains` preserves `antibody_chains`, which the pipeline passes as
        # [heavy, light] - and NOT `sorted(sequences)`. Sorting put 1N8Z's
        # light chain ('A') ahead of its heavy chain ('B') and fed the model a
        # light chain in the heavy slot; `before` is computed from the same
        # order, so the rows stayed self-consistent and an offset sweep could
        # not see it. Measured on the 1N8Z variable domains: CDR wt-argmax
        # 14.9% with the chains swapped against 27.7% with them the right way
        # round.
        order = [c for c in PairedChains if c in sequences] or list(sequences)
        if chain not in order:
            raise KeyError(f"chain {chain!r} is not among {', '.join(order)}")
        # Where the mutated chain starts among the rows the model returns, so a
        # position in the chain becomes a row index. The 4 is AbLang2's own
        # punctuation around each chain (`>`, `|`, `<`) and the leading 1 is this
        # chain's own `<`; see the layout in the class docstring. This was
        # `len + 1` with no leading offset, which read the row before the one it
        # wanted for the heavy chain and three rows early for the light.
        before = 1 + sum(len(sequences[c]) + 4 for c in order[: order.index(chain)])

        # The column-axis guard, once rather than per position: AbLang2's own
        # residue-to-token mapping has to cover the twenty with distinct ids.
        vocabulary = check_vocabulary("ablang2", model.tokenizer.aa_to_token)

        found: dict[int, dict[str, float]] = {}
        for position in positions:
            masked = dict(sequences)
            residues = list(masked[chain])
            residues[position - 1] = "*"  # AbLang2's mask token
            masked[chain] = "".join(residues)
            # A list of two strings, never a joined one. AbLang2 drops the
            # empty `<>` itself, so an absent second chain is "".
            pair = [masked[c] for c in order]
            if len(pair) == 1:
                pair = [pair[0], ""]
            logits = model([pair], mode="likelihood")[0]
            row = logits[before + position - 1]
            total = _logsumexp([float(v) for v in row])
            found[position] = {
                residue: float(row[vocabulary[residue]]) - total
                for residue in AMINO_ACIDS
            }
        return found

    def _load(self):
        try:
            import ablang2
        except ModuleNotFoundError as exc:
            raise MissingDependencyError(
                "ablang2 is not installed. `uv pip install "
                "'design-pipeline[plm]'`, which also brings torch. Its weights "
                "download on first use."
            ) from exc
        return ablang2.pretrained(device=self.device or "cpu")


@dataclass(frozen=True)
class AntiBertyMaskedLM:
    """AntiBERTy, the optional antibody-specific alternative to AbLang2.

    Kept because two antibody priors trained on different data disagree in
    useful ways, and because AntiBERTy installs more easily than AbLang2 on
    some machines.

    Per-position distributions come from the masked LM behind the runner, not
    from `fill_masks`. In antiberty 1.0.0 `fill_masks(sequences)` returns the
    completed strings and accepts no `return_logits`, so the call this adapter
    used to make raised `TypeError: AntiBERTyRunner.fill_masks() got an
    unexpected keyword argument 'return_logits'` and killed any campaign that
    reached this scorer. Its other public methods are no help either:
    `pseudo_log_likelihood` is one scalar per sequence and `embed` returns
    hidden states.

    `AntiBERTyRunner` does expose `.model` (a BERT masked LM, head included)
    and `.tokenizer`, which is all this needs. Two details of that tokenizer,
    both measured on antiberty 1.0.0: `vocab` is a LIST, so residue ids come
    from the `token_to_id` mapping rather than from indexing `vocab` (what this
    adapter used to do, which would have raised even with logits in hand), and
    `[CLS]` occupies row 0, so chain position i is row i. At that offset the
    wild-type residue ranks first at 4 of 6 sampled positions of 1N8Z's heavy
    chain, against 5.67 and 13.67 mean rank for the neighbouring offsets.
    """

    device: str | None = None

    def log_probs(
        self,
        sequences: Mapping[str, str],
        chain: str,
        positions: Sequence[int],
    ) -> Mapping[int, Mapping[str, float]]:

        import torch

        runner = self._load()
        model, tokenizer = runner.model, runner.tokenizer
        # The column-axis guard. `token_to_id` is the model's own mapping;
        # `vocab` next to it is a LIST in antiberty 1.0.0, which is why nothing
        # here indexes that.
        vocabulary = check_vocabulary("antiberty", tokenizer.token_to_id)
        sequence = sequences[chain]
        found: dict[int, dict[str, float]] = {}
        for position in positions:
            residues = list(sequence)
            residues[position - 1] = "_"  # AntiBERTy's mask character
            encoded = tokenizer(["".join(residues)])
            inputs = {
                key: value.to(model.device)
                for key, value in encoded.items()
                if hasattr(value, "to")
            }
            with torch.no_grad():
                output = model(**inputs)
            logits = getattr(output, "prediction_logits", None)
            if logits is None:
                logits = getattr(output, "logits", output[0])
            row = logits[0, position]  # [CLS] is row 0, so position i is row i
            total = _logsumexp([float(v) for v in row])
            found[position] = {
                residue: float(row[vocabulary[residue]]) - total
                for residue in AMINO_ACIDS
            }
        return found

    def _load(self):
        try:
            from antiberty import AntiBERTyRunner
        except ModuleNotFoundError as exc:
            raise MissingDependencyError(
                "antiberty is not installed. `uv pip install "
                "'design-pipeline[plm]'`, which also brings torch."
            ) from exc
        return AntiBERTyRunner()


def _logsumexp(values: Sequence[float]) -> float:
    """Normalizer for a row of logits, in plain Python.

    These two libraries return logits rather than log probabilities, and the
    difference matters: log P(mut) - log P(wt) is unchanged by the normalizer,
    but a column of unnormalized logits is not comparable between positions,
    and `report` correlates the column as a whole.
    """
    import math

    if not values:
        return 0.0
    largest = max(values)
    return largest + math.log(sum(math.exp(v - largest) for v in values))


# The four. ESM-1v first of the general pair: zero-shot variant effect is the
# task it was published for.
esm1v = MaskedMarginalScorer(
    name="esm1v",
    backend=HuggingFaceMaskedLM("facebook/esm1v_t33_650M_UR90S_1"),
)
esm2_650m = MaskedMarginalScorer(
    name="esm2_650m",
    backend=HuggingFaceMaskedLM("facebook/esm2_t33_650M_UR50D"),
)
ablang2 = MaskedMarginalScorer(name="ablang2", backend=AbLang2MaskedLM())
antiberty = MaskedMarginalScorer(name="antiberty", backend=AntiBertyMaskedLM())

SEQUENCE_SCORERS: tuple[str, ...] = ("ablang2", "esm1v", "esm2_650m", "antiberty")
"""The structure-free subset by name. The general selector is the tag:

    free = [n for n, a in REGISTRY.items() if not getattr(a, "requires_structure", True)]

which also picks up `consensus`, a positional prior that needs no coordinates
either.
"""
