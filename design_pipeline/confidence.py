"""How much to trust a predicted complex, and what to do about it.

`get_structure` returns this alongside the paths. It holds three kinds of
evidence, in increasing order of how much it is worth:

1. **What the model says about itself** - ipTM, pTM, mean inter-chain PAE.
   Cheap and routinely misleading. pLDDT is recorded and explicitly marked
   uninformative for this purpose: both chains of an antibody fold well
   individually whether or not the docking is right, so a high pLDDT on an
   antibody-antigen complex is close to a constant.
2. **Seed agreement** - cluster the ensemble by interface RMSD and report the
   dominant fraction. Eighteen of twenty seeds landing on one epitope is
   stronger evidence than any single ipTM, because the documented bottleneck
   is identifying the correct model among those generated, not generating one.
   It is also the only signal here that is not the model grading its own work.
3. **The decoy control** - predict the same antibody against a mismatched
   antigen on identical settings. If the decoy scores comparably, confidence
   is uninformative *for this system*, and the caller has to be able to see
   that rather than infer it.

The asymmetry matters and is stated in every docstring here: **low ipTM
reliably means "do not trust this"; high ipTM does not reliably mean "this is
correct."** Published benchmarks show these models assigning confident scores
to deliberately shuffled, non-cognate antibody-antigen pairs. So the warnings
in this module only ever fire downward. Nothing here will tell you a structure
is good.

`recommend_weights` turns that into a suggestion and nothing more: it returns a
dict and applies nothing. Deciding what to rank on stays the caller's.
"""

from __future__ import annotations

import json
import warnings
from pathlib import Path
from typing import Any, Callable, Mapping, Sequence

# Seed agreement and the decoy control are Category 2 diagnostics: they need
# no ground truth, which is exactly why they are the two confidence signals
# available on a real target.
from . import metrics
from .metrics import diagnostics
from .designs import Design

FILENAME = "confidence.json"

UNINFORMATIVE: tuple[str, ...] = ("plddt",)
"""Recorded, and not to be acted on. See the module docstring."""

WEAK_IPTM = 0.5
"""Below this, the model is telling you not to trust it, and that direction of
the signal is the reliable one."""

SPLIT_SEEDS = 0.5
"""A dominant cluster smaller than this means the ensemble has no consensus:
the predictor generated several different answers and cannot say which."""

SMALL_DECOY_GAP = 0.10
"""An ipTM gap this small between the real antigen and a mismatched one means
the score is not discriminating, whatever its absolute value."""

TRUSTED_STRUCTURE = 1.0
SEQUENCE_WHEN_TRUSTED = 0.5
"""Sequence-only scorers sit below the structure-based ones when the structure
is trustworthy: they know nothing about the interface. They come back to parity
when it is not, because then they are the only columns still meaningful."""


def summarize(
    designs: Sequence[Design],
    *,
    model_name: str,
    experimental: bool = False,
    chains: Sequence[str] | None = None,
    partners: Sequence[str] | None = None,
    decoy_iptm: float | None = None,
    seeds: bool = True,
) -> dict[str, Any]:
    """The flat confidence map for one `get_structure` run.

    Flat on purpose: it goes straight into `metrics.save`, so a notebook can
    read a campaign's confidence without knowing which predictor produced it.

    `experimental` is True for a downloaded structure, where the model-reported
    fields are absent rather than zero - a crystal structure has no ipTM, and
    writing one in would invite it to be compared against a prediction's.
    """
    representative = designs[0]
    scores = dict(representative.scores)
    summary: dict[str, Any] = {
        "model_name": model_name,
        "experimental": bool(experimental),
        "n_seeds": len(designs),
        "iptm": _number(scores.get("iptm")),
        "ptm": _number(scores.get("ptm")),
        "mean_interchain_pae": None if experimental else _interchain_pae(representative),
        # Kept for the record and flagged, because somebody will otherwise
        # read it as a docking confidence.
        "plddt": _number(scores.get("plddt")),
        "plddt_uninformative": True,
        "dominant_cluster_fraction": None,
        "n_clusters": None,
        "decoy_iptm": _number(decoy_iptm),
        "decoy_gap": None,
    }
    summary["decoy_gap"] = _number(
        diagnostics.decoy_gap(summary["iptm"], summary["decoy_iptm"])
    )

    if seeds and len(designs) > 1:
        try:
            fraction, clusters = diagnostics.seed_agreement(
                [d.path for d in designs], chains, partners
            )
        except (ValueError, KeyError) as exc:
            # A single-chain prediction has no interface to cluster on. Not a
            # failure of the run, so it is recorded as absent.
            summary["seed_agreement_error"] = str(exc)[:200]
        else:
            summary["dominant_cluster_fraction"] = fraction
            summary["n_clusters"] = clusters

    warn(summary)
    return summary


def warn(confidence: Mapping[str, Any]) -> list[str]:
    """Warn about every signal pointing *down*. Returns the warnings raised.

    One direction only. There is no branch here that says a structure looks
    good, because none of these numbers supports that conclusion.
    """
    messages = [reason for reason, triggered, _ in _triggered(confidence)]
    for message in messages:
        warnings.warn(
            f"{message}. Low confidence is the reliable direction of this "
            f"signal: treat the structure-based scorer columns as suspect and "
            f"see design_pipeline.confidence.recommend_weights.",
            RuntimeWarning,
            stacklevel=3,
        )
    return messages


RULES: tuple[tuple[str, Callable[[Mapping[str, Any]], bool], float], ...] = (
    (
        "the decoy antigen scored within "
        f"{SMALL_DECOY_GAP} ipTM of the real one, so confidence is "
        "uninformative for this system",
        lambda c: _below(c.get("decoy_gap"), SMALL_DECOY_GAP),
        0.2,
    ),
    (
        f"fewer than {SPLIT_SEEDS:.0%} of the seeds agree on one interface",
        lambda c: _below(c.get("dominant_cluster_fraction"), SPLIT_SEEDS),
        0.3,
    ),
    (
        f"ipTM is below {WEAK_IPTM}",
        lambda c: _below(c.get("iptm"), WEAK_IPTM),
        0.4,
    ),
)
"""An explicit rule table, not a model. Three readable conditions, each with
the weight it suggests for structure-based scorers; the lowest one wins. The
numbers are judgment calls rather than fitted parameters, which is exactly why
they are a visible table that can be argued with and edited."""


def recommend_weights(
    confidence: Mapping[str, Any] | None,
    models: Sequence[str] | None = None,
) -> dict[str, float]:
    """Suggested per-scorer weights. Returns a recommendation; applies nothing.

    Structure-based scorers are downweighted when the evidence says the
    structure may be wrong, and sequence-only scorers come up to parity in the
    same breath - if the complex is a fiction, they are the only columns left
    that mean anything.

    With no confidence at all (a structure from an unknown source) it declines
    to guess and weights everything equally, because a default that quietly
    favoured one class would be the library making the call this function
    exists to hand back.
    """
    names = list(models) if models is not None else _known_scorers()
    if confidence is None:
        return {name: 1.0 for name in names}

    triggered = [weight for _, _, weight in _triggered(confidence)]
    structure_weight = min(triggered) if triggered else TRUSTED_STRUCTURE
    sequence_weight = TRUSTED_STRUCTURE if triggered else SEQUENCE_WHEN_TRUSTED
    return {
        name: structure_weight if _requires_structure(name) else sequence_weight
        for name in names
    }


def weight_reasons(confidence: Mapping[str, Any] | None) -> list[str]:
    """Why `recommend_weights` suggested what it did. Empty means nothing fired."""
    if confidence is None:
        return ["no confidence recorded for this structure"]
    return [reason for reason, _, _ in _triggered(confidence)]


def write(run_dir: Path | str, confidence: Mapping[str, Any]) -> Path:
    """Write `confidence.json` beside the structure, and mirror it into
    `metrics.json` so the notebook and the harness find it the usual way."""
    run_dir = Path(run_dir)
    run_dir.mkdir(parents=True, exist_ok=True)
    path = run_dir / FILENAME
    path.write_text(json.dumps(dict(confidence), indent=1, sort_keys=True))
    metrics.save(
        run_dir,
        {
            f"confidence_{key}": value
            for key, value in confidence.items()
            if isinstance(value, (float, int, str, bool, type(None)))
        },
    )
    return path


def read(run_dir: Path | str) -> dict[str, Any] | None:
    """The recorded confidence, or None if the step never wrote one."""
    path = Path(run_dir) / FILENAME
    if not path.is_file():
        return None
    try:
        return json.loads(path.read_text())
    except json.JSONDecodeError:
        return None


def _triggered(
    confidence: Mapping[str, Any],
) -> list[tuple[str, bool, float]]:
    return [(reason, True, weight) for reason, test, weight in RULES if test(confidence)]


def _below(value: Any, limit: float) -> bool:
    """True only when there is a number and it is below the limit.

    A missing signal is not a warning. Nineteen of the twenty confidence fields
    can legitimately be absent - no seeds, no decoy, a downloaded structure -
    and treating absence as failure would fire on every honest run.
    """
    return isinstance(value, (int, float)) and value == value and float(value) < limit


def _number(value: Any) -> float | None:
    if value is None:
        return None
    try:
        number = float(value)
    except (TypeError, ValueError):
        return None
    return None if number != number else number


def _interchain_pae(prediction: Design) -> float | None:
    """Mean PAE across the interface, from the predictor's own scores file.

    Via the existing `scoring.pae_interaction`, which already averages both
    off-diagonal blocks and already knows where each adapter puts its PAE
    matrix. A second implementation would only be a second convention.
    """
    try:
        from .scoring import pae_interaction

        return _number(pae_interaction(prediction, None))
    except Exception:
        # No PAE matrix, or one whose shape does not match the structure. Not
        # a reason to fail a run that otherwise produced coordinates.
        return None


def _known_scorers() -> list[str]:
    from .adapters._plm import SEQUENCE_SCORERS
    from .adapters._scorers import MATRIX_SCORERS, PER_VARIANT_SCORERS

    return list(MATRIX_SCORERS) + list(PER_VARIANT_SCORERS) + list(SEQUENCE_SCORERS)


def _requires_structure(name: str) -> bool:
    """The registry's own tag, defaulting to True for anything unregistered.

    True is the safe default: an unknown scorer that turns out to need
    coordinates would otherwise keep full weight on a structure nobody trusts.
    """
    from .registry import REGISTRY

    adapter = REGISTRY.get(name)
    if adapter is None:
        return True
    return bool(getattr(adapter, "requires_structure", True))
