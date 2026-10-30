"""The de novo branch: RFdiffusion and MPNN, kept but off the default path.

These were the first two steps of this package when it was a de novo design
pipeline. Affinity maturation does not need them - you do not rebuild a loop to
test a point mutation - so they are no longer part of the public surface.
`diffuse_loops` is the one entry point that reaches them, and what it returns
is an ordinary variant table scored by the ordinary `score_variants`.

A warning about that table, and it is not a small one: a rebuilt CDR differs
from the parent at ten or twenty positions at once. Every scorer in this
package was validated on single and double mutants, so their numbers on these
candidates are extrapolation far outside the edit-distance range where anybody
has measured whether they work. Treat the ranking as a hypothesis generator
and nothing more; the validation harness cannot tell you about this regime
either, because the published affinity datasets do not cover it.
"""

from __future__ import annotations

from pathlib import Path
from typing import Any, Sequence

from . import numbering, variants as _variants
from .adapters._folding import fasta_chains
from .api import RUNS, PathLike, _run_step
from .designs import Design
from .structures import Complex, chain_ids, chain_residue_ids, chain_sequence
from .variants import VariantTable

IMGT_CDR: dict[int, tuple[int, int]] = {
    1: (27, 38),
    2: (56, 65),
    3: (105, 117),
}
"""IMGT CDR definitions, inclusive. The reason `cdr="H3"` can be resolved
without the caller counting residues - and the reason ANARCI is needed for it."""


def make_protein_backbones(
    inputs: PathLike | Sequence[Design] | None = None,
    output_dir: PathLike = RUNS / "backbones",
    *,
    model_name: str = "rfdiffusion",
    input_pdb: PathLike | None = None,
    num_designs: int = 1,
    contigs: str | None = None,
    partial_T: int | None = None,
    seed: int = 0,
    cache: bool = True,
    extra_args: list[str] | None = None,
    **kwargs: Any,
) -> list[Design]:
    """Generate backbones.

    `inputs` is either a path to a target/motif PDB (fresh generation) or a
    `list[Design]` to diffuse around (partial diffusion, one invocation per
    design, each as `inference.input_pdb`). Same function, same adapter.

    `partial_T` is a raw timestep count relative to `diffuser.T`, not a 0-1
    noise fraction: published values like 80 assume the original 200-step
    model, so at today's default T=50 the equivalent is roughly 20. Passing
    designs without `partial_T` is legal and means "generate afresh from each",
    which is rarely what you want.
    """
    designs: list[Design] = []
    if isinstance(inputs, (str, Path)):
        input_pdb = inputs
    elif inputs:
        designs = list(inputs)
    if input_pdb is not None and not Path(input_pdb).is_file():
        raise FileNotFoundError(f"input PDB does not exist: {input_pdb}")

    return _run_step(
        model_name=model_name,
        inputs=designs,
        output_dir=output_dir,
        seed=seed,
        cache=cache,
        extra_args=extra_args,
        kwargs={
            "input_pdb": input_pdb,
            "num_designs": num_designs,
            "contigs": contigs,
            "partial_T": partial_T,
            **kwargs,
        },
    )


def design_sequences(
    backbones: Sequence[Design],
    output_dir: PathLike = RUNS / "sequences",
    *,
    model_name: str = "protein_mpnn",
    num_seqs: int = 1,
    temperature: float = 0.1,
    seed: int = 0,
    cache: bool = True,
    extra_args: list[str] | None = None,
    **kwargs: Any,
) -> list[Design]:
    """Design sequences for backbones. One Design per sampled sequence.

    `model_name` picks the MPNN variant: `protein_mpnn` (the default),
    `ligand_mpnn` or `soluble_mpnn`.

    `temperature` 0.1 is the tool's own default; raise it toward 0.3 for more
    diversity. `seed` is always passed on explicitly - left unset it is random
    upstream, which would poison the cache key.
    """
    return _run_step(
        model_name=model_name,
        inputs=list(backbones),
        output_dir=output_dir,
        seed=seed,
        cache=cache,
        extra_args=extra_args,
        kwargs={"num_seqs": num_seqs, "temperature": temperature, **kwargs},
    )


def predict_structure(
    sequences: Sequence[Design],
    output_dir: PathLike = RUNS / "predictions",
    *,
    model_name: str | list[str] = "chai",
    seed: int = 0,
    cache: bool = True,
    extra_args: list[str] | None = None,
    **kwargs: Any,
) -> list[Design]:
    """Fold designed sequences.

    Reachable for the de novo branch; `get_structure` is the supported way to
    fold an antibody-antigen complex, because it returns the whole ensemble
    rather than only the top-ranked model.

    `model_name` is one model or a list of them. A list folds the same
    sequences with each named model in turn and concatenates the results, so N
    sequences x M models gives N*M predictions: each keeps the sequence it was
    folded from as `parent_id` and names the model that made it in `producer`.
    The models run one after another and the first failure propagates, so a
    partial result is never returned.
    """
    predictions: list[Design] = []
    for name in _model_names(model_name):
        predictions.extend(
            _run_step(
                model_name=name,
                inputs=list(sequences),
                output_dir=output_dir,
                seed=seed,
                cache=cache,
                extra_args=extra_args,
                kwargs=kwargs,
            )
        )
    return predictions


def _model_names(model_name: str | list[str]) -> list[str]:
    """The models to run, validated before any of them starts.

    A plain string is the one-model case and stays indistinguishable from a
    single-element list. Caching is per (model, inputs) already, because the
    model is part of the resolved config, so a repeated name would run once
    and hand back the same predictions twice.
    """
    if isinstance(model_name, str):
        return [model_name]
    names = list(model_name)
    if not names:
        raise ValueError("model_name is an empty list: name at least one model")
    duplicates = sorted({n for n in names if names.count(n) > 1})
    if duplicates:
        listed = ", ".join(repr(d) for d in duplicates)
        raise ValueError(
            f"model_name lists {listed} more than once; every model runs over the "
            "same inputs, so a repeat would only duplicate its predictions"
        )
    return names


def diffuse_loops(
    complex_: Complex,
    cdr: str = "H3",
    num_designs: int = 8,
    output_dir: PathLike | None = None,
    *,
    span: tuple[int, int] | None = None,
    contigs: str | None = None,
    partial_T: int = 20,
    num_seqs: int = 1,
    backbone_model: str = "rfdiffusion",
    sequence_model: str = "protein_mpnn",
    seed: int = 0,
    cache: bool = True,
    extra_args: list[str] | None = None,
    **kwargs: Any,
) -> VariantTable:
    """Rebuild one CDR by partial diffusion, as a variant table.

        loops = diffuse_loops(complex_, cdr="H3", num_designs=100)
        scores = score_variants(complex_, loops, models=["antifold"])

    Not part of the default path. The loop is rebuilt at its existing length so
    the designed sequence can be diffed against the parent position by
    position, which is what makes the result an ordinary variant table and
    scorable by the ordinary scorers. See the module docstring on how far
    outside their validated range those scores are.

    `cdr` is a chain id and an IMGT CDR number, so "H3" and "L1" are the usual
    arguments and resolving them needs ANARCI. Pass `span=(start, end)` - 1
    based, into the chain as given - to say it directly, or `contigs=` to hand
    RFdiffusion its own contig string and bypass both.
    """
    out = Path(output_dir) if output_dir is not None else complex_.output_dir / "loops"
    chain = cdr[:1]
    sequence = chain_sequence(complex_.path, chain)
    first, last = span if span is not None else _cdr_span(complex_.path, cdr)

    backbones = make_protein_backbones(
        [Design(id="parent", path=complex_.path, parent_id=None, producer="get_structure")],
        out / "backbones",
        model_name=backbone_model,
        num_designs=num_designs,
        contigs=contigs or _contigs(complex_.path, chain, first, last),
        partial_T=partial_T,
        seed=seed,
        cache=cache,
        extra_args=extra_args,
        **kwargs,
    )
    designed = design_sequences(
        backbones,
        out / "sequences",
        model_name=sequence_model,
        num_seqs=num_seqs,
        seed=seed,
        cache=cache,
    )

    order = chain_ids(complex_.path)
    if chain not in order:
        raise KeyError(f"chain {chain!r} is not in {complex_.path}")
    index = order.index(chain)
    labels = numbering.imgt_labels(sequence)

    rows: list[dict[str, Any]] = []
    for design in designed:
        _name, chains = fasta_chains(design.path)
        if index >= len(chains):
            raise ValueError(
                f"{design.path} has {len(chains)} chain(s) and chain {chain!r} "
                f"is number {index + 1} of the complex; the backbone step "
                f"dropped a chain, so the diff would compare the wrong one"
            )
        mutations = _diff(sequence, chains[index], chain, design.path)
        if not mutations:
            continue  # the model re-derived the parent loop; not a variant
        rows += _variants._rows_for(mutations, {chain: labels})

    if not rows:
        raise ValueError(
            f"all {len(designed)} designs re-derived the parent sequence at "
            f"{chain}{first}-{last}; raise partial_T or num_seqs"
        )

    import pandas as pd

    return _variants.write_table(
        pd.DataFrame(rows),
        out / "loops.parquet",
        output_dir=complex_.output_dir,
        structure=complex_.path,
        config={
            "cdr": cdr,
            "span": [first, last],
            "num_designs": num_designs,
            "partial_T": partial_T,
            "backbone_model": backbone_model,
            "sequence_model": sequence_model,
            # Recorded so a reader of the table knows these are not point
            # mutants: the scorers were never validated this far out.
            "edit_distance": "loop-scale",
        },
    )


def _cdr_span(structure: Path, cdr: str) -> tuple[int, int]:
    """(first, last) chain positions of an IMGT CDR, 1-based."""
    chain, index = cdr[:1], cdr[1:]
    try:
        low, high = IMGT_CDR[int(index)]
    except (KeyError, ValueError):
        raise ValueError(
            f"{cdr!r} is not a CDR: expected a chain id and 1, 2 or 3, "
            f"e.g. 'H3'"
        ) from None
    labels = numbering.imgt_labels(chain_sequence(structure, chain))
    positions = [
        position
        for position, label in enumerate(labels, start=1)
        if label and low <= int("".join(c for c in label if c.isdigit())) <= high
    ]
    if not positions:
        raise ValueError(
            f"could not locate {cdr} in chain {chain} of {structure}. ANARCI "
            f"{'is not installed' if not numbering.available() else 'found no variable domain'}; "
            f"pass span=(first, last) in chain positions instead."
        )
    return positions[0], positions[-1]


def _contigs(structure: Path, chain: str, first: int, last: int) -> str:
    """An RFdiffusion contig string that rebuilds one span and keeps the rest.

    The rebuilt segment is pinned to its original length (`n-n`), because a
    variant table is a position-by-position diff against the parent and an
    indel would leave nothing to diff. Chains other than the target are kept
    whole as context.
    """
    segments: list[str] = []
    for other in chain_ids(structure):
        residues = chain_residue_ids(structure, other)
        if any(not r.lstrip("-").isdigit() for r in residues):
            raise ValueError(
                f"chain {other} of {structure} has insertion codes "
                f"(e.g. {next(r for r in residues if not r.isdigit())!r}) and "
                f"RFdiffusion contigs are numeric ranges. Renumber the "
                f"structure, or pass contigs= directly."
            )
        if other != chain:
            segments.append(f"{other}{residues[0]}-{residues[-1]}")
            continue
        length = last - first + 1
        parts = []
        if first > 1:
            parts.append(f"{chain}{residues[0]}-{residues[first - 2]}")
        parts.append(f"{length}-{length}")
        if last < len(residues):
            parts.append(f"{chain}{residues[last]}-{residues[-1]}")
        segments.append("/".join(parts))
    return "[" + "/0 ".join(segments) + "]"


def _diff(
    wild: str, designed: str, chain: str, source: Path
) -> list[tuple[str, int, str, str]]:
    """Every position where the designed chain differs from the parent."""
    if len(wild) != len(designed):
        raise ValueError(
            f"{source} designed {len(designed)} residues for chain {chain} and "
            f"the parent has {len(wild)}. diffuse_loops rebuilds the loop at "
            f"its original length precisely so this diff is possible; a "
            f"different length means the contig string changed it."
        )
    return [
        (chain, position, wt, mut)
        for position, (wt, mut) in enumerate(zip(wild, designed), start=1)
        if wt != mut
    ]
