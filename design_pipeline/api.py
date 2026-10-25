"""The three model steps. Scoring and selection live in their own modules.

Every step is the same four moves: resolve the config, hash it with the seed
and the inputs, skip or run the adapter, record what came out. The steps differ
only in their ergonomic argument names and their defaults.
"""

from __future__ import annotations

from dataclasses import replace
from pathlib import Path
from typing import Any, Sequence

from .cache import cache_key, read_manifest, write_manifest
from .designs import Design
from .errors import RunFailedError
from .registry import get_adapter

PathLike = str | Path
RUNS = Path("runs")
"""Default root for output. Pass output_dir to put a campaign somewhere else."""


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
    diversity. `seed` is always passed on explicitly — left unset it is random
    upstream, which would poison the cache key.

    Model-specific arguments go through **kwargs, e.g.
    `checkpoint_protein_mpnn="proteinmpnn_v_48_020.pt"` for the variant trained
    with 0.20 A noise, which tolerates imperfect backbones better.
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
    model_name: str = "colabfold",
    seed: int = 0,
    cache: bool = True,
    extra_args: list[str] | None = None,
    **kwargs: Any,
) -> list[Design]:
    """Fold designed sequences."""
    return _run_step(
        model_name=model_name,
        inputs=list(sequences),
        output_dir=output_dir,
        seed=seed,
        cache=cache,
        extra_args=extra_args,
        kwargs=kwargs,
    )


def _run_step(
    *,
    model_name: str,
    inputs: Sequence[Design],
    output_dir: PathLike,
    seed: int,
    cache: bool,
    extra_args: list[str] | None,
    kwargs: dict[str, Any],
) -> list[Design]:
    adapter = get_adapter(model_name)
    kwargs = {k: v for k, v in kwargs.items() if v is not None}
    resolved = {
        "model_name": model_name,
        "kwargs": kwargs,
        "extra_args": list(extra_args or ()),
    }
    key = cache_key(resolved, inputs, _hashed_files(inputs, kwargs), seed)
    run_dir = Path(output_dir) / key

    if cache:
        hit = read_manifest(run_dir)
        if hit is not None:
            return [replace(d, meta={**d.meta, "cached": True}) for d in hit]

    run_dir.mkdir(parents=True, exist_ok=True)
    produced = adapter.run(
        inputs, run_dir, extra_args=list(extra_args or ()), seed=seed, **kwargs
    )
    outside = [d.path for d in produced if not _within(d.path, run_dir)]
    if outside:
        raise RunFailedError(
            f"{model_name} returned paths outside its run directory: {outside}"
        )

    # Ancestry is propagated here rather than in every adapter: an adapter sees
    # only its own inputs, and scoring needs to reach a backbone two steps up.
    parents = {d.id: d for d in inputs}
    designs = [
        replace(
            d,
            meta={
                **d.meta,
                "seed": seed,
                "cache_key": key,
                "cached": False,
                "ancestors": _ancestors(parents.get(d.parent_id or "")),
            },
        )
        for d in produced
    ]
    write_manifest(run_dir, designs)
    return designs


def _ancestors(parent: Design | None) -> tuple[str, ...]:
    if parent is None:
        return ()
    return (*tuple(parent.meta.get("ancestors") or ()), parent.id)


def _hashed_files(inputs: Sequence[Design], kwargs: dict[str, Any]) -> list[Path]:
    """Kwargs that point at an existing file, e.g. a target PDB.

    Input designs are hashed by the cache itself, by id as well as content.
    """
    missing = [d.path for d in inputs if not d.path.is_file()]
    if missing:
        raise FileNotFoundError(f"input files do not exist: {missing}")
    files = []
    for value in kwargs.values():
        if isinstance(value, (str, Path)) and Path(value).is_file():
            files.append(Path(value))
    return files


def _within(path: Path, parent: Path) -> bool:
    return parent.resolve() in path.resolve().parents
