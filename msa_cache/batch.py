"""Fold all variants of a run in one `boltz predict` call, with the cached alignments.

Kira's Boltz adapter starts one `boltz predict` process per design, and every process reloads the
weights and prepares its input before it folds anything. Boltz accepts a directory of inputs, so one
process can fold them all and pay that cost once. The alignments come from the MsaCache, so no server
is called either.

    cache = MsaCache("runs/msa_parent")
    folded = fold_batch(cache, variants, output_dir="runs/r1/predictions/batch", seed=0)

The result is a list of Designs and a DONE manifest, in the same shape the adapter writes, so scoring
and run_report treat them alike. `design_pipeline` is not edited. The one private name relied on is
BoltzAdapter._collect, which turns Boltz's output files into a Design; it sits in `_collect` below so a
change on the other side breaks one function and its test.
"""

from __future__ import annotations

import os
import shutil
from pathlib import Path
from typing import Sequence

from .cache import CacheError, MsaCache

INPUTS = "inputs"
"""Boltz names its results folder after the input directory: boltz_results_inputs."""
DEFAULT_BOLTZ_CACHE = Path.home() / ".boltz"


def write_batch_inputs(cache: MsaCache, designs: Sequence, inputs_dir: Path | str) -> list[Path]:
    """One Boltz FASTA per design, each chain pointing at its alignment in the cache.

    Everything is checked before anything is written, so a design that is too far from the parent
    fails the batch up front instead of after the first folds.
    """
    from design_pipeline.adapters._folding import fasta_chains

    ids = [d.id for d in designs]
    if len(set(ids)) != len(ids):
        raise CacheError(f"design ids must be unique within a batch, got {sorted(i for i in set(ids) if ids.count(i) > 1)}")
    plan: list[tuple[str, list[str], list[Path]]] = []
    for design in designs:
        if Path(design.id).name != design.id or design.id.startswith("."):
            raise CacheError(f"design id {design.id!r} cannot be used as a file name")
        _, chains = fasta_chains(design.path)
        msas = cache.for_chains(chains)
        for k, msa in enumerate(msas):
            if msa is None:
                raise CacheError(
                    f"design {design.id}: chain {k + 1} does not fit a stretch of the parent's chain within "
                    f"{cache.max_mismatch_fraction:.0%} differing residues. Build a cache for this "
                    "antibody, or fold this design on its own with use_cached_msas(fallback='server')."
                )
        plan.append((design.id, chains, msas))

    inputs_dir = Path(inputs_dir)
    if inputs_dir.exists():
        shutil.rmtree(inputs_dir)
    inputs_dir.mkdir(parents=True)
    written = []
    for design_id, chains, msas in plan:
        path = inputs_dir / f"{design_id}.fasta"
        path.write_text("".join(f">{chr(ord('A') + k)}|protein|{m.resolve()}\n{c}\n" for k, (c, m) in enumerate(zip(chains, msas))))
        written.append(path)
    return written


def _boltz_cache() -> Path:
    env = os.environ.get("BOLTZ_CACHE")
    cache = Path(env) if env else DEFAULT_BOLTZ_CACHE
    if env and not cache.is_absolute():
        raise CacheError(f"BOLTZ_CACHE must be an absolute path, got {env!r}")
    return cache


def _collect(prediction_dir: Path, output_dir: Path, design):
    """Boltz's files for one design as a Design. The only use of the adapter's private method."""
    from design_pipeline.adapters._boltz import BoltzAdapter

    return BoltzAdapter()._collect(prediction_dir, output_dir, design)


def fold_batch(
    cache: MsaCache,
    designs: Sequence,
    output_dir: Path | str,
    *,
    seed: int = 0,
    extra_args: Sequence[str] = (),
    executable: str = "boltz",
) -> list:
    """Fold every design in one `boltz predict` process. Needs a GPU, like any Boltz run.

    The designs are sequence Designs (a FASTA with one record, chains separated by ':'), as the adapter
    takes. Returns one Design per input, with `parent_id` set to the input's id as the adapter does, and
    writes the DONE manifest into `output_dir`. A design that Boltz did not write a prediction for
    fails the whole call rather than returning a short list.
    """
    from design_pipeline.cache import write_manifest
    from design_pipeline.containers import stream
    from design_pipeline.errors import MissingDependencyError, RunFailedError

    if not designs:
        raise CacheError("no designs to fold")
    exe = shutil.which(executable)
    if exe is None:
        raise MissingDependencyError(f"{executable!r} is not on PATH, source the install's env.sh first")

    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    results = output_dir / f"boltz_results_{INPUTS}"
    if results.exists():
        shutil.rmtree(results)
    inputs_dir = output_dir / INPUTS
    write_batch_inputs(cache, designs, inputs_dir)

    stream(
        [
            exe,
            "predict",
            str(inputs_dir),
            "--out_dir",
            str(output_dir),
            "--cache",
            str(_boltz_cache()),
            "--output_format",
            "mmcif",
            "--write_full_pae",
            "--seed",
            str(seed),
            *extra_args,
        ]
    )

    folded = []
    for design in designs:
        prediction_dir = results / "predictions" / design.id
        if not prediction_dir.is_dir():
            raise RunFailedError(f"boltz wrote no prediction for {design.id!r} under {results / 'predictions'}")
        folded.append(_collect(prediction_dir, output_dir, design))
    write_manifest(output_dir, folded)
    return folded
