"""The public API: get a structure, score variants against it, gate the result.

Affinity maturation starts from a molecule that already exists, so the pipeline
is not a chain of generate-design-predict. It is:

    1. get_structure        obtain the complex (download, or predict)
    2. enumerate_variants   every single mutation of a chain  (variants.py)
    3. score_variants       fan out over scorers, join into one table
    4. gate                 drop confident developability failures
    5. select_pairs         combine survivors into pairs       (variants.py)

Steps 3 and 4 are a fan-out and a join, not a link in a chain: one structure
feeds several models in parallel, each writes one column keyed on `variant_id`,
and the columns are joined into one table. Nothing is ever reduced to a single
number here - no combined score, no weights, no learned ranker - because which
column to trust comes out of the validation harness, not out of a guess in a
library.

`score_variants` is called on singles and again on pairs. That reuse is the
point: a pair is an ordinary candidate with two rows, and no code path anywhere
estimates a pair's score by adding its singles' scores together.
"""

from __future__ import annotations

import warnings
from dataclasses import replace
from pathlib import Path
from typing import Any, Mapping, Sequence

from . import confidence as _confidence
from . import variants as _variants
from .adapters._gates import DEFAULT_CHECKS
from .adapters._scorers import DEFAULT_SCORERS
from .adapters._tabular import CANDIDATES, write_candidates
from .cache import cache_key, hash_file, hash_ids, read_manifest, tool_key, write_manifest
from . import metrics as _metrics
from .designs import Design
from .errors import RunFailedError
from .registry import get_adapter
from .structures import Complex, read_fasta, write_fasta
from .variants import VariantTable

PathLike = str | Path
RUNS = Path("runs")
"""Default root for output. Pass output_dir to put a campaign somewhere else."""

STRUCTURE_SOURCES = ("pdb_download", "dummy_complex")
"""Step-1 models that obtain a file rather than predict one. They take no
sequences, so naming one together with sequences is a mistake worth catching
here: the adapter would otherwise complain about a missing pdb_id, which sends
you looking at the wrong argument."""

PER_VARIANT_LIMIT = 2_000
"""Candidates a per-variant scorer will take without an explicit cap. Above
this, `max_variants` is required: a ddG model on a full 4,598-row scan is a
week of compute, and queuing it silently is worse than refusing."""


def get_structure(
    output_dir: PathLike = RUNS,
    *,
    model_name: str = "pdb_download",
    pdb_id: str | None = None,
    antibody_seq: str | Sequence[str] | None = None,
    antigen_seq: str | Sequence[str] | None = None,
    fasta: PathLike | None = None,
    chains: Sequence[str] | None = None,
    chain_ids: Sequence[str] | None = None,
    num_seeds: int | None = None,
    decoy_antigen_seq: str | Sequence[str] | None = None,
    confidence: bool = True,
    seed: int = 0,
    cache: bool = True,
    extra_args: list[str] | None = None,
    **kwargs: Any,
) -> Complex:
    """Step 1. The complex this campaign starts from.

    The default fetches it from the PDB, because experimental coordinates beat
    a prediction every time:

        complex_ = get_structure("runs/her2_v1", pdb_id="1N8Z")

    `chains` is checked against the file, which is worth doing: antibody
    structures label the heavy and light chains inconsistently, and finding out
    now beats finding out after 4,000 variants have been scored at the wrong
    positions.

    When no solved structure exists, a predictor takes the sequences instead:

        complex_ = get_structure(
            "runs/unknown_v1", model_name="alphafold3",
            antibody_seq="EVQLVESGGG...", antigen_seq="TQVCTGTDMK...",
            num_seeds=20,
        )

    Sequences are written to a FASTA - one record, chains joined with ':' - and
    handed to the adapter, so `alphafold3`, `boltz2`, `chai1` and `colabfold`
    all take the same input. Every seed a predictor writes becomes a path on
    the returned `Complex`; `complex_.path` is the tool's own top pick and
    `complex_.paths` is the ensemble behind it, because ranking of
    antibody-antigen complexes is not reliable enough to hide the rest.

    `fasta` takes those sequences from a file you already have instead, in
    whichever layout it is written - one record per chain, or one record with
    the chains joined by ':':

        complex_ = get_structure(
            "runs/unknown_v1", model_name="alphafold3",
            fasta="data/structures/my_complex.fasta",
            chain_ids=("H", "L", "A"), num_seeds=20,
        )

    Chain order is the file's order, and that is the only thing that carries
    meaning: nothing downstream distinguishes an antibody chain from an
    antigen one, so `antibody_seq` + `antigen_seq` is just a readable way of
    writing an ordered list. The file is normalised into the campaign's own
    `inputs/complex.fasta`, which is what gets hashed, so moving or renaming
    the original does not invalidate the run.

    An existing *structure* is `source=` rather than `fasta=`, and skips the
    prediction entirely:

        complex_ = get_structure("runs/local_v1", source="data/structures/x.cif")

    Every run records how much to trust the result on `complex_.confidence`,
    and beside the structure as `confidence.json`: the model's own ipTM and
    pTM, the mean inter-chain PAE, how many of the seeds agree on one
    interface, and - when `decoy_antigen_seq` is given - what the same
    antibody scores against a mismatched antigen on identical settings. That
    last one is the control worth paying for: these predictors assign
    confident scores to non-cognate pairs, so a high ipTM only means something
    if a decoy scores materially lower.

        complex_ = get_structure(
            "runs/unknown_v1", model_name="boltz2",
            antibody_seq="EVQLVESGGG...", antigen_seq="TQVCTGTDMK...",
            decoy_antigen_seq="MKTIIALSYIFCLVFA...",   # an unrelated protein
            num_seeds=20,
        )
        print(complex_.confidence["dominant_cluster_fraction"])

    `confidence=False` skips it, which is worth doing only for a large
    ensemble you have already measured: the seed clustering is pairwise.

    Model-specific arguments pass through as **kwargs, and `extra_args` is the
    raw escape hatch for a flag this package has never heard of.
    """
    output_dir = Path(output_dir)
    inputs: list[Design] = []
    chain_order: list[str] = []
    sequences: list[str] = []
    given = antibody_seq is not None or antigen_seq is not None
    if decoy_antigen_seq is not None and antibody_seq is None:
        raise ValueError(
            "decoy_antigen_seq needs antibody_seq, so the decoy keeps the same "
            "antibody chains and swaps only the antigen. With fasta= the split "
            "between antibody and antigen chains is not knowable from here."
        )
    if given and fasta is not None:
        raise ValueError(
            "pass fasta=, or antibody_seq / antigen_seq, not both: both say "
            "which chains to fold and in what order"
        )
    if given or fasta is not None:
        if fasta is not None:
            sequences = read_fasta(Path(fasta))
        else:
            sequences = _chains(antibody_seq) + _chains(antigen_seq)
            if not sequences:
                raise ValueError("antibody_seq / antigen_seq are empty")
        if model_name in STRUCTURE_SOURCES:
            raise ValueError(
                f"{model_name} obtains a structure rather than predicting one, "
                f"so it takes no sequences. Name a predictor - "
                f"model_name='alphafold3', 'boltz2', 'chai1' or 'colabfold' - "
                f"or pass pdb_id= / source= instead of the sequences."
            )
        written = write_fasta(
            output_dir / "inputs" / "complex.fasta", "complex", sequences
        )
        # A fixed id: the cache key hashes input ids as well as content, and a
        # name that moved would make every run a miss. It is also why `fasta`
        # itself is not passed down: the normalised copy is the input, so the
        # key follows the sequences rather than where they were read from.
        inputs = [
            Design(id="complex", path=written, parent_id=None, producer="get_structure")
        ]
        chain_order = list((chain_ids or ())[: len(sequences)])

    resolved = {
        "pdb_id": pdb_id,
        "chains": chains,
        "chain_ids": chain_ids,
        "num_seeds": num_seeds,
        **kwargs,
    }
    designs = _run_step(
        model_name=model_name,
        inputs=inputs,
        output_dir=output_dir,
        seed=seed,
        cache=cache,
        extra_args=extra_args,
        kwargs=dict(resolved),
    )
    run_dir = Path(output_dir) / str(designs[0].meta.get("cache_key", ""))

    recorded = None
    if confidence:
        # Read it back when it is already there: the seed clustering is
        # pairwise over the ensemble, so recomputing it on every cache hit
        # would make a resumed campaign slower than a fresh one.
        recorded = _confidence.read(run_dir)
        if recorded is None:
            recorded = _confidence.summarize(
                designs,
                model_name=model_name,
                experimental=model_name in STRUCTURE_SOURCES,
                chains=_antibody_chain_ids(chain_ids, antibody_seq),
                decoy_iptm=_decoy_iptm(
                    model_name=model_name,
                    antibody_seq=antibody_seq,
                    decoy_antigen_seq=decoy_antigen_seq,
                    output_dir=output_dir,
                    seed=seed,
                    cache=cache,
                    extra_args=extra_args,
                    resolved=resolved,
                ),
            )
            _confidence.write(run_dir, recorded)

    return Complex(
        paths=tuple(d.path for d in designs),
        output_dir=output_dir,
        confidence=recorded,
        config={
            "model_name": model_name,
            "pdb_id": pdb_id,
            "fasta": str(fasta) if fasta is not None else None,
            "chains": list(chains) if chains else None,
            "chain_order": chain_order or None,
            "num_seeds": num_seeds,
            "seed": seed,
            "cache_key": designs[0].meta.get("cache_key"),
            "kwargs": {k: v for k, v in kwargs.items() if v is not None},
        },
    )


def _antibody_chain_ids(
    chain_ids: Sequence[str] | None, antibody_seq: str | Sequence[str] | None
) -> list[str] | None:
    """Which of the named chain ids are the antibody, for seed clustering.

    Only knowable when the caller split the sequences themselves; otherwise
    None, and `metrics.seed_agreement` falls back to reading H and L as the
    antibody and everything else as the antigen.
    """
    if not chain_ids or antibody_seq is None:
        return None
    count = len(_chains(antibody_seq))
    return list(chain_ids)[:count] or None


def _decoy_iptm(
    *,
    model_name: str,
    antibody_seq: str | Sequence[str] | None,
    decoy_antigen_seq: str | Sequence[str] | None,
    output_dir: Path,
    seed: int,
    cache: bool,
    extra_args: list[str] | None,
    resolved: Mapping[str, Any],
) -> float | None:
    """ipTM for the same antibody against a mismatched antigen.

    Identical settings, identical seed, its own cache entry under
    `<output_dir>/decoy`. The point is not whether the decoy prediction is any
    good - it is that a confidence score which cannot tell the real antigen
    from an unrelated protein is telling you nothing about this system, and
    that has to be visible rather than inferred.
    """
    if decoy_antigen_seq is None:
        return None
    sequences = _chains(antibody_seq) + _chains(decoy_antigen_seq)
    fasta = write_fasta(
        output_dir / "decoy" / "inputs" / "complex.fasta", "decoy", sequences
    )
    designs = _run_step(
        model_name=model_name,
        inputs=[
            Design(id="decoy", path=fasta, parent_id=None, producer="get_structure")
        ],
        output_dir=output_dir / "decoy",
        seed=seed,
        cache=cache,
        extra_args=extra_args,
        kwargs=dict(resolved),
    )
    # .get rather than [...]: Scores raises on a miss, and a predictor that
    # reports no ipTM is a missing control, not a failed run.
    return designs[0].scores.get("iptm")


def _chains(sequences: str | Sequence[str] | None) -> list[str]:
    """One sequence or several, as a list. A ':' separates chains, as in the
    FASTA the folding adapters already read."""
    if sequences is None:
        return []
    if isinstance(sequences, str):
        return [part.strip() for part in sequences.split(":") if part.strip()]
    return [part.strip() for part in sequences if part.strip()]


def score_variants(
    complex_: Complex,
    variants: VariantTable,
    output_dir: PathLike | None = None,
    *,
    models: Sequence[str] = DEFAULT_SCORERS,
    max_variants: int | None = None,
    structure: PathLike | None = None,
    cache: bool = True,
    extra_args: list[str] | None = None,
    **kwargs: Any,
) -> VariantTable:
    """Step 3. One column per model, joined on `variant_id`.

        scores = score_variants(complex_, variants, models=["antifold", "esm_if"])

    Each model is an adapter that writes `variant_id,score[,raw]` and knows
    nothing about the others; this function runs them, caches them
    individually, and joins their columns into `score_<model>` plus
    `score_<model>_raw`. `score_` is normalized at the adapter boundary so
    higher always means better binding whatever the tool's own convention, and
    the raw value is kept so a convention we got backwards is recoverable
    without re-running anything. The columns are never averaged together.

    The two cost classes are explicit rather than hidden. A matrix scorer -
    `antifold`, `esm_if`, `consensus` - sees the whole table, because one
    forward pass over the complex answers every single mutation. A per-variant
    scorer - `gearbind`, `rosetta_ddg`, `foldx` - builds a mutant structure per
    candidate, so `max_variants` caps what it is given and above
    PER_VARIANT_LIMIT candidates an uncapped one raises instead of queueing a
    week of compute. The cap applies only to those: capping the cheap scorers
    would throw away columns that cost nothing.

    Caching is per scorer, on (structure, candidate set, model, kwargs). Adding
    a model to `models` does not invalidate the ones that already ran.

    `**kwargs` and `extra_args` go to every model named, so a flag one model
    rejects is an error for the whole call; the remedy is one call per model,
    which joins into the same table anyway.
    """
    out = Path(output_dir) if output_dir is not None else complex_.output_dir / "scores"
    target = Path(structure) if structure is not None else complex_.path

    import pandas as pd

    frame = variants.frame
    results = _fan_out(
        structure=target,
        frame=frame,
        models=models,
        output_dir=out,
        max_variants=max_variants,
        cache=cache,
        extra_args=extra_args,
        kwargs=kwargs,
    )

    for model, csv_path in results:
        answers = _read_answers(csv_path, model, ("score", "raw"))
        frame = frame.drop(
            columns=[f"score_{model}", f"score_{model}_raw"], errors="ignore"
        )
        frame = frame.merge(
            answers.rename(
                columns={"score": f"score_{model}", "raw": f"score_{model}_raw"}
            ),
            on="variant_id",
            how="left",
        )

    return _variants.write_table(
        frame,
        out / f"scored_{variants.path.stem}.parquet",
        output_dir=variants.output_dir,
        structure=target,
        config={
            "models": list(models),
            "max_variants": max_variants,
            "structure": str(target),
            "from": str(variants.path),
        },
    )


def gate(
    variants: VariantTable,
    checks: Sequence[str] = DEFAULT_CHECKS,
    output_dir: PathLike | None = None,
    *,
    complex_: Complex | None = None,
    structure: PathLike | None = None,
    max_variants: int | None = None,
    cache: bool = True,
    extra_args: list[str] | None = None,
    **kwargs: Any,
) -> VariantTable:
    """Step 4. Drop the candidates that fail a developability check.

        survivors = gate(scores, checks=["liabilities"])

    Each check answers a boolean and a reason per candidate. A failure is
    dropped, never traded off against binding score: a 10-fold affinity gain is
    not worth a new glycosylation sequon, and the only way to be sure the
    trade is not made is to make it impossible.

    The annotated table - every candidate, every `pass_<check>` and
    `reason_<check>`, failures included - is written next to the survivors, so
    a gate that dropped too much can be re-tuned without re-running the tools.

    Checks follow the same two cost classes as scorers, so `max_variants` caps
    the structure-building ones (`thermompnn`, `aggrescan3d`) the same way.
    `liabilities` is pure sequence regex and runs on everything.
    """
    out = Path(output_dir) if output_dir is not None else variants.output_dir / "gated"
    target = Path(structure) if structure is not None else (
        complex_.path if complex_ is not None else variants.structure
    )
    if target is None:
        raise ValueError(
            "gate needs the complex: the checks read the chain sequences and "
            "the coordinates. Pass complex_= (or structure=), or use a table "
            "that came from enumerate_variants, which carries it."
        )

    frame = variants.frame
    results = _fan_out(
        structure=target,
        frame=frame,
        models=checks,
        output_dir=out,
        max_variants=max_variants,
        cache=cache,
        extra_args=extra_args,
        kwargs=kwargs,
    )

    passing = None
    for check, csv_path in results:
        answers = _read_answers(csv_path, check, ("passed", "reason"))
        frame = frame.drop(
            columns=[f"pass_{check}", f"reason_{check}"], errors="ignore"
        )
        frame = frame.merge(
            answers.rename(
                columns={"passed": f"pass_{check}", "reason": f"reason_{check}"}
            ),
            on="variant_id",
            how="left",
        )
        # A candidate the check said nothing about - one a cap left out, or one
        # the tool could not model - is not a confident failure, so it passes.
        verdict = frame[f"pass_{check}"].fillna(True).astype(bool)
        frame[f"pass_{check}"] = verdict
        # A check whose every reason was empty reads back as an all-NaN float
        # column, so the cast comes before the fill.
        frame[f"reason_{check}"] = frame[f"reason_{check}"].astype(object).fillna("")
        passing = verdict if passing is None else (passing & verdict)

    annotated = _variants.write_table(
        frame,
        out / f"annotated_{variants.path.stem}.parquet",
        output_dir=variants.output_dir,
        structure=target,
        config={"checks": list(checks), "from": str(variants.path)},
    )
    survivors = frame if passing is None else frame[passing]
    if not len(survivors):
        warnings.warn(
            f"gate dropped every one of the {variants.candidates} candidates. "
            f"The reasons are in {annotated.path}; the developability "
            f"predictors are weak, so a threshold that strict is usually the "
            f"cause rather than the variants.",
            RuntimeWarning,
            stacklevel=2,
        )
    return _variants.write_table(
        survivors,
        out / f"survivors_{variants.path.stem}.parquet",
        output_dir=variants.output_dir,
        structure=target,
        config={
            "checks": list(checks),
            "annotated": str(annotated.path),
            "from": str(variants.path),
        },
    )


def _fan_out(
    *,
    structure: Path,
    frame: Any,
    models: Sequence[str],
    output_dir: Path,
    max_variants: int | None,
    cache: bool,
    extra_args: list[str] | None,
    kwargs: Mapping[str, Any],
) -> list[tuple[str, Path]]:
    """Run each named adapter over the candidate set; return its result CSV.

    The shared core of `score_variants` and `gate`. Both hand one structure and
    one candidate list to several adapters and join one column back per
    adapter, which is why neither needs a Protocol of its own.
    """
    names = _unique(models)
    if not names:
        raise ValueError("name at least one model")
    adapters = {name: get_adapter(name) for name in names}  # validated up front
    ids = list(dict.fromkeys(frame["variant_id"]))
    structure_hash = hash_file(structure)
    resolved_kwargs = {k: v for k, v in kwargs.items() if v is not None}

    results: list[tuple[str, Path]] = []
    for name in names:
        adapter = adapters[name]
        per_variant = bool(getattr(adapter, "per_variant", False))
        subset = ids
        if per_variant:
            if max_variants is None and len(ids) > PER_VARIANT_LIMIT:
                raise ValueError(
                    f"{name} builds a structure per candidate and there are "
                    f"{len(ids)} of them. Pass max_variants=<n> to run it on a "
                    f"subset (the first n rows of the table, so sort on the "
                    f"column you trust first), or drop {name!r} from this call."
                )
            if max_variants is not None:
                subset = ids[:max_variants]

        resolved = {"kwargs": resolved_kwargs, "extra_args": list(extra_args or ())}
        key = tool_key(
            structure_hash=structure_hash,
            candidate_hash=hash_ids(subset),
            model_name=name,
            resolved=resolved,
        )
        run_dir = output_dir / key

        if cache:
            hit = read_manifest(run_dir)
            if hit is not None:
                results.append((name, hit[0].path))
                continue

        run_dir.mkdir(parents=True, exist_ok=True)
        wanted = set(subset)
        candidates = write_candidates(
            run_dir / CANDIDATES,
            (row._asdict() for row in frame.itertuples() if row.variant_id in wanted),
        )
        produced = adapter.run(
            [
                Design(
                    id=structure.stem,
                    path=structure,
                    parent_id=None,
                    producer="get_structure",
                )
            ],
            run_dir,
            candidates=str(candidates),
            extra_args=list(extra_args or ()),
            **resolved_kwargs,
        )
        if not produced:
            raise RunFailedError(f"{name} wrote no result table under {run_dir}")
        outside = [d.path for d in produced if not _within(d.path, run_dir)]
        if outside:
            raise RunFailedError(
                f"{name} returned paths outside its run directory: {outside}"
            )
        write_manifest(
            run_dir,
            [
                replace(d, meta={**d.meta, "cache_key": key, "candidates": len(subset)})
                for d in produced
            ],
        )
        # Every step that produces numbers leaves a metrics.json, so a campaign
        # can be swept with metrics.load_all() without knowing its layout.
        _metrics.save(
            run_dir,
            {
                "model_name": name,
                "candidates": len(subset),
                "per_variant": per_variant,
                "capped": max_variants is not None and len(subset) < len(ids),
            },
        )
        results.append((name, produced[0].path))
    return results


def _read_answers(path: Path, model: str, columns: Sequence[str]) -> Any:
    """One adapter's result CSV as a frame of variant_id plus `columns`.

    Validated here rather than trusted: a duplicated id would multiply rows on
    the join, and an id nobody asked for means the join key convention has
    drifted, which is worth an error rather than a silently dropped column.
    """
    import pandas as pd

    answers = pd.read_csv(path)
    missing = [c for c in ("variant_id", *columns) if c not in answers.columns]
    if missing:
        raise RunFailedError(
            f"{model} wrote {path} without column(s) {missing}; expected "
            f"variant_id,{','.join(columns)}"
        )
    duplicated = answers["variant_id"].duplicated()
    if duplicated.any():
        repeated = answers.loc[duplicated, "variant_id"].iloc[0]
        raise RunFailedError(
            f"{model} wrote {int(duplicated.sum())} duplicate rows in {path} "
            f"(e.g. {repeated!r}); one row per candidate, not per mutation"
        )
    return answers[["variant_id", *columns]]


def _unique(models: Sequence[str]) -> list[str]:
    names = [models] if isinstance(models, str) else list(models)
    duplicates = sorted({n for n in names if names.count(n) > 1})
    if duplicates:
        listed = ", ".join(repr(d) for d in duplicates)
        raise ValueError(
            f"{listed} named more than once; every model runs over the same "
            f"candidates, so a repeat would only overwrite its own column"
        )
    return names


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
    """Resolve the config, hash it with the inputs, skip or run, record.

    Unchanged from the first version of this package: it drives `get_structure`
    and the de novo branch in `denovo.py`. The scoring fan-out does not go
    through it, because its key is per scorer rather than per step.
    """
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
    if not produced:
        # Without this the empty list reaches `get_structure`, which indexes
        # designs[0] for the cache key and reports an IndexError. The usual
        # cause is a predictor handed no sequences at all.
        detail = (
            "; it was handed no input, so there was nothing to fold"
            if not inputs
            else ""
        )
        raise RunFailedError(
            f"{model_name} returned no designs under {run_dir}{detail}"
        )
    outside = [d.path for d in produced if not _within(d.path, run_dir)]
    if outside:
        raise RunFailedError(
            f"{model_name} returned paths outside its run directory: {outside}"
        )

    # Ancestry is propagated here rather than in every adapter: an adapter sees
    # only its own inputs, and the de novo branch needs to reach a backbone two
    # steps up.
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
