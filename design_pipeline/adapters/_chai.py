"""Chai-1 through its Python API.

In process, unlike the other adapters here, and for one reason: the `chai-lab
fold` command line writes pTM and ipTM to an npz and throws the PAE matrix
away, while `run_inference` returns it. `pae_interaction` is the metric this
campaign ranks on and Chai-1 is the default folding model, so the adapter
imports chai_lab (and with it torch) inside `run` and keeps the whole
prediction. Nothing is imported at module import time.

Apache 2.0 for code and for weights, commercial use included. MSA-free by
default, which is the fast path for screening; `use_msa_server=True` routes to
the public ColabFold MMseqs2 server. The model runs in bfloat16, so pre-Ampere
cards are out, and `require_gpu` says so before any weights download.

Weights land in `$CHAI_DOWNLOADS_DIR` on first use; unset, chai_lab writes them
into a `downloads/` directory inside its own installation, which is rarely what
you want on a cluster.
"""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
from typing import Any, Callable, Sequence

from ..containers import make_id
from ..designs import Design
from ..errors import MissingDependencyError, RunFailedError
from ._folding import (
    cif_to_pdb,
    fasta_chains,
    first_match,
    require_gpu,
    work_dir,
    write_scores_json,
)

PRODUCER = "chai"


def load_run_inference() -> Callable[..., Any]:
    """`chai_lab.chai1.run_inference`, imported on use.

    Importing design_pipeline must not import a model, and this import brings
    torch with it.
    """
    import os

    try:
        from chai_lab.chai1 import run_inference
    except ModuleNotFoundError as exc:  # pragma: no cover - needs the extra
        # The one Env B tool with no path to point at. Every other adapter here
        # shells out, so DESIGN_PIPELINE_<TOOL> is enough to reach it from Env
        # A; this one imports, so the variable can only tell you where the
        # environment is - you have to be running inside it. Said explicitly,
        # because a set variable and a failing import looks like a bug
        # otherwise.
        configured = os.environ.get("DESIGN_PIPELINE_CHAI", "").strip()
        where = (
            f" $DESIGN_PIPELINE_CHAI is set to {configured}, but chai runs in "
            "process: run this from inside that environment rather than "
            "pointing at it from here."
            if configured
            else ""
        )
        raise MissingDependencyError(
            "chai_lab is not installed in this interpreter. Install the "
            "optional extra with `uv pip install -e '.[chai]'` (or `pip "
            "install chai_lab`) and run from inside that environment, or use "
            f"model_name='dummy_chai' to run the loop without a folding model.{where}"
        ) from exc
    return run_inference


def write_chai_fasta(path: Path, sequence: Design) -> Path:
    """One record per chain, in Chai-1's `entity|name=label` header format."""
    name, chains = fasta_chains(sequence.path)
    records = [
        f">protein|name={name}-{chr(ord('A') + index)}\n{chain}\n"
        for index, chain in enumerate(chains)
    ]
    path.write_text("".join(records))
    return path


@dataclass(frozen=True)
class ChaiAdapter:
    """Satisfies the ModelAdapter protocol."""

    extra: str = "chai"

    def run(
        self,
        inputs: Sequence[Design],
        output_dir: Path,
        *,
        extra_args: Sequence[str] | None = None,
        use_msa_server: bool = False,
        num_diffusion_samples: int | None = None,
        seed: int = 0,
        **kwargs: Any,
    ) -> list[Design]:
        require_gpu("chai", self.extra, bfloat16=True)
        if extra_args:
            raise ValueError(
                "chai runs in process and has no command line for extra_args="
                f"{list(extra_args)}. Pass run_inference arguments as keyword "
                "arguments instead, e.g. num_trunk_recycles=5."
            )
        run_inference = load_run_inference()
        # num_diffusion_samples is this package's name for the tool's
        # num_diffn_samples; left unset, the tool's own default stands.
        options = {"num_diffn_samples": num_diffusion_samples, **kwargs}
        options = {key: value for key, value in options.items() if value is not None}

        designs: list[Design] = []
        for sequence in inputs:
            work = work_dir(output_dir, sequence)
            fasta = write_chai_fasta(work / "input.fasta", sequence)
            # run_inference asserts its output directory is empty, so the
            # converted input stays in the parent rather than beside it.
            predictions = work / "prediction"
            candidates = run_inference(
                fasta_file=fasta,
                output_dir=predictions,
                use_msa_server=use_msa_server,
                seed=seed,
                **options,
            )
            designs.append(
                self._best(candidates.sorted(), predictions, output_dir, sequence)
            )

        if not designs:
            raise RunFailedError(f"chai returned no predictions under {output_dir}")
        return designs

    def _best(
        self,
        candidates: Any,
        work: Path,
        output_dir: Path,
        sequence: Design,
    ) -> Design:
        """The highest-ranked candidate, as one Design.

        Chai-1 writes one mmCIF per diffusion sample; `sorted()` puts the best
        aggregate score first and only that one becomes a design, which is how
        the ColabFold adapter treats rank_001 too. The rest stay on disk.
        """
        if not candidates.cif_paths:
            raise RunFailedError(
                f"chai wrote no mmCIF for {sequence.id!r} under {work}"
            )
        cif = Path(candidates.cif_paths[0])
        pdb = cif_to_pdb(cif)
        ranking = candidates.ranking_data[0]
        ptm = float(ranking.ptm_scores.complex_ptm)
        iptm = float(ranking.ptm_scores.interface_ptm)
        # Chai-1 reports pLDDT per token on a 0-1 scale; ColabFold's key is 0-100.
        plddt = [100.0 * float(value) for value in candidates.plddt[0]]
        pae = [[float(value) for value in row] for row in candidates.pae[0]]
        scores_json = write_scores_json(
            work / f"{pdb.stem}_scores.json",
            plddt=plddt,
            pae=pae,
            ptm=ptm,
            iptm=iptm,
        )
        meta: dict[str, Any] = {
            "cif": str(cif),
            "scores_json": str(scores_json),
        }
        # scores.model_idx_N.npz beside pred.model_idx_N.cif: the ranking order
        # does not rename the files, so the index comes off the mmCIF name.
        raw = first_match(work, f"scores.{cif.stem.removeprefix('pred.')}.npz")
        if raw is not None:
            meta["raw_scores"] = str(raw)
        return Design(
            id=make_id(pdb, output_dir),
            path=pdb,
            parent_id=sequence.id,
            producer=PRODUCER,
            scores={
                "plddt": sum(plddt) / len(plddt),
                "ptm": ptm,
                "iptm": iptm,
                # Chai-1's own ranking score: 0.2 pTM + 0.8 ipTM, minus a clash
                # penalty. Prefixed, because it is not a ColabFold key.
                "chai_aggregate_score": float(ranking.aggregate_score),
            },
            meta=meta,
        )
