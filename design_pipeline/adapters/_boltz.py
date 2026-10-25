"""Boltz-2 through the `boltz predict` command line.

A subprocess, as for ColabFold and MPNN: its stdout is the progress report we
want streamed, and it keeps torch out of this process.

MIT licensed, code and weights. Boltz expects MSAs, so `use_msa_server`
defaults to True and routes to the public ColabFold MMseqs2 server; turn it off
and pass precomputed alignments through the FASTA header
(`>A|protein|/path/to/msa.a3m`) with extra_args if you have them.

Weights, the CCD and the molecule archive download on first prediction into
`$BOLTZ_CACHE` (default `~/.boltz`), several GB in total. Boltz 2.x has no
`boltz download` subcommand - `predict` fetches what is missing itself - so
this adapter creates that directory, passes it explicitly and says, once,
before a run that would download.

Boltz-2 also predicts binding affinity for protein-ligand inputs. Those
numbers come through under `boltz_` names and nothing here is built on them.
"""

from __future__ import annotations

import json
import os
import shutil
import sys
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any, Mapping, Sequence

import numpy as np

from ..containers import format_kwargs, make_id, stream
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

PRODUCER = "boltz"
DEFAULT_CACHE = Path.home() / ".boltz"
CHECKPOINT = "boltz2_conf.ckpt"
"""Inside the cache. Its absence is what makes the first run a download."""


def write_boltz_fasta(path: Path, sequence: Design) -> Path:
    """One record per chain, in Boltz's `chain_id|entity_type` header format."""
    _, chains = fasta_chains(sequence.path)
    records = [
        f">{chr(ord('A') + index)}|protein\n{chain}\n"
        for index, chain in enumerate(chains)
    ]
    path.write_text("".join(records))
    return path


def parse_confidence(path: Path) -> dict[str, float]:
    """Boltz's confidence summary, under the keys ColabFold uses.

    `complex_plddt` is a fraction, ColabFold's `plddt` is 0-100. Everything
    past the shared three keeps a `boltz_` prefix.
    """
    try:
        payload = json.loads(path.read_text())
    except (OSError, json.JSONDecodeError):
        return {}
    scores: dict[str, float] = {}
    if payload.get("complex_plddt") is not None:
        scores["plddt"] = 100.0 * float(payload["complex_plddt"])
    for key in ("ptm", "iptm"):
        if payload.get(key) is not None:
            scores[key] = float(payload[key])
    for key in ("confidence_score", "ligand_iptm", "protein_iptm", "complex_iplddt"):
        if payload.get(key) is not None:
            scores[f"boltz_{key}"] = float(payload[key])
    return scores


def parse_affinity(path: Path) -> dict[str, float]:
    """Boltz-2's affinity head, for protein-ligand inputs only."""
    try:
        payload = json.loads(path.read_text())
    except (OSError, json.JSONDecodeError):
        return {}
    return {
        f"boltz_{key}": float(value)
        for key, value in payload.items()
        if isinstance(value, (int, float))
    }


@dataclass(frozen=True)
class BoltzAdapter:
    """Satisfies the ModelAdapter protocol."""

    executable: str = "boltz"
    extra: str = "boltz"
    aliases: Mapping[str, str] = field(default_factory=dict)
    """Boltz's own options are underscored, so this is empty."""

    def run(
        self,
        inputs: Sequence[Design],
        output_dir: Path,
        *,
        extra_args: Sequence[str] | None = None,
        use_msa_server: bool = True,
        seed: int = 0,
        recycling_steps: int | None = None,
        diffusion_samples: int | None = None,
        **kwargs: Any,
    ) -> list[Design]:
        require_gpu("boltz", self.extra)
        exe = self._resolve()
        cache = self._cache()
        args = format_kwargs(
            {
                "seed": seed,
                "recycling_steps": recycling_steps,
                "diffusion_samples": diffusion_samples,
                **kwargs,
            },
            "flag",
            self.aliases,
        )
        if use_msa_server:
            args.append("--use_msa_server")

        designs: list[Design] = []
        for sequence in inputs:
            work = work_dir(output_dir, sequence)
            fasta = write_boltz_fasta(work / f"{sequence.id}.fasta", sequence)
            stream(
                [
                    exe,
                    "predict",
                    str(fasta),
                    "--out_dir",
                    str(work),
                    "--cache",
                    str(cache),
                    "--output_format",
                    "mmcif",
                    # The PAE matrix is not in the confidence summary, and
                    # pae_interaction needs it.
                    "--write_full_pae",
                    *args,
                    *(extra_args or ()),
                ]
            )
            designs.append(self._collect(work, output_dir, sequence))

        if not designs:
            raise RunFailedError(f"boltz returned no predictions under {output_dir}")
        return designs

    def _resolve(self) -> str:
        found = shutil.which(self.executable)
        if found:
            return found
        raise MissingDependencyError(
            f"{self.executable!r} is not on PATH. Install the optional extra "
            "with `uv pip install -e '.[boltz]'` (or `pip install "
            "'boltz[cuda]'`), or use model_name='dummy_boltz' to run the loop "
            "without a folding model."
        )

    def _cache(self) -> Path:
        """Where Boltz keeps its weights, created before the run needs it.

        BOLTZ_CACHE is the tool's own variable and must be absolute; this
        respects it rather than inventing a cache of our own.
        """
        env = os.environ.get("BOLTZ_CACHE")
        cache = Path(env) if env else DEFAULT_CACHE
        if env and not cache.is_absolute():
            raise MissingDependencyError(
                f"BOLTZ_CACHE must be an absolute path, got {env!r}. That is "
                "Boltz's own requirement, and it fails the same way inside the "
                "tool."
            )
        cache.mkdir(parents=True, exist_ok=True)
        if not (cache / CHECKPOINT).is_file():
            # Not an error: predict downloads them. Said once, because several
            # GB arriving silently looks like a hang.
            print(
                f"boltz weights are not in {cache}; the first prediction "
                "downloads them (several GB). Set BOLTZ_CACHE to put them "
                "elsewhere.",
                file=sys.stderr,
                flush=True,
            )
        return cache

    def _collect(self, work: Path, output_dir: Path, sequence: Design) -> Design:
        """The rank-0 prediction, as one Design.

        Boltz writes `boltz_results_<input>/predictions/<record>/<record>_model_N.cif`,
        plus a confidence JSON and the PAE npz beside it. Globbing rather than
        rebuilding that path: the layout has moved between versions.
        """
        cif = first_match(work, "*_model_0.cif", "*.cif")
        if cif is None:
            raise RunFailedError(
                f"boltz exited 0 but wrote no mmCIF for {sequence.id!r} under {work}"
            )
        pdb = cif_to_pdb(cif)
        confidence = first_match(work, "confidence_*_model_0.json", "confidence_*.json")
        scores = parse_confidence(confidence) if confidence else {}
        affinity = first_match(work, "affinity_*.json")
        if affinity is not None:
            scores.update(parse_affinity(affinity))

        meta: dict[str, Any] = {"cif": str(cif)}
        if confidence is not None:
            meta["raw_scores"] = str(confidence)
        if affinity is not None:
            meta["affinity_json"] = str(affinity)
        pae = first_match(work, "pae_*_model_0.npz", "pae_*.npz")
        per_token = first_match(work, "plddt_*_model_0.npz", "plddt_*.npz")
        meta["scores_json"] = str(
            write_scores_json(
                work / f"{pdb.stem}_scores.json",
                # Per-token pLDDT is a fraction here too.
                plddt=(
                    100.0 * np.load(per_token)["plddt"]
                    if per_token is not None
                    else None
                ),
                pae=np.load(pae)["pae"] if pae is not None else None,
                ptm=scores.get("ptm"),
                iptm=scores.get("iptm"),
            )
        )
        return Design(
            id=make_id(pdb, output_dir),
            path=pdb,
            parent_id=sequence.id,
            producer=PRODUCER,
            scores=scores,
            meta=meta,
        )
