"""Protenix through the `protenix predict` command line.

Two things make this adapter different from the other two.

**Its input is JSON, not FASTA.** One job per design:
`[{"name": <design id>, "sequences": [{"proteinChain": {"sequence": ..., "count": 1}}]}]`,
one entry per chain. `write_protenix_json` is that conversion.

**The checkpoint matters.** `checkpoint` names it and the default is the
tool's own, `protenix_base_default_v1.0.0`. The Mini checkpoints are
substantially faster at comparable quality and are the better screening
choice: `protenix_mini_default_v0.5.0` (MSA), and `protenix_mini_esm_v0.5.0` /
`protenix_mini_ism_v0.5.0`, which run from a language-model embedding and need
no MSA at all. Pass `use_default_params=True` with a Mini checkpoint to get
the cycle and step counts it was tuned for (4 and 5, against 10 and 200).

**Licensing depends on the checkpoint and is not settled here.** The pinned
version's own LICENSE file and README say Apache 2.0 for code and model
parameters; the project's GitHub page says the same, and separately that the
`protenix-v2` checkpoint's weights are proprietary and not released under any
open-source licence. Picking a checkpoint is therefore a licence decision as
well as a speed one. See the model table in README.md; nothing here asserts
that any of these weights are cleared for commercial use.

Weights and the CCD data download on first use under `$PROTENIX_ROOT_DIR` (the
tool's own variable, default the home directory).
"""

from __future__ import annotations

import json
import shutil
from dataclasses import dataclass
from pathlib import Path
from typing import Any, Sequence

from ..containers import make_id, stream
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

PRODUCER = "protenix"
BASE_CHECKPOINT = "protenix_base_default_v1.0.0"
"""The tool's default: full size, MSA, 10 cycles and 200 diffusion steps."""
MINI_CHECKPOINT = "protenix_mini_default_v0.5.0"
"""Smaller and much faster at comparable quality; 4 cycles, 5 steps."""


def write_protenix_json(path: Path, sequence: Design) -> Path:
    """The one-job input file Protenix reads, from one sequence design."""
    name, chains = fasta_chains(sequence.path)
    job = {
        "name": name,
        "sequences": [
            {"proteinChain": {"sequence": chain, "count": 1}} for chain in chains
        ],
    }
    path.write_text(json.dumps([job], indent=1))
    return path


def parse_confidence(path: Path) -> dict[str, float]:
    """Protenix's summary confidence file, under the keys ColabFold uses.

    Its `plddt` is already 0-100. Everything past the shared three keeps a
    `protenix_` prefix.
    """
    try:
        payload = json.loads(path.read_text())
    except (OSError, json.JSONDecodeError):
        return {}
    scores: dict[str, float] = {}
    for key in ("plddt", "ptm", "iptm"):
        value = payload.get(key)
        if isinstance(value, (int, float)):
            scores[key] = float(value)
    for key in ("ranking_score", "has_clash", "disorder"):
        value = payload.get(key)
        if isinstance(value, (int, float)):
            scores[f"protenix_{key}"] = float(value)
    return scores


def parse_pae(path: Path) -> list[list[float]] | None:
    """The token-pair PAE matrix out of a full-data file, when one was written."""
    try:
        payload = json.loads(path.read_text())
    except (OSError, json.JSONDecodeError):
        return None
    pae = payload.get("token_pair_pae") or payload.get("pae")
    return pae if isinstance(pae, list) else None


def _args(options: dict[str, Any]) -> list[str]:
    """Protenix options are `type=bool`, not switches.

    `format_kwargs` renders True as a bare flag, which click rejects here, so
    booleans are passed as values instead and None drops out.
    """
    args: list[str] = []
    for key, value in options.items():
        if value is None:
            continue
        args += [f"--{key}", str(value)]
    return args


@dataclass(frozen=True)
class ProtenixAdapter:
    """Satisfies the ModelAdapter protocol."""

    executable: str = "protenix"
    extra: str = "protenix"
    env_var: str | None = "DESIGN_PIPELINE_PROTENIX"
    """The same variable `routes.py` probes. Protenix is not provisioned any
    more - it served de novo design - but the adapter stays, and an adapter that
    stays has to be reachable the same way the others are."""

    def run(
        self,
        inputs: Sequence[Design],
        output_dir: Path,
        *,
        extra_args: Sequence[str] | None = None,
        checkpoint: str = BASE_CHECKPOINT,
        use_msa: bool = True,
        seed: int = 0,
        **kwargs: Any,
    ) -> list[Design]:
        require_gpu("protenix", self.extra)
        exe = self._resolve()
        args = _args(
            {
                "use_msa": use_msa,
                # The token-pair PAE lives in the full-data file, and
                # pae_interaction needs it.
                "need_atom_confidence": True,
                **kwargs,
            }
        )

        designs: list[Design] = []
        for sequence in inputs:
            work = work_dir(output_dir, sequence)
            job = write_protenix_json(work / "input.json", sequence)
            stream(
                [
                    exe,
                    "predict",
                    "--input",
                    str(job),
                    "--out_dir",
                    str(work),
                    "--seeds",
                    str(seed),
                    "--model_name",
                    checkpoint,
                    *args,
                    *(extra_args or ()),
                ]
            )
            designs.append(self._collect(work, output_dir, sequence, checkpoint))

        if not designs:
            raise RunFailedError(f"protenix returned no predictions under {output_dir}")
        return designs

    def _resolve(self) -> str:
        from ..routes import configured_path

        configured = configured_path(self.env_var) if self.env_var else None
        if configured is not None:
            return str(configured)
        found = shutil.which(self.executable)
        if found:
            return found
        raise MissingDependencyError(
            f"{self.executable!r} is not on PATH and "
            f"${self.env_var} is not set. Install the optional extra "
            "with `uv pip install -e '.[protenix]'` (or `pip3 install "
            "protenix`) into its own environment and point "
            f"${self.env_var} at that environment's `protenix`, or use "
            "model_name='dummy_protenix' to run the loop without a folding "
            "model."
        )

    def _collect(
        self, work: Path, output_dir: Path, sequence: Design, checkpoint: str
    ) -> Design:
        """The rank-0 sample, as one Design.

        Protenix writes
        `<dataset>/<name>/seed_<seed>/predictions/<name>_sample_<rank>.cif`
        with the summary confidence beside it. Globbed rather than rebuilt:
        the directory names come from the input and the seed.
        """
        cif = first_match(work, "*_sample_0.cif", "*.cif")
        if cif is None:
            raise RunFailedError(
                f"protenix exited 0 but wrote no mmCIF for {sequence.id!r} "
                f"under {work}"
            )
        pdb = cif_to_pdb(cif)
        confidence = first_match(
            work, "*_summary_confidence_sample_0.json", "*summary_confidence*.json"
        )
        scores = parse_confidence(confidence) if confidence else {}
        full = first_match(work, "*_full_data_sample_0.json", "*full_data*.json")

        meta: dict[str, Any] = {"cif": str(cif), "checkpoint": checkpoint}
        if confidence is not None:
            meta["raw_scores"] = str(confidence)
        if full is not None:
            meta["full_data_json"] = str(full)
        meta["scores_json"] = str(
            write_scores_json(
                work / f"{pdb.stem}_scores.json",
                pae=parse_pae(full) if full is not None else None,
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
