"""ColabFold, via the `colabfold_batch` entry point.

Subprocess rather than the Python API: it is the interface ColabFold documents,
it keeps jax out of this process, and its stdout is the progress report we want
streamed.
"""

from __future__ import annotations

import json
import shutil
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any, Mapping, Sequence

from ..containers import format_kwargs, make_id, stream
from ..designs import Design
from ..errors import MissingDependencyError, RunFailedError

PUBLIC_MSA_SERVER = "https://api.colabfold.com"


def parse_scores(path: Path) -> dict[str, float]:
    """ColabFold's per-prediction scores JSON."""
    try:
        payload = json.loads(path.read_text())
    except (OSError, json.JSONDecodeError):
        return {}
    scores: dict[str, float] = {}
    plddt = payload.get("plddt")
    if plddt:
        scores["plddt"] = sum(plddt) / len(plddt)
    for key in ("ptm", "iptm"):
        if payload.get(key) is not None:
            scores[key] = float(payload[key])
    return scores


def find_scores_json(pdb: Path) -> Path | None:
    """The scores file beside a prediction. Names vary between ColabFold
    versions, so match on the rank rather than on the full name."""
    rank = next((p for p in pdb.stem.split("_") if p.isdigit()), None)
    candidates = sorted(pdb.parent.glob("*scores*.json"))
    if rank is not None:
        ranked = [p for p in candidates if f"rank_{rank}" in p.name]
        if ranked:
            return ranked[0]
    return candidates[0] if candidates else None


@dataclass(frozen=True)
class ColabFoldAdapter:
    executable: str = "colabfold_batch"
    aliases: Mapping[str, str] = field(
        default_factory=lambda: {"recycles": "num-recycle", "seed": "random-seed"}
    )

    def run(
        self,
        inputs: Sequence[Design],
        output_dir: Path,
        *,
        extra_args: Sequence[str] | None = None,
        host_url: str = PUBLIC_MSA_SERVER,
        **kwargs: Any,
    ) -> list[Design]:
        exe = shutil.which(self.executable)
        if exe is None:
            raise MissingDependencyError(
                f"{self.executable!r} is not on PATH. Install the optional "
                "extra with `uv pip install 'design-pipeline[colabfold]'` (or "
                "`pip install colabfold`), or use model_name="
                "'dummy_structures' to run the loop without a folding model."
            )
        args = format_kwargs(kwargs, "flag", self.aliases, dash_flags=True)
        # Explicit, so a local MSA server can be pointed at later.
        args += ["--host-url", host_url]

        designs: list[Design] = []
        for fasta in inputs:
            # One call per FASTA: colabfold_batch also takes a directory, but
            # our inputs can come from more than one.
            before = set(output_dir.rglob("*.pdb"))
            stream(
                [exe, str(fasta.path), str(output_dir), *args, *(extra_args or ())]
            )
            produced = sorted(p for p in output_dir.rglob("*.pdb") if p not in before)
            # Only the best-ranked model becomes a design; the rest stay on
            # disk for inspection.
            ranked = [p for p in produced if "rank_001" in p.name] or produced[:1]
            for pdb in ranked:
                scores_json = find_scores_json(pdb)
                designs.append(
                    Design(
                        id=make_id(pdb, output_dir),
                        path=pdb,
                        parent_id=fasta.id,
                        producer="colabfold",
                        scores=parse_scores(scores_json) if scores_json else {},
                        # Scoring reads the PAE matrix back out of this.
                        meta={"scores_json": str(scores_json)} if scores_json else {},
                    )
                )
        if not designs:
            raise RunFailedError(
                f"{self.executable} wrote no PDB files under {output_dir}"
            )
        return designs
