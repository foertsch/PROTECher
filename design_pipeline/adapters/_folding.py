"""What the co-folding adapters (Chai-1, Boltz-2, Protenix) share.

Three tools, three input formats, three score files, one `Design`. The
pipeline's contract is narrow and the adapters bend to it rather than the other
way round:

* `score_designs` superimposes CA atoms on the backbone a prediction came from,
  so `path` is always a PDB. These models write mmCIF, so `cif_to_pdb` converts
  and the original stays in `meta["cif"]`.
* `select` sorts one score name across a pool that a model list leaves mixed,
  so every adapter reports `plddt` (0-100, the scale ColabFold uses), `ptm` and
  `iptm`, whatever the tool calls them. The tool's own file is kept in
  `meta["raw_scores"]`, and anything beyond the shared three goes into `scores`
  under a producer-prefixed name.
* `pae_interaction` reads the PAE matrix back out of `meta["scores_json"]`,
  which is ColabFold's format, so `write_scores_json` writes that shape from
  whatever numbers the tool gives us.

All three need Linux and a CUDA GPU and download weights on first use, so
`require_gpu` runs before anything is fetched. Each tool keeps its weights
where it keeps them; there is no cache here to share.
"""

from __future__ import annotations

import json
import shutil
import subprocess
import sys
from pathlib import Path
from typing import Any, Iterable, Sequence

from ..designs import Design
from ..errors import MissingGPUError

BFLOAT16_CAPABILITY = 8.0
"""Compute capability where bfloat16 arrives (Ampere). Below it, Chai-1 is out."""


def gpus() -> list[tuple[str, float | None]]:
    """(name, compute capability) for each visible GPU, from nvidia-smi.

    nvidia-smi rather than `torch.cuda.is_available()`: torch belongs to the
    model's environment, and importing it to ask one question would cost
    seconds on every call and pull a GPU stack into a process that otherwise
    only shells out. A driver too old for the compute_cap query still lists its
    cards, and an unknown capability is reported as None rather than guessed.
    """
    exe = shutil.which("nvidia-smi")
    if exe is None:
        return []
    listed = _nvidia_smi(exe, "--query-gpu=name,compute_cap", "--format=csv,noheader")
    if listed is None:
        return [(name, None) for name in _nvidia_smi(exe, "-L") or ()]
    found: list[tuple[str, float | None]] = []
    for line in listed:
        name, _, capability = line.partition(",")
        try:
            found.append((name.strip(), float(capability)))
        except ValueError:
            found.append((name.strip(), None))
    return found


def _nvidia_smi(exe: str, *args: str) -> list[str] | None:
    try:
        done = subprocess.run(
            [exe, *args], capture_output=True, text=True, timeout=60, check=False
        )
    except OSError:
        return None
    if done.returncode != 0:
        return None
    return [line for line in done.stdout.splitlines() if line.strip()]


def require_gpu(model: str, extra: str, *, bfloat16: bool = False) -> None:
    """Fail before the weights download, not deep inside the model.

    `extra` names the optional extra, and every message points at the dummy
    counterpart: wanting to run the loop on a laptop is the common reason to
    meet this error.
    """
    offline = f"To run the loop without a folding model, use model_name='dummy_{extra}'."
    if sys.platform != "linux":
        raise MissingGPUError(
            f"{model} runs on Linux with a CUDA GPU; this is {sys.platform}. "
            f"Run it on the cluster or another Linux host with an NVIDIA card. "
            f"{offline}"
        )
    found = gpus()
    if not found:
        raise MissingGPUError(
            f"{model} needs a CUDA GPU and nvidia-smi reports none. Check that "
            f"this is a GPU node (`nvidia-smi -L`) and that the driver is "
            f"loaded. {offline}"
        )
    if bfloat16:
        known = [capability for _, capability in found if capability is not None]
        if known and max(known) < BFLOAT16_CAPABILITY:
            cards = ", ".join(
                f"{name} (compute {capability})" for name, capability in found
            )
            raise MissingGPUError(
                f"{model} runs in bfloat16, which needs compute capability "
                f"{BFLOAT16_CAPABILITY} or newer (Ampere); this host has "
                f"{cards}. Use model_name='boltz' or 'colabfold' on older "
                f"cards. {offline}"
            )


def cif_to_pdb(cif: Path, pdb: Path | None = None) -> Path:
    """Convert an mmCIF prediction to PDB, keeping the B-factor column.

    Scoring reads CA atoms and their B-factors, which is where these models
    write per-atom pLDDT, so the conversion has to carry that annotation
    across: dropping it would leave `mean_plddt` with nothing to fall back on.
    """
    from biotite.structure.io import pdbx
    from biotite.structure.io.pdb import PDBFile

    target = pdb if pdb is not None else cif.with_suffix(".pdb")
    structure = pdbx.get_structure(
        pdbx.CIFFile.read(str(cif)), model=1, extra_fields=["b_factor", "occupancy"]
    )
    out = PDBFile()
    out.set_structure(structure)
    out.write(str(target))
    return target


def fasta_chains(fasta: Path) -> tuple[str, list[str]]:
    """The record name and its chains.

    The sequence step writes one record per design and separates chains with
    ':' - LigandMPNN's convention, which ColabFold also reads. These three
    tools want the chains listed separately instead, so this is where the one
    record becomes a list.
    """
    records = [r for r in fasta.read_text().split(">") if r.strip()]
    if len(records) != 1:
        raise ValueError(
            f"{fasta} holds {len(records)} FASTA records; the folding adapters "
            "expect one design per file, which is what design_sequences writes. "
            "Split a multi-record FASTA into one Design per record first."
        )
    header, _, body = records[0].partition("\n")
    chains = [chain.strip() for chain in body.replace("\n", "").split(":")]
    chains = [chain for chain in chains if chain]
    if not chains:
        raise ValueError(f"{fasta} has a header but no sequence")
    return header.split()[0], chains


def write_scores_json(
    path: Path,
    *,
    plddt: Sequence[float] | None = None,
    pae: Iterable[Iterable[float]] | None = None,
    **scalars: float | None,
) -> Path:
    """Write another tool's confidence numbers in ColabFold's scores-JSON shape.

    `scoring.pae_interaction` reads this file and nothing else does, so what
    has to match is the `pae` key and the 0-100 pLDDT scale, not the rest of
    ColabFold's payload.
    """
    payload: dict[str, Any] = {
        key: float(value) for key, value in scalars.items() if value is not None
    }
    if plddt is not None:
        payload["plddt"] = [float(v) for v in plddt]
    if pae is not None:
        payload["pae"] = [[float(v) for v in row] for row in pae]
    path.write_text(json.dumps(payload))
    return path


def work_dir(output_dir: Path, design: Design) -> Path:
    """An empty directory per input, inside the run directory.

    One directory per input keeps outputs apart - all three tools name files
    after the input record - and satisfies Chai-1, which refuses to write into
    a directory that is not empty. A leftover from an interrupted attempt is
    removed rather than worked around: the run directory is content-addressed,
    so a half-finished one holds nothing worth keeping.
    """
    work = output_dir / design.id
    if work.exists():
        shutil.rmtree(work)
    work.mkdir(parents=True)
    return work


def first_match(work: Path, *patterns: str) -> Path | None:
    """The first file matching any pattern, patterns tried in order.

    These tools bury outputs under directories named after the input, the
    model and the seed, and the layout moves between versions, so the adapters
    glob for what they need instead of rebuilding the path.
    """
    for pattern in patterns:
        found = sorted(work.rglob(pattern))
        if found:
            return found[0]
    return None
