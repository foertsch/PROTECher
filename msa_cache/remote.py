"""Run fold_batch on a machine that does not hold the data, such as a Modal container.

The side that has the files stages one folder, the remote side folds it, and the results come back with
their absolute paths rewritten to local ones. Nothing here knows about Modal, so it is tested with a
stand-in `boltz`; cloud/modal_fold.py is the thin layer that moves the folders.

    stage("designs/", "data/msa_candidates/trastuzumab_her2", "stage/run1")   # locally
    run_staged("stage/run1", seed=0)                                          # remotely: folds into stage/run1/out
    localize_manifest("out/", remote_out="/work/run1/out")                    # locally, after the download

Only the standard library is imported on the staging side, so it runs in a bare Modal client environment.
"""

from __future__ import annotations

import json
import shutil
from pathlib import Path
from typing import Sequence

from .batch import fold_batch
from .cache import MAX_MISMATCH_FRACTION, CacheError, MsaCache, _best_offset

ORDER = "ORDER"
"""File in designs/ listing the ids in the order they are folded. A design's random draw depends on its place in the batch."""


def read_chains(fasta: Path) -> list[str]:
    """The chains of a one-record FASTA whose chains are separated by ':', the pipeline's sequence files."""
    records = [r for r in Path(fasta).read_text().split(">") if r.strip()]
    if len(records) != 1:
        raise CacheError(f"{fasta} holds {len(records)} FASTA records, expected one design per file")
    chains = [c.strip() for c in records[0].partition("\n")[2].replace("\n", "").split(":") if c.strip()]
    if not chains:
        raise CacheError(f"{fasta} has no sequence")
    return chains


def stage(designs_dir: Path | str, cache_dir: Path | str, dest: Path | str, order: Sequence[str] | None = None) -> list[str]:
    """Collect the design FASTAs and the parent's alignments into dest/designs and dest/cache, ready to upload.

    Every design is checked against the cache here, so one that is too far from the parent fails now with
    its name instead of after the upload. The cache's `variants/` folder is left behind: it is rebuilt remotely
    from the parent's alignment and would only add megabytes. Returns the ids in folding order.
    """
    designs_dir, cache_dir, dest = Path(designs_dir), Path(cache_dir), Path(dest)
    if not (cache_dir / "parent.json").is_file():
        raise CacheError(f"{cache_dir} holds no parent alignments (no parent.json)")
    found = {p.stem: p for p in sorted(designs_dir.glob("*.fa"))}
    if not found:
        raise CacheError(f"no *.fa files in {designs_dir}")
    ids = list(order) if order is not None else list(found)
    if sorted(ids) != sorted(found):
        raise CacheError(f"order {sorted(ids)} does not match the designs {sorted(found)}")

    parent = json.loads((cache_dir / "parent.json").read_text())["chains"]
    for design_id, path in found.items():
        for k, chain in enumerate(read_chains(path)):
            if k >= len(parent) or (chain != parent[k] and _best_offset(chain, parent[k], MAX_MISMATCH_FRACTION) is None):
                raise CacheError(f"design {design_id}: chain {k + 1} does not fit a stretch of the parent's chain in {cache_dir}")

    if dest.exists():
        shutil.rmtree(dest)
    (dest / "designs").mkdir(parents=True)
    (dest / "cache").mkdir()
    for design_id in ids:
        shutil.copy(found[design_id], dest / "designs" / f"{design_id}.fa")
    (dest / "designs" / ORDER).write_text("\n".join(ids) + "\n")
    for name in ["parent.json", *sorted(p.name for p in cache_dir.glob("parent_*.csv"))]:
        shutil.copy(cache_dir / name, dest / "cache" / name)
    return ids


def run_staged(workdir: Path | str, seed: int = 0) -> list:
    """Fold everything staged in workdir with fold_batch, in the staged order, into workdir/out."""
    from design_pipeline import Design

    workdir = Path(workdir)
    ids = (workdir / "designs" / ORDER).read_text().split()
    designs = [Design(id=i, path=workdir / "designs" / f"{i}.fa", parent_id=None, producer="manual") for i in ids]
    return fold_batch(MsaCache(workdir / "cache"), designs, workdir / "out", seed=seed)


def localize_manifest(out_dir: Path | str, remote_out: str) -> int:
    """Point the paths in a downloaded DONE manifest's meta at out_dir instead of the remote folder.

    The folding adapter records absolute paths (the mmCIF, the scores JSON). scoring.pae_interaction
    reads the scores JSON through them, so they have to exist locally. Returns how many were rewritten.
    """
    out_dir = Path(out_dir).resolve()
    marker = out_dir / "DONE"
    records = json.loads(marker.read_text())
    prefix = remote_out.rstrip("/")
    changed = 0
    for record in records:
        for key, value in record.get("meta", {}).items():
            if isinstance(value, str) and (value == prefix or value.startswith(prefix + "/")):
                record["meta"][key] = str(out_dir) + value[len(prefix) :]
                changed += 1
    marker.write_text(json.dumps(records, indent=1))
    return changed
