"""Reuse one set of Boltz-2 alignments (MSAs) for every variant of a parent.

The public alignment server takes about 2 minutes per sequence and is shared and rate limited.
Variants of one parent differ in a handful of positions, so the parent's alignments serve them all.
Only the first row of an alignment, the query itself, changes: it is swapped for the variant's own
sequence. In a test on the trastuzumab scFv this stayed within seed noise of calling the server for
each variant, while folding with no alignment lowered ipTM by 0.08 to 0.11.

    cache = MsaCache("runs/msa_parent")
    if not cache.has_parent():
        build_parent(cache, [scfv, antigen], work_dir="runs/msa_parent_fold")   # one server call
    with use_cached_msas(cache):
        predict_structure(variants, model_name="boltz", use_msa_server=False)   # no server calls
"""

from __future__ import annotations

import csv
import hashlib
import json
import shutil
from contextlib import contextmanager
from datetime import datetime, timezone
from pathlib import Path
from typing import Iterator, Sequence


class CacheError(RuntimeError):
    """The cache cannot serve what was asked, said in words rather than failing inside Boltz."""


def _read(path: Path) -> list[list[str]]:
    with open(path, newline="") as fh:
        return list(csv.reader(fh))


def _write(path: Path, rows: list[list[str]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with open(path, "w", newline="") as fh:
        csv.writer(fh).writerows(rows)


def _sequence_column(rows: list[list[str]]) -> int:
    try:
        return rows[0].index("sequence")
    except ValueError:
        raise CacheError(f"alignment has no 'sequence' column, header is {rows[0]}") from None


class MsaCache:
    """A directory holding the parent's alignments, one csv per chain, and the variants derived from them."""

    def __init__(self, directory: Path | str, max_mismatch_fraction: float = 0.1) -> None:
        self.dir = Path(directory)
        self.dir.mkdir(parents=True, exist_ok=True)
        self.max_mismatch_fraction = max_mismatch_fraction

    @property
    def _meta(self) -> Path:
        return self.dir / "parent.json"

    def has_parent(self) -> bool:
        return self._meta.is_file()

    @property
    def parent_chains(self) -> list[str]:
        if not self.has_parent():
            raise CacheError(f"{self.dir} holds no parent alignments yet: run build_parent or import_boltz_msas")
        return json.loads(self._meta.read_text())["chains"]

    def import_boltz_msas(self, msa_dir: Path | str, chains: Sequence[str]) -> None:
        """Keep the alignments Boltz saved for the parent, found as `<name>_<chain index>.csv`."""
        msa_dir = Path(msa_dir)
        for index, chain in enumerate(chains):
            found = sorted(msa_dir.glob(f"*_{index}.csv"))
            if len(found) != 1:
                raise CacheError(f"expected one *_{index}.csv in {msa_dir}, found {[p.name for p in found]}")
            rows = _read(found[0])
            if len(rows) < 2 or rows[1][_sequence_column(rows)] != chain:
                raise CacheError(f"{found[0].name}: the first alignment row is not chain {index + 1} of the parent")
            shutil.copy(found[0], self.dir / f"parent_{index}.csv")
        self._meta.write_text(
            json.dumps(
                {"chains": list(chains), "created": datetime.now(timezone.utc).strftime("%Y-%m-%d %H:%M UTC")},
                indent=1,
            )
        )

    def alignment_for(self, index: int, chain: str) -> Path | None:
        """The alignment file for one chain of a variant, or None when the chain is not close to the parent's."""
        chains = self.parent_chains
        if index >= len(chains):
            return None
        parent_path = self.dir / f"parent_{index}.csv"
        parent = chains[index]
        if chain == parent:
            return parent_path
        if len(chain) != len(parent):
            return None
        mismatches = sum(a != b for a, b in zip(chain, parent))
        if mismatches > self.max_mismatch_fraction * len(parent):
            return None
        out = self.dir / "variants" / f"{hashlib.sha1(chain.encode()).hexdigest()[:12]}_{index}.csv"
        if not out.is_file():
            rows = _read(parent_path)
            rows[1][_sequence_column(rows)] = chain
            _write(out, rows)
        return out

    def for_chains(self, chains: Sequence[str]) -> list[Path | None]:
        return [self.alignment_for(i, c) for i, c in enumerate(chains)]

    def info(self) -> dict:
        chains = self.parent_chains
        rows = [len(_read(self.dir / f"parent_{i}.csv")) - 1 for i in range(len(chains))]
        return {
            "directory": str(self.dir),
            "chain_lengths": [len(c) for c in chains],
            "alignment_rows": rows,
            "variant_files": len(list((self.dir / "variants").glob("*.csv"))) if (self.dir / "variants").is_dir() else 0,
        }


def build_parent(cache: MsaCache, chains: Sequence[str], work_dir: Path | str, seed: int = 0):
    """Fold the parent once with the public server and keep its alignments. Needs a GPU, like any Boltz run."""
    from design_pipeline import Design, predict_structure

    work_dir = Path(work_dir)
    work_dir.mkdir(parents=True, exist_ok=True)
    fasta = work_dir / "parent.fa"
    fasta.write_text(">parent\n" + ":".join(chains) + "\n")
    parent = Design(id="parent", path=fasta, parent_id=None, producer="manual")
    folded = predict_structure(
        [parent], output_dir=work_dir / "predictions", model_name="boltz", seed=seed, use_msa_server=True
    )
    msa_dirs = [p for p in (work_dir / "predictions").rglob("msa") if p.parent.name.startswith("boltz_results_")]
    if len(msa_dirs) != 1:
        raise CacheError(f"expected one Boltz msa folder under {work_dir}, found {len(msa_dirs)}")
    cache.import_boltz_msas(msa_dirs[0], chains)
    return folded


@contextmanager
def use_cached_msas(cache: MsaCache, fallback: str = "error") -> Iterator[MsaCache]:
    """Make the Boltz adapter point every chain at the cache, for folds run with use_msa_server=False.

    A chain too far from the parent's raises CacheError, so no fold silently runs without an alignment.
    With fallback="server" such a chain is written without a path and Boltz asks the public server for
    it, which needs use_msa_server=True. Kira's package is patched in memory for the duration only.
    """
    if fallback not in ("error", "server"):
        raise ValueError("fallback must be 'error' or 'server'")
    from design_pipeline.adapters import _boltz
    from design_pipeline.adapters._folding import fasta_chains

    original = _boltz.write_boltz_fasta

    def write(path: Path, sequence) -> Path:
        _, chains = fasta_chains(sequence.path)
        records = []
        for k, (chain, msa) in enumerate(zip(chains, cache.for_chains(chains))):
            if msa is None and fallback == "error":
                raise CacheError(
                    f"design {sequence.id}: chain {k + 1} is not within {cache.max_mismatch_fraction:.0%} "
                    "substitutions of the parent's chain of the same length. Use fallback='server' (with "
                    "use_msa_server=True) or build a cache for this antibody."
                )
            header = f">{chr(ord('A') + k)}|protein" + ("" if msa is None else f"|{msa.resolve()}")
            records.append(f"{header}\n{chain}\n")
        path.write_text("".join(records))
        return path

    _boltz.write_boltz_fasta = write
    try:
        yield cache
    finally:
        _boltz.write_boltz_fasta = original
