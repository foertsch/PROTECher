"""Fixture adapters, one per step, so the whole round loop runs on a laptop.

No GPU, no container runtime, no network. They write real (if tiny) PDB,
FASTA and scores-JSON files, so the scoring module is exercised for real
rather than stubbed.
"""

from __future__ import annotations

import hashlib
import json
from dataclasses import dataclass
from pathlib import Path
from typing import Any, Sequence

from ..containers import make_id
from ..designs import Design
from ._folding import cif_to_pdb, work_dir
from ._folding import write_scores_json as write_normalized_scores

# Aliased on import: this module has its own write_scores_json, which writes
# ColabFold's shape directly for DummyStructures, while DummyFold goes through
# the shared normalizer the way the real co-folding adapters do.

CHAIN_A = 10
"""Residues in the target chain."""
CHAIN_B = 12
"""Residues in the binder chain."""
SEQUENCE = "MKTAYIAKQRQISFVKSHFSRQLEERLGLIEVQ"


def _jitter(token: str) -> float:
    """Deterministic pseudo-displacement in [0, 1), from a label."""
    digest = hashlib.sha256(token.encode()).digest()
    return int.from_bytes(digest[:4], "big") / 2**32


def write_pdb(
    path: Path,
    token: str,
    spread: float = 0.0,
    residues: tuple[int, int] = (CHAIN_A, CHAIN_B),
) -> None:
    """A CA-only two-chain structure. `spread` offsets the binder chain so a
    prediction differs from the backbone it came from and sc-RMSD is non-zero.
    """
    lines: list[str] = []
    serial = 1
    for chain, count in zip(("A", "B"), residues):
        for resi in range(1, count + 1):
            shift = spread * _jitter(f"{token}:{chain}{resi}") if chain == "B" else 0.0
            x = 1.5 * resi + shift
            y = 3.0 if chain == "B" else 0.0
            z = 0.0
            plddt = 50.0 + 40.0 * _jitter(f"{token}:b:{chain}{resi}")
            lines.append(
                "ATOM  "
                + f"{serial:5d}"
                + " "
                + " CA "
                + " "
                + f"{'GLY':>3}"
                + " "
                + chain
                + f"{resi:4d}"
                + "    "
                + f"{x:8.3f}{y:8.3f}{z:8.3f}"
                + f"{1.0:6.2f}{plddt:6.2f}"
                + "          C"
            )
            serial += 1
        lines.append("TER")
    lines.append("END")
    path.write_text("\n".join(lines) + "\n")


def write_scores_json(path: Path, token: str) -> dict[str, float]:
    """A ColabFold-shaped scores file, plus the summary scores it implies."""
    n = CHAIN_A + CHAIN_B
    plddt = [50.0 + 40.0 * _jitter(f"{token}:p:{i}") for i in range(n)]
    # Inter-chain PAE deliberately worse than intra-chain, as in a real run.
    pae = [
        [
            (3.0 if (i < CHAIN_A) == (j < CHAIN_A) else 12.0)
            + 4.0 * _jitter(f"{token}:{i}:{j}")
            for j in range(n)
        ]
        for i in range(n)
    ]
    ptm = 0.5 + 0.4 * _jitter(f"{token}:ptm")
    iptm = 0.3 + 0.5 * _jitter(f"{token}:iptm")
    path.write_text(json.dumps({"plddt": plddt, "pae": pae, "ptm": ptm, "iptm": iptm}))
    return {
        "plddt": sum(plddt) / len(plddt),
        "ptm": ptm,
        "iptm": iptm,
    }


def write_cif(pdb: Path, cif: Path) -> Path:
    """The PDB fixture again as mmCIF, which is what the co-folding tools write.

    Here only so `DummyFold` can hand the real `cif_to_pdb` a real mmCIF: the
    conversion is on the path of every prediction from those three models, and
    stubbing it offline would leave it untested.
    """
    from biotite.structure.io import pdbx
    from biotite.structure.io.pdb import PDBFile

    structure = PDBFile.read(str(pdb)).get_structure(
        model=1, extra_fields=["b_factor", "occupancy"]
    )
    out = pdbx.CIFFile()
    pdbx.set_structure(out, structure)
    out.write(str(cif))
    return cif


def write_raw_scores(path: Path, token: str) -> dict[str, Any]:
    """A confidence file in a co-folding tool's own shape.

    A fractional complex pLDDT and a token-pair PAE matrix, i.e. what Boltz and
    Protenix write and what the adapters have to normalize. Deliberately not
    ColabFold's shape: the normalization is the part worth testing.
    """
    n = CHAIN_A + CHAIN_B
    plddt = [0.5 + 0.49 * _jitter(f"{token}:p:{i}") for i in range(n)]
    pae = [
        [
            (3.0 if (i < CHAIN_A) == (j < CHAIN_A) else 12.0)
            + 4.0 * _jitter(f"{token}:{i}:{j}")
            for j in range(n)
        ]
        for i in range(n)
    ]
    payload = {
        "complex_plddt": sum(plddt) / len(plddt),
        "per_token_plddt": plddt,
        "ptm": 0.5 + 0.4 * _jitter(f"{token}:ptm"),
        "iptm": 0.3 + 0.5 * _jitter(f"{token}:iptm"),
        "pae": pae,
    }
    path.write_text(json.dumps(payload))
    return payload


@dataclass(frozen=True)
class DummyFold:
    """Stands in for one of the co-folding adapters: Chai-1, Boltz-2, Protenix.

    Writes an mmCIF and a tool-shaped confidence file into a directory per
    input, then goes through the same `cif_to_pdb` and `write_scores_json` the
    real adapters use. So the offline suite exercises the conversion and the
    score normalization rather than stubbing them, and `producer` differs per
    instance, which is what a mixed pool needs.
    """

    producer: str

    def run(
        self,
        inputs: Sequence[Design],
        output_dir: Path,
        *,
        extra_args: Sequence[str] | None = None,
        seed: int = 0,
        **kwargs: Any,
    ) -> list[Design]:
        designs: list[Design] = []
        for sequence in inputs:
            work = work_dir(output_dir, sequence)
            token = f"{sequence.id}:{seed}"
            # The real tools write mmCIF only, so the PDB this starts from is
            # not left behind to be mistaken for the prediction.
            staging = work / "staging.pdb"
            write_pdb(staging, token=token, spread=1.2)
            cif = write_cif(staging, work / f"{sequence.id}_model_0.cif")
            staging.unlink()

            raw_file = work / "raw_scores.json"
            raw = write_raw_scores(raw_file, token=token)
            pdb = cif_to_pdb(cif)
            plddt = [100.0 * value for value in raw["per_token_plddt"]]
            scores_json = write_normalized_scores(
                work / f"{pdb.stem}_scores.json",
                plddt=plddt,
                pae=raw["pae"],
                ptm=raw["ptm"],
                iptm=raw["iptm"],
            )
            designs.append(
                Design(
                    id=make_id(pdb, output_dir),
                    path=pdb,
                    parent_id=sequence.id,
                    producer=self.producer,
                    scores={
                        "plddt": 100.0 * raw["complex_plddt"],
                        "ptm": raw["ptm"],
                        "iptm": raw["iptm"],
                    },
                    meta={
                        "cif": str(cif),
                        "raw_scores": str(raw_file),
                        "scores_json": str(scores_json),
                        "seed": seed,
                    },
                )
            )
        return designs


@dataclass(frozen=True)
class DummyBackbones:
    """Stands in for RFdiffusion, for both fresh generation and partial diffusion."""

    def run(
        self,
        inputs: Sequence[Design],
        output_dir: Path,
        *,
        extra_args: Sequence[str] | None = None,
        num_designs: int = 1,
        seed: int = 0,
        **kwargs: Any,
    ) -> list[Design]:
        parents: Sequence[Design | None] = inputs or [None]
        designs: list[Design] = []
        # Numbered from a flat counter, like RFdiffusion's fixed output_prefix:
        # names stay bounded instead of compounding the parent's name every
        # round. Ancestry lives in parent_id, which is the point of having it.
        counter = 0
        for parent in parents:
            for _ in range(int(num_designs)):
                path = output_dir / f"design_{counter}.pdb"
                counter += 1
                write_pdb(path, token=f"{path.stem}:{seed}")
                designs.append(
                    Design(
                        id=make_id(path, output_dir),
                        path=path,
                        parent_id=None if parent is None else parent.id,
                        producer="dummy_backbones",
                        meta={"seed": seed},
                    )
                )
        return designs


@dataclass(frozen=True)
class DummySequences:
    """Stands in for the MPNN step: one design per sampled sequence.

    Mirrors what `ligandmpnn` actually writes — the same score names and the
    same model_type/temperature in meta — so offline tests exercise the shape
    the real adapter produces.
    """

    def run(
        self,
        inputs: Sequence[Design],
        output_dir: Path,
        *,
        extra_args: Sequence[str] | None = None,
        num_seqs: int = 1,
        temperature: float = 0.1,
        seed: int = 0,
        **kwargs: Any,
    ) -> list[Design]:
        designs: list[Design] = []
        for backbone in inputs:
            for i in range(int(num_seqs)):
                # Named after the input structure, as the real tool does.
                path = output_dir / f"{backbone.path.stem}_sample{i}.fa"
                design_id = make_id(path, output_dir)
                # The record name carries the run key, so two sequence arms over
                # the same backbones cannot collide in one fold run.
                path.write_text(f">{design_id}\n{SEQUENCE}\n")
                token = f"{path.stem}:{seed}"
                designs.append(
                    Design(
                        id=design_id,
                        path=path,
                        parent_id=backbone.id,
                        producer="dummy_sequences",
                        scores={
                            "overall_confidence": 0.5 + 0.5 * _jitter(token),
                            "seq_rec": 0.3 + 0.4 * _jitter(f"{token}:r"),
                        },
                        meta={
                            "seed": seed,
                            "model_type": "dummy_sequences",
                            "temperature": temperature,
                        },
                    )
                )
        return designs


@dataclass(frozen=True)
class DummyStructures:
    """Stands in for ColabFold: a predicted PDB plus a scores JSON."""

    def run(
        self,
        inputs: Sequence[Design],
        output_dir: Path,
        *,
        extra_args: Sequence[str] | None = None,
        seed: int = 0,
        **kwargs: Any,
    ) -> list[Design]:
        designs: list[Design] = []
        for sequence in inputs:
            token = f"{sequence.id}:{seed}"
            # ColabFold names outputs after the FASTA record, which is the
            # design id; two sequence arms can therefore share one fold run.
            pdb = output_dir / f"{sequence.id}_unrelaxed_rank_001.pdb"
            scores_json = output_dir / f"{sequence.id}_scores_rank_001.json"
            write_pdb(pdb, token=token, spread=1.2)
            scores = write_scores_json(scores_json, token=token)
            designs.append(
                Design(
                    id=make_id(pdb, output_dir),
                    path=pdb,
                    parent_id=sequence.id,
                    producer="dummy_structures",
                    scores=scores,
                    meta={"scores_json": str(scores_json), "seed": seed},
                )
            )
        return designs
