"""MPNN sequence design through the pip-installable `ligandmpnn` package.

One codebase serves ProteinMPNN, LigandMPNN and SolubleMPNN via `--model_type`,
the model parameters ship in the wheel, and it is MIT licensed — so this step
needs no container and no weight download.

Subprocess rather than importing the package, for two reasons: its stdout is
the progress report we want streamed, and its pins (torch==2.2.1,
numpy==1.23.5, Python <3.12) conflict with ours (biotite needs numpy>=2.1), so
in practice it lives in its own environment. Point `executable` at that
environment's `ligandmpnn` and nothing else changes.

The tool writes a multi-record FASTA per input under `seqs/`. We split it into
one file per sampled sequence, so 100 backbones x 8 sequences is 800 Designs
with 800 paths and ColabFold can be handed them one at a time.
"""

from __future__ import annotations

import shutil
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any, Mapping, Sequence

from ..containers import format_kwargs, make_id, stream
from ..designs import Design
from ..errors import MissingDependencyError, MissingLigandError, RunFailedError

MODEL_TYPES = ("protein_mpnn", "ligand_mpnn", "soluble_mpnn")

# Residues that are het atoms but not a ligand worth conditioning on.
SOLVENT = {"HOH", "WAT", "DOD"}


def header_scores(header: str) -> dict[str, float]:
    """Numeric `key=value` fields from a FASTA header.

    LigandMPNN writes `overall_confidence`, `ligand_confidence` and `seq_rec`
    on sampled records. Non-numeric fields (`use_ligand_context=True`,
    `model_path=...`) are skipped rather than coerced.
    """
    scores: dict[str, float] = {}
    for part in header.split(","):
        key, sep, value = part.partition("=")
        if not sep:
            continue
        try:
            scores[key.strip()] = float(value)
        except ValueError:
            continue
    return scores


def het_residues(pdb: Path) -> set[str]:
    """Non-solvent het residue names in a PDB."""
    names: set[str] = set()
    for line in pdb.read_text(errors="ignore").splitlines():
        if line.startswith("HETATM"):
            name = line[17:20].strip()
            if name not in SOLVENT:
                names.add(name)
    return names


@dataclass(frozen=True)
class MPNNAdapter:
    """Satisfies the ModelAdapter protocol. One class, one instance per model_type."""

    model_type: str = "protein_mpnn"
    executable: str = "ligandmpnn"
    """Console script from the `ligandmpnn` wheel. Override with an absolute
    path to run it out of a separate environment."""
    env_var: str | None = "DESIGN_PIPELINE_MPNN"
    """The same variable `routes.py` probes, and the one all three model types
    share - they are one install. Not provisioned any more (de novo design), but
    reached the same way as everything else in Env B while the adapter exists."""
    aliases: Mapping[str, str] = field(default_factory=dict)
    """LigandMPNN's own flags are underscored, so this is mostly empty."""

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
        exe = self._resolve()
        args = self._args(num_seqs, temperature, kwargs)

        designs: list[Design] = []
        seq_dir = output_dir / "seqs"
        for backbone in inputs:
            if self.model_type == "ligand_mpnn":
                self._require_ligand(backbone)
            before = set(seq_dir.glob("*.fa")) | set(seq_dir.glob("*.fasta"))
            stream(
                [
                    exe,
                    "--model_type",
                    self.model_type,
                    "--seed",
                    str(seed),
                    "--pdb_path",
                    str(backbone.path),
                    "--out_folder",
                    str(output_dir),
                    *args,
                    *(extra_args or ()),
                ]
            )
            produced = sorted(
                p
                for p in list(seq_dir.glob("*.fa")) + list(seq_dir.glob("*.fasta"))
                if p not in before
            )
            for source in produced:
                designs += self._split(source, output_dir, backbone, temperature)

        if not designs:
            raise RunFailedError(
                f"{self.executable} exited 0 but wrote no sequences under {seq_dir}"
            )
        return designs

    def _resolve(self) -> str:
        from ..routes import configured_path

        # Asked first, for the same reason the explicit executable= below is
        # honoured: the pins mean this almost never lives on our PATH.
        configured = configured_path(self.env_var) if self.env_var else None
        if configured is not None:
            return str(configured)
        found = shutil.which(self.executable)
        if found:
            return found
        # An absolute path into another environment is the supported way to
        # work around the dependency pins.
        if Path(self.executable).is_file():
            return str(Path(self.executable))
        raise MissingDependencyError(
            f"{self.executable!r} is not on PATH and ${self.env_var} is not "
            "set. Install the optional extra "
            "with `uv pip install -e '.[mpnn]'` (or `pip install ligandmpnn`). "
            "Its pins (torch==2.2.1, numpy==1.23.5, Python <3.12) conflict "
            "with this package's, so the usual fix is a separate environment: "
            "create one, install ligandmpnn there, and point the adapter at "
            "it with dataclasses.replace(protein_mpnn, executable='/path/to/"
            "that/env/bin/ligandmpnn'). To run the pipeline without a "
            "sequence model, use model_name='dummy_sequences'."
        )

    def _args(
        self, num_seqs: int, temperature: float, kwargs: Mapping[str, Any]
    ) -> list[str]:
        # num_seqs is batch_size x number_of_batches. One batch is the simple
        # reading; pass either flag yourself when memory is tight (on CPU the
        # tool wants batch_size 1) and num_seqs steps aside.
        batching: dict[str, Any] = {}
        if "batch_size" not in kwargs and "number_of_batches" not in kwargs:
            batching = {"batch_size": num_seqs, "number_of_batches": 1}
        return format_kwargs(
            {"temperature": temperature, **batching, **kwargs},
            "flag",
            self.aliases,
        )

    def _require_ligand(self, backbone: Design) -> None:
        """LigandMPNN without a ligand is just unconditioned design.

        Base RFdiffusion does not carry ligands into its output, so this would
        otherwise fail silently: plausible sequences, and the ligand context
        that motivated choosing this model entirely absent.
        """
        if not het_residues(backbone.path):
            raise MissingLigandError(
                f"model_type 'ligand_mpnn' needs ligand context, but "
                f"{backbone.path} has no non-solvent HETATM records "
                f"(design {backbone.id!r}, from {backbone.producer}). "
                "RFdiffusion output does not carry ligands through, so the "
                "ligand has to come from the input structure. Use "
                "model_name='protein_mpnn' or 'soluble_mpnn' for ligand-free "
                "backbones."
            )

    def _split(
        self, source: Path, output_dir: Path, backbone: Design, temperature: float
    ) -> list[Design]:
        """One Design per sampled sequence, `path` on its own FASTA."""
        designs: list[Design] = []
        for record in source.read_text().split(">"):
            header, _, body = record.partition("\n")
            if not header.strip() or not body.strip():
                continue
            scores = header_scores(header)
            # Only sampled records carry a confidence; the first record is the
            # input sequence.
            if "overall_confidence" not in scores:
                continue
            index = int(scores.get("id", len(designs) + 1))
            path = source.with_name(f"{source.stem}_seq{index}.fa")
            design_id = make_id(path, output_dir)
            # The record name, not the file name, is what the folding step names
            # its outputs after. It carries the run key so that two sequence
            # arms over the same backbones cannot collide in one fold run.
            path.write_text(f">{design_id}\n{body.strip()}\n")
            meta: dict[str, Any] = {
                "model_type": self.model_type,
                "temperature": temperature,
                "fasta": str(source),
            }
            packed = self._packed_pdb(output_dir, source.stem, index)
            if packed is not None:
                meta["packed_pdb"] = str(packed)
            designs.append(
                Design(
                    id=design_id,
                    path=path,
                    parent_id=backbone.id,
                    producer=self.model_type,
                    scores=scores,
                    meta=meta,
                )
            )
        return designs

    @staticmethod
    def _packed_pdb(output_dir: Path, stem: str, index: int) -> Path | None:
        """The side-chain-packed structure, when --pack_side_chains was on.

        The directory name is not pinned down in the tool's docs, so this looks
        for a packed/ directory and then for a *_packed* file, and records
        nothing if neither is there rather than guessing.
        """
        for candidate in sorted((output_dir / "packed").glob(f"{stem}*.pdb")):
            return candidate
        for candidate in sorted(output_dir.rglob(f"{stem}*packed*.pdb")):
            return candidate
        return None


# One class, one instance per model_type: the registry maps three names onto
# these, so model_type never has to travel through the generic kwargs path.
protein_mpnn = MPNNAdapter(model_type="protein_mpnn")
ligand_mpnn = MPNNAdapter(model_type="ligand_mpnn")
soluble_mpnn = MPNNAdapter(model_type="soluble_mpnn")
