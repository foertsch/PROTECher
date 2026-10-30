"""AlphaFold3, the fallback for step 1 when no solved structure exists.

Two things make this more than a `ContainerAdapter` entry. The input is a JSON
file rather than flags, so it is written here from the FASTA the API hands over;
and `num_seeds` has to come back as an *ensemble*. Ranking of antibody-antigen
complexes is unreliable - the top-ranked seed is regularly not the closest to
the crystal structure - so every seed becomes a path on the `Complex` and
downstream code can score across them instead of trusting the ranking.

The image is the user's own build: DeepMind distributes the code and the
weights separately, and the weights are not redistributable.
"""

from __future__ import annotations

import json
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any, Mapping, Sequence

from ..containers import ContainerAdapter, make_id
from ..routes import require
from ..designs import Design
from ..errors import RunFailedError
from ._folding import cif_to_pdb, fasta_chains

DEFAULT_IDS = ("H", "L", "A", "B", "C", "D")
"""Chain ids written into the JSON when the caller does not name them: the
antibody chains first, then the antigen. AlphaFold3 honours them, so a complex
predicted this way can be addressed as chain='H' the same way a downloaded one
can. The other predictors label chains themselves."""


def write_input_json(
    path: Path,
    name: str,
    chains: Sequence[str],
    *,
    chain_ids: Sequence[str] | None = None,
    num_seeds: int = 1,
    seed: int = 0,
) -> Path:
    """AlphaFold3's input JSON for one complex.

    Seeds are consecutive from `seed`, so the ensemble is reproducible and the
    cache key over (config, seed) actually describes what came out.
    """
    ids = list(chain_ids or DEFAULT_IDS)
    if len(ids) < len(chains):
        raise ValueError(
            f"{len(chains)} chains to predict and only {len(ids)} chain ids; "
            f"pass chain_ids=(...) with one per chain"
        )
    payload = {
        "name": name,
        "modelSeeds": [seed + i for i in range(max(1, int(num_seeds)))],
        "sequences": [
            {"protein": {"id": ids[i], "sequence": sequence}}
            for i, sequence in enumerate(chains)
        ],
        "dialect": "alphafold3",
        "version": 1,
    }
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(payload, indent=1))
    return path


@dataclass(frozen=True)
class AlphaFold3Adapter:
    """Satisfies the ModelAdapter Protocol; the container work is delegated."""

    name: str = "alphafold3"
    image: str | None = None
    """No default: the weights come from DeepMind on request and the image is
    your own build, so there is no name this package could ship. `run` raises
    with that route until you set one."""
    model_dir: str = "/root/models"
    db_dir: str = "/root/public_databases"
    aliases: Mapping[str, str] = field(default_factory=dict)

    def run(
        self,
        inputs: Sequence[Design],
        output_dir: Path,
        *,
        num_seeds: int = 1,
        chain_ids: Sequence[str] | None = None,
        seed: int = 0,
        extra_args: Sequence[str] | None = None,
        **kwargs: Any,
    ) -> list[Design]:
        if not inputs:
            raise ValueError(
                "alphafold3 needs the sequences to fold; get_structure writes "
                "them to a FASTA from antibody_seq / antigen_seq, or from the "
                "one you pass as fasta="
            )
        if self.image is None:
            require(self.name)
        fasta = inputs[0]
        name, chains = fasta_chains(fasta.path)
        json_path = write_input_json(
            output_dir / "input.json",
            name,
            chains,
            chain_ids=chain_ids,
            num_seeds=num_seeds,
            seed=seed,
        )
        container = ContainerAdapter(
            image=self.image,
            command=("python", "run_alphafold.py"),
            output_template=("--output_dir", "{path}"),
            fixed_args=(
                f"--model_dir={self.model_dir}",
                f"--db_dir={self.db_dir}",
            ),
            aliases=self.aliases,
            path_kwargs=("json_path",),
            output_globs=("*.cif",),
            # The JSON carries the input, so there is one invocation whatever
            # happens; the FASTA is not passed on the command line.
            per_input=False,
            gpu=True,
        )
        produced = container.run(
            [], output_dir, json_path=str(json_path), extra_args=extra_args, **kwargs
        )
        cifs = sorted(
            (d.path for d in produced),
            # The top-level ranked model first, then the per-seed directories,
            # so Complex.path is the tool's own pick and Complex.paths is the
            # whole ensemble behind it.
            key=lambda p: (len(p.relative_to(output_dir).parts), p.name),
        )
        if not cifs:
            raise RunFailedError(f"{self.image} wrote no mmCIF under {output_dir}")

        designs: list[Design] = []
        for cif in cifs:
            pdb = cif_to_pdb(cif)
            designs.append(
                Design(
                    id=make_id(pdb, output_dir),
                    path=pdb,
                    parent_id=fasta.id,
                    producer=self.name,
                    meta={"cif": str(cif), "seed": seed, "num_seeds": int(num_seeds)},
                )
            )
        return designs


alphafold3 = AlphaFold3Adapter()
