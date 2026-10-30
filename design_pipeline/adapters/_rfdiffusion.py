"""RFdiffusion: backbone generation and partial diffusion, same adapter.

Weights ship inside the image. Args are Hydra-style key=value, so the kwargs
mapping is mostly f"{k}={v}" with a small alias table for the ergonomic names.

On `partial_T`: it is a count of timesteps, relative to `diffuser.T`. Published
values like 80 assume the original 200-step model; at today's default T=50 the
equivalent is roughly 20. It is exposed raw and NOT converted to a 0-1 noise
fraction, because people copy numbers straight out of papers and a silent
rescaling would be worse than the arithmetic.
"""

from __future__ import annotations

from pathlib import Path

from ..containers import ContainerAdapter
from ..routes import RFDIFFUSION_IMAGE

ALIASES = {
    "input_pdb": "inference.input_pdb",
    "num_designs": "inference.num_designs",
    "contigs": "contigmap.contigs",
    "output_prefix": "inference.output_prefix",
    "hotspot_res": "ppi.hotspot_res",
    "partial_T": "diffuser.partial_T",
    "T": "diffuser.T",
    "seed": "inference.seed",
    "ckpt_override_path": "inference.ckpt_override_path",
    "symmetry": "inference.symmetry",
}


def collect(path: Path, run_dir: Path) -> list[tuple[Path, dict[str, float]]]:
    """One design per PDB, with whatever the sibling .trb holds."""
    return [(path, parse_trb(path.with_suffix(".trb")))]


def parse_trb(trb: Path) -> dict[str, float]:
    """Scores from RFdiffusion's trajectory file.

    The .trb is a pickled dict of numpy arrays whose contents vary by version
    and by inference mode, so this is best-effort: we take the final pLDDT if
    one is there and skip the rest. Anything missing simply means no score.
    """
    if not trb.is_file():
        return {}
    try:
        import numpy as np

        data = np.load(trb, allow_pickle=True)
        payload = data.item() if hasattr(data, "item") else dict(data)
        plddt = payload.get("plddt")
        if plddt is None:
            return {}
        arr = np.asarray(plddt, dtype=float)
        # Shape is (timesteps, residues) during diffusion; the last step is the
        # one that describes the structure that was written out.
        final = arr[-1] if arr.ndim > 1 else arr
        return {"rfdiffusion_plddt": float(final.mean())}
    except Exception:
        # A .trb we cannot read is not a failed run.
        return {}


# The one published image of the ten this package used to name, pinned by
# digest in routes.py rather than tracking a tag that can move under a
# campaign without changing the config.
rfdiffusion = ContainerAdapter(
    image=RFDIFFUSION_IMAGE,
    arg_style="hydra",
    # Partial diffusion feeds each input design in as inference.input_pdb;
    # fresh generation passes the target the same way, through path_kwargs.
    input_template=("inference.input_pdb={path}",),
    output_template=("inference.output_prefix={path}/design",),
    fixed_args=("inference.model_directory_path=/app/RFdiffusion/models",),
    aliases=ALIASES,
    path_kwargs=("input_pdb",),
    output_globs=("*.pdb",),
    collect=collect,
)
