# Installing and running `design_pipeline`

Three model steps plus scoring, each swappable by name. The core package is
small and pure Python; the models run in containers or as an optional extra, so
you install only what you intend to run.

Order matters: do step 1 and 2 on your laptop first. They need no GPU, no
Docker and no network, and they tell you the package works before you start
fighting with CUDA.

---

## 1. Python environment

Python 3.11 or newer. We use [uv](https://docs.astral.sh/uv/).

```bash
uv venv --python 3.11
uv pip install -e '.[dev,analysis]'
```

Activate it with `source .venv/bin/activate` (Linux/macOS) or
`.venv\Scripts\activate` (Windows PowerShell). You can also skip activation and
call `.venv/bin/python` / `.venv\Scripts\python.exe` directly.

Without uv:

```bash
python3.11 -m venv .venv
.venv/bin/pip install -e '.[dev,analysis]'
```

What you just installed:

| | |
|---|---|
| core | `numpy`, `biotite` — needed by scoring, both light |
| `[dev]` | `pytest`, `pandas` |
| `[analysis]` | `pandas`, for `to_dataframe` |

No torch, no jax, nothing CUDA. Importing `design_pipeline` imports no model.
RFdiffusion runs in a container and adds **no** Python dependency. The two
model extras, `[mpnn]` and `[colabfold]`, are heavy and are covered in section
4 — install them only if you intend to run those models locally.

---

## 2. Check it works, with no GPU and no Docker

The `dummy_*` models write small fixture files and plausible scores. They exist
so the whole loop is testable on a laptop. Run the test suite:

```bash
.venv/bin/python -m pytest -q
```

68 tests, about 20 seconds, fully offline. Then run the real round loop against
fixtures:

```bash
.venv/bin/python examples/refine.py --dummy
```

That generates backbones, designs sequences, "folds" them, scores, selects the
best 8 and re-diffuses — three rounds — and prints the top 10 as a dataframe.
Run it twice: the second run resumes from cache and does no work.

If both of those pass, the package is fine and anything that breaks later is
the environment, not the code.

---

## 3. Container runtime (for RFdiffusion and ProteinMPNN)

You need **either** apptainer or docker. The adapter detects which is present
and prefers apptainer, because that is the one shared clusters allow.

- **apptainer** (clusters, Linux): [install docs](https://apptainer.org/docs/admin/main/installation.html).
  On a cluster, try `module load apptainer` or `module avail apptainer` first.
- **docker** (workstation): [install docs](https://docs.docker.com/get-docker/).
  On Windows, Docker Desktop with the WSL2 backend. Apptainer is Linux-only.

For GPU access:

```bash
nvidia-smi                 # the driver must be there first
```

- docker additionally needs the [NVIDIA Container Toolkit](https://docs.nvidia.com/datacenter/cloud-native/container-toolkit/latest/install-guide.html);
  the adapter passes `--gpus all`.
- apptainer passes `--nv` and needs no extra package.

Check the GPU is visible *inside* a container before blaming the pipeline:

```bash
docker run --rm --gpus all nvidia/cuda:12.4.0-base-ubuntu22.04 nvidia-smi
apptainer run --nv docker://nvidia/cuda:12.4.0-base-ubuntu22.04 nvidia-smi
```

If no runtime is on `PATH`, the pipeline raises `MissingRuntimeError` naming
both options rather than failing somewhere deep in a subprocess.

---

## 4. Model images and weights

### RFdiffusion — weights ship inside the image

The adapter uses image `rosettacommons/rfdiffusion` and passes
`inference.model_directory_path=/app/RFdiffusion/models`, which is where the
weights live inside that image. Nothing to download separately.

```bash
docker pull rosettacommons/rfdiffusion
# apptainer pulls on first use and caches it; to pre-pull:
apptainer pull docker://rosettacommons/rfdiffusion
```

> **Verify this image name before the event.** It comes from our brief, not
> from a pull we have done. RFdiffusion's repository ships a `Dockerfile` you
> are expected to build yourself, and the weights are downloaded by a script
> rather than baked in. If the pull above 404s, build it locally and point the
> adapter at your image — see *Overriding an adapter* below.

Set a big cache dir if `$HOME` is small, which it usually is on a cluster:

```bash
export APPTAINER_CACHEDIR=/scratch/$USER/apptainer
```

### MPNN sequence design — pip, no container, parameters bundled

Sequence design uses the upstream [`ligandmpnn`](https://pypi.org/project/ligandmpnn/)
package (MIT, model parameters in the wheel). One codebase serves three models,
selected by `--model_type`, and they are registered as three names:

| `model_name` | `--model_type` | For |
|---|---|---|
| `protein_mpnn` | `protein_mpnn` | the default |
| `soluble_mpnn` | `soluble_mpnn` | sequences biased away from aggregation |
| `ligand_mpnn` | `ligand_mpnn` | backbones that carry a ligand (see the guard below) |

```bash
uv pip install -e '.[mpnn]'
ligandmpnn --help
```

> **Expect this extra not to resolve in the same environment.** `ligandmpnn`
> 0.1.2 pins `numpy==1.23.5`, `torch==2.2.1` and Python `<3.12`, while biotite
> needs `numpy>=2.1`. The adapter shells out to the `ligandmpnn` console
> script rather than importing it, so the supported fix is a second
> environment:
>
> ```bash
> uv venv .venv-mpnn --python 3.11
> uv pip install --python .venv-mpnn ligandmpnn
> ```
>
> then point the adapter at it, once, in your script:
>
> ```python
> import dataclasses
> from design_pipeline import REGISTRY
> from design_pipeline.adapters import protein_mpnn
>
> REGISTRY["protein_mpnn"] = dataclasses.replace(
>     protein_mpnn, executable=".venv-mpnn/bin/ligandmpnn"   # Scripts/ on Windows
> )
> ```

Checkpoint variants are reachable as ordinary kwargs. Higher training noise
tolerates imperfect backbones better, which matters for diffused ones:

```python
design_sequences(
    backbones, num_seqs=8, temperature=0.2,
    checkpoint_protein_mpnn="proteinmpnn_v_48_020.pt",   # 0.20 A noise
)                                                        # _010 = 0.10, _002 = 0.02
```

**The ligand guard.** `ligand_mpnn` raises `MissingLigandError` if the input PDB
has no non-solvent `HETATM` records. This is deliberate: base RFdiffusion does
not carry ligands into its output, so without the check LigandMPNN silently
degrades to unconditioned design — plausible sequences, with the ligand context
that motivated choosing the model entirely absent. The ligand has to come in on
the input structure; there is no transplanting helper.

### ColabFold — optional Python extra, weights download on first run

```bash
uv pip install -e '.[colabfold]'
```

This one is heavy (it brings jax), which is exactly why it is an extra. The
adapter shells out to `colabfold_batch`, so after installing, check:

```bash
colabfold_batch --help
```

- **Weights**: the AlphaFold2 parameters download automatically on the first
  prediction and are cached; budget a few GB and a slow first run. See
  `colabfold_batch --help` for the flag that relocates that cache in your
  version, and pass it through as shown below.
- **MSAs**: the public MSA server `https://api.colabfold.com` is the default.
  It is rate-limited and it means your sequences leave the machine. To point at
  a local server later, pass `host_url`:

  ```python
  predict_structure(seqs, host_url="http://localhost:8080")
  ```

If `colabfold_batch` is not on `PATH`, the pipeline raises
`MissingDependencyError` naming the extra to install.

### Overriding an adapter

Adapters are frozen dataclasses in a plain dict registry, so changing an image
or a weights path is a few lines in your own script — never an edit to this
package.

```python
import dataclasses
from design_pipeline import REGISTRY
from design_pipeline.adapters import rfdiffusion

REGISTRY["rfdiffusion"] = dataclasses.replace(rfdiffusion, image="rfdiffusion:local")
```

If you built the image yourself and the weights are on the host rather than
inside it, route the weights directory in as a path argument: anything named in
`path_kwargs` gets bind-mounted and rewritten to its in-container path.

```python
from design_pipeline.adapters._rfdiffusion import ALIASES

REGISTRY["rfdiffusion"] = dataclasses.replace(
    rfdiffusion,
    image="rfdiffusion:local",
    fixed_args=(),                      # drop the in-image weights path
    path_kwargs=("input_pdb", "model_directory_path"),
    aliases={**ALIASES, "model_directory_path": "inference.model_directory_path"},
)

make_protein_backbones(
    "inputs/target.pdb",
    model_directory_path="/shared/models/rfdiffusion",
    contigs="[10-40/A163-181/10-40]",
    num_designs=100,
)
```

Note that the *parent* of each path argument is what gets bind-mounted, so keep
weights in a directory you are happy to expose.

---

## 5. Running the pipeline

### The three steps

```python
from design_pipeline import (
    make_protein_backbones, design_sequences, predict_structure,
    score_designs, select, to_dataframe,
)

backbones = make_protein_backbones(
    "inputs/target.pdb",
    output_dir="runs/v1/bb",
    num_designs=100,
    contigs="[10-40/A163-181/10-40]",
    hotspot_res=["A59", "A83", "A91"],
    seed=0,
)

seqs  = design_sequences(backbones, "runs/v1/seqs", num_seqs=8)
preds = predict_structure(seqs, "runs/v1/fold")
scored = score_designs(preds, reference=backbones)

best = select(scored, top_k=10, by="pae_interaction")
print(to_dataframe(scored).sort_values("pae_interaction").head(20))
```

Every function returns a plain `list[Design]`. Each `Design` holds an `id`, a
`path` on disk, the `parent_id` it came from, its `producer`, a `scores` dict
and free-form `meta`. Nothing is held in memory between steps.

### Any model flag, without touching this package

Keyword arguments are mapped to the model's own CLI flags — Hydra-style
`key=value` for RFdiffusion, `--flag value` for the others. A small alias table
covers the ergonomic names above; anything else passes straight through. For
something the mapping does not reach, use `extra_args`:

```python
make_protein_backbones(
    "inputs/target.pdb",
    num_designs=100,
    extra_args=["denoiser.noise_scale_ca=0.5"],
)
```

### Swapping a model

```python
seqs = design_sequences(backbones, "runs/v1/seqs", model_name="dummy_sequences")
```

Registered: `rfdiffusion`, `protein_mpnn`, `ligand_mpnn`, `soluble_mpnn`,
`colabfold`, `dummy_backbones`, `dummy_sequences`, `dummy_structures`. An
unknown name lists the known ones.

Running several sequence-design variants in one round needs no API support —
the steps return plain lists, so the arms concatenate and scoring sees them
together. See [examples/fanout.py](examples/fanout.py).

### Iterative refinement

The round loop is a `for` statement in your own script, not a library feature.
See [examples/refine.py](examples/refine.py):

```python
for rnd in range(n_rounds):
    seqs      = design_sequences(backbones, num_seqs=8)
    preds     = predict_structure(seqs)
    scored    = score_designs(preds, reference=backbones)
    survivors = select(scored, top_k=8, by="pae_interaction")
    backbones = make_protein_backbones(survivors, partial_T=20, num_designs=12, seed=rnd)
```

`make_protein_backbones` takes either a target path (fresh generation) or a
`list[Design]` (partial diffusion around them). `partial_T` is a raw timestep
count relative to `diffuser.T`: published values like 80 assume the original
200-step model, so at today's default `T=50` the equivalent is roughly 20.

Run it:

```bash
.venv/bin/python examples/refine.py            # real models
.venv/bin/python examples/refine.py --dummy    # fixtures, offline
```

---

## 6. Output layout, caching and seeds

Each step writes to `output_dir/<key>/`, where `<key>` hashes the resolved
config, the seed, and the content and ids of the input files:

```
runs/v1/bb/2e3eb49f1a0c4b77/
    design_0.pdb
    design_1.pdb
    DONE                     <- the design records: paths, scores, parent ids
```

- A directory with a `DONE` marker is skipped and its records returned, so
  rerunning a script resumes where it stopped. That is the whole reason we do
  not need a workflow engine.
- `cache=False` on any step forces a rerun.
- **The seed is part of the key.** Asking for another 100 designs with the same
  config and a new `seed` really runs them, instead of handing back the cached
  100.
- Deleting output files turns the entry back into a miss.
- `runs/` is git-ignored. Model weights and structures over a few MB do not go
  in git either — see the working rules in [README.md](README.md).

Subprocess stdout and stderr stream straight through, so a running job is
visible and the command being run is echoed to stderr.

---

## 7. Troubleshooting

| What you see | What it means |
|---|---|
| `MissingRuntimeError` | Neither apptainer nor docker on `PATH`. Install one, or use the `dummy_*` models. |
| `MissingDependencyError` | `colabfold_batch` or `ligandmpnn` is missing. The message names the extra; for MPNN see the separate-environment note in section 4. |
| `MissingLigandError` | `ligand_mpnn` on a backbone with no ligand. Use `protein_mpnn` or `soluble_mpnn`, or supply a structure that carries the ligand. |
| `UnknownModelError` | Typo in `model_name`; the message lists the registered ones. |
| `RunFailedError: exit code N` | The model itself failed. Its output was streamed above the traceback — read that, not the traceback. |
| `RunFailedError: ... wrote nothing matching` | The container ran and exited 0 but produced no output files. Usually a wrong output path or a silently skipped input. |
| Docker pull 404 | The image name needs verifying; build locally and override the adapter (section 4). |
| CUDA / no GPU inside the container | Test with the `nvidia-smi` container commands in section 3 before suspecting the pipeline. |
| A step you expected to run was skipped | Cache hit. Change the seed, or pass `cache=False`. |
