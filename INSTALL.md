# Installing and running `design_pipeline`

Affinity maturation in five functions, each model swappable by name. The core
package is small and pure Python; the models run in containers or as an
optional extra, so you install only what you intend to run.

Order matters: do step 1 and 2 on your laptop first. They need no GPU, no
Docker and no network, and they tell you the package works before you start
fighting with CUDA.

---

## 1. Python environment

Python 3.11 or newer. We use [uv](https://docs.astral.sh/uv/).

```bash
uv venv --python 3.11
uv pip install -e '.[dev]'
```

Activate it with `source .venv/bin/activate` (Linux/macOS) or
`.venv\Scripts\activate` (Windows PowerShell). You can also skip activation and
call `.venv/bin/python` / `.venv\Scripts\python.exe` directly.

Without uv:

```bash
python3.11 -m venv .venv
.venv/bin/pip install -e '.[dev]'
```

What you just installed:

| | |
|---|---|
| core | `numpy`, `biotite` for structures; `pandas` + `pyarrow` for the variant table. All light. |
| `[dev]` | `pytest` |
| `[numbering]` | `anarci`, for the `imgt` column |

No torch, no jax, nothing CUDA. Importing `design_pipeline` imports no model.
The container adapters — every scorer, every gate check, RFdiffusion — add
**no** Python dependency at all. The model extras `[mpnn]` and `[colabfold]`
are heavy and are covered in section 4; install them only if you intend to run
those models locally.

`[numbering]` is optional and needs HMMER on `PATH`. Without it the `imgt`
column is empty and a warning says so once; nothing else changes, because every
join in this package is on chain positions. IMGT labels matter when you want a
result to carry across antibodies.

---

## 2. Check it works, with no GPU and no Docker

The `dummy_*` models write small fixture files and plausible scores. They exist
so the whole loop is testable on a laptop. Run the test suite:

```bash
.venv/bin/python -m pytest -q
```

197 tests, about 30 seconds, fully offline. Then run the real pipeline against
fixtures:

```bash
.venv/bin/python examples/maturation.py --dummy
```

That takes a fixture complex, enumerates every single mutation of its heavy
chain, scores them, gates the sequence liabilities, builds 20 pairs at least
10 A apart and scores those. Run it twice: the second run resumes from cache
and does no work. The validation harness runs offline too, on the two fixture
systems that ship with the package:

```bash
.venv/bin/python -c "from design_pipeline.harness import run_systems, report; print(report(run_systems(models=['dummy_scorer'])))"
```

If both of those pass, the package is fine and anything that breaks later is
the environment, not the code.

---

## 3. Container runtime (for the scorers, the gate checks and RFdiffusion)

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

### Chai-1, Boltz-2, Protenix — one extra each, weights download on first run

```bash
uv pip install -e '.[chai]'      # the default folding model
uv pip install -e '.[boltz]'     # Boltz-2
uv pip install -e '.[protenix]'  # Protenix, pinned to 2.0.0
```

All three need **Linux and a CUDA GPU**, and Chai-1 additionally needs bfloat16
(compute capability 8.0 or newer, i.e. Ampere and later). The adapters check
this up front and raise `MissingGPUError` naming what is missing, rather than
failing somewhere inside the model. Licences and MSA requirements are in the
model table in [README.md](README.md) — read it before using Protenix for
anything that ships.

Each tool keeps its own weights cache; there is deliberately no shared one:

| Model | Cache | Notes |
|---|---|---|
| `chai` | `$CHAI_DOWNLOADS_DIR` | Unset, chai_lab writes into its own installation directory. Set it. |
| `boltz` | `$BOLTZ_CACHE` (default `~/.boltz`) | Must be absolute. Weights, CCD and the molecule archive are several GB; `boltz predict` fetches what is missing on the first run, and the adapter says so before it starts. |
| `protenix` | `$PROTENIX_ROOT_DIR` (default `$HOME`) | Checkpoints and CCD data download on first use. |

Two things to know about the arguments:

- Chai-1 runs **in process** (the adapter imports `chai_lab`, and with it
  torch), because its command line throws away the PAE matrix that
  `pae_interaction` needs. It therefore has no `extra_args`; pass
  `run_inference` arguments as keywords, e.g.
  `get_structure(out, model_name="chai1", antibody_seq=..., antigen_seq=...,
  num_trunk_recycles=5)`.
- Protenix picks the checkpoint: `checkpoint="protenix_mini_default_v0.5.0"`
  is much faster at comparable quality, and the `protenix_mini_esm_v0.5.0` /
  `protenix_mini_ism_v0.5.0` variants need no MSA. Add
  `use_default_params=True` to get the cycle and step counts a Mini checkpoint
  was tuned for.

To run step 1 with no GPU and no network, use `dummy_complex`; the de novo
branch has the fixture folders `dummy_chai`, `dummy_boltz` and
`dummy_protenix`.

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
  get_structure(out, model_name="colabfold", antibody_seq=..., antigen_seq=...,
                host_url="http://localhost:8080")
  ```

If `colabfold_batch` is not on `PATH`, the pipeline raises
`MissingDependencyError` naming the extra to install.

### Scorer and gate images

Every scorer except `consensus` and every gate check except `liabilities` runs
in a container, which is why they add no Python dependency. Each adapter is
plain data, so the image name, the flag the candidates CSV arrives on, and the
column the result is read from are all overridable without touching this
package:

```python
from dataclasses import replace
from design_pipeline import REGISTRY
from design_pipeline.adapters import _scorers

REGISTRY["antifold"] = replace(
    _scorers.antifold,
    container=replace(_scorers.antifold.container, image="ghcr.io/ourlab/antifold:2.1"),
    value_column="log_odds",
)
```

Or per call, when it is only the column names that moved:

```python
score_variants(complex_, variants, models=["antifold"],
               id_column="mutation", value_column="llr")
```

The images in the registry are the published ones where a published one exists
and the lab's own build where it does not (`foldx` needs a licence,
`alphafold3` needs weights DeepMind distributes separately). Verify the name
before assuming a pull will work; a 404 means build it locally and override as
above.

The adapter contract is small enough to wrap anything: the container is handed
a structure, a `candidates.csv` of `variant_id,chain,position,wt,mut`, and an
output directory, and has to write a CSV with one row per `variant_id`. A
candidate with several rows is a multi-mutation candidate and must be scored as
one mutant — never as the sum of its parts.

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

### The five steps

```python
from design_pipeline import (
    get_structure, enumerate_variants, score_variants, gate, select_pairs,
)

complex_ = get_structure("runs/her2_v1", pdb_id="1N8Z", chains=["A", "B", "C"])
variants = enumerate_variants(complex_, chain="B")
scores = score_variants(
    complex_, variants, "runs/her2_v1/scores",
    models=["antifold", "esm_if", "consensus"],
)
survivors = gate(scores, checks=["liabilities"])
pairs = select_pairs(survivors, max_pairs=20, min_distance=10.0,
                     by="score_antifold")
pair_scores = score_variants(complex_, pairs, models=["antifold", "gearbind"])
```

`get_structure` returns a `Complex`: the structure files it wrote plus the
resolved config. The other four return a `VariantTable`: a parquet file of
`variant_id, chain, position, imgt, wt, mut` plus one `score_<model>` column
per scorer. Both are frozen, both hold paths, and nothing is held in memory
between steps — read a table with `.frame` when you want the pandas object.

A pair is two rows sharing one `variant_id`, in the same schema, which is why
`score_variants` is simply called again on it. There is no code path anywhere
that adds two single scores together to estimate a pair, and `select_pairs`
drops the score columns rather than carry the singles' numbers forward.

When the target has no solved structure, step 1 predicts one instead:

```python
complex_ = get_structure(
    "runs/unknown_v1", model_name="alphafold3",
    antibody_seq="EVQLVESGGG...", antigen_seq="TQVCTGTDMK...",
    num_seeds=20,
)
print(len(complex_.paths), complex_.path)   # every seed, not only the top one
```

### Costs, caps and the cache

`score_variants` caches per scorer, on (structure content, candidate set, model
name, kwargs), so adding a model to `models=[...]` does not invalidate the ones
that already ran. The matrix scorers see the whole table — one forward pass over
the complex answers every single mutation. The per-variant scorers build and
relax a structure per candidate, so they are capped:

```python
# gearbind on 4,598 candidates is a week of compute, so this raises
score_variants(complex_, variants, models=["gearbind"])

# rank on the column you decided to trust, then cap
best = survivors.frame.sort_values("score_antifold", ascending=False)
score_variants(complex_, survivors, models=["gearbind"], max_variants=200)
```

### The validation harness

```python
from design_pipeline.harness import run_systems, report

results = run_systems(
    systems=["1n8z_her2", "4fqi_flu"],
    models=["antifold", "esm_if", "gearbind", "consensus"],
    output_dir="runs/harness",
)
print(report(results))        # Spearman and precision@10, per scorer per system
```

Structures download; the measured affinities you supply once, as a CSV of
`variant_id,measured_affinity` under `data/ground_truth/` (or wherever
`$DESIGN_PIPELINE_GROUND_TRUTH` points). Sources are FLAb, AbBiBench and
AB-Bind — `design_pipeline.harness.systems` names them per system, and so does
the error you get for a missing file.

Variant ids are 1-based into the chain as given in the structure file, so a set
numbered in IMGT or Kabat needs converting first. `enumerate_variants(only=...)`
checks every wild-type residue against the structure and refuses the rows that
disagree, which is the one guard against a correlation computed at the wrong
positions.

### Any model flag, without touching this package

Keyword arguments are mapped to the model's own CLI flags — Hydra-style
`key=value` for RFdiffusion, `--flag value` for the others. A small alias table
covers the ergonomic names above; anything else passes straight through. For
something the mapping does not reach, use `extra_args`:

```python
score_variants(complex_, variants, models=["antifold"],
               extra_args=["--batch_size", "16"])
```

`**kwargs` and `extra_args` go to every model named in one call, so a flag one
model rejects is an error for the whole call; the remedy is one call per model,
which joins into the same table anyway.

### Swapping a model

```python
score_variants(complex_, variants, models=["dummy_scorer"])
gate(scores, checks=["dummy_check"])
```

The registry is a plain dict, so adding your own model is one assignment and no
entry points:

```python
from design_pipeline import REGISTRY
REGISTRY["my_scorer"] = MyScorer()     # run(inputs, output_dir, **kwargs)
```

Registered names, by stage:

| Stage | Names |
|---|---|
| structure | `pdb_download`, `alphafold3`, `boltz2`, `chai1`, `chai`, `boltz`, `protenix`, `colabfold` |
| scorers | `antifold`, `esm_if`, `abmpnn`, `proteinmpnn_score`, `consensus`, `gearbind`, `rosetta_ddg`, `foldx` |
| gate | `liabilities`, `thermompnn`, `aggrescan3d`, `netsolp` |
| de novo | `rfdiffusion`, `protein_mpnn`, `ligand_mpnn`, `soluble_mpnn` |
| fixtures | `dummy_complex`, `dummy_scorer`, `dummy_per_variant`, `dummy_check`, and the de novo `dummy_*` set |

An unknown name lists the known ones.

### The de novo branch

RFdiffusion and the MPNN variants are kept, demoted off the default path. The
one entry point is `diffuse_loops`, which rebuilds a CDR at its original length
and returns an ordinary variant table:

```python
from design_pipeline import diffuse_loops

loops = diffuse_loops(complex_, cdr="H3", num_designs=100)
scores = score_variants(complex_, loops, models=["antifold"])
```

Those variants differ from the parent at ten or twenty positions at once, far
outside the edit-distance range where any of these scorers was validated, and
the harness cannot tell you about that regime either, because the published
affinity datasets do not cover it. Treat the ranking as a hypothesis generator.

The old three-step API is still there under `design_pipeline.denovo`
(`make_protein_backbones`, `design_sequences`, `predict_structure`), with
`score_designs` in `design_pipeline.scoring` and `select` in
`design_pipeline.designs`. See [examples/refine.py](examples/refine.py) and
[examples/fanout.py](examples/fanout.py):

```bash
.venv/bin/python examples/refine.py --dummy    # fixtures, offline
```

---

## 6. Output layout, caching and seeds

Each run writes to `output_dir/<key>/`. For `get_structure` and the de novo
branch, `<key>` hashes the resolved config, the seed, and the content and ids
of the input files. For `score_variants` and `gate` it hashes (structure
content, candidate set, model name, kwargs) — per scorer rather than per step,
which is what keeps a finished column alive when the scorer set changes:

```
runs/her2_v1/4f1c90ab77e2d305/
    1n8z.pdb
    DONE                      <- the records: paths, scores, parent ids
runs/her2_v1/variants/B_singles.parquet
runs/her2_v1/scores/
    9ad3e7c1b4028f66/         <- antifold over these candidates
        candidates.csv
        scores.csv
        DONE
    c70b1e55d9a34f82/         <- esm_if over the same ones
    scored_B_singles.parquet  <- the joined table
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
| `MissingDependencyError` | A model's package or entry point is missing (`colabfold_batch`, `ligandmpnn`, `boltz`, `protenix`, `chai_lab`). The message names the extra; for MPNN see the separate-environment note in section 4. |
| `MissingGPUError` | `chai`, `boltz` or `protenix` on a non-Linux host, a node with no CUDA GPU, or (Chai-1 only) a card older than Ampere. Move to a GPU node, pick another model, or use the matching `dummy_*` model. |
| `MissingLigandError` | `ligand_mpnn` on a backbone with no ligand. Use `protein_mpnn` or `soluble_mpnn`, or supply a structure that carries the ligand. |
| `UnknownModelError` | Typo in a `model_name`, a `models=[...]` entry or a `checks=[...]` entry; the message lists the registered ones. |
| `MissingDependencyError: ... OAS` | `consensus` needs a frequency table. Point `frequencies=` at a CSV of `imgt,aa,frequency`, set `DESIGN_PIPELINE_OAS_TABLE`, or drop it from `models=[...]`. |
| `ValueError: ... max_variants` | A per-variant scorer on an uncapped table of thousands. Rank on the column you trust, then pass `max_variants=n`. |
| `ValueError: ... does not match` | A variant id whose wild-type residue is not what the structure holds at that position: the numbering differs (IMGT or Kabat rather than chain positions), or it is the wrong chain. |
| `gate dropped every ... candidates` | A threshold too strict for a weak predictor. The reasons are in the annotated table the warning names. |
| `RunFailedError: exit code N` | The model itself failed. Its output was streamed above the traceback — read that, not the traceback. |
| `RunFailedError: ... wrote nothing matching` | The container ran and exited 0 but produced no output files. Usually a wrong output path or a silently skipped input. |
| Docker pull 404 | The image name needs verifying; build locally and override the adapter (section 4). |
| CUDA / no GPU inside the container | Test with the `nvidia-smi` container commands in section 3 before suspecting the pipeline. |
| A step you expected to run was skipped | Cache hit. Change the seed, or pass `cache=False`. |
