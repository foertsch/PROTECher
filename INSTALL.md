# Installing and running `design_pipeline`

Affinity maturation in five functions, each model swappable by name. The core
package is small and pure Python; the models install into one of two
environments, or into the single container that genuinely needs one, so you
install only what you intend to run.

Order matters: do step 1 and 2 on your laptop first. They need no GPU, no
Docker and no network, and they tell you the package works before you start
fighting with CUDA.

At any point, in any environment:

```bash
python -m design_pipeline.preflight
```

prints one row per adapter — step, which environment it belongs to, whether it
resolves here, and what was probed — and exits non-zero when a model the
default campaign needs is missing. Run it in each environment before the event.
It imports nothing for real and pulls nothing, so it is safe on a login node.

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
That last sentence stays true once `[plm]` is installed too — the model imports
are inside the adapters that need them — and `[plm]` is what the default scorer
needs, so the `core` that `provision/core.sh` builds includes it with a CPU
torch. On a laptop you can leave it out; `ablang2` is then the one row
`preflight` reports unresolved.
The container adapters — every scorer, every gate check, RFdiffusion — add
**no** Python dependency at all. The model extras `[mpnn]` and `[colabfold]`
are heavy and are covered in section 4; install them only if you intend to run
those models locally.

`[numbering]` is optional and needs HMMER on `PATH`. Without it the `imgt`
column is empty and a warning says so once; nothing else changes, because every
join in this package is on chain positions. IMGT labels matter when you want a
result to carry across antibodies.

### Where each model goes

Models do not all fit together, so they do not all go in one place.
`design_pipeline/routes.py` is the machine-readable version of this table — the
errors and `preflight` both read it — and `provision/manifest.toml` is the
machine-readable version of the environments themselves.

| | holds | why |
|---|---|---|
| **`core`** | this package, `ablang2`, `antiberty`, `transformers` + `torch`, `fair-esm`, `anarci` | the sequence-only scorers and the whole table/gating path. Build and test this first. |
| **Env B** | `antifold`, `thermompnn`, and the co-folding models (`chai`, `boltz`, `protenix`, `colabfold`) | each pins torch against the others, and AntiFold pins `biotite==0.38` |
| **container** | `gearbind` | TorchDrug pins Python ≤3.10 and an old torch: the one conflict a second environment cannot fix |
| **container** | `rfdiffusion` | it already had a published image; pinned by digest |

"Env B" is one label for several environments, and the dry run has since settled
how many. Their pins forced five, each with its own script and its own key in
[provision/manifest.toml](provision/manifest.toml): `antifold` (Python 3.13),
`thermompnn` (3.10), `boltz`, `chai` and `esm_env`. AntiFold and ThermoMPNN do
not share one. Nothing in the code assumes a grouping either way — each model is
reached by its own path, which is what `executable=` is for — so the split cost
nothing but disk.

Adapters in another environment are reached by path, not by import — the
pipeline shells out:

```python
score_variants(complex_, variants, models=["antifold"],
               executable="/opt/envs/antifold/bin/python")
# or, once, for the whole session:
export DESIGN_PIPELINE_ANTIFOLD=/opt/envs/antifold/bin/python
export DESIGN_PIPELINE_THERMOMPNN=/opt/envs/thermompnn/bin/python
export THERMOMPNN_SRC=/opt/src/ThermoMPNN
```

`antifold` is run through the repo's `antifold_score` module, and its table route
needs a second variable: `ANTIFOLD_SRC`, the AntiFold source folder that holds
`models/model.pt` (or `antifold_src=` as a keyword argument). Without it the run
stops with `--route table needs --table, or --antifold_src (or ANTIFOLD_SRC)`.
`preflight` still lists `antifold` as resolved and adds `(ANTIFOLD_SRC not set)`
to its row. `antifold` scores single mutations only; for pairs and larger sets
use `antifold_seq`.

`thermompnn` needs the same pair of variables, for the same reason and then
one more: the repo sets `package = false` and ships no console script, so
`DESIGN_PIPELINE_THERMOMPNN` is that environment's **python**, not a binary,
and `THERMOMPNN_SRC` (or `src=`) is the clone, the folder holding
`analysis/custom_inference.py`. That script site-saturates one chain per call
and takes no mutation list, so the adapter runs it once per chain the
candidates touch and joins the result back on; a multi-mutation candidate gets
no value and therefore passes the gate unmeasured.

#### AntiFold cannot go in `core`

`pip install antifold` pins `torch==2.3.1`, `torch-geometric==2.4.0`,
`biotite==0.38`, `biopython==1.83` and `numpy==1.26.*`. This package needs
`biotite>=2.0`, because biotite 2.0 changed `superimpose()` and the scoring
relies on the new signature. Checked on 2026-10-30:

```
$ pip install --dry-run "biotite>=2.0" "numpy>=1.26" antifold
ERROR: Cannot install antifold==0.2.0, antifold==0.2.1, antifold==0.2.2,
antifold==0.2.3, antifold==0.3.0, antifold==0.3.1 and biotite>=2.0 because
these package versions have conflicting dependencies.
ERROR: ResolutionImpossible
```

That is every published version of AntiFold. So it gets its own environment and
is invoked by path. It is no longer in `DEFAULT_SCORERS` (2026-10-31, see
[CHANGELOG.md](CHANGELOG.md)), so a default `score_variants` call does not need
Env B at all; `models=["antifold"]` or `models=["antifold_seq"]` does, and the
adapters and routes are unchanged.

#### AntiFold 0.4.0 from GitHub, and `antifold_seq`

The PyPI releases are not the newest. The GitHub source, v0.4.0 (2026-09-27),
needs Python 3.13 or later but does not pin biotite or numpy: its requirements
are `biotite>=1.0`, `numpy>=1.26`, `torch>=2.3` and `pandas>=2.0`. On 2026-10-30
it installed from a clone into a Python 3.13 environment next to biotite 2.0.0,
numpy 2.5.4 and pandas 3.1.0, and ran on 1N8Z; `python -m
design_pipeline.preflight` run inside that environment reported the default
campaign resolved. So where `core` can be Python 3.13, AntiFold may not need an
environment of its own. That is one environment on one machine, checked with
the preflight and not through a full campaign.

`antifold_seq` needs that version. It is the whole-sequence counterpart of
`antifold`: AntiFold's own `score_sequence_in_complex`, one forward pass per
mutated chain, so a candidate with several mutations is scored as a whole and
not as a sum of singles. It disagrees with the matrix lookup on a single
mutation (20 trastuzumab singles: correlation 0.81), so it is a separate
column. It reads the same `DESIGN_PIPELINE_ANTIFOLD` and takes the chain ids as
keyword arguments:

```python
score_variants(complex_, variants, models=["antifold_seq"], heavy="B", light="A",
               weights="/path/to/antifold/models/model.pt")
```

The weights (567 MB) download on first use into AntiFold's own folder; pass
`weights=` to reuse a download. It is a per-variant scorer, about 0.7 s a
candidate on a CPU, so an uncapped table above `PER_VARIANT_LIMIT` is refused;
cap it with `max_variants=`.

#### ESM-IF, resolved empirically

Two separate questions, two different answers. The pure-Python part resolves
into `core` with no conflict:

```
$ pip install --dry-run "biotite>=2.0" fair-esm torch-geometric
Would install biotite-2.0.0 fair-esm-2.0.0 torch-geometric-2.8.1 numpy-2.5.4
  scipy-1.18.1 networkx-3.7 ... (no conflict)
```

`torch-scatter` does not:

```
$ pip install --dry-run "biotite>=2.0" fair-esm torch-geometric torch-scatter
  File "<string>", line 8, in <module>
ModuleNotFoundError: No module named 'torch'
error: subprocess-exited-with-error
Getting requirements to build wheel did not run successfully.
```

It ships no wheel for current torch, so pip builds it from the sdist, and its
`setup.py` imports `torch` at build time inside an isolated build environment
where torch is not present. It cannot be installed in the same `pip install` as
torch. The sequence that works is two steps, and needs a C++ compiler (plus the
CUDA toolkit for GPU support):

```bash
pip install torch                     # first, on its own
pip install torch-scatter --no-build-isolation
```

Note also what this does *not* establish: `fair-esm` installing is not ESM-IF
running. ESM-IF is a Python API with no command line, so this package has no
verified way to invoke it, and the `esm_if` adapter says exactly that rather
than shelling out to a command nobody has checked.

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

## 3. Container runtime (for RFdiffusion and GearBind only)

Two adapters use a container. Everything else installs into `core` or one of the
per-model environments in [provision/manifest.toml](provision/manifest.toml), so
you can skip this section entirely for a scoring-and-gating campaign.

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

## 4. Model routes and weights

### RFdiffusion — the one published image, pinned by digest

Checked on Docker Hub on 2026-10-30: `rosettacommons/rfdiffusion` exists and
has exactly one tag, `latest`, 8.1 GB, pushed 2025-06-25. The adapter pins the
digest that tag pointed at rather than the tag itself, because a tag that moves
under a campaign changes the scores without changing the config, and the cache
key would not notice. The pin lives in `design_pipeline/routes.py` as
`RFDIFFUSION_IMAGE`.

The adapter passes `inference.model_directory_path=/app/RFdiffusion/models`,
which is where the weights live inside that image. Nothing to download
separately.

```bash
docker pull rosettacommons/rfdiffusion@sha256:b88dc5d7ccd166557fb9b922bbaeb15bdf8d0568907c50405a8bf282e67254fe
# apptainer pulls on first use and caches it; to pre-pull:
apptainer pull docker://rosettacommons/rfdiffusion@sha256:b88dc5d7ccd166557fb9b922bbaeb15bdf8d0568907c50405a8bf282e67254fe
```

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

### Scorers and gate checks — routes, not images

**This package ships no container image names except RFdiffusion's.** An
earlier version of this file claimed the registry held "the published ones
where a published one exists". That was wrong: ten adapters named ten images,
and when they were checked against the Docker Hub API on 2026-10-30, eight of
the ten returned 404. The names were GitHub `org/repo` slugs used as registry
paths. A default that 404s is worse than no default, because it fails inside
`docker pull` with a registry error instead of telling you what to install.

What each one actually needs, and whether we have checked it:

| adapter | where | route | checked |
|---|---|---|---|
| `antifold` | Env B | `pip install antifold` (PyPI 0.3.1, BSD-3; weights download on first run) | yes |
| `esm_if` | `core` | `pip install fair-esm` + `torch-geometric` + `torch-scatter` — but no CLI exists to invoke | route yes, invocation no |
| `ablang2` | `core` | `pip install ablang2` | yes |
| `esm1v`, `esm2_650m` | `core` | `pip install 'design-pipeline[plm]'` | yes |
| `antiberty` | `core` | `pip install antiberty` | no |
| `consensus`, `liabilities` | — | bundled; `consensus` needs an OAS table that has never been built, so it is out of `DEFAULT_SCORERS` ([docs/provisioning.md](docs/provisioning.md)) | yes |
| `thermompnn` | Env B | `git clone https://github.com/Kuhlman-Lab/ThermoMPNN`, then `uv sync --extra cu121` against its own `pyproject.toml` (or `mamba env create -f environment.yaml`); needs `THERMOMPNN_SRC` | route yes, invocation no |
| `gearbind` | container | build from `https://github.com/DeepGraphLearning/GearBind`; no published image | yes |
| `rosetta_ddg` | manual | licence-gated. An official `rosettacommons/rosetta` image does exist (23 tags) but is not configured here | yes |
| `foldx` | manual | licence-gated; your own image | yes |
| `alphafold3` | manual | weights from DeepMind on request; your own image | yes |
| `aggrescan3d`, `netsolp`, `abmpnn`, `proteinmpnn_score` | — | **no route established.** These carried invented image names and nothing replaces them | no |

AbMPNN has no public repository and its weights are request-gated, so there is
nothing to install at all today.

An adapter whose tool is not present raises `MissingRuntimeError` quoting its
route verbatim. An adapter with no route says that instead of guessing one —
and `preflight` marks an unchecked route as such, so the table never implies
more than was actually verified.

To wire one up once you have built it, each adapter is plain data, so the image
name, the flag the candidates CSV arrives on, and the column the result is read
from are all overridable without touching this package:

```python
from dataclasses import replace
from design_pipeline import REGISTRY
from design_pipeline.adapters import _scorers

REGISTRY["gearbind"] = replace(
    _scorers.gearbind,
    tool=_scorers.matrix_container("ghcr.io/ourlab/gearbind:1.0"),
    value_column="ddg_pred",
)
```

Or per call, when it is only the column names that moved:

```python
score_variants(complex_, variants, models=["antifold"],
               id_column="mutation", value_column="llr")
```

The registry holds no image names of its own, so there is nothing here to
verify before a pull: the only image this package names is RFdiffusion's,
pinned by digest. Every other container is one you build and set yourself, by
the override above.

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
# a pair has two mutations: `antifold`'s lookup has no answer, `antifold_seq`
# scores each pair as a whole and needs the chain ids
pair_scores = score_variants(complex_, pairs, models=["antifold_seq"],
                             heavy="B", light="A")
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
