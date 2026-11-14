# Provisioning

How each environment this project needs is built, what forces it to be separate,
and what cannot be installed at all because it is licence-gated or was never
built. Written 2026-11-06.

The companion document is [cluster_gpu_node_rebuild.md](cluster_gpu_node_rebuild.md),
which records what was *measured* on `fgcz-r-023` before its scratch was
cleared. This document is the other half: what to *run* to get it back. Where
the two disagree, the rebuild doc is the observation and this one is the
intention.

Three things hold the pieces together:

| | |
|---|---|
| [provision/manifest.toml](../provision/manifest.toml) | the machine-readable inventory: per environment a Python version and a freeze file, per source install a git URL and commit, plus the weights and the prerequisites |
| `provision/*.sh` | one script per environment, written against the *constraint* that forces the separation rather than against exact pins |
| [design_pipeline/routes.py](../design_pipeline/routes.py) | what the code believes about where each tool lives. `tests/test_route_drift.py` fails CI when the two disagree |

## Inventory

**Required.** Without all four there is no final ranking.

| Env | GPU | Python | What for |
|---|---|---|---|
| `core` | no | 3.12 | this package, the tests, the sequence-only scorers, every table and plot. CPU torch, for the scorers and nothing else |
| `foldx` | no | n/a | the physics ddG baseline. A licensed binary, not a pip install |
| `boltz` | **yes** | 3.12 | Boltz-2. The folds the final ranking's ipSAE and interface metrics are computed from |
| `antifold` | **yes** | 3.13 | AntiFold, for the evidence funnel's veto, plus AbLang2 in the same environment |

**Optional.** Each buys one thing and costs an environment; none is on the path
to a ranking.

| Env | GPU | Python | What for |
|---|---|---|---|
| `esm_env` | yes | 3.12 | ESM-1v and ESM-IF1 — the scorer calibration, not the screen |
| `thermompnn` | yes | 3.10 | the stability gate |
| `dockq` | no | 3.12 | DockQ, for the ground-truth half of `design_pipeline.metrics` |
| `chai` | yes | 3.12 | Chai-1, as a second opinion against Boltz-2. See the caveat below |

### The default scorer is in `core`, and has to be

`DEFAULT_SCORERS` is `("ablang2",)`, so `ablang2` is in preflight's
`DEFAULT_CAMPAIGN` and decides its exit status. [routes.py](../design_pipeline/routes.py)
puts it in `core` and declares no `env_var` — `tests/test_route_drift.py` forbids
one for a `core` route — and `AbLang2MaskedLM._load` imports `ablang2` in
process. Those three facts together mean `core` has to hold it: there is no path
by which it can be reached elsewhere, the way `antifold` and `boltz` are.

This document used to say `core` had "no torch at all", and
`provision/core.sh` asserted it with `uv pip list | grep torch` while
`provision/antifold.sh` installed `ablang2` into the *AntiFold* environment. The
scorer was therefore in a venv nothing reached it from, and preflight in a
freshly provisioned `core` exited 1 on the only row it gates on. INSTALL.md's
`core` row — "this package, `ablang2`, `antiberty`, `transformers` + `torch`,
`fair-esm`, `anarci`" — was right all along; the script and this file were the
drift.

`core.sh` now installs the `plm` extra with `--torch-backend=cpu`, pinning
`ablang2==0.2.1` to the version [freeze_antifold.txt](cluster_freezes/freeze_antifold.txt)
measured so the two copies cannot diverge, and asserts the constraint as what it
actually is: a fresh interpreter importing `design_pipeline` must not pull in
`torch` or `jax`. That is a property of the lazy imports in
`adapters/_plm.py`, not of the package list, and `tests/test_steps.py` asks it
the same way.

`antifold.sh` keeps its own `ablang2`, because
[freeze_antifold.txt](cluster_freezes/freeze_antifold.txt) records it and the
freeze is the measurement. `core` is where a campaign reaches it.

**Dropped.** `protenix` and the LigandMPNN environment are not provisioned.
Both served de novo design, which this project is not doing. Their adapters,
routes and extras stay in the package — `routes.py` describes what the code can
reach, `manifest.toml` describes what we build, and those are different
questions — but no script here builds them, and nothing in the required path
calls them.

Chai-1 is kept in the manifest, with its freeze, as a drop-in alternative to
Boltz-2. Nothing requires it. Note the asymmetry: the Boltz adapter shells out
to the `boltz` binary, so it can be driven from `core` by path, while the Chai
adapter calls `run_inference` **in process** (the CLI drops the PAE matrix, and
`pae_interaction` needs it). Chai therefore has to be driven from inside its own
environment; `DESIGN_PIPELINE_CHAI` is a preflight probe only and does not make
the adapter reachable from `core`. That is the one place where the "same variable
the adapter reads" rule in `routes.Route` cannot hold today.

## What forces each separation

Written as constraints, because a constraint survives a version bump and a pin
does not. The scripts below install against these, not against the freezes.

| Env | The conflict | Why `core` cannot absorb it |
|---|---|---|
| `boltz` | `boltz[cuda]` brings torch 2.15 plus a mixed CUDA 12/13 `nvidia-*` wheel set | `core`'s torch is the CPU build the sequence-only scorers run on; this is a different resolve, and importing `design_pipeline` must still import no model |
| `antifold` | needs **Python >= 3.13**; `core` is 3.12 | a different interpreter is a different environment, full stop |
| `antifold` (PyPI route, not used) | `pip install antifold` pins `biotite==0.38` against this package's `biotite>=2.0` | `ResolutionImpossible`, recorded verbatim in [INSTALL.md](../INSTALL.md) |
| `thermompnn` | its `environment.yaml` is Python 3.10 with `pytorch-cuda 11.7` | two majors of CUDA behind everything else here |
| `chai` | `chai_lab` pins torch 2.6 against boltz's 2.15 | the two folding models pin torch against each other, so they cannot share |
| `esm_env` | `transformers` 5.x plus its own CUDA torch | `core`'s torch is the CPU build; ESM-1v on a GPU is its own resolve |
| `dockq` | brings biopython, which nothing in this package imports | would be dead weight in `core`, and DockQ is wrapped not reimplemented |
| `foldx` | a licensed binary | not a Python package at all |
| `gearbind` (container, not an env) | TorchDrug pins Python <= 3.10 and an old torch | the one conflict a second *environment* cannot fix, since nothing in it can share a process with us |

### A note on `biotite==0.38`

The constraint that *forces* AntiFold out of `core` is `biotite==0.38`, which
every published PyPI release of AntiFold pins. The install we actually use is
**not** that one: it is the GitHub source at 0.4.0 or later, which drops the
biotite pin and sits happily next to `biotite 2.0.0` — which is what
[freeze_antifold.txt](cluster_freezes/freeze_antifold.txt) measures (biotite
2.0.0, numpy 2.5.4, Python 3.13.14).

So both facts are true and they are about different routes. `biotite==0.38` is
recorded in `manifest.toml` as the separation constraint, under the PyPI route
it belongs to. The source route's separation driver is Python 3.13. Installing
`biotite==0.38` into the environment we actually build would *break* it.

## The CUDA variant is a parameter

No freeze file in this repository carries a `+cuXXX` local-version tag, and none
records an `--index-url`. CUDA arrives entirely as ordinary `nvidia-*` wheel
dependencies of torch, which means a replay resolves torch from whatever index
is default and the CUDA build is implicit.

Every GPU script therefore takes the variant explicitly:

```bash
provision/boltz.sh --cuda cu128
```

`--cuda` accepts `auto` (the default; uv detects the driver), `cpu`, or a
`cuXXX` moniker, and is passed to uv as `--torch-backend`. The CPU-only
environments (`core`, `dockq`, `foldx`) accept the flag and ignore it, so a
driver script can pass one value to all of them.

### The Boltz exception

`provision/boltz.sh` does **not** install against constraints. It replays
[freeze_venv-boltz.txt](cluster_freezes/freeze_venv-boltz.txt) verbatim, because
that environment carries a mixed CUDA wheel set — `nvidia-cublas==13.4.1.3`
alongside `nvidia-cublas-cu12==12.9.2.10`, and the same doubling for
`nvidia-cuda-nvrtc` — almost certainly `cuequivariance-ops-*` pulling cu12 next
to torch 2.15's cu13. It worked on the L40S. A fresh resolve is not guaranteed
to reproduce that pairing, and this is the environment whose numbers went into
the final ranking, so it is the one place where the exact set is worth more than
a readable constraint.

`--cuda` is accepted and refused with a message on this script, rather than
silently ignored: a replay has its CUDA variant baked in.

The freezes are `uv pip freeze`, not `pip freeze` — every one begins with uv's
`Using Python ... environment at: ...` banner, which pip never emits. Three
further lines are not installable: `-e file:///scratch/foertsch/protecher` in
`freeze_venv.txt` and `freeze_venv-chai.txt`, and
`antifold @ file:///scratch/.../antifold/src` in `freeze_antifold.txt`. The
filter recorded in [cluster_gpu_node_rebuild.md](cluster_gpu_node_rebuild.md) —

```bash
grep -vE '^(-e|Using)|@ file' docs/cluster_freezes/freeze_venv-boltz.txt
```

— handles the banner and all three, and is what the scripts use.

## FoldX

FoldX is a licensed download and nothing here can fetch it. What the script does
is validate an unpacked copy, because every FoldX failure we have had was a
layout problem rather than an install problem.

`provision/foldx.sh --zip ~/Downloads/foldx5_Linux_0.zip` unpacks into
`tools/foldx/` (gitignored) so that the binary and `molecules/` end up
**siblings**, then checks three things:

1. the binary exists and is executable (`chmod +x` if not),
2. `molecules/` is beside it, not nested a level down — FoldX reads its rotamer
   data by relative path and fails obscurely when it is not,
3. `--version` runs and exits 0.

It prints the `export FOLDX_BINARY=...` line to add to the environment, and the
two further variables the adapter reads (`DESIGN_PIPELINE_FOLDX_PDB`,
`DESIGN_PIPELINE_FOLDX`).

### The expiry warning

Academic FoldX licences are time-limited and the build is stamped with its
expiry in the filename: the copy in `~/tools/foldx/` is `foldx_20261231`, i.e.
31 December 2026. `design_pipeline.preflight` now parses that date out of
`FOLDX_BINARY` and warns when it is inside 30 days or already past. It is a
warning, not a failure: the date is a filename convention, not something the
binary reports, so a renamed binary must not break preflight. A rename simply
means no warning.

## Prerequisites that no script can satisfy

| Prerequisite | State | What it blocks |
|---|---|---|
| FoldX licensed download (`foldx5_Linux_0.zip`) | **held.** `~/tools/foldx/` on the cluster's shared home, which is the only path that survived the scratch clear | the `foldx` scorer |
| OAS frequency table (`DESIGN_PIPELINE_OAS_TABLE`) | **does not exist, and never did** | the `consensus` scorer |

The OAS table was never committed (`git log --all --diff-filter=A` over the
whole history matches nothing), is absent from the node inventory's survive and
scratch tables, and `DESIGN_PIPELINE_OAS_TABLE` is not among the variables
[cluster_gpu_node_rebuild.md](cluster_gpu_node_rebuild.md) records as having
been set. Nor does `consensus` appear in
[data/track_a_final/ranking.csv](../data/track_a_final/ranking.csv). Nothing was
lost; it was simply never built.

Accordingly `consensus` is **out of `DEFAULT_SCORERS`**, which is now
`("ablang2",)`. The adapter, the route and the registry entry are unchanged, so
`models=["consensus"]` still works the moment a table exists.

### Open task: build the OAS table

Not done, and it is its own piece of work:

1. download a paired-sequence slice from the Observed Antibody Space
   (`https://opig.stats.ox.ac.uk/webapps/oas/`),
2. number it to IMGT with ANARCI or anarcii,
3. count per `(imgt, aa)` and normalise per position,
4. write `imgt,aa,frequency` as CSV.

The reader is `_frequency_table` in
[design_pipeline/adapters/_scorers.py](../design_pipeline/adapters/_scorers.py):
it keys on `imgt` when that column is present and falls back to `position`,
which only means anything for the antibody the table was built from. The scorer
is single-mutant only by design — combining two positional marginals is the
additivity shortcut under another name, so a double mutant gets an empty cell.

## Source installs with unknown commits

Two environments install from a git tree and **neither commit is known**.

| Env | URL | Commit |
|---|---|---|
| `antifold` | `https://github.com/oxpig/AntiFold` | `UNKNOWN` |
| `thermompnn` | `https://github.com/Kuhlman-Lab/ThermoMPNN` | `UNKNOWN` |

The node's copies had no `.git` directory, so neither origin nor commit was
recorded before the scratch was cleared. The existing results came from an
unidentified tree in both cases. This is recorded as `commit = "UNKNOWN"` in
`manifest.toml` rather than guessed, and the scripts take `--commit` so that the
next build can be pinned even though the last one cannot be reproduced.

For AntiFold the version is partly constrained — 0.4.0 or later, per
[INSTALL.md](../INSTALL.md), which is what the Python 3.13 requirement and the
absent biotite pin both point at. For ThermoMPNN nothing narrows it beyond "a
tree whose `environment.yaml` is Python 3.10 with `pytorch-cuda 11.7`".

## Weights

No script downloads these; each tool fetches its own on first use, into the
location its variable names. Sizes measured on the node.

| Artifact | Size | Env var | Source |
|---|---|---|---|
| Boltz-2 (`boltz2_conf.ckpt`, `boltz2_aff.ckpt`, `mols.tar`) | 6.3 GB | `BOLTZ_CACHE` | first `boltz predict` |
| Chai-1 (`models_v2/`, traced ESM-2, conformers) | 6.9 GB | `CHAI_DOWNLOADS_DIR` | first `run_inference` |
| AntiFold `model.pt` | 567 MB | `ANTIFOLD_SRC` (`src/models/`) | first use; `weights=` reuses a download |
| AbLang2 paired weights | 179 MB | — (inside whichever venv imports it) | first use |
| ESM-IF1 `esm_if1_gvp4_t16_142M_UR50.pt` | 142 M params | `TORCH_HOME` | `torch.hub` |
| ESM-1v (`facebook/esm1v_t33_650M_UR90S_1`) | 650 M params | `HF_HOME` | Hugging Face |

`BOLTZ_CACHE` must be absolute — that is Boltz's own requirement, and the
adapter raises on a relative value rather than letting it fail inside the tool.

## The route contract, and the test that enforces it

A tool in another environment cannot be imported from this one. So for any route
that does not live in `core`, the only probe that can succeed from outside it is
"has the path been configured" — `importlib.util.find_spec` asks *this*
interpreter about a package that is deliberately not in it and always answers
no. The reasoning was already written in
[routes.py](../design_pipeline/routes.py) for AntiFold; it now applies to every
such route.

Nine routes previously declared `module` or `executable` only and so were
probed, wrongly, against `core`'s interpreter: `chai`, `chai1`, `boltz`,
`boltz2`, `protenix`, `colabfold`, `protein_mpnn`, `ligand_mpnn`,
`soluble_mpnn`. Concretely, `boltz` ran the entire final screen and preflight on
the node would still have printed `import boltz failed`. Each now declares an
`env_var`:

| Variable | Points at | Routes |
|---|---|---|
| `DESIGN_PIPELINE_ANTIFOLD` | the AntiFold environment's `python` | `antifold`, `antifold_seq` |
| `DESIGN_PIPELINE_THERMOMPNN` | the ThermoMPNN environment's entry point | `thermompnn` |
| `DESIGN_PIPELINE_BOLTZ` | that environment's `boltz` binary | `boltz`, `boltz2` |
| `DESIGN_PIPELINE_CHAI` | the Chai environment's `python` (probe only — see above) | `chai`, `chai1` |
| `DESIGN_PIPELINE_PROTENIX` | that environment's `protenix` binary | `protenix` |
| `DESIGN_PIPELINE_COLABFOLD` | that environment's `colabfold_batch` | `colabfold` |
| `DESIGN_PIPELINE_MPNN` | that environment's `ligandmpnn` | `protein_mpnn`, `ligand_mpnn`, `soluble_mpnn` |
| `DESIGN_PIPELINE_FOLDX` | the python that runs `adapters/_foldx_score.py` | `foldx` |

`tests/test_route_drift.py` is what keeps this from rotting: every route whose
environment is not `core` and not bundled must declare an `env_var`, and setting
that variable must make `resolve()` return true and cite it. It fails on the
pre-change tree for exactly those nine routes.

`env` and `env_var` answer different questions, which is why `foldx` can be
`env=MANUAL` and still carry one: `env` says where the tool *comes from*
(licence-gated, you obtain it yourself), `env_var` says how it is *reached* (by
path, like Env B). `Route.env_var`'s docstring said "for an Env B route", which
was the inaccuracy; it now says "reached by path".

`foldx`'s `needs_env` also now lists all three variables the adapter actually
reads — `FOLDX_BINARY`, `DESIGN_PIPELINE_FOLDX_PDB` and `DESIGN_PIPELINE_FOLDX`.
It listed only the first, so preflight never reported the other two as unset
even though [_foldx_score.py](../design_pipeline/adapters/_foldx_score.py) reads
`DESIGN_PIPELINE_FOLDX_PDB`.

## Container runtime

`detect_runtime()` tried `apptainer` then `docker` and took the first on PATH,
with no way to say which you meant. On a cluster apptainer is the only one
allowed, so preferring it is right there; on a workstation with both installed
you got apptainer silently.

`DESIGN_PIPELINE_RUNTIME` now overrides it:

```bash
DESIGN_PIPELINE_RUNTIME=docker python -m design_pipeline.preflight
```

Set it to `apptainer` or `docker`. An unrecognised value is an error naming the
two that work, and a named runtime that is not on PATH is an error too rather
than a silent fall back to the other one — if you asked for docker, getting
apptainer is not a favour.

This does not make GearBind reachable. `gearbind` is a `ContainerScorer` with
`tool=None` and its route declares no image, so it resolves to "not configured"
whatever the runtime. Docker makes the *build* easy — TorchDrug on Python 3.10
is exactly what a container is for — but you still have to build from
`https://github.com/DeepGraphLearning/GearBind` and set `REGISTRY['gearbind']`.

## Housekeeping: the per-job repo copies

The node's scratch held roughly 60 `protecher_*` job directories of 1 to 12 MB,
each a copy of the whole repo tree plus a log, listed as "not worth keeping" in
[cluster_gpu_node_rebuild.md](cluster_gpu_node_rebuild.md). Together with
`uv-cache` (28 GB) and `apptainer/cache` (15 GB) they were the bulk of the
reclaimable space.

`provision/clean_scratch.sh` removes them. It is **dry-run by default** and
prints what it would delete; `--apply` performs it. It refuses to touch the
directories the rebuild doc lists as worth keeping, matches only the per-job
name pattern, and will not run against a path that is not a scratch root.

Note that the node's `/scratch/foertsch` was already cleared on 2026-11-02, so
on today's node there is nothing for it to find. The script exists for the next
node, and nothing in this repository can reach the cluster anyway.

## Running it

```bash
provision/core.sh
```

```bash
provision/foldx.sh --zip ~/Downloads/foldx5_Linux_0.zip
```

```bash
provision/boltz.sh
```

```bash
provision/antifold.sh --cuda cu128
```

Then, from `core`, with the variables from the scripts' closing `export` lines in
the environment:

```bash
python -m design_pipeline.preflight
```

Every required row should read `yes`. A row marked `B` that does not resolve
inside `core` means its `DESIGN_PIPELINE_*` variable is unset, not that the tool
is missing.
