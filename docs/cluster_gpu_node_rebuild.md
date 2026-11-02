# GPU node (fgcz-r-023) setup, and how to rebuild it

Written 2026-11-02, before the node's `/scratch/foertsch` was cleared. It records what was installed for the Track A work, with the versions measured on the node that day, and a rebuild order. The exact package lists (`uv pip freeze`) are in [cluster_freezes/](cluster_freezes/); the cluster's shared home also keeps them, with the layout notes, in `~/protecher_inventory/`, which is not removed with the node's scratch.

What this is not: a record of the commands first typed. The versions are measured; the install commands are reconstructed from them and from [INSTALL.md](../INSTALL.md), and the sections below say where that matters.

## The node

| | |
|---|---|
| Host | `fgcz-r-023`, 2 x NVIDIA L40S (driver 590.48.01), Debian 12, kernel 6.12.57 |
| Local disk | `/scratch`, 7.0 TB, **node-local** (the login node `fgcz-c-043` has its own `/scratch`) |
| Tools on the node | apptainer 1.4.2, uv 0.11.21 (`/usr/bin/uv`) |
| Access | login `ssh -o HostKeyAlias=fgcz-c-043 foertsch@130.60.193.43`; jobs by `sbatch -p GPU_L40S -w fgcz-r-023` |
| Shared across nodes | home (`/home/foertsch`) only |

Everything below lives under `/scratch/foertsch/` unless it says `~`.

## What survives clearing the node

| Item | Where | Notes |
|---|---|---|
| Code | GitHub `foertsch/PROTECher`, main at `9fe5c1c` or later | the node's `protecher/` was a clone at `01d85d2` |
| FoldX | `~/tools/foldx/` | `foldx_20261231` (87 MB binary), `foldx5_Linux_0.zip`, `molecules/`. The zip is the licensed download; keep it |
| Parent alignments | `~/protecher_msa_candidates/` | `trastuzumab_her2/parent_0.csv`, `parent_1.csv` (4.1 and 3.7 MB, same size as the node's), the other 5 antibodies, `structures/` (1N8Z, 1S78, 1YY9, 5GGS, 5X8L, 9PI9 mmCIF) |
| Package lists | `docs/cluster_freezes/`, and `~/protecher_inventory/` | exact `uv pip freeze` of every environment |
| Scores | the repo, `data/track_a_final/` | rows, ranking and picks of the final screen |
| Modal | volumes `protecher-work`, `boltz-weights` | fold outputs of the Modal slices; weights volume used by `cloud/modal_fold.py` |

## What exists only on the node scratch (copy first if wanted)

| Path | Size | What |
|---|---|---|
| `protecher_screen2`, `screen3`, `screen4_s1`, `screen4_s2` and the small `protecher_screen*` job dirs | about 1.0 GB | raw Boltz-2 folds of the cluster slices (mmCIF, PAE `.npz`, confidence JSON). The scores are in the repo; the structures are not |
| `protecher/runs/` | 612 MB | earlier pipeline runs and smoke tests |
| `protecher_track_a_run`, `protecher_loop_arm` | 141 MB, 168 MB | evidence-funnel run output; the exploratory RFdiffusion loop arm |
| `protecher_pr23_foldx_results`, `protecher_foldx_calib`, `protecher_results`, `esm1v_scan` | under 3 MB each | FoldX run, calibration and ESM-1v scan tables |
| `protecher_tools/ipsae/` | 0.1 MB + `work/` | `ipsae.py` (Dunbrack **v4**) and `ipsae_reference_results.csv`, the reference values the in-process scorer was checked against on 30 Boltz-2 complexes. `design_pipeline/ipsae.py` in the repo is v3 |
| `protecher_cache/msa_candidates/trastuzumab_her2/variants/` | part of 548 MB | per-variant alignments derived from the parent's; regenerable by `msa_cache` |

Not worth keeping: `protecher_tools/uv-cache` (28 GB), `apptainer/cache` (15 GB), `~60` small per-job `protecher_*` dirs of 1 to 12 MB (each is a copy of the repo tree plus a log), `protecher_cache/_*_test` dirs, `.cache/pip`, `.cache/pixi`.

## Components

`T` = `/scratch/foertsch/protecher_tools`, `C` = `/scratch/foertsch/protecher_cache`.

| Component | Path | Size | Versions (measured) | Weights, source |
|---|---|---|---|---|
| **Env A** (this package, tests, scoring) | `protecher/.venv` | 0.6 GB | Python 3.12.16; biotite 2.0.0, numpy 2.5.4, pandas 3.1.0, scipy 1.18.1, pyarrow 26.0.0, matplotlib 3.11.2, pytest 9.2.0; the package editable | none |
| **Boltz-2** | `protecher/.venv-boltz`, link `T/bin/boltz` | 7.1 GB | Python 3.12.16; boltz 2.2.1, torch 2.15.0, cuequivariance 0.13.0 (+ `-torch`, `-ops-*`), numpy 1.26.4, pandas 3.0.6, rdkit 2026.9.1 | `C/boltz`: `boltz2_conf.ckpt` 2.3 GB, `boltz2_aff.ckpt` 2.1 GB, `mols.tar` 1.9 GB. `boltz predict` downloads what is missing on first run ([INSTALL.md](../INSTALL.md)); `BOLTZ_CACHE` points here |
| **Chai-1** (comparison only, not used for the final screen) | `protecher/.venv-chai` | 6.2 GB | Python 3.12.16; chai-lab 0.6.1, torch 2.6.0, numpy 1.26.4, pandas 2.3.3 | `C/chai`: `models_v2/` (trunk 633 MB, diffusion 477 MB, small heads), `esm/traced_sdpa_esm2_t36_3B_UR50D_fp16.pt` 5.7 GB, `conformers_v1.apkl` 125 MB; `CHAI_DOWNLOADS_DIR` |
| **AntiFold + AbLang2** | `T/antifold/{src,venv}` | 6.8 GB | Python 3.13; AntiFold from the GitHub source (v0.4.0 per INSTALL.md; the node copy has no `.git`), torch 2.15.0, torch-geometric 2.8.1, biotite 2.0.0, numpy 2.5.4, pandas 3.1.0, ablang2 0.2.1, anarcii 2.0.8 | `src/models/model.pt` 567 MB (AntiFold, downloaded on first use; INSTALL.md says pass `weights=`); AbLang2 paired weights 179 MB inside the venv's `ablang2/model-weights-ablang2-paired/` |
| **AbLang2 as a package** | `T/ablang2_pkgs` | 178 MB | ablang2 0.2.1, beartype 0.23.1, einops 0.8.2, rotary_embedding_torch 0.9.1 (a `pip --target` dir) | as above |
| **ESM-1v and ESM-IF** | `T/esm_env` | 5.4 GB | Python 3.12.16; torch 2.15.0, transformers 5.22.1, huggingface-hub 2.4.0, numpy 2.5.4 | ESM-1v: `C/hf` (`facebook/esm1v_t33_650M_UR90S_1`, `HF_HOME`). ESM-IF1: `T/torch/hub/checkpoints/esm_if1_gvp4_t16_142M_UR50.pt` |
| **LigandMPNN** | `T/mpnn-env` | 5.2 GB | Python 3.11.18; ligandmpnn 0.1.2, torch 2.2.1, numpy 1.23.5, biopython 1.79, prody 2.4.1, scipy 1.12.0, ml-collections 0.1.1 | parameters bundled in the pip package ([INSTALL.md](../INSTALL.md)) |
| **ThermoMPNN** | `T/thermompnn/src` | 5.5 GB | source tree with `environment.yaml` (Python 3.10, pytorch-cuda 11.7); no `.git`, so the origin and commit were not recorded | data and checkpoints inside the tree |
| **DockQ** | `T/dockq` | 93 MB | Python 3.12.13; dockq 2.1.3, biopython 1.88, numpy 1.26.4 | none |
| **Notebook runner** | `T/nbrun_pkgs` | 83 MB | ipykernel 7.4.0, jupyter_client 8.10.0, nbclient 0.11.0 (a `pip --target` dir) | none |
| **RFdiffusion** | `apptainer/images/rfdiffusion.sif` | 7.8 GB | image `rosettacommons/rfdiffusion@sha256:b88dc5d7ccd166557fb9b922bbaeb15bdf8d0568907c50405a8bf282e67254fe`, apptainer 1.4.2, built 2026-10-31 | inside the image (`/app/RFdiffusion/models`) |
| Test images | `apptainer/images/alpine.sif`, `cuda-base.sif` | 4 MB, 90 MB | `alpine`; `nvidia/cuda:12.4.1-base-ubuntu22.04` | none; made to test the runtime |
| **uv, Pythons** | `T/bin/{uv,uvx}`, `T/python/` (CPython 3.10.22, 3.11.18, 3.12.16), `T/uv-python` (3.12.13) | 0.6 GB | uv in `T/bin` is the 2026-10-27 download; the system uv is 0.11.21 | |
| **Structures** | `C/structures` | 11 MB | 1N8Z, 1S78, 1YY9, 5GGS, 5X8L, 9PI9 mmCIF | RCSB |
| **MSA cache** | `C/msa_candidates` | 548 MB | 6 antibodies; trastuzumab/HER2 is the one used | made with `msa_cache` against `https://api.colabfold.com` (a remote server call, the only one) |

The shared environment variables were set by `T/env.sh`:

```bash
export PATH=/scratch/foertsch/protecher_tools/bin:$PATH
export UV_CACHE_DIR=/scratch/foertsch/protecher_tools/uv-cache
export UV_PYTHON_INSTALL_DIR=/scratch/foertsch/protecher_tools/python
export BOLTZ_CACHE=/scratch/foertsch/protecher_cache/boltz
export CHAI_DOWNLOADS_DIR=/scratch/foertsch/protecher_cache/chai
export HF_HOME=/scratch/foertsch/protecher_cache/hf
```

Further variables the code reads: `DESIGN_PIPELINE_ANTIFOLD` (the AntiFold venv's python), `ANTIFOLD_SRC` (`T/antifold/src`), `FOLDX_BINARY` (`~/tools/foldx/foldx_20261231`), `DESIGN_PIPELINE_PDB_CACHE`, `APPTAINER_CACHEDIR`. Installs and caches stay under `/scratch/foertsch`, not in `$HOME`.

## Rebuild order

Do this on `fgcz-r-023` through `sbatch`, since `/scratch` is node-local; a login-node install does not exist for the GPU node.

1. **Folders and `env.sh`.** Recreate `T` and `C`, write `T/env.sh` as above, `source` it.
2. **Repo.** `git clone git@github.com:foertsch/PROTECher.git /scratch/foertsch/protecher`.
3. **Env A.** `uv venv --python 3.12 .venv && uv pip install -e .` plus `pytest`; compare with `cluster_freezes/freeze_venv.txt`. Check with `pytest tests -q`; the suite was 771 passed, 1 skipped on 2026-11-02.
4. **Boltz-2.** `uv venv --python 3.12 .venv-boltz && uv pip install -e '.[boltz]'` (the repo's extra), or `grep -vE '^(-e|Using)|@ file' docs/cluster_freezes/freeze_venv-boltz.txt | uv pip install -r -` (the filter drops the editable-install and header lines the freeze files carry) to get the measured versions. Link `T/bin/boltz`. Fetch the weights by running one `boltz predict` with `BOLTZ_CACHE` set. The `boltz-weights` Modal volume probably holds a copy; its contents were not checked.
5. **Parent alignments.** `cp -r ~/protecher_msa_candidates/* $C/msa_candidates/` and the structures to `$C/structures/`. Nothing else is needed to fold the parent and its variants, because the variant alignments are derived from the parent's (`python -m msa_cache`); no server call is made after that.
6. **AntiFold and AbLang2** (only to rerun the evidence funnel's veto). Python 3.13 venv, install AntiFold from the GitHub source at v0.4.0, `uv pip install ablang2==0.2.1`; compare with `freeze_antifold.txt`. Export `DESIGN_PIPELINE_ANTIFOLD` and `ANTIFOLD_SRC`.
7. **FoldX.** Nothing to rebuild; set `FOLDX_BINARY=~/tools/foldx/foldx_20261231`. The licence is the zip in `~/tools/foldx`.
8. **Optional**, only if those models are needed again: Chai-1 (`.venv-chai`, weights about 7 GB on first run), ESM environment (`freeze_esm_env.txt`), LigandMPNN (`freeze_mpnn-env.txt`), ThermoMPNN (source origin not recorded; find it by its `environment.yaml`), DockQ (`dockq==2.1.3`), RFdiffusion (`APPTAINER_CACHEDIR=/scratch/$USER/apptainer apptainer pull docker://rosettacommons/rfdiffusion@sha256:b88dc5d7...`, 8 GB).
9. **Smoke test.** `python examples/track_a_screen.py fold` on 2 candidates at 1 seed; the screen ran at 10 to 19 s a fold on this GPU.

Folding the final designs again needs only steps 1 to 5. Steps 6 to 8 are for the earlier funnel and the scorer calibration.

## Not measured or not known

- The commands originally run for each environment. The freezes are exact, the commands above are reconstructed.
- The origin and commit of the node's ThermoMPNN tree and of the AntiFold source (neither has `.git`).
- The Modal volume sizes.
- Whether each model's download URL is unchanged; weights were fetched on 2026-10-30 and 2026-10-31.
