# PROTECher

Team repo for the Breast Cancer Awareness Hackathon (HerCode, part of Swiss {ai} Weeks), Zürich, 31 Oct 2026, 9:00–19:00. We are aiming at the Adaptyv protein design challenge.

Event page: https://luma.com/4ogriih0

## Team

| Name | Background |
|---|---|
| Arion Förtsch | Bioinformatics, R, cluster work |
| Kira | ML engineer (Berlin), language and graph models, protein track at the last Future of Health hackathon |
| Maria | MSc medical physics and computational science, UZH |
| Iason | MSc physics, ETH; simulations and computer vision |
| Alixe | Life sciences engineering, EPFL; biosensing and neuroscience |

## Layout

```
docs/             challenge notes and filtering strategy
design_pipeline/  backbones, sequences, folding, scoring, selection (Kira)
examples/         runnable scripts, e.g. the iterative refinement loop
filters/          scoring and filtering code
notebooks/        exploration
data/             local only, git-ignored
```

How to install the environment, the container images and the model weights, and
how to run the pipeline: [INSTALL.md](INSTALL.md).

## Structure prediction models

Four co-folding models are registered for `predict_structure`, one optional
extra each, plus a fixture counterpart per model that runs offline with no GPU.
`model_name` also takes a list, which folds every sequence with every model
named.

| `model_name` | Extra | Code licence | Weights licence | MSA | Hardware |
|---|---|---|---|---|---|
| `chai` (default) | `[chai]` | Apache 2.0 | Apache 2.0, commercial use stated | optional, off by default; `use_msa_server=True` for the public MMseqs2 server | Linux + CUDA, bfloat16 (compute capability 8.0+, i.e. Ampere or newer) |
| `boltz` | `[boltz]` | MIT | MIT | expected; `use_msa_server=True` is the default | Linux + CUDA |
| `protenix` | `[protenix]` | Apache 2.0 | **unresolved, see below** | on by default; the `protenix_mini_esm_*` / `protenix_mini_ism_*` checkpoints run MSA-free | Linux + CUDA |
| `colabfold` | `[colabfold]` | MIT | AlphaFold2 parameters, CC BY 4.0 (attribution) | public MMseqs2 server by default | GPU strongly recommended |
| `dummy_chai`, `dummy_boltz`, `dummy_protenix`, `dummy_structures` | none | — | — | no | anything |

Licences above are what each project's own bundled `LICENSE` and README say,
read on 2026-10-25: Chai-1 "is released under an Apache 2.0 License (both code
and model weights)", Boltz "All the code and weights are provided under MIT
license", Protenix as below. None of this is legal advice — check again before
anything ships.

Every adapter reports `plddt` (0-100), `ptm` and `iptm` whatever its tool calls
them, writes its prediction as a PDB with the original mmCIF in `meta["cif"]`,
and keeps the tool's own score file in `meta["raw_scores"]`. Model-specific
numbers come through prefixed, e.g. `chai_aggregate_score`,
`boltz_affinity_pred_value`.

**Protenix weights: do not assume commercial use is cleared.** For the pinned
version (`protenix==2.0.0`) the wheel's own `LICENSE` is Apache 2.0 and its
README reads "The Protenix project including both code and model parameters is
released under the Apache 2.0 License. It is free for both academic research
and commercial use." The project's GitHub page says the same for the project,
*and separately* states that the weights of the `protenix-v2` checkpoint "are
proprietary and confidential information of the rights holder" and "are not
released under any open-source license". So the licence depends on the
checkpoint you actually run, the default (`protenix_base_default_v1.0.0`) and
Mini checkpoints are covered by the project statement, and `protenix-v2` is
not. Checked 2026-10-25; re-read both before relying on either.

**ColabFold.** Kept and still working. The AlphaFold2 parameters it uses are
CC BY 4.0 per DeepMind's repository, i.e. attribution rather than
non-commercial — worth knowing if a licence audit is what pushed us to the
other three.

## Working rules

- Branch, open a PR, get one review. Keep `main` runnable.
- No data, structures over a few MB, or model weights in git. Share those through the team chat or a shared drive.
- Update `docs/challenge_notes.md` as soon as the organisers announce the target, metric and limits.

## Status

Challenge details are unknown until the event. Preparation notes are in [docs/challenge_notes.md](docs/challenge_notes.md) and [docs/filtering_strategy.md](docs/filtering_strategy.md).
