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
design_pipeline/  the pipeline: structure, variants, scorers, gate, harness (Kira)
examples/         runnable scripts, e.g. examples/maturation.py
filters/          scoring and filtering code
notebooks/        exploration
data/             local only, git-ignored
```

How to install the environment, the container images and the model weights, and
how to run the pipeline: [INSTALL.md](INSTALL.md).

## The pipeline

Antibody affinity maturation, not de novo design: it starts from a molecule
that already exists. Five functions, everything on disk, no DAG and no
scheduler.

```python
from design_pipeline import (
    get_structure, enumerate_variants, score_variants, gate, select_pairs,
)

complex_ = get_structure("runs/her2_v1", pdb_id="1N8Z", chains=["A", "B", "C"])
variants = enumerate_variants(complex_, chain="B")              # 4,598 rows
scores = score_variants(complex_, variants,
                        models=["antifold", "esm_if", "consensus"])
survivors = gate(scores, checks=["liabilities"])
pairs = select_pairs(survivors, max_pairs=20, min_distance=10.0)
pair_scores = score_variants(complex_, pairs, models=["antifold", "gearbind"])
```

Run it offline, with no GPU, no Docker and no network:

```bash
python examples/maturation.py --dummy
```

On the GPU box, starting from sequences rather than an accession, in a notebook:
[notebooks/03_end_to_end.ipynb](notebooks/03_end_to_end.ipynb) takes the antibody
and antigen FASTA, folds the parent over several seeds, runs the scan, the
scorers, the gate and the pairs, refolds the designs as a structure veto and
writes the K submission sequences. It has no fixture mode — it runs `preflight`
on the adapters it is about to use and stops if one of them is missing.

Three things about this that are deliberate:

- **Scoring is a fan-out and a join, not a chain.** One structure feeds several
  models in parallel; each writes one `score_<model>` column keyed on
  `variant_id`; the columns are joined into one table. The cache is per scorer,
  so adding GearBind does not invalidate the AntiFold run that already
  finished.
- **Nothing is averaged into a combined score.** No weights, no learned ranker.
  Ranking is the caller's business, and the right weights come out of the
  harness below rather than out of a guess in the library.
- **A pair is an ordinary candidate.** `score_variants` is called again on the
  pairs, and no code path anywhere estimates a pair's score by adding its
  singles' scores together. `select_pairs` drops the score columns for exactly
  that reason.

### Scorers

Two cost classes, and the API makes the split explicit rather than hiding it.
A matrix scorer runs once on the complex and answers every single mutation from
one forward pass. A per-variant scorer builds and relaxes a mutant structure
per candidate, so `score_variants(..., max_variants=n)` caps it and an uncapped
one refuses a table of more than a couple of thousand rows.

| `models=[...]` | Cost | What it is |
|---|---|---|
| `antifold` | matrix | antibody-specific inverse folding; out of the default set since 2026-10-31, see [CHANGELOG.md](CHANGELOG.md) |
| `antifold_seq` | per variant | the same model over a whole mutated chain, so a combination is scored as one mutant |
| `esm_if` | matrix | general inverse folding |
| `abmpnn`, `proteinmpnn_score` | matrix | same shape, lower priority |
| `consensus` | matrix | positional consensus from a cached OAS table; no model, no GPU |
| `gearbind` | per variant | supervised ddG, the one with wet-lab validation |
| `rosetta_ddg` | per variant | Rosetta Flex ddG; slow, errors uncorrelated with the learned models |
| `foldx` | per variant | cheaper physics baseline, licence permitting |
| `dummy_scorer`, `dummy_per_variant` | — | deterministic fixtures, offline |

With no `models=`, `score_variants` runs `ablang2` and `consensus`
(`DEFAULT_SCORERS`), neither of which needs a structure-model environment.

`score_<model>` is normalized at the adapter boundary so higher always means
better binding; `score_<model>_raw` keeps the tool's own number, so a
convention we got backwards is recoverable without re-running anything.

### Gate

Each check answers a boolean and a reason. Failures are dropped, never traded
off against binding score, and the annotated table keeps every failure and its
reason so a gate that dropped too much can be re-tuned without re-running the
tools. Gate only on confident failures: the developability predictors are weak,
and the reliable signal is a hard sequence liability.

| `checks=[...]` | Cost | What it looks for |
|---|---|---|
| `liabilities` | free | new N-glycosylation sequon (N-x-S/T), deamidation (NG/NS/NT), isomerization (DG/DS/DD), fragmentation (DP), introduced cysteine |
| `netsolp` | sequence | predicted solubility |
| `thermompnn`, `aggrescan3d` | per variant | folding stability, aggregation |

### Validation harness

The single most important part. On the day there is no feedback, so the only
way to know which scorer to trust is to have measured it beforehand on systems
where the answers are published.

```python
from design_pipeline.harness import run_systems, report

results = run_systems(
    systems=["1n8z_her2", "4fqi_flu"],
    models=["antifold", "esm_if", "gearbind", "consensus"],
    output_dir="runs/harness",
)
print(report(results))        # Spearman and precision@10, per scorer per system
```

This is the ordinary pipeline in a `for` loop with a correlation on the end,
not a second code path. Ground truth comes from
[FLAb](https://github.com/Graylab/FLAb), AbBiBench and
[AB-Bind](https://github.com/sarahsirin/AB-Bind-Database); none of it is
redistributable, so the real systems read it from `data/ground_truth/` and two
fixture systems ship with the package so the harness runs in CI against the
dummy scorer. `report` breaks results out per scorer, never for an ensemble
only — the ensemble number hides the thing you came to learn.

## Structure prediction models

`get_structure` downloads by accession by default, because experimental
coordinates beat a prediction every time. When no solved complex exists, these
predict one from `antibody_seq` and `antigen_seq` instead — one optional extra
each, plus a fixture counterpart per model that runs offline with no GPU.
`alphafold3` takes `num_seeds` and writes every seed onto `complex_.paths`,
because the ranking of antibody-antigen complexes is not reliable enough to
hide the rest of the ensemble.

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

RFdiffusion and the MPNN variants are kept, demoted to the optional branch
behind `diffuse_loops(complex_, cdr="H3")`, which rebuilds one CDR and returns
an ordinary variant table. Those variants differ from the parent at ten or
twenty positions at once, far outside the edit-distance range where any of the
scorers was validated.

Every folding adapter reports `plddt` (0-100), `ptm` and `iptm` whatever its
tool calls them, writes its prediction as a PDB with the original mmCIF in `meta["cif"]`,
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

Challenge details are unknown until the event. Preparation notes are in [docs/challenge_notes.md](docs/challenge_notes.md) and [docs/filtering_strategy.md](docs/filtering_strategy.md). The chain that produced the submitted Track A designs is in [docs/track_a_final_pipeline.md](docs/track_a_final_pipeline.md).
