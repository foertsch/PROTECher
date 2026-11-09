# Track A final pipeline: what was run for the 15 submitted designs

Status: written on 2026-11-02, after the submission of 2026-10-31. It records the chain that produced the 15 scFv designs, so it can be rerun. The challenge judged the pooled Boltz-2 ipSAE of the designs, so the last steps rank on ipSAE. The earlier plan is in [track_a_filtering.md](track_a_filtering.md); this is what replaced its folding and ranking stages.

## Chain

| # | Step | Script | Input | Output |
|---|---|---|---|---|
| 1 | Candidate sets | Kira's candidates notebook (PR 23, FoldX tiers), `examples/track_a_run.py design` | trastuzumab/HER2 1N8Z, the evidence table | `candidates.csv`: 432 sets of at least 5 substitutions, 400 from the FoldX tiers (Kira, `k001`-`k400`) and 32 from the evidence funnel (`o01`-`o32`) |
| 2 | Fold | `examples/track_a_screen.py fold` (cluster L40S) or `cloud/modal_fold.py::fold_dir` (Modal) | `candidates.fasta`, the parent's cached alignments | Boltz-2 folds of each candidate and of the parent, one batch per slice and seed |
| 3 | Score | `examples/track_a_chem_rank.py score` | the fold folders | `rows_*.csv`: ipSAE (max and min), ipTM, pLDDT and the charge-contact terms per fold |
| 4 | Rank | `examples/track_a_chem_rank.py rank --chem-weight 0.15` | `rows_*.csv`, `candidates.csv` | `ranking.csv` |
| 5 | Pick | `examples/track_a_final_pick.py` | the same | `final15.csv`, `final15.parent.json` |
| 6 | Export | `examples/track_a_export.py` | the picked variant ids and the parent | scFv sequences and `submission.csv` |

Only the 115 candidates `k001`-`k115` went through steps 2 to 5; the 32 of ours were not folded. Each was folded at 3 seeds (345 candidate folds), and the parent was folded in every batch (12 folds). The inputs and outputs of steps 1 to 5 are in `data/track_a_final/`.

```
python examples/track_a_screen.py fold --fasta candidates.fasta --subset 0-39 --seeds 1 --out screen/fold_a
python examples/track_a_chem_rank.py score --root screen/fold_a --out rows_a.csv --workers 8
python examples/track_a_chem_rank.py rank --rows data/track_a_final/rows_*.csv --candidates data/track_a_final/candidates.csv --out ranking.csv --chem-weight 0.15
python examples/track_a_final_pick.py --rows data/track_a_final/rows_*.csv --candidates data/track_a_final/candidates.csv --out final15.csv
python examples/track_a_chem_rank.py audit --ranking data/track_a_final/ranking.csv
```

The last command reproduces `data/track_a_final/final15.csv` from the committed rows files.

## Folding

`boltz predict` was run with the model's defaults (`boltz2`, 3 recycling steps, 200 sampling steps, 1 diffusion sample, no potentials) and `--output_format mmcif --write_full_pae --seed N --cache`. The scFv is chain A, as a protein record whose alignment comes from the parent's cache, so no MSA server is called; the antigen is HER2 residues 480-607 (128 residues). The Boltz version was 2.2.1 in the repository's environment and was not checked again on the cluster.

`track_a_screen.py` folds a whole slice in one process (`msa_cache.batch.fold_batch`), so the weights load once: 10 to 19 s a fold on an L40S, against about 55 s when each design starts its own `boltz predict`. `modal_fold.py::fold_dir --seeds N` runs seeds 0 to N-1 in N containers at the same time. Several processes must not share one input folder; two runs that did lost one slice of candidates, and it was folded again.

The parent was folded in every batch and every candidate is scored against the parent of its own batch and seed, because the poses (and so the charge-contact terms) depend on the seed.

## Scoring

ipSAE is `design_pipeline.scoring.ipsae`, which follows Dunbrack's `ipsae.py` v4 (d0 per aligned residue, `d0 = max(1, 1.24 (max(n, 26) - 15)^(1/3) - 1.8)`, PAE cutoff 15 Å). It reproduces the script's value within 5e-7 on 30 Boltz-2 complexes of the scFv and HER2. Only the PAE cutoff matters for the number; the distance cutoff of the script only changes the counts it reports. `ipsae_ab` aligns on the scFv and scores the antigen, `ipsae_ba` the reverse, `ipsae_max` and `ipsae_min` are the larger and smaller; the ranking uses `max`.

The reference the scorer is checked against lives on the GPU node, not in this repository: `protecher_tools/ipsae/` carries the **v4** script and `ipsae_reference_results.csv`, the 30 complexes ([cluster_gpu_node_rebuild.md](cluster_gpu_node_rebuild.md)). To re-run the comparison by hand there, copy one fold's PAE and structure into an empty folder (the script writes its output table beside its inputs) and run:

```
python protecher_tools/ipsae/ipsae.py pae_<name>_model_0.npz <name>_model_0.cif 15 10
```

Three things were dropped from the repository once nothing used them, and all three are in git history:

| Dropped | What it was | Where it is |
|---|---|---|
| `vendor/ipsae_v3.py` (2026-11-09) | a verbatim copy of the **v3** script, which floors d0 at 1.04 rather than 1.0 and cannot be imported as a function - a command-line script only, which is why the in-process scorer was written separately | `git show 1fb9362:design_pipeline/ipsae.py` (it was added as `design_pipeline/ipsae.py`; the commit subject calls it v4, the file header says version 3, and the 1.04 floor agrees with the header) |
| `design_pipeline/ipsae_score.py` (2026-11-08) | a subprocess wrapper that ran the script and parsed its output table | `git show 1fb9362:design_pipeline/ipsae_score.py` |
| `design_pipeline/rank_ipsae.py` (2026-11-08) | a min-direction ranking layer, superseded by `track_a_chem_rank.py` | `git show 1fb9362:design_pipeline/rank_ipsae.py` |

Nothing imported, called or tested any of them; the only references were in these docs. If the organisers' script ever needs parsing, recover the wrapper from history rather than writing it again. A submission scored by the organisers may use another variant (direction, d0, cutoff); that is not known.

Chemistry is `charge_contact_report` (PR 24), computed per fold against the parent's fold of the same seed: `contact` is the salt-bridge ramp summed over the added cationic substitutions, `coulomb` the parent's minus the design's screened electrostatic energy of the mutated charge sites against the antigen (larger is more favourable). Both are averaged only over the seeds whose pose keeps the parent's epitope (`eligible`). A candidate's chemistry is `chem_z = z(contact) + z(coulomb)` over the candidates that have an on-site seed, and 0 for one that has none. The motivation is a judge's remark that the epitope is rich in Asp and Glu; the spec is [chem_score_spec.md](chem_score_spec.md).

## Ranking and the pick

`rank_weighted` ranks by `score = z(ipSAE max, mean over seeds) + 0.15 * chem_z`. The weight is deliberately small: chemistry reorders near ties and cannot make up a large ipSAE gap. `rank` also reports `in_band` against the parent's mean minus one SD; the pick uses two.

### What the score is made of, and what it is not

The score has exactly three independent inputs, named as `SCORE_INPUTS` in `track_a_chem_rank.py` beside the formula:

| Column | In the score | Nominal weight |
|---|---|---|
| `ipsae_max` | yes, as `z(ipsae_max)` | 1.0 |
| `contact` | yes, inside `chem_z` | 0.15 |
| `coulomb` | yes, inside `chem_z` | 0.15 |
| `ipsae_min`, `iptm`, `in_band`, `n_on_site`, `para_net`, `chem` | no | 0 |

So the composite is `1 * z(ipsae_max) + 0.15 * z(contact) + 0.15 * z(coulomb)` and nothing else. The columns in the lower row are written to `ranking.csv` and `final15.csv` to be read, not scored: `chem` is the band-rule chemistry of `rank_rows` and only `chem_z` reaches the number, `in_band` is a flag derived from `ipsae_max`, and `iptm`, `n_on_site` and `para_net` are context.

**`ipsae_min` is excluded deliberately** (recorded 2026-11-09). It is close to a fixed multiple of `ipsae_max` - about `0.79 * ipsae_max`, leaving a residual SD of 0.023 against the column's own SD of 0.241, so roughly 99% of it is already in `ipsae_max`, and the two order the candidates alike (Spearman about 0.995). Scoring both would have put weight 2 on one co-folding signal against 0.15 each on the two physics terms. It stays as a column, since the same fold yields it for free and it is worth seeing when the two directions disagree, but it must not enter `score`. `tests/test_track_a_chem_rank.py::test_ipsae_min_is_reported_but_never_enters_the_score` perturbs the column and asserts the scores do not move.

A pair of score inputs that measure the same thing is the same bug, so the check is permanent rather than a one-off reading:

- `check_score_correlations(records, threshold=0.9)` computes Pearson r and Spearman rho over `SCORE_INPUTS`, pairwise-complete, and warns on stderr for every pair that reaches 0.9 in either. `rank` calls it on every ranking it writes, so a redundant column cannot be added silently; `--correlation-limit` moves the threshold.
- `python examples/track_a_chem_rank.py audit --ranking ranking.csv` prints both correlation matrices over `AUDIT_COLUMNS` (the three score inputs plus every reported column), the flagged pairs, and each input's share of the score's variance - which is what a nominal weight is actually worth once the inputs correlate. `--drop <column>` prints the shares with a column removed.
- `notebooks/06_score_correlation_audit.ipynb` is the same functions with the `ipsae_min` regression above, the pair-by-pair table, and a re-ranking of the committed `ranking.csv` under each candidate drop, to see whether a drop would move the picks.

The audit over all of `AUDIT_COLUMNS` is expected to flag pairs that are arithmetic rather than findings: `chem` and `chem_z` are the same `z(contact) + z(coulomb)` over different populations (`chem` over the in-band on-site candidates, `chem_z` over all on-site ones with 0 for a candidate with no on-site seed), and `in_band` is a threshold on `ipsae_max`. Only a flagged pair whose two members are both in `SCORE_INPUTS` is a double-weighting.

`track_a_final_pick.py` takes the first 10 by score, then 5 more by `chem_z` from the candidates whose mean ipSAE is at least the parent's mean minus 2 SD and which differ from every earlier pick by at least two substitutions. The parent's ipSAE (max) over its 12 folds is 0.755 with SD 0.091, so the floor is 0.573.

| Design | Candidate | Group | ipSAE max | chem_z | Substitutions |
|---|---|---|---|---|---|
| 01 | k085 | ipSAE + chemistry | 0.801 | 0.98 | A:T94S+B:D31R+B:Y52R+B:R59F+B:Y105W |
| 02 | k063 | ipSAE + chemistry | 0.802 | 0.93 | A:D28W+A:Y49W+A:T94S+B:Y52R+B:D102W |
| 03 | k006 | ipSAE + chemistry | 0.686 | 3.37 | A:D28W+A:S50R+A:Y92W+B:Y52R+B:D102W |
| 04 | k103 | ipSAE + chemistry | 0.780 | 0.36 | A:S50R+A:Y92W+B:D31R+B:Y57W+B:M107R |
| 05 | k105 | ipSAE + chemistry | 0.782 | -0.40 | A:D28W+A:Y49W+A:T94S+B:Y57W+B:G101F |
| 06 | k031 | ipSAE + chemistry | 0.762 | -0.30 | A:D28W+A:Y49W+A:Y92W+B:Y57W+B:Y105W |
| 07 | k090 | ipSAE + chemistry | 0.686 | 2.11 | A:T94S+B:Y52R+B:R59F+B:D102W+B:M107R |
| 08 | k056 | ipSAE + chemistry | 0.776 | -0.86 | A:Y49W+A:Y92W+B:Y33W+B:Y57W+B:D102W |
| 09 | k025 | ipSAE + chemistry | 0.746 | 0.07 | A:D28W+A:Y49W+A:Y92W+B:Y33W+B:D102W |
| 10 | k023 | ipSAE + chemistry | 0.786 | -1.47 | A:D28W+A:Y49W+A:Y92W+B:Y33W+B:R59F |
| 11 | k098 | chemistry | 0.610 | 2.45 | A:D28W+A:T31S+A:Y49W+B:D31R+B:Y105W |
| 12 | k067 | chemistry | 0.597 | 1.66 | A:D28W+A:S50R+A:T94S+B:Y33W+B:Y57W |
| 13 | k028 | chemistry | 0.689 | 0.17 | A:D28W+A:S50R+A:Y92W+B:Y57W+B:G101F |
| 14 | k062 | chemistry | 0.688 | -0.12 | A:D28W+A:F53Y+A:Y92W+B:Y57W+B:D102W |
| 15 | k035 | chemistry | 0.636 | -0.30 | A:D28W+A:F53Y+A:Y92W+B:Y33W+B:Y105W |

Substitutions are in the crystal's numbering (heavy chain B, light chain A). The design numbers are the order of the table, which is the order of the submitted file.

## What this does not show

- Seed noise is large. The parent's 12 folds span 0.58 to 0.85 (mean 0.755, SD 0.091), and its per-seed means were 0.69, 0.76 and 0.82. The best candidates, k085 and k063, average 0.80 over three seeds, inside that range, so the screen does not show that any design beats the parent.
- ipSAE from one predictor says how confident Boltz-2 is in the interface, not that the variant binds tighter. No score we computed ranked affinity in the public data we calibrated on, FoldX included. Nothing here was measured.
- The chemistry weight of 0.15 was set to keep the judge's idea in the ranking, not fitted to anything.
- An AlphaFold2 cross-check of the top 10 against the parent was planned and is not done.
