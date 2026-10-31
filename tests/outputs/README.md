# AntiFold run results (trastuzumab, 1N8Z)

AntiFold scores from the Adaptyv Track A work (improve the trastuzumab scFv against HER2, parent PDB 1N8Z).
`tests/test_antifold_outputs.py` reads the files here and asserts the headline numbers, so the numbers in
`docs/track_a_filtering.md` can be recomputed from the repo. The files are results copied from the cluster node; the test
runs nothing. About 0.2 MB. No FoldX output is included (its licence does not allow sharing it).

"(not checked)" marks a statement taken from an earlier session's note or a repo document and not verified here.

## What these results are

A calibration of AntiFold, not a set of designs. AntiFold scored trastuzumab variants whose binding had already been
measured by others (AbBiBench, SKEMPI), to see whether its score tracks binding. No candidate designs for submission are
in this folder.

AntiFold scores how likely an antibody chain's sequence is, given the backbone of the whole antibody (heavy and light chain
from 1N8Z, no antigen). A variant that AntiFold finds more likely than the wild type gets a higher score. The question is
whether that agrees with measured binding.

| Data | n | Spearman of AntiFold score with measured value | Same, mutation count only (no model) | By mutation count |
|---|---|---|---|---|
| AbBiBench 1N8Z file | 419 | 0.19 | 0.26 | 1 to 3: 0.69 (n 23); 4 to 6: -0.09 (n 86); 7 to 13: 0.15 (n 310) |
| AbBiBench 4D5 file | 2,079 | 0.31 | 0.23 | 4 to 6: 0.63 (n 645); 7 to 9: 0.07 (n 1,429) |
| SKEMPI 1N8Z rows, against minus ddG | 8 | 0.60 | not computed | 6 single-mutation rows: 0.71 |

Unrounded: 0.1945, 0.2568, 0.6947, -0.0946, 0.1515 (1N8Z file); 0.3141, 0.2318, 0.6263, 0.0747 (4D5 file); 0.595 and
0.714 (SKEMPI). The test asserts each at a tolerance of 0.005.

How to read it:

- **Spearman** is the rank correlation between AntiFold's score and the measured value over the variants. +1 means the
  same order, 0 means unrelated.
- **Mutation count only** is the same correlation for minus the number of mutations: a "scorer" that only penalises
  change. AntiFold has to beat it to add anything. It does on the 4D5 file (0.31 against 0.23) and not on the 1N8Z file
  (0.19 against 0.26).
- **By mutation count** compares variants with similar numbers of mutations. Below 7 mutations the result is mixed:
  0.69 for 1 to 3 (1N8Z file, only 23 variants), -0.09 for 4 to 6 (1N8Z file) and 0.63 for 4 to 6 (4D5 file). At 7 or
  more it is weak in both files (0.15 and 0.07). `docs/track_a_filtering.md` draws the conclusion that AntiFold works as
  a filter and not as a ranking.
- **SKEMPI** has 8 rows, so its correlations are anecdotes. Positive ddG means weaker binding than the wild type, so the
  test correlates the score with minus ddG.

## How AntiFold was run

- **Tool:** `oxpig/AntiFold`, source clone on node fgcz-r-023 at commit `9918344067f34a6b026f2634fdba212fc3e6e986`
  (2026-09-27, "ANARCII (#35)"), version 0.4.0. Read from the clone and its venv on 2026-10-31.
- **Weights:** `src/models/model.pt`, 566,838,063 bytes, SHA-256
  `d5c442fa0372c28f4d0026d2f551b6f8ba7e7a127cb6837813a88093ed233e9e` (computed 2026-10-31; the file was written
  2026-10-30 21:03, before the first run below). Its origin as release asset `weights-v1/model.pt` is from an earlier
  session's report (not checked).
- **Glue:** `antifold_score` in this repo (AbBiBench) and `design_pipeline.score_variants(...,
  models=["antifold_seq"])` (SKEMPI). The repo commit the node copy was built from was not recorded.
- **Input:** heavy and light chain only, no antigen. Structure files on the node: `1N8Z.pdb` (A light, B heavy, C HER2)
  and `1N8Z.cif`. Public PDB entry, not copied here.
- **Score (`score`):** change in total log-likelihood against the wild type, in nats, summed over the chains the
  candidate touches. 0 means as likely as the unmutated 1N8Z sequence, positive means more likely, higher is better.
  Several mutations in one chain are scored together, not as a sum of singles. A candidate with mutations in both
  chains gets the sum of the two chain scores. (From `antifold_score/scorer.py` and
  `design_pipeline/adapters/_antifold_seq.py`.)

## abbibench_trastuzumab/

Trastuzumab variants scored on 1N8Z (heavy B, light A), date 2026-10-30.

- **Source:** AbBiBench, "AbBiBench: A Benchmark for Antibody Binding Affinity Maturation and Design" (arXiv
  2506.04235), dataset https://huggingface.co/datasets/AbBibench/Antibody_Binding_Benchmark_Dataset, CC BY 4.0
  (licence read on the dataset page on 2026-10-31). Files used: `1n8z_benchmarking_data.csv` and
  `4d5_her2_benchmarking_data_trimmed.csv` (names as downloaded by the earlier session, not re-checked). The raw files
  are not copied.
- **What the two files are:** the 1N8Z file has 419 variants, all with heavy-chain mutations only (3,204 substitutions in
  total, 1 to 13 per variant). The 4D5 file has 2,080 variants, 2,076 of them with light-chain mutations as well
  (7,503 heavy and 7,118 light substitutions; 3 to 9 per scored variant). The 1N8Z variants sit at heavy-chain
  positions 29 to 33, 52 to 55, 57 and 97 to 109, the 4D5 variants at heavy 100, 101, 102 and 106 and light 91 to 94 and
  97 (positions 1-based into the 1N8Z chains). `docs/track_a_filtering.md` says the files hold no framework variants, so
  they calibrate AntiFold on loop positions only (not checked against a CDR definition).
- **Ids:** `n8z_NNNN` and `d5_NNNN` are our ids: row NNNN (1-based, in file order) of the two AbBiBench files. The
  mutations are the positions where the AbBiBench heavy and light sequence differs from 1N8Z chains B (first 120
  residues) and A (first 107), after cutting the extra residues AbBiBench adds (`A` after VH, `RT` after VL).
- **Commands** (node, AntiFold environment, `AF=/scratch/foertsch/protecher_tools/antifold`, default CPU device). The
  `cand_*.csv` inputs (one row per mutation: `variant_id,chain,position,wt,mut`) are not kept in the repo:

```bash
python -m antifold_score --pdb $AF/1N8Z.pdb --heavy B --light A --antifold_src $AF/src --mutations_csv cand_n8z.csv --out_dir out_n8z_complex
python -m antifold_score --pdb $AF/1N8Z.pdb --heavy B --light A --antifold_src $AF/src --mutations_csv cand_d5.csv  --out_dir out_d5_complex
python -m antifold_score --pdb $AF/1N8Z.pdb --heavy B --light A --antifold_src $AF/src --mutations_csv cand_singles.csv --route table --out_dir out_singles_table
```

| File | What it is | Rows | Columns |
|---|---|---|---|
| `scores_n8z_complex.csv` | AntiFold's own function for scoring a whole mutated chain given the backbone (`complex` route, the default), 1N8Z file | 419 | `variant_id`; `score` (nats, see above; range -42.3 to 6.5, median -7.5); `raw` (total log-likelihood of the whole mutant, nats); `note` (empty here, otherwise says why a candidate got no score) |
| `scores_d5_complex.csv` | Same, 4D5 file | 2,079 | same (`score` range -39.2 to 3.4, median -22.1) |
| `scores_singles_table.csv` | `table` route: AntiFold's per-residue table for the wild-type structure, exact for single mutations only. The 236 distinct single mutations that occur in either file, id `<chain><position><wt>><mut>` (e.g. `A91H>F`). `score` is (wild-type mean minus mutant mean) times the table rows, in nats; `raw` is the mutant's mean negative log-likelihood (lower fits better). Not used by the headline numbers | 236 | same (`score` range -13.2 to 3.1) |
| `meta_n8z.csv` | Mutation counts and measured value, 1N8Z file | 419 | `variant_id`; `n_heavy`, `n_light` (substitutions in each chain against 1N8Z); `binding_score` (below) |
| `meta_d5.csv` | Same, 4D5 file | 2,080 | same |

`binding_score` is AbBiBench's measured value, copied unmodified, higher meaning stronger binding. **It is not on the same
scale in the two files, so do not compare values across them.** In the 1N8Z file it runs from 5.6 to 9.25 and the
AbBiBench paper lists it as -log Kd, from Shanehsazzadeh et al. 2023 (read on the arXiv v2 page on 2026-10-31 through a
page-summarising tool, not read in full, units not checked). In the 4D5 file it runs from -3.73 to 5.12; the paper
page I read did not describe that file, so its source and meaning are not documented here (not checked).

`meta_d5.csv` has one row more than the scores: `d5_0991` has no mutation (the unmutated sequence, `binding_score`
2.166), so it was not a candidate and has no score.

## skempi_1n8z/

The eight SKEMPI 2.0 rows of the 1N8Z block (ids 5176 to 5183), scored with `antifold_seq` on `1N8Z.cif` (heavy B, light
A, GPU), one job on 2026-10-31: three VL single mutations (N30S, H91F, T94S), three VH singles (D102W, Y105F, Y109V), a
VL triple (N30S, H91F, Y92W) and the seven-mutation combination of all of them. SKEMPI:
https://life.bsc.es/pid/skempi2/database/details/5183, CC BY 4.0 (per `docs/research/08_trastuzumab_variant_datasets.md`,
not checked); reference field PMID 12206766, identified as Gerstner 2002 by an earlier agent (not checked). SPR at
310 K, wild-type KD 0.35 nM; chain A is VL, chain B is VH. The score column was supplied by the session that ran the
job; its output is not kept here.

| File | What it is | Rows | Columns |
|---|---|---|---|
| `skempi_1n8z_antifold_seq.csv` | KD, ddG and score for the eight rows | 8 | `id` (SKEMPI id); `label` (mutations, quoted because the last one contains commas); `kd_nM` (measured KD, nM); `ddg_exp_kcal_mol` (= RT ln(KD / 0.35 nM), RT 0.616 kcal/mol, two decimals, positive is weaker than the wild type); `score_antifold_seq` (nats, same definition as `score` above, range -5.9 to 3.8) |

## What cannot be recovered from these files

- Which mutations each AbBiBench variant has. The candidate tables are not kept, and the ids depend on the row order of the
  AbBiBench files at the time (dataset version not recorded). The rule above regenerates them from the dataset.
- The exact commands, since they read `cand_*.csv` files that are not kept. The route, structure, chains and flags are
  above.
- The repo commit the node's `antifold_score` and `design_pipeline` were built from.
- The SKEMPI job output, the per-residue AntiFold tables and the raw job logs.

## Not included

- FoldX output of any kind (licence).
- The cetuximab scFv scores (a different antibody and target, scored on 1YY9). The cetuximab numbers in
  `docs/track_a_filtering.md` therefore cannot be recomputed from this folder.
- The combine-tuning sets (300 variants each, `antifold_seq` on 1N8Z, 2026-10-31): they stay on the node at
  `/scratch/foertsch/protecher_combine_tune/runs/scored_30`, `scored_45` and `scored_60` (node-local scratch of
  fgcz-r-023, file `scored_sets.parquet`, about 23 KB each).
