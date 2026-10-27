# 01. Filtering de novo binders: which in-silico metrics predicted wet-lab binding

Research date 2026-10-24. Source: subagent report, merged and edited. "Rec" marks a number that the research agent recomputed from a public table named in the text. It is not a published value. "UNVERIFIED" marks anything that rests on a secondary source or that could not be opened.

## Bottom line

1. No single in-silico metric predicts binding well in Adaptyv data. AUROCs on tested designs fall between 0.55 and 0.83 with wide intervals. Scores enrich a pool but do not pick winners reliably. The AUROC values in this file are recomputed by the research session from public tables (marked Rec), not published results. Do not quote them to a jury as published.
2. The best single score in the literature is the interface-restricted ipSAE_min from an AF3-class model. It gives AUROC 0.79 and average precision 0.35 (base rate 10.7%) over 3,669 designs on 15 targets (Rec from the Overath et al. table), and AUROC 0.83 on the 302 Adaptyv EGFR designs in that table (16 binders, interval 0.71 to 0.92).
3. Boltz-2 ipSAE, which Adaptyv now selects on, gave AUROC 0.63 on the Nipah competition (1,028 designs) and 0.74 on the TREM2 hackathon (110 Proteinbase rows, 46 binders). Both are Rec from Proteinbase, and both pools were already cut on ipSAE, so the values are understated. The 0.74 is for the `boltz2_ipsae` column. The stricter `boltz2_min_ipsae` column gives 0.62 on the same rows (re-run 2026-10-26), so which ipSAE variant is ranked matters. The TREM2 table has 110 rows and 46 binders, while the organiser's post reports 100 tested and 37 bound. The extra rows are unexplained.
4. ipTM and AF2 iPAE rank worse than ipSAE (EGFR round 2: AUROC 0.64 and 0.61) but work as gates: iPAE at most 10 A, binder pLDDT at least 80, ipTM above 0.5.
5. ESM2 pseudo-log-likelihood (PLL) tracks similarity to known binders. Per-residue AUROC is 0.73 on all EGFR designs and 0.43 among de novo ones. The unnormalised version was gamed by trimming sequences. Rosetta interface terms add modest signal (AUROC about 0.65). No FoldX benchmark on binder success was found.
6. Averaging two or three co-folding models helps a little. Requiring all models to agree raises precision and loses recall, so rank by the mean and keep the AND rule for the last 10 slots.
7. Absolute cut-offs do not transfer between targets. Overath's F1-optimal ipSAE_min thresholds span 0.5 to 0.8 with precision from 0.1 to 1.0. Treat any cut as a rank inside our own pool.
8. The EGFR and Proteinbase datasets are public under open licences. No dataset contains HER2.

## 1. Checking the numbers in docs/challenge_notes.md

Primary source for this section: the EGFR paper, https://www.biorxiv.org/content/10.1101/2025.04.17.648362.full.pdf (CC BY 4.0).

- Round 1: 726 designs submitted, 201 tested (top 100 by AF2-Multimer iPAE plus 101 picked by Adaptyv), 146 expressed, 5 binders (5/201 = 2.5%). Two more binders were disqualified for similarity to known sequences. The notes are right.
- Round 2: 1,131 submitted, 400 tested (100 ranked plus 300 hand-picked), 378 expressed, 53 binders (13%). Totals over both rounds: 1,857 submitted, 601 characterised by BLI, best optimised binder 1.21 nM, best de novo binder 82 nM (BindCraft). The Adaptyv blog post for round 2 gives 91.5 nM for the best de novo binder, so treat 82 nM as the paper's number and check before quoting (see file 05).
- Ranking: round 1 used iPAE alone. Round 2 used the rank average of unnormalised ESM2-650M PLL, iPAE and ipTM (https://github.com/adaptyvbio/competition_metrics). The notes omit that most tested designs were hand-picked.
- Later format, confirmed for the muni x Adaptyv TREM2 hackathon (one day, 16 teams): at most 10 designs per team, 141 received, top 100 by Boltz-2 ipSAE_d0res (PAE cutoff 15, distance cutoff 15) tested (https://www.adaptyvbio.com/blog/agents-vs-humans). The Nipah competition tested 600 designs by collection-average ipSAE, 400 by expert choice and 200 by community vote (https://www.adaptyvbio.com/blog/nipah-submissions). Whether HerCode copies either is unannounced.
- "1 to 5 designs per designer" appears only in a secondary write-up and in no primary source we could open. UNVERIFIED.

## 2. EGFR rounds 1 and 2

The paper reports, for the 400 round-2 designs, ROC AUC for binding of 0.72 (esm_pll), 0.66 (pLDDT) and 0.64 (ipTM). For high expression the AUCs are 0.58 (ipTM), 0.55 (pLDDT) and 0.44 (esm_pll). The authors conclude that the filter metrics cannot accurately predict binding and find no significant correlation between KD and any ranking metric.

Rec from https://github.com/adaptyvbio/egfr2024_post_competition (all_submissions.csv; 60 binders among 523 labelled tested designs, round 2 alone 53 of 377):

| Metric | R1+R2 AUROC | R2 AUROC (95% CI) |
|---|---|---|
| AF2-M iPAE (lower is better) | 0.59 | 0.61 (0.52 to 0.70) |
| AF2-M ipTM | 0.62 | 0.64 (0.56 to 0.72) |
| Binder pLDDT | 0.66 | 0.65 (0.58 to 0.73) |
| ESM2 PLL, unnormalised | 0.51 | 0.55 (0.45 to 0.64) |
| ESM2 PLL per residue | 0.72 | 0.73 (0.67 to 0.80) |
| Rosetta interface dG / interface H-bonds | 0.66 / 0.69 | 0.65 / 0.68 |

- The paper's 0.72 matches the per-residue value, so the paper's PLL is most likely length-normalised.
- Selection arms (Rec). Round 1: the top 100 by iPAE held 3 binders and the 101 hand-picks held 4, so iPAE-only ranking gave no enrichment. Round 2: the rank-average top 100 held 24 binders against 29 of 300 hand-picks (24% vs 9.7%), but 83 of those 100 were variants of known binders and only 4 were de novo.
- Hit rate by the designers' own category (Rec): diversified known binders 32/184, optimised binders 13/88, hallucination 8/78 (10%), de novo 2/163 (1.2%). These categories are what designers declared in `all_submissions.csv`. Adaptyv's blog and file 05 use another split (BindCraft de novo 9%, 6 of 65), and file 03 counts 7 of 8 de novo binders as BindCraft. The three figures cover different populations. BindCraft designs probably sit under "hallucination" in the CSV, which would explain why that row has 10% and the "de novo" row 1.2% (inference, not confirmed). For planning, use 1% to 10% and expect the higher end only for BindCraft-style hallucination.
- Among the 215 de novo and hallucination designs (10 binders, Rec): AUROC 0.80 for ipTM, 0.77 for interface pLDDT, 0.73 for pLDDT, 0.70 for iPAE, 0.78 for interface H-bonds, 0.38 for unnormalised PLL.

## 3. ipSAE

- Dunbrack (https://www.biorxiv.org/content/10.1101/2025.02.10.637595, code https://github.com/DunbrackLab/IPSAE) defines ipSAE and tests it on 40 PDB heterodimers against 70 mismatched pairs. PAE cutoffs of 10 or 15 A look most suitable. It has no wet-lab binder validation and gives no binder threshold.
- Overath et al. (https://www.biorxiv.org/content/10.1101/2025.08.14.670059v1) compare 3,766 designs on 15 targets (436 binders), re-predicted with AF2 initial guess, ColabFold, AF3 and Boltz-1. AF3 ipSAE_min was the best single feature, with 1.4 times the average precision of AF2-IG ipAE. The pooled F1-optimal cut is 0.61. Per-target optima mostly lie at 0.5 to 0.8 with precision from 0.1 to 1.0. Their second recipe pre-filters on input shape complementarity above 0.62 and binder RMSD below 3.73 A, then takes the top K by ipSAE_min. Some of their data had been pre-selected on ipAE, which deflates ipAE.
- Rec from their Zenodo table (pooled AUROC / AP): AF3 ipSAE_min 0.787 (CI 0.761 to 0.814) / 0.350; AF3 ipTM 0.751 / 0.229; Boltz-1 ipSAE_min 0.768 / 0.324; ColabFold ipSAE_min 0.748 / 0.285; AF2-IG pae_interaction 0.721 / 0.218. The recomputed AP ratio to ipAE is 1.6, where the authors report 1.4. At AF3 ipSAE_min above 0.61: precision 0.38 and recall 0.52 pooled, precision 0.26 and recall 0.36 on EGFR.
- Adaptyv EGFR subset of that table (302 designs, 16 binders), AUROC: AF3 ipSAE_min 0.83, ColabFold ipSAE_min 0.87, Boltz-1 ipSAE_min 0.76, AF3 ipTM 0.69, AF2-IG iPAE 0.67, binder pLDDT 0.51.
- Boltz-2 ipSAE on later Adaptyv data (Rec from https://storage.proteinbase.com/proteinbase_latest.csv, snapshot 2026-10-12; pools were pre-selected on this score, so AUROC is understated):
  - Nipah (aggregate statistics only, 1,028 labelled designs, 101 binders = 9.8%): AUROC 0.63 (0.58 to 0.69), ipTM 0.60. The top 100 by ipSAE held 16 binders. ipSAE of at least 0.8 gave 14.2% precision (1.45 times enrichment).
  - TREM2 (110 Proteinbase rows, 46 binders; the organiser's post says 141 received, top 100 tested, 37 bound, and the difference is unexplained): AUROC 0.74 (0.64 to 0.83) for `boltz2_ipsae`, 0.62 for `boltz2_min_ipsae`, ipTM 0.69. Reproduced on 2026-10-26. Hit rate by ipSAE band: 17/21 at 0.9 or more, 22/53 from 0.8 to 0.9, 6/22 from 0.7 to 0.8, 1/14 below 0.7.
  - Promera reports AUROC 0.70 on Nipah for its own ipSAE variant and 0.66 for the next method (916 designs, 94 binders; https://www.biorxiv.org/content/10.64898/2026.06.07.729267v1.full).
- For nanobody or scFv formats, a 2026 benchmark found ipSAE rose when CDRs were mutated and called it the worst confidence measure in its alanine scan (https://www.frontiersin.org/journals/bioinformatics/articles/10.3389/fbinf.2026.1858891/full). This matters if the challenge allows antibody formats.

## 4. ipTM, iPAE and pLDDT

- Bennett et al. (https://pmc.ncbi.nlm.nih.gov/articles/PMC10163288/): AF2 initial-guess pAE_interaction showed a sharp rise in success below 10. Filtering at pAE_interaction below 10 and complex RMSD below 5 A (about 2.3% of designs pass) gave 8-fold and 30-fold higher experimental success than Rosetta-only filtering on LTK and IL-2Ralpha site 1. The ROC curves are in supplementary figures only.
- BindCraft (https://www.biorxiv.org/content/10.1101/2024.09.30.615802v2, filters at https://github.com/martinpacesa/BindCraft/blob/main/settings_filters/default_filters.json) gates on pLDDT above 0.8, i_pTM above 0.5, i_pAE below 0.35 (normalised by 31, about 10.9 A), shape complementarity above 0.55 per model (0.6 on the average), binder RMSD below 3.5 A, surface hydrophobicity below 35%, dG below 0, at least 7 interface residues and 3 H-bonds. The authors rank by i_pTM and call it a good binary predictor of binding and a poor predictor of affinity. They order the top 10 of 20 from at least 100 passing designs.
- EGFR round 2 (Rec): iPAE of at most 10 A gave 24 binders in 105 designs (22.9% vs 14.0% base rate, recall 45%). None of 5 designs under 7.5 A bound. ipTM of at least 0.9 gave 30 in 160 (18.8%). pLDDT of at least 90 gave 11 in 40 (27.5%) and at least 85 gave 20.7%. Round 1 (Rec): 1 of 20 designs under 10 A bound, and the five valid binders had iPAE 9.2 to 21.4. Adaptyv's own post doubts iPAE as an affinity proxy (https://www.adaptyvbio.com/blog/po102).
- Binder pLDDT AUROC is 0.67 pooled in Overath's data but 0.51 on its Adaptyv subset (Rec). The nanobody benchmark found pLDDT best for spotting out-of-distribution sequences and weaker for picking binders.

## 5. ESM PLL, Rosetta and FoldX

- PLL: Adaptyv added it in round 2 as an expression proxy and said it may favour natural-looking designs (https://github.com/adaptyvbio/competition_metrics). It did not predict expression (AUROC 0.44). Designers trimmed binders to raise it, and no top-100 sequence exceeded 100 residues (https://www.adaptyvbio.com/blog/po103). It works as a similarity-to-known-binders signal and falls below chance among de novo designs (numbers in section 2). If the organisers use unnormalised PLL again, short sequences are favoured, which matters for the length choice.
- Rosetta: Overath found that Rosetta dG, dSASA and shape complementarity add information to ipSAE, although extra features helped only some targets. Input-structure terms alone are weak (Rec pooled AUROC 0.54 for shape complementarity, 0.58 for dG, 0.70 for dG on EGFR). On EGFR, binder per-residue Rosetta total energy gives AUROC 0.62 when inverted, EvoEF2 0.61 and BUDE 0.56 (Rec from the DE-STRESS table in the round-2 repo). Bennett et al. call Rosetta ddG an effective filter that AF2 initial guess beats.
- FoldX: no benchmark on binder success was found. Report it, do not cut on it.

## 6. Cross-model agreement

- Overath: AF3 beat Boltz-1 and ColabFold, and pooling features from all models did not improve median AP. Rec on their table: the rank mean of AF3 and Boltz-1 ipSAE_min gives AUROC 0.794 against 0.786 for AF3 alone. Requiring all three models above 0.61 gives precision 0.415 and recall 0.35 (AF3 alone: 0.377 and 0.52). On the 293 Adaptyv EGFR designs, 6 of 12 bound under the three-model rule against 6 of 30 for AF3 alone (tiny n).
- Anthropic technical report (https://www-cdn.anthropic.com/30bf50e22a01388bb29bf077ee3f244531594b7a.pdf, 2026-08-18): on the Overath data (3,532 designs, 13 targets, 391 binders) macro-AP was 0.55 for AF3, 0.61 to 0.62 for ESMFold2 variants and 0.66 for a z-score ensemble of three predictors (+0.11 over AF3, interval 0.02 to 0.26). On Claude's own 1,320 designs the ensemble gave mean AP 0.52 against 0.31 by chance. The hit rate was 49% for rank 1, 44% for the top 5, 39% for the top 10 and 28% over all 30. Scores tracked affinity weakly (Spearman 0.16), and failed targets scored almost as high as successful ones.
- ProtDBench (https://arxiv.org/abs/2605.04118, Cao data, 9 targets, 1,385 binders): AF2-IG, Boltz-1, Boltz-2, Chai-1, ColabFold and Protenix recover largely different binders at top-1% ipTM, and recall falls sharply as more verifiers must agree.

## 7. Practical points from later Adaptyv rounds

- Hit rates on Proteinbase competition cards (https://proteinbase.com/competitions): EGFR R1 2.5% of 201, EGFR R2 14% of 400 (probably counts the two controls, 55 of 400; the paper says 13%), Nipah 9.3% of 1,200 (best 0.37 nM), GEM x Adaptyv RBX1 2.8% of 322 (best 23.7 nM; selection rule not found). The ongoing Anthropic x Adaptyv competition says selection will weigh several in-silico metrics and will not hinge on one (https://proteinbase.com/competitions/anthropic-adaptyv-2026).
- Reproducibility: Adaptyv found that local reruns differ with hardware and seed, and that ipSAE variance was uneven across designs (Nipah post above). One team's best design scored 0.7896 locally and 0.858 officially (https://seqera.io/blog/building-a-scalable-nextflow-protein-design-pipeline-for-adaptyv-s-nipah-binder/, secondary). Run Adaptyv's pipeline (https://github.com/adaptyvbio/nipah_ipsae_pipeline, MIT, Colab, precomputed MSA) with several seeds, rank on the median across seeds and apply floors to the minimum. Taking the best seed would inflate the local score relative to Adaptyv's single official run, which is the gap this paragraph describes. Keep a margin above any cut. The known-binder prescreen was updated on 2026-10-20 (https://github.com/adaptyvbio/binder-prescreen). Run it before submitting.
- HER2: none of the 35 target slugs in Proteinbase is HER2. One small academic study designed HER2 domain IV minibinders on the trastuzumab epitope (RFdiffusion, ProteinMPNN, AF2 initial guess; target and numbers re-read on the PMC page by the merging session). Filtering at pAE_interaction below 10 and pLDDT above 70 left 31 of about 5,000 designs, 25 were tested by E. coli surface display, 2 bound, and the best had Kd 19 nM by ITC (https://pmc.ncbi.nlm.nih.gov/articles/PMC12650588/). The paper's abstract and results disagree on how this affinity compares with trastuzumab, so quote the 19 nM only. Reading: expect de novo hit rates of 1% to 10% on a new receptor (EGFR de novo 1.2%, EGFR hallucination 10%, Nipah 9.3%), and higher on easy targets after ipSAE selection (TREM2: 37 of 100).

## 8. Recommended thresholds

Rows marked (judgment) are a suggestion without direct data. Sources are in section 10.

| Metric / model | Use | Value | Evidence |
|---|---|---|---|
| Boltz-2 ipSAE (the column Adaptyv scores; direction not documented, see section 10; 3+ seeds) | Primary rank and floor | Rank on the median over seeds; floor 0.5 on the minimum over seeds; aim at 0.7 to 0.8 or higher; rank inside the pool | F1-optimal 0.61 [S5]; TREM2 hit rate 17/21 at 0.9+ (Rec) [S9]; Adaptyv selects on it [S2, S4] |
| ipTM, any model | Floor, tie-break | Above 0.5; prefer 0.8+ | BindCraft gate [S7]; AUROC 0.64 EGFR R2 [S1], 0.60 Nipah and 0.69 TREM2 (Rec) [S9] |
| AF2-IG iPAE (only if run) | Gate | At most 10 A | [S8]; EGFR R2 23% vs 14% (Rec) [S11] |
| Binder pLDDT | Fold gate | At least 80; tie-break 85 to 90 | [S7]; AUROC 0.66 [S1]; 90+ gives 27.5% vs 14% (Rec) [S11] |
| Binder RMSD, design vs co-fold | Gate | Below 3.5 A (BindCraft), 3.73 A (Overath) | [S5, S7]; Anthropic: adds nothing to ipSAE ranking [S10] |
| Shape complementarity, design model | Gate | At least 0.6 on average; 0.62 | [S5, S7]; EGFR AUROC 0.60 to 0.69 (Rec) [S11] |
| Rosetta dG, dSASA, interface residues, H-bonds | Gate | dG below 0; dSASA at least 1; 7+ residues; 3+ H-bonds | [S7]; EGFR AUROC 0.65 to 0.69 (Rec) [S11] |
| Cross-model | Rank by z-score mean of 2 or 3 models; AND rule for the last 10 slots | All models above 0.61 | Rec on [S16]; [S10]; [S12] |
| ESM2 PLL per residue | Outlier flag (judgment: bottom 5 to 10%) | Never use the unnormalised value | [S1]; below chance among de novo designs (Rec) [S11] |
| FoldX, total Rosetta energy | Report only | No cut | No supporting benchmark; DE-STRESS energy AUROC about 0.6 (Rec) |

## 9. Dataset availability

| Dataset | URL | Format and size | Licence |
|---|---|---|---|
| EGFR round 1 | https://github.com/adaptyvbio/egfr_competition_1 | results/result_summary.csv, 201 designs plus a control (plddt, pae_interaction, kd; no ipTM or PLL); replicate_summary.csv 461 rows; raw kinetics 131 MB zip | Code Apache 2, data ODbL |
| EGFR round 2 | https://github.com/adaptyvbio/egfr_competition_2 | result_summary.csv, 400 designs plus 2 controls (binding true/false/unknown = 53/325/22, pae_interaction, esm_pll, iptm, plddt); DE-STRESS CSVs (70 columns); structure_predictions.zip 38 MB; embeddings.zip 70 MB; package.zip 406 MB | Same |
| All 1,857 submissions | https://github.com/adaptyvbio/egfr2024_post_competition, adaptyv/data/all_submissions.csv (7.5 MB) and .parquet (1.9 MB) | Git LFS: raw GitHub serves a 132-byte pointer, so fetch https://media.githubusercontent.com/media/adaptyvbio/egfr2024_post_competition/main/adaptyv/data/all_submissions.csv. Includes iPAE, ipTM, PLL, Rosetta interface terms, selection arm, binding label | No licence file |
| Proteinbase (EGFR, Nipah, TREM2, RBX1, PD-L1 and more) | https://huggingface.co/datasets/adaptyvbio/proteinbase, https://storage.proteinbase.com/proteinbase_latest.csv | 4,069 proteins, 6 columns, nested JSON "evaluations" (KD, kon, koff, binding, expression, Boltz-2 ipSAE, ipTM, pLDDT, LIS); parquet 3.7 MB, CSV 46 MB; Boltz-2 metrics for 2,194 rows, none for EGFR or RBX1 | ODC-By |
| Overath meta-analysis | https://zenodo.org/records/15722219, code https://github.com/DigBioLab/de_novo_binder_scoring | final_dataset.csv 82 MB, 312 columns, 3,676 rows (3,669 labelled, 394 binders; the paper's abstract says 3,766 designs and 436 binders, reason unchecked, so quote 3,669 only as the size of the table we analysed); tarballs 0.24 to 3.8 GB | CC BY 4.0, MIT code |
| Anthropic binder campaign | https://huggingface.co/datasets/Anthropic/claude-protein-binder-design | 1,440 designs, 1,320 with binding data, ten predictors x 5 seeds with PAE matrices; about 79 GB | CC BY 4.0, MIT scripts |

Not public or not found: any HER2 data; scores for Nipah or TREM2 designs that were submitted but not tested; HerCode scoring rules; the Anthropic x Adaptyv selection prompt (promised after the event); the RBX1 selection rule.

Suggested use: download all_submissions.csv (via the media URL) and the Proteinbase CSV before 31 Oct, and fit the stage weights of the filter on them. Expect noisy weights, because the EGFR positives number 60 and the TREM2 positives 46.

## 10. What could not be verified

- The HerCode metric, target, length limit and format. The event page only says "a protein design challenge by Adaptyv" (re-read 2026-10-24). Everything else is inference from the TREM2 and Nipah formats.
- AUROC and AP values that exist only in figures (Bennett, Overath, ProtDBench, per-predictor results in the Anthropic report). They were recomputed from released tables where possible. The recomputed Overath AP ratio (1.6) differs from the authors' (1.4), and the released CSV is smaller than the abstract states.
- What Proteinbase's `boltz2_ipsae` column is exactly. The download page does not say, and neither do the TREM2 post, the Nipah pipeline README or Dunbrack's README (which defines asym and max per chain pair and no min). Proteinbase also carries `boltz2_min_ipsae`. On the 110 TREM2 rows `boltz2_ipsae` is at or above the min column in every row, so it is probably the max or one fixed direction (inference). The check is to rescore public TREM2 designs and compare. designMethod is blank for 36% of rows, so per-method hit rates are unreliable.
- Secondary claims, all UNVERIFIED: "1 to 5 designs per designer"; per-method Nipah hit rates for BoltzGen, BindCraft and RFdiffusion from a community summary; one team's report that its local ipSAE did not correlate with the official score (that blog returned HTTP 403).
- The recomputed AUROCs depend on label definitions that differ by source (Overath KD below 10 uM, Adaptyv's own labels), on pre-selection by the same metrics, and on small positive counts (10 to 16 in several subsets).
- Viral-target rows in Proteinbase were used only as aggregate statistics. No pathogen sequences were collected.

## 11. Sources

| ID | URL | Supports | Type |
|---|---|---|---|
| S1 | https://www.biorxiv.org/content/10.1101/2025.04.17.648362.full.pdf | EGFR paper: counts, ranking, AUROC tables | Primary |
| S2 | https://www.adaptyvbio.com/blog/agents-vs-humans | TREM2 hackathon format and results | Primary (vendor) |
| S3 | https://github.com/adaptyvbio/egfr_competition_1, https://github.com/adaptyvbio/egfr_competition_2, https://github.com/adaptyvbio/competition_metrics | EGFR data, ranking formula, licences | Primary |
| S4 | https://www.adaptyvbio.com/blog/nipah-submissions | Nipah selection rule, ipSAE variance | Primary (vendor) |
| S5 | https://www.biorxiv.org/content/10.1101/2025.08.14.670059v1 | Overath meta-analysis and thresholds | Primary |
| S6 | https://www.biorxiv.org/content/10.1101/2025.02.10.637595 | Dunbrack ipSAE definition | Primary |
| S7 | https://www.biorxiv.org/content/10.1101/2024.09.30.615802v2, https://github.com/martinpacesa/BindCraft | BindCraft filters | Primary |
| S8 | https://pmc.ncbi.nlm.nih.gov/articles/PMC10163288/ | Bennett iPAE below 10 rule | Primary |
| S9 | https://storage.proteinbase.com/proteinbase_latest.csv, https://proteinbase.com/download | Proteinbase data; Nipah and TREM2 recomputation | Primary data |
| S10 | https://www-cdn.anthropic.com/30bf50e22a01388bb29bf077ee3f244531594b7a.pdf | Anthropic campaign and ensemble benchmark | Primary (technical report) |
| S11 | https://github.com/adaptyvbio/egfr2024_post_competition | all_submissions.csv | Primary data |
| S12 | https://arxiv.org/abs/2605.04118 | ProtDBench verifier agreement | Primary (preprint) |
| S13 | https://www.biorxiv.org/content/10.64898/2026.06.07.729267v1.full | Promera Nipah AUROC | Primary (preprint, vendor method) |
| S14 | https://www.frontiersin.org/journals/bioinformatics/articles/10.3389/fbinf.2026.1858891/full | Nanobody co-folding benchmark | Primary |
| S15 | https://pmc.ncbi.nlm.nih.gov/articles/PMC12650588/ | HER2 domain IV minibinder study | Primary |
| S16 | https://zenodo.org/records/15722219 | Overath dataset | Primary data |
| S17 | https://huggingface.co/datasets/adaptyvbio/proteinbase, https://huggingface.co/datasets/Anthropic/claude-protein-binder-design | Dataset cards, sizes, licences | Primary |
| S18 | https://proteinbase.com/competitions, https://proteinbase.com/competitions/anthropic-adaptyv-2026, https://proteinbase.com/competitions/adaptyv-egfr-binder2 | Competition cards and rules | Primary (vendor) |
| S19 | https://www.adaptyvbio.com/blog/po102, https://www.adaptyvbio.com/blog/po103 | Round 1 and 2 commentary, PLL length bias | Primary (vendor) |
| S20 | https://github.com/adaptyvbio/nipah_ipsae_pipeline, https://github.com/adaptyvbio/binder-prescreen | Scoring pipeline, prescreen | Primary |
| S21 | https://seqera.io/blog/building-a-scalable-nextflow-protein-design-pipeline-for-adaptyv-s-nipah-binder/ | Local vs official ipSAE gap | Secondary (participant) |
| S22 | https://scour.ing/@minezone/p/https://blog.booleanbiotech.com/adaptyv-binder-design-competition | "1 to 5 designs per designer" | Secondary, UNVERIFIED |
