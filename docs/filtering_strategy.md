# Filtering strategy

This is the single filter specification for the team. It replaces the first draft of 2026-10-23 and applies the changes proposed in `research/SUMMARY.md` after the review of 2026-10-26. Evidence and sources are in `research/01_binder_filtering.md` and `research/04_adaptyv_rules_and_assay.md`. Terms are defined in the glossary at the top of `research/SUMMARY.md`.

Target, metric, length limit and submission cap of the 31 Oct challenge are not announced. Every number below is a starting value taken from earlier Adaptyv events (mostly the TREM2 hackathon, 10 designs per team, Boltz-2 ipSAE ranking). Items marked **(decision)** need the team. Items marked **(unverified)** rest on one secondary source or on an assumption.

## 1. Job of the filter

Generation produces hundreds to thousands of candidates, and we can submit at most K (K = 10 if the TREM2 format repeats). No score ranks binders well: AUROC is 0.6 to 0.8 on the best public data. The filter therefore has three tasks, in this order:

1. Remove designs that fail for reasons we can predict: they do not fold, do not express, copy a known binder, or sit on a patch that is blocked in reality (glycans).
2. Rank the survivors by the announced metric, or by Boltz-2 ipSAE if none is announced.
3. Pick K designs that differ from each other in epitope and in method, because organisers have rewarded difference and ten near-copies fail together.

## 2. Metrics and where they come from

One row per (design, model, seed). "Exists" means it is in `design_pipeline` today (checked 2026-10-26).

| Key | Meaning | Status |
|---|---|---|
| `plddt` | Mean pLDDT of the prediction, 0 to 100 (BindCraft states its gate as 0.8, which is 80 here) | Exists, reported by the adapters |
| `ptm`, `iptm` | Predicted TM score and interface TM score | Exists, reported by the adapters |
| `pae_interaction` | Mean predicted aligned error across the interface, in A (lower is better) | Exists in `scoring.py`, needs the PAE matrix in the scores file |
| `sc_rmsd` | CA RMSD of the prediction against the backbone it came from | Exists, returns None when the CA counts differ, so a complex needs a binder-only variant (to add) |
| `ipsae_ab`, `ipsae_ba`, `ipsae_max`, `ipsae_min`, and `ipsae` | ipSAE for each chain direction, their max and min, and `ipsae`, the rank key (see below) | To add, from Dunbrack's `ipsae.py` or Adaptyv's `nipah_ipsae_pipeline` |
| sequence checks | Length, alphabet, cysteine count, repeats, hydrophobic runs, alanine fraction | To add, sequence only, no GPU |
| `contacts` | Target residues within 4.5 A of the binder in the predicted complex | To add |
| shape complementarity, Rosetta dG | Interface quality terms | To add, optional, report only |

**Which ipSAE is the rank key (unverified).** Dunbrack's README reports ipSAE per ordered chain pair ("asym") and the maximum of the two ("max"). The TREM2 post and the Nipah pipeline README say "ipSAE" and "ipSAE_d0res" and do not say which direction. Proteinbase stores two Boltz-2 columns for the TREM2 designs: `boltz2_ipsae` and `boltz2_min_ipsae`. On the 110 TREM2 rows (re-run 2026-10-26), `boltz2_ipsae` is higher than or equal to the min column in every row, so it is probably the max or one fixed direction, and it separates binders much better: AUROC 0.74 against 0.62 for the min. We therefore use `ipsae` as the name of the rank key and set it to `ipsae_max` until the check below says otherwise. Record all variants. On the day, use whatever the organisers' scoring script outputs.

Check to do before the event: score 10 to 20 public TREM2 designs from Proteinbase with our Boltz-2 and ipSAE setup and compare with their `boltz2_ipsae` column. Agreement within a few hundredths confirms the variant. A systematic gap says our settings (seed, cutoffs, direction) differ.

## 3. Stages

Gates remove a design. Ranks order what is left. A metric that is missing for a design neither passes nor fails it: the gate is skipped and the design gets a flag `skipped:<gate>`. Note that `select(by=...)` in `designs.py` raises a KeyError on a missing score, so gates must be applied before selection and the pool for each `by` must contain only designs that have it.

### Stage 0: format (sequence only)

Keep a design if all hold: 20 natural amino acids only; no His, FLAG, Strep or GFP11 tag, no signal peptide; length within the announced limit, or 50 to 120 residues until one is announced; sequence unique in the pool. Adaptyv appends its own 64-residue tag, so submit the bare binder.

Length: the public EGFR round-2 table has no expression failure under 100 residues (210 designs) and the first failures at 119 residues. The Anthropic campaign expressed 95% at 50 to 120 residues. Prefer 100 residues or fewer where the generator allows it and treat 100 to 120 as acceptable.

### Stage 1: structure confidence

Per design, aggregate over seeds for the primary model (Boltz-2): rank value = median over seeds, gate value = minimum over seeds. In `design_pipeline` terms, a design is the sequence that was folded (`parent_id` of the prediction), the model is `producer`, and the seed is `meta_seed` after `to_dataframe`. Each seed is a separate `predict_structure(..., seed=n)` call, and the seed is part of the cache key, so new seeds really rerun. Adaptyv scores a design once, and local reruns differ with hardware and seed (one team saw 0.79 locally and 0.86 officially, secondary report), so the gate has to hold in every seed we ran. This replaces "best of 3 or more seeds" from the first research summary **(decision, confirm with Kira)**. Run at least 3 seeds.

Gates, all required:

- `ipsae` at least 0.5 in every seed
- `iptm` at least 0.5 in every seed
- `plddt` at least 80 in every seed
- `pae_interaction` at most 10 A, only if an AF2-family model was run
- binder RMSD, design against co-fold, below 3.5 A, only where computable
- shape complementarity at least 0.6, only where computed (not every generator emits it, so some designs will have no value)

Rank by median `ipsae` (or the announced metric). Absolute cuts do not transfer between targets (Overath's per-target optima for AF3 `ipsae_min` run from 0.5 to 0.8), so read the floors as a way to drop the bottom of the pool and read the rest as ranks.

### Stage 2: second model

Refold the survivors with a model outside the generator's own loop: Chai-1 or AF2 for BoltzGen designs, Boltz-2 or Chai-1 for BindCraft designs. Compute `ipsae` there with the same seeds rule.

Cross-model score = mean of the two models' z-scores of median `ipsae`, with each z taken over all designs that passed stage 1 in that model. A design is **all-agree** if both models give median `ipsae` of at least 0.6 (0.61 is Overath's pooled optimum for AF3 `ipsae_min`, applied here to Boltz-2 and Chai-1 by assumption, unverified). At final selection, all-agree designs come first, then the rest by cross-model score. The AND rule raises precision and lowers recall in Overath's data, and with only K slots precision matters more.

### Stage 3: epitope and glycans (manual, 30 to 50 designs)

Take the top 30 to 50 by stage 2 into manual review (Maria and Alixe). For each design compute `contacts`, then check:

- The binder sits on a functional or accessible face of the target, not on a buried patch or on the construct's tag or linker.
- No glycosylated Asn within 5 A of the interface (starting value). The binder is made in E. coli without glycans, while the target is a glycosylated protein from mammalian cells, so a model that ignores glycans overrates patches next to them. HER2 sites per file 02: N68, N124, N187, N259, N530, N571, N629. Reject unless the reviewer overrides, and write the reason down.
- Residue numbers are converted. 1N8Z, 1S78, 3MZW and 3BE1 number HER2 from the mature protein (UniProt minus 22), other entries use UniProt numbering.

Epitope cluster: two designs share an epitope if the Jaccard overlap of their `contacts` sets is at least 0.4 (starting value). Cluster greedily in rank order.

### Stage 4: developability (sequence only, cheap)

Apply the checklist in `research/04_adaptyv_rules_and_assay.md` section 4. In short, flag or drop designs with: an odd cysteine count (a precaution from Adaptyv's docs, the EGFR table did not show an effect), a hydrophobic run of 6 or more, repeats or homopolymer runs, alanine fraction at or above 0.3, high GRAVY, or a low NetSolP or SoluProt score. Prefer acidic, Glu- and Lys-rich surfaces as a weak prior. Use ESM pseudo-log-likelihood only per residue, only to flag the bottom 5 to 10% of the pool. The unnormalised version favours short sequences and was gamed in the EGFR rounds.

Similarity to known binders: Adaptyv has used a 10-residue rule and a 25% edit-distance rule. Build a HER2 reference set (trastuzumab, pertuzumab, affibody ZHER2:342, DARPins, nanobodies, from PLAbDab, Thera-SAbDab and patents) and drop any design below the strictest of the rules announced. `binder-prescreen` ships no HER2 set.

### Stage 5: pick K

Sort by (all-agree first, then cross-model score). Walk down the list and take a design unless its epitope cluster already holds `ceil(K/3)` designs. Reserve 2 or 3 of the K slots for a different method from the main generator, taking the best design of that method even if it ranks lower. Aim for at least 3 epitope clusters.

If fewer than K designs pass every gate, do not pad. Relax one gate at a time in this order: shape complementarity, binder RMSD, `plddt` 80 to 75, `ipsae` 0.5 to 0.4. Record the gate relaxed for each design, and have the reviewers read every design that needed it **(decision)**.

### Stage 6: specificity (report only)

Cofold each of the top designs against the closest paralogs (for HER2: EGFR, HER3, HER4) and report `ipsae` on each. Flag designs whose paralog value reaches their on-target value. It is a column, not a gate, unless the organisers announce a specificity goal (the May FGFR2 track had one).

## 4. What changes when the metric is announced

| Situation | Action |
|---|---|
| Metric is Boltz-2 ipSAE (TREM2 format) | Use the organisers' definition and script as the rank key. Our gates stay. Keep a margin above any official cut. |
| Another metric (iPAE, ipTM, a rank average) | Swap the rank key in stages 1 and 2 for the announced metric, computed with the organisers' tool. Gates stay, because they guard failure modes and do not rank. |
| Nothing announced by the end of the kickoff | Rank by median `ipsae` as above, record every metric, and tell the team we guessed. |
| A scoring script is provided | Run it on 20 designs in the first hour and compare with our local values. Fix any systematic offset before the first selection. |
| Target is not HER2 | Stages 0, 1, 2, 4 and 5 are unchanged. For stage 3, spend the first 30 minutes on the supplied structure: trim to the construct, check the numbering offset, list glycosylation sites from UniProt, choose hotspots. File 02 lists a recipe and six alternative targets. |

## 5. Calibration

- EGFR data (`all_submissions.csv`, 60 binders among 523 labelled tested designs) can calibrate expression and developability weights and the AF2-based metrics (iPAE, ipTM, pLDDT). Its rows carry no Boltz-2 ipSAE.
- Boltz-2 ipSAE can be checked only on the TREM2 and Nipah rows of Proteinbase. Both pools were already cut on ipSAE, which understates AUROC. The TREM2 table has 110 rows and 46 binders, while the organiser's post reports 100 tested and 37 bound, and the extra rows are unexplained (see file 01, section 3).
- With 40 to 60 positives a fitted weight is noisy. Default: fixed gates plus rank by the announced metric, with no fitted weights. Decide by 29 Oct whether a fit is worth a day **(decision)**.

## 6. Where this lives in the code

- Metrics go into `METRICS` in `design_pipeline/scoring.py`, one function of (prediction, reference) that returns None when it does not apply.
- Stages 0 to 2 and 5 are table operations over `to_dataframe(designs)`. `filters/` is empty today and is the proposed home **(decision, Kira)**.
- `select()` sorts ascending by default. Use `ascending=False` for ipSAE, ipTM and pLDDT, and `ascending=True` for `pae_interaction` and RMSD.

## 7. Open decisions

1. Median rank with minimum gate across seeds instead of best of seeds (stage 1).
2. Whether to fit weights on public data (section 5), by 29 Oct.
3. Where the filter code lives, and the names of the ipSAE columns. The rank-key variant is settled by the TREM2 check in section 2.
4. The gate relaxation order and the epitope similarity cutoff (0.4), both starting values.
5. Who builds the HER2 known-binder reference set, and by when.

## What is already implemented

Written as of 2026-10-30, in `design_pipeline` rather than in `filters/`, because the pipeline turned out to be affinity maturation from a known complex rather than de novo generation:

- **Stage 4, developability** is the `gate` step. `liabilities` is pure sequence regex (new N-glycosylation sequon, deamidation, isomerization, fragmentation, introduced cysteine) and costs nothing; `netsolp`, `thermompnn` and `aggrescan3d` are container checks. Failures are dropped, never weighted against binding score, and the annotated table keeps every reason so a threshold can be re-tuned without re-running anything.
- **Stage 1, the metric cut** is `score_variants`: one `score_<model>` column per scorer, joined on `variant_id`, never averaged. The open decision below — hard thresholds versus a weighted score — is deliberately left open in the code: there is no combined column and no default weights.
- **Stage 2, cross-model agreement** is what several scorer columns on one table give you. Two inverse-folding models plus a supervised ddG model is three opinions on the same candidate.
- **Calibration** is `design_pipeline.harness`: `run_systems` + `report` give Spearman and precision@10 per scorer per system against published affinity data (FLAb, AbBiBench, AB-Bind). That is where the stage weights should come from, rather than from a guess. It is the answer to this document's calibration section, and it works on measured antibody systems whether or not the Adaptyv data turns out to be downloadable.
- **Stage 5, diversity** is partly `select_pairs`, which spreads pairs across the top of the ranking and enforces a minimum CA-CA separation rather than returning 20 variations on one position.

Stages 3 (epitope) and 6 (specificity) are still manual, as this document says they should be.
