# Filtering strategy (draft)

The generation step will produce far more candidates than we can submit, so the filter decides the outcome. Which stages carry the most weight depends on the judging metric, which we do not know yet.

## Funnel

1. **Metric cut.** Use whatever the judges announce (ipSAE, ipTM or iPAE), plus binder pLDDT. Until then, record all of them for every candidate.
2. **Cross-model agreement.** Re-predict each candidate complex with a second co-folding model (Boltz or Chai). Keep designs where both models place the binder in the same pose.
3. **Epitope check.** The binder should sit on a functional epitope, not on a glycosylated or buried patch. This stage needs biological judgment (Maria, Alixe).
4. **Developability.**
   - no unpaired cysteines
   - no large hydrophobic patches
   - sensible pI and net charge
   - no low-complexity or repeated sequence
   - within the length limit
   - ESM pseudo-log-likelihood as an expressibility proxy
5. **Diversity.** Cluster the survivors and submit one or two per cluster. Ten near-copies of one design waste the submission cap, and novelty of approach may be scored.
6. **Specificity.** Where the target has paralogues (HER2, EGFR, HER3), check that the binder does not sit on a conserved face.

## Calibration

If the Adaptyv EGFR competition data is public, we can check which in-silico metrics predicted binding in the lab and set stage weights from that. About 2.5% of validated designs bound in round 1, so the signal is thin and we should expect noisy weights. Not yet checked whether the dataset is downloadable.

## Open decisions

- Hard thresholds versus a weighted score. Hard thresholds are easier to explain to a jury; a weighted score uses the metrics better.
- How many candidates to carry into stage 3, since epitope review is manual and slow.
- Whether to spend Iason's physics-based rescoring (FoldX, Rosetta, short MD) on all survivors or only the top 20.

## Code layout

Each stage lives in `filters/` as one function that takes a table of candidates and returns the same table with a pass flag and the metric values. Stages can then be reordered or dropped once the real metric is known.

## What is already implemented

Written as of 2026-10-30, in `design_pipeline` rather than in `filters/`, because the pipeline turned out to be affinity maturation from a known complex rather than de novo generation:

- **Stage 4, developability** is the `gate` step. `liabilities` is pure sequence regex (new N-glycosylation sequon, deamidation, isomerization, fragmentation, introduced cysteine) and costs nothing; `netsolp`, `thermompnn` and `aggrescan3d` are container checks. Failures are dropped, never weighted against binding score, and the annotated table keeps every reason so a threshold can be re-tuned without re-running anything.
- **Stage 1, the metric cut** is `score_variants`: one `score_<model>` column per scorer, joined on `variant_id`, never averaged. The open decision below — hard thresholds versus a weighted score — is deliberately left open in the code: there is no combined column and no default weights.
- **Stage 2, cross-model agreement** is what several scorer columns on one table give you. Two inverse-folding models plus a supervised ddG model is three opinions on the same candidate.
- **Calibration** is `design_pipeline.harness`: `run_systems` + `report` give Spearman and precision@10 per scorer per system against published affinity data (FLAb, AbBiBench, AB-Bind). That is where the stage weights should come from, rather than from a guess. It is the answer to this document's calibration section, and it works on measured antibody systems whether or not the Adaptyv data turns out to be downloadable.
- **Stage 5, diversity** is partly `select_pairs`, which spreads pairs across the top of the ranking and enforces a minimum CA-CA separation rather than returning 20 variations on one position.

Stages 3 (epitope) and 6 (specificity) are still manual, as this document says they should be.
