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
