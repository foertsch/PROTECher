# Research summary (2026-10-24, corrected 2026-10-26 after review)

Details, numbers and URLs are in files 01 to 05 in this folder. The filter itself is specified in `../filtering_strategy.md`, which is the single source for it. Target, metric, length limit and compute are unannounced (the Luma page names only "a protein design challenge by Adaptyv"). Everything that depends on them is inference from earlier Adaptyv events.

## Glossary

- **Co-folding model:** predicts the binder and the target together as one complex (Boltz-2, Chai-1, AlphaFold). Gives the scores below.
- **pLDDT:** per-residue confidence of a prediction, 0 to 100. Used as a mean per design.
- **ipTM, pTM:** interface and overall TM-score estimates of a predicted complex, 0 to 1, higher is better.
- **iPAE (pae_interaction):** mean predicted aligned error across the binder-target interface, in A, lower is better.
- **ipSAE:** an interface score from Dunbrack that is computed per ordered chain pair from the PAE matrix and counts only residues whose predicted error to the partner is below a cutoff. 0 to 1, higher is better. Variants differ in how the length scale d0 is set (the `d0res` variant sets it from the number of residues passing the cutoff) and in direction: A to B, B to A, their max, or their min (min is Overath's convention). Which variant Adaptyv ranks on is not documented (file 01, section 10).
- **PAE cutoff, distance cutoff:** the ipSAE settings used in the TREM2 hackathon, both 15 A.
- **AUROC:** the chance that a random binder scores higher than a random non-binder. 0.5 is chance, 1.0 is perfect.
- **AP:** average precision, the area under the precision-recall curve. Chance level equals the share of binders.
- **PLL:** ESM2 pseudo-log-likelihood, a language-model estimate of how natural a sequence looks. "Per residue" divides by length.
- **Rec:** a number recomputed by the research session from a public table. It is not a published value.
- **KD:** dissociation constant, lower is tighter binding. **BLI, SPR:** the two binding assays Adaptyv uses.
- **Hit rate:** share of tested designs that bind.
- **Gate and rank:** a gate removes designs below a threshold. A rank orders what is left.

## Top recommendations

1. **Plan for Boltz-2 ipSAE as the ranking metric.** The one-day TREM2 hackathon (Feb 2026, 16 teams) allowed 10 designs per team, ranked by Boltz-2 ipSAE_d0res (PAE and distance cutoff 15) and sent the top 100 of 141 to the lab (89 expressed, 37 bound). Script this now with Adaptyv's pipeline (github.com/adaptyvbio/nipah_ipsae_pipeline) and run 3 or more seeds. Rank on the median over seeds and apply floors to the minimum over seeds, because Adaptyv scores each design once and one participant saw 0.79 locally against 0.86 officially (secondary report). Before the event, rescore public TREM2 designs from Proteinbase and compare with their `boltz2_ipsae` column to learn which ipSAE variant Adaptyv uses.
2. **Treat scores as ranks.** ipSAE is the best single score in the comparisons we found (AUROC 0.63 to 0.83 depending on data set, recomputed by the research session from public tables, not published values), and no metric predicts binding well. De novo hit rates on new receptors ran from about 1% to 10%, depending on how designs are categorised (file 01, section 2).
3. **Design for expression.** In the public EGFR round-2 table, no design under 100 residues failed to express (210 designs) and the first failures are at 119 residues. The Anthropic campaign expressed 95% at 50 to 120 residues. Prefer 100 residues or fewer, allow up to 120. Use SolubleMPNN, no unpaired cysteines (an Adaptyv docs precaution, not an EGFR finding), no repeats, and length-normalised ESM likelihood. Adaptyv appends a 64-residue tag, so submit untagged sequences.
4. **Generate with BoltzGen first and BindCraft second.** Both are MIT-licensed and have wet-lab records on receptors. Download weights and build containers this week, and dry-run on HER2 and EGFR on the real hardware before 30 Oct.
5. **Prepare HER2 now.** Convert numbering (1N8Z uses mature numbering, UniProt is +22), pick two or three sites, and build a known-binder reference set, because Adaptyv has used 10-residue and 25%-edit-distance similarity rules.
6. **Submit a spread.** Organisers reward difference (Adaptyv hand-picks for diversity and method novelty, Berlin picked finalists for workflow novelty). Cover several epitopes and methods, with 2 or 3 slots for a different approach.
7. **Calibrate with care.** EGFR rounds 1 and 2 (60 binders in 523 labelled designs) can calibrate expression and AF2-based metrics but carry no Boltz-2 ipSAE. Only the TREM2 and Nipah rows of Proteinbase carry it. The Overath et al. table has 3,669 labelled designs on 15 targets (the paper reports 3,766). No HER2 data exist. Decide by 29 Oct whether a fit on 40 to 60 positives is worth the time.

## Likely targets (tiers, file 02)

The first version of this table gave percentages. They had no stated method and were dropped.

| Tier | Target | Why |
|---|---|---|
| Likely | HER2 | Iconic breast cancer target, many ECD complexes among 62 RCSB entries, trastuzumab, pertuzumab, affibody, DARPin and nanobody controls. 1N8Z and 1S78 verified. Not in past Adaptyv rounds, but the same "new target" argument applies to TROP2, Nectin-4 and HER3 |
| Possible | EGFR | Adaptyv's best-characterised target, live in the Anthropic x Adaptyv series |
| Possible | TROP2 | Approved TNBC drug target, 1.56 A Fab structure (9PI9) |
| Long tail | Nectin-4, HER3, PD-L1, FGFR2 | Nectin-4 has few glycans. HER3 is hard to design against. PD-L1 is easiest but has a weak breast link. HerCode's May 2026 track used FGFR2 |

## Filtering

The operational filter (stages, thresholds, missing-value rule, what to do when fewer than K designs pass, what changes when the metric is announced) is in `../filtering_strategy.md`. The changes proposed here earlier have been applied there, with these corrections from the review: median over seeds with a minimum gate, instead of best of seeds; the rank key is Adaptyv's `boltz2_ipsae` column and not the min over both directions; and calibration is split by what each data set can support.

## Still open (needs the team, no owner yet)

1. A day-by-day plan for 26 to 30 Oct with a named owner per item.
2. The GPUs available to the team on the day. File 03's throughput budget depends on it.
3. The HER2 known-binder reference set (half a day, a disqualification risk if skipped).
4. A dry run of the submission mechanics: channel, method note, format check, duplicate check.
5. A rule for when to switch to hosted compute (Tamarind, Rowan, PXDesign server) if local compute fails on the day.
6. Expected funnel numbers worked out in advance: designs generated, surviving each stage, and the plan if only three survive.
7. A pitch owner. File 05, section 3 has a 3-minute pitch checklist.

## Questions for the organisers

1. Which target, and which construct or PDB chain is supplied?
2. Which ranking metric, computed by which model and cutoffs, per design or per submission set? Is a scoring script provided?
3. How many designs per team, and how many go to the lab (about 100 overall in TREM2)?
4. Maximum length, and does it include the 64-residue tag Adaptyv appends? How are tags and the start methionine handled?
5. Is wet-lab testing guaranteed, for how many designs, and when do results arrive?
6. Which similarity rule to known binders applies?
7. Which tools are allowed, and what compute is provided (GPU type, credits, internet access)?
8. How is the jury scored (metric, novelty, biology, pitch), and what is the pitch format?
9. How does this event relate to Swiss {ai} Weeks and to the Anthropic x Adaptyv competition (challenge 5 opens 26 Oct)? One press report gives Swiss {ai} Weeks 2026 as 1 Sep to 4 Oct, which conflicts with the Luma page.

Unverified: the target, metric and all rules of this event. Each file ends with its own list of unverified items.
