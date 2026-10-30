# Track A filtering plan

Status: agreed with Arion on 2026-10-30, not yet reviewed by the rest of the team. It covers Track A only (improve an announced scFv with at least 5 substitutions, 15 designs, one blind round). `filtering_strategy.md` still holds the de novo plan and is left unchanged. Terms are in `track_a_glossary.md`. Evidence is in `research/04`, `06` and `07`, and in Adaptyv's po101 blog post.

## What the filter is for

No score we can compute predicts KD. In Adaptyv's public data ipTM, pLDDT and ESM likelihood sit near AUROC 0.5 to 0.6 for affinity (07). The funnel therefore removes likely losers and builds designs, and the selection spreads risk across mechanisms. The jury question "why these 15" is answered by one line per design: its mutations, their source and the effect we predict.

## Funnel

| # | Stage | Rule or threshold | Basis |
|---|---|---|---|
| 1 | Rules | Substitutions only, at least 5 per design, all inside VH or VL. Freeze the conserved cysteine pairs, the linker, the VH-first orientation and the C-terminal end. No new cysteine, no indels, no tag in the sequence. Length unchanged. | Rule (Track A text, 04, 06) |
| 2 | Position map | One class per residue, in this order: paratope (any heavy atom within 5 Å of the antigen), other CDR (approximate Kabat), second shell (framework within 8 Å of the antigen or 4.5 Å of a paratope residue), VH-VL interface (within 4.5 Å of the other domain), framework. Structure is the crystal if the antibody is known (1N8Z for trastuzumab), otherwise a Boltz-2 prediction. | Structural fact. On 1N8Z: 21 paratope, 40 other CDR, 12 second shell, 17 interface, 137 framework. |
| 3 | Single scan | Score all single substitutions (about 4,200: 223 positions times 19) with ESM-1b, five ESM-1v models, ESM-2, and ESM-IF1, AntiFold and ProteinMPNN with the antigen in the input. Keep a substitution if its likelihood ratio exceeds 1 in at least 2 of the 6 sequence models, one model ranks its residue first at that site, and the inverse-folding ratio with the antigen exceeds 1. A scorer leaves the consensus if its AUROC for improved against not improved on the 57 cetuximab scFv variants is at or below 0.5. That test is too small to pick a best scorer. | Hie 2024 rule (k ranged 1 to 4, Fab and IgG). Shanker 2024 (ESM-IF1, IgG). The antigen condition and the drop rule are our judgment. Aim for about 300 survivors (06). |
| 4 | Liabilities and energy | Drop new unpaired cysteines, new NG and DG motifs and new exposed hydrophobic patches. Rank the rest by ΔΔG from FoldX or Rosetta, electrostatic term first. Keep about 60 tagged mutations. The ΔΔG cut is set after we see the score distribution. | Judgment. ΔΔG agreement is weak to moderate (AB-Bind r 0.16 to 0.45, 06). |
| 5 | Build sets | Draw sets of 5 to 10 from the shortlist and score each set jointly with a language model. For each mutation count from 5 to 10 keep the best-scoring sets that meet the spread rules (a proximal frontier, as in PEX). Spread rules: each mutation appears in 2 to 4 designs, no CDR-H3 mutation in more than 2 designs, designs differ by at least 3 substitutions, no two mutations within 4 Å in one design unless scored together, isoelectric point shift under 1. | Judgment. Combinations interact (Lippow: 4-fold where 8.4-fold was expected, 06). The PEX frontier comes from Adaptyv's po101 survey (GFP and AAV data, no antibodies). |
| 6 | Structure veto | Fold the parent and the finalists with Boltz-2 over at least 5 seeds. Reject a variant if ipTM or ipSAE falls below the parent's mean minus 2 SD of seed noise, or if the paratope Cα atoms move more than 1.5 Å against the parent's prediction after superposing the antigen. The veto never ranks. Calibrate first on about 40 damaging and 40 framework variants. If the damaging variants stay inside the parent's noise band, the veto is advisory and the pitch says so. | Thresholds are 06's own, not published. Scores did not rank affinity in public data (07). |
| 7 | Pick the 15 | Fill the arms below, then write the one-line rationale per design. | Judgment |

## How the 15 are picked

| Arm | Designs (default) | Substitutions | Where | Basis |
|---|---|---|---|---|
| Consensus framework | 8 | 5 to 10 | Framework and second shell where language models and germline agree | Cradle's cetuximab winner changed 10 framework positions (07) |
| Stability and interface | 4 | 5 to 7 | VH-VL interface and CDR-adjacent surface, including VH K30D | 06, from Dudgeon 2012. VL S52D sits in CDR-L2 and would use a loop slot. |
| Loop periphery | 3 | 1 to 2 CDR changes plus 3 to 4 framework changes | Peripheral CDR positions, never a CDR-H3 hot spot | 06, 07 |

Mutations per design: 5 to 6 for most, 7 to 8 for some, 9 to 10 for two, with the sets taken from the frontier of stage 5. The frontier plot (joint score against mutation count) goes on a slide for the jury.

Loop designs are a dial, set on the day from the announced parent. A parent below 1 nM that expresses well has little framework headroom, so use 5 or 6 and shrink the consensus arm. A weak or poorly expressing parent gets 2 or 3. If expression is not announced, judge from KD alone or keep 3. The stability arm stays at 4.

Selection is greedy: highest joint score first, subject to the spread rules, tie broken toward fewer CDR changes. Each design's line carries every mutation with its class and source (literature, model consensus or structure) and the predicted change in kon, koff, KD or expression.

Expectation for the pitch (07, inference): about half of framework-only variants bind and few beat a strong parent. In the cetuximab data, Gleghorn's loop-light variants gave 10 binders and 5 improved out of 17, and variants of Cradle's winner improved in 1 of 18.

## Precomputed today

Done: a 227-row position table for 1N8Z (scratch copy, not in the repo yet), the recheck of 07's CDR split from the Proteinbase table, and the check of Hie's rule against the paper.

To run, in priority order:
1. Score the 57 cetuximab scFv variants (public Adaptyv data) with every scorer, and see which separates the 18 improved from the rest. This is the only same-format calibration we have.
2. Run stages 2 and 3 for trastuzumab (1N8Z), pertuzumab (1S78), sacituzumab (9PI9) and cetuximab (1YY9).
3. Calibrate the veto: parent over 5 seeds plus about 80 variants, about 3 hours on one L40S.
4. Build a full 15 for trastuzumab. On the day the same script reruns on the announced sequence, because the reformatted scFv may differ from the crystal.

Code to write: apply and count mutations on a parent, ESM likelihood scanner, inverse-folding scorer with antigen, position classifier (the 1N8Z script is the start), set builder with the spread rules and the loop-count parameter, veto script.

## Open decisions

1. Second-shell definition. The antigen-distance rule in 06 gives 2 residues on 1N8Z and ours gives 12. Either way the framework arm rests mainly on sequence signals.
2. FoldX or Rosetta licence. Without one, stage 4 uses liability filters and inverse-folding scores only.
3. A second structure model. Chai-1 gave ipTM 0.16 to 0.18 on the trastuzumab input where Boltz-2 gave 0.94. Use ColabFold or Protenix for a second opinion, or none.
4. Scorer weights. Equal under the consensus rule. A scorer drops out if it is at or below chance on the cetuximab test, and a clear winner there may earn a higher weight.
5. Veto thresholds and a floor for the SD, after calibration.
6. Who writes which scorer today.
7. Mutation-count quota. The frontier shows what each extra mutation costs in joint score. Keep the fixed quota above, or set it from where the frontier flattens.

## Questions for the organisers

Which orientation and linker, and is a start Met present? Do linker mutations count toward the 5? Are novelty rules applied to variants? Is the parent's expression announced with its KD? Do lab results reach the jury, and when?
