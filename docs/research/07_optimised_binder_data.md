# 07. Optimised and mutated binders: what the public Adaptyv data show

Research date 2026-10-30. "Rec" marks a number recomputed in this session from the public tables in section 8. "Pub" marks a published value. "UNVERIFIED" marks anything not opened or not reproducible. Counts depend on the thresholds in section 5.

## Bottom line

1. Only one public parent has our format: the cetuximab scFv (241 aa, VH-(G4S)3-VL). Rec: 57 variants, 43 tested, 27 expressed, 24 bound, 18 with a lower KD than the parent (6.64 nM, the mean of replicates at 3.33 and 9.94 nM). Pub: the paper and Cradle quote 9.94 nM.
2. Mutation site decided the outcome. Rec: 13 of 13 tested variants without CDR mutations bound and beat the parent. Variants with 3 or more CDR mutations bound in 2 of 13 (1.8 nM and 10 µM). Only Cradle and the Gleghorn lab made CDR-free variants, so method and region are confounded.
3. Ten framework mutations were enough. Rec: Cradle's 12 sequences (10 mutations each, 25 positions, no CDR change) all expressed and bound at 0.37 to 4.9 nM. Pub: the winner was 1.21 nM, 8.2-fold better than the control.
4. No score separated improved from non-improved variants. Rec AUROC for binding on 107 EGFR variants (37 bound): ipTM 0.59, pLDDT 0.70, ESM2 per-residue PLL 0.59. Boltz-2 ipSAE gave 0.59 and 0.55 on 450 ProteinMPNN redesigns of AlphaProteo binders (PD-L1, VEGF) and 0.51 and 0.53 when ranking the binders by KD.
5. Expect about half of 15 framework-only variants to bind and few to beat a strong parent. Rec: 3 to 11 extra mutations on Cradle's winner gave 1 improved of 18 (parent by inference). An 80 pM PD-L1 binder lost affinity in 20 of 20 pLM variants.

## 1. Parents and variant counts (Rec)

EGFR universe: 1,857 submissions plus 205 Proteinbase-only designs, 2,062 unique sequences; 834 tested, 737 expressed, 129 bound.

| Parent (reference source) | Variants | Tested | Expressed | Bound | Authors or methods |
|---|---|---|---|---|---|
| Cetuximab scFv (control row) | 57 | 43 | 27 | 24 | Cradle 12, Gleghorn lab 20, nhaas 10, Eureka 5, junhee.park 5, corrochano 4, kerrding 1 |
| Human EGF (control row) | 100 | 89 | 87 | 48 | many; ProtRL is pLM-generated, so "variant" is loose |
| Human TGFa (UniProt) | 32 | 23 | 23 | 4 | Aurelia Bustos and others |
| BTC, EREG, EPGN (UniProt) | 8 | 5 | 5 | 1 | |
| EGFR nanobodies 7D12, EgA1 (PDB), at most 20 mutations | 32 | 14 | 9 | 5 | ahmedsameh, chrisxushaoyong, cchallacombe, tim |
| scFv near panitumumab (PDB) | 11 | 2 | 2 | 0 | julian.faust crossbreds |
| **All known-parent variants** | **240** | **176** | **153** | **82** | |
| All other EGFR designs | 1,821 | 657 | 583 | 46 | |

With competition-time labels only, variants bound 28 of 98 (29%) and others 32 of 503 (6%). Follow-up measurements in Proteinbase (Cradle, Gleghorn lab, ProtRL) lift the variant rate to 47%. Designer categories (Rec): optimised 13 of 88 bound, diversified 32 of 184. My references explain only 39 of 135 optimised and 67 of 382 diversified submissions.

Non-EGFR parent sets (human targets only):
- DeepMind SynthID arm: 15 AlphaProteo parents per target with 255 ProteinMPNN redesigns each on PD-L1 and VEGF110 (3 of 15 VEGF parents do not bind). The RBD arm is aggregate only: 255 variants, 250 expressed, 213 bound.
- Gleghorn lab "dsm-synteract" on PD-L1: 20 variants nearest to one AlphaProteo parent (80 pM), 19 of them 6 to 29 mutations away.
- Pro-1 on FGF-1 (label fgf-r1): wild type (15.2 nM) plus 17 point-mutant variants (1 to 7 mutations) and 2 tandem fusions.

## 2. Outcomes versus the parent (Rec)

Improved: KD below the parent's. Matched: up to 3x weaker. Weaker: over 3x. Lost: expressed, no binding.

| Set (parent KD) | n | Expressed | Bound | Improved | Matched | Weaker | Lost | Not expressed |
|---|---|---|---|---|---|---|---|---|
| Cetuximab scFv (6.64 nM) | 43 | 27 | 24 | 18 | 0 | 6 | 3 | 16 |
| Cradle | 12 | 12 | 12 | 12 | 0 | 0 | 0 | 0 |
| Gleghorn lab | 20 | 12 | 11 | 6 | 0 | 5 | 1 | 8 |
| Eureka, corrochano, nhaas | 11 | 3 | 1 | 0 | 0 | 1 | 2 | 8 |
| EGF-derived (759 nM) | 88 | 86 | 47 | 33 | 8 | 6 | 39 | 2 |
| FGF-1 Pro-1 (15.2 nM) | 17 | 14 | 5 | 1 | 2 | 2 | 9 | 3 |
| PD-L1 ProteinMPNN | 252 | 252 | 175 | 70 | 32 | 68 | 77 | 0 |
| VEGF ProteinMPNN (binding parents) | 202 | 201 | 119 | 47 | 21 | 50 | 82 | 1 |
| PD-L1 Gleghorn lab pLM (80 pM) | 20 | 20 | 15 | 0 | 0 | 15 | 5 | 0 |

Against 9.94 nM the cetuximab row has 18 improved, 2 matched, 4 weaker. Fifteen of the 24 binders are at least 2x better than 6.64 nM. The control expressed "low" and Cradle reports high expression for the winner, so part of the gain may be more active protein (inference).

| Mutations | Cetuximab scFv: tested / bound / improved | Pooled parent-measured sets: n / bound / improved |
|---|---|---|
| 1 to 2 | none | 2 / 2 / 1 |
| 3 to 5 | none | 11 / 2 / 0 |
| 6 to 10 | 13 / 12 / 12 | 23 / 19 / 12 |
| 11 to 20 | 28 / 12 / 6 | 65 / 38 / 12 |
| 21 to 40 | none | 215 / 137 / 55 |
| over 40 | 2 / 0 / 0 | 218 / 140 / 56 |

The EGFR rules required at least ten mutations from known binders (Pub, proteinbase.com/competitions), so no cetuximab variant has fewer than 9. The 1 to 5 rows come from 13 Pro-1 point mutants and the over-20 rows from ProteinMPNN redesigns, so the pooled column is a loose prior only.

| CDR mutations (cetuximab variants) | Tested | Bound | Improved |
|---|---|---|---|
| 0 | 13 | 13 | 13 |
| 1 to 2 | 15 | 9 | 4 |
| 3 to 5 | 5 | 1 | 1 |
| over 5 | 8 | 1 | 0 |

Nanobodies (mutations from the nearest of 7D12, EgA1, EgB4, 34E5): 0 to 5 mutations 4 of 7 tested bound, 11 to 20 mutations 1 of 7, 21 to 45 mutations 6 of 38. No nanobody parent KD was measured, so no improved share exists.

## 3. Which scores predicted anything (Rec)

AUROC with bootstrap 95% interval. Non-expressers count as negatives for EGFR. "Variants" are the 107 tested known-parent EGFR variants with AF2-multimer scores (37 bound).

| Score | EGFR variants | All 629 tested EGFR (70 bound) | PD-L1 redesigns (249) | VEGF redesigns (201) |
|---|---|---|---|---|
| ipTM | 0.59 (0.47 to 0.70) | 0.64 (0.57 to 0.69) | 0.57 (Boltz-2) | 0.55 (Boltz-2) |
| pLDDT | 0.70 (0.60 to 0.79) | 0.69 (0.64 to 0.76) | 0.54 | 0.56 |
| iPAE (lower better) | 0.59 (0.47 to 0.69) | 0.63 (0.57 to 0.69) | not applicable | not applicable |
| Boltz-2 ipSAE | not available | not available | 0.59 (0.51 to 0.67) | 0.55 (0.47 to 0.63) |
| ESM2 PLL, unnormalised | 0.49 (0.38 to 0.60) | 0.49 (0.41 to 0.56) | not run | not run |
| ESM2 PLL per residue | 0.59 (0.48 to 0.69) | 0.71 (0.65 to 0.76) | not run | not run |
| Rosetta interface dG | 0.51 (0.40 to 0.62) | 0.64 (0.58 to 0.72) | not run | not run |
| Fewer mutations | 0.64 (0.52 to 0.76) | not applicable | 0.45 (0.38 to 0.53) | 0.47 (0.38 to 0.55) |

- Per-residue PLL and dG lose most of their signal inside variants, since both track similarity to known binders.
- Among SynthID binders, ipSAE did not rank by KD (0.51 on PD-L1, 70 improved vs 100 not; 0.53 on VEGF). The change in ipSAE from the parent separated binders from non-binders at 0.67 (0.59 to 0.73) on PD-L1 and 0.57 on VEGF, and improved from non-improved at 0.44 and 0.53.
- Cetuximab variants alone (21 scored, 11 bound): iPAE 0.89, fewer mutations 0.92. Cradle's 10 designs drive both, so they are too small to use.
- Fewer CDR mutations: AUROC 0.88 (0.77 to 0.96) for binding and 0.89 (0.74 to 0.99) for improved vs not among 24 binders; 0.75 (0.58 to 0.90) without Cradle (29 variants).
- DE-STRESS (81 round-2 variants, 25 bound, 62 columns): the largest deviation (rosetta_fa_sol_per_aa, 0.16) exceeds a label-shuffle null (95th percentile 0.22). It probably tracks family composition, which I did not test. Nanobody scores (12 tested, 3 bound) are too small to read.

## 4. Who did well and how

| Entry | Declared method | Mutations (Rec) | Outcome |
|---|---|---|---|
| Cradle, round 2 winner | Zero-shot ML optimisation of the scFv with CDRs blocked (Pub: Cradle blog, po104). The paper cites evolutionary data (MMseqs2 on uniref30 and colabfold_envdb) | 10 per design, 33 unique at 25 positions, all framework: heavy FR3 (70 to 88), light FR1 to FR3 (137 to 156, 173 to 183, 210 to 214); up to 19 apart | 1 of 10 tested in the competition (1.21 nM). Pub: 11 follow-ups all would have won, 3 sub-nM, best 339 pM. Rec: 12 of 12 bound, 0.37 to 4.9 nM |
| Gleghorn lab (Synthyra) "dsm-synteract" | Protein language model (details UNVERIFIED) | 12 to 20 from cetuximab, 3 to 11 from Cradle's winner; all 11 binders carry at least 9 of Cradle's 10 mutations (inference: the public winner was the parent) | 11 of 20 bound, best 0.86 nM. Of 18 within 11 mutations of the winner: 1 improved, 5 within 3x, 5 weaker, 7 not expressed or lost |
| Chris Xu (Tihumab), nanobody, rank 2 | CDR grafting onto human germline frameworks (Pub: po104) | 35 edits from EgA1; 7D12-like set 12 to 15 | 5.18 and 18.0 nM; 7D12-like set 1 of 3 bound (238 nM) |
| jvogt4, rank 4 | AntiBERTy proposals ranked with DeepAB (Pub: paper) | Parent unidentified; 50 to 55 mismatches to nearest patent hit | 5 of 5 scFvs bound: 80, 149, 151, 401 nM, 10 µM |
| Aurelia Bustos | TGFa and EGF redesign | 10 to 25 edits | 8 of 8 tested bound, 51.6 nM to 2.9 µM |
| ProtRL (Ferruz lab) | RL-tuned pLM, EGF-like | best at 7 from EGF | 35 of 87 bound, best 11.6 nM |

K5Q, N70S, N73T, N88D, R179K and K183R occur in all 11 Gleghorn lab binders. The Cradle follow-up post (Pub) says the best binders were not the closest to the winner.

Anthropic x Adaptyv (proteinbase.com): five single-design challenges (EGFR, TNF-a, GLP1R, TfR1, KRAS G12D pMHC). None names an antibody, scFv or optimisation track, no results are posted, and selection "will not rely on a single in silico metric".

## 5. Method note

- References: cetuximab scFv and human EGF from the Adaptyv control rows; EGFR ligands from UniProt; 7D12, EgA1, EgB4, 34E5 and the VH/VL of cetuximab, necitumumab, panitumumab, nimotuzumab, matuzumab from RCSB FASTA.
- Alignment: own standard-library Python (banded global edit distance, substitution and gap cost 1), run with `python3 -I`. A variant is at most 30 edits (scFv), 25 (ligands) or 20 mutations (nanobodies). Mutation count is substitutions plus internal indels; terminal tags and a 6-residue linker extension are not counted. Check: my alignment recovers the 10 mutations in the Cradle winner's name, and the 12 Cradle sequences give 33 unique mutations at 25 positions with a maximum pairwise distance of 19, as in the Cradle follow-up post.
- CDRs: ANARCI is not installed, so approximate Kabat by motif (H1 from Cys22+9 to before WVRQ, H2 from W+14, H3 between CAR and WGQG, L1 after Cys23, L2 seven residues from W+15, L3 between the last Cys and FGxG). For cetuximab the H3 (ALTYYDYEFAY) matches paper Table 3; the other five CDRs were not checked against a source.
- Labels: competition label if tested then, else Proteinbase; Proteinbase KD is the arithmetic mean of replicates. Expressed means the competition label is not "Unknown", or any positive split-GFP call in Proteinbase.
- SynthID: each variant goes to its nearest parent and is kept if the second-nearest is at least 5 edits further (3 PD-L1 and 2 VEGF variants dropped).
- AUROC: Mann-Whitney, 500 or 1,000 bootstrap resamples, sign flipped for iPAE, dG and mutation count.

## 6. What could not be verified

- The hackathon antibody is unknown; no scFv other than cetuximab is covered.
- The Gleghorn lab method, the SynthID design and Pro-1 were not read in primary papers. Pro-1's FGF-1 target comes from a search summary (https://adaptyvbio.substack.com/p/protein-designer-spotlight-can-a). The DSM link to Cradle's winner rests on shared mutations only.
- Cetuximab's KD is 3.33 or 9.94 nM by replicate, so "improved" below 3x is within noise. Cradle's 339 pM is Pub only (my best mean is 0.37 nM).
- Most declared "optimised" submissions have no identified parent here.
- No Boltz-2 scores exist for EGFR in Proteinbase, so ipSAE on scFv variants is untested.
- 6 labels differ between all_submissions.csv and Proteinbase (I used the competition label); 8 tested designs are missing from Proteinbase.
- Nipah antibody-format entries (aggregate only, crude motif detector): scFv-like 5 of 25 bound, VHH-like 1 of 96; parents unknown.
- The Gleghorn lab, Pro-1 and PD-L1 pLM rows have about 20 variants each, too few for intervals.

## 7. Implication for 15 variants with at least 5 mutations each

Put at least ten slots on framework-only variants with 5 to 10 mutations: 13 of 13 CDR-free cetuximab variants bound and beat the parent, while 2 of 13 with 3 or more CDR mutations bound. Choose positions where the parent differs from its family or germline consensus, as Cradle did, and keep the CDRs intact; spend at most two or three slots on 1 to 2 CDR mutations, which bound in 9 of 15 and were over 3x weaker in 5 of those 9. Do not rank by ipTM, pLDDT, PLL or ipSAE for affinity (AUROC 0.5 to 0.6, no signal for KD order); use them to drop variants predicted to unfold, and rank by consensus match and position diversity. Binding near 50% is plausible, with an improved KD in roughly one of five (6 of 31 non-Cradle cetuximab variants), and much lower if the parent is already strong and well expressed (1 of 18 on Cradle's winner, 0 of 20 on an 80 pM PD-L1 parent). Gains may be largest if the announced scFv expresses poorly, as cetuximab did (inference). The evidence is one antibody on EGFR plus a few miniprotein targets, so treat these as priors.

## 8. Sources

- https://storage.proteinbase.com/proteinbase_latest.csv (46,310,408 bytes, downloaded 2026-10-30, latest measurement 2026-10-10)
- https://media.githubusercontent.com/media/adaptyvbio/egfr2024_post_competition/main/adaptyv/data/all_submissions.csv
- https://github.com/adaptyvbio/egfr_competition_1 and https://github.com/adaptyvbio/egfr_competition_2 (result_summary, replicate_summary, destress_binder_only)
- https://www.biorxiv.org/content/10.1101/2025.04.17.648362.full.pdf (sections 3.2, 4.2, 4.3, Table 3)
- https://www.cradle.bio/blog/adaptyv-protein-design-competition and https://www.cradle.bio/blog/adaptyv2
- https://www.adaptyvbio.com/blog/po104/ and https://www.adaptyvbio.com/blog/po103
- https://proteinbase.com/competitions and https://proteinbase.com/competitions/anthropic-adaptyv-2026
- RCSB FASTA for 4KRL, 4KRN, 7OM4, 8UKV, 1YY9, 6B3S, 5SX4, 3GKW, 3C08 (https://www.rcsb.org/fasta/entry/ID); UniProt P01133, P01135, P15514, Q99075, P35070, O14944, Q6UW88
- https://arxiv.org/abs/2506.08293 (DSM, search summary only, UNVERIFIED)
