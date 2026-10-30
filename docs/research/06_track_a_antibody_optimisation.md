# Track A: improving an existing breast cancer scFv in one blind round

Written 2026-10-30. Source tags: (r) I read the full text or raw data; (s) a fetch tool summarised the page for me. UNVERIFIED marks secondary or unopened sources. Each number comes from the source named beside it.

## Bottom line

1. **Framework edits are the safest blind move.** Adaptyv's public EGFR round-2 table holds 7 scFv variants of its cetuximab control with 10 to 14 mutations. The one with 10 framework-only mutations (Cradle) bound at 1.21 nM with high expression. The control gave 6.64 nM in the table (9.94 nM in the paper) with low expression. The six variants that edited CDRs did not improve: four did not express, one bound weakly, one did not bind [S10, S11].
2. **Language-model consensus and antigen-aware inverse folding have the best no-data record.** ESM-1b plus ESM-1v consensus: 14 to 71% of single substitutions improved Kd at least 1.1-fold (Hie). ESM-IF1 with the antigen: about a third of singles improved neutralisation (Shanker). Both built multi-mutants from round-1 measurements, so neither shows blind combination [S3, S4].
3. **Blind multi-mutation designs mostly keep binding when they avoid the paratope.** LAffAb on 6G08: 10 designs with 6 mutations, 8 bound, 2 improved up to 30-fold (preprint). AbLIFT on G6: 18 designs with 4 to 10 mutations, 7 matched or beat the parent. LAffAb on cetuximab edited CDRs and 0 of 12 kept binding [S7, S8].
4. **Trastuzumab has one literature-backed stability module.** VH K30D plus VL S52D (Kabat) gave a scFv with wild-type-like HER2 affinity and no heat aggregation. It supplies 2 of the at least 5 mutations [S12].
5. **Use co-folding scores as vetoes.** I found no benchmark of Boltz-2 or Chai-1 ranking antibody point mutants by affinity. The Boltz-2 affinity head takes small-molecule ligands only [S27, S28].
6. **Spread the 15 over four groups:** stability and expression (4), consensus affinity (5), physics-based (3), bolder framework-weighted (3). Section 5 gives the rules.
7. **A 242-residue chain with a 15-residue linker is exactly trastuzumab** (VH 120, VL 107, no Met). Pertuzumab and cetuximab give 241, sacituzumab 243, pembrolizumab 246, atezolizumab 240 [S31, my arithmetic].
8. **Today:** precompute single-site scans for trastuzumab, pertuzumab and cetuximab. Test the scorers on trastuzumab CDR-H3 data and measure co-folding seed noise on the same run.

## 1. Methods with wet-lab evidence

| Method, study, format | Tested | Improved | Typical gain | Mutations per variant | Combining |
|---|---|---|---|---|---|
| ESM-1b + 5 ESM-1v, Hie 2024, Fab/IgG, 7 antibodies (r) [S3] | 122 variants, at most 20 per antibody, 2 rounds | 14–71% of round-1 singles; 32 of 76 recommended singles | up to 7-fold on 4 mature antibodies, up to 160-fold on 3 unmatured | 1, then 2 or more | Round-2 combinations of round-1 hits mostly improved; Tm rose in 21 of 31 affinity-enhancing variants |
| ESM-IF1 + antigen, Shanker 2024, IgG (r) [S4] | 20 singles on each of 2 antibodies | about 1/3 of singles (IC50) | up to 25-fold neutralisation, 37-fold apparent affinity | 1 to 8 | All 16 multi-mutants beat parent; chosen from round 1 |
| Electrostatics-ranked Rosetta, Lippow 2007, yeast scFv (s) [S5] | cetuximab 5 singles; D44.1 11 singles | cetuximab 3 of 5; D44.1 6 of 9 electrostatics picks, 0 of 2 total-energy picks | cetuximab triple 9.5-fold (490 to 52 pM); D44.1 quadruple 102-fold | 3–4 | Roughly additive; two more cetuximab mutations cut the gain to 3-fold |
| VL-VH interface design AbLIFT, Warszawski 2019 (r) [S7] | G6: 18 designs; h492.1: 20 | G6: 7 matched or beat parent; h492.1: 20 of 20 expressed (parent undetectable) | G6: 270 to 60 pM, Tm +19 °C, IgG yield about 10-fold | 4–10 | Blind |
| Rosetta + PSSM design LAffAb, 2025 preprint (s) [S8] | 6G08: 10 designs; cetuximab 12 | 6G08: 8 bound, 2 improved; cetuximab 0 | 6G08 up to 30-fold | 6 | Blind |
| Zero-shot framework design, Cradle, Adaptyv scFv (r data, s blog) [S9, S11] | 1 in the competition, 11 more in follow-up (blog) | 12 of 12 strong binders per blog, 3 sub-nM | winner 5.5 to 8.2-fold | 10, framework only | Blind |
| AntiBERTy proposals ranked by DeepAB, EGFR round 2 (r) [S10] | 5 | 5 of 5 binders | KD not given | not given | not given |

Clark 2006 improved an anti-VLA1 Fab 10-fold with a quadruple mutant from more than 80 tested variants (abstract) [S6]. I found no prospective wet-lab antibody optimisation with ProteinMPNN, AbMPNN or AntiFold. Retrospectively (s) [S24], the median rank percentile of improved variants in Hie's data was AntiFold 80, ProteinMPNN 73, ESM-IF1 57, AbMPNN 55. On the anti-lysozyme scan Spearman was AntiFold 0.418, ESM-IF1 0.334, AbMPNN 0.322, ProteinMPNN 0.301, ESM-2 0.264. AbBiBench ranks inverse folding best on average across about 15 models (s) [S23]. ProteinMPNN raised expression and Tm of myoglobin and TEV protease [S32]. AbLang and Sapiens did worse than general models in Hie's retrospective test. I found no wet-lab study for AbLang2 or IgLM. Methods that need measured data (DyAb, EVOLVEpro, Li 2023 with 28.7-fold, UNVERIFIED) do not apply [S1, S33].

Epistasis: Lippow's L31E plus L93A gave 4-fold where 8.4-fold was expected. In trastuzumab CDR-H3 yeast scFv data (Mason 2021, reanalysed by Chinery) the binder fraction fell from 98% at one mutation to 17% at ten (s) [S21, S22].

## 2. Which positions

- **Framework and second shell carry about half the gains.** 36 of 76 recommended singles and 18 of 32 improving singles in Hie were framework. Shanker found gains at the interface and far from it. Cradle mutated framework only [S3, S4, S11].
- **Most mutable CDR sites do not touch antigen.** In LAffAb, 15 of 18 mutable D44.1 CDR positions did not contact it. AbLIFT's scan found enhancing mutations at the VL-VH interface core [S7, S8].
- **Distal sites can move affinity and stability together.** VL position 83 (F83A in G6.31) did both by changing VH/VL dynamics [S13, S4]. Trastuzumab VL has F at 83 by my count; untested there.
- **Paratope.** Lippow's three cetuximab hits sit in CDR-L1 and L3 and were picked by the electrostatic term. For trastuzumab, a model trained on Fab data marks H55, H107, H111, H112, H113 as intolerant and H59, H62, H110 as tolerant (numbering probably IMGT, my inference; preprint) [S5, S20].

Expression and stability mutation types:
- Asp or Glu at CDR-adjacent sites (VH 28, 30–33, 35; VL 24, 49–53, 56) gave non-aggregating, well-expressed domains and were retrofitted onto Herceptin [S12].
- Polar replacement of a hydrophobic residue exposed at the former CH1/CL surface raised scFv periplasmic yield 25-fold (VH V84D; numbering not in the abstract) [S14].
- Consensus, hydrophilicity and new hydrogen bonds: 14 designed mutations in scFvD1.3, 75% stabilising, additive, +9.1 kcal/mol combined [S16].
- 4D5 frameworks fold well: grafting CDRs onto them rescued soluble expression [S15].

The scFv format exposes the former CH1/CL surface and leaves the VL-VH interface as the only domain contact. A 15- or 25-residue linker gave a stable monomer for one scFv, though dimers can form by domain swapping and persist [S17]. Hie and Shanker measured Fab or IgG, so gains may not transfer. Keep the linker and orientation, and leave every cysteine alone.

## 3. Known mutations for the likely antibodies

PDB IDs checked in RCSB: 1N8Z (trastuzumab Fab-HER2), 1S78 (pertuzumab), 6OGE (HER2 with both Fabs), 6BGT (trastuzumab Fab mutant), 1YY9 (cetuximab), 6ARU (cetuximab Fab mutant), 9PI9 (sacituzumab-TROP2, 1.56 Å), 5GGS (pembrolizumab-PD-1), 5X8L (atezolizumab-PD-L1) [S31].

| Antibody | Mutation | Source, format, numbering | Effect |
|---|---|---|---|
| Trastuzumab | VH K30D, VL S52D | Dudgeon 2012; scFv; Kabat (s) [S12] | Wild-type-like HER2 affinity by SPR; no aggregation at 80–85 °C. KD values are in a figure I did not see. 6BGT carries exactly these two (my sequence diff) |
| Trastuzumab | one Trp-for-Asp change at the most variable CDR site | Gerstner 2002; Fab, phage (abstract) [S18] | About 3-fold tighter; position not in the abstract, UNVERIFIED |
| Trastuzumab | framework variants of the humanisation | Carter 1992 (abstract) [S19] | humAb4D5-8 binds 250-fold tighter than humAb4D5-1, so framework sites near CDRs are risky |
| Pertuzumab | VL Y96A; VH T30A, A49G, L69F | Patent US 10,000,576, Kabat, Fab (UNVERIFIED, search summary; scanned PDF unreadable) [S34] | Reported affinity gain, stability unchanged |
| Cetuximab | VL S26D, T31E, N93A | Lippow 2007; yeast scFv; chain plus position, residues match 1YY9 (s) [S5] | Singles 1.4, 2.4, 3.5-fold; triple 9.5-fold. VH N56A and T61E did not help |
| Cetuximab | VH K5Q N70S K71R N73T S87T N88D; VL R45K K49R E79D S80P | Cradle winner (r) [S11]; VL numbers are scFv chain positions 179, 183, 213, 214 minus 134 | 1.21 nM against 6.64 or 9.94 nM; high expression |
| Cetuximab | VH S28D, N31D; VL S52D, S56D | Crystal 6ARU, "To be published", Christie and Christ (Christ is Dudgeon's senior author); my sequence diff | Complex crystallised; KD effect unknown to me |
| Sacituzumab | none found | | 9PI9 gives the structure |

I found no primary affinity-mutation list for pertuzumab. Kang 2019 selected histidine variants that bind less at acidic pH on purpose. Pertuzumab shares 4D5-type frameworks with trastuzumab. Its VL has S52 (SASYRYT) and its VH position 30 is T, not K. Papers number by Kabat, IMGT or chain position, so map every literature mutation by sequence alignment.

## 4. Ranking designs in silico

- **Positive retrospective result.** Protenix ipTM correlated at r = −0.61 and Deepconformer at −0.57 with affinity of HyHEL10 heavy-chain variants (higher ipTM, tighter binding), pooled over three antigen systems. For the flexible system r was −0.10; Chai-1 with a crystal template reached −0.59 (r) [S29]. Chai-1 ipTM fell for three of four deleterious minibinder interface mutants. The fourth kept its ipTM and some binding (r) [S30].
- **Negative result.** In a benchmark of Boltz-2, Chai-1, IntFold and AF3 on nanobodies, AF3 ipSAE rose as CDR3 residues were replaced by alanine and was the worst measure; AF3 pLDDT fell as it should. The paper has no scFv data, and its figure caption and alt text disagree on the ipSAE trend. Best median AUC for true binders against baits was 0.90 (IntFold Rosetta dG); Boltz-2's lead was not significant (r) [S26].
- **Competition data.** On 400 EGFR round-2 designs, ipTM gave AUC 0.64 for binding and 0.58 for expression; ESM2 PLL gave 0.72 and 0.44 (r) [S10]. Unnormalised PLL favoured short sequences.
- **Boltz-2 affinity head.** The docs allow one small molecule of at most 128 atoms and state that Boltz "only supports the computation of affinity of small molecules to protein targets" (s). It does not apply to scFv-antigen binding. Chai-1's README lists no affinity output (s) [S27, S28].
- **ipSAE thresholds** in Adaptyv's skill target true-versus-false complexes. The skill gives no mutant-level evidence [S2].
- **ΔΔG predictors.** On AB-Bind, correlations were r = 0.16–0.45. FoldX reached AUC 0.88 for improved against weakened binders and put 42% of the 80 most improved binders in its top 5% (abstract) [S25]. A 2025 study concluded that machine-learned antibody ΔΔG needs orders of magnitude more data (search summary, UNVERIFIED) [S35]. In Lippow the electrostatic term beat total energy.
- **False positives to expect.** Language models pull toward natural sequence, and four of the seven Adaptyv cetuximab variants did not express. ESM-IF1 misses gains that need a backbone shift [S4]. Co-folding confidence can stay high or rise while binding drops. Flexible paratopes break rank order.

## 5. One-day plan

| Time | Step and tools | Output |
|---|---|---|
| 0:00–0:45 | Align the scFv to the six reference VH/VL pairs here. Get orientation, linker, Met. Number with ANARCI or abnumber. Label each position: CDR, paratope (within 5 Å of antigen), second shell (5–8 Å), VL-VH interface, framework | 242-row position table |
| 0:45–1:15 | Boltz-2 and Chai-1, 5 seeds, with and without a template; compare the pose with the crystal | Pose check; seed-noise floor |
| 1:15–3:00 | Score about 4,300 substitutions at non-Cys, non-linker sites: ESM-1b, 5 ESM-1v, ESM-2; ESM-IF1 and AntiFold with antigen; ProteinMPNN scoring. Keep a substitution if LLR is above 0 in at least 2 of 6 sequence models and it is the top residue in one (Hie's rule), and a structure model agrees | Ranked table |
| 3:00–4:00 | FoldX or Rosetta ΔΔG (fold and binding) on the top 300, electrostatic term first. Drop new unpaired Cys, new NG/DG motifs, new hydrophobic patches | About 60 tagged mutations |
| 4:00–5:30 | Build the 15 designs. Co-folding veto: reject if ipTM or ipSAE falls below the wild-type mean minus 2 SD of seed noise, or the interface moves more than 1.5 Å (my thresholds) | 15 sequences |
| 5:30–8:00 | Checks (at least 5 mutations, pairwise distance of 3 or more, pI shift under 1), submission file, slides, buffer | Submission and talk |

Portfolio for trastuzumab (n = 15):

| Group | n | Mutations | Where | Source | Risk |
|---|---|---|---|---|---|
| A stability, expression | 4 | 5–7 | Framework, VL-VH interface, CDR-adjacent D/E sites, at least 8 Å from antigen | K30D + S52D in 2 designs; consensus and interface repacking in the rest | Low for binding |
| B consensus affinity | 5 | 5–6 | Framework, second shell, at most 1 CDR change, no H3 hot spot | 2 of 3 model families agree | Medium |
| C physics | 3 | 5–6 | Contact shell outside hot spots, plus 2–3 from group A | Lippow-type polar-to-hydrophobic or added charge | High |
| D bold | 3 | 8–10 | Union of A and B, framework-weighted | Cradle pattern | Medium to high |

Rules: each mutation appears in 2 to 4 designs so the talk can attribute outcomes, no CDR-H3 mutation in more than 2 designs, no two mutations within 4 Å in one design unless co-scored. My estimate, not a published rate: 9 to 13 bind, 1 to 3 beat the parent, several express better. The BLI series is 30 to 1000 nM (team doc 04), so KD below 1 nM may not resolve; koff and expression then carry the outcome (my inference).

If the antibody is not trastuzumab:
- **Pertuzumab:** same groups; S52D is the only literature carry-over.
- **Cetuximab:** reuse the Cradle pattern (public in S11), Lippow's triple as group C, no CDR combinations. The 2024 control expressed poorly, so group A grows to 5. If the lead equals that control, the Cradle mutations are public.
- **Sacituzumab:** no literature. TROP2 is dimeric in 9PI9, so avidity can flatter KD. Rely on groups B and D.
- **Pembrolizumab, atezolizumab:** 5GGS and 5X8L exist. Neither matches 242 residues.

Prepare today:
1. Install and test ESM-1v (5 models), ESM-2, ESM-IF1, AntiFold, ProteinMPNN scoring, Boltz-2, Chai-1. The team pipeline has no ESM or AntiFold step.
2. Run the full scan on the 1N8Z, 1S78 and 1YY9 antibody chains.
3. Score the 419 trastuzumab heavy-chain variants in AbBiBench (dataset location UNVERIFIED) to set scorer weights [S23].
4. Write the seed-noise and pose-check script.
5. Ask organisers about orientation, Met, linker mutations and novelty rules.

## What could not be verified

- Cradle's follow-up 12 sequences and KDs. The two Cradle blog pages disagree on mutation counts; the CSV confirms 10 for the winner. Control KD: 6.64 nM (CSV) against 9.94 nM (paper).
- Position of Gerstner's Trp-for-Asp change; Dudgeon's Herceptin KD values; the pertuzumab patent mutations; Li 2023 and Hummer 2025 numbers; AbBiBench per-model Spearman values (figures only).
- No benchmark of Boltz-2 or Chai-1 ranking antibody point mutants was found, and no wet-lab evidence for AbLang2 or IgLM.
- LAffAb is an unreviewed preprint; I read about 78% of it through a summarising tool, as with Lippow, Dudgeon (Herceptin part), AbBiBench, Chinery, AntiFold, the Boltz docs, the Chai README and both Cradle blogs.
- Whether the antibody is trastuzumab. Every trastuzumab module is conditional on it.

## Sources

S1 https://www.adaptyvbio.com/blog/po101 (s)
S2 https://github.com/adaptyvbio/protein-design-skills (skills esm, ipsae, germinal read raw) (r)
S3 Hie 2024 https://pmc.ncbi.nlm.nih.gov/articles/PMC10869273/ (r)
S4 Shanker 2024 https://pubmed.ncbi.nlm.nih.gov/38963838/ ; PDF https://static1.squarespace.com/static/562d385fe4b0e2b30b0c79b4/t/668c6729471412128dd04ce0/1720477482498/science.adk8946.pdf (r)
S5 Lippow 2007 https://pmc.ncbi.nlm.nih.gov/articles/PMC2803018 (s)
S6 Clark 2006 https://pmc.ncbi.nlm.nih.gov/articles/PMC2242497 (abstract)
S7 Warszawski 2019 https://pmc.ncbi.nlm.nih.gov/articles/PMC6728052/ (r)
S8 LAffAb https://www.biorxiv.org/content/10.1101/2025.11.26.690765v1.full (s)
S9 Cradle https://cradle.bio/blog/adaptyv2 ; https://cradle.bio/blog/adaptyv-protein-design-competition (s)
S10 EGFR competition paper https://www.biorxiv.org/content/10.1101/2025.04.17.648362 (r)
S11 https://github.com/adaptyvbio/egfr_competition_2 (results/result_summary.csv) (r)
S12 Dudgeon 2012 https://pmc.ncbi.nlm.nih.gov/articles/PMC3390889/ (s)
S13 Koenig 2017 https://pubmed.ncbi.nlm.nih.gov/28057863/ (abstract)
S14 Nieba 1997 https://pubmed.ncbi.nlm.nih.gov/9194169/ (abstract)
S15 Jung and Plückthun 1997 https://pubmed.ncbi.nlm.nih.gov/9415446/ (abstract)
S16 Monsellier and Bedouelle 2006 https://pubmed.ncbi.nlm.nih.gov/16926023/ (abstract)
S17 Arndt 1998 https://pubmed.ncbi.nlm.nih.gov/9737871/ (abstract)
S18 Gerstner 2002 https://pubmed.ncbi.nlm.nih.gov/12206766/ (abstract)
S19 Carter 1992 https://pmc.ncbi.nlm.nih.gov/articles/PMC49066/ (abstract)
S20 Bachas 2022 https://www.biorxiv.org/content/10.1101/2022.08.16.504181v1 (r)
S21 Mason 2021 https://pubmed.ncbi.nlm.nih.gov/33859386/ (abstract)
S22 Chinery 2024 https://www.biorxiv.org/content/10.1101/2024.03.26.586756v1.full (s)
S23 AbBiBench https://arxiv.org/abs/2506.04235 (s)
S24 AntiFold https://pmc.ncbi.nlm.nih.gov/articles/PMC11961221/ (s)
S25 AB-Bind https://pmc.ncbi.nlm.nih.gov/articles/PMC4815335/ (abstract)
S26 Nanobody benchmark https://www.frontiersin.org/journals/bioinformatics/articles/10.3389/fbinf.2026.1858891/full (r)
S27 Boltz docs https://github.com/jwohlwend/boltz/blob/main/docs/prediction.md (s)
S28 Chai-1 https://github.com/chaidiscovery/chai-lab (s)
S29 Xu 2026 https://pmc.ncbi.nlm.nih.gov/articles/PMC13679318/ (r)
S30 Broske 2026 https://pmc.ncbi.nlm.nih.gov/articles/PMC13494018/ (r)
S31 RCSB entries and sequences, https://data.rcsb.org/rest/v1/core/entry/1N8Z and the same path for each ID (r)
S32 Sumida 2024 https://pmc.ncbi.nlm.nih.gov/articles/PMC10811672/ (abstract)
S33 DyAb https://pmc.ncbi.nlm.nih.gov/articles/PMC13531962/ (abstract)
S34 Patent https://image-ppubs.uspto.gov/dirsearch-public/print/downloadPdf/10000576 (UNVERIFIED)
S35 Hummer 2025 https://pmc.ncbi.nlm.nih.gov/articles/PMC12374841/ (search summary, UNVERIFIED)
