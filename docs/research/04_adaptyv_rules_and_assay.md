# 04. Adaptyv submission rules and assay

Research date 2026-10-24. Source: subagent report, merged and edited. "UNVERIFIED" marks anything not confirmed in a primary source. Source codes S1 to S42 are listed in section 7.

## Bottom line

1. No public page states the rules of the 31 Oct challenge. The Luma page has one sentence about it (re-read 2026-10-24), and HerCode, Hackathon Hub and Adaptyv pages add nothing [S1, S2, S3].
2. The nearest precedent for submission rules is the Anthropic x Adaptyv competition running now (28 Sep to 8 Nov 2026): single chain, 10 to 250 residues per chain, no tags in the submitted sequence, a ranked CSV with name, sequence and molecule_class, at most 20 designs per participant, de novo only, free testing, results public on Proteinbase [S14, S15].
3. The nearest precedent for a one-day format is the muni/BioArena x Adaptyv TREM2 hackathon: at most 10 designs per team, ranked by Boltz-2 ipSAE, top 100 of 141 sent to the lab, 89 expressed and 37 bound [S17, S18].
4. Adaptyv expresses designs cell-free (E. coli by default) and appends a 64-residue C-terminal tag (split-GFP11 plus Twin-Strep). It captures the binder on a Strep-Tactin sensor and flows the target over it, by BLI (Gator Bio) or SPR (Carterra) [S24 to S27].
5. A design is a binder if it gives a fittable curve with a KD, or an association shift of at least 300% over the negative control. Screening tiers are Strong below 50 nM, Medium 50 nM to 1 uM and Weak above 1 uM [S15, S28, S29].
6. Expression success was 73% in EGFR round 1, 95% in round 2, 89% in TREM2 and 95% for the Claude designs. In the public round-2 data none of the 210 designs under 100 residues failed to express (the research agent's count) [S6, S11, S18, S20].
7. Of the three claims in the notes, the 10-per-team ipSAE format and the 700/200/5 EGFR figures hold. "1 to 5 designs per designer" appears nowhere we could read.
8. Ask the organisers at kickoff for the length limit, designs per team, ranking metric, how tags and the start Met are handled, and who pays.

## 1. Submission rules per competition or round

| Event | Length | Tags | Expression | Designs per team | Format | Ranking | Tested | Src |
|---|---|---|---|---|---|---|---|---|
| EGFR round 1 (Jul to Sep 2024) | under 200 aa | not stated | cell-free (host not named) | unlimited; leaderboard showed each person's top 10 | web form, auto-scored | AF2-M iPAE | 201 of 726 (100 top + 101 hand-picked) | S6, S7 |
| EGFR round 2 (Oct to Nov 2024) | under 250 aa | not stated | cell-free | unlimited | web form | Borda count of unnormalised ESM2 PLL, iPAE, ipTM | 400 of 1,131 (100 top + 300 hand-picked, at least one per designer) | S6, S8, S10 |
| GEM x Adaptyv RBX1 (15 Feb to 26 Mar 2026) | 250 aa or fewer | tags and signal peptides not allowed | not stated; BLI | up to 100, ranked | ranked CSV plus 2-page method PDF | own ranking; slots = 300 / number of teams; de novo, at least 25% edit distance to UniRef50 | 322 listed (rules page says 300) | S16, S13 |
| muni/BioArena x Adaptyv TREM2 hackathon (28 Feb 2026) | not stated | not stated | cell-free (lab standard, not stated on page); BLI | up to 10 | ranked by ipSAE | Boltz-2 ipSAE_d0res (PAE and distance cutoff 15); at least 10 aa from known therapeutics | top 100 of 141 | S17, S18 |
| Anthropic x Adaptyv (28 Sep to 8 Nov 2026, 5 challenges) | 10 to 250 aa per chain | none in sequence; own C-terminal Twin-Strep captures the design, tags are off limits as epitope | cell-free | up to 20 (Tracks 2, 3); 20 to 40 (Track 1) | CSV: name, sequence, molecule_class; Fab as VH:VL; row order = rank | Tracks 2, 3: Claude-based selection, no single metric; Track 1: top 20 passing filters; de novo, zero-shot, novelty level 3 or higher; winner set by lab result | about 1,500 per challenge | S14, S15 |
| Foundry service (docs) | 8 to 700 aa (expression page says about 50 to 700) | none to add; Adaptyv appends tag at C-terminus (N-terminus on request) | cell-free E. coli default, yeast alternative, CHO not available | no cap; minimum 20 (expression) or 24 (binding) per target | CSV or FASTA in portal, or API | you choose | your list; duplicates rejected | S24, S25, S30, S37 |

Cost, deadlines, IP and tools:

- Competitions test at no cost to participants. Anthropic and Adaptyv fund the 2026 competition, Twist part-funded DNA in EGFR round 2, and Track 2 gets $500 of Modal credit. EGFR winners got 10 free assays. There is no cash prize in the 2026 competition. RBX1 paid $1,000 and $100 [S6, S14, S16].
- Foundry: binding screening from $99 per protein (includes expression). The configurator showed $169 per sequence for 60 proteins in duplicate, with assay depth not visible in the text. Target reagent costs extra. Releasing results on Proteinbase gives 20% off. Turnaround is 3 to 4 weeks [S39, S31].
- 2026 windows: each challenge is open two weeks and closes 23:59 AoE. Challenge 5 opens Oct 26 and closes Nov 8. Lab results are due early January [S14].
- IP: in competitions you keep ownership, but sequences and results are published under ODC-BY, and the terms advise filing patents before submitting (the RBX1 page says ODC-ODbL). In the Foundry service your data is confidential and yours, with a licence for Adaptyv's internal models [S14, S16, S31].
- Tools and AI: no tool restrictions beyond licences. Commercial tools need your own licence. Embedded instructions in a submission can disqualify it. Every design is biosecurity-screened before DNA synthesis [S14, S15].

Claims in docs/challenge_notes.md:

- Up to 10 per team, ranked by ipSAE, 100 to the wet lab: verified for the TREM2 hackathon [S17, S18]. Nipah used average ipSAE per collection of designs [S23].
- EGFR round 1, 700+ submitted, 200 validated, 5 bound: verified (726, 201, 5 binders = 2.5%, plus 2 later disqualified) [S6]. The Substack says 90 designers and the paper says about 50 [S12, S6].
- "1 to 5 designs per designer": UNVERIFIED. The cited write-up does not say it [S40]. What the primary sources support is "at least one protein from each designer" in the round-2 hand-picked set [S8].

## 2. Assay and hit definitions

- Expression: cell-free, PURE-type E. coli system with DsbC and glutathione for disulfides. The DNA is a codon-optimised linear template, and expression runs 8 h at 37 C. Expression is read by split-GFP against a standard curve and reported as expressed, low or not detected [S32, S34, S36, S26].
- Capture: Twin-Strep on Strep-Tactin XT. The design is the ligand and the target is the analyte, so avidity arises only from the target (dimeric or Fc-tagged targets look tighter). Targets must not carry Twin-Strep or biotin [S25, S35].
- BLI (Gator Bio Pro): 25 C, 5 Hz, half-log series of 4 concentrations from 1000 to 30 nM, 220 s association, 240 s dissociation, duplicates. SPR (Carterra LSA XT): single-cycle, 7 concentrations in the docs example, 5 concentrations (10 nM to 1 uM) with 5 min association in the Claude benchmark [S26, S27, S20].
- Fitting: global 1:1 Langmuir, with bivalent-analyte and two-state fits tried. An apparent KD is reported where those fit better. The quantifiable range is about 0.1 nM to 10 uM [S26, S28].
- Hit call: True if a fittable curve with a KD exists, or if the association shift is at least 300% over the negative control. False if expressed but at or below control. Unknown if no or low expression or a marginal signal. Screening tiers: Strong below 50 nM, Medium (called Moderate in one guide) 50 nM to 1 uM, Weak above 1 uM [S15, S26, S28, S29].
- No KD cut-off defines "binder" in the docs. In EGFR round 1 a 22.9 uM design counted as a binder [S6]. Hit definitions elsewhere differ (RFdiffusion: 50% of control signal; BindCraft: SPR signal at 10 uM) [S6].
- Failed expression: the binding call is Unknown (no or low expression) and the expression call is "not detected". The TREM2 post lists 11 not expressed and 52 not binding. The Claude report counts an unexpressed design as a tested non-binder [S18, S20, S36].
- Rounds 1 and 2 used BLI with at least 2 replicates, with EGF and cetuximab controls. The 2024 paper does not give the instrument, the concentration series or the three "binding thresholds" of its Figure 8 [S6].

## 3. Developability and expression findings

Expression rates: EGFR R1 146/201 (73%); R2 378/400 (95%) [S6]. TREM2 89/100 [S18]. Claude benchmark: 61 of 1,296 not expressed (95% expressed) [S20]. Nipah and RBX1: expression rate not found (UNVERIFIED); hit rates 9.3% and 2.8% [S13]. Hit rates for context: R1 2.5%, R2 13%, TREM2 37%, Claude 26.8% [S6, S18, S19].

What the EGFR paper reports for expression (high expression versus everything else, binder-only models, ROC AUC) [S6]:

- Glu fraction 0.77, Lys fraction 0.73, Rosetta solvation per residue 0.68, backbone H-bonds per residue 0.68, helix fraction 0.67, ipTM 0.58, pLDDT 0.55, ESM PLL 0.44 (not predictive).
- A DNA-feature classifier reached F1 0.93 and 90% accuracy. Sequence length was its top feature (7.13% accuracy gain), followed by non-unique 6-mers (repeats).
- Round-1 analysis: lower Aggrescan3D, higher NetSolP usability and smaller size went with expression.
- For binding, ESM PLL had the best AUC (0.72), then pLDDT 0.66 and ipTM 0.64. Unnormalised PLL favours short sequences. The paper reports no pI, cysteine or instability-index results.

Other statements:

- Round-2 expression rose from 73% to 95% after SolubleMPNN, NetSolP and SoluProt came into use, but the paper says causality is not established [S6].
- All ten designs from the BindCraft authors expressed highly. Diversified binders expressed at 92%. Five GB1-fused composites all expressed and none bound [S6, S8]. The top 100 by metric gave 24% binders, the hand-picked set 10% [S8].
- In TREM2, 19% of agent designs and 7% of human designs failed to express. Agent designs carried median alanine fractions of 0.27 to 0.46 (corpus 23.3%), but the post draws no link to failure [S18].
- The Claude designs (50 to 120 residues, 1,133 of 1,320 sequenced with SolubleMPNN) expressed at 95% [S20].
- Adaptyv's docs list high hydrophobicity, repeats, unusual composition and size above 50 kDa as expression risks, and advise removing unpaired cysteines [S32, S33].

The research agent's own check on the public round-2 table (402 rows: 281 high, 97 medium, 2 low, 22 none; scripts in the session scratchpad) [S11]. It reproduces the paper's Glu AUC (0.77) and Lys AUC (0.72), which serves as a sanity check.

- All 22 non-expressers are 119 to 127 aa (9, VHH-like) or 239 to 250 aa (13, scFv-like). None of 210 designs under 100 aa failed. 22 of 192 designs at or above 100 aa did (11.5%).
- Within the 100+ aa group, lower computed pI (median 5.0 for expressed vs 7.9 for none) and higher Glu and Lys fractions went with expression (AUC 0.24 for pI, 0.82 Glu, 0.79 Lys). Format and charge are confounded and there are only 22 failures.
- Unnormalised ESM PLL separates expressed from none overall (AUC 0.74), but only through length. Within the 100+ aa group it is 0.48 [0.36 to 0.59]. Odd cysteine counts did not differ (3 of 22 vs 47 of 380).
- Adaptyv's round-1 blog said ESM2 likelihood tracks expression. The paper's AUC of 0.44 does not support that [S7, S6].

## 4. Checklist of sequence-level constraints before submitting

1. Alphabet: 20 natural amino acids, uppercase, no X, B, Z, U or O. Chains are joined by ":" [S37].
2. Single chain, 10 to 250 aa, unless the organisers say otherwise. With no stated limit, stay at or below 200 (round 1). Prefer 100 aa or fewer, since the public EGFR round-2 table has no expression failure below 100 aa and the first failures at 119 aa. 50 to 120 aa is the range of the Anthropic campaign, which expressed 95% [S6, S15, S20, S11].
3. No His, FLAG, Strep, GFP11, signal peptide or extra linker. Adaptyv appends `GGGSGGGSRDHMVLHEYVNAAGITGGGSGGGSGGSAWSHPQFEKGGGSGGGSGGSAWSHPQFEK`. Keep the binder C-terminus free of the interface and never design against the tag [S15, S25, S16].
4. Even cysteine count and no free cysteine. Paired disulfides are supported [S32, S33].
5. Run NetSolP and SoluProt. Design sequences with SolubleMPNN or similar [S33, S20].
6. No hydrophobic run of 6 or more, low mean Kyte-Doolittle hydropathy (a screen in the 2026 selection prompt), low GRAVY [S14, S38].
7. No low-complexity or repeated segments (non-unique 6-mers) and no homopolymer runs. Alanine fraction under 0.3 (Adaptyv's own BoltzGen filter) [S6, S32, S38].
8. Prefer acidic, Glu- and Lys-rich surfaces over basic or near-neutral pI. Treat this as a weak prior [S6, S11].
9. Normalise ESM PLL by length, or skip it for expression [S6, S9].
10. Novelty: at least 10 residues from known binders (EGFR, TREM2), 25% edit distance (RBX1), novelty level 3 or higher (2026). HER2 has many public binders, and binder-prescreen ships no HER2 set (only TNF-a, GLP1R, TfR1), so build a reference set from PLAbDab, Thera-SAbDab and patents [S21, S22].
11. Unique sequences. The API rejects duplicates [S37]. Order the CSV by your own rank.
12. Run ipSAE with several seeds. Boltz-2 scores varied with hardware and seed in the Nipah round [S42].
13. If controls are allowed, add a known binder as scFv, Fab or small scaffold (no IgG). Controls count towards the total, and plates fit 90 or fewer, or 120 or more [S24, S30].
14. Write a short method note, and assume any tested design and its result will be published openly [S14].

## 5. What could not be verified

- Any rule of the 31 Oct challenge: length, designs per team, metric, tags, who pays, deadline, whether wet-lab testing happens at all. HerCode's May hackathon (FGFR2 track) shows no Adaptyv role [S4, S41]. The Catalyst hackathon page says only that winners get validation [S5].
- Whether Anthropic x Adaptyv challenge 5 (opens Oct 26, target unannounced) is related to the hackathon. This is inference only.
- Tags, instrument, concentrations and exact thresholds for the 2024 EGFR rounds. Whether the 64-residue tag counts toward length limits. Start-Met handling (the template example nanobody has no Met).
- Expression rates for Nipah and RBX1. Nipah rules beyond aggregate numbers (the research on that pathogen target was stopped on purpose).
- Source conflicts: docs length (8 vs 50 minimum), minimum order (20 vs 24), Medium vs Moderate; TREM2 dates (Proteinbase 28 Feb vs blog 13 May); R1 designer count (90 vs about 50) and expression (73% vs 76% in a blog); R2 binders (53 in the paper, 55 flagged true in the CSV); the Proteinbase R2 page text looks copied from round 1; Anthropic Track 1 tested count (about 15 in terms, 20 on pages).
- Whether HER2 is in Adaptyv's target catalogue (page not readable by this agent; see file 02). X and LinkedIn posts were blocked. The Overath et al. meta-analysis, the Twist page and the Hugging Face Claude dataset were not opened for this file.

## 6. Reading for kickoff

Questions that follow from this file are collected in SUMMARY.md. The ones about tags and the start Met matter because the appended tag is 64 residues long. If the organisers apply a length limit to the full construct, a 120-residue design becomes 184 residues.

## 7. Sources

S1 https://luma.com/4ogriih0 : event page, one mention of the Adaptyv challenge; primary, no rules.
S2 https://www.hercode.ch/events : lists event, no details; primary.
S3 https://hackathonhub.eu/hackathons/healthtech : aggregator listing; secondary.
S4 https://bsse.ethz.ch/news-and-events/d-bsse-news/2026/06/24-hours-of-protein-design.html : student account of the HerCode May hackathon; secondary.
S5 https://luma.com/4gca12wa : Catalyst hackathon, Adaptyv track; primary, no rules.
S6 https://www.biorxiv.org/content/10.1101/2025.04.17.648362v2.full.pdf : EGFR paper, full text; primary.
S7 https://www.adaptyvbio.com/blog/po102/ : round 1 results and expression; primary.
S8 https://www.adaptyvbio.com/blog/po104/ : round 2 selection, expression, AUCs; primary.
S9 https://www.adaptyvbio.com/blog/po103/ : length bias in round 2; primary.
S10 https://github.com/adaptyvbio/competition_metrics : round 2 metric; primary.
S11 https://github.com/adaptyvbio/egfr_competition_2 : round 2 data (results/result_summary.csv); primary.
S12 https://adaptyvbio.substack.com/p/designing-better-cancer-drugs-using : round 1 headline numbers; primary.
S13 https://proteinbase.com/competitions : tested counts and hit rates; primary.
S14 https://proteinbase.com/competitions/anthropic-adaptyv-2026 (also /terms, /selection) : 2026 rules, terms, selection prompt; primary.
S15 https://proteinbase.com/competitions/anthropic-adaptyv-2026/challenges/tfr1 (egfr, tnf-alpha, glp1r same) : requirements, FAQ, hit rule; primary.
S16 https://proteinbase.com/competitions/gem-adaptyv-rbx1 and https://www.gembio.ai/gemadaptyv2026 : RBX1 rules; primary.
S17 https://proteinbase.com/competitions/bioarena-adaptyv-trem2 : hackathon rules; primary.
S18 https://www.adaptyvbio.com/blog/agents-vs-humans : TREM2 results; primary.
S19 https://www.adaptyvbio.com/blog/anthropic-1 : Claude benchmark, 26.8%; primary.
S20 https://www-cdn.anthropic.com/30bf50e22a01388bb29bf077ee3f244531594b7a.pdf : Claude report, assay and expression; primary.
S21 https://www.adaptyvbio.com/blog/novelty : novelty scale; primary.
S22 https://github.com/adaptyvbio/binder-prescreen : prior-art screen and CSV template; primary.
S23 https://github.com/adaptyvbio/nipah_ipsae_pipeline : average ipSAE per collection (aggregate rule only); primary.
S24 https://docs.adaptyvbio.com/docs/experiment-setup/supported-protein-formats.md : length, formats, systems; primary.
S25 https://docs.adaptyvbio.com/wiki/technologies/construct-tag-design.md : appended tag; primary.
S26 https://docs.adaptyvbio.com/wiki/methods/bli-methods.md : BLI method; primary.
S27 https://docs.adaptyvbio.com/wiki/methods/spr-methods.md : SPR method; primary.
S28 https://docs.adaptyvbio.com/docs/experiment-types/binding.md : tiers, labels; primary.
S29 https://docs.adaptyvbio.com/guides/binding-characterization.md : calls, fit quality; primary.
S30 https://docs.adaptyvbio.com/docs/experiment-setup/creating-experiments.md : upload, minimums, controls; primary.
S31 https://docs.adaptyvbio.com/docs/faq/overview.md : cost, timing, data ownership; primary.
S32 https://docs.adaptyvbio.com/wiki/technologies/cell-free-protein-expression.md : system and failure causes; primary.
S33 https://docs.adaptyvbio.com/guides/improve-protein-expression.md : advice; primary.
S34 https://docs.adaptyvbio.com/wiki/technologies/split-GFP.md : expression readout; primary.
S35 https://docs.adaptyvbio.com/docs/experiment-setup/target-proteins.md : target rules; primary.
S36 https://docs.adaptyvbio.com/docs/experiment-types/expression.md : expression calls; primary.
S37 https://devs.adaptyvbio.com/api/v1/openapi.json : alphabet, duplicates; primary.
S38 https://github.com/adaptyvbio/protein-design-skills : Adaptyv-written QC heuristics (skills/protein-qc), not outcome data; primary.
S39 https://www.adaptyvbio.com/services/binding : prices; primary.
S40 https://blog.booleanbiotech.com/adaptyv-binder-design-competition : checked, no per-designer rule; secondary.
S41 https://futureofhealthhack.hercode.ch/ : HerCode May hackathon site, tracks and partners, no Adaptyv mention; primary.
S42 https://www.adaptyvbio.com/blog/nipah-submissions : ipSAE reproducibility (aggregate only); primary.
