# 05. Past winners and jury criteria

Research date 2026-10-24. Source: subagent report, merged and edited. The merging session re-read the Luma event page and the D-BSSE write-up of the HerCode May hackathon. "UNVERIFIED" marks anything from a secondary source, a social-media post that could not be opened, or a team's own claim without an organiser page.

## Bottom line

1. In every Adaptyv event we could read, the in-silico metric decided which designs reached the lab, and wet-lab KD decided who won. The metric itself predicted binding poorly: EGFR round 2 gave ROC-AUC near 0.5 for ipTM, iPAE and unnormalised ESM2 likelihood, and TREM2 gave ipSAE AUC 0.65 (the figure from the Adaptyv/muni write-ups; file 01 has recomputed values of 0.74 on a 110-design subset).
2. Hit rates varied from about 1% to 38% with target and execution. Plan for a low rate.
3. The hackathon template Adaptyv now uses is up to 10 designs per team, ranked by Boltz-2 ipSAE, top 100 to the lab (TREM2, 28 Feb 2026).
4. Winners used different tools (Mosaic, BindCraft, PXDesign, optimised antibodies), so no single pipeline wins.
5. Where organisers looked beyond the metric, they named method novelty and design diversity (Adaptyv), or workflow novelty and spinout potential (Berlin).
6. Two claims in the repo notes needed correction (done in `docs/challenge_notes.md` on 2026-10-26). "1 to 5 designs per designer" is in no source we could open. The May 2026 HerCode protein track was a Lemna Bio track with no public results, no Adaptyv involvement and no reported wet-lab step.
7. The recurring failures are time lost to environment setup, expression failures, similarity rules and near-duplicate submissions.

## 1. Events

### 1.1 Adaptyv EGFR round 1 (submissions to 18 Aug 2024, results Sep 2024)

- Format: online, unlimited submissions, about one month. Rules: single chain, under 200 aa, at least 10 aa from known binders. Only a participant's top 10 appeared on the leaderboard. [bioRxiv v1](https://www.biorxiv.org/content/10.1101/2025.04.17.648362v1.full), [Adaptyv Substack](https://adaptyvbio.substack.com/p/crowdsourcing-better-cancer-drugs)
- Scoring and selection: 726 submissions. 100 were ranked by AF2 iPAE and 101 hand-picked "to increase diversity" across methods and to "test novel approaches". [bioRxiv v1](https://www.biorxiv.org/content/10.1101/2025.04.17.648362v1.full)
- Outcome: 73% expressed (146/201). 5 binders (2.5%). Two further binders were disqualified for similarity to existing therapeutics. [Adaptyv blog](https://www.adaptyvbio.com/blog/po102/)
- Winners by KD:
  - First: Martin Pacesa and Lennart Nickel, BindCraft (AF2 hallucination), KD 4.91e-7 M.
  - Second: Khondamir Rustamov, AF2Seq plus SolubleMPNN, iPAE rank 54.
  - Third: Adrian Tripp and Sigrid Kaltenbrunner, ProtFlow (RFdiffusion, LigandMPNN, ESMFold and ColabFold filters, 3 cycles), iPAE rank 89. [Adaptyv blog](https://www.adaptyvbio.com/blog/po102/)
- Leaderboard topper versus lab winner: Alex Naka held the top iPAE spots with model-based active learning, yet is listed as a special mention and not as a KD winner. [Boolean Biotech](https://blog.booleanbiotech.com/adaptyv-binder-design-competition), [Adaptyv blog](https://www.adaptyvbio.com/blog/po102/)
- Jury: none. Prizes went by KD (10 free assays for the top 3). [bioRxiv v1](https://www.biorxiv.org/content/10.1101/2025.04.17.648362v1.full)

### 1.2 Adaptyv EGFR round 2 (4 Oct to 4 Nov 2024, results post 8 Dec 2024)

- Format: in-person kickoffs in five cities (Lausanne among them) on 12 Oct, then online submission. Limit: under 250 aa. [Adaptyv Substack](https://adaptyvbio.substack.com/p/designing-better-cancer-drugs-using), [Proteinbase](https://proteinbase.com/competitions/adaptyv-egfr-binder2)
- Selection: 1,131 submissions from about 130 designers. The top 100 came from a Borda count of ESM2 likelihood, iPAE and ipTM, plus 300 hand-picked designs, and every designer had at least one design tested. The Proteinbase page describes the picks as 100 by pAE and 100 by "rule of cool" and lists 200 tested, which conflicts with the 400 reported elsewhere. [bioRxiv v1](https://www.biorxiv.org/content/10.1101/2025.04.17.648362v1.full), [Proteinbase](https://proteinbase.com/competitions/adaptyv-egfr-binder2)
- Outcome: 378 of 400 expressed (95%). 53 bound (13 to 14%). 30 of 130 participants had a binder. Hit rate by strategy: diversified known binders 28%, optimised binders 18%, BindCraft de novo 9% (6 of 65; Adaptyv's blog categories, which differ from the designer-declared categories in file 01, section 2, so do not compare the hit rates directly). Top-100 picks bound 24% of the time, Adaptyv's diversity picks 10%. [Adaptyv blog](https://www.adaptyvbio.com/blog/po104/)
- Winners:
  - Cradle: cetuximab framework variant, 1.21 nM.
  - Chris Xu: nanobody by CDR grafting.
  - Aurelia Bustos: TGFa-based, RFdiffusion plus ProteinMPNN, 8 of 10 bound.
  - Best de novo: Nickel and Pacesa, BindCraft, 91.5 nM per the Adaptyv post (the bioRxiv paper gives 82 nM; see file 01). [Adaptyv blog](https://www.adaptyvbio.com/blog/po104/)
- Pitfalls named by the organisers: unnormalised ESM2 likelihood favoured very short sequences, and Cradle's entry first ranked 300 to 600 for that reason. Some diversified designs drifted too far from natural sequences to express. [Adaptyv blog](https://www.adaptyvbio.com/blog/po104/), [Cradle](https://www.cradle.bio/blog/adaptyv-protein-design-competition)
- Cradle spent about 30 minutes on setup and ran overnight. It submitted 10 designs, and only 1 was tested. [Cradle](https://www.cradle.bio/blog/adaptyv-protein-design-competition)

### 1.3 Bits to Binders (UT Austin, BioML Society; kickoff 24 Aug 2024, deadline 30 Sep 2024)

- Format: 5 weeks, 28 teams, 80-aa binders to CD20 as CAR-T binding domains, about 500 designs asked per team. Sequences were released to the public domain. Adaptyv is credited for binding affinity. [BioML Society](https://www.biomlsociety.org/challenge), [bioRxiv](https://www.biorxiv.org/content/10.64898/2026.03.03.709355v2.full)
- Winners: Perez Lab Gators, Amigo Acids and Schoeder Lab. The highest round-1 hit rate was Nucleate UK (38%). Top teams all used RFdiffusion or Chroma plus ProteinMPNN or SolubleMPNN. Perez Lab Gators was the only team to fold designs inside the full CAR. [BioML Society](https://www.biomlsociety.org/challenge), [bioRxiv](https://www.biorxiv.org/content/10.64898/2026.03.03.709355v2.full)
- Rewarded: cytotoxicity, cytokine production and expansion, each normalised to controls. ipTM, pLDDT and ESM-2 likelihood had little predictive power. Synthesis failures clustered in repetitive, high-GC sequences from 5 teams. [bioRxiv](https://www.biorxiv.org/content/10.64898/2026.03.03.709355v2.full)
- The BAGEL-CAR team (414 designs, 38.4% hit rate) says refolding with ESMFold would have discarded every enriched design. [bioRxiv](https://www.biorxiv.org/content/10.64898/2026.09.02.748856v1)

### 1.4 Align to Innovate Protein Engineering Tournament (pilot, 2023)

Enzyme prediction and design, not binders. 28 teams. The generative round was scored on the best passing variant and the median of passing variants. TUM Rostlab won the generative round. Top designs sat within 1 to 4 mutations of the best known variant. The organisers say criteria should match real engineering goals. [PMC](https://pmc.ncbi.nlm.nih.gov/articles/PMC12517241)

### 1.5 Adaptyv x bioArena x muni TREM2 hackathon, San Francisco, 28 Feb 2026

- Format: one day. Hack time 10:45 to 16:00, submissions 16:00 to 17:00. 16 teams (the Adaptyv post says 16), including 6 autonomous LLM agents. Up to 10 designs per team, ranked by Boltz-2 ipSAE_d0res, top 100 of 141 to the lab, top 3 teams by experimental binding. [Luma](https://luma.com/a6t92ohv), [Adaptyv blog](https://www.adaptyvbio.com/blog/agents-vs-humans), [muni](https://www.muni.bio/research/agents-vs-humans-trem2)
- Pipeline numbers: 4,167 binders designed before in-silico filtering, 141 submitted, 100 tested. [muni](https://www.muni.bio/research/agents-vs-humans-trem2)
- Outcome: 37 of 100 bound, 89% expressed. Humans had a 38.5% hit rate and agents 34.3%. Agents had 19% expression failures against 7% for humans, and agent designs were more similar to each other. [Adaptyv blog](https://www.adaptyvbio.com/blog/agents-vs-humans)
- Winners: MRAZS with Mosaic (1.11 and 1.3 nM), BBraiNSEY with BindCraft (1.9 nM), then PXDesign designs from GPT-5.2 (3.64 nM) and Grok 4.1 Fast. All six agents picked PXDesign. [Proteinbase](https://proteinbase.com/collections/adaptyv-x-muni-hackathon-ai-agents-vs-humans)
- ipSAE cut-off 0.547 (Adaptyv default). The top 10 by ipSAE had a 70% hit rate, but ipSAE as a ranker had AUC 0.65. [muni](https://www.muni.bio/research/agents-vs-humans-trem2)
- Jury remarks: none. The only criterion was lab KD and expression.

### 1.6 Berlin Bio x AI Hackathon (Nucleate Germany), 27 to 28 Feb 2026, Adaptyv protein track on 15-PGDH

- Format: 24 hours (Fri 16:00 to Sat 16:00), 76 participants, 17 teams. Public pitch to 200+ guests and a named jury of 6 from academia, industry and a foundation. [biohack.berlin](https://www.biohack.berlin/), [B-TU](https://www.b-tu.de/en/news/article/30704-default-f84d0a64438424f8f64a9b7deb76f660)
- Selection: 6 protein-track teams. A panel picked finalists on "workflow novelty and spinout potential", and 3 finalist teams sent sequences. 117 designs were tested, 91% expressed and 1 bound weakly (1.7 uM, about 0.9%). Methods: BindCraft, BoltzGen, RFdiffusion. [Proteinbase](https://proteinbase.com/collections/berlin-bio-x-adaptyv-15-pgdh-binder-design-competition)
- Winner: the page lists "In Silico" for Protein Design, but its layout is ambiguous. Adaptyv prize: $10K lab credits. [biohack.berlin](https://www.biohack.berlin/)

### 1.7 GEM x Adaptyv RBX1 (15 Feb to 26 Mar 2026, results 26 Apr 2026)

- Rules: de novo only, up to 250 aa, at most 75% identity to known proteins, up to 100 ranked sequences per team plus a 1 to 2 page method PDF. Adaptyv tests each team's top (300 / number of teams) by the team's own ranking. Prize $1,000. [GEM](https://www.gembio.ai/gemadaptyv2026)
- Outcome (team-reported): 1 "Strong" binder among 322 tested designs, next best 85 nM, 20 of 21 designs of one team expressed. Mandrake Bio claims the win (26 nM) with ORBIT (BindCraft framework, Protenix, SolubleMPNN, 12,000+ designs generated, 21 tested, manual "vibe check"). This is self-reported, and no organiser confirmation was found. [Mandrake Bio](https://research.mandrake.bio/p/we-dont-even-design-binders)

### 1.8 HerCode Future of Health Hackathon, ETH, 8 to 10 May 2026

- Format: the official page says 36 hours, 400 applications, 110 accepted, three tracks (Genomics with SOPHiA GENETICS, Drug Discovery with Lemna Bio, AI in Healthcare with USZ) and 13+ jurors and mentors. No scoring, prizes or winners are published. [HerCode](https://futureofhealthhack.hercode.ch/)
- The protein track: a D-BSSE team write-up calls it "Protein Design" and describes a task to design a binder that blocks FGFR2 from its growth factor while avoiding FGFR1. The team focused on the ligand-binding pocket, where only a few residues differ from FGFR1. A Lemna Bio coach's LinkedIn post (seen on a reposting page) says he coached Lemna's "protein design track". No Adaptyv role and no wet-lab step is stated anywhere. [D-BSSE](https://bsse.ethz.ch/news-and-events/d-bsse-news/2026/06/24-hours-of-protein-design.html), [LinkedIn repost page](https://ch.linkedin.com/in/moustafa-houmani)
- That team (self-report): RFdiffusion pipeline plus BoltzGen in parallel. Hours were lost early on models that would not run. They set up runs around midnight, slept, and picked top candidates and wrote slides the next morning. Two members gave a 3-minute pitch. The team placed top 3 in the track and won best pitch overall (self-reported). No thresholds, design counts or wet-lab step are reported. [D-BSSE](https://bsse.ethz.ch/news-and-events/d-bsse-news/2026/06/24-hours-of-protein-design.html)
- Jury: "experts in protein design, startup founders and industry leaders". SOPHiA says judges weighed "innovation, technical execution, impact and more" (general, not track-specific). Gene Hackers won SOPHiA's own track. [D-BSSE](https://bsse.ethz.ch/news-and-events/d-bsse-news/2026/06/24-hours-of-protein-design.html), [SOPHiA GENETICS](https://careers.sophiagenetics.com/blogs/sg-updates-blog/building-the-future-of-health-one-hack-at-a-time)
- Relevance: FGFR2 is a cancer-relevant receptor. The paralogue-avoidance goal (FGFR2 over FGFR1) matches stage 6 of our funnel. If the same host reuses the format, expect a specificity objective.

### 1.9 BioIncubate Catalyst, ETH Student Project House Höngg, 30 to 31 May 2026

- Format: 24 hours, free, CHF 1,600 in cash prizes, launch event 12 May, pitches to a jury of scientists, engineers and ETH faculty. Track 2 "Binder Design and Lead Optimization" (Adaptyv Bio x BioIncubate); winners get wet-lab validation. Compute was provided. No scoring criteria, targets or results are on the organiser pages. [Luma](https://luma.com/4gca12wa), [BioIncubate](https://www.bioincubate.ch/catalyst)
- Reported winner (UNVERIFIED, from LinkedIn snippets that could not be opened): "Basel Binder Designers" (Oliver Nagl, Nathan Bowman, Sofia Filippini, Christoph Luh) won the binding track on an Adaptyv challenge, lead optimisation of cetuximab. First prize reportedly included 50 designs tested at Adaptyv. Nagl and Bowman are also on the FGFR2 team above. [LinkedIn (Nagl)](https://www.linkedin.com/in/oliver-nagl-41a40a1b0/)
- Cofolding track: a mentor's post names Apheris judges. No winners are given. [LinkedIn (Krause)](https://de.linkedin.com/in/falkokrause)
- No official recap was found.

### 1.10 Others

- Nipah de novo portion (benign approach summary only): Escalante Bio's Mosaic pipeline used no inverse folding and no second-model filter. They submitted 10 designs and 9 bound, and say an epitope picked by the optimiser beat a hand-chosen one. [Escalante Bio](https://blog.escalante.bio/winning-the-de-novo-portion-of-the-adaptyv-nipah-binder-competition/)
- Anthropic x Adaptyv (live, 28 Sep to 8 Nov 2026): selection by a Claude workflow that looks at method novelty, design diversity, several confidence metrics rather than one, and several designs per selected submission. Selection is not winning, because lab results decide. [Proteinbase](https://proteinbase.com/competitions/anthropic-adaptyv-2026), [selection page](https://proteinbase.com/competitions/anthropic-adaptyv-2026/selection)
- QBI UCSF hackathon 2026 (48 hours, 19 teams): winners were not binder projects. Two teams built target-preparation tools for BindCraft, RFdiffusion, BoltzGen and Germinal. [QBI](https://qbi.ucsf.edu/hackathon2026)
- Upcoming: Protein Design Hackathon Basel, 21 to 22 Nov 2026, Adaptyv binder challenge, pitches to a jury, criteria unstated. [protein-design.fun](https://protein-design.fun/)

## 2. Cross-event patterns

1. Metrics select and labs decide. iPAE against KD had Pearson r of -0.25 on PDBbind (round 1). Round-2 AUCs were near 0.5. ipSAE had AUC 0.65 on TREM2. ipTM and pLDDT were weak in Bits to Binders. [Adaptyv blog](https://www.adaptyvbio.com/blog/po102/), [Adaptyv blog](https://www.adaptyvbio.com/blog/po104/), [muni](https://www.muni.bio/research/agents-vs-humans-trem2), [bioRxiv](https://www.biorxiv.org/content/10.64898/2026.03.03.709355v2.full)
2. Hit rates span a wide range: 2.5% (EGFR 1), 13 to 14% (EGFR 2), 37% (TREM2), about 0.9% (Berlin 15-PGDH), 1 strong binder in 322 (RBX1, team-reported). Same tools, different targets.
3. Tool diversity among winners: BindCraft (EGFR 1, EGFR 2 de novo, TREM2 third), Mosaic (TREM2 first, Nipah), PXDesign (agents), optimised antibodies (EGFR 2). In EGFR 2, starting from an existing binder gave the best hit rates and de novo hallucination gave 9%.
4. Organisers reward difference. EGFR had hand-picked slots for diversity and novel approaches. The Anthropic x Adaptyv selection prompt names method novelty and diversity. Berlin picked finalists for workflow novelty and spinout potential.
5. Expression and similarity rules remove designs. Round 1 expressed 73% and round 2 95% after SolubleMPNN and a naturalness term. Agent designs failed to express at 19%. Two round-1 binders were disqualified for similarity. RBX1 capped identity at 75%.
6. Time is short. TREM2 gave about 5 hours of hacking. The FGFR2 team lost hours to setup. Agent runs lost designs because humans ran jobs longer. [Adaptyv blog](https://www.adaptyvbio.com/blog/agents-vs-humans)
7. Open data is normal: public-domain release (Bits to Binders), published tested designs including negatives (Anthropic x Adaptyv), a method PDF (GEM).

## 3. Lessons and a pitch checklist (inference)

- Before 31 Oct, run the full pipeline once end to end on a public target (TREM2 and EGFR data are on Proteinbase), with Boltz-2 ipSAE_d0res (PAE and distance cutoffs 15) already scripted. Record ipSAE, ipTM, iPAE, pLDDT and ESM likelihood for every design, since the metric is unknown.
- Fix the environment problems between 24 and 30 Oct, not on the day. The FGFR2 team's only named loss was setup time.
- Assume a cap near 10 designs per team and a lab cut near 100 across teams. Submit a spread of methods, epitopes and clusters, so a miss on one epitope does not sink the set. Hold 2 to 3 slots for a deliberately different approach (for example an optimised known binder, if the rules permit it).
- Filter for expression: SolubleMPNN or equivalent, no repeats, no excess cysteines, length-normalised scores. Check identity to known HER2 binders against whatever similarity rule is announced (10 aa and 75% identity are the precedents).
- Time plan for 10 hours (9:00 to 19:00): generation running by 11:00, selection frozen by 16:00, slides and any method note done by 18:00, submit with margin. TREM2 allowed one hour for final submission. Ours is unknown.
- Pitch (3 minutes, 2 speakers, the FGFR2 team's format): target and epitope rationale with clinical relevance; funnel numbers (generated, filtered, submitted); why these designs differ from each other; the main risk (expression, selectivity among paralogues); code link and seeds for reproducibility.
- Pitfalls to avoid: ranking only by the announced metric; ten near-copies; short sequences that game an unnormalised score; refolding with a single model that disagrees with the scorer; a wrong target construct (read the supplied PDB or chain carefully).

## 4. What could not be verified

- "1 to 5 designs per designer" is absent from the bioRxiv paper (both halves were read through a summariser), the Adaptyv posts and the Boolean Biotech page. The scour.ing mirror cited in the notes did not render. Only "hand-picked for diversity and novel approaches" and "each round-2 designer had at least one design tested" are supported. UNVERIFIED.
- HerCode protein track: no official scoring, results or winner list. The team says 24 h and the organiser says 36 h. Best pitch is also claimed by another team in a LinkedIn snippet that could not be opened. UNVERIFIED.
- Catalyst winners and prizes: LinkedIn only, no recap.
- Berlin protein-track winner and criteria. The official RBX1 winner. Round-1 and round-2 best de novo KD (82 nM in the notes, 91.5 nM in the Adaptyv post).
- Nipah details: one Adaptyv page was blocked by the fetch tool's safety filter and was not pursued. A 9.6% hit rate appears only in a search snippet.
- No documented 2025 to 2026 protein binder hackathon recaps from EPFL, Hugging Face, Ginkgo, Nucleate (other than Berlin), Synbio or iGEM-adjacent events turned up. No other Swiss {ai} Weeks protein hackathon was found.
- LinkedIn and X returned errors for most pages. WebFetch summaries come from a small model and were re-checked only where noted (the D-BSSE page and the Luma page were re-read by the merging session).

## 5. Sources

Primary (organiser or team):

- https://luma.com/4ogriih0 : our event page. Adaptyv challenge named, no criteria.
- https://www.biorxiv.org/content/10.1101/2025.04.17.648362v1.full : EGFR rules, selection, AUCs.
- https://www.adaptyvbio.com/blog/po102/ : EGFR 1 winners, tools, expression.
- https://www.adaptyvbio.com/blog/po104/ : EGFR 2 winners, hit rates by method.
- https://adaptyvbio.substack.com/p/crowdsourcing-better-cancer-drugs : EGFR 1 call and deadlines.
- https://adaptyvbio.substack.com/p/designing-better-cancer-drugs-using : EGFR 1 stats, round 2 kickoffs.
- https://proteinbase.com/competitions/adaptyv-egfr-binder2 : EGFR 2 stages.
- https://www.cradle.bio/blog/adaptyv-protein-design-competition : Cradle process.
- https://www.biomlsociety.org/challenge : Bits to Binders format, winners.
- https://www.biorxiv.org/content/10.64898/2026.03.03.709355v2.full : Bits to Binders analysis.
- https://www.biorxiv.org/content/10.64898/2026.09.02.748856v1 : BAGEL-CAR.
- https://pmc.ncbi.nlm.nih.gov/articles/PMC12517241 : Align tournament.
- https://luma.com/a6t92ohv : TREM2 hackathon schedule and rules.
- https://www.adaptyvbio.com/blog/agents-vs-humans : TREM2 results.
- https://www.muni.bio/research/agents-vs-humans-trem2 : TREM2 pipeline counts, ipSAE cut.
- https://proteinbase.com/collections/adaptyv-x-muni-hackathon-ai-agents-vs-humans : TREM2 leaderboard.
- https://www.biohack.berlin/ : Berlin format, jury, prizes.
- https://proteinbase.com/collections/berlin-bio-x-adaptyv-15-pgdh-binder-design-competition : Berlin results, panel criteria.
- https://www.gembio.ai/gemadaptyv2026 : RBX1 rules.
- https://research.mandrake.bio/p/we-dont-even-design-binders : RBX1 team method (self-report).
- https://futureofhealthhack.hercode.ch/ : HerCode tracks, jury, scale.
- https://bsse.ethz.ch/news-and-events/d-bsse-news/2026/06/24-hours-of-protein-design.html : FGFR2 team (self-report).
- https://luma.com/4gca12wa and https://www.bioincubate.ch/catalyst : Catalyst format.
- https://blog.escalante.bio/winning-the-de-novo-portion-of-the-adaptyv-nipah-binder-competition/ : Mosaic approach.
- https://proteinbase.com/competitions/anthropic-adaptyv-2026 and /selection : current selection method.
- https://qbi.ucsf.edu/hackathon2026 and https://protein-design.fun/ : other hackathons.

Secondary or unverified:

- https://blog.booleanbiotech.com/adaptyv-binder-design-competition : iPAE leaderboard view. It does not contain the 1 to 5 claim.
- https://careers.sophiagenetics.com/blogs/sg-updates-blog/building-the-future-of-health-one-hack-at-a-time : sponsor blog, one-phrase mention of protein design, judging terms.
- https://ch.linkedin.com/in/moustafa-houmani : reposts mentioning Lemna's protein design track.
- https://www.linkedin.com/in/oliver-nagl-41a40a1b0/ and https://de.linkedin.com/in/falkokrause : Catalyst claims (snippets or partial). UNVERIFIED.
- https://www.b-tu.de/en/news/article/30704-default-f84d0a64438424f8f64a9b7deb76f660 : Berlin preview.
- https://www.adaptyvbio.com/blog/anthropic-1 : Adaptyv case study of Claude designs (context only).
- https://intuitionlabs.ai/articles/anthropic-adaptyv-protein-design-competition-guide : third-party guide, September 2026.
