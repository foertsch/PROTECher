# 03. Generation tools for the 31 Oct binder challenge

Research date 2026-10-24. Source: subagent report, merged and edited. Numbers in square brackets point to the source list in section 7. "Authors" means the numbers come from the tool's own paper or company. "Secondary" means the agent read a summary or mirror, not the primary. "UNVERIFIED" marks anything that could not be confirmed.

## Bottom line

1. Run BoltzGen as the main generator and BindCraft as the second one, then add a ProteinMPNN resampling pass on the survivors.
2. Both are MIT-licensed, install with pip or a conda script, and have wet-lab records on cell-surface receptors. BindCraft gave the EGFR round 1 winner and 7 of the 8 de novo binders in round 2 [49]. BoltzGen was the most used tool in a recent Adaptyv competition on a viral target [50] and refolds with Boltz-2, the model Adaptyv used to rank designs there.
3. GPU throughput limits the day more than tool choice. BoltzGen costs about 30 to 125 s per design (secondary timings), so the 60,000 designs per target in its paper equal 500 to 2,100 GPU-hours. One hackathon day yields hundreds to low thousands of designs.
4. All published hit rates are author numbers from much deeper sampling. No independent wet-lab head-to-head of BindCraft, BoltzGen and RFdiffusion3 was found. The one in-silico benchmark (ProtDBench) comes from the PXDesign group.
5. HER2 evidence is thin. A 2025 RFdiffusion study got 2 of 25 constructs binding (19 nM) [53]. One BindCraft binder retargeted an AAV to HER2 without a reported affinity [4]. No HER2 campaign was found for BoltzGen, Chai-2, Latent-X or PXDesign.
6. Chai-2 (partners only) and AlphaProteo (no release) are unavailable. BoltzProt-1 is hosted and paid. RFdiffusion3 has no wet-lab protein-binder data. Germinal, mBER and RFantibody design antibody formats.
7. Spend about 5 hours generating (filtering runs on a rolling basis from hour 2 to 3) and 2 hours on final filtering, with GPU time split roughly 60/40 between BoltzGen and BindCraft.
8. Before 31 Oct, download all weights, build one container per tool, check the NVIDIA driver (BoltzGen main needs a CUDA 13 driver [11]), install PyRosetta for BindCraft, and run a full dry run on HER2 (PDB 1N8Z, see file 02 for the ID check) and EGFR on the real hardware.

## 1. Comparison table

| Tool / version | Reported hit rate (source, who) | Time per design | GPU need | Licence | Install | HER2 evidence |
|---|---|---|---|---|---|---|
| BindCraft, no tagged release; last commit 2026-09-21 [1, 5] | 10 to 100% per target, mean 46.3% over 12 targets (PD-1 13/53, CD45 4/16, PD-L1 7/9). Authors, preprint [4] | About 5 min per 250-residue trajectory on H100; about 100 trajectories (easy target) to 1,000+ (hard) for enough accepted designs [3]; 300 to 3,000 trajectories for about 100 accepted [4] | 32 GB fits about 550 residues, 80 GB about 950 [3]; about 30 GB per 500 residues [4] | MIT [2]; PyRosetta needs a paid licence for commercial use [1]; AF2 weights CC BY 4.0 [7] | conda script, PyRosetta from a RosettaCommons server [6]; Colab works but slow [1, 3] | One HER2 binder (HER2_b1) retargeted an AAV; no Kd or hit rate [4] |
| BoltzGen, PyPI 0.3.2 (2026-05-01), GitHub main newer [8, 12] | Preprint v2 (16 Jun 2026): screening hits on 6/10 novel targets (nanobodies) and 5/10 (miniproteins), 15 tested per target out of 60,000 designs; 4/5 and 3/5 on easier targets. Authors, Adaptyv SPR/BLI [10]. Hit calls were later tightened [16] | About 125 s at about 200 aa (Tamarind, GPU unnamed) [13]; 30 to 60 s on L40S [14]. Both secondary | Min 24 GB, 48 GB advised (secondary) [14]; README gives no VRAM figure [8] | MIT including weights [9, 10] | pip or Docker; about 6 GB weights, main adds 33 GB ESMFold2 and needs a CUDA 13 driver [8, 11] | None found. HER2 has 755 PDB complexes, so it is in-distribution [10] |
| RFdiffusion v1 + ProteinMPNN + AF2 filter | 19% overall, 95 designs per target, 5 targets (authors, Nature 2023) [25]. Independent (30 designs, 6 targets): 2 of 5 Strep-tag binders worked in Western blots, 0 of 15 for three growth factors [26]. PD-L1 minibinders 15/17 bound (2026 study) [54] | Not stated; default 50 steps [24] | Not stated [24] | BSD [24] | conda env with SE3 (CUDA 11.1 yml) or Docker [24] | 2/25 expressed constructs bound domain IV, Kd 19 nM by ITC [53] |
| RFdiffusion3 (Foundry 1.0.0) [20] | No wet-lab protein-binder rate in the preprint; in silico beats RFD1 on 4/5 targets; DNA binder 1/5 (authors) [22]. NVIDIA's phage display comparison: 0.07% specific hits [44] | About 10x faster than RFD2 on A6000 [22] | Not stated [20, 23] | BSD-3 covering weights [21] | `pip install rc-foundry`, official Docker [20] | None found |
| ProteinMPNN / LigandMPNN | Sequence design only. Soluble-MPNN use was credited with better expression in EGFR round 2 [49] | Not stated; cheap relative to folding | Not stated [29, 30] | MIT, parameters included [29, 30] | Clone + pip [29, 30] | n/a |
| Chai-1 / Chai-2 | Chai-2 (authors): antibodies 16% over 52 targets with up to 20 designs each, 26 targets with a hit; miniproteins 68% [32, 33] | n/a | Chai-1: A100/H100 80 GB or L40S 48 GB advised [31] | Chai-1 Apache-2.0 for code and weights [31]; Chai-2 closed | Chai-1 pip; Chai-2 only for selected partners [32] | None found |
| Germinal [34] | 4 to 22% per target (PD-L1, IL-3, IL-20, BHRF1; 43 to 101 designs each), authors [35] | About 1 passing design per H100-hour; 2 to 8 min per iteration [34] | 40 GB+ [34] | Apache-2.0 [36]; PyRosetta and IgLM non-commercial [34] | conda/Docker, manual AF-Multimer params [34] | None (antibody format) |
| mBER [37] | 65 of 145 tested targets (45%) with significant on-design hits; per-target hit rate median about 0.4%, up to 8%; best epitope 8/21 (authors) [38] | Not stated | 32 GB advised, under 16 GB for small targets [37] | MIT code [37] | conda, about 9 GB weights [37] | None seen |
| RFantibody [28] | Not quantified; VHH found in some 95-design sets, about 10k designs may be needed [28] | Not stated | Not stated | MIT [28] | Docker + uv [28] | None |
| PXDesign [39] | 17 to 82% nM hits on 6 of 7 targets; PD-L1 8/11 (authors) [39, 40] | Not in text [40] | No VRAM figure; L20, A800 used [39, 40] | Apache-2.0 code; weights licence unstated [39] | Docker or conda; free web server [39] | None |
| Proteina-Complexa (Mar 2026) [42] | Vendor: PDGFR 122/191 (63.5%); target-specific hits 1.87% vs BoltzGen 0.60% (1.44% with MPNN) on 83 targets [41, 44] | 15.6 s per sample [41] | A100 tested [42] | Code Apache-2.0; weights NVIDIA Open Model License [43] | uv or Docker, Ubuntu 22.04+ [42] | Not mentioned |
| Mosaic [48] | README claims PD-L1 8/10, IL7Ra 7/10, Adaptyv event wins (self-reported, unchecked) | Not stated | JAX GPU | MIT [48] | uv | None |
| Latent-X / X2 [45, 46] | Macrocycles above 90%, 30 to 100 designs per target (authors) [45] | n/a | Hosted | Web platform, early access | Web UI | None |
| AlphaProteo [47] | BHRF1 88% (authors); other targets not given on the page | n/a | n/a | No release found | n/a | Not a target |

## 2. Tool notes

**BindCraft.** One version exists on `main`, and the repo has no releases [1]. Commits since June 2025 are install fixes, a PyRosetta API change and two bug fixes in April 2026 [5]. ProtDBench used the June 2025 commit [52]. The original numbers use 6 to 53 designs per target [4]. The authors say to run until at least 100 designs pass the default filters, then order 5 to 20, and that ipTM predicts binding but not affinity [1, 3]. The Nature version (doi 10.1038/s41586-025-09429-6) was paywalled, and a search-result summary gives 10 to 100% in its abstract. A blog post mentions a "pre-release BindCraft2" [58]. No repository or paper was found (UNVERIFIED).

**BoltzGen.** Pipeline steps are design, inverse folding (SolubleMPNN default), Boltz-2 refolding, analysis and diversity-aware filtering. `--num_designs` of 10,000 to 60,000 is advised for a real run and 50 for a test [8]. The authors tested 15 designs per target out of 60,000 [10]. They use "screening hit" for any SPR or BLI signal and "confirmed binder" for clean sensorgrams [10], and Boltz says some Adaptyv classifications disagreed with a second CRO [16]. Older secondary summaries claiming nanomolar binders for 6 of 9 targets [18] predate that change. Avoid binder lengths 73 to 76, where the model collapses to ubiquitin [10]. BoltzGen main now scores with ESMFold2, which needs a CUDA 13 driver, and the docs describe no switch to turn it off [11]. PyPI 0.3.2 dates from May [12], and whether it already includes ESMFold2 is unconfirmed.

**BoltzProt-1.** Announced 16 Jun 2026 as a BoltzGen successor with a protein-protein affinity scorer. Boltz claims about three times the nanobody hit rate of BoltzGen on ten hard targets, with no counts [17]. It is served through Boltz Lab and the API and is not marked open source. The $100 academic launch credit covered the first two weeks only [17].

**RFdiffusion family.** v1 is BSD, needs the SE3 environment, and the README says to filter at pae_interaction below 10 and to use ProteinMPNN afterwards [24]. RFdiffusion2 targets enzymes and small molecules and has no protein-binder workflow [27]. RFdiffusion3 (in Foundry with RF3, ProteinMPNN, LigandMPNN) is BSD-3 including weights [21], but its preprint reports wet-lab data only for DNA binders and enzymes [22]. The tutorial recommends running MPNN after RFD3 and refolding with RF3 [23]. Wet-lab proof that the RFdiffusion pipeline works on cancer receptors exists for PD-L1 (15 of 17 minibinders, 2 to 230 nM), with lower yields for B7-H3 and B7-H4 [54].

**ProteinMPNN / LigandMPNN.** Both MIT. ProteinMPNN has soluble-protein weights, and LigandMPNN adds ligand-aware, membrane and side-chain packing models [29, 30]. LigandMPNN gives no advantage for a protein-only target.

**Chai.** Chai-1 is an open co-folding model and a candidate second predictor for the filter. Chai-2 goes to selected "academic and industry partners" under a responsible deployment policy, with an early-access sign-up [32]. Treat it as unavailable.

**AlphaProteo, Latent-X.** AlphaProteo has no public release [47]. Latent-X has a web platform with early-access sign-up [45, 46]. A free tier and the output terms could not be confirmed.

**Germinal, mBER, RFantibody.** All three design VHH or scFv formats, so use them only if the challenge asks for antibodies. Germinal needs 40 GB+ and PyRosetta [34], and its authors quote 200 to 500 H100-hours for 200 to 400 passing designs [35]. mBER's per-target median hit rate is 0.4%, rising only under strict ipTM filtering [38]. RFantibody's README says no reliable filter exists and about 10k designs may be needed [28].

**PXDesign, Proteina-Complexa, Mosaic.** PXDesign (ByteDance) has a free web server and Apache-2.0 code [39]. In its paper PD-L1 gave 8 of 11 tested designs, with six more excluded for low expression [40]. It tops ProtDBench, whose authors work at ByteDance [52]. Complexa is the newest (March 2026) and the numbers are NVIDIA's. The weights licence allows commercial use but forbids bypassing safety guardrails [43]. Mosaic gradient-optimises sequences through Boltz-2 and won several Adaptyv-run events by its own account [48]. Its README warns that it often needs hand-tuning and can produce designs that fail simple in-silico tests [48].

**EvoBind2, BindMaster.** EvoBind2 designs peptides from sequence and does not suit a protein-binder challenge (only secondary pages were seen). Nothing was found on a tool called BindMaster.

## 3. HER2 evidence

- RFdiffusion + ProteinMPNN + AF2 screening: 5,000 backbones, 31 passed, 25 cloned, 2 bound HER2 domain IV (the trastuzumab epitope, template 1N8Z), best Kd 19 nM by ITC [53]. The abstract of this MDPI Biomolecules paper claims higher affinity than existing drugs, while its results report lower affinity than the trastuzumab Fab. The merging session re-read the page and confirms the inconsistency.
- BindCraft: one HER2 binder (HER2_b1) gave enhanced cell specificity when inserted into an AAV capsid. The BoltzGen paper lists the BindCraft HER2 target as residues 23 to 268 of 1N8Z [4, 10]. 1N8Z chain C is the full HER2 ectodomain (607 residues, checked against the RCSB API), so that range covers domain I and most of domain II. Whether the BindCraft preprint meant UniProt or PDB numbering (they differ by 22) and why it chose that region is unclear (UNVERIFIED). File 02 has the numbering trap.
- Antibody format: Absci's zero-shot HCDR designs gave 421 SPR-characterised binders, 71 at low nM, but that is antibody CDR redesign with trastuzumab as comparator and not a de novo miniprotein [55].
- Nothing else that was opened tested HER2.

## 4. Hosted options and compute

- Tamarind has a free tier, with no model list or limits on the pricing page [60].
- Rowan runs BoltzGen with a free-tier cap of 10 designs per run (2025 post) [19].
- PXDesign has a free server [39].
- The TREM2 hackathon page lists compute support from Rowan, OpenRouter and Modal [57]. Adaptyv's skill pages describe `modal_bindcraft.py` and `modal_boltzgen.py` wrappers on an L40S [14, 15]. Ask the organisers whether compute will be credits rather than GPUs.

## 5. Head-to-head evidence

- ProtDBench (in silico, single A100, 24 h budget, AF2-IG-Easy filter, ByteDance authors). Successful backbones on PD-L1 were PXDesign 1247, BoltzGen 549, RFdiffusion3 542 and BindCraft 231. On IL-7RA they were 695, 279, 182 and 89. Diffusion tools produced more passing backbones, and BindCraft had higher per-sequence success [52]. HER2 is not a target. Scale by 10/24 for a 10-hour day.
- NVIDIA's compute-matched phage display (83 targets, about 32 GPU-hours per target): Complexa 1.87% specific hits, BoltzGen with MPNN 1.44%, BoltzGen native 0.60%, RFdiffusion3 0.07%, BindCraft 232 specific hits (no rate given) [44]. This is a competitor's paper.
- Adaptyv's skill pages quote $1.80 (BoltzGen) and $2.90 (BindCraft) per accepted design over 7 targets [14, 15]. Adaptyv's primary post was not found.

What Adaptyv competitions say about tools:

- EGFR (2024 to 2025): RFdiffusion was the most common de novo pipeline. BindCraft won round 1 (491 nM from 5 designs) and made 7 of 8 de novo binders in round 2 (best 82 nM from 10 designs; the categories differ from the CSV split in file 01, section 2, so do not compare the hit rates directly). Cradle's cetuximab optimisation won overall at 1.21 nM. The ranking metrics (iPAE, ipTM, ESM2) showed no significant correlation with KD [49].
- Viral-target competition (late 2025, method-level facts only): BoltzGen "emerged as the most popular design tool", RFdiffusion and BindCraft were widely used, and 600 of 1,200 tested designs were picked by Boltz-2 ipSAE. The top team combined BindCraft with TIMED sequence design and Boltz-2 evaluation [50].
- TREM2 one-day hackathon (Feb 2026): the top binder (1.11 nM) came from a Mosaic design. Humans bound 25/65 and agents 12/35 (secondary, and file 05 gives other counts) [48, 56]. Ranking used ipSAE [57].
- Anthropic's August 2026 study put typical campaigns at 10 to 15% and reported 26.8% over 1,320 designs. It names no tools [51]. Anthropic ran it, so treat the field-norm figure as theirs.

## 6. Recommended setup for one 10-hour day

**Tools.** BoltzGen plus BindCraft. Neither needs more than a pip or conda install, and their model families differ (Boltz-2 versus AF2), which gives the cross-model check in the filtering plan a real second opinion. Add Mosaic only if Kira has spare days, since it needs hand-tuning per target. Skip RFdiffusion3 (no wet-lab binder data), Complexa (new, NVIDIA licence terms, little community use) and all closed tools. Keep the PXDesign web server as a free extra batch if time allows.

**Sequence design.** BoltzGen and BindCraft already run SolubleMPNN and MPNN-sol internally. After filtering, resample the top 50 to 100 backbones with soluble ProteinMPNN (8 sequences each, as in ProtDBench [52]), refold with Boltz-2 and keep a sequence only if the metric holds. NVIDIA's data suggest MPNN on BoltzGen backbones helps (0.60% to 1.44% specific hits) [44]. Skip LigandMPNN.

**Budget for a 5-hour generation window.** Rates come from the timings above. They are estimates, and an H100 is likely faster than the L40S behind most of them.

| Hardware | BoltzGen | BindCraft | Notes |
|---|---|---|---|
| 1 x 24 GB | 3 h: about 100 to 300 designs, `--budget` 30 | 2 h: about 15 to 25 trajectories | Trim the HER2 construct hard; by extrapolation of [3], 24 GB fits about 400 residues. Use free PXDesign or Rowan batches for extra candidates |
| 4 x 80 GB | 2 GPUs: about 300 to 1,000 designs, `--devices 2` | 1 GPU: about 60 trajectories | 1 GPU stays free for filtering from hour 2 |
| 8 x 80 GB | 4 GPUs: about 600 to 2,000 designs | 3 GPUs: about 180 trajectories | Run separate BindCraft processes per GPU writing to one folder [3] |

BindCraft yields few accepted designs per GPU-hour: 3 to 30 trajectories per accepted design [4], so expect roughly 2 to 20 accepted designs from 60 trajectories.

**Schedule (09:00 to 19:00).**

- 09:00 to 10:15: target prep (trim the supplied structure to the epitope region, pick one or two epitopes, fill the prepared BoltzGen YAML with `binding_types` and the BindCraft JSON with hotspots). Run a 5-design smoke test per tool.
- 10:15 to 10:30: launch production runs.
- 10:30 to 15:30: generation. From 12:00 a filtering worker scores finished batches (Boltz-2 metrics, second predictor, developability).
- 15:30 to 17:30: final funnel, ProteinMPNN resampling, clustering, epitope review by Maria and Alixe, Iason's physics rescoring on the top 20.
- 17:30 to 18:30: choose the submission set (probably up to 10) and write the method note.
- 18:30 to 19:00: buffer for upload problems.

**Precompute checklist (24 to 30 Oct).**

1. Ask the organisers (or the cluster admin) for GPU type, VRAM, driver version, internet access and any Modal or Rowan credits. Decide today between BoltzGen from PyPI 0.3.2 and from main.
2. Download and cache everything on a persistent volume or inside images: BoltzGen weights (about 6 GB, plus about 33 GB ESMFold2 if on main), AF2 params (about 5.3 GB) for BindCraft, Boltz-2 and Chai-1 weights for the filter, soluble ProteinMPNN weights.
3. Build one Apptainer or Docker image per tool with weights baked in. Install PyRosetta for BindCraft (academic licence) and test the import.
4. Run the BindCraft PD-L1 example and the BoltzGen example (`--num_designs 10 --budget 2`), then 50 designs of each on HER2 domain IV (PDB 1N8Z) and on EGFR. Record seconds per design and accepted counts on the actual GPU. These replace the estimates above.
5. Write the YAML and JSON templates, the trimming script, and the batch runner that sends finished designs to the filter.
6. Create accounts on Tamarind, Rowan and the PXDesign server as fallbacks.
7. Freeze on 30 Oct: one final smoke test on the real machine and no new installs afterwards.

Two cautions for filtering:

- HER2 has many antibody complexes in the PDB (62 entries carry UniProt P04626, re-checked on RCSB on 2026-10-26; an earlier draft quoted 755 entries from the BoltzGen preprint [10], which could not be reproduced and is dropped), so Boltz-2 and AF2 may assign inflated confidence to antibody-like poses on the trastuzumab or pertuzumab epitopes (inference). The caution is qualitative and holds at either count.
- If the judge is Boltz-2 ipSAE, tools that optimise through Boltz-2 will score well in silico. One Adaptyv participant deliberately avoided direct ipSAE optimisation [50].

## 7. What could not be verified

- BoltzGen seconds per design from the paper (Figure 14 values are images). The timings above are Tamarind and Adaptyv mirrors, with no GPU named for the first.
- Per-target design counts for BoltzGen (supplement F.2) and for the RFantibody paper (bioRxiv returned 403).
- Whether PyPI BoltzGen 0.3.2 contains the mandatory ESMFold2 scoring.
- Any wet-lab protein-binder rate for RFdiffusion3, and its GPU memory.
- BindCraft2 (single blog mention [58]) and BindMaster (nothing found).
- Chai-2 access beyond the "partners" statement. A Lilly TuneLab free-trial claim appeared only in a search snippet.
- Latent-X free tier and mini-binder range of 10 to 64% (secondary), and AlphaProteo access status.
- Adaptyv's method-level results for the viral-target competition (the results post was not opened; a skill-page mirror lists very different hit rates per tool, which were not used).
- The claim that all six TREM2 agents chose PXDesign (file 05 cites a Proteinbase page for it, here it rests on a search summary), and the TREM2 hit counts (conflicting secondary sources).
- Nature version of BindCraft (paywalled). The bioRxiv v2 text was used.
- HER2 results for BoltzGen, Chai-2, Latent-X, PXDesign, Germinal and mBER: none found, which is absence of evidence only.
- No independent matched wet-lab comparison of the main generators exists in what was found.
- The challenge's metric, compute, sequence limit and tool rules.
- PXDesign weights licence.

## 8. Sources

[1] https://github.com/martinpacesa/BindCraft : README, MIT, PyRosetta commercial licence note, install, 32 GB Colab note. Primary.
[2] https://raw.githubusercontent.com/martinpacesa/BindCraft/main/LICENSE : MIT, 2024 Martin Pacesa. Primary.
[3] https://github.com/martinpacesa/BindCraft/wiki/De-novo-binder-design-with-BindCraft : timing, memory per residue, trajectory counts. Primary (author wiki).
[4] https://www.biorxiv.org/content/10.1101/2024.09.30.615802v2.full : BindCraft preprint, per-target counts, mean 46.3%, 30 GB, HER2 AAV binder. Primary, authors.
[5] https://api.github.com/repos/martinpacesa/BindCraft/commits?per_page=20 : commit history to Sep 2026. Primary.
[6] https://raw.githubusercontent.com/martinpacesa/BindCraft/main/install_bindcraft.sh : PyRosetta and AF2 parameter sources. Primary.
[7] https://raw.githubusercontent.com/google-deepmind/alphafold/main/README.md : AF2 parameters CC BY 4.0. Primary.
[8] https://github.com/HannesStark/boltzgen : README, install, protocols, design counts, ESMFold2. Primary.
[9] https://raw.githubusercontent.com/HannesStark/boltzgen/main/LICENSE : MIT, 2025 Hannes Stärk. Primary.
[10] https://www.biorxiv.org/content/10.1101/2025.11.20.689494v2 : BoltzGen v2 (16 Jun 2026), wet-lab results, MIT including weights, ubiquitin issue, Table 10. Primary, authors.
[11] https://raw.githubusercontent.com/HannesStark/boltzgen/main/docs/esmfold2.md : ESMFold2 scoring and CUDA 13 driver. Primary.
[12] https://pypi.org/pypi/boltzgen/json : release dates, latest 0.3.2. Primary.
[13] https://docs.tamarind.bio/tools/boltzgen : per-step runtimes, design counts. Secondary.
[14] https://smithery.ai/skills/adaptyvbio/boltzgen : Adaptyv skill mirror, L40S timings, 24 GB minimum, $1.80 per accepted design. Secondary.
[15] https://raw.githubusercontent.com/adaptyvbio/protein-design-skills/HEAD/skills/bindcraft/SKILL.md : BindCraft L40S timings, $2.90 per accepted design. Secondary.
[16] https://boltz.com/experimental-protocol : Boltz note on hit classification and Adaptyv disagreement. Primary (company).
[17] https://boltz.com/boltzmol-boltzprot-api : BoltzProt-1 announcement, API credits. Primary (company).
[18] https://www.tamarind.bio/blog/boltzgen-validated-de-novo-binder-design-for-diverse-targets : earlier BoltzGen numbers (6/9). Secondary.
[19] https://rowansci.substack.com/p/breaking-boltzgen-now-live-on-rowan : Rowan limits, GPU-day cost. Secondary.
[20] https://github.com/RosettaCommons/foundry : Foundry (RFD3, RF3, MPNN), install, Docker. Primary.
[21] https://raw.githubusercontent.com/RosettaCommons/foundry/production/LICENSE.md : BSD-3 covering code and weights. Primary.
[22] https://www.biorxiv.org/content/10.1101/2025.09.18.676967v2.full : RFdiffusion3 preprint. Primary, authors.
[23] https://rosettacommons.github.io/foundry/models/rfd3/tutorials/ppi_design_tutorial.html : RFD3 binder tutorial. Primary.
[24] https://raw.githubusercontent.com/RosettaCommons/RFdiffusion/main/README.md : BSD, install, filters. Primary.
[25] https://www.ebi.ac.uk/europepmc/webservices/rest/PMC10468394/fullTextXML : RFdiffusion Nature 2023 full text (19%). Primary, authors.
[26] https://www.biorxiv.org/content/10.1101/2025.02.07.636769v1.full : independent low-success RFdiffusion study (30 designs). Primary, independent.
[27] https://github.com/RosettaCommons/RFdiffusion2 : scope and licence. Primary.
[28] https://github.com/RosettaCommons/RFantibody : MIT, install, hit-rate caveat. Primary.
[29] https://raw.githubusercontent.com/dauparas/ProteinMPNN/main/LICENSE : MIT. Primary.
[30] https://github.com/dauparas/LigandMPNN : MIT for code and parameters, model variants. Primary.
[31] https://raw.githubusercontent.com/chaidiscovery/chai-lab/main/README.md : Chai-1 Apache-2.0, GPU advice. Primary.
[32] https://www.chaidiscovery.com/news/introducing-chai-2 : Chai-2 results and partner access. Primary, company.
[33] https://www.biorxiv.org/content/10.1101/2025.07.05.663018.full.pdf : Chai-2 preprint. Primary, authors.
[34] https://github.com/SantiagoMille/germinal : README, Apache-2.0, GPU, PyRosetta/IgLM terms, run times. Primary.
[35] https://www.biorxiv.org/content/10.1101/2025.09.19.677421v1.full : Germinal per-target results. Primary, authors.
[36] https://raw.githubusercontent.com/SantiagoMille/germinal/main/LICENSE : Apache-2.0. Primary.
[37] https://github.com/manifoldbio/mber-open : mBER README, MIT, GPU, weights. Primary.
[38] https://www.biorxiv.org/content/10.1101/2025.09.26.678877v1.full : mBER paper. Primary, authors.
[39] https://github.com/bytedance/PXDesign : README, Apache-2.0, install, free server, 17 to 82%. Primary.
[40] https://www.biorxiv.org/content/10.1101/2025.08.15.670450v2.full : PXDesign paper. Primary, authors.
[41] https://research.nvidia.com/labs/genair/proteina-complexa : Complexa results page. Primary, vendor.
[42] https://github.com/NVIDIA-Digital-Bio/Proteina-Complexa : install, A100, workflows. Primary.
[43] https://raw.githubusercontent.com/NVIDIA-BioNeMo/Proteina-Complexa/dev/licenses/license_weights.txt : NVIDIA Open Model License (code Apache-2.0 in licenses/license_code.txt). Primary.
[44] https://www.biorxiv.org/content/10.64898/2026.09.12.751118v1.full : Complexa preprint, compute-matched comparison. Primary, vendor authors.
[45] https://arxiv.org/abs/2507.19375 : Latent-X report abstract. Primary, authors.
[46] https://www.latentlabs.com : product page. Primary.
[47] https://deepmind.google/discover/blog/alphaproteo-generates-novel-proteins-for-biology-and-health-research/ : AlphaProteo results and limits. Primary.
[48] https://raw.githubusercontent.com/escalante-bio/mosaic/main/README.md : Mosaic, MIT, claims. Primary (self-reported).
[49] https://www.biorxiv.org/content/10.1101/2025.04.17.648362.full.pdf : Adaptyv EGFR competition paper. Primary.
[50] https://adaptyvbio.substack.com/p/what-happened-in-the-nipah-protein : tool popularity and top-team pipelines in a viral-target competition (method-level only). Primary (Adaptyv).
[51] https://www.anthropic.com/research/Claude-accelerates-protein-design and https://adaptyvbio.substack.com/p/case-study-benchmarking-claudes-protein : field norm 10 to 15%, 26.8%. Primary, Anthropic and Adaptyv.
[52] https://arxiv.org/html/2605.04118v1 : ProtDBench, in silico, ByteDance authors. Primary.
[53] https://www.ebi.ac.uk/europepmc/webservices/rest/PMC12650588/fullTextXML : HER2 domain IV minibinders (Biomolecules 2025). Primary.
[54] https://www.ebi.ac.uk/europepmc/webservices/rest/PMC13494018/fullTextXML : RFdiffusion minibinders for B7-family cancer receptors (2026). Primary.
[55] https://www.biorxiv.org/content/10.1101/2023.01.08.523187v4 : Absci zero-shot HER2 antibody CDR designs (421 SPR binders). Primary, antibody and not minibinder.
[56] https://scalevise.com/resources/claude-protein-binder-design-trem2-results/ : TREM2 hackathon counts. Secondary.
[57] https://luma.com/a6t92ohv : TREM2 hackathon page, ipSAE ranking, compute sponsors. Primary.
[58] https://blog.escalante.bio/winning-the-de-novo-portion-of-the-adaptyv-nipah-binder-competition/ : mentions a pre-release BindCraft2. Primary, one sentence used.
[59] https://luma.com/4ogriih0 : event page, no tool or compute rules. Primary.
[60] https://www.tamarind.bio/pricing : free and academic tiers, no model list. Primary.
