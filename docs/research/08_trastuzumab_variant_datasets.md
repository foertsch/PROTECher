# 08. Public trastuzumab and anti-HER2 variant datasets

Research date 2026-10-31. "Chk" marks something we opened or recomputed in this session. "Agent" marks a claim from the research subagent that we did not open ourselves. "UNVERIFIED" marks anything neither of us opened. Link checks are HTTP status only, made on 2026-10-31, and say nothing about content.

## Bottom line

1. No public dataset gives per-variant measurements for trastuzumab or 4D5 with framework or VH-VL interface mutations at any scale. Agent: every large set mutates CDR-H3 only (Mason 2021, Chinery 2024, Absci 2023). The only measured non-H3 variants are a few CDR-L1 and L3 rows in SKEMPI and about 10 BLI variants of a trastuzumab scFv in the AlphaBind paper.
2. Proteinbase has no HER2 entry. Agent: 45 targets, none is HER2, and its antibody-format data are cetuximab scFv rounds only.
3. The best-documented mutation data are the eight SKEMPI rows of the 1N8Z block (Gerstner 2002). Chk: wild type 0.35 nM by SPR at 310 K. VH D102W is the only clear gain, 3.2-fold tighter (0.11 nM). VH Y105F is 4-fold weaker. The seven-mutation combination is near wild type (0.38 nM), so D102W's gain does not carry into the combination.
4. The large sets suit scorer calibration, not mutation mining. Agent: Chinery has about 524,000 variants with coarse labels (high, medium or low HER2 binding), Mason has about 39,000 with binary labels, Absci has about 421 SPR binders. All three are CDR-H3 only, so like the AbBiBench files already scored they can only calibrate the loop arm.
5. Numbering needs no conversion. Chk: the SKEMPI wild-type letters match the parent sequence at plain sequential positions in each chain (section 2), so VL position plus 135 gives the scFv position.
6. Framework evidence for trastuzumab is thin. Agent: AlphaBind (A-Alpha Bio, 2025) reports about 10 BLI variants of a trastuzumab scFv with edits in FR2 and the CDRH2-FR3-CDRH3 window, parent 7.43 pM by KinExA against our 0.82 nM by BLI. Whether the sequences are public is UNVERIFIED.

## 1. Datasets

| Dataset | Link (check) | Variants | Measured | Mutations | Framework | Licence |
|---|---|---|---|---|---|---|
| Chinery 2024, HER2-aff-large | https://zenodo.org/records/10831512 (403 to a script, UNVERIFIED) | 524,346 after deduplication (Agent) | Yeast display sorted into high, medium and low HER2 binding, NGS (https://www.biorxiv.org/content/10.1101/2024.03.26.586756v1.full, 206) | Up to 10 or more. Binder fraction 98% at 1 edit, 17% at 10 (Agent) | No, CDR-H3 only, IMGT 107-116 (B99-108 by the agent's mapping) | CC BY 4.0 (Agent). data.zip 2.54 GB, raw FASTQ 4-5 GB each |
| Mason 2021 | https://github.com/dahjan/DMS_opt/tree/master/data (206) | 38,733, of which 11,277 binders, per Chinery (Agent) | Mammalian display and FACS, binary labels | Multi-site in 10 CDR-H3 residues, distribution not seen | No (some CDR-L3 files) | Not confirmed. LICENSE.md is linked, raw file 404 (Agent) |
| Mason mirror | https://huggingface.co/datasets/alchemab/her2-binding-prediction (200) | 22,779 rows (Agent) | As above | As above | No | None stated |
| Absci 2023 | https://github.com/AbSciBio/unlocking-de-novo-antibody-design (206) | About 421 SPR binders, plus 13 candidates with Fab and mAb KD (Agent) | SPR KD, ELISA, ADCC, developability | Edit-distance column present, distribution not read | No, HCDR3 only, IMGT 105-117 (B97-109) | BSD (Agent) |
| SKEMPI 2.0, 1N8Z block | https://life.bsc.es/pid/skempi2/database/details/5183 (200), IDs 5176 to 5183; full file https://life.bsc.es/pid/skempi2/database/download/skempi_v2.csv (200, not downloaded) | 8 rows (Chk) | SPR, 310 K, wild type 0.35 nM (Chk) | 1, 3 or 7 (Chk) | No: CDR-L1, L3 and H3 only | CC BY 4.0 (Agent) |
| AlphaBind 2025 | https://pmc.ncbi.nlm.nih.gov/articles/PMC12296056/ (200); code https://github.com/A-Alpha-Bio/alphabind (206) | About 10 BLI variants of a trastuzumab scFv (Agent) | KinExA and BLI, parent 7.43 pM | 2 to 11 across the campaign (Agent) | Yes: CDRH2-FR3-CDRH3 window with many FR2 edits (Agent) | Paper CC BY-NC 4.0, code MIT (Agent). Tables S5 and S7. Sequences UNVERIFIED |
| FLAb | https://github.com/Graylab/FLAb (206) | Aggregator | Its trastuzumab content is the Absci files only (Agent) | n/a | No | Per study |
| AB-Bind | https://github.com/sarahsirin/AB-Bind-Database (206) | Unknown | `1N8Z.pdb` is in the file list (Agent). Whether its experimental table holds trastuzumab rows is UNVERIFIED | Unknown | Unknown | None stated |
| Proteinbase | https://proteinbase.com/download (200); https://storage.proteinbase.com/proteinbase_latest.csv (206) | 45 targets, none HER2 (Agent) | Cetuximab scFv rounds only | n/a | n/a | ODC-BY |

Other notes (all Agent):

- BindingGYM's 4D5-HER2 set (2,080 variants, per the ALLM-Ab paper, PMC12606632) is the AbBiBench file already scored in `docs/track_a_filtering.md`.
- A Hugging Face dataset `w139700701/Trastuzumab` (232 MB, no card, no licence, mismatched CSVs including `mHER_H3_AgNeg.csv`) looks like a copy of the Mason data. Not checked further.
- Chinery's 700 BLI-validated designs: the v1 preprint says validation is ongoing and no data are deposited. A later version may have them.
- The Anthropic x Adaptyv 2026 page the agent read (https://proteinbase.com/competitions/anthropic-adaptyv-2026, 200) lists five challenges (EGFR, TNF-a, GLP1R, TfR1, KRAS G12D pMHC) and no HER2 or antibody track. That differs from the challenge slides photographed on 2026-10-31. Which page applies is unresolved.

## 2. The SKEMPI 1N8Z block (Chk)

Read from the eight detail pages on 2026-10-31. The reference field is PMID 12206766. That this is the Gerstner 2002 paper in `06` is the agent's identification and was not checked. All rows: SPR, 310 K, wild type 3.5E-10 M. SKEMPI numbers positions by the PDB chain (N30 of chain A is VL, D102 of chain B is VH). The scFv position of a VL residue is the chain A number plus 135. Chk: with plain sequential numbering, the parent sequence has D at VH 102, Y at VH 105 and 109, N at VL 30, H at VL 91, Y at VL 92 and T at VL 94, and the scFv has N, H, Y and T at 165, 226, 227 and 229 (and S at 187 for VL S52). These match the wild-type letters in SKEMPI. The location label is SKEMPI's own (COR core, SUP support, RIM rim, INT interior). Its definition is general knowledge, not checked here.

| SKEMPI ID | Mutation (chain, residue) | scFv position | KD (nM) | Against wild type | Location |
|---|---|---|---|---|---|
| 5176 | VL N30S | 165 | 0.39 | 1.1-fold weaker | COR |
| 5177 | VL H91F | 226 | 0.37 | 1.1-fold weaker | SUP |
| 5178 | VL N30S, H91F, Y92W | 165, 226, 227 | 0.31 | 1.1-fold tighter | COR, SUP, COR |
| 5179 | VL T94S | 229 | 0.59 | 1.7-fold weaker | COR |
| 5180 | VH D102W | 102 | 0.11 | 3.2-fold tighter | RIM |
| 5181 | VH Y105F | 105 | 1.4 | 4.0-fold weaker | COR |
| 5182 | VH Y109V | 109 | 0.51 | 1.5-fold weaker | INT |
| 5183 | VL N30S, H91F, Y92W, T94S and VH D102W, Y105F, Y109V | 165, 226, 227, 229, 102, 105, 109 | 0.38 | 1.1-fold weaker | mixed |

Reading (ours, inference):

- D102W is the only single gain, sits at the rim and is in CDR-H3, so it is a loop-arm candidate. It sits in the seven-mutation set, which is near wild type, so a gain here is not additive with the others.
- Y105F is 4-fold weaker, in line with VH Y105 being one of the nine FoldX alanine hotspots (`docs/track_a_filtering.md`). VL N30 is also a FoldX hotspot, yet N30S is neutral. D102 is not a hotspot there.
- The assay is not ours: wild type 0.35 nM by SPR at 310 K against our parent's 0.82 nM by BLI as an scFv. The assay format of these SKEMPI rows is not recorded in the fields read. Effect sizes may not transfer.
- The Dudgeon 2012 stability pair (VH K30D, VL S52D) is not in SKEMPI. Its positions map directly: VH 30, and VL 52 is scFv 187. The effect is not rechecked here.

## 3. Use for Track A (ours, inference)

- Calibrating scorers: Chinery, Mason and Absci are the only large sets with multi-mutation labels. They are CDR-H3 only, so they test the loop arm, like AbBiBench. Chinery's labels are coarse, so use it for binder against non-binder separation and not for ranking.
- Mutation mining: the SKEMPI block is small but clean. AlphaBind is the only framework lead and needs its supplementary tables.
- Not downloaded: Chinery's `data.zip` is 2.54 GB, so a download needs a decision and a place under `/scratch/foertsch`.

## 4. Not found or not verified (Agent)

- No Kd table for humAb4D5 framework variants beyond Carter 1992 (already in `06`).
- The Mason licence, the contents and columns of Chinery's `data.zip`, and the mutation-count distributions for Absci and Mason were not opened.
- Whether the AlphaBind trastuzumab sequences are public, and which FR2 edits they made.
- How many 1N8Z rows SKEMPI holds in total. IDs 5176 to 5183 were read here (Chk). The full CSV would show the rest.
- SAbDab and Thera-SAbDab were not checked (structures, not per-variant measurements).
- Abstract only: Vajdos 2002 (Fab2C4, the pertuzumab lineage, not 4D5), Bostrom 2009 (bH1), Mason 2021.
- Studied other antibodies, not trastuzumab: Koenig 2017 (anti-VEGF G6.31), Kiyoshi 2014 (11K2 and MCP-1), Warren 2025 (4B5, peptide variants only).
- Two broad literature searches (trastuzumab deep mutational scan, affinity maturation) were refused by a tool-side classifier and were not rephrased, so some papers may be missing.
