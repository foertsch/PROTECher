# 02. Likely targets: breast-cancer-relevant extracellular proteins

Research date 2026-10-24. Source: subagent report, merged and edited. Resolutions are in angstrom (A). Residue numbers are UniProt numbering (signal peptide included) unless marked "PDB numbering". "Computed" means the research agent measured heavy-atom contacts (4.5 A cutoff) on the RCSB coordinate files. No paper stated those lists, so treat them as a starting point for hotspot choice and not as published epitopes. "UNVERIFIED" marks anything from a secondary source or that could not be opened.

PDB verification: every ID below was pulled from RCSB's own UniProt cross-reference search by the research agent (RCSB Data API). The merging session re-checked these IDs against the RCSB REST API: 1N8Z, 1S78, 3MZW, 5MY6, 4HRN, 9PI9, 1YY9, 4FRW and 7PEE (titles, methods and resolutions match). The other IDs rest on the agent's check.

## Bottom line

1. No public page names the hackathon target. The Luma page (re-read 2026-10-24) says the event has "a protein design challenge by Adaptyv" and the agent reads its description as "real challenges related to breast cancer". hercode.ch, ai-weeks.ch, the Adaptyv docs and blog, and Proteinbase are silent.
2. Both IDs in the notes are correct. 1N8Z is the HER2 ectodomain with trastuzumab Fab (X-ray, 2.52 A, HER2 chain C has 607 residues, UniProt 23 to 629). 1S78 is the HER2 ectodomain with pertuzumab Fab (X-ray, 3.25 A). The domain IV (trastuzumab) and domain II (pertuzumab) assignments hold.
3. The 15 to 20% figure for HER2-positive breast cancer holds for current usage (a 2020 expert review, PMID 32139271). The 2003 Nature structure paper gave 20 to 30%.
4. Every candidate target below has an entry in Adaptyv's catalogue file, so the catalogue does not separate them. Of the 38 Proteinbase entries in that file, only EGFR and PD-L1 are on our list. HER2, HER3, TROP2 and Nectin-4 do not appear there.
5. Tiers (the research agent's judgment, not measured; section 9 gives the reasons): likely is HER2. Possible are EGFR and TROP2. The long tail is Nectin-4, HER3, PD-L1, FGFR2 and everything else.
6. FGFR2 deserves a place on the list because the previous HerCode protein design track (May 2026) asked for an FGFR2 binder.
7. For HER2 the best-documented design sites are domain IV (trastuzumab site, 579 to 627), domain II (pertuzumab site, 258 to 337), domain I (109 to 186) and a domain III/IV site (394 to 456, 507 to 539) used by an affibody with 22 pM affinity.
8. Numbering trap: 1N8Z, 1S78, 3MZW and 3BE1 number HER2 from the mature protein (UniProt minus 22), while 5MY6, 6LBX and the 2025 cryo-EM entries use UniProt numbering. Always convert before setting hotspots.

## 1. HER2 (ERBB2, UniProt P04626)

### 1.1 Structures

62 RCSB entries carry P04626 (count re-run against the RCSB search API on 2026-10-26: 62). The table shows the ECD complexes among them. All verified.

| ID | Title (short) | Method, resolution | Chains | Partner |
|---|---|---|---|---|
| 1N8Z | HER2 ECD with Herceptin Fab | X-ray 2.52 A | HER2 C (607 aa, UniProt 23 to 629); Fab L A, H B | trastuzumab Fab; NAG, SO4 |
| 1S78 | ErbB2-pertuzumab complex | X-ray 3.25 A | HER2 A, B; Fab L C, E, H D, F | pertuzumab Fab; NAG |
| 3MZW | HER2 ECD with affibody ZHER2:342 | X-ray 2.90 A | HER2 A; affibody B | affibody (58 aa) |
| 2A91 | ErbB2 domains 1 to 3 | X-ray 2.50 A | A | none |
| 5MY6 | HER2-nanobody complex | X-ray 2.25 A | HER2 A; Nb B | nanobody 2Rs15d |
| 7QVK | NM-02 with HER2 ECD | X-ray 3.10 A | HER2 AAA; NM-02 BBB | humanised sdAb |
| 4HRN | DARPin H10-2-G with HER2 domain IV | X-ray 2.65 A | HER2 C, D (103 aa); DARPin A, B | DARPin |
| 4HRL / 4HRM | DARPin 9_29 / 9_26 with HER2 domain I | X-ray 2.55 / 3.20 A | HER2 C / A, C; DARPin A / B, D | DARPins |
| 6LBX | HER2 domain IV with repebody Rb-H2 | X-ray 2.03 A | HER2 B; repebody A | computationally guided repebody |
| 6OGE, 8PWH | HER2 with trastuzumab and pertuzumab Fabs | cryo-EM 4.36 / 3.17 A | HER2 A or E plus 4 Fab chains | both Fabs |
| 7MN5 | HER2/HER3/NRG1b ECD | cryo-EM 2.93 A | HER2-MBP fusion B; HER3 A; NRG1 H | HER3, ligand |
| 9MTE, 9MTX | HER2 ECD with nanobodies 27A05/47D05, 29E09 | cryo-EM 2.66 / 3.06 A | HER2 C / B | VHHs (Mol Cancer Ther 2026) |

### 1.2 Epitopes and reference binders

Domain limits from Pfam in InterPro: I about 23 to 189, II about 190 to 343, III about 344 to 510, IV about 511 to 652.

- Trastuzumab (1N8Z, 3BE1): computed contacts 579 to 583, 591 to 595, 601, 613 to 627 in domain IV (PDB numbering 557 to 561, 569 to 573, 579, 591 to 605). Cho 2003 (PMID 12610629) states it binds the juxtamembrane region. Margetuximab targets the same epitope (Fc engineered), and no separate HER2-margetuximab entry was found.
- Pertuzumab (1S78): computed contacts 258, 267 to 270, 274 to 279, 290 to 292, 307 to 319, 333 to 337 plus two outliers at 150 and 178 (domain II). Franklin 2004 (PMID 15093539) says "near the center of domain II".
- Affibody ZHER2:342 (3MZW): computed contacts 394 to 401, 417 to 424, 449 to 452, 456, 507 to 509, 523 to 527, 538 to 539. The paper (PMID 20696930) reports KD 22 pM and an epitope distant from trastuzumab and pertuzumab. Fab37 (3N85) hits the same face (395 to 401, 417 to 426, 449 to 452, 503 to 509, 523 to 527).
- DARPins: 4HRN (domain IV, computed 534, 543 to 558, 571 to 577, 585) and 4HRL (domain I, computed 109 to 112, 154 to 157, 165 to 170, 178 to 186). G3 binds subdomain IV and 9.26 binds subdomain I (Stuber 2021, PMC8217238). No RCSB entry is named G3, and 4HRN ("H10-2-G") may be the same lineage (UNVERIFIED). Zahnd 2006 (PMID 16963452): several DARPins share the trastuzumab epitope, none the pertuzumab epitope, best Kd 7.3 nM.
- Nanobodies: 2Rs15d (5MY6) recognises domain I (computed 109 to 112, 156 to 157, 170 to 183, 217), Kd 4.74 nM on cells, and does not compete with trastuzumab or pertuzumab. NM-02 (7QVK) also has a non-overlapping epitope. Domain I antibodies chA21 (3H3B) and HuA21 (6J71) contact 157 to 166 and 185 to 197.
- Repebody Rb-H2 (6LBX): domain IV, computed 533, 542 to 548, 558, 572 to 577, 582 to 586, 591.
- De novo: an RFdiffusion plus ProteinMPNN plus AF2 study made domain IV minibinders (Biomolecules 2025, https://pmc.ncbi.nlm.nih.gov/articles/PMC12650588/). Of 25 cloned candidates (31 passed filters), 2 bound on E. coli surface display, and the lead has Kd 19 nM by ITC (called 0_703_6 in the abstract and 2_703_6 in the text). The paper cites 1N8Z as template but describes removing a repebody, so the real template was probably 6LBX (UNVERIFIED). Its abstract and results disagree on how the affinity compares with trastuzumab. The BindCraft preprint uses a HER2 binder for AAV retargeting, and no HER2 success rate was found in it. The Chai-2 preprint text does not name HER2. No other wet-lab-tested de novo HER2 binder was found.

### 1.3 Glycans and constructs

- UniProt: ECD 23 to 652. N-glycans at N68, N124, N187, N259, N530, N571, N629. 25 annotated disulfides.
- Glycans modelled in structures: N187 and N259 (1N8Z), N68 (2A91, 3MZW), N530 and N571 (1S78).
- The 4HRN domain IV fragment carries N8D and N49D (UniProt N530D, N571D), so the DARPin face at 534 to 577 sits next to two glycans. The pertuzumab contact at 258 abuts N259. The trastuzumab epitope ends at 627, two residues before N629.
- Constructs: structure entries use ECD residues 23 to 629 (1N8Z), the full ECD, domain IV fragments, or HER2-Fc fusions (8VQD, 8VQE).
- Vendor ECDs are HEK293 products with C-terminal His tags (Sino 10004-H08H, residues 23 to 652, reseller listing, UNVERIFIED). Sino's domain IV fragment (489 to 630) is biotinylated, which Adaptyv's docs forbid for targets. Adaptyv takes Fc-tagged targets (with an avidity warning) and buys from Acro, Sino or Kactus.

### 1.4 Adaptyv

The catalogue file lists "Isoform 1 of Receptor tyrosine-protein kinase erbB-2 P04626-1" plus rat, mouse and macaque forms. No Adaptyv competition found used HER2 (BenchBB targets: EGFR, IL-7Ra, PD-L1, BBF-14, BHRF1, MBP, Cas9). No page names HER2 for this event.

### 1.5 Difficulty: medium

The ECD is rigid and extended (Cho 2003: fixed, ligand-activated-like conformation). The main epitopes are not glycan-covered, but the domain IV and II sites border glycans. One wet-lab-tested de novo domain IV binder exists (2 of 25, 19 nM). HER2 has the best controls of any target (trastuzumab Fab, affibody, DARPin, nanobody).

## 2. EGFR (P00533)

Structures (all verified):

| ID | Title (short) | Method, res | Chains | Partner |
|---|---|---|---|---|
| 1IVO | EGF with EGFR ECD (Proteinbase reference PDB) | X-ray 3.30 A | EGFR A, B; EGF C, D | EGF |
| 1YY9 | EGFR ECD with cetuximab Fab | X-ray 2.60 A | EGFR A; L C, H D | cetuximab |
| 3B2V | EGFR ECD with IMC-11F8 (necitumumab) Fab | X-ray 3.30 A | A; H, L | necitumumab |
| 5SX4 | panitumumab Fab with domain III | X-ray 2.80 A | M, N; H, I, J, L | panitumumab |
| 4KRL | nanobody 7D12 with domain III | X-ray 2.85 A | A; B | VHH |
| 3QWQ | ECD with adnectin | X-ray 2.75 A | A; B | adnectin |
| 7SYD | extracellular module with EGF | cryo-EM 3.10 A | A, B | EGF |

- Computed epitopes (UniProt = mature + 24): cetuximab 373 to 377, 406 to 408, 432 to 442, 462 to 467, 489 to 497. Panitumumab 372 to 373, 377, 408, 432 to 433, 441 to 444, 462 to 467, 488 to 495. Necitumumab 372 to 373, 408, 432 to 436, 441 to 442, 462 to 467, 488 to 497. All sit in domain III and overlap the EGF site. 7D12: 349, 370 to 382, 408, 433, 442. Adnectin: 36 to 42, 69, 93, 123 to 127, 149 to 152 (domain I). The competition paper lists S440 and G441 as round-2 hotspots, which fall inside the cetuximab/EGF domain III patch under either numbering.
- Adaptyv data: 826 tested and 134 binders (Proteinbase). Best de novo 82 nM and best optimised 1.21 nM against cetuximab at 9.94 nM (bioRxiv).
- ECD 25 to 645. 13 N-glycans (56, 73, 128, 175, 196, 352, 361, 413, 444, 528, 568, 603, 623). 25 disulfides. The bioRxiv paper calls the ECD "easily produced in human cells". The exact Adaptyv construct (tag, Fc) is not stated in the text read.
- In the catalogue and on Proteinbase (Human EGFR). Anthropic x Adaptyv Challenge 1 targeted EGFR (28 Sep to 7 Oct 2026).
- Difficulty: lowest risk (hit rate rose from about 3% to 13% across rounds), but 13 glycans and a heavily mined epitope.

## 3. HER3 (ERBB3, P21860)

- Structures (all verified): 1M6B (HER3 ECD, X-ray 2.60 A, chains A, B); 4LEO (RG7116 Fab, 2.64 A; HER3 C); 5CUS (KTN3379 Fab, 3.20 A; HER3 A to D); 4P59 (MOR09825/LJM716 Fab, 3.40 A; HER3 A); 5O4O (MF3178 Fab, 3.40 A; HER3 C); 3P11 (DL11 Fab, 3.70 A; HER3 A); 7D85 (ISU104 Fab, 2.50 A; HER3 A, D, domain 3); 7BHE / 7BHF (DARPin D5 with domain 4, 2.30 / 2.00 A); 7MN5 (HER2/HER3 cryo-EM).
- Epitopes (computed unless noted): RG7116 domain 1, 36, 67, 91, 111 to 126, 148 to 151, 172 to 181. LJM716 domains 2 and 4 (paper), 266 to 277, 315, 571, 582 to 584, 596 to 602, 609 to 615. KTN3379 310 to 311, 328 to 338, 356 to 358, 419 to 427, 448 to 452, 475 to 477, 481. MF3178 400, 424 to 433, 453 to 455. ISU104 domain 3 (NRG1 site), 335 to 341, 400, 421 to 433, 453 to 455, 483 to 486. DARPin D5 domain 4, residues 568 to 577 (stated in the abstract). Designed binders: computational plus phage-display Affilin proteins (Int J Mol Sci 2025). No de novo minibinder paper was found.
- ECD 20 to 643. 10 N-glycans (126, 250, 353, 408, 414, 437, 469, 522, 566, 616). 23 disulfides. Pfam: L1 56 to 166, CR1 182 to 332, L2 353 to 472, CR2 499 to 629. 4LEO shows glycans at N250, N353, N408, N414, N437, N616.
- In catalogue (P21860, plus mouse). Not in BenchBB or Proteinbase.
- Difficulty: medium to high. The ECD switches between tethered and open states, it carries 10 glycans, and no de novo precedent exists.

## 4. TROP2 (TACSTD2, P09758)

- Structures (verified): 7PEE (ECD, X-ray 2.81 A, chains A to D); 7E5M (trans assembly, 3.20 A, A, B); 7E5N (cis assembly, 3.20 A, A to D); 9PI9 (sacituzumab Fab with Trop2 dimer, X-ray 1.56 A, Trop2 E, F, Fab A to D, Structure 2026). NMR entries of the cytoplasmic tail (2MAE, 2MVK, 2MVL) are irrelevant.
- Sacituzumab epitope: Q237 to Q252 (iScience abstract), computed in 9PI9 as 32, 37 to 38, 66, 174 to 185, 237 to 250. The ECD packs as cis and trans dimers. The paralogue EpCAM (P16422) needs a specificity check. Sacituzumab govitecan is approved for TNBC (iScience abstract). No other binder structures and no de novo binder were found.
- ECD 27 to 274. N-glycans N33, N120, N168, N208. Tg1 domain 70 to 145. 3 annotated disulfides. NAG at N168 is modelled in 7PEE and 9PI9, six residues before the 174 to 185 epitope segment.
- In catalogue (P09758). Absent from BenchBB and Proteinbase.
- Difficulty: medium. Small ECD, few clinical epitopes, a dimer, and glycans at the antibody site.

## 5. Nectin-4 (NECTIN4, Q96NY8)

- Structures (verified): 4FRW (D1-D2 fragment, X-ray 3.50 A, chains A to F); 4JJH (D1, 2.25 A, A, B); 9KKJ (D1 with 9MW2821 Fab, cryo-EM 3.26 A; Nectin-4 A, Fab H, L); 9UCL (IgV with nanobody-based "Trimbody", cryo-EM 2.43 A; Nectin-4 G, H, I). 4GJT also contains Nectin-4, but its partner is a viral glycoprotein, so it was not analysed or used.
- 9MW2821 binds the front beta-sheet of D1. Mutagenesis shows Q77, E78, H83 and E95 matter (abstract). Computed contacts 64, 68, 76 to 97, 116, 128. That face also forms the homodimer. Bicyclic peptides for imaging TNBC exist (PMID 42915462). A phase I/II of 9MW2821 included 20 TNBC patients (PMID 40288679). The enfortumab epitope is UNVERIFIED.
- ECD 32 to 349 (V 32 to 144, C2 148 to 237, C2 248 to 331). UniProt annotates one N-glycan (N281, in the membrane-proximal C2). The V domain has none. 3 disulfides.
- In catalogue ("Nectin-4 Q96NY8.1").
- Difficulty: the lowest glycan burden, a small IgV and a flat face. Only 4 usable entries and no de novo precedent.

## 6. PD-L1 (CD274, Q9NZQ7)

- Structures (all X-ray, verified): 3BIK (PD-1/PD-L1, 2.65 A); 4ZQK (human PD-1/PD-L1, 2.45 A); 5X8L (atezolizumab, 3.10 A); 5X8M (durvalumab, 2.66 A); 5GRJ (avelumab, 3.21 A); 5GGT (BMS-936559 Fab, 2.80 A); 5JDS (nanobody, 1.70 A); 6PV9 (macrocyclic peptide, 2.00 A). 79 entries in total.
- PD-1 face (computed from 3BIK): 18 to 26, 45 to 58, 63 to 66, 113 to 125. Atezolizumab 49 to 69, 76, 111 to 125. Durvalumab 23 to 26, 56 to 61, 111 to 115, 121 to 127. Avelumab 54 to 69, 73 to 78, 113 to 117, 121 to 125. Nanobody 54 to 69, 73, 113 to 115, 119 to 125. All overlap the PD-1 face.
- ECD 19 to 238 (V 19 to 127, C2 133 to 225). N-glycans N35, N192, N200, N219. 2 disulfides. Chai-2 reports picomolar PD-L1 miniproteins (preprint).
- In catalogue, on Proteinbase (Human PD-L1), and in BenchBB.
- Difficulty: lowest, heavily benchmarked. The breast link is indirect and UNVERIFIED here.

## 7. FGFR2 (P21802)

- Structures (X-ray, verified): 1DJS (D2-D3 with FGF1, 2.40 A); 1NUN (FGF10-FGFR2b, 2.90 A); 2FDB (FGF8b-FGFR2c, 2.28 A); 1E0O, 1EV2, 3OJ2; 4WV1 (Fab 2B.1.3 with D2, 2.36 A).
- HerCode's May 2026 track set an FGFR2 ligand-blocking binder task with FGFR1 avoidance (ETH D-BSSE news, https://bsse.ethz.ch/news-and-events/d-bsse-news/2026/06/24-hours-of-protein-design.html).
- ECD 22 to 377. Ig domains 25 to 125, 154 to 247, 256 to 358. 8 N-glycans. Paralogue FGFR1.
- In catalogue (P21802, isoform 3).
- Difficulty: medium. The FGFR1 paralogue is a specificity trap.

## 8. Long tail

All PDB IDs verified by the research agent.

| Target | ECD, glycans (UniProt) | PDB | Notes |
|---|---|---|---|
| FRalpha (P15328) | mature 25 to 234; N69, N161, N201; 8 disulfides | 4KM6 (1.55 A), 4KM7 (1.80), 4KMX (2.20), 4LRH (2.80, folate), 5IZQ (3.60) | no antibody complex in RCSB under this accession; folate pocket buried; epitopes of farletuzumab or mirvetuximab UNVERIFIED |
| B7-H3 (Q5ZPR3) | 29 to 466, 4 Ig domains; N104, N189, N215, N322, N407, N433 | 9LY5 (cryo-EM 2.98, Fab 20G5), 9LY6 (3.67), 9LME (X-ray 2.40, nanobody T3CL11) | T3CL11 binds the IgV face opposite the two glycans (abstract); all 2026 entries |
| CD47 (Q08722) | 19 to 141; N23, N34, N50, N73, N111 | 2JJS (1.85, SIRPa), 2VSC (1.90), 4KJY (1.93, FD6), 5IWL (2.80, 5F9 diabody), 5TZU (2.10, B6H12.2) | SIRPa/antibody face 21 to 26, 45 to 57, 115 to 124 includes the N23 and N50 positions (computed) |
| c-MET (P08581) | 25 to 932; 11 N-glycans; Sema 27 to 515 | 1SHY (3.22), 4K3J (2.80, onartuzumab), 7MO7 (cryo-EM 4.80) | large multi-domain ECD; onartuzumab contacts 326 to 339, 368 to 369, 426, 446 to 449, 469 |
| CDH3 (P22223) | 108 to 654, 5 EC domains; N200, N566 | 6ZTR (2.10, Fab CQY684), 5JYL (2.55, scFv TSP7), 4ZML (1.85) | Ca2+ dependent; CQY684 contacts 124 to 125, 139 to 142, 151 to 164, 168 to 172, 190 to 191 |
| MUC1 (P15941) | repeat region, O-glycosylated | 1SM3 (1.95, SM3 with peptide TSAPDTRPAPGST), 5T78 (2.20), 6BSC (SEA domain, 1.30) | no folded ECD, glycan-dependent epitopes; poor BLI target |
| LIV-1 (Q13433) | 29 to 325; N67, N241, N266, N283 | none | catalogue has Q13433-1; ladiratuzumab vedotin was discontinued after early TNBC responses (PMID 41030069) |

## 9. Likelihood tiers

The research agent first gave percentages (HER2 about 50%, EGFR 15%, TROP2 8%). They rest on no stated method, so they were replaced after review by three tiers with the reason for each. The order inside a tier is not a ranking.

**Likely**

- **HER2.** The iconic breast cancer target. It has positive controls with four non-overlapping epitopes (domains I, II, IV and III/IV), many ECD complexes among its 62 RCSB entries and a vendor ECD, and a "new target" fits an event with no Adaptyv HER2 history. Against it: no document names it, and the same "new target" argument applies to TROP2, Nectin-4 and HER3.

**Possible**

- **EGFR.** Adaptyv's best-characterised target (826 tested on Proteinbase, 1,800+ designs in the competition paper), with cetuximab as control. It is also the first problem in the live Anthropic x Adaptyv series, so material and assay are already in use. Against it: reuse, and the breast link (UNVERIFIED here).
- **TROP2.** Approved TNBC ADC target, with a 1.56 A Fab complex published in 2026. No Adaptyv history.

**Long tail**

- **Nectin-4.** Glycan-light IgV domain and TNBC clinical activity (PMID 40288679).
- **HER3.** HER2 partner with many Fab structures, hard to design against.
- **PD-L1.** Easiest design target and already in Adaptyv's data, with a weak breast link.
- **FGFR2.** HerCode precedent, but organisers may avoid repeating it.
- **Everything else** (FRalpha, B7-H3, CD47, c-MET, CDH3, MUC1, LIV-1, anything unlisted).

The tiers are the research agent's judgment from the evidence above. They are not measured.

## 10. Evidence on the actual hackathon target

- Luma page: host Zeynep Savaseri, no target, rules or sponsors. hercode.ch lists the event and shows no Adaptyv logo.
- ai-weeks.ch/events (3.9 MB page, searched): other 31 Oct events appear, this one does not. A Swiss press report gives Swiss {ai} Weeks 2026 as 1 Sep to 4 Oct (secondary), which conflicts with "officially part of Swiss {ai} Weeks" on 31 Oct. This is worth a question at kickoff, because it affects who sets the rules.
- Adaptyv docs, blog (latest post 25 Aug 2026) and Proteinbase: no HerCode, HER2 or breast cancer mention. x.com (HTTP 402), nitter (403) and LinkedIn profile pages (redirect) were not readable, and search engines returned no relevant post.
- HerCode precedent: Future of Health Hackathon (8 to 10 May 2026), protein design track on FGFR2, challenge written with Lemna Bio, computational only, no Adaptyv mention.
- Adaptyv precedents: TREM2 hackathon with muni (early 2026; Boltz-2 ipSAE_d0res, 10 designs per team, top 100 to lab, from the Adaptyv blog). Berlin 15-PGDH track (search snippet, UNVERIFIED here; file 05 has the Proteinbase page). Catalyst Zurich (30 to 31 May, target not named). Anthropic x Adaptyv: EGFR, TNF-a, GLP1R, TfR1 so far. Challenge 5 (26 Oct to 8 Nov) is unannounced and may reveal a target on Monday, with no evidence it is linked to our event.
- Adaptyv constraints: constructs 8 to 700 aa, cell-free E. coli or yeast, no glycosylation, no full IgG, 20 designs minimum per target, about 150 ug target per run, targets bought from Acro, Sino or Kactus, soluble ECDs only.

Practical consequence of the constraints (inference by the merging session): the binder is made in cell-free E. coli and is not glycosylated, while the target ECD is a glycosylated HEK product. Designs against glycan-adjacent patches (HER2 near N259, N530 and N571, TROP2 near N168) can look good on an unglycosylated model and lose affinity against the glycosylated target.

## 11. What could not be verified

- Event target, metric, length cap, compute. Any X, LinkedIn or Instagram post (blocked).
- The live Adaptyv catalogue (login wall). The agent read the embedding file behind targets.adaptyvbio.com. A listing means orderable, not run.
- Residue lists marked "computed". The Cho 2003 and Franklin 2004 papers are not open access.
- Whether 4HRN is the G3 DARPin. G3 about 90 pM (patent and secondary snippets only).
- "BindCraft succeeds 10% on HER2": no support in the preprint text read.
- Vendor construct details (reseller pages). Epitopes of enfortumab, farletuzumab, magrolimab. EGFR and PD-L1 breast relevance. A 10 to 250 aa length cap for the Anthropic x Adaptyv series (search summary only; file 04 read the primary page).
- Proteinbase targets page: only the first 20 entries were read. The catalogue CSV lists 38 Proteinbase entries (no HER2, HER3, TROP2, Nectin-4, FGFR2).

## 12. Sources

Primary unless marked.

- https://luma.com/4ogriih0 : event text, no target.
- https://www.hercode.ch : event listing. https://futureofhealthhack.hercode.ch/ : May 2026 partners (Lemna Bio, no Adaptyv).
- https://ai-weeks.ch/events : event list without this hackathon. https://www.netzwoche.ch/news/2026-09-03/swiss-ai-weeks-gehen-in-die-zweite-runde : 2026 dates (secondary).
- https://bsse.ethz.ch/news-and-events/d-bsse-news/2026/06/24-hours-of-protein-design.html : FGFR2 track (news).
- https://data.rcsb.org/graphql, https://data.rcsb.org/rest/v1/core/entry/1N8Z, https://search.rcsb.org/rcsbsearch/v2/query, https://files.rcsb.org/download/<ID>.pdb : all PDB data and contacts.
- https://rest.uniprot.org/uniprotkb/<ACC>.json : topology, glycans, disulfides. https://www.ebi.ac.uk/interpro/api/entry/pfam/protein/uniprot/P04626 : Pfam domains.
- PubMed abstracts, https://pubmed.ncbi.nlm.nih.gov/<PMID>/ : 12610629 (Cho 2003), 15093539 (Franklin 2004), 20696930 (affibody), 24095059 (DARPins), 28751451 (2Rs15d), 16963452 (Zahnd 2006), 34992090 (margetuximab), 32139271 (15 to 20%), 41314216 and 34693228 (TROP2), 41101504 (Nectin-4 Fab), 40288679 and 42915462 (Nectin-4 clinical, imaging), 23780344, 26460020, 23928993, 34196609, 33782100, 40429825 (HER3), 42525538 and 41735305 (B7-H3), 23882082 (onartuzumab), 33879555 (CDH3), 41030069 (LIV-1).
- https://pmc.ncbi.nlm.nih.gov/articles/PMC12650588/ : de novo HER2 minibinders. https://pmc.ncbi.nlm.nih.gov/articles/PMC8217238/ : G3 and 9.26 subdomains. https://pmc.ncbi.nlm.nih.gov/articles/PMC10351726/ : NM-02.
- https://www.biorxiv.org/content/10.1101/2024.09.30.615802v2.full : BindCraft. https://www.biorxiv.org/content/10.1101/2025.07.05.663018v1.full : Chai-2. https://www.biorxiv.org/content/10.1101/2025.04.17.648362v2.full : Adaptyv EGFR competition.
- https://docs.adaptyvbio.com/docs/experiment-setup/target-proteins.md, https://docs.adaptyvbio.com/docs/experiment-setup/supported-protein-formats.md, https://docs.adaptyvbio.com/docs/faq/overview.md, https://docs.adaptyvbio.com/llms.txt : constraints.
- https://target-catalog-embeddings.t3.tigrisfiles.io/target_catalog_embedding_2d.csv : catalogue listing (7,455 rows).
- https://proteinbase.com/competitions/anthropic-adaptyv-2026, https://proteinbase.com/targets, https://proteinbase.com/targets/egfr, https://www.adaptyvbio.com/blog/benchbb, https://www.adaptyvbio.com/blog/agents-vs-humans, https://www.adaptyvbio.com/blog : Adaptyv history.
- https://luma.com/4gca12wa : Catalyst track (no target).
- https://www.labome.com/product/Sino-Biological/10004-H08H4-B.html : Sino domain IV fragment (secondary).
