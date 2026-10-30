# Track A glossary

Terms used in the Track A plan: choose 15 scFv variants of an announced breast cancer antibody, each with at least 5 substitutions, and filter them. Terms already defined in `research/SUMMARY.md` (pLDDT, ipTM, iPAE, ipSAE, PAE cutoff, AUROC, PLL, KD, BLI, SPR, hit rate, gate and rank) are not repeated. Facts and numbers come from `research/04`, `06` and `07` and keep their caveats. "06's own" marks a threshold that file proposes without published support.

## The antibody and its parts

- **scFv (single-chain variable fragment):** The two binding domains of an antibody, VH and VL, joined by a linker into one chain. It has no constant domains and no Fc. Adaptyv's layout is VH, a 15-residue (G4S)3 linker, then VL (07, for its 241-residue cetuximab scFv). Trastuzumab in the same layout gives 242 residues (VH 120, VL 107), which is why it is our leading guess (06).
- **VH, VL:** The heavy-chain and light-chain variable domains. Each folds as a compact sandwich of beta sheets, with three loops at one end.
- **Antigen, epitope, paratope:** The antigen is the target protein (HER2 if the antibody is trastuzumab). The epitope is the antigen surface the antibody touches. The paratope is the antibody side of that contact. In our position table a residue is paratope if it lies within 5 Å of the antigen (06's rule).
- **CDR (complementarity-determining region):** The six loops that form the binding site: H1, H2, H3 on VH and L1, L2, L3 on VL. H3 is usually the longest and most variable. Changes here can alter binding a lot, and in Adaptyv's cetuximab scFv data variants with whole CDRs replaced mostly lost binding (07, rechecked against the Proteinbase table).
- **Framework (FR):** The stretches of each domain between the CDRs. They form the scaffold, so changes here mostly act on folding, stability and expression. Cradle's winning variant changed only framework positions (07).
- **Second shell:** Residues behind the paratope, 5 to 8 Å from it in 06, with no direct contact to the antigen. They hold the loops in shape, so changing them can alter binding indirectly. We count them with the framework slots.
- **VH-VL interface:** The surface where the two domains pack. In an scFv it is the only contact between the domains (06). It affects how the halves pair and how stable and well expressed the molecule is.
- **Linker:** The (G4S)3 stretch joining VH to VL. Frozen. Dimers can form by domain swapping, and 06 advises keeping the linker and the VH-first orientation unchanged.
- **Disulfide bond, conserved cysteines:** Each domain has one internal disulfide between two conserved cysteines. Frozen. A new unpaired cysteine is a known expression risk (04, section 4).
- **Twin-Strep tag:** The 64-residue C-terminal tag Adaptyv appends. The scFv is captured on a Strep-Tactin sensor through it, and the antigen flows over (04). We submit the sequence without the tag and leave the C-terminal end alone.
- **Numbering (Kabat, IMGT):** Conventions that give the same site the same number across antibodies. Papers use different schemes, and the numbers differ from the position in the chain. Map every literature mutation by sequence alignment, as 06 says. ANARCI and abnumber assign the numbers.
- **Germline, consensus:** The germline is the gene segments an antibody derives from. The consensus is the most common residue at a position among related antibodies. A position where the parent departs from both is a candidate site, which is Cradle's approach as 07 reads it.
- **Parent, variant, substitution:** The parent (or wild type) is the announced scFv. A variant is the parent with substitutions. A substitution is written as original residue, position, new residue, for example VH K30D (Kabat, from 06). The mutation count of a variant is its number of substitutions from the parent, and Track A needs at least 5.

## Assay terms

- **kon, koff:** The association and dissociation rate constants. KD is koff divided by kon. A slower koff (the antibody lets go more slowly) is the usual route to a lower KD.
- **Sensorgram:** The binding signal over time as the antigen flows over the captured scFv. Adaptyv fits it to get the two rate constants and KD.
- **Quantifiable range:** About 0.1 nM to 10 µM (04). The BLI concentration series runs from 1000 to 30 nM, so KD below 1 nM may not resolve and koff and expression carry more weight (06).
- **Expression call:** Expressed, low or not detected, read from split-GFP (04). A variant that does not express gives no binding result.

## Scores and models

- **Language model (ESM-1b, ESM-1v, ESM-2):** A model trained on protein sequences that gives each amino acid a probability at each position, given the rest of the sequence. It scores how natural a substitution looks.
- **Zero-shot:** Scoring without training on measured data for this antibody. Adaptyv's po101 blog advises zero-shot language-model scores when the budget is small.
- **Log-likelihood ratio (LLR):** The log of the probability of the new residue divided by the probability of the parent residue at that position. Above zero means the model finds the substitution at least as natural as the parent.
- **Masked marginal:** A way to compute an LLR for every substitution at a site in one pass: hide the residue, read the model's probabilities for all 20 amino acids. It makes scoring all single substitutions cheap.
- **Consensus rule (Hie):** Keep a substitution when its likelihood ratio exceeds 1 (their α = 1) in at least k of six models (ESM-1b plus five ESM-1v), and when at least one model ranks it highest at that site. k was set between 1 and 4 per antibody. 06 uses 2 as the common value. Measured on Fab and IgG, not scFv.
- **Inverse folding (ESM-IF1, AntiFold, ProteinMPNN):** Models that take a backbone structure and score how well each sequence fits it. With the antigen in the input they see the interface as well. AntiFold is trained on antibodies.
- **ΔΔG (FoldX, Rosetta):** The predicted change in stability or binding energy from a substitution. Negative means more favourable. 06 reports weak to moderate agreement with measurements on AB-Bind (abstract-level source).
- **Co-folding model (Boltz-2):** Predicts the structure of the scFv with the antigen together. On the known trastuzumab scFv plus a HER2 fragment, Boltz-2 gave ipTM 0.94, so a known complex scores high before any mutation.
- **Seed noise:** The spread of a score when the same input is folded with different random seeds. A veto threshold must sit outside it. 06 proposes the parent's mean minus 2 SD (06's own).
- **RMSD, pose check:** RMSD is the average distance in Å between matching atoms after superposition. The pose check compares the predicted binding position and loop shapes of a variant with the parent's. 06 proposes rejecting a movement above 1.5 Å (06's own).
- **Veto:** A gate that removes a variant when its structure score falls clearly below the parent's, and ranks nothing. Co-folding scores did not rank affinity in the public data (07) and AF3 ipSAE rose as CDR3 was replaced by alanine in a nanobody benchmark (06), so we use them only this way.

## Design terms

- **Epistasis:** The effect of a combination differs from the sum of the single effects. In Lippow's data two helpful substitutions gave 4-fold where 8.4-fold was expected (06, from a summary).
- **Hamming distance, edit distance:** The number of positions at which two sequences of equal length differ. Edit distance also counts insertions and deletions. We use it for the 5-substitution minimum against the parent and for a minimum distance between designs (06 suggests 3).
- **Frozen positions:** Positions we never mutate: conserved cysteines, the linker and the C-terminal end.
- **Headroom:** How much room the parent leaves to improve. A tight, well-expressed parent has little. Variants of Cradle's winner improved in 1 of 18 cases (07).
- **Arm, portfolio:** The 15 designs together are the portfolio. An arm is the group built by one rule, such as framework and second-shell changes, or loop-periphery changes.
- **Loop slots:** The number of designs that carry 1 to 2 CDR changes. Default 3, set on the day from the announced KD and expression (agreed in discussion, 2026-10-30).
