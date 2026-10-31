# Glossary

Short definitions for the Track A brief, merged from the team repo (`docs/track_a_glossary.md` and `docs/research/SUMMARY.md`) on 2026-10-30. If a definition is wrong, fix it here and in the repo.

## The antibody and its parts

- **scFv (single-chain variable fragment):** The two binding domains of an antibody, VH and VL, joined by a linker into one chain. It has no constant domains and no Fc. In this challenge: VH, a 15-residue (G4S)3 linker, then VL, 242 residues in total.
- **VH, VL:** The heavy-chain and light-chain variable domains. Each folds as a sandwich of beta sheets, with three loops at one end.
- **Antigen, epitope, paratope:** The antigen is the target protein (HER2 if the antibody is trastuzumab). The epitope is the part of its surface the antibody touches. The paratope is the antibody's side of that contact.
- **CDR (complementarity-determining region):** The six loops that form the binding site: H1, H2, H3 on VH and L1, L2, L3 on VL. H3 is usually the longest and most variable. Changes here can alter binding a lot. In Adaptyv's public cetuximab scFv data, most variants with edited CDRs lost binding or did not express.
- **Framework:** The stretches of each domain between the CDRs. They form the scaffold, so changes here mostly act on folding, stability and expression. Cradle's winning scFv variant in Adaptyv's EGFR round changed only framework positions.
- **Second shell:** Residues just behind the paratope. They do not touch the antigen but hold the loops in shape, so changing them can alter binding indirectly.
- **VH-VL interface:** The surface where the two domains pack. In an scFv it is the only contact between them. It affects how stable the molecule is and how well it expresses.
- **Linker, disulfide, Twin-Strep tag:** The linker joins VH to VL. Each domain has one internal disulfide bond between two conserved cysteines. Adaptyv adds a 64-residue Twin-Strep tag at the C-terminal end, which captures the scFv on the sensor. We keep the linker and the cysteines unchanged and submit the sequence without the tag.
- **Numbering (Kabat, IMGT):** Conventions that give the same site the same number across antibodies. Papers use different schemes, and the numbers differ from the position in the chain. Map every literature mutation by sequence alignment.
- **Germline, consensus:** The germline is the gene segments an antibody derives from. The consensus is the most common residue at a position among related antibodies. A position where the parent departs from both is a candidate site.
- **Parent, variant, substitution:** The parent (wild type) is the announced scFv. A variant is the parent with substitutions. A substitution is written as original residue, position, new residue, for example VH K30D. Track A needs at least 5 substitutions per variant.

## Binding and the lab assay

- **KD:** The dissociation constant. Lower means tighter binding. It equals koff divided by kon.
- **kon, koff:** The association and dissociation rate constants. A slower koff (the antibody lets go more slowly) is the usual route to a lower KD.
- **BLI, SPR:** The two binding assays Adaptyv uses. The scFv is captured on a sensor and the antigen flows over it.
- **Sensorgram:** The binding signal over time. Adaptyv fits it to get kon, koff and KD.
- **Quantifiable range:** About 0.1 nM to 10 µM. A KD below 1 nM may not resolve well, so koff and expression carry more weight.
- **Expression call:** Expressed, low or not detected. A variant that does not express gives no binding result.
- **Hit rate:** The share of tested designs that bind.

## Scores and models

- **Co-folding model (Boltz-2, Chai-1, AlphaFold):** Predicts the structure of the scFv and the antigen together as one complex. On the known trastuzumab scFv with a HER2 fragment, Boltz-2 gave ipTM 0.94, so a known complex scores high before any mutation.
- **pLDDT:** Per-residue confidence of a predicted structure, 0 to 100, averaged over the design.
- **ipTM, pTM:** Interface and overall TM-score estimates of a predicted complex, 0 to 1, higher is better.
- **PAE, iPAE:** Predicted aligned error, the model's expected position error in Å. iPAE is the mean across the interface, lower is better.
- **ipSAE:** An interface score computed from the PAE matrix that counts only residues whose predicted error to the partner is below a cutoff. 0 to 1, higher is better. It has several variants.
- **AUROC:** The chance that a random binder scores higher than a random non-binder. 0.5 is chance, 1.0 is perfect.
- **Language model (ESM-1b, ESM-1v, ESM-2):** A model trained on protein sequences that gives each amino acid a probability at each position, given the rest of the sequence. It scores how natural a substitution looks. Zero-shot means scoring without any measured data for this antibody.
- **Log-likelihood ratio (LLR):** The log of the probability of the new residue divided by the probability of the parent residue at that position. Above zero means the model finds the substitution at least as natural as the parent.
- **Inverse folding (ProteinMPNN, ESM-IF1, AntiFold):** Models that take a backbone structure and score how well each sequence fits it. AntiFold is trained on antibodies.
- **ΔΔG (FoldX, Rosetta):** The predicted change in stability or binding energy from a substitution. Negative means more favourable. Agreement with measurements is weak to moderate.
- **Seed noise:** The spread of a score when the same input is folded with different random seeds. A cut-off has to sit outside it.
- **RMSD:** The average distance in Å between matching atoms of two structures after superposition. Used to check that a variant keeps the parent's binding pose and loop shapes.
- **Gate, rank, veto:** A gate removes designs below a threshold. A rank orders what is left. A veto is a gate that only removes variants whose structure score falls clearly below the parent's and ranks nothing. Co-folding scores did not rank affinity well in the public data, so veto is the safer use.

## Design terms

- **Epistasis:** The effect of a combination of substitutions differs from the sum of the single effects.
- **Hamming distance, edit distance:** The number of positions at which two sequences of equal length differ. Edit distance also counts insertions and deletions.
- **Frozen positions:** Positions we never mutate, such as the conserved cysteines, the linker and the C-terminal end.
