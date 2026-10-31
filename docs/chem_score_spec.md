# Chemistry score for Track A (spec, not implemented)

Status: soft, rank-only design agreed with Arion on 2026-10-31. The defaults for the four open questions below are mine and need a yes. Nothing here is built or committed.

## Why

The judge's note: the HER2 epitope is rich in Asp/Glu, so the paratope should offer positive groups, and contacts beyond about 4 Å (2.5 Å around O plus 1.5 Å around N on the whiteboard) do not count. ipSAE does not separate the AH06-family designs (parity with the parent, seed noise about 0.1), so a chemistry score ranks them.

Read-only check on the 14 designs (`feasibility_tables.txt` in the shared notes folder, script `chem_feasibility.py` there):

- Epitope (heavy atom within 5 Å of the scFv): 1N8Z has 19 residues, 5 acidic (E558 D560 D570 D596 E598), 2 Lys, net -3. Boltz parent poses: 6 acidic, 2 Lys, net -4. HER2 residues present in 1N8Z sum to net -12; the 480-607 fragment is -6.
- L:D28K reaches 2.1 to 3.1 Å from E598 (in 1N8Z the D28 carboxylate is 4.5 Å from E598, acid against acid). H:D31K and D31R sit 5 to 8 Å from D585 and have no antigen charge within 10 Å in 1N8Z. A hard 4 Å rule would reject them.
- Added Lys/Arg within 4.0 Å of a carboxylate O: 0 or 1 per design, so a count has almost no resolution.
- Some seed-0 poses bind elsewhere (epitope overlap with the parent 0.10 to 0.26 for designs 01, 02, 07, 10). The AH16 designs are off-site in most poses.

## Definition

Inputs per design: the substitution list (crystal numbering), the Boltz-2 model of each seed, the parent's model of the same seed.

1. **Pose filter.** A seed counts only if the Jaccard overlap of its epitope (HER2 residues within 5 Å of the scFv) with the parent's epitope from the same seed is at least 0.5. No eligible seed gives NaN and a flag, never zero.
2. **Contact term C.** For each added Lys or Arg: `s(d)` with d the shortest distance from its charged N (NZ; NE, NH1, NH2) to any HER2 carboxylate O (OD1, OD2, OE1, OE2). `s = 1` for d up to 3.5 Å, falling linearly to 0 at 4.0 Å. His counts 0.5 (ND1, NE2). Substitutions that add no charge, or sit more than 4 Å away, contribute 0 and are never penalised. C is the sum over substitutions.
3. **Charge term Q.** Screened Coulomb energy of the changed groups against the antigen's charged groups, as `E = sum q_i q_j exp(-r/lambda) / r`, with lambda 8 Å and charge sites at the carboxylate C (Asp, Glu, minus one) and NZ or CZ (Lys, Arg, plus one). `Q = E_parent - E_design`, so more favourable is larger. This keeps D31K and D31R in play through steering. This is our own formula, a few lines of numpy; the established alternative (APBS) is heavier than the time allows. The antigen charges come from the 128-residue fragment in the Boltz model only.
4. **Aggregate.** Mean of C and of Q over eligible seeds. Also report the best seed, the number of eligible seeds, and the paratope net charge of the design as a flag (no gate).
5. **Combine with ipSAE.** Only designs within the ipSAE noise band of the parent (mean at least parent minus 2 SD of the seed SD) are ranked by chemistry: `chem = z(C) + z(Q)` over the candidate pool. Outside the band ipSAE decides. The score never vetoes.

## Fit in the code

`design_pipeline/scoring.py` holds `ipsae()` and the METRICS registry. The new metric is `charge_contact`, built on biotite (`CellList`, `distance`), registered the same way, with the per-seed aggregation in the Track A summarize step of `examples/track_a_run.py`. Tests: synthetic two-chain structures (a Lys N placed 3.0, 3.8 and 4.5 Å from a carboxylate O; the parent against itself gives C and Q of 0; an off-site pose is excluded), plus a regression test on one stored Boltz model.

## Open questions with my defaults

1. 4.0 Å is an N to O distance, not a limit on side-chain reach. Default: distance.
2. Charge term is the screened Coulomb above, not FoldX or APBS. Default: Coulomb.
3. Paratope net charge is reported as a flag, with no cap. Default: flag.
4. Pose filter threshold 0.5. Default: 0.5.
