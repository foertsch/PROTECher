# Article numbers

Every number we would put in writing, with the file it comes from. Claim, value, source,
population, n. No prose, no interpretation: the argument lives elsewhere, this is the
ledger it has to agree with.

Verification state, 2026-11-25: every DESIGNS row in sections 8 to 10, the alanine-scan
SD row in section 9 and the ipSAE rows in section 7 were recomputed from the files named
beside them. Corrections from that pass are marked in bold in place. Sections 2 to 6
(BENCHMARK and CONTROLS) are read straight out of the Part 2 to 4 output CSVs and were
checked value by value against them. The remaining unverified claims are in section 13.

The three files this ledger is checked against, all written from the repository and
none of them recomputing a score:

| file | what it is | written by |
|---|---|---|
| `runs/negative_controls/verdict_states.csv` | the scorer verdict table, one of four states per test | `examples/verdict_table.py` |
| `data/publication/` | the publication CSVs - benchmark, negatives, qualification, plus a README naming each block's estimand. **BENCHMARK only**: no CONTROLS row and no DESIGNS row | `examples/export_publication.py` |
| this file | the claim ledger, including every DESIGNS figure absent from the export | by hand, against the files named in each row |

Populations, used as labels throughout and never mixed in one row:

| label | what it is | size |
|---|---|---|
| BENCHMARK | molecules others measured; scorer test data only | 185 variants: ADAPT 22 (ddG), Moon 7 (KD), Chinery 156 fitted + 97 unfittable. Direction set 21 improvements + 8 losses = 29, of which 13 are singles. Within-position set: 3 positions, 9 improvements, 126 pairs |
| CONTROLS | deliberately-bad variants we generated | 86 |
| DESIGNS | ours | 4,313 singles -> 432 sets -> all 432 FoldX-scored -> 15 reported |

Test states, used in every verdict row: **PASS** = CI entirely above 0.5; **FAIL** = CI
entirely below; **UNPROVEN** = CI spans 0.5; **N/T** = the test is impossible for that
scorer.

## The verdict table

Source: `runs/negative_controls/verdict_states.csv`, written by
`examples/verdict_table.py`, which re-labels the AUROCs and intervals already in the
Part 2 to 4 output CSVs and scores nothing. The notebook's own `verdict.csv` is left in
place as the record of what that run printed; it carries the old two-state labels and
should not be quoted.

| scorer | within position (126 pairs) | state | direction | n | state | clean kinetics | n | state | manufactured fail | wt_argmax |
|---|---|---|---|---|---|---|---|---|---|---|
| `foldx` | 0.8095 (0.7143–0.8889) | **PASS** | 0.7917 (0.5833–0.9583) | 29 | **PASS** | 0.7162 (0.6507–0.7799) | 253 | **PASS** | 8.6% | — |
| `antifold_seq` | 0.2778 (0.1270–0.4365) | **FAIL** | 0.5595 (0.3155–0.8036) | 29 | UNPROVEN | 0.6100 (0.5375–0.6795) | 253 | **PASS** | 2.2% | — |
| `antifold` | 0.3492 (0.2619–0.4286) | **FAIL** | 0.7000 (0.4000–0.9667) | 13 | UNPROVEN | — | 14 | N/T | 3.3% | 66.7% |
| `antiberty` | 0.3175 (0.1667–0.4683) | **FAIL** | 0.3667 (0.0000–0.7667) | 13 | UNPROVEN | — | 14 | N/T | 28.6% | 19.0% |
| `esm1v` | 0.3333 (0.1984–0.4526) | **FAIL** | 0.2333 (0.0000–0.6000) | 13 | UNPROVEN | — | 14 | N/T | 11.0% | 14.3% |
| `esm2_650m` | 0.2619 (0.1667–0.3413) | **FAIL** | 0.4000 (0.0000–0.8000) | 13 | UNPROVEN | — | 14 | N/T | 17.6% | 23.8% |
| `ablang2` | 0.2460 (0.1429–0.3333) | **FAIL** | 0.6667 (0.0000–1.0000) | 13 | UNPROVEN | — | 14 | N/T | 29.7% | 9.5% |
| `blosum62` | 0.5357 (0.4330–0.6390) | UNPROVEN | — | — | N/T | — | — | N/T | — | — |
| `ipsae_max` | — | N/T | 0.4762 (0.2262–0.7321) | 29 | UNPROVEN | — | — | N/T | — | — |
| `consensus` | — | N/T | — | — | N/T | — | — | N/T | — | — |
| `gearbind` | — | N/T | — | — | N/T | — | — | N/T | — | — |

**Every FAIL in the project is a within-position FAIL.** Six scorers, 126 pairs, each
interval entirely below 0.5. No other test in this project produces a FAIL, and that is
a fact about panel width as much as about the scorers: within-position is the only panel
wide enough to put an interval clear of 0.5 on either side.

**The direction test cannot resolve the five matrix scorers.** They score singles only,
so their panel is 10 improvements against 3 losses, **n = 13**, against 21 and 8
(n = 29) for `foldx`, `antifold_seq` and `ipsae_max`. Four of the five intervals touch
0.000 and `ablang2`'s runs 0.000 to 1.000. The smallest p the design can reach is
**0.0035** (`direction/auroc_overall.csv`, `p_floor`). All five are UNPROVEN on this
test and a point estimate from it carries no information - including `antifold`'s 0.700,
which is not evidence of anything.

**`ipsae_max` is UNPROVEN on direction, not FAIL.** 0.4762, CI 0.2262 to 0.7321, n = 29:
the interval spans 0.5. Its exclusion from the ranking rests on **ICC 0.298** over three
seeds on 116 designs (§7), which is a precision argument and is not in this table.

**`blosum62` is UNPROVEN, and that is the conservativeness control's result.** 0.5357,
CI 0.433 to 0.639. A substitution matrix is at chance on this test, so exchangeability
does not explain the six FAILs - the measured improvements are mid-table substitutions,
not conservative ones. The interval was printed by Part 4 and not written to
`conservativeness.csv`; §2 says where it was recovered from.

**The floor checks carry no state.** A failure rate above 25% or a `wt_argmax` below 60%
says the column is unusable for either claim; it is not an interval against 0.5 and
cannot be PASS, FAIL or UNPROVEN. `wt_argmax` here is Part 4's paratope-alphabet rate,
not the by-region rates in §6.

Why a row has no state, where N/T needs qualifying:

- `foldx` — nan
- `antifold_seq` — nan
- `antifold` — kinetics: the Chinery panel returns 14 fitted wells and 0 unfittable for the matrix scorers - one label, so there is no AUROC to have
- `antiberty` — kinetics: the Chinery panel returns 14 fitted wells and 0 unfittable for the matrix scorers - one label, so there is no AUROC to have
- `esm1v` — kinetics: the Chinery panel returns 14 fitted wells and 0 unfittable for the matrix scorers - one label, so there is no AUROC to have
- `esm2_650m` — kinetics: the Chinery panel returns 14 fitted wells and 0 unfittable for the matrix scorers - one label, so there is no AUROC to have
- `ablang2` — kinetics: the Chinery panel returns 14 fitted wells and 0 unfittable for the matrix scorers - one label, so there is no AUROC to have
- `blosum62` — within position: interval computed in the 2026-11-22 run and printed, not written to conservativeness.csv - recovered here from the executed notebook's own output. Re-running Part 4 now persists it. Not tested on direction or kinetics: it ran only as the conservativeness reference
- `ipsae_max` — within position and kinetics: never scored on either panel. Direction is the only test it has, and it is UNPROVEN. Its exclusion from the ranking rests on ICC 0.298 over 3 seeds (116 designs), not on this row
- `consensus` — no route on any box: needs an OAS frequency table that does not exist in this repository. Never run, so nothing is shown either way
- `gearbind` — no route on this box: container image not provisioned. Never run

---

## 1. Lead result: the matching confound

| claim | value | source file | population | n |
|---|---|---|---|---|
| esm1v AUROC, unmatched | 0.875 | `runs/negative_controls/03_negative_controls.executed.ipynb` | BENCHMARK | 126 pairs |
| esm1v AUROC, position-matched | **0.333** | `runs/negative_controls/manufactured/position_matched_auroc.csv` | BENCHMARK | 126 pairs |
| esm1v matched CI | 0.198 to 0.453 | `runs/negative_controls/manufactured/position_matched_auroc.csv` | BENCHMARK | 126 pairs |
| antiberty AUROC, unmatched | 0.806 | `runs/negative_controls/03_negative_controls.executed.ipynb` | BENCHMARK | 126 pairs |

The unmatched number is a position effect. Matching within position removes it and the
scorer lands below chance.

---

## 2. Within-position test — every FAIL in the project originates here

Source: `runs/negative_controls/manufactured/position_matched_auroc.csv`. Population
BENCHMARK, n = 126 pairs for every row.

| scorer | AUROC | CI | state |
|---|---|---|---|
| foldx | 0.8095 | 0.7143 to 0.8889 | **PASS** |
| blosum62 | 0.5357 | **0.433 to 0.639** | **UNPROVEN** - spans 0.5, which is the result: exchangeability does not explain the FAILs |
| antifold | 0.3492 | 0.2619 to 0.4286 | **FAIL** |
| esm1v | 0.3333 | 0.1984 to 0.4526 | **FAIL** |
| antiberty | 0.3175 | 0.1667 to 0.4683 | **FAIL** |
| antifold_seq | 0.2778 | 0.1270 to 0.4365 | **FAIL** |
| esm2_650m | 0.2619 | 0.1667 to 0.3413 | **FAIL** |
| ablang2 | 0.2460 | 0.1429 to 0.3333 | **FAIL** |

blosum62's AUROC is in `runs/negative_controls/manufactured/conservativeness.csv`.
**Its interval is not in that file and was never missing**: Part 4 computed it with
`stratified_ci` and printed it - "BLOSUM62 AS A SCORER on the position-matched test:
0.536 (0.433 to 0.639) over 126 pairs", cell 30 output of
`runs/negative_controls/03_negative_controls.executed.ipynb` - and the write to
conservativeness.csv dropped it. Notebook 03 now persists `ci_low` and `ci_high` to
that file, so the next Part 4 run ends the workaround. **Do not write that it was
never computed.**

Conservativeness control, same file: `rho_vs_blosum` is **-0.2063** for foldx and
**+0.4247** for esm2_650m — the scorer that passes is the one that disagrees with a
substitution matrix.

---

## 3. Direction test — n differs by scorer, and no FAIL originates here

Source: `runs/negative_controls/direction/auroc_overall.csv`. Population BENCHMARK.

| scorer | AUROC | CI | n improvements | n losses | state |
|---|---|---|---|---|---|
| foldx | 0.7917 | 0.5833 to 0.9583 | 21 | 8 | **PASS** |
| antifold | 0.7000 | 0.4000 to 0.9667 | 10 | 3 | UNPROVEN |
| ablang2 | 0.6667 | 0.0000 to 1.0000 | 10 | 3 | UNPROVEN |
| antifold_seq | 0.5595 | 0.3155 to 0.8036 | 21 | 8 | UNPROVEN |
| esm2_650m | 0.4000 | 0.0000 to 0.8000 | 10 | 3 | UNPROVEN |
| antiberty | 0.3667 | 0.0000 to 0.7667 | 10 | 3 | UNPROVEN |
| esm1v | 0.2333 | 0.0000 to 0.6000 | 10 | 3 | UNPROVEN |

**The five matrix scorers see 10 improvements and 3 losses, not 21 and 8**, because they
score singles only. Four of their CIs touch 0.000 and ablang2's spans 0.000 to 1.000:
**this test cannot resolve the matrix scorers at all, and a point estimate from it carries
no information.** The lowest p-value reachable at n = 13 is 0.0035
(`direction/auroc_overall.csv`, `p_floor`).

---

## 4. Clean kinetics

Source: `runs/negative_controls/chinery/auroc_clean.csv`. Population BENCHMARK.

| scorer | AUROC | CI | n fitted | n unfittable | state |
|---|---|---|---|---|---|
| foldx | 0.7162 | 0.6507 to 0.7799 | 156 | 97 | **PASS** |
| antifold_seq | 0.6100 | 0.5375 to 0.6795 | 156 | 97 | **PASS** |
| antifold | — | — | 14 | 0 | **N/T** |
| ablang2 | — | — | 14 | 0 | **N/T** |
| esm1v | — | — | 14 | 0 | **N/T** |
| esm2_650m | — | — | 14 | 0 | **N/T** |
| antiberty | — | — | 14 | 0 | **N/T** |

---

## 5. CONTROLS — beats_positives leads, calls_improving is secondary

Source: `runs/negative_controls/manufactured/failures.csv`. Population CONTROLS.
`beats_positives` is the share of deliberate breaks a scorer ranks above a real
improvement. `calls_improving` is the share it calls an improvement outright.

n-weighted over the five break groups (charge_reversal 16, geometric_hotspot_alanine 10,
glycine_for_buried_bulky 20, mirror 22, proline_in_strand 20; n = 88):

| scorer | beats_positives | calls_improving |
|---|---|---|
| **foldx** | **0.00 in every group** | 7 of 88 = **7.95%** |

Per-group `beats_positives` for foldx is 0.0 at charge_reversal, geometric_hotspot_alanine,
glycine_for_buried_bulky, mirror and proline_in_strand — it never ranks a deliberate break
above a real improvement, in any group.

Worst `beats_positives` by group, for contrast, same file:

| scorer | group | beats_positives | n |
|---|---|---|---|
| esm2_650m | hotspot_alanine | 1.00 | 10 |
| esm2_650m | geometric_hotspot_alanine | 1.00 | 10 |
| esm1v | geometric_hotspot_alanine | 0.90 | 10 |
| esm1v | hotspot_alanine | 0.90 | 10 |
| antiberty | geometric_hotspot_alanine | 0.90 | 10 |
| antifold_seq | mirror | 0.9091 | 22 |

### surface_nearmiss is scored by distance from zero, not by sign

A near-miss at 10 A should score about zero, so "above zero" is a coin flip by
construction and `calls_improving` on this group measures nothing.

| scorer | median | calls_improving | n |
|---|---|---|---|
| **foldx** | **+0.0087** | 0.50 | 24 |
| antifold | -8.4572 | 0.0417 | 24 |
| esm1v | -6.7665 | 0.0000 | 24 |
| antiberty | -6.8922 | 0.0000 | 24 |
| esm2_650m | -5.0337 | 0.1250 | 24 |
| ablang2 | -5.5974 | 0.1250 | 24 |
| antifold_seq | -4.3196 | 0.2500 | 24 |

FoldX's median of **+0.0087 kcal/mol is the correct answer** — it returns approximately
zero for a substitution 10 A from the interface. Its `calls_improving` of exactly 0.50 is
that correctness, not a failure. Every other scorer returns a large negative number for a
mutation that should do nothing.

---

## 6. wt_argmax by region

Source: `runs/negative_controls/wt_argmax_region/wt_argmax_by_region.csv`. Population
BENCHMARK positions. Chance = **5%** (1 of 20 residues).

| scorer | framework rate | framework n | CDR rate | CDR n | verdict in file |
|---|---|---|---|---|---|
| antifold | 0.9833 (177/180) | 180 | 0.7660 (36/47) | 47 | column works |
| esm1v | 0.8889 (160/180) | 180 | 0.3830 (18/47) | 47 | column works |
| antiberty | 0.8389 (151/180) | 180 | 0.4894 (23/47) | 47 | column works |
| esm2_650m | 0.7611 (137/180) | 180 | 0.2553 (12/47) | 47 | column works |
| ablang2 | 0.4889 (88/180) | 180 | 0.2766 (13/47) | 47 | inconclusive |

Every scorer is above chance in CDRs (`cdr_vs_chance` in the file). Every one is far
weaker in CDRs than in framework — the region where the designs are.

---

## 7. ipSAE — excluded on precision, not on the AUROC

| claim | value | source file | population | n |
|---|---|---|---|---|
| ICC(1) | **0.298** | `runs/06_score_correlation_audit.executed.ipynb` | DESIGNS x 3 seeds | 116 designs |
| between-design SD | 0.1591 | same | DESIGNS | 116 |
| seed SD | 0.2442 | same | DESIGNS | 116 |
| SE of the 3-seed mean | 0.1410 | same | DESIGNS | 116 |
| **SE as a share of between-design SD** | **88.6%** | derived from the two rows above | DESIGNS | 116 |
| designs constant across seeds | 0 of 116 | same | DESIGNS | 116 |
| seed-to-seed rank correlation | rho 0.31 | `examples/track_a_chem_rank.py` (comment) | DESIGNS | 116 |
| seed top-15 overlap | 3 to 6 of 15 | `examples/track_a_chem_rank.py` (comment) | DESIGNS | 116 |
| direction-panel AUROC | **0.4762**, CI **0.2262 to 0.7321** | `runs/negative_controls/ipsae/direction_ipsae.csv`, reproduced 2026-11-25 | BENCHMARK | 29 (21 pos, 8 neg) |
| direction-panel state | **UNPROVEN** (CI spans 0.5) | — | BENCHMARK | 29 |

Reproduction method for the row above: the notebook's own `auroc` (average ranks) and
`auroc_ci` (stratified bootstrap, `BOOTSTRAP_SEED = 20261103`, 10,000 resamples, 95%)
applied to `score_ipsae` against `label` in `direction_ipsae.csv`. The brief's CI low of
**0.220** does not reproduce exactly; **0.2262** does. The notebook shares one RNG across
every bootstrap call, so its draw depends on call order - quote 0.226, or re-run Part 2.

**ipSAE's exclusion from ranking rests on ICC 0.298, not on the AUROC.** The AUROC is
UNPROVEN, which is not evidence against it. A confidence measure that does not order point
mutants by affinity is not failing at its own job.

Independent recomputation on the direction panel, a different and smaller population:

| claim | value | source file | population | n |
|---|---|---|---|---|
| ICC(1,1) on `ipsae_max` | 0.2651 | computed from `runs/negative_controls/ipsae/folds.csv` | BENCHMARK direction panel x 3 seeds | 29 designs, 87 rows |
| SE of the 3-seed mean | 0.0317 | same | BENCHMARK direction panel | 29 |
| SE as a share of between-design SD | 69.3% | same | BENCHMARK direction panel | 29 |

Do not quote 0.298 and 0.2651 as the same number. They are two populations.

---

## 8. The 432 rerank

Source: `runs/rerank/ranking_qualified_v1_432.csv`, scores from
`runs/track_a_foldx/foldx.csv`, profile `qualified_v1` = `z(-foldx_ddg)`
(`examples/track_a_chem_rank.py`, `COULOMB_WEIGHT = 0.0`).

| claim | value | population | n |
|---|---|---|---|
| designs FoldX-scored | 432 of 432, 0 failed | DESIGNS | 432 |
| of the 317 previously unscored, better than the old rank-15 cutoff | **0** | DESIGNS | 317 |
| old rank-15 cutoff | k034 at -2.8299 kcal/mol | DESIGNS | 1 |
| best of the 317 | k128 at -2.7614, overall rank 17, short by **0.069** | DESIGNS | 1 |
| new top 15 vs old top 15 | **15 of 15, identical order** | DESIGNS | 15 |
| of the 32 "ours" in the new top 15 | **0** | DESIGNS | 32 |
| best "ours" rank | 71 (o13) of 432 | DESIGNS | 32 |
| "ours" rank range | 71 to 176, median **119.5** | DESIGNS | 32 |
| kira median rank | **232.5** | DESIGNS | 400 |
| ddg > 0 (weaker than wild type) | 261 of 432 | DESIGNS | 432 |
| ddg > 0, the 115 | 29 of 115 | DESIGNS | 115 |
| ddg > 0, the 317 | 232 of 317 | DESIGNS | 317 |
| ddg > 0, "ours" | 2 of 32 | DESIGNS | 32 |
| Spearman(priority_or_additive, measured ddg) | **-0.3219**, p 4.5e-4 | DESIGNS, the selected 115, all source=kira | 115 |

The Spearman is computed inside the selected 115, so it is range-restricted and
understates the proxy's coarse discrimination. Across the selection boundary the same
proxy separated well: median ddg **-1.344** for the 115 against **+1.721** for the 317.

### Predicted improvement rates

Source: `data/track_a_final/candidates.csv` joined to `runs/track_a_foldx/foldx.csv`,
`ddg < 0`. Partition is exact: 32 + 115 + 285 = 432, no overlap.

| group | predicted to improve | rate | population | n |
|---|---|---|---|---|
| "ours" | 30 of 32 | **94%** | DESIGNS | 32 |
| kira, selected for the first FoldX cut | 86 of 115 | **75%** | DESIGNS | 115 |
| kira, not selected | 55 of 285 | **19%** | DESIGNS | 285 |

**PREDICTED, by FoldX. Not measured.** Nothing in this project measures the affinity of
any design.

---

## 9. FoldX precision on the design population

Source: `runs/track_a_foldx/foldx_sd_sample.csv` (n_runs = 3, 60 candidates) and
`runs/track_a_foldx/foldx_runs20_top50.csv` (n_runs = 20, top 50). Scripts in
`examples/foldx_precision/`.

This is a result about the candidate population, not about FoldX. FoldX's documented
accuracy against experimental ddG is roughly 0.8 to 1.0 kcal/mol RMSE, so a run-to-run SD
of 0.56 is smaller than its known error, and BuildModel's stochasticity is the reason
`numberOfRuns` exists. Precision and accuracy are separate: the AUROC rows in sections 2
to 4 measure whether the mean tracks measured ddG and they stand. A scorer can be valid
and imprecise at once.

| claim | value | source file | population | n |
|---|---|---|---|---|
| median per-candidate SD, 3 runs | **0.5643** | `foldx_sd_sample.csv` | DESIGNS | 60 |
| p90 per-candidate SD, 3 runs | **1.0928** (the brief's 1.1340 is `method='higher'`; linear interpolation gives 1.0928) | `foldx_sd_sample.csv` | DESIGNS | 60 |
| per-candidate SD range, 3 runs | 0.0327 to 1.8000, mean 0.5993 | same | DESIGNS | 60 |
| candidates with SD = 0, 3 runs | 0 of 60 | same | DESIGNS | 60 |
| median per-candidate SD, 20 runs | 0.6470 | `foldx_runs20_top50.csv` | DESIGNS top 50 | 50 |
| median SE of the mean, 3 runs | 0.326 | derived, SD/sqrt(3) | DESIGNS | 60 |
| median SE of the mean, 20 runs | **0.1447** | `foldx_runs20_top50.csv` | DESIGNS top 50 | 50 |
| median SD, single alanine mutants | **0.0000** (mean 0.0404, p90 **0.0702**, exactly zero in **151 of 210**) | `trastuzumab_foldx_alanine_scan.csv`, `ddg_bind_sd`, `status == ok` (227 rows, 17 `wt_is_ala`) | BENCHMARK positions | 210 |

**The single-mutant precision estimate is off by more than an order of magnitude for 5 to
7 substitution designs: 0.0000 median against 0.5643.** Mechanism: more substitutions give
BuildModel more to rearrange. Same edit-distance problem that bit the scorers in section 2.

### The designs are packed tighter than the scorer's resolution

| claim | value | source file | population | n |
|---|---|---|---|---|
| distinct positions the 432 move | **23** (chain A 10, chain B 13) | `runs/track_a_foldx/candidates_all432.csv` | DESIGNS | 432 |
| ddG span, top 22 | **1.556** kcal/mol | `ranking_qualified_v1_432.csv` | DESIGNS | 22 |
| ddG span, top 30 | 1.802 | same | DESIGNS | 30 |
| ddG span, top 50 | 2.402 | same | DESIGNS | 50 |
| ddG span, all 432 | 12.115 | same | DESIGNS | 432 |
| ddG span, ADAPT | **2.790** (-2.040 to +0.750) | `data/benchmarks/trastuzumab_mutation_evidence.csv`, `ddg_reported_kcal_mol` | BENCHMARK | 22 |
| designs vs benchmark at matched n = 22 | **1.79x more tightly packed** | derived | — | 22 vs 22 |
| median adjacent gap, top 30 | 0.0363 | `ranking_qualified_v1_432.csv` | DESIGNS | 30 |
| SE of a 3-run mean as a share of the top-30 span | 18% | derived | DESIGNS | 30 |

The practice this supports: **measure a scorer's precision on your own candidate
population and compare it to the spread you are trying to resolve.** One short run. The
60-candidate n=3 measurement took 1,027 s and the 50-candidate n=20 measurement 4,570 s.

---

## 10. What the 15 reported designs are

| claim | value | source file | population | n |
|---|---|---|---|---|
| designs submitted | 15 | `runs/rerank/ranking_qualified_v1.csv` | DESIGNS | 15 |
| band they were drawn from, 3 runs, median SD | **39 designs, ranks 5 to 43** | derived from `foldx_sd_sample.csv` + `ranking_qualified_v1_432.csv` | DESIGNS | 432 |
| band width, 3 runs | 1.7264 kcal/mol | same | DESIGNS | 39 |
| designs within 1 SD of the rank-15 cutoff, 3 runs | 26 | same | DESIGNS | 432 |
| designs within 2 SD of the rank-15 cutoff, 3 runs | 51 | same | DESIGNS | 432 |
| adjacent ranks separable only if SD < | **0.0227** | derived | DESIGNS | 30 |

15 designs were submitted, drawn from a band of 39 (ranks 5 to 43) that FoldX at 3 runs
cannot distinguish. Per-candidate SD median 0.5643 kcal/mol; SE of a 3-run mean 0.326
against a top-30 spread of 1.802.

### Revised at 20 runs

Source: `runs/track_a_foldx/foldx_runs20_top50.csv`. Every test uses each candidate's own
SE: resolvable when `|mean_i - mean_j| >= 1.96 * sqrt(se_i^2 + se_j^2)`.

| claim | value | population | n |
|---|---|---|---|
| **Spearman(3-run order, 20-run order)** | **+0.7567**, p 2.1e-10 | DESIGNS top 50 | 50 |
| Kendall tau | +0.5739 | DESIGNS top 50 | 50 |
| ranks that changed | 46 of 50 | DESIGNS top 50 | 50 |
| top-15 membership change | 11 of 15 kept | DESIGNS top 50 | 15 |
| entered the top 15 | k027, k047, k069, k128 | DESIGNS | 4 |
| left the top 15 | k028, k032, k034, k052 | DESIGNS | 4 |
| **the old rank-15 cutoff, k034, at 20 runs** | **rank 21** | DESIGNS top 50 | 1 |
| largest moves | k039 20->49, k041 19->42, k032 10->30 | DESIGNS top 50 | 3 |
| tied band around k034, own SEs | **21 of 50**, ranks 11 to 35, 0.904 kcal/mol wide | DESIGNS top 50 | 50 |
| pairs resolvable | 856 of 1,225 = **69.9%** | DESIGNS top 50 | 50 |
| resolvable tiers, uncorrected | **8** - notebook 07's own clustering; **not reproduced** by chaining adjacent pairs at z = 1.96 on each candidate's own SE, which gives 1 tier (every adjacent gap in the top 50 is smaller than 1.96*sqrt(se_i^2+se_j^2)). Quote 8 only with 07's definition beside it, or quote the pair count instead | DESIGNS top 50 | 50 |
| resolvable tiers, Bonferroni z = 4.103 | 5, top tier 9 - same caveat, adjacent chaining gives 1 | DESIGNS top 50 | 50 |

**Top tier, 5 designs**, all source=kira, and **all five carry `A:S50R`** (checked against the `variant_id` column of `foldx_runs20_top50.csv`):

| name | mean ddG | SD | SE | variant |
|---|---|---|---|---|
| k003 | -3.850 | 0.647 | 0.145 | A:D28W+A:S50R+A:Y92W+B:Y52R+B:Y105W |
| k006 | -3.839 | 0.775 | 0.173 | A:D28W+A:S50R+A:Y92W+B:Y52R+B:D102W |
| k049 | -3.599 | 0.681 | 0.152 | A:S50R+A:Y92W+B:Y52R+B:G101F+B:M107R |
| k010 | -3.571 | 0.422 | 0.094 | A:D28W+A:S50R+B:D31R+B:Y52R+B:Y105W |
| k001 | -3.498 | 0.977 | 0.219 | A:S50R+A:Y92W+B:Y52R+B:R59F+B:Y105W |

### Selection bias in the reported values

| claim | value | source file | population | n |
|---|---|---|---|---|
| mean shift, 3-run mean to 20-run mean | **+0.2565** kcal/mol worse | `foldx_runs20_top50.csv` | DESIGNS top 50 | 50 |
| candidates that got worse | 38 of 50 | same | DESIGNS top 50 | 50 |
| corr(SD, shift) | **+0.699** Pearson, +0.659 Spearman | same | DESIGNS top 50 | 50 |
| mean shift, SD < 0.4 | -0.046 | same | DESIGNS top 50 | 11 |
| mean shift, SD > 0.9 | **+0.762** | same | DESIGNS top 50 | 11 |
| best ddG, 3 runs vs 20 runs | **-4.182 -> -3.850** | same | DESIGNS | 1 |

Selecting the top 50 on a 3-run mean promoted candidates whose three runs happened to land
low, and those are the high-variance ones; well-determined candidates did not move. The
reported best ddG was optimistic by 0.33 kcal/mol by construction of the selection.

---

## 11. D98W

| claim | value | source file | population | n |
|---|---|---|---|---|
| affinity improvement | **4.3x** (KD 0.11 nM against the Moon parent KD of row `TZM-014`) | `data/benchmarks/trastuzumab_mutation_evidence.csv`, row **`TZM-019`**, `affinity_fold_reported` | BENCHMARK | 1 |
| NCI-N87 activity ratio | **0.48** | same file, row **`TZM-023`**, `metric = relative_antiproliferative_activity`, `units = IC50(WT)/IC50(variant)` | BENCHMARK | 1 |
| SK-BR-3 activity ratio, same variant | 0.89 | same file, row `TZM-022` | BENCHMARK | 1 |

A 4.3x affinity gain against an activity ratio below 1 in the same molecule, on both
cell lines. The two readouts disagree in direction.

**Two values that both read 0.48 are easy to collide. Never put both in one
sentence.** Moon's wild-type parent KD is **0.48 nM** (`TZM-014`), in nanomolar, and
it is what every Moon fold change is measured against. D98W's NCI-N87 relative
antiproliferative activity is **0.48** (`TZM-023`), a unitless IC50 ratio on a cell
line. The affinity claim is `TZM-019`: KD 0.11 nM against that parent, reported as
4.3x.

---

## 12. The liability rule merge

Source for every row: `examples/liability_union_check.py`, which reads the chain sequences
and the disulfide map off `data/structures/1N8Z_fv_trimmed_Repair.pdb` exactly as
`_gates.LiabilityCheck.run` does, so it measures the shipping rule. Position sanity: the
`wt` column agrees with the structure for 4,313 of 4,313 singles.

| claim | value | population | n |
|---|---|---|---|
| flagged by the merged rule | **511** (11.8%) | DESIGNS singles | 4,313 |
| flagged by the narrower copy | 372 | DESIGNS singles | 4,313 |
| disagreements | **145** (142 merged only, 3 copy only) | DESIGNS singles | 4,313 |
| motifs the copy lacked | isomerization 82, deamidation 40, fragmentation 21 | DESIGNS singles | 4,313 |
| merged rule by reason | free_cysteine 223, isomerization 83, disulfide_break 76, deamidation 68, n_glycosylation 53, fragmentation 21 | DESIGNS singles | 4,313 |
| variable-domain cysteines, all paired | A 23/88, B 22/96 | structure | 4 |
| 29-panel failures | **1 of 29**: `A:T93I+A:P95N+A:P96I`, reason `A:n_glycosylation@95`, and it is a **loss** | BENCHMARK | 29 |
| measured improvements vetoed | **0 of 21** | BENCHMARK | 21 |
| the 432 design sets failing | **0** | DESIGNS | 432 |
| the shipped 115 failing | 0 | DESIGNS | 115 |

The union does not re-filter the pool: the merge is a correctness fix, not a change to what
was selected.

**The merged rule has a blind spot the narrower copy did not have.** Three singles are
flagged by the copy and not by the merge - `A:T31G`, `B:T78G`, `B:S85G`. Each turns `NT` or
`NS` into `NG` at an asparagine that already carried a motif. `new_motifs` compares
occurrences by start position, and the merged deamidation class `N[GST]` matched before the
substitution as well as after, so the position does not move and nothing is reported. `NG`
is the high-risk deamidation motif and `NT`/`NS` are weaker, so this is a real increase in
liability at a site the rule is already watching. Widening a motif class while matching
occurrence-wise hides changes in severity within the class.

These numbers do not reproduce the ones carried in the brief and in `_gates.py`'s own
docstring - 437 vs 352 with 233 disagreements, split 158 motifs and 75 cysteine calls.
Neither total matches, and the mechanism does not either: all four variable-domain cysteines
are paired, so the union's `disulfide_break` catches every cysteine-count change the copy
catches and the copy makes **no** cysteine calls the union misses. The docstring's claim of
75 is not reproducible on this structure and this population.

---

## 13. Untraceable

Carried from the brief, not reproduced from a file in this repository. Each needs a source
or removal before publication.

| claim | value as given | status |
|---|---|---|
| hydrophobic-exposure gate | **not repeated here, by design** | **moved out of this ledger 2026-11-26.** Four values (an AUROC in the wrong direction, two group means, an RSA correlation, three area figures) were computed 2026-11-23 in scratch scripts (`q3b.py`, `q3c.py`) under the laptop's Windows temp directory, during a diagnostics-only pass that was instructed to commit nothing; that directory has since been cleaned, so none of them can be re-checked against a file. They live in `docs/filtering_strategy.md`, marked NOT REPRODUCED, where their only job is to stop the gate being rebuilt. **Nothing in the article cites them, and they are not restated in this ledger so that nothing here can be mistaken for a traceable number.** |
| ipSAE direction-panel AUROC | 0.476, CI 0.220 to 0.732 | **resolved 2026-11-25**: reproduced from `runs/negative_controls/ipsae/direction_ipsae.csv` as 0.4762, CI 0.2262 to 0.7321, now in §7. Only the brief's CI low (0.220) is wrong |
| blosum62 within-position CI | 0.433 to 0.639 | **resolved 2026-11-26**: printed by Part 4 of the 2026-11-22 run and recovered from cell 30 of the executed notebook. UNPROVEN rests on the interval, not on the point estimate. `conservativeness.csv` carries it after the next Part 4 run |
| CONTROLS size | 86 variants | does not reconcile with `failures.csv`, whose break groups sum to 88 per scorer and 112 with the neutral groups |
| "the four bugs" as a set of four | — | **the count is untraceable**: nothing in the repository enumerates four. Five fixes with a before/after are traceable and are listed in §14; a sixth has no number. Pick four from there or drop the count |
| D98W NCI-N87 ratio | 0.48 | **resolved 2026-11-26**: row `TZM-023` of the evidence file. See §11 for the collision with Moon's 0.48 nM parent KD |
| liability union count over the singles | 437 vs 352, 233 disagreements, split 158 motifs / 75 cysteine | **contradicted by measurement**, see §12: measured 511 vs 372 with 145 disagreements and 0 cysteine-only calls. `_gates.py`'s docstring carries the unreproducible figures and should be corrected |

---

## 14. The bug fixes, with before and after

Traceable from the working tree of `final_cleanup` (uncommitted) and from `CHANGELOG.md`.
Five, not four. No number in sections 1 to 11 changes because of any of them except the
FoldX SD rows, which exist only because of the first.

| # | bug | before | after | source |
|---|---|---|---|---|
| 1 | FoldX's per-run spread was computed and thrown away | no precision column; two candidates 0.001 kcal/mol apart were indistinguishable from two 1.0 apart | `ddg_sd` written per candidate, empty at `FOLDX_RUNS=1`; measured median **0.5643** at n=3 over 60 candidates | `design_pipeline/adapters/_foldx_score.py`, `score_candidate` |
| 2 | the FoldX total-energy column was named for the wrong quantity | `ddg_fold` | `ddg_folding`, with `read_scored()` rewriting an older table in place under the current header | same file, `FIELDS` / `RENAMED` |
| 3 | two copies of the liability rules, the narrower one shipping | the copy flagged **372** of 4,313 singles; the docstring claimed 437 vs 352 with 233 disagreements, split 158 motifs / 75 cysteine | the merged rule flags **511**, disagreeing on **145** (142 merge-only, 3 copy-only), with **0** cysteine-only calls; the docstring's figures do not reproduce on this structure | `design_pipeline/adapters/_gates.py`, `liability_reasons`; measured by `examples/liability_union_check.py` |
| 4 | `new_motifs` compares occurrences by start position, so a severity change inside one character class is silent | `A:T31G`, `B:T78G`, `B:S85G` turn `NT`/`NS` into `NG` at an already-matched asparagine and return `[]` | documented as a known blind spot, not fixed; 3 of 4,313 singles, none of them in the 432 | `_gates.py`, `new_motifs` docstring |
| 5 | nine Env B routes were probed against the Env A interpreter | preflight reported `import boltz failed` while `boltz` ran the entire final screen out of `.venv-boltz`; 9 routes affected | each route declares the `env_var` its adapter reads, guarded by `tests/test_route_drift.py` against a reconstruction of the pre-change table | `CHANGELOG.md`, 2026-11-06; `design_pipeline/routes.py` |

A sixth, with no number attached: `AntiBertyMaskedLM.log_probs` used a bare `import torch`
and raised `ModuleNotFoundError` instead of this package's `MissingDependencyError`
(`design_pipeline/adapters/_plm.py`, `_torch`). Counting it would make six.

---

## 15. Provenance of the FoldX numbers

| file | what produced it |
|---|---|
| `runs/track_a_foldx/foldx.csv` | `design_pipeline/adapters/_foldx_score.py`, 432 rows, `FOLDX_RUNS=3`. **Header on disk is the pre-migration one: `variant_id, ddg, ddg_fold, n_runs` - there is no `ddg_sd` column at all, not an empty one.** `read_scored()` rewrites it to `variant_id, ddg, ddg_sd, ddg_folding, n_runs` on the next run and leaves `ddg_sd` empty for all 432 migrated rows, because a past run's spread cannot be recovered. **Every real SD in this ledger comes from `foldx_sd_sample.csv` (n=3, 60 candidates) or `foldx_runs20_top50.csv` (n=20, top 50)**, never from `foldx.csv`. |
| `runs/track_a_foldx/candidates_all432.csv` | `examples/foldx_precision/expand_432.py` |
| `runs/rerank/ranking_qualified_v1_432.csv` | `examples/foldx_precision/rank_432.py` |
| `runs/track_a_foldx/foldx_sd_sample.csv` | `examples/foldx_precision/measure_sd.py`, `FOLDX_RUNS=3` |
| `runs/track_a_foldx/foldx_runs20_top50.csv` | `examples/foldx_precision/measure_runs20.py`, `FOLDX_RUNS=20` |
| tier, band and cluster figures | `notebooks/07_foldx_precision_and_rerank.ipynb` |
| the within-1-SD and within-2-SD counts | `examples/foldx_precision/sd_report.py` |
| the rank distribution of the 32 "ours" | `examples/foldx_precision/summary_432.py` |
| run logs | `examples/foldx_precision/logs/` |
| the whole run, end to end | `notebooks/07_foldx_precision_and_rerank.ipynb`, executed copy at `runs/07_foldx_precision_and_rerank.executed.ipynb` |

Every figure in sections 8 to 10 is re-derived from the output files by notebook 07, so
executing it checks this ledger. Its FoldX steps are guarded and skip when their output is
already on disk; `FORCE_FOLDX = True` overrides that and costs hours.

Structure for every FoldX number above: `data/structures/1N8Z_fv_trimmed_Repair.pdb`,
sha256 `79ce088822d539a99517770ccb2d4b089ad25778a5529fb6d0e37fe92ddbaba9`, identical to the
committed blob and to the structure the original 115 were scored against.
