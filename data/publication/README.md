# Publication CSVs - BENCHMARK only

Written by `examples/export_publication.py`. Three files, each a view of a source table
in this repository with clearer column names and nothing recomputed.

**No design sequence of ours is in these files, and no manufactured control is.** No
affinity of any design in this project was measured, so no design appears in a table
other people will read as data, and the design sequences are the competition entry. The
CONTROLS population (deliberately-bad variants we generated) is excluded for the same
reason in reverse: those labels encode our own chemistry assumptions, not a measurement.
Summary statistics computed *over* the designs are a different thing and are included -
see `qualification.csv` below.

## The three estimands, which are not interchangeable

A fold change from one block does not mean the same thing as a fold change from another.
Pooling them is the single easiest mistake to make with these files.

| block | estimand | format | comparator | n |
|---|---|---|---|---|
| `adapt` | **monovalent ddG of binding**, Fab, SPR against a batch-matched parent | Fab | Herceptin Fab, batch-matched | 22 variants carry `ddg_kcal_mol` (26 variants, 50 rows in all) |
| `moon` | **IgG KD**, SPR, against a matched wild-type parent measured at **0.48 nM** in the same study | IgG | hu4D5, measured alongside | 7 variants with KD, plus the WT reference (8 variants, 39 rows) |
| `chinery` | **IgG KD, anti-human-Fc captured, 1:1 fit**. **No matched wild-type** and the block is **outcome-selected**: 156 clean fits out of 720 submitted designs | IgG | trastuzumab, not co-measured per plate | 720 `variant_outcome` rows plus 40 `assay_control` (760) |
| `other_literature` | single-study context, one estimand each - read `assay`, `metric` and `comparator` per row before using one | mixed | mixed | 80 rows over 12 studies |

Consequences worth stating once:

* **`chinery` is a lower bound on any correlation.** Both the weakest and the
  unmeasurable candidates are gone, the dynamic range is truncated, and the truncation
  is size-dependent: a larger design is likelier to fail to fit at all, so the
  large-substitution bins are the survivors of a harsher filter.
* **`moon`'s fold changes are against a parent in the same experiment**, which `chinery`
  has no equivalent of. Only `moon` supports "x-fold better than wild type" directly.
* **`adapt` is the only ddG block.** Every kcal/mol in this export comes from it.
* **Two values that both read 0.48, which must never share a sentence.** Moon's wild-type parent
  KD is 0.48 nM (`TZM-014`), in nanomolar, and it is what every Moon fold change is
  measured against. D98W's NCI-N87 relative antiproliferative activity is 0.48
  (`TZM-023`), a unitless IC50 ratio on a cell line. The affinity claim is `TZM-019`:
  KD 0.11 nM against that parent, reported as 4.3x.

## `I.C.` and `N/A` are literal

In `value`, `kd_molar`, `kon_per_M_per_s`, `koff_per_s` and `source_fit_residual_sd`:

* `I.C.` - the instrument returned binding with an inconclusive 1:1 fit (102 rows each).
* `N/A` - the well returned no number at all (454 rows each).

Neither is a zero and neither is a missing value to impute. `N/A` is consistent with a
design that does not bind, with one that never expressed, and with a sensor that did not
regenerate; the source does not distinguish them. A blank cell in these files is a blank
cell in the source.

## `numbering_scheme`

Added here, because a position is meaningless without the scheme it indexes into. It is
a controlled vocabulary over the eight free-text values of the source's
`numbering_source`, which is kept beside it as `numbering_source_text`:

| `numbering_scheme` | means |
|---|---|
| `1n8z_pdb` | residue numbering of PDB 1N8Z |
| `kabat` | Kabat numbering as the source gave it |
| `mature_domain_sequential` | 1-based from the start of each mature variable domain, derived from the supplied VH/VL |
| `fc_source` | the source's own Fc numbering; a hinge mutation, contextual |
| `not_assigned` | a full-sequence loop replacement with no unique coordinate |
| `not_stated` | the source gave none |

`mapping_basis` says how a row was tied to 1N8Z, and `not mapped` (89 rows) means it was
not. Do not read a position off a `not mapped` row.

## The files

### `benchmark.csv`

One row per measurement, not per variant: a variant with KD, kon and koff is three rows.
Every study in `data/benchmarks/trastuzumab_mutation_evidence.csv`, with `block` added
so the three estimands above can be selected for, and `in_direction_panel` marking the
29 variants of the direction test (21 measured improvements, 8 measured losses).
`record_type` separates `variant_outcome` from `assay_control`;
`evidence_type` marks the one `computational_only` ADAPT row, which is not a measurement.

### `negatives.csv`

One row per substitution, so a 7-substitution design is 7 rows: 4,551 rows over **700
designs**. `label_class` counted by design, not by row:

| `label_class` | designs | what it is |
|---|---|---|
| `negative_caveated` | 434 | Chinery `N/A` wells. Ambiguous by construction - excluded from every AUROC in this project |
| `positive` | 156 | a clean fitted KD |
| `negative_clean` | 97 | bound but could not be fitted |
| `negative_affinity` | 8 | measured weaker than the parent, from the literature |
| `negative_developability` | 3 | a developability liability |
| `negative_function` | 2 | a functional loss |

`L:R66G` is the weak one among the 8: the source reports 4.8 +/- 0.7 nM against
2.0 +/- 0.4 nM **with no significance test** and describes the affinity as similar.

### `qualification.csv`

One row per candidate score column, with whether it entered the ranking and why not.
`panel_state` uses the four states of `runs/negative_controls/verdict_states.csv`:
**PASS** (95% interval entirely above 0.5), **FAIL** (entirely below), **UNPROVEN**
(spans it), **N/T** (not testable).

Two things to read carefully:

* **`ipsae_max` is UNPROVEN, not FAIL and not untested.** The source table said "no
  panel test"; the panel is in `runs/negative_controls/ipsae/direction_ipsae.csv` and
  gives AUROC 0.476, CI 0.226 to 0.732 on 21 against 8. The interval spans 0.5, so
  nothing was shown either way. **ipSAE's exclusion from the ranking rests on its ICC
  of 0.298, a precision argument, not on this AUROC.**
* **The ICC columns are summary statistics over the design population** (115 or 116
  designs, 3 Boltz seeds each), and they belong here. What is withheld from this export
  is the competition entry - the design sequences and their per-design scores - not
  results computed over them. A reliability coefficient is a result: it says a column is
  mostly seed noise, which is why the column is not in the ranking. `foldx` has no ICC
  by construction:
  one fixed repaired structure, no seed, and its run-to-run spread is averaged away by
  `FOLDX_RUNS` and was not retained at the time.

## Related files not exported here

* `trastuzumab_foldx_alanine_scan.csv` (repository root) is the per-position FoldX
  alanine scan of the scFv. **It uses the older column name `ddg_fold`** for the
  folding total-energy change; the scorer writes that quantity as `ddg_folding` since
  2026-11-24. Same quantity, older header, and it is not the binding ddG - that is
  `ddg_bind` in the same file.
* `runs/negative_controls/verdict_states.csv` is the scorer verdict table.
* `docs/article_numbers.md` is the claim-by-claim ledger, including every DESIGNS
  figure deliberately absent from this export.
