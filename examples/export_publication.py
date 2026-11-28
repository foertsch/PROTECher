"""The publication CSVs: the measured benchmark, its labels, and what qualified.

BENCHMARK only. Molecules other people measured, the labels this project put on them,
and the per-column qualification table. **No CONTROLS row and no DESIGNS row is
exported** - the manufactured negatives encode our own assumptions and the 432 designs
have no measured affinity, so neither belongs in a table other people will read as
data. The design-population figures live in `docs/article_numbers.md`.

Three files, written to `data/publication/` with a README naming each block's estimand:

    benchmark.csv      one row per measurement, every study in the evidence file
    negatives.csv      one row per substitution, the labelled direction and Chinery sets
    qualification.csv  one row per score column, with its panel state and its ICC

What this does to the source tables, and nothing else:

* renames columns whose names do not say what they hold (`ddg_reported_kcal_mol` ->
  `ddg_kcal_mol`, `kon_Minv_sinv` -> `kon_per_M_per_s`, `raw_KD_M` -> `kd_molar`);
* adds `numbering_scheme`, a controlled vocabulary over the eight free-text values in
  `numbering_source`, which is kept beside it as `numbering_source_text`. A position
  without its scheme is not a position;
* adds `block`, so the three estimands in the README can be selected for;
* adds `in_direction_panel`, from `runs/negative_controls/direction_set.csv`, so the
  29-variant panel is reconstructable from the export alone;
* relabels the qualification table's panel verdicts into the four states used
  everywhere else (PASS / FAIL / UNPROVEN / N/T) and fills in the ipSAE row, which
  said "no panel test" while its panel sat in `runs/negative_controls/ipsae/`.

**`I.C.` and `N/A` are left exactly as they are.** They are the instrument's own words:
`I.C.` is an inconclusive fit and `N/A` is an empty well, and the source does not say
which of "did not bind", "did not express" or "sensor did not regenerate" produced
either. Turning them into blanks or zeros would invent a measurement. Every read and
every write here uses `keep_default_na=False` for that reason, so a blank in the output
is a blank in the source and nothing else.
"""

from __future__ import annotations

import argparse
from pathlib import Path

import pandas as pd

EVIDENCE = Path("data/benchmarks/trastuzumab_mutation_evidence.csv")
NEGATIVES = Path("data/benchmarks/trastuzumab_negatives.csv")
QUALIFICATION = Path("data/track_a_final/qualification.csv")
DIRECTION_SET = Path("runs/negative_controls/direction_set.csv")

# The three blocks with an estimand of their own, named in the README. Everything else
# in the evidence file is single-study literature context and is marked as such rather
# than dropped: eight of the direction panel's negatives come from it.
BLOCKS = {"adapt2017": "adapt", "moon2016": "moon", "chinery2026": "chinery"}
OTHER_BLOCK = "other_literature"

# `numbering_source` is free text with eight values. A position means nothing without
# the scheme it is an index into, so the scheme becomes a column of its own with a
# vocabulary an unfamiliar reader can act on. An unmapped value is an error, not a
# fallback: silently calling an unknown scheme "not_stated" is how an off-by-one ships.
NUMBERING_SCHEME = {
    "1N8Z PDB numbering": "1n8z_pdb",
    "Kabat": "kabat",
    "Kabat; incomplete clone sequences": "kabat",
    "mature-domain sequential; derived from supplied VH/VL": "mature_domain_sequential",
    "source labels; consistent with mature-domain sequential positions": "mature_domain_sequential",
    "full-sequence loop replacement; no unique deletion coordinate assigned": "not_assigned",
    "Fc source numbering; hinge mutation contextual": "fc_source",
    "": "not_stated",
}

EVIDENCE_RENAMES = {
    "source_index": "source_row_index",
    "data_source_index": "source_table_index",
    "numbering_source": "numbering_source_text",
    "ddg_reported_kcal_mol": "ddg_kcal_mol",
    "affinity_improvement_calculated": "affinity_fold_calculated",
    "affinity_improvement_reported": "affinity_fold_reported",
    "raw_KD_M": "kd_molar",
    "kon_Minv_sinv": "kon_per_M_per_s",
    "koff_sinv": "koff_per_s",
    "source_capture_RU": "source_capture_ru",
    "target_or_cell": "target_or_cell_line",
    "full_VH_sequence": "vh_sequence",
    "full_VL_sequence": "vl_sequence",
    "source_CDRH3_design_sequence": "cdrh3_design_sequence",
}

NEGATIVES_RENAMES = {
    "source_file": "source_block",
    "wt": "wt_residue",
    "mut": "mutant_residue",
}

QUALIFICATION_RENAMES = {
    "column": "score_column",
    "qualified": "used_in_ranking",
    "tested_on_negatives": "tested_on_negative_panel",
    "panel_verdict": "panel_state",
    "icc": "icc_all_designs",
    "icc_varying": "icc_designs_that_vary",
    "icc_used": "icc_reported",
    "n_designs": "n_designs_for_icc",
    "seeds": "n_seeds_for_icc",
    "why": "reason",
}

# The four states, applied to the qualification table so it agrees with
# `runs/negative_controls/verdict_states.csv` instead of carrying its own words.
PANEL_STATE = {"pass": "PASS", "fail": "FAIL", "unproven": "UNPROVEN", "": "N/T"}

# ipSAE's panel test exists and the qualification table said it did not. The row is
# filled in from the file that holds it; the exclusion still rests on the ICC.
IPSAE_PANEL = {
    "panel_state": "UNPROVEN",
    "panel_detail": "direction AUROC 0.476 (CI 0.226-0.732), n = 29; interval spans 0.5",
}


def read(path: Path) -> pd.DataFrame:
    """Everything as text, with `I.C.` and `N/A` left as the strings they are."""
    return pd.read_csv(path, dtype=str, keep_default_na=False, low_memory=False)


def scheme_of(value: str) -> str:
    if value not in NUMBERING_SCHEME:
        raise KeyError(
            f"numbering_source {value!r} has no entry in NUMBERING_SCHEME. Add one "
            f"rather than letting a position ship without the scheme it indexes into."
        )
    return NUMBERING_SCHEME[value]


def benchmark(evidence: Path, direction_set: Path) -> pd.DataFrame:
    frame = read(evidence)
    panel = set(read(direction_set)["evidence_id"]) if direction_set.is_file() else set()
    frame["numbering_scheme"] = frame["numbering_source"].map(scheme_of)
    frame["block"] = frame["study"].map(BLOCKS).fillna(OTHER_BLOCK)
    frame["in_direction_panel"] = frame["evidence_id"].isin(panel).map(
        {True: "yes", False: "no"})
    frame = frame.rename(columns=EVIDENCE_RENAMES)
    lead = ["evidence_id", "block", "study", "year", "variant", "format",
            "numbering_scheme", "numbering_source_text", "mapping_basis",
            "substitutions_mature_domain_1based", "metric", "value", "units",
            "parent_value", "ddg_kcal_mol", "ddg_sd_kcal_mol", "affinity_fold_reported",
            "outcome", "in_direction_panel"]
    rest = [c for c in frame.columns if c not in lead]
    return frame[lead + rest]


def negatives(path: Path, evidence: Path) -> pd.DataFrame:
    frame = read(path)
    # The scheme a negatives row's position is in is the scheme of the evidence row it
    # came from, so it is joined rather than assumed. Two joins are needed: the 31
    # literature rows carry an `evidence_id` and the 4,520 Chinery rows carry none,
    # only a `variant`, which resolves 687 of 687 against the evidence file and gives
    # one scheme per variant. Falling back to "not_stated" would have labelled every
    # Chinery position as scheme-less when the scheme is known.
    source = read(evidence)[["evidence_id", "variant", "numbering_source"]]
    by_id = (source[source["evidence_id"] != ""]
             .drop_duplicates("evidence_id").set_index("evidence_id")["numbering_source"])
    by_variant = (source[source["variant"] != ""]
                  .drop_duplicates("variant").set_index("variant")["numbering_source"])
    frame["numbering_source"] = (frame["evidence_id"].map(by_id)
                                 .fillna(frame["variant"].map(by_variant))
                                 .fillna(""))
    frame["numbering_scheme"] = frame["numbering_source"].map(scheme_of)
    blank = frame[(frame["numbering_scheme"] == "not_stated")
                  & frame["study"].isin(BLOCKS)]
    if len(blank):
        raise ValueError(
            f"{len(blank)} negatives row(s) from a named block resolved to no "
            f"numbering scheme: {sorted(set(blank['variant']))[:5]}. A position without "
            f"its scheme should not ship."
        )
    frame["block"] = frame["study"].map(BLOCKS).fillna(OTHER_BLOCK)
    frame = frame.rename(columns={**NEGATIVES_RENAMES,
                                  "numbering_source": "numbering_source_text"})
    lead = ["variant", "evidence_id", "block", "study", "label_class", "label_reason",
            "n_mutations", "chain", "position", "wt_residue", "mutant_residue",
            "numbering_scheme", "numbering_source_text", "format", "source_block"]
    rest = [c for c in frame.columns if c not in lead]
    return frame[lead + rest]


def qualification(path: Path) -> pd.DataFrame:
    frame = read(path)
    frame["panel_verdict"] = frame["panel_verdict"].str.strip().str.lower().map(PANEL_STATE)
    if frame["panel_verdict"].isna().any():
        unknown = frame.loc[frame["panel_verdict"].isna(), "column"].tolist()
        raise ValueError(f"panel_verdict not in the four states for: {unknown}")
    for key, value in IPSAE_PANEL.items():
        column = {"panel_state": "panel_verdict", "panel_detail": "panel_detail"}[key]
        frame.loc[frame["column"] == "ipsae_max", column] = value
    return frame.rename(columns=QUALIFICATION_RENAMES)


README = """# Publication CSVs - BENCHMARK only

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
"""


def main(argv=None) -> int:
    parser = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    parser.add_argument("--out-dir", type=Path, default=Path("data/publication"))
    parser.add_argument("--evidence", type=Path, default=EVIDENCE)
    parser.add_argument("--negatives", type=Path, default=NEGATIVES)
    parser.add_argument("--qualification", type=Path, default=QUALIFICATION)
    parser.add_argument("--direction-set", type=Path, default=DIRECTION_SET)
    args = parser.parse_args(argv)

    args.out_dir.mkdir(parents=True, exist_ok=True)
    tables = {
        "benchmark.csv": benchmark(args.evidence, args.direction_set),
        "negatives.csv": negatives(args.negatives, args.evidence),
        "qualification.csv": qualification(args.qualification),
    }
    for name, frame in tables.items():
        path = args.out_dir / name
        frame.to_csv(path, index=False)
        print(f"{name:<20}{len(frame):>6} rows{len(frame.columns):>5} cols  -> {path}")

    readme = args.out_dir / "README.md"
    readme.write_text(README, encoding="utf-8", newline="\n")
    print(f"{'README.md':<20}{'':>6}      {'':>5}      -> {readme}")

    bench = tables["benchmark.csv"]
    print("\nblocks:", bench["block"].value_counts().to_dict())
    print("direction panel rows:", int((bench["in_direction_panel"] == "yes").sum()))
    literal = [c for c in ("value", "kd_molar", "kon_per_M_per_s", "koff_per_s",
                           "source_fit_residual_sd") if c in bench.columns]
    for column in literal:
        counts = bench[column].value_counts()
        print(f"  {column:<24}I.C. {counts.get('I.C.', 0):>4}   "
              f"N/A {counts.get('N/A', 0):>4}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
