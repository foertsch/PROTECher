# Candidate expression and stability risk screening

This is a local screening report for substitutions of the organisers' 242-aa
trastuzumab scFv. It identifies sequence-rule violations and descriptive risk
signals. It does **not** predict expression yield, melting temperature, folding
free energy or binding affinity. Absence of warnings is not experimental validation.

## Run

Run commands from the `PROTECher` repository root. Relative input and output paths
are resolved against that root, even when launching the script from elsewhere.

```bash
# Existing candidate set; no downloads, model inference or GPU.
.venv/bin/python hackday/05_developability.py

# Candidate-only FASTA or CSV; the official parent is inserted automatically.
.venv/bin/python hackday/05_developability.py --input data/designs/new_candidates.csv

# Sequence checks without reading structure or AntiFold files.
.venv/bin/python hackday/05_developability.py --sequence-only

# Directly consume track_a_export.py's default submission CSV (no header).
.venv/bin/python hackday/05_developability.py \
  --input submission/submission.csv --csv-no-header --sequence-only

# Attach externally calculated ipSAE scores, without ranking or gating on them.
.venv/bin/python hackday/05_developability.py \
  --input data/designs/new_candidates.csv \
  --ipsae data/designs/ipsae.csv \
  --out data/designs/new_candidates_review.csv
```

The legacy `--fasta` flag is an alias of `--input`. Input format is CSV for a
`.csv` suffix and FASTA otherwise. CSV headers must include `name,sequence`.
For the team's headerless submission CSV, explicitly add `--csv-no-header`;
that mode requires exactly two columns and rejects a header. Header detection
is deliberately not guessed, so the first candidate cannot silently disappear.
The exporter's `scfv_for_folding.fasta` can also be read directly.
FASTA headers are used in full as identifiers. Whitespace is removed from
sequences and lowercase amino acids are normalized to uppercase. IDs must be
unique and nonempty; preserve the same IDs across the team's tools.

`PARENT`, `parent`, `trastuzumab_scFv_parent` and `trastuzumab_scfv_parent` are reference
aliases. If present, they must exactly match the official sequence. `--parent`
changes the reference row's label, not the reference sequence. An unmodified
sequence under any other name is treated as a candidate with zero substitutions.

The optional `developability` extra in `pyproject.toml` declares Biopython:

```bash
python -m pip install -e '.[developability]'
```

Use the project's environment for installation. This workspace's `.venv` already
has Biopython 1.88 installed. If Biopython is absent on another machine, screening
continues with the instability metric explicitly unavailable. NumPy is required;
Biotite is needed only for structural patches. Neither missing optional data nor
missing Biopython is replaced by a numeric zero.

## Rules and reference mapping

The fixed reference in the script is the exact sequence supplied by the organisers.
All reported positions are 1-based positions in the full scFv:

| Region | Positions |
|---|---|
| VH | 1–120 |
| Linker, `GGGGSGGGGSGGGGS` | 121–135 |
| VL | 136–242 |
| Protected cysteines | 22, 96, 158, 223 |

Candidates must have at least five substitutions and preserve the linker and
protected cysteines. Newly introduced cysteines are reported separately as risks,
including when total cysteine count has stayed unchanged. No maximum number of
submissions or mutations is invented by this script.

Positional comparison currently supports equal-length, substitution-only designs.
Length-changing designs retain their sequence-level metrics but get
`review_required`; mutation count, linker/cysteine mapping, new-motif comparison,
structure scores and AntiFold comparison are not inferred using a truncated zip.
They need alignment and clarification of the applicable design rules.

Malformed files, duplicate IDs and an incorrect reference stop the run before
writing results. Invalid individual sequences are retained as `invalid_input`
rows with numeric scores left blank, so other candidates can still be reviewed.

## Metrics and limitations

- **pI and charge at pH 7.4:** fixed-pKa sequence estimates, retaining the existing
  script's method. A pI change of at least 1.0 produces a heuristic warning.
  This is not an expression-system-specific solubility model.
- **GRAVY and hydrophobic runs:** sequence descriptors, with differences from
  the parent. A longer run of `AVILMFWC` residues generates a warning.
- **Instability index:** the complete Guruprasad calculation provided by
  [Biopython ProteinAnalysis](https://biopython.org/docs/latest/api/Bio.SeqUtils.ProtParam.html).
  The previous incomplete weight table has been removed. This descriptive index
  is not a validated scFv rejection criterion; no cutoff is used.
- **New liability motifs:** NG/NS/NT, DG/DS/DD, DP and N-X-S/T with X other than P.
  Potential glycosylation now produces both `new_liability` and
  `new_glycosylation_motif`. A motif is a sequence warning, not proof of a chemical
  event; expression host and structural context affect interpretation.
- **PSH/PPC/PNC:** retained column names for custom hydrophobic/positive/negative
  patch heuristics. These are **not the published TAP metrics**. `--tol` (default
  1.0) is an absolute increment in each heuristic's own units, not a standard
  deviation, calibrated clinical threshold or success probability.
- **SFvCSP:** the existing product of estimated VH and VL charges, with each domain
  treated separately. This is a descriptor, not the published structural TAP
  charge-symmetry implementation or a model of the connected scFv's termini.
- **Structure:** default `data/structures/1N8Z.pdb`, chain B residues 1–120 and
  chain A residues 1–107. Coverage, amino acid identities, insertion codes and
  finite coordinates/exposures are checked before scoring. This path expects that
  exact numbering, not an arbitrary Boltz output. All patch metrics reuse parent
  coordinates and exposure, without rebuilding mutant side chains or modelling
  the linker. Mutations at relative exposure below 0.25 get `approx_buried`.
- **AntiFold:** default `data/scores/antifold/1N8Z_BAC.csv`. Cached reference residue
  identities and unique mapped positions are verified. `antifold_d` sums mutant
  minus parent log probabilities at changed positions. Higher values mean greater
  preference under this additive approximation; they are not energies. This does
  not model interactions among mutations or perform a new full-variant inference.
  If even one changed position is unscored, the whole result is unavailable.

There is no combined score, automatic candidate rejection, calibrated expression
predictor or ranking by ipSAE. All candidates remain in the report.

## Outputs

Default outputs are `data/designs/developability_review.csv` and
`data/designs/developability_review.metadata.json`. The previous
`data/designs/developability.csv` is preserved. Rerunning overwrites the selected
output pair; choose `--out` for separate snapshots. Input paths cannot also be
used as output paths.

Each row contains the original normalized sequence, a sequence SHA-256, mutation
list, rule checks, metric values, parent-relative deltas, flags and availability
statuses. Empty numeric cells mean **unavailable**, not zero. All reference deltas
are zero when the corresponding metric is available.

`review_notes` explains the warnings in plain English. `liability_details` names
the substitutions overlapping each newly introduced motif, using full-scFv
numbering. For example, K30S creates `glyc:N28`, and K30D creates `isomer:DD@30`.
These identify the local sequence change, not a confirmed chemical reaction.

Identical candidate sequences under different IDs are retained and flagged as
`duplicate_sequence`; `duplicate_sequence_ids` lists the matching IDs. This helps
coordinate redundant GPU runs without removing intentional controls or seed
replicates. It is not a biological failure or an additional competition rule.
The automatically inserted reference is excluded from candidate duplicate checks.

| Column | Interpretation |
|---|---|
| `input_status` | Valid sequence or invalid amino acids/empty sequence |
| `rules_status` | `reference`, `pass`, `fail` or `review_required` |
| `rule_issues` | Specific sequence-rule violations or unsupported length change |
| `flags` | Potential biological risks and buried-mutation approximation warning |
| `review_notes` | Plain-English explanation of rules, warnings and missing evidence |
| `liability_details` | New motifs and the substitutions at their positions |
| `duplicate_sequence_ids` | Other candidate IDs with the same normalized sequence |
| `data_warnings` | Missing/unusable resources or skipped positional calculations |
| `*_status` | Whether the corresponding metric was available and its method |
| `assessment` | `reference`, `invalid_input`, `rule_violation`, `needs_review` or `no_warnings_detected` |

`no_warnings_detected` only means this limited set of checks raised no warning.
Missing ipSAE is recorded separately and does not change the developability
assessment. Negative AntiFold deltas are reported without an arbitrary cutoff.

The metadata records input/script hashes, reference sequence, resource paths,
package versions, parameters, timestamp and methodological limitations.

## Integrating teammates' results

Provide new candidates as FASTA or `name,sequence` CSV. They can omit the parent.

Provide ipSAE as `design,sequence,ipsae` CSV, with one externally selected or
aggregated value per design. The sequence is required to catch stale scores after
a design has changed. Duplicate/unknown IDs, mismatched sequences and non-finite
or out-of-range scores stop the import. Partial coverage is allowed; unscored
designs have blank ipSAE cells. Use the report's reference label if supplying a
parent score. This script does not compute ipSAE or aggregate seeds.

### Needed now, before GPU runs

The only required new input for initial screening is the candidate sequence list
with stable, unique IDs. Provisional candidates are enough. The official parent
and sequence constraints are already configured; no GPU, ipSAE or new structures
are required. Please confirm whether the 15 local variants are the current team
set, or provide the replacement list. Give untagged VH-linker-VL scFv sequences;
describe tags/signal peptides separately rather than silently adding them.

Helpful context, if known (not a blocker): which candidates are intended for the
next GPU batch, the rationale for the mutations, and any intentional duplicate
controls. Expression host/system, tags, secretion format, buffer/pH/temperature
and the relevant assay endpoint help interpret warnings; unknown is acceptable.
Existing rationale may already be in the team's research notes.

### Later, after cofolding

- Send ipSAE scores keyed to the same IDs **and exact sequences**, plus the
  metric/aggregation convention, seeds, Boltz-2 settings and HER2 domain-IV
  sequence/boundaries used. Include the parent control if it was run. These
  results can be attached through `--ipsae`; they are not required for screening.
- Provide candidate structures and chain/residue mappings if candidate-specific
  structural checks are wanted. That adapter is not implemented yet.
- If whole-candidate AntiFold scores already exist, retain their method and
  antigen-inclusion metadata. They need a separate importer and must not replace
  the additive proxy column silently or be mixed with it in one ranking. The
  new ADAPT notebooks do not establish reliable affinity ranking for our designs.
- Experimental expression/stability results can support later calibration;
  they are not a prerequisite for the initial report.

## Verification and current local results

```bash
.venv/bin/python -m pytest -q tests/test_hackday_developability.py
```

Checks cover parent identity, zero deltas, the published Biopython numerical
example, rule boundaries, cysteine relocation, glycosylation, invalid inputs,
length changes, missing/partial scores, reference mismatch, score joins and CLI
outputs. They require no network or model inference.
Integration tests feed the actual team exporter outputs (FASTA, CSV with a header
and CSV without a header) into screening and check numbering across VH and VL.
Regression checks also cover duplicate sequences and mutation-level explanations.

The initial revised run on the 15 existing variants found all 15 compliant with
the checked sequence rules; 11 had warnings and 4 had none under these heuristics.
These are not 11 failures and 4 experimentally validated successes. No ipSAE
scores were supplied. The corrected parent instability index is 43.41 (the old
partial-table value was 21.2); old and new instability columns are not comparable.
