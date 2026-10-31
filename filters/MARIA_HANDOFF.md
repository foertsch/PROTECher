# Maria: model comparison, integration and filtering

**Working branch:** `maria/results-integration-handoff`
**Based on:** your existing `maria/score-comparison` work (`0d229c3`), with current
`main` (`d9e56d4`) merged locally. Neither existing branch was changed remotely.

## What your part does

Bring model outputs and Iason's developability report together for the **same
candidate sequences**, make differences and missing evidence visible, and apply
the filtering policy agreed with Alixe and the model owners. This tool does not
generate antibody designs or run models. It supports the final selection; it does
not claim that model agreement proves binding or successful expression.

There are two distinct stages:

| Script | Input and purpose |
|---|---|
| `filters/inspect_scores.py` (your existing work) | An already combined mutation-level pipeline CSV/Parquet; one row per variant, model ranks, same-site ranks for singles, correlations. It assumes normalized higher-is-better `score_*` columns. |
| `filters/filter_variants.py` (your existing work) | Stage-3 consensus **for single substitutions only**, using the existing TOML configuration. It is not the final 5+-mutation design filter. |
| `filters/compare_candidates.py` (new) | The authoritative whole-candidate sequences, separate model CSVs and an optional Iason report; identity-checked integration, seed aggregation, comparison, explicit filters and an optional diverse shortlist. |

The existing scripts and TOML file have been preserved. Their model list and
calibration assumptions have not been approved by this implementation. In the
latest `CHANGELOG.md`, default scorers changed to `ablang2` and `consensus`;
AntiFold remains available but is no longer a default ranking tool. Do **not**
apply the historical six-sequence/four-inverse-model consensus configuration to
whatever columns happen to arrive, or transfer a small-study AUROC to a different
format without review. Confirm the actual model set and policy with the team.

## Completed and usable now

- Join independent model outputs to a candidate manifest without dropping
  unscored candidates or multiplying rows accidentally.
- Verify exact sequences or sequence hashes. Unknown IDs, stale sequences,
  duplicate seeds and mixed declared protocols cause clear errors.
- Support a frozen run identity map for model files that only carry mutation IDs
  or model-specific IDs.
- Aggregate repeated predictions **within one source/protocol** using median,
  count, minimum, maximum and sample standard deviation. Missing seeds remain
  visible; insufficient replicate coverage is not treated as a complete score.
- Preserve each score's direction and source namespace. Report per-metric ranks,
  within-mutation-count ranks, and deltas from an explicit parent control.
- Compare models using direction-adjusted Spearman correlations with overlap
  counts, both overall and within mutation-count groups. This is descriptive
  agreement, not validation against experimental ground truth.
- Import Iason's report: sequence-rule failures can reject a design, biological
  warnings request review, and absent required evidence means INCOMPLETE.
- Apply explicitly configured numerical gates. There are no invented production
  cutoffs or averaged scores across models.
- Optionally select a deterministic shortlist using agreed score priorities and
  minimum sequence distance. Exact duplicates cannot both enter the shortlist.
  The algorithm reports a shortfall rather than silently relaxing the constraints.
- Produce CSVs, a readable Markdown summary and JSON provenance (input hashes,
  actual configuration, source protocols and script hash).
- Provide a runnable **synthetic** demo and automated regression tests. The demo
  is not a candidate proposal, experiment or calibration result.

## Get the branch

From your own clean checkout (preserve any uncommitted work first):

```bash
git fetch origin
git switch --track origin/maria/results-integration-handoff
```

If that local branch already exists, use `git switch maria/results-integration-handoff`
and then `git pull --ff-only`. Do not reset your existing branch or overwrite
uncommitted work. This branch already contains your original reporting scripts.
No GPU packages are needed for comparison; the repository's usual pandas
dependency is sufficient. Tests use pytest.

## Run the demonstration now

Run from the repository root. Replace `.venv/bin/python` with your Python
environment's interpreter if necessary.

```bash
.venv/bin/python filters/compare_candidates.py \
  --candidates filters/demo/candidates.csv \
  --config filters/demo/comparison.json \
  --developability filters/demo/developability.csv \
  --out runs/maria_demo

.venv/bin/python -m pytest -q tests/test_candidate_comparison.py
.venv/bin/python -m pytest -q tests/test_maria_legacy_filters.py
```

Expected demo: one ELIGIBLE, one REVIEW, one INCOMPLETE and one REJECT, plus the
parent. The requested shortlist of two contains only one because the defaults do
not promote a warning-bearing or incomplete candidate automatically. Demo
thresholds are arbitrary illustrative values, **not recommendations**.

Verification completed: **139 tests passed** across the new integration, the
preserved single-substitution scripts, existing combination/reporting code, scFv
mapping and export tests. The older `run_report` plotting tests need matplotlib;
it was installed in the local environment for that check. The new comparison
script itself only needs pandas and the standard library.

An additional local-format check joined Iason's existing report successfully:
15 provisional candidates remained INCOMPLETE because the configured model result
files were absent. No final candidate ranking was inferred from that check.

## Use real candidates when available

```bash
.venv/bin/python filters/compare_candidates.py \
  --candidates /path/to/approved_candidates.fasta \
  --config filters/candidate_comparison.json \
  --developability /path/to/developability_review.csv \
  --out runs/maria_candidate_review
```

The supplied configuration is an editable starting point, not a committed model
selection. Missing source files produce an INCOMPLETE report, so preparation can
start before model runs. Once actual files arrive, replace the `SET_...` protocol
placeholders, point at the real files, and agree which metrics are required.

Paths inside the JSON are relative to that JSON's directory. CLI paths are
relative to your current working directory. Output names are fixed within `--out`;
use a new output directory for each candidate version. Inputs cannot be overwritten.

### Candidate identity

The manifest is FASTA or CSV with `name,sequence`. Include the exact parent and
set `parent_id` to its ID (or omit it when no parent is available). No parent is
silently inferred. IDs must be unique; sequences are normalized to uppercase with
whitespace removed. The entire FASTA header is the ID. Keep the same convention
across tools. Manifest validation checks syntax, not scFv design rules; those
come from Iason's report. Duplicate sequences under different IDs are reported.

The 15 sequences previously found under Iason's local `data/designs` are
**provisional local inputs, not the team's agreed shortlist**. They are not
included in this branch. Obtain confirmation or the replacement set before
presenting biological conclusions.

### Model input formats

Each JSON source has a unique `name`, CSV `path`, `protocol` description,
`id_column` (default `design`) and a list of metrics. Each metric declares
`name`, source `column`, `better` (`higher`, `lower`, or `descriptive`), and optional
`required` (default false) and `units`. Output keys are `source.metric`.

Choose **exactly one** identity method per source:

1. `sequence_column`: each prediction row carries the exact candidate sequence.
2. `hash_column`: each row carries SHA-256 of the normalized candidate sequence.
3. `identity_file`: a CSV with `source_id,design,sequence`, saved **with the model
   run**. This maps raw `variant_id` or cofolding complex IDs back to manifest IDs.
   Do not manufacture this mapping from the latest candidate file for an old run;
   that would conceal stale inputs. The tool verifies declared identities, not the
   authenticity of the producer's history.

For repeated predictions, set `replicate_column` (e.g. `seed`) and
`expected_replicates`. If the model emits several predictions per seed, give each
observation a unique compound ID such as `seed0:model1`, and define the expected
count accordingly; do not duplicate the same candidate/replicate key. Scores
from distinct checkpoints, antigen constructs or scoring methods must be separate
sources. Optional `protocol_column` requires every row to match the source's
declared protocol. Blank score cells mean unavailable; malformed/non-finite values
are rejected rather than silently interpreted as successful observations.

No raw-model format is guessed. `score_variants` CSV outputs can use a frozen
identity map and `id_column: "variant_id"`. A Boltz summary with `complex,seed,...`
can use `id_column: "complex"` and its own map. ipTM is **not** ipSAE: only declare
an ipSAE metric when that score has actually been calculated. Future raw JSON or
Parquet adapters require representative files and a tested mapping.

### Iason's report

Pass the CSV to `--developability`. Required columns are `design,sequence,
input_status,rules_status,flags,data_warnings`; optional `rule_issues` and
`review_notes` are carried into the output. Unknown or duplicate IDs and sequence
mismatches fail the join. The parent row can be aliased to `parent_id` only when
explicitly marked as a reference and its exact sequence matches. Include the
parent in your manifest when using Iason's default report.

This importer does not recalculate his metrics or promote his additive AntiFold
proxy to a whole-candidate model score. His package does not need to be installed
on your branch: the CSV is the interface.

### Filtering and selection

Decision precedence: REFERENCE; known rule/gate failure → REJECT; missing required
evidence → INCOMPLETE; warnings → REVIEW; otherwise ELIGIBLE. Every reason is
retained even when another reason takes precedence. ELIGIBLE only means the
configured checks did not block selection, not experimentally validated quality.

`gates` defaults to `[]`. An enabled gate explicitly gives `metric`, `operator`
(`>=` or `<=`), `threshold` and `action` (`review` or `reject`). An absent or
under-replicated gate score means INCOMPLETE, never zero or pass. Missing optional
metrics remain visible without blocking other work.

`shortlist` defaults to `null`. After policy agreement, set `size`, ordered
`sort_by` metric keys, `min_sequence_distance` (>=1) and optionally
`include_review` (default false). This is lexicographic priority, not an ensemble
weighted average. ipSAE is permitted as a gate/tie-breaker, not the first sorting
metric. Sorting does not prove affinity. The greedy diversity rule uses Hamming
distance and does not mix unequal sequence lengths in one shortlist; indel-aware
selection and campaign-specific arm quotas are not implemented.

## Deliverables

| Output | Purpose |
|---|---|
| `combined.csv` | One row per manifest candidate; all scores, availability, ranks, deltas, warnings and decisions |
| `observations.csv` | Individual seed/run metric values with verified sequence hash and protocol |
| `eligible.csv` | ELIGIBLE candidates only; not an automatically approved submission |
| `shortlist.csv` | Explicitly requested selection, otherwise empty with headers |
| `model_agreement.csv` | Direction-adjusted correlations, overlap counts and mutation-count strata |
| `summary.json` | Configuration, source status and reproducibility hashes |
| `report.md` | Readable counts, side-by-side comparison and limitations for the demo |

## What remains, who supplies it, and when

| Dependency | Owner/input needed | When and next action |
|---|---|---|
| Candidate identity/version | Design owners: final or provisional manifest with stable IDs; flag intentional controls | As soon as a version is ready: validate it and freeze run identity maps before scoring |
| Model set and score semantics | Arion/model owners: model/checkpoint, exact score definition/direction, full-candidate vs single-mutation method, output columns and run input mapping | Before configuring real sources; latest defaults are a starting point, not team approval |
| Priority and cutoffs | Alixe + team: required metrics, gate actions/thresholds, seed policy, shortlist size and diversity requirements | Before declaring a final shortlist; until then compare and report missing evidence |
| Sequence/developability results | Iason: report for the exact manifest version | Can arrive before GPU; rerun the join whenever sequences change |
| Model outputs | Model owners: CSV/Parquet/JSON examples, with exact input identity and seed/run IDs | As each run finishes: validate/adapt format and rerun the report; do not wait for every model to start integration |
| Boltz/ipSAE | Cofolding owner: structures, ipSAE calculation definition/cutoffs, HER2 domain-IV boundaries, seeds and parent controls | After GPU runs and ipSAE calculation: attach scores with matching identities/protocols |
| Calibration evidence | Kira + model owners: comparable labelled measurements and assay metadata | Before claiming validated thresholds/weights; otherwise keep decisions explicitly heuristic |
| Final choice/demo | Alixe + design owners + Maria: inspect REJECT/REVIEW/INCOMPLETE reasons and confirm selection | When enough evidence is available; demo can be rehearsed now on labelled synthetic inputs |

No dates or completion estimates are invented for other people's runs. The
integration code, demo and tests are ready now; final biological ranking waits
on the inputs and decisions above. No models or GPU jobs have been started, and
no sequences have been submitted.

## Prompt for Maria's ChatGPT/Codex

> Continue the PROTECher results-integration task on branch
> `maria/results-integration-handoff`. Read `filters/MARIA_HANDOFF.md` first.
> Preserve my original `inspect_scores.py` / single-substitution filtering work;
> use `compare_candidates.py` for complete candidate designs. Run the synthetic
> demo and regression tests before adapting inputs. Inspect the real candidate
> manifest, model outputs and run identity maps when they arrive; do not invent
> missing scores or assume IDs alone prove sequence identity. Confirm model
> directions/protocols with the model owners and filtering priorities with Alixe.
> Keep ipSAE as structural evidence, distinguish whole-candidate AntiFold from
> additive proxies, and retain incomplete candidates with reasons. Generate the
> combined report and only produce a final shortlist once the policy is agreed.
> Do not change main or other people's branches. Keep work on my branch and ask
> me before any merge, submission or new GPU job.
