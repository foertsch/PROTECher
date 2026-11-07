# Changelog

What changed in `design_pipeline` and why, newest first. A change that moves a
default, a scale or a numbering belongs here with the number that moved it; a
refactor that changes no result does not.

## 2026-11-06

### `DEFAULT_SCORERS` is now `("ablang2",)`

`consensus` is out of the default set. Not because it ranked badly - it was
never run at all - but because it cannot run: it needs an OAS frequency table
that `DESIGN_PIPELINE_OAS_TABLE` points at, and that table does not exist and
never did.

Three independent checks, none of which found it:

| Where it would be | Result |
|---|---|
| the repository's history | `git log --all --diff-filter=A --name-only` matches no `oas`/`frequen` path on any branch |
| the GPU node inventory | absent from both tables in `docs/cluster_gpu_node_rebuild.md`, and `DESIGN_PIPELINE_OAS_TABLE` is not among the variables it records as set |
| the results | no `consensus` column in `data/track_a_final/ranking.csv` |

So no number in this repository was ever produced with it, and a default that
raises `MissingDependencyError` on every fresh checkout is not a default. The
adapter, the route and the registry entry are unchanged: `models=["consensus"]`
runs exactly as before the moment a table exists, and
[docs/provisioning.md](docs/provisioning.md) records the four steps to build
one as an open task.

`ablang2` is the one default left. It installs into Env A, needs no GPU, and
needs no structure model, which keeps the default campaign that
`python -m design_pipeline.preflight` gates on inside one environment.

### Nine Env B routes were probed against the wrong interpreter

No result changes. Recorded because it means every preflight table printed
before today understated what was installed.

`chai`, `chai1`, `boltz`, `boltz2`, `protenix`, `colabfold`, `protein_mpnn`,
`ligand_mpnn` and `soluble_mpnn` declared `module` or `executable` only, so
`resolve()` probed them with `importlib.util.find_spec` or `shutil.which`
against Env A - asking this interpreter about packages that are deliberately
not in it. The reasoning against doing that was already written in `routes.py`
for `antifold` and had not been applied to the routes that came from
`pyproject.toml`.

The concrete effect: `boltz` ran the entire final screen out of `.venv-boltz`
and preflight in Env A reported `import boltz failed` the whole time. Each
route now declares the `env_var` its adapter reads, `tests/test_route_drift.py`
fails CI if a new one does not, and that test also runs its check against a
reconstruction of the pre-change table so the rule cannot pass vacuously.

`foldx`'s `needs_env` now lists all three variables the adapter reads rather
than one, so preflight reports `DESIGN_PIPELINE_FOLDX_PDB` - the one with no
sensible default - when it is unset.

### `DESIGN_PIPELINE_RUNTIME` chooses the container runtime

`detect_runtime()` took the first of apptainer, docker on PATH. That is right on
a cluster, where apptainer is the only one permitted, and wrong on a workstation
with both, where apptainer won silently. The variable now overrides it, and a
runtime it names that is not on PATH is an error rather than a fall back to the
other one.

## 2026-10-31

### `DEFAULT_SCORERS` was `("ablang2", "consensus")`

`antifold` is out of the default set. Its adapter, `antifold_seq`'s adapter and
both routes are unchanged: `models=["antifold"]` or `models=["antifold_seq"]`
runs either exactly as before, and `antifold` stays in `MATRIX_SCORERS` and in
the registry. Only what a `score_variants` call with no `models=` touches has
changed, which also takes Env B (and a GPU) off the default campaign that
`python -m design_pipeline.preflight` reports on.

The reason is that on the one measured ddG panel we have for trastuzumab - the
22 ADAPT Fab variants (Vivcharuk et al. 2017), 11 singles and 11 combinations
at four sites, SPR against one batch-matched parent - AntiFold does not rank:

| Spearman against ddG | Fab only | HER2 (chain C) in the input |
|---|---|---|
| `antifold` matrix lookup, 11 singles | +0.291 (R^2 0.066) | **+0.045** (R^2 0.002) |
| `antifold_seq`, 7 doubles, within a fixed count | +0.179 | **+0.357** |
| `antifold_seq`, 4 triples, within a fixed count | 0.000 | **0.000** |
| `antifold_seq`, 11 combinations pooled | -0.655 | -0.609 |
| mutation count alone, same 11 | +0.837 | +0.837 |

The right-hand column is the condition closest to a binding campaign, and it is
where the structure-based scorer has least to say: +0.045 over the singles, and
R^2 0.002 on the fit of ddG.

The pooled figure for the combinations is the one to distrust in both columns.
`antifold_seq` reaches -0.655 over all 11, but minus the mutation count reaches
+0.837 on those same variants, so the pooled number is a correlation with how
far a variant sits from the parent rather than with binding. Holding the count
fixed is what removes that, and it is what leaves the within-size rows above.

Two things this is *not*. It is not an index offset: both runs matched
AntiFold's `pdb_res` column to the chain at all 434 antibody positions, and
`python -m antifold_score.check_wt` re-runs that as a round trip the notebook
cannot make of itself - the matrix carries AntiFold's own wild-type column
rather than the PDB's - alongside a shift profile, the label mapping and the
column order. And it is not a verdict on AntiFold as a filter - the AbBiBench
calibration in `tests/outputs/README.md` (419 and 2,079 variants) is what
`docs/track_a_filtering.md` reads as "a filter and not a ranking", and that
conclusion stands. A default that costs a second environment and a GPU has to
earn more than +0.045.

Provenance: two runs, both in the repo, differing only in whether AntiFold saw
the antigen.

- **`notebooks/03_adapt_calibration.ipynb`** (fgcz-r-023, 2026-10-31, committed
  in `e14d035`) is Fab only: `run_antifold` passes `--heavy_chain` and
  `--light_chain` and nothing else, and `antifold_seq`'s `ComplexScorer` loads
  the same two chains, so nothing in it saw HER2. AntiFold table `1N8Z_BA.csv`,
  434 rows.
- **`notebooks/03_adapt_calibration_antigen.ipynb`** (PR #17, `6af3749`) is the
  same panel, the same numbering gate and the same two scorers with chain C in
  AntiFold's input. Table `1N8Z_BAC.csv`, 1,015 rows. Its Part 6 re-scores
  `antifold_seq` without the antigen inside that same run as a control and
  reproduces -0.655 exactly, so the two notebooks are comparable and the
  difference between the columns above is the antigen and nothing else.

Why the doubles move from +0.179 to +0.357 while the triples sit at 0.000 in
both: the antigen shifts `antifold_seq` by mean -0.921 nats (sd 1.150, largest
3.023) over the 22 variants, and the two score columns still agree at Spearman
+0.972 - so the ordering barely changes. On 7 points Spearman is
1 - (sum d^2)/56, which moves in steps of 0.018, and +0.179 (sum 46) to +0.357
(sum 36) is about one adjacent swap. **Neither doubles figure is a result**; at
n=7 and n=4 the within-size numbers are what the pooled -0.655 looks like once
the mutation-count confound (+0.837 on the same variants) is taken out, which
is to say: nothing survives. The decision rests on the singles at +0.045, not
on these two.
