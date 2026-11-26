# FoldX precision on the design population

What these scripts measured, and why they exist: `_foldx_score.py` reduced each
candidate's `FOLDX_RUNS` models to their mean and deleted the per-run values, so nothing
on disk said how well-determined any ddG was. These scripts keep the per-run values and
report the spread. The wrapper now carries a `ddg_sd` column, so new runs no longer need
them for that purpose — they stay as the provenance of the numbers in
[docs/article_numbers.md](../../docs/article_numbers.md) sections 9 and 10.

Each imports `run_foldx` and `parse_interaction_energy` from
`design_pipeline/adapters/_foldx_score.py`, so the numbers come from exactly the calls that
produced `runs/track_a_foldx/foldx.csv`. Validation: the remeasured 3-run means matched the
stored `ddg` for all 60 sampled candidates, maximum difference 0.0000.

| script | produced | settings |
|---|---|---|
| `expand_432.py` | `runs/track_a_foldx/candidates_all432.csv` | expands the 432 variant ids to one row per mutation through `design_pipeline.variants.parse_id` |
| `rank_432.py` | `runs/rerank/ranking_qualified_v1_432.csv` | `qualified_v1` = `z(-foldx_ddg)`, population SD as `statistics.pstdev` does, matching the package's own `_z` |
| `measure_sd.py` | `runs/track_a_foldx/foldx_sd_sample.csv` | `FOLDX_RUNS=3`, 60 candidates: ranks 1-30 plus 30 evenly spaced over 31-432 |
| `measure_runs20.py` | `runs/track_a_foldx/foldx_runs20_top50.csv` | `FOLDX_RUNS=20`, top 50. Sorted best-first after the run and given a `rank_n20` column; the script itself appends rows as the thread pool finishes them, which is arrival order and carries no information |
| `sd_report.py` | the n=3 band figures | median and p90 SD, designs within 1 and 2 SD of the rank-15 cutoff |
| `summary_432.py` | the rerank distribution figures | per-group ddG summaries and rank distributions |

`logs/` holds the run logs, with per-candidate completion timestamps. Wall times on a
24-core box: 60 candidates at 3 runs, 8 workers, 1,027 s; 50 candidates at 20 runs, 12
workers, 4,570 s.

`tiers.py`, `boundary.py` and `spread.py` were deleted on 2026-11-25. Each printed one
table, and notebook 07 recomputes every figure they produced from the same CSVs — verified
against its stored outputs before they were removed. `sd_report.py` and `summary_432.py`
stayed because they are the only source of two numbers the ledger cites: the counts within
1 and 2 SD of the rank-15 cutoff, and the rank distribution of the 32 "ours".

One note on reading the tier figures in notebook 07: the count is definition-dependent.
Pairwise "not resolvable" is not transitive, so its connected components chain the whole
list into one blob — the greedy leader definition is the one to read, and a
Bonferroni-corrected variant is printed alongside it.

These are measurements, not pipeline steps. Nothing in `design_pipeline` imports them.
