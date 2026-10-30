# Run msa_batch_seeds_01: batch-vs-solo-seeds

**15/15 designs pass all gates · best ipTM 0.93 · parent 0.93 · earlier best run msa_cache_01: 7/7**

Created 2026-10-30 13:44 UTC · models boltz · seeds n/a · code `34bf9ba` · GPU job 411598

ok = reaches the target, ! = passes the gate but misses the target, ✗ = fails the gate. Thresholds are in `thresholds.toml`.

## Scorecard

| Metric | Median | Range | Parent | Best earlier median | ok / warn / fail |
|---|---|---|---|---|---|
| ipTM | 0.92 | 0.89 to 0.93 | 0.93 | 0.94 | 15 / 0 / 0 |
| pLDDT | 92.1 | 91.4 to 92.3 | 92.0 | 91.9 | 15 / 0 / 0 |
| pTM | 0.94 | 0.92 to 0.95 | 0.94 | 0.95 | 15 / 0 / 0 |

## Designs

| Design | Model | Group | Mutations | ipTM | pLDDT | pTM | Δ ipTM | Result |
|---|---|---|---|---|---|---|---|---|
| v4_model_0-batch_s2 | boltz | test_b | VH:S25T+VH:A88S+VL:Q3E+VL:V19I+VL:T69S | 0.93 ✓ | 92.1 ✓ | 0.95 ✓ | +0.01 | pass |
| wt_model_0-batch_s1 (parent) | boltz | parent | - | 0.93 ✓ | 92.1 ✓ | 0.95 ✓ | +0.01 | pass |
| v1_model_0-batch_s1 | boltz | test_a | VH:S17T+VH:T69S+VH:S71T+VH:D90E+VL:T85S | 0.93 ✓ | 91.7 ✓ | 0.94 ✓ | +0.00 | pass |
| wt_model_0-batch_s2 (parent) | boltz | parent | - | 0.93 ✓ | 91.9 ✓ | 0.94 ✓ | -0.00 | pass |
| v2_model_0-batch_s2 | boltz | test_a | VH:Q13E+VH:S21T+VH:A23S+VH:S75T+VL:S14T | 0.93 ✓ | 92.2 ✓ | 0.94 ✓ | -0.00 | pass |
| v2_model_0-batch_s3 | boltz | test_a | VH:Q13E+VH:S21T+VH:A23S+VH:S75T+VL:S14T | 0.92 ✓ | 91.4 ✓ | 0.94 ✓ | -0.00 | pass |
| v1_model_0-batch_s2 | boltz | test_a | VH:S17T+VH:T69S+VH:S71T+VH:D90E+VL:T85S | 0.92 ✓ | 91.6 ✓ | 0.94 ✓ | -0.00 | pass |
| v3_model_0-batch_s3 | boltz | test_b | VH:I70V+VH:A88S+VL:V15I+VL:S67T+VL:E81D | 0.92 ✓ | 92.1 ✓ | 0.94 ✓ | -0.01 | pass |
| v4_model_0-batch_s1 | boltz | test_b | VH:S25T+VH:A88S+VL:Q3E+VL:V19I+VL:T69S | 0.92 ✓ | 92.0 ✓ | 0.94 ✓ | -0.01 | pass |
| v1_model_0-batch_s3 | boltz | test_a | VH:S17T+VH:T69S+VH:S71T+VH:D90E+VL:T85S | 0.92 ✓ | 92.1 ✓ | 0.94 ✓ | -0.01 | pass |
| v3_model_0-batch_s2 | boltz | test_b | VH:I70V+VH:A88S+VL:V15I+VL:S67T+VL:E81D | 0.92 ✓ | 92.0 ✓ | 0.94 ✓ | -0.01 | pass |
| wt_model_0-batch_s3 (parent) | boltz | parent | - | 0.92 ✓ | 92.1 ✓ | 0.94 ✓ | -0.01 | pass |
| v2_model_0-batch_s1 | boltz | test_a | VH:Q13E+VH:S21T+VH:A23S+VH:S75T+VL:S14T | 0.92 ✓ | 92.3 ✓ | 0.94 ✓ | -0.01 | pass |
| v4_model_0-batch_s3 | boltz | test_b | VH:S25T+VH:A88S+VL:Q3E+VL:V19I+VL:T69S | 0.90 ✓ | 91.8 ✓ | 0.92 ✓ | -0.03 | pass |
| v3_model_0-batch_s1 | boltz | test_b | VH:I70V+VH:A88S+VL:V15I+VL:S67T+VL:E81D | 0.89 ✓ | 92.2 ✓ | 0.92 ✓ | -0.04 | pass |

Showing 15 of 15 designs, passing first, best ipTM first.

## Two metrics at once

![scatter](scatter.png)

## Compared with earlier runs

| Run | Created | Models | Designs | Pass | ipTM median | ipTM best |
|---|---|---|---|---|---|---|
| test_01 | 2026-10-30 09:06 UTC | boltz | 5 | 5 (100%) | 0.94 | 0.95 |
| msa_cache_01 | 2026-10-30 12:03 UTC | boltz | 7 | 7 (100%) | 0.93 | 0.94 |
| msa_helper_01 | 2026-10-30 12:15 UTC | boltz | 5 | 5 (100%) | 0.93 | 0.94 |
| msa_batch_01 | 2026-10-30 13:28 UTC | boltz | 5 | 5 (100%) | 0.93 | 0.94 |
| msa_batch_seeds_01 (this run) | 2026-10-30 13:44 UTC | boltz | 15 | 15 (100%) | 0.92 | 0.93 |

![compare](compare.png)

## Groups

test_a: 6, test_b: 6, parent: 3

## Files

`designs.csv` (one row per design with every score and grade), `summary.json` (numbers for the leaderboard).
