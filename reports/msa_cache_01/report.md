# Run msa_cache_01: cached-msa

**7/7 designs pass all gates · best ipTM 0.94 · parent 0.92 · earlier best run test_01: 5/5**

Created 2026-10-30 12:03 UTC · models boltz · seeds 0, 1, 2 · code `34bf9ba` · GPU job 411385

ok = reaches the target, ! = passes the gate but misses the target, ✗ = fails the gate. Thresholds are in `thresholds.toml`.

## Scorecard

| Metric | Median | Range | Parent | Best earlier median | ok / warn / fail |
|---|---|---|---|---|---|
| ipTM | 0.93 | 0.92 to 0.94 | 0.92 | 0.94 | 7 / 0 / 0 |
| pLDDT | 91.7 | 91.4 to 92.5 | 92.5 | 91.8 | 7 / 0 / 0 |
| pTM | 0.95 | 0.93 to 0.95 | 0.94 | 0.95 | 7 / 0 / 0 |

## Designs

| Design | Model | Group | Mutations | ipTM | pLDDT | pTM | Δ ipTM | Result |
|---|---|---|---|---|---|---|---|---|
| v1_model_0-e7f0fe38 | boltz | test_a | VH:S17T+VH:T69S+VH:S71T+VH:D90E+VL:T85S | 0.94 ✓ | 91.6 ✓ | 0.95 ✓ | +0.02 | pass |
| wt_model_0-e7f0fe38 | boltz | parent | - | 0.94 ✓ | 91.9 ✓ | 0.95 ✓ | +0.02 | pass |
| v4_model_0-e7f0fe38 | boltz | test_b | VH:S25T+VH:A88S+VL:Q3E+VL:V19I+VL:T69S | 0.93 ✓ | 91.7 ✓ | 0.95 ✓ | +0.01 | pass |
| wt_model_0-d4897021 | boltz | parent | - | 0.93 ✓ | 92.1 ✓ | 0.95 ✓ | +0.01 | pass |
| v2_model_0-e7f0fe38 | boltz | test_a | VH:Q13E+VH:S21T+VH:A23S+VH:S75T+VL:S14T | 0.92 ✓ | 91.7 ✓ | 0.94 ✓ | +0.00 | pass |
| wt_model_0-62ff411b (parent) | boltz | parent | - | 0.92 ✓ | 92.5 ✓ | 0.94 ✓ | +0.00 | pass |
| v3_model_0-e7f0fe38 | boltz | test_b | VH:I70V+VH:A88S+VL:V15I+VL:S67T+VL:E81D | 0.92 ✓ | 91.4 ✓ | 0.93 ✓ | -0.01 | pass |

Showing 7 of 7 designs, passing first, best ipTM first.

## Two metrics at once

![scatter](scatter.png)

## Compared with earlier runs

| Run | Created | Models | Designs | Pass | ipTM median | ipTM best |
|---|---|---|---|---|---|---|
| test_01 | 2026-10-30 09:06 UTC | boltz | 5 | 5 (100%) | 0.94 | 0.95 |
| msa_cache_01 (this run) | 2026-10-30 12:03 UTC | boltz | 7 | 7 (100%) | 0.93 | 0.94 |

![compare](compare.png)

## Groups

parent: 3, test_a: 2, test_b: 2

## Files

`designs.csv` (one row per design with every score and grade), `summary.json` (numbers for the leaderboard).
