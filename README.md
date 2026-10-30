# Results

Reports from GPU runs, pushed automatically after each successful job. This branch holds results only.
The code is on `main`.

- `LEADERBOARD.md` compares every run, one row per run, and lists the best designs across runs.
- `reports/<run_id>/report.md` is the report for one run. Start with its first line (the verdict) and the scorecard.
- `reports/<run_id>/designs.csv` has every design with all scores and grades. `summary.json` has the numbers the leaderboard reads.

Thresholds that turn scores into ok, warn and fail are in `thresholds.toml` on `main`.
