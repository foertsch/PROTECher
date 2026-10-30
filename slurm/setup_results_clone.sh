#!/usr/bin/env bash
# One-off, run on the GPU node (node-local scratch): create the clone that report jobs publish from,
# and the `results` branch on GitHub if it does not exist. The branch holds only reports/ and
# LEADERBOARD.md, so it cannot break the code on main.
#   ssh login '. ~/.config/protecher/cluster.env; sbatch -p $PROTECHER_GPU_PARTITION -w $PROTECHER_GPU_NODE ...' < slurm/setup_results_clone.sh
set -euo pipefail
. "${PROTECHER_CONFIG:-$HOME/.config/protecher/cluster.env}"
: "${PROTECHER_SCRATCH:?PROTECHER_SCRATCH is not set, copy slurm/cluster.env.example to ~/.config/protecher/cluster.env}"
RESULTS=${RESULTS:-$PROTECHER_SCRATCH/protecher_results}
URL=${REPO_URL:-git@github.com:foertsch/PROTECher.git}

if [ -d "$RESULTS/.git" ]; then echo "already set up: $RESULTS"; exit 0; fi
git clone -q --no-checkout "$URL" "$RESULTS"
cd "$RESULTS"
git config user.name "Arion Förtsch"
git config user.email "68943888+foertsch@users.noreply.github.com"

if git ls-remote --exit-code --heads origin results >/dev/null 2>&1; then
  git checkout -q -b results origin/results
  echo "results branch exists on GitHub, checked out"
else
  git checkout -q --orphan results
  git rm -rq --cached . 2>/dev/null || true
  find . -mindepth 1 -maxdepth 1 ! -name .git -exec rm -rf {} +
  cat > README.md <<'README'
# Results

Reports from GPU runs, pushed automatically after each successful job. This branch holds results only.
The code is on `main`.

- `LEADERBOARD.md` compares every run, one row per run, and lists the best designs across runs.
- `reports/<run_id>/report.md` is the report for one run. Start with its first line (the verdict) and the scorecard.
- `reports/<run_id>/designs.csv` has every design with all scores and grades. `summary.json` has the numbers the leaderboard reads.

Thresholds that turn scores into ok, warn and fail are in `thresholds.toml` on `main`.
README
  git add README.md
  git commit -q -m "Start the results branch

Reports and the leaderboard from GPU runs are pushed here by slurm/report_job.sh.

Co-Authored-By: Claude Sonnet 5.5 <noreply@anthropic.com>"
  git push -q -u origin results
  echo "created and pushed the results branch"
fi
git log --oneline -3
