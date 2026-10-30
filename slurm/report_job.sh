#!/usr/bin/env bash
# Slurm job body: build one run's report and push it to the `results` branch.
# Submitted by slurm/submit_with_report.sh with --dependency=afterok, so it only runs after a
# successful GPU job, on the same node (run outputs live on that node's local /scratch).
#
# Environment (set by the submit script):
#   RUN_ID, RUN_DIR   the run to report on
#   REPORT_ARGS       extra arguments for `python -m run_report build` (e.g. "--parent wt --label test")
#   GPU_JOB_ID        recorded in the report
# Optional overrides: CODE, RESULTS, TOOLS
set -euo pipefail

CODE=${CODE:-/scratch/foertsch/protecher}
RESULTS=${RESULTS:-/scratch/foertsch/protecher_results}
TOOLS=${TOOLS:-/scratch/foertsch/protecher_tools}
: "${RUN_ID:?RUN_ID is not set}" "${RUN_DIR:?RUN_DIR is not set}"

. "$TOOLS/env.sh"
export PYTHONPATH="$CODE"
PY="$CODE/.venv/bin/python"
[ -d "$RUN_DIR" ] || { echo "STOP: $RUN_DIR does not exist on this node"; exit 1; }
[ -d "$RESULTS/.git" ] || { echo "STOP: $RESULTS is not set up, run slurm/setup_results_clone.sh once"; exit 1; }

exec 9>"$RESULTS/.lock"
flock 9  # one report job at a time on this node, so pushes do not race each other

cd "$RESULTS"
for attempt in 1 2 3; do
  git fetch -q origin results
  git reset -q --hard origin/results
  OUT="reports/$RUN_ID"
  if [ -e "$OUT" ]; then echo "STOP: report $RUN_ID already exists on the results branch"; exit 1; fi

  # shellcheck disable=SC2086
  "$PY" -m run_report build "$RUN_DIR" --run-id "$RUN_ID" --out "$OUT" --history reports \
        --code-dir "$CODE" --thresholds "$CODE/thresholds.toml" ${REPORT_ARGS:-} | tee /tmp/run_report_$$.log
  VERDICT=$(sed -n 's/^VERDICT: //p' /tmp/run_report_$$.log | tail -1)
  rm -f /tmp/run_report_$$.log
  "$PY" -m run_report leaderboard reports --out LEADERBOARD.md

  BIG=$(find reports LEADERBOARD.md -type f -size +2M | head -3)
  if [ -n "$BIG" ]; then echo "STOP: files over 2 MB, not committing: $BIG"; exit 3; fi

  git add reports LEADERBOARD.md
  git commit -q -m "report $RUN_ID: ${VERDICT:-no verdict}"
  if git push -q origin HEAD:results; then
    echo "pushed report $RUN_ID to the results branch ($(git rev-parse --short HEAD))"
    exit 0
  fi
  echo "push rejected (attempt $attempt), fetching and rebuilding"
  sleep $((attempt * 3))
done
echo "STOP: could not push after 3 attempts; the report is in $RESULTS/reports/$RUN_ID only"
exit 4
