#!/usr/bin/env bash
# Run this on the laptop. Submits a GPU job and, after it succeeds, a report job that publishes to
# the `results` branch. Both jobs are pinned to the node that holds the install (node-local scratch).
#
#   slurm/submit_with_report.sh RUN_ID GPU_JOB_SCRIPT [-- REPORT_ARGS...]
#
# The GPU job script runs with RUN_ID and RUN_DIR (default $CODE/runs/$RUN_ID) in its environment and
# must write its design_pipeline output under $RUN_DIR (pass output_dir=... to the steps).
# Example:  slurm/submit_with_report.sh boltz_wt_01 ./my_gpu_job.sh -- --parent wt --label "parent only"
#
# REPORT_ARGS go through sbatch --export, so they must not contain commas.
# Environment: CLUSTER_SSH (default foertsch@fgcz-c-043), NODE (default fgcz-r-023), PARTITION (GPU_L40S),
#              GPU_TIME (default 02:00:00), CODE (default /scratch/foertsch/protecher)
set -euo pipefail

CLUSTER_SSH=${CLUSTER_SSH:-foertsch@fgcz-c-043}
NODE=${NODE:-fgcz-r-023}
PARTITION=${PARTITION:-GPU_L40S}
GPU_TIME=${GPU_TIME:-02:00:00}
CODE=${CODE:-/scratch/foertsch/protecher}
HERE=$(cd "$(dirname "$0")" && pwd)

[ $# -ge 2 ] || { sed -n '2,12p' "$0"; exit 1; }
RUN_ID=$1; GPU_SCRIPT=$2; shift 2
[ "${1:-}" = "--" ] && shift
REPORT_ARGS="$*"
[[ "$RUN_ID" =~ ^[A-Za-z0-9_.-]+$ ]] || { echo "RUN_ID may only use letters, digits, . _ -"; exit 1; }
[ -f "$GPU_SCRIPT" ] || { echo "no such file: $GPU_SCRIPT"; exit 1; }
RUN_DIR=${RUN_DIR:-$CODE/runs/$RUN_ID}

LOGS='$HOME/protecher_logs'
ssh -o BatchMode=yes "$CLUSTER_SSH" 'mkdir -p "$HOME/protecher_logs"'

GPU_JOB=$(ssh -o BatchMode=yes "$CLUSTER_SSH" \
  "sbatch --parsable -p $PARTITION -w $NODE --gres=gpu:1 --cpus-per-task=8 --mem=32G -t $GPU_TIME \
     -J run_$RUN_ID -o $LOGS/run_${RUN_ID}_%j.out --export=ALL,RUN_ID=$RUN_ID,RUN_DIR=$RUN_DIR" < "$GPU_SCRIPT")

# afterok: the report job runs only if the GPU job finished with exit code 0.
REPORT_JOB=$(ssh -o BatchMode=yes "$CLUSTER_SSH" \
  "sbatch --parsable -p $PARTITION -w $NODE --cpus-per-task=2 --mem=8G -t 00:20:00 \
     --dependency=afterok:$GPU_JOB --kill-on-invalid-dep=yes \
     -J report_$RUN_ID -o $LOGS/report_${RUN_ID}_%j.out \
     --export=ALL,RUN_ID=$RUN_ID,RUN_DIR=$RUN_DIR,GPU_JOB_ID=$GPU_JOB,REPORT_ARGS='$REPORT_ARGS'" < "$HERE/report_job.sh")

echo "run $RUN_ID: GPU job $GPU_JOB, then report job $REPORT_JOB (only if the GPU job succeeds)"
echo "logs: ~/protecher_logs/ on the cluster; results: the 'results' branch on GitHub"
