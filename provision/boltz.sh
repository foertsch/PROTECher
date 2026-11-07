#!/usr/bin/env bash
# Boltz-2: the one environment replayed from its freeze, verbatim.
#
# Every other script here installs against the constraint that forces the
# separation, because a constraint survives a version bump and a pin does not.
# This one does the opposite, on purpose.
#
# freeze_venv-boltz.txt holds nvidia-cublas==13.4.1.3 AND
# nvidia-cublas-cu12==12.9.2.10, and the same doubling for nvidia-cuda-nvrtc:
# a mixed CUDA 12/13 wheel set, almost certainly cuequivariance-ops-* pulling
# cu12 alongside torch 2.15's cu13. It worked on the L40S. A fresh resolve is not
# guaranteed to reproduce that pairing, and this is the environment whose numbers
# went into the final ranking - so here the exact set is worth more than a
# readable constraint.

. "$(dirname "$0")/_common.sh"

usage() {
  cat >&2 <<'USAGE'
provision/boltz.sh [options]

Replays docs/cluster_freezes/freeze_venv-boltz.txt verbatim into .venv-boltz.

  --cuda is REFUSED here, not ignored: a verbatim replay has its CUDA variant
  baked into the wheel set, and accepting the flag silently would imply
  otherwise. Use provision/chai.sh, or install '.[boltz]' by hand, if you need
  a different variant and are prepared to resolve it fresh.
USAGE
  usage_common
}

# Parsed before the refusal so that --help still works and an unknown flag still
# fails the same way as everywhere else.
parse_args "$@"
[ "$CUDA" = "auto" ] || die "--cuda=$CUDA was given, but this script replays freeze_venv-boltz.txt verbatim and its CUDA variant is fixed by the wheel set (a mixed cu12/cu13 one - see the header). To choose a variant you have to resolve fresh: uv pip install -e '$REPO_ROOT[boltz]' --torch-backend=$CUDA, and then the versions are no longer the ones the final ranking came from."

need_uv

PY="$(python_version 3.12)"
VENV="$PREFIX/.venv-boltz"
FREEZE="$REPO_ROOT/docs/cluster_freezes/freeze_venv-boltz.txt"

step "Boltz-2: Python $PY at $VENV"
uv venv --python "$PY" "$VENV"

step "replaying $(basename "$FREEZE")"
say "$(freeze_filter "$FREEZE" | wc -l | tr -d ' ') installable lines (the banner and any -e / @ file lines are filtered)"
# Piped rather than written to a temporary file: the filter is part of the
# record, and a temp file invites someone editing it and wondering why the
# versions moved.
freeze_filter "$FREEZE" | VIRTUAL_ENV="$VENV" uv pip install -r -

step "this package, editable, so the adapters are importable here too"
# The freeze's `-e file:///scratch/foertsch/protecher` line is what the filter
# drops; this is that line, pointed at this clone instead.
VIRTUAL_ENV="$VENV" uv pip install -e "$REPO_ROOT"

step "checking torch sees the GPU"
# The point of this environment. A CPU-only torch build resolves without
# complaint and then every fold is a silent 100x slower, so it is checked here
# rather than discovered in a screen.
"$VENV/bin/python" - <<'PY' || die "torch in $VENV cannot see a GPU. On the cluster this script has to run through sbatch on the GPU node, not on the login node."
import sys
import torch

print(f"torch {torch.__version__}, cuda {torch.version.cuda}", file=sys.stderr)
if not torch.cuda.is_available():
    print("torch.cuda.is_available() is False", file=sys.stderr)
    raise SystemExit(1)
print(f"{torch.cuda.device_count()} device(s): {torch.cuda.get_device_name(0)}", file=sys.stderr)
PY

step "the boltz entry point"
[ -x "$VENV/bin/boltz" ] || die "no boltz entry point in $VENV/bin"
say "$VENV/bin/boltz"

report "Boltz-2" \
  "DESIGN_PIPELINE_BOLTZ=$VENV/bin/boltz" \
  "BOLTZ_CACHE=${BOLTZ_CACHE:-<an absolute path with 6.3 GB free>}"
say ""
say "BOLTZ_CACHE must be absolute - that is Boltz's own requirement, and the"
say "adapter raises on a relative value rather than letting it fail inside the"
say "tool. The weights (6.3 GB) download on the first prediction."
