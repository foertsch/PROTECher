#!/usr/bin/env bash
# Chai-1: optional, a second opinion against Boltz-2.
#
# Nothing requires this. It is here because the freeze exists and because having
# a second co-folding model is cheap insurance against a Boltz-specific
# artefact - not because any result depends on it.
#
# The separation constraint is torch: chai_lab pins 2.6 against boltz's 2.15, so
# the two folding models cannot share an environment with each other, let alone
# with Env A.
#
# Read the caveat below before wiring it into anything.

. "$(dirname "$0")/_common.sh"

usage() {
  cat >&2 <<'USAGE'
provision/chai.sh [options]

Builds the Chai-1 environment into .venv-chai, against the chai extra.

  --cuda <v>  the torch variant. Needs a GPU with compute capability 8.0 or
              newer: Chai-1 runs in bfloat16.
USAGE
  usage_common
}

parse_args "$@"
need_uv

PY="$(python_version 3.12)"
VENV="$PREFIX/.venv-chai"

step "Chai-1: Python $PY at $VENV"
uv venv --python "$PY" "$VENV"

step "installing the chai extra, $(torch_backend_args)"
# Against the extra, not the freeze: unlike Boltz this environment's numbers are
# in no result we keep, so a fresh resolve costs nothing. freeze_venv-chai.txt
# stays in docs/cluster_freezes as the record of what was measured
# (chai-lab 0.6.1, torch 2.6.0).
VIRTUAL_ENV="$VENV" uv pip install -e "$REPO_ROOT[chai]" "$(torch_backend_args)"

step "checking torch sees a bfloat16-capable GPU"
"$VENV/bin/python" - <<'PY' || die "torch in $VENV cannot see a suitable GPU. Run this on the GPU node; on pre-Ampere cards use boltz instead."
import sys
import torch

print(f"torch {torch.__version__}, cuda {torch.version.cuda}", file=sys.stderr)
if not torch.cuda.is_available():
    print("torch.cuda.is_available() is False", file=sys.stderr)
    raise SystemExit(1)
major, minor = torch.cuda.get_device_capability(0)
capability = float(f"{major}.{minor}")
print(f"{torch.cuda.get_device_name(0)}, compute {capability}", file=sys.stderr)
# Chai-1 runs in bfloat16, which arrives with Ampere (8.0). The adapter checks
# this too, via nvidia-smi; checking it at install means not finding out later.
if capability < 8.0:
    print("compute capability below 8.0: Chai-1 needs bfloat16", file=sys.stderr)
    raise SystemExit(1)
PY

report "Chai-1" \
  "DESIGN_PIPELINE_CHAI=$VENV/bin/python" \
  "CHAI_DOWNLOADS_DIR=${CHAI_DOWNLOADS_DIR:-<an absolute path with 6.9 GB free>}"
say ""
say "Caveat, and it is the reason DESIGN_PIPELINE_CHAI is unlike the others:"
say "the Chai adapter calls run_inference IN PROCESS, because the Chai CLI drops"
say "the PAE matrix that pae_interaction needs. So that variable is a preflight"
say "probe and nothing more - it lets preflight report chai from Env A, and it"
say "does not make chai runnable from Env A. Drive Chai from inside $VENV."
