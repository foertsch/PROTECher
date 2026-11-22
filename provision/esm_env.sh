#!/usr/bin/env bash
# ESM-1v and ESM-IF1: optional, for the scorer calibration rather than the screen.
#
# The separation constraint is the same negative one that defines Env A: this
# brings torch and transformers, and Env A has neither on purpose. There is no
# conflict between this and anything else here - it just cannot live in A.
#
# ESM-IF1 is the awkward half. torch-scatter ships no wheel for current torch and
# builds from source, which needs a compiler and --no-build-isolation; that is
# why esm_if is not in DEFAULT_SCORERS and why this script treats it as optional
# within an already-optional environment.

. "$(dirname "$0")/_common.sh"

usage() {
  cat >&2 <<'USAGE'
provision/esm_env.sh [options]

Builds the ESM environment into tools/esm_env: transformers for ESM-1v, and
fair-esm plus the torch-geometric stack for ESM-IF1.

  --cuda <v>  the torch variant.
USAGE
  usage_common
}

parse_args "$@"
need_uv

PY="$(python_version 3.12)"
VENV="$TOOLS/esm_env"

step "ESM: Python $PY at $VENV"
uv venv --python "$PY" "$VENV"

step "ESM-1v: transformers and torch, $(torch_backend_args)"
# The plm extra, which is also how Env A would get these if it were allowed to.
VIRTUAL_ENV="$VENV" uv pip install -e "$REPO_ROOT[plm]" "$(torch_backend_args)"

step "ESM-IF1: fair-esm and the torch-geometric stack"
VIRTUAL_ENV="$VENV" uv pip install fair-esm torch-geometric "$(torch_backend_args)"
# Separated out and allowed to fail: torch-scatter builds from source against
# the torch that is already installed, so it needs that torch present (hence a
# second command) and --no-build-isolation (hence the flag), and it needs a
# compiler, which a login node may not have. Without it ESM-1v still works and
# only esm_if is unavailable.
if ! VIRTUAL_ENV="$VENV" uv pip install torch-scatter --no-build-isolation; then
  say ""
  say "note: torch-scatter did not build, so esm_if is not available in this"
  say "environment. ESM-1v is unaffected. torch-scatter ships no wheel for"
  say "current torch and compiles against the installed one: it needs a"
  say "compiler and the CUDA toolkit headers. INSTALL.md records what the"
  say "resolver said when this was last attempted."
fi

step "checking torch sees the GPU"
"$VENV/bin/python" - <<'PY' || die "torch in $VENV cannot see a GPU. Run this on the GPU node."
import sys
import torch

print(f"torch {torch.__version__}, cuda {torch.version.cuda}", file=sys.stderr)
if not torch.cuda.is_available():
    print("torch.cuda.is_available() is False", file=sys.stderr)
    raise SystemExit(1)
print(f"{torch.cuda.device_count()} device(s): {torch.cuda.get_device_name(0)}", file=sys.stderr)
PY

report "ESM" \
  "HF_HOME=${HF_HOME:-<an absolute path for the ESM-1v download>}" \
  "TORCH_HOME=${TORCH_HOME:-<an absolute path for the ESM-IF1 checkpoint>}"
say ""
say "Both models download on first use: ESM-1v (650 M params) from Hugging Face"
say "into HF_HOME, ESM-IF1 (142 M params) through torch.hub into TORCH_HOME."
say "This environment has no DESIGN_PIPELINE_* variable of its own: esm_if and"
say "the PLM scorers are Env A routes, reached by import, so run them from here."
