# Shared by every provision/*.sh. Sourced, not run.
#
# Holds the three things all of them need: the argument parsing (so --cuda means
# the same everywhere), the freeze-file filter, and the closing report that
# prints the export lines the environment needs. Nothing here installs anything.

# shellcheck shell=bash

set -euo pipefail

# Where environments go. The node's layout was protecher/.venv* for the venvs
# this repo owns and protecher_tools/<name> for the separate tools, so these two
# defaults reproduce it from inside a clone.
REPO_ROOT="${REPO_ROOT:-$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)}"
PREFIX="${PROTECHER_PREFIX:-$REPO_ROOT}"
TOOLS="${PROTECHER_TOOLS:-$REPO_ROOT/tools}"

CUDA="auto"
COMMIT=""
ZIP=""
APPLY=0
PYTHON_OVERRIDE=""

say() { printf '%s\n' "$*" >&2; }
step() { printf '\n== %s\n' "$*" >&2; }
die() { printf 'error: %s\n' "$*" >&2; exit 1; }

need() {
  command -v "$1" >/dev/null 2>&1 || die "$1 is not on PATH. $2"
}

need_uv() {
  need uv "Install it from https://docs.astral.sh/uv/ or use the system one (/usr/bin/uv on fgcz-r-023)."
}

usage_common() {
  cat >&2 <<'USAGE'
Common options:
  --cuda <variant>   auto (default), cpu, or a cuXXX moniker such as cu128.
                     Passed to uv as --torch-backend. Ignored by the CPU-only
                     environments, so one value can be passed to all of them.
  --python <version>  override the Python version this environment is built on.
  --prefix <dir>     where venvs this repo owns are created (default: repo root).
  --tools <dir>      where separate tool environments go (default: <repo>/tools).
  -h, --help         this message.
USAGE
}

parse_args() {
  while [ $# -gt 0 ]; do
    case "$1" in
      --cuda) CUDA="${2:?--cuda needs a value}"; shift 2 ;;
      --cuda=*) CUDA="${1#*=}"; shift ;;
      --commit) COMMIT="${2:?--commit needs a value}"; shift 2 ;;
      --commit=*) COMMIT="${1#*=}"; shift ;;
      --zip) ZIP="${2:?--zip needs a value}"; shift 2 ;;
      --zip=*) ZIP="${1#*=}"; shift ;;
      --python) PYTHON_OVERRIDE="${2:?--python needs a value}"; shift 2 ;;
      --python=*) PYTHON_OVERRIDE="${1#*=}"; shift ;;
      --prefix) PREFIX="${2:?--prefix needs a value}"; shift 2 ;;
      --prefix=*) PREFIX="${1#*=}"; shift ;;
      --tools) TOOLS="${2:?--tools needs a value}"; shift 2 ;;
      --tools=*) TOOLS="${1#*=}"; shift ;;
      --apply) APPLY=1; shift ;;
      -h|--help) usage; exit 0 ;;
      *) die "unknown option $1 (try --help)" ;;
    esac
  done
  case "$CUDA" in
    auto|cpu|cu[0-9][0-9][0-9]) ;;
    *) die "--cuda takes auto, cpu, or a cuXXX moniker such as cu128; got '$CUDA'" ;;
  esac
}

python_version() {
  # $1 is the version the manifest records for this environment.
  printf '%s' "${PYTHON_OVERRIDE:-$1}"
}

torch_backend_args() {
  # uv resolves the CUDA build through --torch-backend. No freeze in this repo
  # records an index URL or a +cuXXX tag, so the variant is never implied by the
  # pins and always has to be said out loud. 'auto' lets uv read the driver.
  printf '%s' "--torch-backend=$CUDA"
}

freeze_filter() {
  # Drop the three kinds of line a `uv pip freeze` carries that cannot be
  # installed from: uv's "Using Python ... environment at:" banner, an editable
  # install of this repo, and a non-editable `name @ file:///...` local path.
  # This is the filter recorded in docs/cluster_gpu_node_rebuild.md.
  [ -f "$1" ] || die "freeze file $1 not found"
  grep -vE '^(-e|Using)|@ file' "$1"
}

report() {
  # $1 the environment name, then any number of 'VAR=value' lines to export.
  step "$1 is built"
  shift
  if [ $# -gt 0 ]; then
    say "Add these to the environment (T/env.sh on the cluster):"
    for pair in "$@"; do
      say "  export $pair"
    done
  fi
  say ""
  say "Then, from Env A with those set: python -m design_pipeline.preflight"
}
