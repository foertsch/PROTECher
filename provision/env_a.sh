#!/usr/bin/env bash
# Env A: this package, the tests, the sequence-only scorers. CPU, no torch.
#
# The one environment with no separation problem to solve, because it is the one
# everything else is separated FROM. Its defining constraint is negative:
# importing design_pipeline must not import a model, so nothing that brings
# torch or jax belongs here.

. "$(dirname "$0")/_common.sh"

usage() {
  cat >&2 <<'USAGE'
provision/env_a.sh [options]

Builds Env A: an editable install of this package plus pytest, at the Python
version the manifest records (3.12). CPU only - --cuda is accepted and ignored,
so a driver script can pass one value to every environment.
USAGE
  usage_common
}

parse_args "$@"
need_uv

PY="$(python_version 3.12)"
VENV="$PREFIX/.venv"

step "Env A: Python $PY at $VENV"
uv venv --python "$PY" "$VENV"

step "installing this package, editable, plus the dev extra"
# Against pyproject, not against freeze_venv.txt. The freeze is the measurement;
# the dependency list is the intent, and it is deliberately short and light.
VIRTUAL_ENV="$VENV" uv pip install -e "$REPO_ROOT[dev]"

step "optional extras that stay CPU-only"
# numbering needs HMMER on PATH for ANARCI, so it is attempted and not required:
# without it the imgt column is empty and a warning says so, nothing else
# changes.
VIRTUAL_ENV="$VENV" uv pip install -e "$REPO_ROOT[numbering]" \
  || say "note: the numbering extra did not install (ANARCI needs HMMER on PATH). The imgt column will be empty."
VIRTUAL_ENV="$VENV" uv pip install -e "$REPO_ROOT[notebooks]"

step "checking that no torch came in"
# The constraint, asserted rather than assumed. plm/boltz/chai all bring torch,
# and any of them landing here by accident is the failure this guards.
if VIRTUAL_ENV="$VENV" uv pip list 2>/dev/null | grep -qiE '^(torch|jax)[[:space:]]'; then
  die "torch or jax is in Env A. That breaks the one rule this environment has: importing design_pipeline must not import a model. Check which extra pulled it in."
fi
say "none, as intended"

step "running the suite"
# 771 passed, 1 skipped on 2026-11-02. A different count is not necessarily
# wrong, but it is worth knowing about before anything else is built on top.
"$VENV/bin/pytest" "$REPO_ROOT/tests" -q || die "the suite did not pass in a fresh Env A"

report "Env A" "PATH=$VENV/bin:\$PATH"
