#!/usr/bin/env bash
# DockQ: optional, CPU, tiny.
#
# The mildest separation in the set and worth stating anyway: DockQ brings
# biopython, and nothing in this package imports biopython. Putting it in Env A
# would add a dependency that only one wrapped function needs, so it gets its own
# 93 MB environment instead. `design_pipeline.metrics` imports DockQ lazily, so
# interface_rmsd, epitope_recall and every diagnostic keep working without it.

. "$(dirname "$0")/_common.sh"

usage() {
  cat >&2 <<'USAGE'
provision/dockq.sh [options]

Builds the DockQ environment into tools/dockq, against the metrics extra.

CPU only - --cuda is accepted and ignored.
USAGE
  usage_common
}

parse_args "$@"
need_uv

PY="$(python_version 3.12)"
VENV="$TOOLS/dockq"

step "DockQ: Python $PY at $VENV"
uv venv --python "$PY" "$VENV"

step "installing the metrics extra"
VIRTUAL_ENV="$VENV" uv pip install -e "$REPO_ROOT[metrics]"

step "checking the entry point"
[ -x "$VENV/bin/DockQ" ] || say "note: no DockQ entry point in $VENV/bin; the package is imported rather than shelled out to, so this is not fatal."
"$VENV/bin/python" -c "import DockQ; print('DockQ imports')" >&2 \
  || die "DockQ does not import in $VENV"

report "DockQ" "PATH=$VENV/bin:\$PATH"
say ""
say "No DESIGN_PIPELINE_* variable: design_pipeline.metrics imports DockQ in"
say "process, so run the metrics half from this environment. There is no route"
say "entry for it either - it is not an adapter."
