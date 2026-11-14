#!/usr/bin/env bash
# core: this package, the tests, the sequence-only scorers. CPU.
#
# The one environment with no separation problem to solve, because it is the one
# everything else is separated FROM. Its defining constraint is negative:
# importing design_pipeline must not import a model. That is a statement about
# the imports, not about the package list - the lazy imports in
# adapters/_plm.py are what keep it true - so CPU torch is allowed here, and
# has to be: the default scorer needs it.

. "$(dirname "$0")/_common.sh"

usage() {
  cat >&2 <<'USAGE'
provision/core.sh [options]

Builds core: an editable install of this package, pytest, and the plm extra
that holds the sequence-only scorers, at the Python version the manifest
records (3.12). CPU only - torch comes in as a CPU build and --cuda is accepted
and ignored, so a driver script can pass one value to every environment.
USAGE
  usage_common
}

parse_args "$@"
need_uv

PY="$(python_version 3.12)"
VENV="$PREFIX/.venv"

step "core: Python $PY at $VENV"
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

step "the plm extra: the sequence-only scorers"
# Not optional, despite being an extra. DEFAULT_SCORERS is ("ablang2",), so
# ablang2 is in preflight's DEFAULT_CAMPAIGN and decides its exit status;
# routes.py puts it in core with no env_var - test_route_drift forbids one for
# a core route - and AbLang2MaskedLM._load imports `ablang2` in process. So
# there is no path by which core can reach it elsewhere: without this step
# preflight exits 1 on the only scorer it gates on. INSTALL.md's core row has
# always said this environment holds ablang2, antiberty, transformers + torch
# and fair-esm; it was this script that did not.
#
# CPU torch deliberately, so the inventory's "core: GPU no" stays true.
# AbLang2MaskedLM._load defaults to device="cpu", and a masked-marginal pass is
# one forward per position over a few hundred positions, not a fold.
#
# ablang2 is pinned to the version freeze_antifold.txt measured, in the same
# resolve rather than a second one, so that this copy and the one antifold.sh
# installs cannot drift apart and quietly change the scores.
VIRTUAL_ENV="$VENV" uv pip install -e "$REPO_ROOT[plm]" "ablang2==0.2.1" --torch-backend=cpu

step "checking that importing this package still imports no model"
# The constraint, asserted rather than assumed - but asserted as what it
# actually is. This was `uv pip list | grep torch`, which is a proxy for the
# rule and now the wrong one: torch is here on purpose. The rule is that
# importing design_pipeline must not import a model, and a fresh interpreter is
# the only honest way to ask - inside pytest the answer depends on what ran
# first. tests/test_steps.py asks the same question the same way.
VIRTUAL_ENV="$VENV" "$VENV/bin/python" -c '
import sys
import design_pipeline  # noqa: F401
leaked = sorted({m.split(".")[0] for m in sys.modules} & {"torch", "jax"})
if leaked:
    raise SystemExit("importing design_pipeline pulled in " + ", ".join(leaked))
' || die "importing design_pipeline imports a model. That is core's one rule, and a new top-level import somewhere broke it - the model imports belong inside the adapter that needs them."
say "it does not"

step "checking the default scorer is actually reachable"
# The row preflight's exit status hangs on. Checked here rather than left to
# preflight, because this is the environment that can fix it.
VIRTUAL_ENV="$VENV" "$VENV/bin/python" -c 'import ablang2'   || die "ablang2 did not install into core. DEFAULT_SCORERS is (\"ablang2\",), so preflight will exit 1 and the default campaign has no scorer."
say "ablang2 imports"

step "running the suite"
# 771 passed, 1 skipped on 2026-11-02. A different count is not necessarily
# wrong, but it is worth knowing about before anything else is built on top.
"$VENV/bin/pytest" "$REPO_ROOT/tests" -q || die "the suite did not pass in a fresh core"

report "core" "PATH=$VENV/bin:\$PATH"
