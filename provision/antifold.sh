#!/usr/bin/env bash
# AntiFold (plus AbLang2, which shares the environment).
#
# Installed from the GitHub source, NOT from PyPI, and the distinction matters:
#
#   PyPI    every published release pins biotite==0.38 against this package's
#           biotite>=2.0. pip calls that ResolutionImpossible (verbatim in
#           INSTALL.md). This is the constraint that forces a separate
#           environment at all - and it is NOT a pin to apply to the one we
#           build. Installing biotite==0.38 here would break it.
#   source  0.4.0 and later need Python >= 3.13 but pin neither biotite nor
#           numpy, so they install next to biotite 2.0.0. That is what
#           freeze_antifold.txt measures, and Python 3.13 is this environment's
#           actual separation driver, since Env A is 3.12.
#
# The commit is UNKNOWN: the node's copy had no .git, so the tree the existing
# AntiFold results came from cannot be identified. --commit pins the next one.

. "$(dirname "$0")/_common.sh"

usage() {
  cat >&2 <<'USAGE'
provision/antifold.sh [options]

Builds the AntiFold environment: Python 3.13, AntiFold from the GitHub source,
AbLang2 alongside it.

  --commit <sha>  pin the AntiFold tree. Default: the default branch, because
                  the commit the existing results came from was never recorded
                  (manifest.toml says UNKNOWN, and means it).
  --cuda <v>      the torch variant. This environment needs a GPU.
USAGE
  usage_common
}

parse_args "$@"
need_uv
need git "Needed to clone the AntiFold source; there is no usable PyPI release."

PY="$(python_version 3.13)"
ENVDIR="$TOOLS/antifold"
VENV="$ENVDIR/venv"
SRC="$ENVDIR/src"
URL="https://github.com/oxpig/AntiFold"

step "AntiFold: Python $PY at $VENV, source at $SRC"
mkdir -p "$ENVDIR"

step "cloning $URL"
if [ -d "$SRC/.git" ]; then
  say "$SRC is already a clone; fetching"
  git -C "$SRC" fetch --all --tags
elif [ -e "$SRC" ]; then
  die "$SRC exists and is not a git clone. That is how the node's copy ended up unidentifiable - move it aside rather than installing over it."
else
  git clone "$URL" "$SRC"
fi

if [ -n "$COMMIT" ]; then
  step "pinning to $COMMIT"
  git -C "$SRC" checkout --detach "$COMMIT"
else
  say ""
  say "note: no --commit given, so this takes the default branch."
  say "The tree the existing AntiFold results came from is UNKNOWN - the node's"
  say "copy had no .git - so this build is not reproducing it, it is replacing"
  say "it. Record the commit below in provision/manifest.toml."
fi
RESOLVED="$(git -C "$SRC" rev-parse HEAD)"

step "creating the venv"
uv venv --python "$PY" "$VENV"

step "installing the source tree, editable, $(torch_backend_args)"
# Editable so that ANTIFOLD_SRC and the installed package are the same tree:
# the weights land in src/models/model.pt and the adapter looks for them there.
VIRTUAL_ENV="$VENV" uv pip install -e "$SRC" "$(torch_backend_args)"

step "checking the biotite that landed"
# The whole point of the source route. If biotite 0.38 is here, something
# resolved the PyPI package instead and this environment is the wrong one.
got="$(VIRTUAL_ENV="$VENV" uv pip show biotite 2>/dev/null | awk '/^Version:/{print $2}')"
say "biotite $got"
case "$got" in
  0.38*) die "biotite $got is the PyPI AntiFold's pin, which means the source install did not take. The source route (0.4.0+) does not pin biotite and should give you 2.x. Check what $SRC/pyproject.toml or setup.py actually requires." ;;
  "") die "biotite is not installed in $VENV; the AntiFold install did not complete" ;;
esac

step "AbLang2, in the same environment"
# Same environment deliberately: it is a second sequence-only scorer with the
# same torch, and a third environment for one 179 MB model is not worth it.
VIRTUAL_ENV="$VENV" uv pip install "ablang2==0.2.1" "$(torch_backend_args)"

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

say ""
say "AntiFold tree now at: $RESOLVED"
say "Put that in provision/manifest.toml under [env.antifold.source], replacing UNKNOWN."

report "AntiFold" \
  "DESIGN_PIPELINE_ANTIFOLD=$VENV/bin/python" \
  "ANTIFOLD_SRC=$SRC"
say ""
say "model.pt (567 MB) downloads into $SRC/models on first use; pass weights= to"
say "reuse a download. ANTIFOLD_SRC is what the repo's antifold_score module reads."
