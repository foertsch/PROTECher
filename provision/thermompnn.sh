#!/usr/bin/env bash
# ThermoMPNN: optional, the stability gate.
#
# The odd one out twice over. Its separation constraint is the harshest here -
# environment.yaml is Python 3.10 with pytorch-cuda 11.7, two CUDA majors behind
# everything else - and it is the only environment built by mamba from a yaml
# rather than by uv from an index, because that yaml is what the tree ships and
# reproducing it by hand would be guesswork.
#
# The commit is UNKNOWN, as for AntiFold: the node's tree had no .git. Unlike
# AntiFold nothing even narrows the version, so --commit is the only way to get
# a build you can name later.

. "$(dirname "$0")/_common.sh"

usage() {
  cat >&2 <<'USAGE'
provision/thermompnn.sh [options]

Clones ThermoMPNN and builds its own environment.yaml with mamba.

  --commit <sha>  pin the tree. Default: the default branch. The commit the
                  existing results came from was never recorded.
  --cuda is accepted and IGNORED: the yaml pins pytorch-cuda 11.7 itself, and
  overriding that is how you get a tree that imports and does not run.
USAGE
  usage_common
}

parse_args "$@"
need git "Needed to clone the ThermoMPNN source."

MAMBA=""
# mamba or conda only: both have `env create -f`. micromamba does not, so it is
# not a drop-in here and silently falling back to it would fail further in.
for candidate in mamba conda; do
  if command -v "$candidate" >/dev/null 2>&1; then MAMBA="$candidate"; break; fi
done
[ -n "$MAMBA" ] || die "neither mamba nor conda is on PATH. ThermoMPNN ships an environment.yaml and that is the only recorded way to build it; uv cannot read a conda yaml, and micromamba has no 'env create'."

if [ "$CUDA" != "auto" ]; then
  say "note: --cuda=$CUDA ignored. environment.yaml pins pytorch-cuda 11.7 itself."
fi

ENVDIR="$TOOLS/thermompnn"
SRC="$ENVDIR/src"
ENVNAME="${THERMOMPNN_ENV:-thermompnn}"
URL="https://github.com/Kuhlman-Lab/ThermoMPNN"

step "ThermoMPNN source at $SRC"
mkdir -p "$ENVDIR"
if [ -d "$SRC/.git" ]; then
  say "already a clone; fetching"
  git -C "$SRC" fetch --all --tags
elif [ -e "$SRC" ]; then
  die "$SRC exists and is not a git clone. Move it aside - an unidentifiable tree is exactly the problem this script is trying not to repeat."
else
  git clone "$URL" "$SRC"
fi

if [ -n "$COMMIT" ]; then
  step "pinning to $COMMIT"
  git -C "$SRC" checkout --detach "$COMMIT"
else
  say ""
  say "note: no --commit given. The tree the existing ThermoMPNN results came"
  say "from is UNKNOWN and nothing narrows it, so this is a replacement rather"
  say "than a reproduction. Record the commit below in provision/manifest.toml."
fi
RESOLVED="$(git -C "$SRC" rev-parse HEAD)"

YAML="$SRC/environment.yaml"
[ -f "$YAML" ] || YAML="$SRC/environment.yml"
[ -f "$YAML" ] || die "no environment.yaml in $SRC. The tree's layout has moved; find its install instructions before going further."

step "building $ENVNAME from $(basename "$YAML") with $MAMBA"
"$MAMBA" env create -n "$ENVNAME" -f "$YAML" -y || die "$MAMBA could not build $ENVNAME from $YAML"

step "checking torch sees the GPU"
# Its own README warns the yaml can resolve a CPU-only torch build. That
# resolves without complaint and then the gate runs at a crawl, so it is checked
# here rather than discovered later.
if ! "$MAMBA" run -n "$ENVNAME" python -c "
import sys, torch
print(f'torch {torch.__version__}, cuda {torch.version.cuda}', file=sys.stderr)
sys.exit(0 if torch.cuda.is_available() else 1)
"; then
  die "torch in $ENVNAME cannot see a GPU. Its README warns environment.yaml can resolve a CPU-only build - that is what has happened. Rebuild with the pytorch channel ahead of defaults, or run this on the GPU node if that is the problem."
fi

ENTRY="$("$MAMBA" run -n "$ENVNAME" python -c 'import sys; print(sys.executable)' 2>/dev/null | tr -d '\r')"

say ""
say "ThermoMPNN tree now at: $RESOLVED"
say "Put that in provision/manifest.toml under [env.thermompnn.source], replacing UNKNOWN."

report "ThermoMPNN" "DESIGN_PIPELINE_THERMOMPNN=$ENTRY"
say ""
say "The adapter runs that interpreter against the tree in $SRC."
