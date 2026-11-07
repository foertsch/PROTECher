#!/usr/bin/env bash
# Remove the per-job repo copies from a node's scratch.
#
# Each sbatch job copied the whole repo tree plus a log into its own
# /scratch/$USER/protecher_<something> directory. About 60 of them accumulated at
# 1 to 12 MB each, listed as "not worth keeping" in
# docs/cluster_gpu_node_rebuild.md.
#
# DRY RUN BY DEFAULT. Nothing is deleted without --apply.
#
# The job directories look exactly like the directories that matter - same
# prefix, same parent - so this refuses by name rather than guessing: KEEP below
# is the list from the rebuild doc's "what exists only on the node scratch"
# table, and anything on it is never a candidate however well it matches.
#
# Note: /scratch/foertsch on fgcz-r-023 was already cleared on 2026-11-02, so on
# that node this finds nothing. It is here for the next one.

. "$(dirname "$0")/_common.sh"

usage() {
  cat >&2 <<'USAGE'
provision/clean_scratch.sh [--apply] [<scratch root>]

Lists the per-job repo copies under <scratch root> (default /scratch/$USER) and,
with --apply, deletes them.

  --apply   actually delete. Without it, nothing is removed and the candidates
            are printed with their sizes.

Never a candidate, whatever it matches: the screen output directories, the
pipeline runs, the funnel and loop-arm output, the FoldX result tables, the
ipSAE reference, the tool environments and the caches. Those are named
explicitly below.
USAGE
  usage_common
}

# A job directory is a copy of the repo tree plus a log. Everything in KEEP
# either holds results that are not in the repository or is a built environment.
KEEP_EXACT=(
  protecher
  protecher_tools
  protecher_cache
  protecher_msa_candidates
  protecher_inventory
  apptainer
  tools
)
KEEP_PREFIX=(
  protecher_screen
  protecher_track_a_run
  protecher_loop_arm
  protecher_pr23_foldx_results
  protecher_foldx_calib
  protecher_results
  esm1v_scan
)

parse_args_scratch() {
  ROOT=""
  while [ $# -gt 0 ]; do
    case "$1" in
      --apply) APPLY=1; shift ;;
      -h|--help) usage; exit 0 ;;
      -*) die "unknown option $1 (try --help)" ;;
      *) [ -z "$ROOT" ] || die "give at most one scratch root"; ROOT="$1"; shift ;;
    esac
  done
  ROOT="${ROOT:-/scratch/${USER:-$(id -un)}}"
}

parse_args_scratch "$@"

[ -d "$ROOT" ] || die "$ROOT is not a directory. Give the scratch root as an argument; on fgcz-r-023 it is /scratch/\$USER, and it is node-local so this has to run on the node itself."

# Refuse to run anywhere that is not a scratch root. The cost of a wrong
# argument here is someone's results, so the guard is on the path shape and not
# on what happens to be inside it.
case "$ROOT" in
  */scratch|*/scratch/*) ;;
  *) die "$ROOT does not look like a scratch path (no /scratch component). Refusing: this script deletes directories, and a mistyped root is not recoverable." ;;
esac
[ "$ROOT" != "/" ] || die "refusing to run against /"

keep() {
  local name="$1" k
  for k in "${KEEP_EXACT[@]}"; do
    [ "$name" = "$k" ] && return 0
  done
  for k in "${KEEP_PREFIX[@]}"; do
    case "$name" in "$k"*) return 0 ;; esac
  done
  return 1
}

step "scanning $ROOT"
candidates=()
while IFS= read -r path; do
  name="$(basename "$path")"
  if keep "$name"; then
    continue
  fi
  # A per-job copy, identified by content rather than by name alone: it holds
  # the package directory, which a results directory does not.
  if [ -d "$path/design_pipeline" ] || [ -f "$path/pyproject.toml" ]; then
    candidates+=("$path")
  fi
done < <(find "$ROOT" -mindepth 1 -maxdepth 1 -type d -name 'protecher_*' | sort)

if [ "${#candidates[@]}" -eq 0 ]; then
  say "no per-job repo copies found under $ROOT."
  say ""
  say "If this is fgcz-r-023: its /scratch/foertsch was cleared on 2026-11-02,"
  say "so there is nothing left to find. Nothing to do."
  exit 0
fi

step "${#candidates[@]} per-job repo copies"
total=0
for path in "${candidates[@]}"; do
  size="$(du -sk "$path" 2>/dev/null | cut -f1)"
  total=$((total + ${size:-0}))
  printf '  %8s KB  %s\n' "${size:-?}" "$path" >&2
done
say ""
say "total: $((total / 1024)) MB"

if [ "$APPLY" -ne 1 ]; then
  say ""
  say "DRY RUN. Nothing deleted. Re-run with --apply to remove these."
  say ""
  say "Not listed, and never touched by this script: ${KEEP_EXACT[*]} and"
  say "anything starting ${KEEP_PREFIX[*]}."
  exit 0
fi

step "deleting"
for path in "${candidates[@]}"; do
  say "  rm -rf $path"
  rm -rf -- "$path"
done
say ""
say "reclaimed about $((total / 1024)) MB."
say ""
say "The two big caches are deliberately not in scope here, because they are"
say "caches rather than job output and you may want them warm: "
say "  protecher_tools/uv-cache  (was 28 GB)"
say "  apptainer/cache           (was 15 GB)"
