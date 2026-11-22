#!/usr/bin/env bash
# FoldX: unpack a licensed download and validate the layout.
#
# Nothing here can fetch FoldX - it is a licence-gated download from
# foldxsuite.crg.eu. What this does is unpack a zip you already hold into a
# layout that works, and then check it, because every FoldX failure we have had
# was a layout problem rather than an install problem: the binary not executable,
# or molecules/ one directory too deep so the rotamer data cannot be found by the
# relative path FoldX uses.

. "$(dirname "$0")/_common.sh"

usage() {
  cat >&2 <<'USAGE'
provision/foldx.sh --zip <path to the FoldX Linux zip> [options]

Get the zip from https://foldxsuite.crg.eu/ - register first (the academic
licence is free), then download FoldX 5 for Linux. The filename moves between
releases (foldx5_Linux_0.zip, foldx_Linux_1.zip, ...), which is why the binary
is located by glob rather than by name.

Unpacks the licensed FoldX zip into tools/foldx/ (gitignored) so the binary and
molecules/ are siblings, then checks that the binary is executable, molecules/
is beside it, and --help runs (the only invocation FoldX 5 exits 0 on
without a structure).

  --zip <path>   the licensed download. Required. Nothing here can fetch it.

CPU only - --cuda is accepted and ignored.
USAGE
  usage_common
}

parse_args "$@"

DEST="$TOOLS/foldx"

[ -n "$ZIP" ] || { usage; die "--zip is required: FoldX is a licensed download and this script cannot fetch it. Register at https://foldxsuite.crg.eu/ and download FoldX 5 for Linux. On the cluster the zip is already in ~/tools/foldx/foldx5_Linux_0.zip."; }
[ -f "$ZIP" ] || die "$ZIP does not exist"
need unzip "Install it, or unpack the zip by hand into $DEST."

step "unpacking $ZIP into $DEST"
mkdir -p "$DEST"
# -o so a re-run is idempotent; -q because the listing is long and the checks
# below are what actually matter.
unzip -o -q "$ZIP" -d "$DEST"

step "flattening, if the zip had a top-level directory"
# FoldX zips have shipped both ways across releases. The binary and molecules/
# have to end up siblings directly under $DEST, so if everything landed one
# level down, lift it.
if [ ! -e "$DEST/molecules" ]; then
  nested="$(find "$DEST" -mindepth 2 -maxdepth 2 -type d -name molecules -print -quit || true)"
  if [ -n "$nested" ]; then
    inner="$(dirname "$nested")"
    say "found molecules/ under $inner; lifting its contents to $DEST"
    # shellcheck disable=SC2086
    find "$inner" -mindepth 1 -maxdepth 1 -exec mv -t "$DEST" {} +
    rmdir "$inner" 2>/dev/null || true
  fi
fi

step "locating the binary"
# Named foldx_<expiry>, foldx, or foldx_<platform> depending on the release.
BINARY="$(find "$DEST" -maxdepth 1 -type f -name 'foldx*' ! -name '*.zip' ! -name '*.txt' ! -name '*.pdf' -print | sort | head -1 || true)"
[ -n "$BINARY" ] || die "no foldx binary directly under $DEST. Unpack the zip by hand so the binary and molecules/ are siblings there."
say "$BINARY"

# --------------------------------------------------------- the three checks

step "check 1/3: the binary is executable"
if [ ! -x "$BINARY" ]; then
  say "not executable; fixing with chmod +x"
  chmod +x "$BINARY"
fi
[ -x "$BINARY" ] || die "$BINARY is still not executable"
say "ok"

step "check 2/3: molecules/ is beside the binary"
# Not nested. FoldX reads its rotamer data by a path relative to the binary, and
# gets this wrong with an error that names neither the file nor the directory.
[ -d "$DEST/molecules" ] || die "molecules/ is not in $DEST, beside the binary. FoldX reads its rotamer data from there by relative path and fails obscurely without it. Check what the zip actually contained."
count="$(find "$DEST/molecules" -maxdepth 1 -type f | wc -l | tr -d ' ')"
[ "$count" -gt 0 ] || die "$DEST/molecules exists but is empty"
say "ok, $count files"

step "check 3/3: the binary runs"
# Run from $DEST, which is how the adapter runs it, so that a relative-path
# problem surfaces here rather than mid-campaign.
#
# --help, because it is the only invocation that exits 0 without a PDB. Measured
# on foldx_20270131: `--version` is an option its parser says REQUIRES an
# argument (255), `-version` is 255, and a bare run is 255 too - "No pdbs for
# the run found". This check used to use --version and therefore failed on a
# working binary. --help still catches what it is here for: a binary built for
# another libc, or one that cannot start at all. Whether FoldX can actually
# score is a question only a real PDB answers, and that is the campaign's job.
if ! ( cd "$DEST" && "$BINARY" --help >/dev/null 2>&1 ); then
  die "$BINARY --help did not exit 0 when run from $DEST. Run it by hand there to see what it says; a binary built for another libc lands here."
fi
say "ok"

# ------------------------------------------------------------- the expiry

step "licence expiry"
# The date is in the filename and nowhere else - the binary does not report it.
# preflight warns inside 30 days; this is the same read, said once at install.
stamp="$(printf '%s' "$(basename "$BINARY")" | grep -oE '20[0-9]{2}(0[1-9]|1[0-2])(0[1-9]|[12][0-9]|3[01])' | head -1 || true)"
if [ -n "$stamp" ]; then
  say "expires ${stamp:0:4}-${stamp:4:2}-${stamp:6:2}, read from the filename"
  say "An academic renewal is free but takes a few days. preflight warns inside 30 days."
else
  say "no date in the filename, so nothing to check. preflight will not warn either."
fi

report "FoldX" \
  "FOLDX_BINARY=$BINARY" \
  "DESIGN_PIPELINE_FOLDX=\$(command -v python)" \
  "DESIGN_PIPELINE_FOLDX_PDB=<the trimmed, repaired Fv plus antigen>"
say ""
say "DESIGN_PIPELINE_FOLDX_PDB has no sensible default: it is the structure every"
say "ddG is computed against, so it has to be named. preflight reports it unset."
