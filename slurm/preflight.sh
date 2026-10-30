#!/usr/bin/env bash
# Preflight for the hackathon morning: one OK, FAIL or SKIP line per check, at most 40 lines.
#
# From the laptop, in the repo root:
#   ssh login 'bash -s' < slurm/preflight.sh
# The exit status is the number of FAIL lines.
#
# Login node: free GPUs on the GPU node and, if set, the second GPU node, with who holds them, and your own queue.
# Node (one CPU-only job on the GPU node, no GPU, starts at once): model weights, alignment caches, structures, free
# /scratch, apptainer, outbound network, the Python environments, and the AntiFold environment if it was built.
# /scratch is node-local, so these node checks must run as a job there, not on the login node.
#
# Names and paths come from ~/.config/protecher/cluster.env on the cluster (template: slurm/cluster.env.example).
#
# Everything sits in main() so bash has parsed the whole script before any command can read the ssh stdin.

main() {
  local out="" fails=0 jid log state n i node cfg alloc free nstate holders
  log="$HOME/protecher_logs/preflight_$$.out"
  . "${PROTECHER_CONFIG:-$HOME/.config/protecher/cluster.env}" 2>/dev/null
  if [ -z "${PROTECHER_SCRATCH:-}" ] || [ -z "${PROTECHER_GPU_NODE:-}" ] || [ -z "${PROTECHER_GPU_PARTITION:-}" ]; then
    echo "FAIL cluster config: PROTECHER_SCRATCH, PROTECHER_GPU_NODE or PROTECHER_GPU_PARTITION not set (see slurm/cluster.env.example)"
    return 1
  fi
  mkdir -p "$HOME/protecher_logs"
  emit() { out+="$1"$'\n'; case "$1" in FAIL*) fails=$((fails + 1));; esac; }

  # Submit the node job first so it runs while the login-node checks are made.
  cat > "$log.sh" <<'NODE'
#!/usr/bin/env bash
. "${PROTECHER_CONFIG:-$HOME/.config/protecher/cluster.env}"
ok() { echo "OK   $*"; }
fail() { echo "FAIL $*"; }
skip() { echo "SKIP $*"; }
size() { stat -c %s "$1" 2>/dev/null || echo 0; }
gb() { awk -v b="$1" 'BEGIN { printf "%.1fG", b / 1e9 }'; }
S=$PROTECHER_SCRATCH

B=$S/protecher_cache/boltz
miss=""
for f in boltz2_conf.ckpt boltz2_aff.ckpt mols.tar; do [ "$(size $B/$f)" -gt 1000000000 ] || miss+=" $f"; done
[ -d $B/mols ] || miss+=" mols/"
[ -z "$miss" ] && ok "boltz weights: conf $(gb $(size $B/boltz2_conf.ckpt)) aff $(gb $(size $B/boltz2_aff.ckpt)) mols.tar, mols/ present" \
               || fail "boltz weights missing or short in $B:$miss"
[ -x $S/protecher_tools/bin/boltz ] && ok "boltz binary $S/protecher_tools/bin/boltz" || fail "boltz binary missing"

C=$S/protecher_cache/chai
miss=""
for f in models_v2/trunk.pt models_v2/diffusion_module.pt models_v2/confidence_head.pt conformers_v1.apkl; do
  [ "$(size $C/$f)" -gt 1000000 ] || miss+=" $f"; done
ls $C/esm/*.pt >/dev/null 2>&1 || miss+=" esm/*.pt"
[ -z "$miss" ] && ok "chai-1 weights present ($(du -sh $C | cut -f1))" || fail "chai-1 weights missing in $C:$miss"

M=$S/protecher_cache/msa_candidates
miss=""
for n in trastuzumab_her2 pertuzumab_her2 sacituzumab_trop2 cetuximab_egfr pembrolizumab_pd1 atezolizumab_pdl1; do
  [ -s $M/$n/parent.json ] && [ -s $M/$n/parent_0.csv ] && [ -s $M/$n/parent_1.csv ] || miss+=" $n"; done
[ -z "$miss" ] && ok "alignment caches: 6 of 6 candidates in $M ($(du -sh $M | cut -f1))" || fail "alignment cache missing for:$miss"

P=$S/protecher_cache/structures
miss=""
for id in 1N8Z 1S78 9PI9 1YY9 5GGS 5X8L; do [ -s $P/$id.cif ] || miss+=" $id"; done
[ -z "$miss" ] && ok "structures: 6 of 6 candidate PDB files in $P" || fail "structures missing:$miss"

avail=$(df -B1G --output=avail /scratch | tail -1 | tr -dc '0-9')
[ "$avail" -ge 200 ] && ok "scratch on $PROTECHER_GPU_NODE: ${avail} GB free" || fail "scratch on $PROTECHER_GPU_NODE: only ${avail} GB free"

export APPTAINER_CACHEDIR=$S/apptainer/cache APPTAINER_TMPDIR=$S/apptainer/tmp
if [ -f $S/apptainer/images/alpine.sif ] && apptainer exec $S/apptainer/images/alpine.sif true 2>/dev/null; then
  ok "apptainer $(apptainer --version | awk '{print $NF}') runs alpine.sif (APPTAINER_CACHEDIR=$APPTAINER_CACHEDIR)"
else fail "apptainer cannot run $S/apptainer/images/alpine.sif"; fi

net() { # name url, any HTTP answer except 000 counts: Docker Hub answers 401 on /v2/
  code=$(curl -s -o /dev/null -m 10 -w '%{http_code}' "$2" 2>/dev/null)
  case "$code" in 000|"") fail "network $1: no answer from $2";; *) ok "network $1: HTTP $code";; esac
}
net RCSB https://files.rcsb.org/download/1N8Z.cif
net HuggingFace https://huggingface.co
net PyPI https://pypi.org/simple/pip/
net GitHub https://github.com
net DockerHub https://registry-1.docker.io/v2/

cd $S/protecher
r=$(.venv/bin/python -c "import numpy, pandas, pyarrow, biotite, design_pipeline, msa_cache, run_report; print('numpy', numpy.__version__, 'pandas', pandas.__version__, 'pyarrow', pyarrow.__version__, 'biotite', biotite.__version__)" 2>&1 | tail -1)
case "$r" in numpy*) ok "core env .venv imports: $r";; *) fail "core env .venv imports: $r";; esac
r=$(.venv-boltz/bin/python -c "import boltz, torch; from importlib.metadata import version; print('boltz', version('boltz'), 'torch', torch.__version__)" 2>&1 | tail -1)
case "$r" in boltz*) ok "boltz env imports: $r";; *) fail "boltz env imports: $r";; esac
r=$(.venv-chai/bin/python -c "import chai_lab, torch; print('chai_lab ok, torch', torch.__version__)" 2>&1 | tail -1)
case "$r" in chai_lab*) ok "chai env imports: $r";; *) fail "chai env imports: $r";; esac

A=$S/protecher_tools/antifold
if [ -x $A/venv/bin/python ]; then
  r=$($A/venv/bin/python -c "import antifold.main; print('antifold imports')" 2>&1 | tail -1)
  w=$(size $A/src/models/model.pt)
  if [ "$r" = "antifold imports" ] && [ "$w" -gt 10000000 ]; then ok "antifold env: imports, weights $(gb $w)"
  else fail "antifold env: $r, weights $(gb $w)"; fi
else skip "antifold env not built ($A/venv)"; fi
NODE
  jid=$(sbatch --parsable -p "$PROTECHER_GPU_PARTITION" -w "$PROTECHER_GPU_NODE" --cpus-per-task=2 --mem=8G -t 00:05:00 \
        -J protecher_preflight -o "$log" "$log.sh" </dev/null 2>/dev/null) || jid=""
  rm -f "$log.sh"

  # Login-node checks.
  for node in "$PROTECHER_GPU_NODE" ${PROTECHER_SECOND_GPU_NODE:-}; do
    cfg=$(scontrol show node -o "$node" 2>/dev/null)
    if [ -z "$cfg" ]; then emit "FAIL gpu $node: node not known to Slurm"; continue; fi
    n=$(grep -o ' Gres=[^ ]*' <<<"$cfg" | grep -o 'gpu[^,]*' | head -1 | sed 's/.*://')
    alloc=$(grep -o 'AllocTRES=[^ ]*' <<<"$cfg" | grep -oE 'gres/gpu(:[A-Za-z0-9_]+)?=[0-9]+' | head -1 | sed 's/.*=//')
    nstate=$(grep -o ' State=[^ ]*' <<<"$cfg" | head -1 | sed 's/.*=//')
    n=${n:-0}; alloc=${alloc:-0}; free=$((n - alloc))
    holders=$(squeue -w "$node" -h -t R -o '%u %b' </dev/null 2>/dev/null | grep gpu | awk '{print $1}' | sort | uniq -c \
              | awk '{printf "%s%s x%s", s, $2, $1; s=", "}')
    case "$nstate" in *DOWN*|*DRAIN*|*FAIL*) emit "FAIL gpu $node: state $nstate"; continue;; esac
    if [ "$free" -ge 1 ]; then emit "OK   gpu $node: $free of $n free${holders:+ (running: $holders)}"
    else emit "FAIL gpu $node: 0 of $n free, held by ${holders:-unknown}"; fi
  done

  mine=$(squeue -u "$USER" -h -o '%i %T %M %j' </dev/null 2>/dev/null | grep -v ' protecher_preflight$' || true)
  n=$(grep -c . <<<"$mine" || true)
  if [ "$n" -eq 0 ]; then emit "OK   queue $USER: no other jobs"
  else
    emit "OK   queue $USER: $n job(s) besides this check"
    out+="$(head -4 <<<"$mine" | sed 's/^/     /')"$'\n'
  fi

  # Node job result.
  if [ -z "$jid" ]; then
    emit "FAIL node job: sbatch refused the job"
  else
    for i in $(seq 1 60); do
      state=$(squeue -j "$jid" -h -o %T </dev/null 2>/dev/null)
      [ -z "$state" ] && break
      sleep 2
    done
    if [ -n "$state" ]; then
      scancel "$jid" 2>/dev/null
      emit "FAIL node job $jid still $state after 120 s on $PROTECHER_GPU_NODE, cancelled"
    elif [ -s "$log" ]; then
      out+="$(grep -E '^(OK|FAIL|SKIP) ' "$log")"$'\n'
      fails=$((fails + $(grep -c '^FAIL ' "$log")))
    else
      emit "FAIL node job $jid left no output ($(sacct -j "$jid" -n -X -o State | head -1 | tr -d ' '))"
    fi
    rm -f "$log"
  fi

  { echo "PROTECher preflight $(date '+%Y-%m-%d %H:%M:%S') on $(hostname): $fails FAIL"; printf '%s' "$out"; } | cut -c1-200 | head -40
  return "$fails"
}

main "$@"
exit $?
