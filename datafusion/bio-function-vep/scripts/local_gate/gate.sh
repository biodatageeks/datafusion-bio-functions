#!/usr/bin/env bash
# Full local gate: build a baseline and a candidate, run every arm, print the verdict.
#
# Usage: gate.sh <baseline-functions-checkout> <candidate-functions-checkout> <gate-dir>
#                --target-pct HG002=<p>,ClinVar=<p> [--bgzf-identity]
#                [--candidate-env ENV=VALUE ...]
#
# Per dataset (hg002, then clinvar), with the merged cache:
#   1. candidate warm-up, 8 workers, annotate only (discarded)
#   2. baseline, candidate, candidate, baseline, 8 workers, annotate only
#   3. candidate parity against the VEP golden, 8 workers
#   4. candidate at 4 workers and at 1 worker, baseline at 1 worker, annotate only
#   5. with --bgzf-identity: candidate, 8 workers, BGZF output
# and then report.py over all arms, written to <gate-dir>/verdict.txt.
#
# Resumable: finished chromosomes and up-to-date builds are skipped.
# GATE_HG002_CHROMS and GATE_CLINVAR_CHROMS restrict the chromosomes, for a rehearsal of
# the driver only; a verdict from a restricted run is not a gate result.
# Exit status: that of report.py (0 pass, 1 fail), or 3 when a run stopped on low disk.
set -euo pipefail

usage() {
  sed -n '2,7p' "${BASH_SOURCE[0]}" | sed 's/^# \{0,1\}//' >&2
  exit 2
}
[ $# -ge 3 ] || usage
BASELINE_SRC=$1 CANDIDATE_SRC=$2 GATE=$3
shift 3
TARGETS="" BGZF_IDENTITY="" CANDIDATE_ENV=()
while [ $# -gt 0 ]; do
  case $1 in
    --target-pct)
      [ $# -ge 2 ] || usage
      TARGETS=$2
      shift
      ;;
    --bgzf-identity) BGZF_IDENTITY=1 ;;
    --candidate-env)
      shift
      while [ $# -gt 0 ] && [[ $1 == *=* ]] && [[ $1 != --* ]]; do
        CANDIDATE_ENV+=("$1")
        shift
      done
      continue
      ;;
    *) usage ;;
  esac
  shift
done
[ -n "$TARGETS" ] || usage

HERE=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd -P)
mkdir -p "$GATE"
GATE=$(cd "$GATE" && pwd -P)
log() { echo "[$(date +%H:%M:%S)] $*"; }

# Build unless the venv already holds a clean build of the checkout's HEAD.
build() { # <functions-checkout> <venv-dir>
  local head built
  head=$(git -C "$1" rev-parse HEAD)
  built=$(sed -n 's/.*"functions_sha": "\([0-9a-f]*\)".*/\1/p' "$2/GATE_IDENTITY.json" 2> /dev/null || true)
  if [ "$built" = "$head" ] && [ -z "$(git -C "$1" status --porcelain)" ] \
    && grep -q '"functions_dirty": false' "$2/GATE_IDENTITY.json"; then
    log "build $(basename "$2"): up to date at ${head:0:7}"
    return
  fi
  log "build $(basename "$2") from $1 at ${head:0:7}"
  "$HERE/build_venv.sh" "$1" "$2" > "$2.identity.log"
}
build "$BASELINE_SRC" "$GATE/venv-baseline"
build "$CANDIDATE_SRC" "$GATE/venv-candidate"

arm() { # <venv name> <arm name> <dataset> <workers> <mode> [run_arm.sh arguments...]
  local venv=$1 name=$2 dataset=$3 workers=$4 mode=$5 code=0
  shift 5
  local chroms=${GATE_HG002_CHROMS:-}
  if [ "$dataset" = clinvar ]; then chroms=${GATE_CLINVAR_CHROMS:-}; fi
  if [ -n "$chroms" ]; then set -- --chroms "$chroms" "$@"; fi
  log "$dataset $name: $venv, $workers workers, $mode $*"
  "$HERE/run_arm.sh" "$GATE/venv-$venv" "$GATE/$dataset/$name" "$dataset" "$workers" "$mode" "$@" \
    >> "$GATE/$dataset.$name.log" 2>&1 || code=$?
  if [ "$code" -ge 2 ]; then
    log "$dataset $name stopped with exit $code; see $GATE/$dataset.$name.log"
    exit "$code"
  fi
  if [ "$code" -ne 0 ]; then log "$dataset $name finished with failures"; fi
}
CAND=(${CANDIDATE_ENV[@]+"${CANDIDATE_ENV[@]}"})

for dataset in hg002 clinvar; do
  arm candidate warmup "$dataset" 8 time ${CAND[@]+"${CAND[@]}"}
  arm baseline baseline-r1 "$dataset" 8 time
  arm candidate candidate-r1 "$dataset" 8 time ${CAND[@]+"${CAND[@]}"}
  arm candidate candidate-r2 "$dataset" 8 time ${CAND[@]+"${CAND[@]}"}
  arm baseline baseline-r2 "$dataset" 8 time
  arm candidate candidate-parity "$dataset" 8 parity ${CAND[@]+"${CAND[@]}"}
  arm candidate candidate-w4 "$dataset" 4 time ${CAND[@]+"${CAND[@]}"}
  arm candidate candidate-w1 "$dataset" 1 time ${CAND[@]+"${CAND[@]}"}
  arm baseline baseline-w1 "$dataset" 1 time
  if [ -n "$BGZF_IDENTITY" ]; then
    arm candidate candidate-bgzf "$dataset" 8 time --bgzf ${CAND[@]+"${CAND[@]}"}
  fi
done

dirs() { # <arm name> -> the arm's directory for both datasets
  printf '%s\n' "$GATE/hg002/$1" "$GATE/clinvar/$1"
}
IDENTITY=($(dirs candidate-w4))
if [ -n "$BGZF_IDENTITY" ]; then IDENTITY+=($(dirs candidate-bgzf)); fi
code=0
python3 "$HERE/report.py" \
  --baseline $(dirs baseline-r1) $(dirs baseline-r2) \
  --candidate $(dirs candidate-r1) $(dirs candidate-r2) \
  --parity $(dirs candidate-parity) \
  --identity "${IDENTITY[@]}" \
  --baseline-1w $(dirs baseline-w1) --candidate-1w $(dirs candidate-w1) \
  --target-pct "$TARGETS" > "$GATE/verdict.txt" || code=$?
cat "$GATE/verdict.txt"
exit "$code"
