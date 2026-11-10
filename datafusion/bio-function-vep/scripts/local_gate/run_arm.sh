#!/usr/bin/env bash
# One arm of a local gate: run vepyr's e2e runner once per chromosome with the merged cache.
#
# Usage: run_arm.sh <venv-dir> <arm-dir> <hg002|clinvar> <workers> <parity|time>
#                   [--bgzf] [--chroms "<c1> <c2> ..."] [ENV=VALUE ...]
#
#   parity  annotate and compare with the VEP golden by md5 (strict and canonical for
#           HG002, strict-clin-sig for ClinVar)
#   time    annotate only
#
# Each chromosome appends one row to <arm-dir>/contigs.tsv and its VCFs are deleted before
# the next one starts. A chromosome that already has a row with exit 0 is skipped, so an
# interrupted arm can be resumed with the same command. An arm directory belongs to one
# build and one set of arguments: resuming it with another is refused (exit 2).
#
# Environment (optional): GATE_ROOT (default $HOME/research/benchmarks/local-gate),
# DATA_VEPYR_DIR (default $HOME/workspace/data_vepyr). Every VEP_* and MIMALLOC_* variable
# of the calling shell is dropped; only the ENV=VALUE arguments reach the run.
set -euo pipefail

usage() {
  sed -n '2,6p' "${BASH_SOURCE[0]}" | sed 's/^# \{0,1\}//' >&2
  exit 2
}
[ $# -ge 5 ] || usage
VENV=$1 ARM=$2 DATASET=$3 WORKERS=$4 MODE=$5
shift 5
BGZF="" CHROMS="" EXTRA_ENV=()
while [ $# -gt 0 ]; do
  case $1 in
    --bgzf) BGZF=1 ;;
    --chroms)
      [ $# -ge 2 ] || usage
      CHROMS=$2
      shift
      ;;
    *=*) EXTRA_ENV+=("$1") ;;
    *) usage ;;
  esac
  shift
done
[[ $WORKERS =~ ^[1-9][0-9]*$ ]] || usage

GATE_ROOT=${GATE_ROOT:-$HOME/research/benchmarks/local-gate}
DATA_VEPYR_DIR=${DATA_VEPYR_DIR:-$HOME/workspace/data_vepyr}
SCRIPTS="$GATE_ROOT/vepyr-build/e2e-testing/scripts"
MIN_FREE_GIB=20 MAX_BG_CPU=200

case $DATASET in
  hg002)
    PROFILE=merged SUFFIX=merged MD5_MODE=both
    INPUT_ARGS=(--vcf "$DATA_VEPYR_DIR/input/HG002_normalized.vcf.gz" --no-normalize)
    DEFAULT_CHROMS=$(printf 'chr%s ' {1..22})
    ;;
  clinvar)
    PROFILE=clinvar_merged SUFFIX=clinvar_merged MD5_MODE=strict-clin-sig
    INPUT_ARGS=(--no-normalize)
    DEFAULT_CHROMS="$(printf '%s ' {1..22}) X Y MT"
    ;;
  *) usage ;;
esac
case $MODE in
  parity) MODE_ARGS=(--comparison-mode md5 --md5-mode "$MD5_MODE") ;;
  time) MODE_ARGS=(--skip-compare) ;;
  *) usage ;;
esac
[ -x "$VENV/bin/python" ] || { echo "error: no python in $VENV" >&2; exit 2; }
[ -f "$SCRIPTS/run_comparison.py" ] || { echo "error: $SCRIPTS/run_comparison.py missing; run build_venv.sh first" >&2; exit 2; }
[ -d "$DATA_VEPYR_DIR/cache/116_GRCh38_merged" ] || { echo "error: no merged cache under $DATA_VEPYR_DIR" >&2; exit 2; }

mkdir -p "$ARM/logs"
ARM=$(cd "$ARM" && pwd -P)
TSV="$ARM/contigs.tsv"
[ -f "$TSV" ] || printf 'dataset\tcontig\tworkers\tmode\tcompression\texit\twall_s\trss_gib\ttime_s\trecords\tbody_md5\tverdict\tenv\tbg_cpu\n' > "$TSV"
COMPRESSION=plain EXT=.vcf
if [ -n "$BGZF" ]; then COMPRESSION=bgzf EXT=.vcf.gz; fi
ENV_TEXT=${EXTRA_ENV[*]:-}

# An arm records the build and the arguments of its first run; a later run must match
# them, or its old rows would be reported as results of the new build.
"$VENV/bin/python" -P - "$VENV" "$ARM" "$DATASET" "$WORKERS" "$MODE" "$COMPRESSION" "$ENV_TEXT" << 'EOF'
import datetime
import json
import os
import sys

venv, arm, dataset, workers, mode, compression, env = sys.argv[1:]
identity_path = os.path.join(venv, "GATE_IDENTITY.json")
identity = json.load(open(identity_path)) if os.path.exists(identity_path) else None
current = {
    "extension_sha256": (identity or {}).get("extension_sha256"),
    "dataset": dataset,
    "workers": int(workers),
    "mode": mode,
    "compression": compression,
    "env": env,
}
path = os.path.join(arm, "arm.json")
if os.path.exists(path):
    recorded = json.load(open(path))
    previous = dict(recorded, extension_sha256=(recorded.get("identity") or {}).get("extension_sha256"))
    differing = [key for key, value in current.items() if previous.get(key) != value]
    if differing:
        print(
            f"error: {arm} holds results of another build or other arguments "
            f"(differs in: {', '.join(differing)}); use a new arm directory",
            file=sys.stderr,
        )
        sys.exit(2)
else:
    json.dump(
        {
            "venv": venv,
            "identity": identity,
            "dataset": dataset,
            "workers": int(workers),
            "mode": mode,
            "compression": compression,
            "env": env,
            "started": datetime.datetime.now().astimezone().isoformat(timespec="seconds"),
        },
        open(path, "w"),
        indent=2,
    )
EOF

# Inherited switches must not leak into a measured run.
UNSET=(-u CONDA_PREFIX -u VIRTUAL_ENV -u PYTHONPATH)
while IFS='=' read -r name _; do
  case $name in VEP_* | MIMALLOC_*) UNSET+=(-u "$name") ;; esac
done < <(env)

free_gib() { df -g "$ARM" | awk 'NR == 2 { print $4 }'; }
bg_cpu() { ps -A -o %cpu= | awk '{ s += $1 } END { printf "%d", s }'; }
wait_quiet() {
  local n=0
  while [ "$(bg_cpu)" -gt "$MAX_BG_CPU" ] && [ "$n" -lt 90 ]; do
    sleep 2
    n=$((n + 1))
  done
}

failed=0
for CONTIG in ${CHROMS:-$DEFAULT_CHROMS}; do
  if awk -F'\t' -v d="$DATASET" -v c="$CONTIG" '$1 == d && $2 == c && $6 == "0" { found = 1 } END { exit !found }' "$TSV"; then
    continue
  fi
  # The volume returns the space of files deleted a moment ago (the previous chromosome's
  # output and shards) with a delay, so a low reading is re-checked for two minutes.
  free=$(free_gib) n=0
  while [ "$free" -lt "$MIN_FREE_GIB" ] && [ "$n" -lt 60 ]; do
    sleep 2
    n=$((n + 1))
    free=$(free_gib)
  done
  if [ "$free" -lt "$MIN_FREE_GIB" ]; then
    echo "STOP: $free GiB free on the volume of $ARM, need $MIN_FREE_GIB" >&2
    exit 3
  fi
  wait_quiet
  # Recorded with the row: the wait gives up after three minutes, and report.py fails a
  # timing run that started on a busy machine.
  bg=$(bg_cpu)

  N=${CONTIG#chr}
  WORK="$ARM/results/116/fast_chr$N"
  REPORT="$ARM/reports/fast_chr${N}_${SUFFIX}_116_report.json"
  OUT="$WORK/vepyr_parquet_chr${N}_${SUFFIX}${EXT}"
  LOG="$ARM/logs/$DATASET.$CONTIG.log"
  rm -f "${REPORT:?}"

  code=0
  (
    cd "$SCRIPTS"
    /usr/bin/time -l env "${UNSET[@]}" DATA_VEPYR_DIR="$DATA_VEPYR_DIR" ${EXTRA_ENV[@]+"${EXTRA_ENV[@]}"} \
      "$VENV/bin/python" -P run_comparison.py --release 116 --profile "$PROFILE" \
      --chroms "$CONTIG" --workers "$WORKERS" --output-dir "$ARM" "${INPUT_ARGS[@]}" \
      "${MODE_ARGS[@]}" ${BGZF:+--bgzf}
  ) > "$LOG" 2>&1 || code=$?

  wall=$(awk '$2 == "real" && $4 == "user" && $6 == "sys" { v = $1 } END { print (v == "" ? "NA" : v) }' "$LOG")
  rss=$(awk '/maximum resident set size/ { v = $1 } END { if (v == "") print "NA"; else printf "%.3f", v / 1073741824 }' "$LOG")
  # Annotation time, record count and md5 verdict come from the runner's report: its log
  # words a concordant mode as MATCH and prints PASS only in the closing summary.
  read -r time_s records verdict < <("$VENV/bin/python" -P - "$REPORT" "$MODE" << 'EOF'
import json
import sys

time_s = records = None
verdict = "-" if sys.argv[2] == "time" else "NONE"
try:
    report = json.load(open(sys.argv[1]))
    time_s = report["annotation"].get("time_s")
    records = report["annotation"].get("output_variants")
    modes = report.get("md5") or {}
    if sys.argv[2] == "parity" and modes:
        concordant = all(m.get("body_match") is True and m.get("count_match") is True for m in modes.values())
        verdict = "PASS" if concordant else "FAIL"
except (OSError, ValueError, KeyError, AttributeError):
    pass
print("NA" if time_s is None else time_s, "NA" if records is None else records, verdict)
EOF
  )
  body_md5=NA
  if [ -f "$OUT" ]; then
    if [ -n "$BGZF" ]; then
      body_md5=$({ bgzip -dc "$OUT" | LC_ALL=C grep -v '^#' || true; } | md5 -q)
    else
      body_md5=$({ LC_ALL=C grep -v '^#' "$OUT" || true; } | md5 -q)
    fi
  fi
  printf '%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n' \
    "$DATASET" "$CONTIG" "$WORKERS" "$MODE" "$COMPRESSION" "$code" "$wall" "$rss" "$time_s" "$records" \
    "$body_md5" "$verdict" "$ENV_TEXT" "$bg" >> "$TSV"
  echo "$DATASET $CONTIG w$WORKERS $MODE $COMPRESSION: exit $code, annotate $time_s s, wall $wall s, rss $rss GiB, $records records, $verdict"

  # Whole-genome outputs do not fit on disk together: drop this chromosome's VCFs now.
  rm -f "${WORK:?}"/vepyr_* "${WORK:?}"/vep_*
  if [ "$code" -ne 0 ] || [ "$verdict" = FAIL ] || [ "$verdict" = NONE ]; then failed=1; fi
done
exit "$failed"
