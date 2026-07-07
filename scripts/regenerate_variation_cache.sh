#!/usr/bin/env bash
set -euo pipefail

CACHE_ROOT="/Users/mwiewior/workspace/data_vepyr/homo_sapiens_merged/115_GRCh38"
OUTPUT_DIR="/Users/mwiewior/workspace/data_vepyr/115_GRCh38_merged"
ENTITY_DIR="$OUTPUT_DIR/variation.lance"
BIN="./target/release/examples/build_lance_variation_chrom"
KIND="python scripts/variation_cache_schema_kind.py"
LOG="$OUTPUT_DIR/regen_progress.log"
PARTITIONS=8
FORCE="${FORCE:-0}"

[ -x "$BIN" ] || { echo "binary not built: $BIN (run cargo build first)" >&2; exit 1; }

if [ "$#" -gt 0 ]; then
  CONTIGS=("$@")
else
  CONTIGS=()
  while IFS= read -r name; do
    CONTIGS+=("$name")
  done < <(find "$ENTITY_DIR" -mindepth 1 -maxdepth 1 -name '*.lance' -type d \
    -exec basename {} .lance \; | sort)
fi

echo "=== regen start $(date) contigs=${#CONTIGS[@]} force=$FORCE ===" | tee -a "$LOG"
fail=0
for contig in "${CONTIGS[@]}"; do
  ds="$ENTITY_DIR/$contig.lance"
  if [ "$FORCE" != "1" ] && [ "$($KIND "$ds" 2>/dev/null || true)" = "new" ]; then
    echo "$contig SKIP already-new" | tee -a "$LOG"
    continue
  fi
  start=$(date +%s)
  if "$BIN" --cache-root "$CACHE_ROOT" --output-dir "$OUTPUT_DIR" \
       --chrom "$contig" --cache-source-type merged \
       --partitions "$PARTITIONS" --overwrite >>"$LOG" 2>&1; then
    elapsed=$(( $(date +%s) - start ))
    rows=$($KIND "$ds" >/dev/null 2>&1 && python -c "import lance;print(lance.dataset('$ds').count_rows())")
    echo "$contig OK rows=$rows elapsed=${elapsed}s" | tee -a "$LOG"
  else
    echo "$contig FAIL (see $LOG)" | tee -a "$LOG"
    fail=$((fail+1))
  fi
done
echo "=== regen done $(date) failures=$fail ===" | tee -a "$LOG"
[ "$fail" -eq 0 ]
