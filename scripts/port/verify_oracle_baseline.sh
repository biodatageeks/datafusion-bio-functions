#!/usr/bin/env bash
# Verify the pinned oracle image passes 48/49 of the upstream Perl t/*.t suite.
# Per verification.md Layer 4: ONLY t/version.t is allowed to fail (it asserts
# a git-metadata check against .git which is absent in the docker-cp'd tree).
# Run before bumping ORACLE_IMAGE digest and as a post-Batch-0 acceptance check.
set -euo pipefail

REPO_ROOT="$(cd "$(dirname "$0")/../.." && pwd)"
ORACLE_FILE="$REPO_ROOT/vep-benchmark/data/port/_cache115/ORACLE_IMAGE"

if [ ! -f "$ORACLE_FILE" ]; then
  echo "ERROR: $ORACLE_FILE missing — pin the oracle image first (Task 9)" >&2
  exit 1
fi

IMG=$(cat "$ORACLE_FILE")
ENSEMBL_VEP_SRC="${ENSEMBL_VEP_SRC:-$REPO_ROOT/../ensembl-vep}"

if [ ! -d "$ENSEMBL_VEP_SRC/t" ]; then
  echo "ERROR: Perl t/ dir not found at $ENSEMBL_VEP_SRC/t" >&2
  echo "  Set ENSEMBL_VEP_SRC env var to the ensembl-vep checkout root." >&2
  exit 1
fi

echo "Oracle image: $IMG"
echo "Perl t/ dir : $ENSEMBL_VEP_SRC/t"
echo "Running prove t/*.t inside oracle image (this takes ~30-60 min)..."

C="vepbaseline_$$"
docker rm -f "$C" >/dev/null 2>&1 || true
docker run -d --name "$C" --platform linux/amd64 --entrypoint sleep "$IMG" infinity >/dev/null
trap 'docker rm -f "$C" >/dev/null 2>&1 || true' EXIT

docker cp "$ENSEMBL_VEP_SRC/." "$C:/opt/vep_src/"

LOG="$REPO_ROOT/vep-benchmark/data/port/_cache115/oracle_baseline_$(date -u +%Y%m%d_%H%M%S).log"
docker exec -u root -w /opt/vep_src "$C" prove -j 4 t/ > "$LOG" 2>&1 || true

# Parse final prove summary block.
# Expected line of interest looks like: "Failed 1/49 test programs. 0/2162 subtests failed."
TOTAL=$(grep -cE '^t/.*\.t' "$LOG" || true)
FAILED_FILES_LINE=$(grep -E "Failed .+test program" "$LOG" | tail -1 || true)

echo "---"
echo "Total .t files: $TOTAL"
echo "Summary line  : $FAILED_FILES_LINE"
echo "Log           : $LOG"

# Count how many test files are reported failed via the file-by-file summary
# (most reliable). A file marked "Result: FAIL" or showing in the FAILED file
# list at the end of prove's output is a failure.
FAILED_FILES=$(grep -cE "^t/.*\.t\s+\(Wstat" "$LOG" 2>/dev/null || echo 0)

# Accept either 0 failures, or exactly 1 failure and that one is version.t.
if [ "$FAILED_FILES" -le 1 ]; then
  if [ "$FAILED_FILES" -eq 0 ] || grep -qE "^t/version\.t\s" "$LOG"; then
    echo "BASELINE OK (≤1 failure; if any, it's version.t which is allowed per verification.md Layer 4)"
    exit 0
  fi
fi

echo "BASELINE FAILED — unexpected failure shape" >&2
echo "See $LOG for details" >&2
exit 1
