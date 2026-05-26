#!/usr/bin/env bash
# Generate a golden annotated VCF for a port case via real Ensembl VEP 115.
#
# Usage: gen_golden.sh <case>
#   Reads:  datafusion/bio-function-vep/tests/data/port/<case>/input.vcf
#   Writes: datafusion/bio-function-vep/tests/data/port/<case>/golden.vcf
#
# Per verification.md:
#   Layer 3: header embeds image digest, VEP version, cache hash,
#            generation timestamp, Perl-test prereq status (5 lines).
#   Layer 8: runs `prove t/<root>.t` inside the oracle image BEFORE
#            generating; aborts if the Perl test fails.
# Determinism: runs the annotation twice, byte-compares stripped output.
#
# IMPORTANT: this script mounts the NATIVE Ensembl cache (Storable/Sereal),
# NOT the vepyr parquet cache. Perl VEP cannot read parquet. The native
# cache slice lives at _cache115/native_cache/.
set -euo pipefail

if [ $# -ne 1 ]; then
  echo "Usage: $0 <case>" >&2
  exit 2
fi
CASE="$1"

REPO_ROOT="$(cd "$(dirname "$0")/../.." && pwd)"
CASE_DIR="$REPO_ROOT/datafusion/bio-function-vep/tests/data/port/$CASE"
CACHE_DIR="$REPO_ROOT/vep-benchmark/data/port/_cache115"
ENSEMBL_VEP_SRC="${ENSEMBL_VEP_SRC:-$REPO_ROOT/../ensembl-vep}"

[ -f "$CASE_DIR/input.vcf" ]                     || { echo "ERROR: missing $CASE_DIR/input.vcf" >&2; exit 1; }
[ -d "$CACHE_DIR/native_cache/homo_sapiens/115_GRCh38" ] || {
  echo "ERROR: missing native v115 cache slice at $CACHE_DIR/native_cache/homo_sapiens/115_GRCh38" >&2
  echo "  Run ./scripts/port/build_v115_parquet.sh first." >&2
  exit 1
}
[ -f "$CACHE_DIR/reference.fa" ]                 || { echo "ERROR: missing $CACHE_DIR/reference.fa" >&2; exit 1; }
[ -f "$CACHE_DIR/ORACLE_IMAGE" ]                 || { echo "ERROR: missing $CACHE_DIR/ORACLE_IMAGE — pin the oracle image first (Task 9)" >&2; exit 1; }
IMG=$(cat "$CACHE_DIR/ORACLE_IMAGE")

# Map case name to its Perl .t root for Layer 8 prereq.
# Most case names map 1:1; explicit overrides where they differ.
CASE_PERL_ROOT="$CASE"
case "$CASE" in
  cache_transcript)         CASE_PERL_ROOT="AnnotationSource_Cache_Transcript" ;;
  cache_variation)          CASE_PERL_ROOT="AnnotationSource_Cache_Variation"  ;;
  cache_regfeat|regfeat)    CASE_PERL_ROOT="AnnotationSource_Cache_RegFeat"    ;;
  parser_vcf*)              CASE_PERL_ROOT="Parser_VCF"                        ;;
  test_smoke)               CASE_PERL_ROOT=""                                  ;;
  *) echo "WARNING: no Perl-test mapping for case '$CASE'; skipping Layer 8 prereq" >&2 ;;
esac

C="vepgolden_$$"
docker rm -f "$C" >/dev/null 2>&1 || true
docker run -d --name "$C" --platform linux/amd64 --entrypoint sleep "$IMG" infinity >/dev/null
trap 'docker rm -f "$C" >/dev/null 2>&1 || true' EXIT

# Stage native cache, fasta, input, AND ensembl-vep src for Layer 8 prereq.
# Bind mounts are unavailable on this Mac (Docker Desktop sharing not enabled
# for /Users/wojtek/Documents), so use docker cp throughout.
docker cp "$CACHE_DIR/native_cache/."           "$C:/vep_cache/"
docker cp "$CACHE_DIR/reference.fa"             "$C:/work_reference.fa"
docker cp "$CACHE_DIR/reference.fa.fai"         "$C:/work_reference.fa.fai" 2>/dev/null || true
docker cp "$CASE_DIR/input.vcf"                 "$C:/work_input.vcf"
if [ -n "$CASE_PERL_ROOT" ] && [ -d "$ENSEMBL_VEP_SRC" ]; then
  docker cp "$ENSEMBL_VEP_SRC/."                "$C:/opt/vep_src/"
fi

# === Layer 8 prereq — Perl-run prerequisite ===
if [ -n "$CASE_PERL_ROOT" ] && [ -d "$ENSEMBL_VEP_SRC" ]; then
  echo "Layer 8 prereq: prove t/${CASE_PERL_ROOT}.t inside oracle image..."
  if ! docker exec -u root -w /opt/vep_src "$C" prove -v "t/${CASE_PERL_ROOT}.t" > /tmp/layer8_$$.log 2>&1; then
    echo "Layer 8 FAILED — Perl test t/${CASE_PERL_ROOT}.t does not pass in oracle image" >&2
    tail -50 /tmp/layer8_$$.log >&2
    exit 1
  fi
  echo "Layer 8 PASS"
  PERL_PASSED="true"
else
  PERL_PASSED="false"
fi

# === Generate golden VCF (twice, for determinism check) ===
generate() {
  local out="$1"
  docker exec -u root "$C" bash -c "
    set -e
    vep --offline --cache --dir_cache /vep_cache \\
        --species homo_sapiens --assembly GRCh38 --cache_version 115 \\
        --input_file /work_input.vcf --output_file '$out' \\
        --vcf --force_overwrite --no_stats \\
        --everything \\
        --fasta /work_reference.fa
  "
}

generate /work_golden_1.vcf
generate /work_golden_2.vcf

# Strip per-run header lines before comparing. Two known sources of nondeterminism
# in VEP's output:
#   - ##VEP=...time="YYYY-MM-DD HH:MM:SS"...   (wall-clock timestamp)
#   - ##VEP-command-line=...                    (echoes the output_file path)
# Strip both via the ^##VEP prefix. Use -u root for write access to /;
# default container user is non-root.
docker exec -u root "$C" bash -c '
  grep -vE "^##VEP" /work_golden_1.vcf > /work_golden_1_stripped.vcf
  grep -vE "^##VEP" /work_golden_2.vcf > /work_golden_2_stripped.vcf
  diff -q /work_golden_1_stripped.vcf /work_golden_2_stripped.vcf
' || { echo "Determinism check FAILED — two consecutive VEP runs produced different output" >&2; exit 1; }
echo "Determinism check PASS"

# === Embed Layer 3 header annotations ===
CACHE_HASH=$(find "$CACHE_DIR/native_cache" -type f -exec sha256sum {} \; 2>/dev/null | sort | sha256sum | awk '{print $1}')
TIMESTAMP=$(date -u +%Y-%m-%dT%H:%M:%SZ)

# Insert Layer 3 headers RIGHT AFTER ##fileformat= so the VCF stays
# spec-conformant (##fileformat= must be the first line). awk reads the
# original golden line-by-line and emits the 5 ##vepyr_port_oracle_*
# lines once, immediately after the first ##fileformat= line.
docker exec -u root "$C" bash -c "
  awk -v img='$IMG' -v hash='$CACHE_HASH' -v ts='$TIMESTAMP' -v perl='$PERL_PASSED' '
    /^##fileformat=/ && !done {
      print
      print \"##vepyr_port_oracle_image=\" img
      print \"##vepyr_port_oracle_vep_version=115\"
      print \"##vepyr_port_oracle_cache_hash=\" hash
      print \"##vepyr_port_oracle_generated_at=\" ts
      print \"##vepyr_port_oracle_perl_test_passed=\" perl
      done=1
      next
    }
    { print }
  ' /work_golden_1.vcf > /work_golden_final.vcf
"

docker cp "$C:/work_golden_final.vcf" "$CASE_DIR/golden.vcf"
echo "GOLDEN -> $CASE_DIR/golden.vcf"
echo "  oracle:     $IMG"
echo "  cache hash: $CACHE_HASH"
echo "  perl prereq: $PERL_PASSED"
