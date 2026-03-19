#!/usr/bin/env bash
set -euo pipefail

# ── Configuration ─────────────────────────────────────────────────────
DATA_DIR="/data/vep"
REPO_DIR="$HOME/datafusion-bio-functions"
RESULTS_DIR="${DATA_DIR}/results"
VCF="${DATA_DIR}/vcf/HG002_chr1.vcf.gz"
FASTA="${DATA_DIR}/Homo_sapiens.GRCh38.dna.primary_assembly.fa"
CHR1_CACHE="${DATA_DIR}/chr1-vep"
ENSEMBL_CACHE="${DATA_DIR}/homo_sapiens/115_GRCh38"
VARIATION_PARQUET="${CHR1_CACHE}/115_GRCh38_variation_1_vep.parquet"
VEP_IMAGE="ensemblorg/ensembl-vep:release_115.2"

TIMESTAMP=$(date +%Y%m%d_%H%M%S)
RESULT_FILE="${RESULTS_DIR}/benchmark_${TIMESTAMP}.txt"

mkdir -p "${RESULTS_DIR}"

echo "=== VEP Benchmark on $(hostname) ===" | tee "${RESULT_FILE}"
echo "Date: $(date -u +%Y-%m-%dT%H:%M:%SZ)" | tee -a "${RESULT_FILE}"
echo "Machine: $(uname -a)" | tee -a "${RESULT_FILE}"
echo "CPUs: $(nproc)" | tee -a "${RESULT_FILE}"
echo "RAM: $(free -h | awk '/Mem:/{print $2}')" | tee -a "${RESULT_FILE}"
echo "" | tee -a "${RESULT_FILE}"

# ── 1. Ensembl VEP (Docker) ──────────────────────────────────────────
echo "=== Running Ensembl VEP (Docker) ===" | tee -a "${RESULT_FILE}"

ENSEMBL_OUT="${RESULTS_DIR}/ensembl_chr1_${TIMESTAMP}.vcf"
ENSEMBL_START=$(date +%s%N)

docker run --rm \
  --user "$(id -u):$(id -g)" \
  -v "${DATA_DIR}:/opt/vep/.vep" \
  -v "${DATA_DIR}/vcf:/work:ro" \
  -v "$(dirname "${FASTA}"):/fasta" \
  -v "${RESULTS_DIR}:/output" \
  "${VEP_IMAGE}" \
  vep \
    --input_file /work/HG002_chr1.vcf.gz \
    --output_file /output/ensembl_output.vcf \
    --vcf \
    --offline \
    --cache \
    --dir /opt/vep/.vep \
    --assembly GRCh38 \
    --everything \
    --fasta /fasta/$(basename "${FASTA}") \
    --force_overwrite \
    --no_stats

ENSEMBL_END=$(date +%s%N)
ENSEMBL_MS=$(( (ENSEMBL_END - ENSEMBL_START) / 1000000 ))

echo "Ensembl VEP time: ${ENSEMBL_MS}ms ($(echo "scale=1; ${ENSEMBL_MS}/1000" | bc)s)" | tee -a "${RESULT_FILE}"
echo "" | tee -a "${RESULT_FILE}"

# Move output
if [ -f "${RESULTS_DIR}/ensembl_output.vcf" ]; then
  mv "${RESULTS_DIR}/ensembl_output.vcf" "${ENSEMBL_OUT}"
  echo "Ensembl output: ${ENSEMBL_OUT}" | tee -a "${RESULT_FILE}"
fi

# ── 2. Native datafusion-bio VEP ─────────────────────────────────────
echo "=== Running native datafusion-bio VEP ===" | tee -a "${RESULT_FILE}"

NATIVE_OUT="${RESULTS_DIR}/native_chr1_${TIMESTAMP}.vcf"

source "$HOME/.cargo/env"

VEP_PROFILE=1 cargo run \
  -p datafusion-bio-function-vep \
  --release \
  --manifest-path="${REPO_DIR}/Cargo.toml" \
  --example profile_annotation -- \
  "${VCF}" \
  "${CHR1_CACHE}" \
  320000 \
  --everything \
  --reference-fasta-path="${FASTA}" \
  --output="${NATIVE_OUT}" \
  2>&1 | tee -a "${RESULT_FILE}"

echo "" | tee -a "${RESULT_FILE}"

# ── 3. Summary ───────────────────────────────────────────────────────
echo "=== Summary ===" | tee -a "${RESULT_FILE}"
echo "Ensembl VEP (Docker): ${ENSEMBL_MS}ms" | tee -a "${RESULT_FILE}"
echo "Results saved to: ${RESULT_FILE}" | tee -a "${RESULT_FILE}"

# Count output lines
if [ -f "${ENSEMBL_OUT}" ]; then
  ENSEMBL_LINES=$(grep -cv '^#' "${ENSEMBL_OUT}" || true)
  echo "Ensembl output variants: ${ENSEMBL_LINES}" | tee -a "${RESULT_FILE}"
fi
if [ -f "${NATIVE_OUT}" ]; then
  NATIVE_LINES=$(grep -cv '^#' "${NATIVE_OUT}" || true)
  echo "Native output variants: ${NATIVE_LINES}" | tee -a "${RESULT_FILE}"
fi

echo ""
echo "Done! Results in ${RESULT_FILE}"
