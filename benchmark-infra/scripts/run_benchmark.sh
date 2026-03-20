#!/usr/bin/env bash
set -euo pipefail

# Usage: ./run_benchmark.sh [ensembl|parquet|fjall|all] [chr1|chr22|wgs] [buffer_size] [fork]
#
# Backend:
#   ensembl  — Ensembl VEP via Docker only
#   parquet  — Native VEP with Parquet backend only
#   fjall    — Native VEP with Fjall backend only
#   all      — All three (default)
#
# Dataset:
#   chr1     — HG002 chr1 only, 319K variants (default)
#   chr22    — HG002 chr22 only, ~50K variants (quick test)
#   wgs      — HG002 full genome (chr1-22), ~4M variants
#
# Buffer size:
#   Ensembl VEP --buffer_size (default: 5000)
#
# Fork:
#   Ensembl VEP --fork (default: nproc, all available cores)

# ── Configuration ─────────────────────────────────────────────────────
DATA_DIR="/data/vep"
REPO_DIR="$HOME/datafusion-bio-functions"
RESULTS_DIR="${DATA_DIR}/results"
FASTA="${DATA_DIR}/Homo_sapiens.GRCh38.dna.primary_assembly.fa"
VEP_IMAGE="ensemblorg/ensembl-vep:release_115.2"

MODE="${1:-all}"
DATASET="${2:-chr1}"
BUFFER_SIZE="${3:-5000}"
FORK="${4:-$(nproc)}"

# ── Dataset selection ─────────────────────────────────────────────────
case "${DATASET}" in
  chr1)
    VCF="${DATA_DIR}/vcf/HG002_chr1.vcf.gz"
    CACHE_DIR="${DATA_DIR}/chr1-vep"
    FJALL_CACHE="${CACHE_DIR}/fjall_variation_cache"
    SAMPLE_LIMIT=320000
    ;;
  chr22)
    VCF="${DATA_DIR}/vcf/HG002_chr22.vcf.gz"
    CACHE_DIR="${DATA_DIR}/chr22-vep"
    FJALL_CACHE="${CACHE_DIR}/fjall_variation_cache"
    SAMPLE_LIMIT=0
    ;;
  wgs)
    VCF="${DATA_DIR}/vcf/HG002_GRCh38_1_22_v4.2.1_benchmark.vcf.gz"
    CACHE_DIR="${DATA_DIR}/wgs-vep"
    FJALL_CACHE="${CACHE_DIR}/fjall_variation_cache"
    SAMPLE_LIMIT=0
    ;;
  *)
    echo "Unknown dataset: ${DATASET}. Use chr1, chr22, or wgs."
    exit 1
    ;;
esac

if [ ! -f "${VCF}" ]; then
  echo "ERROR: VCF not found: ${VCF}"
  exit 1
fi

ENSEMBL_OUT="${RESULTS_DIR}/ensembl_${DATASET}.vcf"

TIMESTAMP=$(date +%Y%m%d_%H%M%S)
RESULT_FILE="${RESULTS_DIR}/benchmark_${MODE}_${DATASET}_${TIMESTAMP}.txt"

mkdir -p "${RESULTS_DIR}"

source "$HOME/.cargo/env"

echo "=== VEP Benchmark on $(hostname) [mode=${MODE}, dataset=${DATASET}] ===" | tee "${RESULT_FILE}"
echo "Date: $(date -u +%Y-%m-%dT%H:%M:%SZ)" | tee -a "${RESULT_FILE}"
echo "Machine: $(uname -a)" | tee -a "${RESULT_FILE}"
echo "CPUs: $(nproc)" | tee -a "${RESULT_FILE}"
echo "RAM: $(free -h | awk '/Mem:/{print $2}')" | tee -a "${RESULT_FILE}"
echo "VCF: ${VCF}" | tee -a "${RESULT_FILE}"
echo "Cache: ${CACHE_DIR}" | tee -a "${RESULT_FILE}"
echo "Sample limit: ${SAMPLE_LIMIT} (0=all)" | tee -a "${RESULT_FILE}"
echo "Ensembl buffer_size: ${BUFFER_SIZE}" | tee -a "${RESULT_FILE}"
echo "Ensembl fork: ${FORK}" | tee -a "${RESULT_FILE}"
echo "" | tee -a "${RESULT_FILE}"

# ── Ensembl VEP (Docker) ────────────────────────────────────────────
if [ "$MODE" = "ensembl" ] || [ "$MODE" = "all" ]; then
  if [ -f "${ENSEMBL_OUT}" ]; then
    echo "=== Ensembl VEP: SKIPPED (output exists: ${ENSEMBL_OUT}) ===" | tee -a "${RESULT_FILE}"
  else
    echo "=== Running Ensembl VEP (Docker) [${DATASET}] ===" | tee -a "${RESULT_FILE}"

    ENSEMBL_START=$(date +%s%N)

    docker run --rm \
      --user "$(id -u):$(id -g)" \
      --tmpfs /tmp:exec \
      -v "${DATA_DIR}/homo_sapiens:/opt/vep/.vep/homo_sapiens:ro" \
      -v "${DATA_DIR}/vcf:/work:ro" \
      -v "$(dirname "${FASTA}"):/fasta" \
      -v "${RESULTS_DIR}:/output" \
      "${VEP_IMAGE}" \
      vep \
        --input_file "/work/$(basename "${VCF}")" \
        --output_file "/output/ensembl_${DATASET}.vcf" \
        --vcf \
        --offline \
        --cache \
        --dir /opt/vep/.vep \
        --assembly GRCh38 \
        --everything \
        --fasta "/fasta/$(basename "${FASTA}")" \
        --fork "${FORK}" \
        --buffer_size "${BUFFER_SIZE}" \
        --no_check_variants_order \
        --force_overwrite \
        --no_stats

    ENSEMBL_END=$(date +%s%N)
    ENSEMBL_MS=$(( (ENSEMBL_END - ENSEMBL_START) / 1000000 ))

    echo "Ensembl VEP time: ${ENSEMBL_MS}ms ($(echo "scale=1; ${ENSEMBL_MS}/1000" | bc)s)" | tee -a "${RESULT_FILE}"
  fi
  echo "" | tee -a "${RESULT_FILE}"
fi

# ── Native VEP — Parquet backend ────────────────────────────────────
if [ "$MODE" = "parquet" ] || [ "$MODE" = "all" ]; then
  if [ -d "${CACHE_DIR}" ]; then
    echo "=== Running native VEP (Parquet) [${DATASET}] ===" | tee -a "${RESULT_FILE}"

    PARQUET_OUT="${RESULTS_DIR}/native_parquet_${DATASET}_${TIMESTAMP}.vcf"

    VEP_PROFILE=1 "${REPO_DIR}/target/release/examples/profile_annotation" \
      "${VCF}" \
      "${CACHE_DIR}" \
      "${SAMPLE_LIMIT}" \
      --everything \
      --reference-fasta-path="${FASTA}" \
      --output="${PARQUET_OUT}" \
      2>&1 | tee -a "${RESULT_FILE}"

    echo "" | tee -a "${RESULT_FILE}"
  else
    echo "=== Parquet backend: SKIPPED (no cache at ${CACHE_DIR}) ===" | tee -a "${RESULT_FILE}"
    echo "" | tee -a "${RESULT_FILE}"
  fi
fi

# ── Native VEP — Fjall backend ──────────────────────────────────────
if [ "$MODE" = "fjall" ] || [ "$MODE" = "all" ]; then
  if [ -d "${FJALL_CACHE}/keyspaces" ]; then
    echo "=== Running native VEP (Fjall) [${DATASET}] ===" | tee -a "${RESULT_FILE}"

    FJALL_OUT="${RESULTS_DIR}/native_fjall_${DATASET}_${TIMESTAMP}.vcf"

    VEP_PROFILE=1 "${REPO_DIR}/target/release/examples/profile_annotation" \
      "${VCF}" \
      "${FJALL_CACHE}" \
      "${SAMPLE_LIMIT}" \
      --everything \
      --reference-fasta-path="${FASTA}" \
      --output="${FJALL_OUT}" \
      2>&1 | tee -a "${RESULT_FILE}"

    echo "" | tee -a "${RESULT_FILE}"
  else
    echo "=== Fjall backend: SKIPPED (no cache at ${FJALL_CACHE}) ===" | tee -a "${RESULT_FILE}"
    echo "Build it with:" | tee -a "${RESULT_FILE}"
    echo "  ${REPO_DIR}/target/release/examples/load_cache_full <variation_parquet> ${FJALL_CACHE} 4" | tee -a "${RESULT_FILE}"
    echo "" | tee -a "${RESULT_FILE}"
  fi
fi

# ── Summary ──────────────────────────────────────────────────────────
echo "=== Summary ===" | tee -a "${RESULT_FILE}"

for f in "${ENSEMBL_OUT}" "${PARQUET_OUT:-}" "${FJALL_OUT:-}"; do
  if [ -n "$f" ] && [ -f "$f" ]; then
    LINES=$(grep -cv '^#' "$f" || true)
    echo "$(basename "$f"): ${LINES} variants" | tee -a "${RESULT_FILE}"
  fi
done

echo ""
echo "Done! Results in ${RESULT_FILE}"
