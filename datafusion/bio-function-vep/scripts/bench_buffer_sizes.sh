#!/usr/bin/env bash
set -euo pipefail

repo_root="$(git rev-parse --show-toplevel)"

input="${INPUT:-/Users/mwiewior/workspace/data_vepyr/HG002_GRCh38_1_22_v4.2.1_benchmark.vcf.gz}"
cache="${CACHE:-/Users/mwiewior/workspace/data_vepyr/115_GRCh38_merged}"
fasta="${FASTA:-/Users/mwiewior/workspace/data_vepyr/Homo_sapiens.GRCh38.dna.primary_assembly.fa}"
backend="${BACKEND:-fjall}"
forks="${FORKS:-8}"
contig_parallelism="${CONTIG_PARALLELISM:-1}"
buffer_sizes="${BUFFER_SIZES:-5000,10000,20000,50000}"
limit="${LIMIT:-none}"
sample_interval="${SAMPLE_INTERVAL:-1}"
chunked_buffer_lookup="${CHUNKED_BUFFER_LOOKUP:-0}"
output_dir="${OUTPUT_DIR:-/tmp/vep_buffer_size_sweep}"
write_outputs="${WRITE_OUTPUTS:-0}"
build="${BUILD:-0}"

binary="${BINARY:-${repo_root}/target/release/examples/bench_annotate_vcf}"

if [[ "${build}" == "1" ]]; then
  (
    cd "${repo_root}"
    RUSTFLAGS="${RUSTFLAGS:--C target-cpu=native}" \
      cargo build --release -p datafusion-bio-function-vep --example bench_annotate_vcf --features kv-cache
  )
fi

if [[ ! -x "${binary}" ]]; then
  echo "benchmark binary not found: ${binary}" >&2
  echo "run with BUILD=1 or build it first with native release flags" >&2
  exit 1
fi

mkdir -p "${output_dir}"
results="${output_dir}/buffer_size_sweep.csv"
echo "buffer_size,forks,contig_parallelism,limit,chunked_buffer_lookup,rows,elapsed_s,rows_per_s,avg_cpu_pct,avg_rss_mb,max_rss_mb,log_path" > "${results}"

now_secs() {
  perl -MTime::HiRes=time -e 'printf "%.6f\n", time'
}

for buffer_size in ${buffer_sizes//,/ }; do
  log_path="${output_dir}/buffer_${buffer_size}.log"
  samples_path="${output_dir}/buffer_${buffer_size}.samples"
  rm -f "${log_path}" "${samples_path}"

  if [[ "${write_outputs}" == "1" ]]; then
    output_path="${output_dir}/buffer_${buffer_size}.vcf"
    rm -f "${output_path}"
  else
    output_path="/dev/null"
  fi

  cmd=(
    "${binary}"
    --input "${input}"
    --cache "${cache}"
    --output "${output_path}"
    --backend "${backend}"
    --everything
    --extended-probes
    --reference-fasta "${fasta}"
    --forks "${forks}"
    --contig-parallelism "${contig_parallelism}"
    --buffer-size "${buffer_size}"
    --compression none
    --no-progress
  )
  if [[ "${chunked_buffer_lookup}" == "1" || "${chunked_buffer_lookup}" == "true" ]]; then
    cmd+=(--chunked-buffer-lookup)
  fi
  if [[ "${limit}" != "none" ]]; then
    cmd+=(--limit "${limit}")
  fi

  echo "running buffer_size=${buffer_size} forks=${forks} contig_parallelism=${contig_parallelism} limit=${limit} chunked_buffer_lookup=${chunked_buffer_lookup}" >&2
  start="$(now_secs)"
  "${cmd[@]}" >"${log_path}" 2>&1 &
  pid=$!

  while kill -0 "${pid}" 2>/dev/null; do
    ps -p "${pid}" -o %cpu= -o rss= | awk 'NF == 2 { print $1, $2 }' >> "${samples_path}" || true
    sleep "${sample_interval}"
  done

  status=0
  wait "${pid}" || status=$?
  end="$(now_secs)"
  if [[ "${status}" != "0" ]]; then
    echo "run failed for buffer_size=${buffer_size}; see ${log_path}" >&2
    exit "${status}"
  fi

  elapsed="$(awk -v s="${start}" -v e="${end}" 'BEGIN { printf "%.3f", e - s }')"
  rows="$(awk '/rows:/ { value=$2 } END { gsub(",", "", value); print value + 0 }' "${log_path}")"
  tool_elapsed="$(awk '/time:/ { value=$2 } END { sub(/s$/, "", value); print value + 0 }' "${log_path}")"
  rows_per_s="$(awk -v r="${rows}" -v t="${tool_elapsed}" 'BEGIN { if (t > 0) printf "%.1f", r / t; else print "0.0" }')"
  read -r avg_cpu avg_rss max_rss < <(
    awk '
      NF == 2 {
        cpu += $1;
        rss += $2;
        if ($2 > max_rss) max_rss = $2;
        n++;
      }
      END {
        if (n == 0) {
          printf "0.0 0.0 0.0\n";
        } else {
          printf "%.1f %.1f %.1f\n", cpu / n, (rss / n) / 1024.0, max_rss / 1024.0;
        }
      }
    ' "${samples_path}"
  )

  echo "${buffer_size},${forks},${contig_parallelism},${limit},${chunked_buffer_lookup},${rows},${elapsed},${rows_per_s},${avg_cpu},${avg_rss},${max_rss},${log_path}" >> "${results}"
done

echo "${results}"
