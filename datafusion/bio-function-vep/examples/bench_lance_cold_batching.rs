//! Compare Lance cold lookup scan stats for the same chr1 cold position set
//! split into many small batches vs fewer larger batches.

use std::collections::{HashMap, HashSet};
use std::env;
use std::fmt::Write as _;
use std::fs::File;
use std::io::{BufRead, BufReader};
use std::path::{Path, PathBuf};
use std::sync::{Arc, Mutex};
use std::time::{Duration, Instant};

use arrow_array::Array;
use datafusion_bio_function_vep::allele::{vcf_to_vep_allele, vcf_to_vep_input_allele};
use datafusion_bio_function_vep::annotate_provider::cache_lookup_column_names;
use datafusion_bio_function_vep::kv_cache::key_encoding::chrom_to_code;
use datafusion_bio_function_vep::kv_cache::position_index::{
    PositionIndex, find_position_index_file,
};
use datafusion_bio_function_vep::kv_cache::variant_bloom_index::{
    VariantBloomIndex, find_variant_bloom_index_file,
};
use datafusion_bio_function_vep::warm_cache::chrom_cache::projection_columns_for_cache;
use datafusion_bio_function_vep::warm_cache::key::{
    position_key_from_code, variant_key_from_position,
};
use flate2::read::MultiGzDecoder;
use futures::TryStreamExt;
use lance::Dataset;
use lance::dataset::scanner::{ExecutionSummaryCounts, MaterializationStyle};

type DynError = Box<dyn std::error::Error + Send + Sync>;
type BenchResult<T> = std::result::Result<T, DynError>;

#[derive(Debug)]
struct Workload {
    chrom: String,
    variants: usize,
    probe_attempts: usize,
    ordered_positions: Vec<u64>,
    variant_keys_by_position: HashMap<u64, Vec<u64>>,
    attempts_by_position: HashMap<u64, usize>,
}

#[derive(Debug)]
struct GateStats {
    warm_miss_unique: usize,
    position_index_unique: usize,
    bloom_unique: usize,
    bloom_attempts: usize,
    bloom_positions: Vec<u64>,
}

#[derive(Default, Debug, Clone)]
struct ScanIoStats {
    iops: usize,
    requests: usize,
    bytes_read: usize,
    indices_loaded: usize,
    parts_loaded: usize,
    index_comparisons: usize,
    fragments_scanned: usize,
    ranges_scanned: usize,
    rows_scanned: usize,
}

impl ScanIoStats {
    fn add_execution(&mut self, counts: &ExecutionSummaryCounts) {
        self.iops += counts.iops;
        self.requests += counts.requests;
        self.bytes_read += counts.bytes_read;
        self.indices_loaded += counts.indices_loaded;
        self.parts_loaded += counts.parts_loaded;
        self.index_comparisons += counts.index_comparisons;
        self.fragments_scanned += counts
            .all_counts
            .get("fragments_scanned")
            .copied()
            .unwrap_or_default();
        self.ranges_scanned += counts
            .all_counts
            .get("ranges_scanned")
            .copied()
            .unwrap_or_default();
        self.rows_scanned += counts
            .all_counts
            .get("rows_scanned")
            .copied()
            .unwrap_or_default();
    }
}

#[derive(Debug)]
struct BatchPlanResult {
    name: String,
    query_batch_size: usize,
    scans: usize,
    filter_keys: usize,
    rows: usize,
    result_batches: usize,
    selected_positions: usize,
    elapsed: Duration,
    io: ScanIoStats,
}

#[tokio::main(flavor = "multi_thread")]
async fn main() -> BenchResult<()> {
    let args = env::args().collect::<Vec<_>>();
    if args.len() < 4 || args.len() > 5 {
        eprintln!(
            "usage: {} <cache_root> <input.vcf.gz> <merged_lance_chr_dataset> [report.md]",
            args[0]
        );
        std::process::exit(2);
    }

    let cache_root = PathBuf::from(&args[1]);
    let vcf_path = PathBuf::from(&args[2]);
    let dataset_path = PathBuf::from(&args[3]);
    let report_path = args.get(4).map(PathBuf::from);
    let chrom = env::var("CHROM").unwrap_or_else(|_| "chr1".to_string());
    let cold_scan_batch_size = env_usize("COLD_SCAN_BATCH_SIZE", 2_000);
    let small_batch_size = env_usize("SMALL_QUERY_BATCH_SIZE", 238);
    let large_batch_size = env_usize("LARGE_QUERY_BATCH_SIZE", 2_000);
    let sort_positions = env_bool("SORT_COLD_POSITIONS", true);
    let materialization_style = MaterializationStyle::AllLate;

    println!(
        "cache_root={}\nvcf={}\ndataset={}\nchrom={chrom}\ncold_scan_batch_size={cold_scan_batch_size}\nsmall_query_batch_size={small_batch_size}\nlarge_query_batch_size={large_batch_size}\nsort_positions={sort_positions}",
        cache_root.display(),
        vcf_path.display(),
        dataset_path.display(),
    );

    let parse_started = Instant::now();
    let workload = load_vcf_workload(&vcf_path, &chrom)?;
    println!(
        "parsed variants={} probe_attempts={} unique_probe_positions={} parse_s={:.3}",
        workload.variants,
        workload.probe_attempts,
        workload.ordered_positions.len(),
        parse_started.elapsed().as_secs_f64()
    );

    let dataset = Dataset::open(dataset_path.to_string_lossy().as_ref()).await?;
    let warm_started = Instant::now();
    let warm_positions = load_warm_positions(&dataset).await?;
    println!(
        "warm_positions={} warm_position_scan_s={:.3}",
        warm_positions.len(),
        warm_started.elapsed().as_secs_f64()
    );

    let position_index_path =
        find_position_index_file(cache_root.join("variation.position_index"), &chrom)
            .ok_or_else(|| format!("missing position index for {chrom}"))?;
    let bloom_index_path =
        find_variant_bloom_index_file(cache_root.join("variation.variant_bloom_index"), &chrom)
            .ok_or_else(|| format!("missing variant bloom index for {chrom}"))?;
    let position_index = PositionIndex::read_from_path(&position_index_path)?;
    let bloom_index = VariantBloomIndex::read_from_path(&bloom_index_path)?;

    let mut gate = build_gate(
        &workload,
        &warm_positions,
        &position_index,
        &bloom_index,
        true,
    );
    if sort_positions {
        gate.bloom_positions.sort_unstable();
    }
    println!(
        "gate warm_miss_unique={} position_index_unique={} bloom_unique={} bloom_attempts={}",
        gate.warm_miss_unique, gate.position_index_unique, gate.bloom_unique, gate.bloom_attempts
    );

    let projection = projected_everything_columns(&dataset);
    println!(
        "projection_cols={} projection={}",
        projection.len(),
        projection.join(",")
    );

    let results = vec![
        scan_plan(
            &dataset,
            "small_chunks",
            &gate.bloom_positions,
            small_batch_size,
            cold_scan_batch_size,
            &projection,
            materialization_style.clone(),
        )
        .await?,
        scan_plan(
            &dataset,
            "large_chunks",
            &gate.bloom_positions,
            large_batch_size,
            cold_scan_batch_size,
            &projection,
            materialization_style,
        )
        .await?,
    ];

    let report = render_report(
        &cache_root,
        &vcf_path,
        &dataset_path,
        &workload,
        &gate,
        &projection,
        &results,
    );
    print!("{report}");
    if let Some(path) = report_path {
        if let Some(parent) = path.parent() {
            std::fs::create_dir_all(parent)?;
        }
        std::fs::write(&path, report)?;
        println!("wrote_report={}", path.display());
    }

    Ok(())
}

async fn load_warm_positions(dataset: &Dataset) -> BenchResult<HashSet<u64>> {
    let mut scanner = dataset.scan();
    scanner
        .filter("tier = 0")?
        .project(&["position_key"])?
        .batch_size(65_536)
        .materialization_style(MaterializationStyle::AllLate);
    let mut stream = scanner.try_into_stream().await?;
    let mut positions = HashSet::new();
    while let Some(batch) = stream.try_next().await? {
        let index = batch.schema().index_of("position_key")?;
        let array = batch
            .column(index)
            .as_any()
            .downcast_ref::<arrow_array::UInt64Array>()
            .ok_or("position_key must be UInt64")?;
        for row in 0..array.len() {
            if !array.is_null(row) {
                positions.insert(array.value(row));
            }
        }
    }
    Ok(positions)
}

async fn scan_plan(
    dataset: &Dataset,
    name: &str,
    positions: &[u64],
    query_batch_size: usize,
    scan_batch_size: usize,
    projection: &[String],
    materialization_style: MaterializationStyle,
) -> BenchResult<BatchPlanResult> {
    let started = Instant::now();
    let mut rows = 0usize;
    let mut result_batches = 0usize;
    let mut selected_positions = HashSet::new();
    let mut io = ScanIoStats::default();
    let mut scans = 0usize;
    let mut filter_keys = 0usize;
    let projection_refs = projection.iter().map(String::as_str).collect::<Vec<_>>();

    for chunk in positions.chunks(query_batch_size.max(1)) {
        scans += 1;
        filter_keys += chunk.len();
        let filter = position_filter(chunk);
        let stats_slot = Arc::new(Mutex::new(None::<ExecutionSummaryCounts>));
        let stats_sink = Arc::clone(&stats_slot);
        let mut scanner = dataset.scan();
        scanner
            .filter(&filter)?
            .project(&projection_refs)?
            .batch_size(scan_batch_size)
            .materialization_style(materialization_style.clone())
            .scan_stats_callback(Arc::new(move |stats| {
                *stats_sink.lock().expect("stats lock poisoned") = Some(stats.clone());
            }));
        let mut stream = scanner.try_into_stream().await?;
        while let Some(batch) = stream.try_next().await? {
            rows += batch.num_rows();
            result_batches += 1;
            append_position_keys(&batch, &mut selected_positions)?;
        }
        if let Some(stats) = stats_slot.lock().expect("stats lock poisoned").as_ref() {
            io.add_execution(stats);
        }
    }

    Ok(BatchPlanResult {
        name: name.to_string(),
        query_batch_size,
        scans,
        filter_keys,
        rows,
        result_batches,
        selected_positions: selected_positions.len(),
        elapsed: started.elapsed(),
        io,
    })
}

fn projected_everything_columns(dataset: &Dataset) -> Vec<String> {
    let available = dataset
        .schema()
        .fields
        .iter()
        .map(|field| field.name.clone())
        .collect::<HashSet<_>>();
    let mut preferred = vec!["consequence_types", "most_severe_consequence"];
    for column in cache_lookup_column_names() {
        if !preferred.contains(&column) {
            preferred.push(column);
        }
    }
    let cache_columns = preferred
        .into_iter()
        .filter(|column| available.contains(*column))
        .map(str::to_string)
        .collect::<Vec<_>>();
    projection_columns_for_cache(&cache_columns, true)
        .into_iter()
        .filter(|column| available.contains(column))
        .collect()
}

fn render_report(
    cache_root: &Path,
    vcf_path: &Path,
    dataset_path: &Path,
    workload: &Workload,
    gate: &GateStats,
    projection: &[String],
    results: &[BatchPlanResult],
) -> String {
    let mut out = String::new();
    let _ = writeln!(out, "# Lance Cold Batching Side Benchmark");
    let _ = writeln!(out);
    let _ = writeln!(out, "- Cache root: `{}`", cache_root.display());
    let _ = writeln!(out, "- VCF: `{}`", vcf_path.display());
    let _ = writeln!(out, "- Lance dataset: `{}`", dataset_path.display());
    let _ = writeln!(out, "- Chrom: `{}`", workload.chrom);
    let _ = writeln!(out, "- Parsed variants: `{}`", workload.variants);
    let _ = writeln!(out, "- Raw probe attempts: `{}`", workload.probe_attempts);
    let _ = writeln!(
        out,
        "- Bloom-admitted unique positions: `{}`",
        gate.bloom_unique
    );
    let _ = writeln!(out, "- Bloom-admitted attempts: `{}`", gate.bloom_attempts);
    let _ = writeln!(out, "- Projected columns: `{}`", projection.len());
    let _ = writeln!(out);
    let _ = writeln!(
        out,
        "| plan | query batch | scans | filter keys | rows | result batches | selected positions | seconds | bytes read | requests | ranges | fragments | parts loaded | index comparisons |"
    );
    let _ = writeln!(
        out,
        "|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|"
    );
    for result in results {
        let _ = writeln!(
            out,
            "| {} | {} | {} | {} | {} | {} | {} | {:.3} | {} | {} | {} | {} | {} | {} |",
            result.name,
            result.query_batch_size,
            result.scans,
            result.filter_keys,
            result.rows,
            result.result_batches,
            result.selected_positions,
            result.elapsed.as_secs_f64(),
            result.io.bytes_read,
            result.io.requests,
            result.io.ranges_scanned,
            result.io.fragments_scanned,
            result.io.parts_loaded,
            result.io.index_comparisons
        );
    }
    out
}

fn build_gate(
    workload: &Workload,
    warm_positions: &HashSet<u64>,
    position_index: &PositionIndex,
    bloom_index: &VariantBloomIndex,
    collect_colocated: bool,
) -> GateStats {
    let mut warm_miss_unique = 0usize;
    let mut position_index_unique = 0usize;
    let mut bloom_unique = 0usize;
    let mut bloom_attempts = 0usize;
    let mut bloom_positions = Vec::new();

    for &position_key in &workload.ordered_positions {
        if warm_positions.contains(&position_key) {
            continue;
        }
        warm_miss_unique += 1;
        if !position_index.contains_position_key(position_key) {
            continue;
        }
        position_index_unique += 1;
        let variant_keys = workload
            .variant_keys_by_position
            .get(&position_key)
            .map(Vec::as_slice)
            .unwrap_or(&[]);
        let variant_admitted = bloom_index.contains_any_variant_keys(variant_keys.iter().copied());
        let colocated_admitted = collect_colocated
            && if bloom_index.supports_position_fallback_keys() {
                bloom_index.contains_position_fallback_key(position_key)
            } else {
                true
            };
        if variant_admitted || colocated_admitted {
            bloom_unique += 1;
            bloom_attempts += workload
                .attempts_by_position
                .get(&position_key)
                .copied()
                .unwrap_or(1);
            bloom_positions.push(position_key);
        }
    }

    GateStats {
        warm_miss_unique,
        position_index_unique,
        bloom_unique,
        bloom_attempts,
        bloom_positions,
    }
}

fn load_vcf_workload(path: &Path, chrom: &str) -> BenchResult<Workload> {
    let file = File::open(path)?;
    let reader = BufReader::new(MultiGzDecoder::new(file));
    let chrom_bare = chrom.strip_prefix("chr").unwrap_or(chrom);
    let chrom_code = chrom_to_code(chrom_bare);
    let mut seen_target = false;
    let mut variants = 0usize;
    let mut probe_attempts = 0usize;
    let mut ordered_positions = Vec::new();
    let mut seen_positions = HashSet::new();
    let mut variant_keys_by_position = HashMap::<u64, Vec<u64>>::new();
    let mut attempts_by_position = HashMap::<u64, usize>::new();

    for line in reader.lines() {
        let line = line?;
        if line.starts_with('#') {
            continue;
        }
        let mut fields = line.split('\t');
        let Some(row_chrom) = fields.next() else {
            continue;
        };
        let row_chrom_bare = row_chrom.strip_prefix("chr").unwrap_or(row_chrom);
        if row_chrom_bare != chrom_bare {
            if seen_target {
                break;
            }
            continue;
        }
        seen_target = true;

        let Some(pos_raw) = fields.next() else {
            continue;
        };
        let _id = fields.next();
        let Some(ref_allele) = fields.next() else {
            continue;
        };
        let Some(alt_allele) = fields.next() else {
            continue;
        };
        let pos = pos_raw.parse::<i64>()?;
        let end = pos + ref_allele.len() as i64 - 1;

        for probe_start in build_probe_starts(pos, end, ref_allele, alt_allele, true) {
            let position_key = position_key_from_code(chrom_code, probe_start)?;
            probe_attempts += 1;
            *attempts_by_position.entry(position_key).or_default() += 1;
            if seen_positions.insert(position_key) {
                ordered_positions.push(position_key);
            }

            let keys = variant_keys_by_position.entry(position_key).or_default();
            for key in warm_variant_key_candidates(position_key, ref_allele, alt_allele) {
                push_unique(keys, key);
            }
        }
        variants += 1;
    }

    Ok(Workload {
        chrom: chrom.to_string(),
        variants,
        probe_attempts,
        ordered_positions,
        variant_keys_by_position,
        attempts_by_position,
    })
}

fn position_filter(keys: &[u64]) -> String {
    let mut filter = String::with_capacity(keys.len() * 20 + 40);
    filter.push_str("tier = 1 AND position_key IN (");
    for (idx, key) in keys.iter().enumerate() {
        if idx > 0 {
            filter.push(',');
        }
        let _ = write!(&mut filter, "{key}");
    }
    filter.push(')');
    filter
}

fn append_position_keys(
    batch: &arrow_array::RecordBatch,
    out: &mut HashSet<u64>,
) -> BenchResult<()> {
    let index = batch.schema().index_of("position_key")?;
    let array = batch
        .column(index)
        .as_any()
        .downcast_ref::<arrow_array::UInt64Array>()
        .ok_or("position_key must be UInt64")?;
    for row in 0..array.len() {
        if !array.is_null(row) {
            out.insert(array.value(row));
        }
    }
    Ok(())
}

fn warm_variant_key_candidates(position_key: u64, vcf_ref: &str, vcf_alt: &str) -> Vec<u64> {
    let mut keys = Vec::with_capacity(4);
    for alt in vcf_alt.split(['|', ',']).filter(|alt| !alt.is_empty()) {
        let (vep_ref, vep_alt) = vcf_to_vep_allele(vcf_ref, alt);
        push_unique(
            &mut keys,
            variant_key_from_position(position_key, &vep_ref, &vep_alt),
        );
        push_unique(
            &mut keys,
            variant_key_from_position(position_key, vcf_ref, alt),
        );
        push_unique(
            &mut keys,
            variant_key_from_position(position_key, vcf_ref, &vep_alt),
        );
        push_unique(
            &mut keys,
            variant_key_from_position(position_key, &vep_ref, alt),
        );
    }
    keys
}

fn build_probe_starts(
    norm_start_i64: i64,
    norm_end_i64: i64,
    vcf_ref: &str,
    vcf_alt: &str,
    extended_probes: bool,
) -> Vec<i64> {
    let mut probe_starts = Vec::with_capacity(6);
    push_unique(&mut probe_starts, norm_start_i64);

    if !extended_probes {
        return probe_starts;
    }

    if norm_start_i64 == norm_end_i64 {
        push_unique(&mut probe_starts, norm_start_i64.saturating_add(1));
    } else {
        push_unique(&mut probe_starts, norm_end_i64);
    }

    for alt in vcf_alt.split(['|', ',']).filter(|alt| !alt.is_empty()) {
        let (_, _, input_start) = vcf_to_vep_input_allele(norm_start_i64, vcf_ref, alt);
        push_unique(&mut probe_starts, input_start);

        let shift_usize = common_prefix_len(vcf_ref, alt);
        if shift_usize == 0 {
            continue;
        }
        if let Some(shifted_start) = norm_start_i64.checked_add(shift_usize as i64) {
            push_unique(&mut probe_starts, shifted_start);
        }
    }

    for alt in vcf_alt.split(['|', ',']).filter(|alt| !alt.is_empty()) {
        let (ref_event_len, alt_event_len) = canonical_event_lengths(vcf_ref, alt);
        if ref_event_len == 0 || alt_event_len != 0 {
            continue;
        }
        let del_len = ref_event_len as i64;
        let max_shift = del_len.min(32);
        for base_start in [norm_start_i64, norm_start_i64.saturating_sub(1)] {
            for shift in 0..=max_shift {
                let Some(candidate_start) = base_start.checked_add(shift) else {
                    continue;
                };
                let Some(candidate_end) = candidate_start.checked_add(del_len - 1) else {
                    continue;
                };
                if candidate_start > norm_end_i64 || candidate_end < norm_start_i64 {
                    continue;
                }
                push_unique(&mut probe_starts, candidate_start);
            }
        }
    }

    probe_starts
}

fn push_unique<T: Copy + PartialEq>(values: &mut Vec<T>, value: T) {
    if !values.contains(&value) {
        values.push(value);
    }
}

#[inline]
fn common_prefix_len(left: &str, right: &str) -> usize {
    left.as_bytes()
        .iter()
        .zip(right.as_bytes().iter())
        .take_while(|(a, b)| a == b)
        .count()
}

#[inline]
fn canonical_event_lengths(ref_allele: &str, alt_allele: &str) -> (usize, usize) {
    let ref_bytes = ref_allele.as_bytes();
    let alt_bytes = alt_allele.as_bytes();
    let mut ref_start = 0usize;
    let mut alt_start = 0usize;
    while ref_start < ref_bytes.len()
        && alt_start < alt_bytes.len()
        && ref_bytes[ref_start] == alt_bytes[alt_start]
    {
        ref_start += 1;
        alt_start += 1;
    }
    let mut ref_end = ref_bytes.len();
    let mut alt_end = alt_bytes.len();
    while ref_end > ref_start
        && alt_end > alt_start
        && ref_bytes[ref_end - 1] == alt_bytes[alt_end - 1]
    {
        ref_end -= 1;
        alt_end -= 1;
    }
    (ref_end - ref_start, alt_end - alt_start)
}

fn env_usize(name: &str, default: usize) -> usize {
    env::var(name)
        .ok()
        .and_then(|value| value.parse::<usize>().ok())
        .unwrap_or(default)
}

fn env_bool(name: &str, default: bool) -> bool {
    env::var(name)
        .ok()
        .map(|value| {
            matches!(
                value.to_ascii_lowercase().as_str(),
                "1" | "true" | "yes" | "on"
            )
        })
        .unwrap_or(default)
}
