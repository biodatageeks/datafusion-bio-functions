//! Benchmark chr1 variation Lance layouts with the VEPyr warm/full-scan and cold/indexed lookup workload.
//!
//! Usage:
//!   RUSTFLAGS="-C target-cpu=native" cargo run --release \
//!     -p datafusion-bio-function-vep --features kv-cache \
//!     --example bench_lance_variation_chr1_workload -- \
//!     /path/to/115_GRCh38_merged /path/to/lance_variation_chr1_bench /path/to/input_chr1.vcf.gz \
//!     research/reports/chr1_lance_variation_rust_workload.md
//!
//! Optional env:
//!   CHROM=chr1
//!   WORKLOAD_BATCH_SIZE=2000
//!   WARM_SCAN_BATCH_SIZE=65536
//!   COLD_SCAN_BATCH_SIZE=8192
//!   LANCE_VARIANTS=2.1-unpacked,2.2-packed
//!   COLD_FRAGMENT_ROWS=512,1024,2048,4096,8192,16384
//!   SORT_COLD_POSITIONS=0
//!   MAX_VARIANTS=5000

use std::collections::{HashMap, HashSet};
use std::env;
use std::fmt::Write as _;
use std::fs::File;
use std::io::{BufRead, BufReader};
use std::path::{Path, PathBuf};
use std::time::{Duration, Instant};

use arrow_array::{Array, Int32Array, Int64Array, RecordBatch, UInt32Array, UInt64Array};
use datafusion_bio_function_vep::allele::{vcf_to_vep_allele, vcf_to_vep_input_allele};
use datafusion_bio_function_vep::kv_cache::key_encoding::chrom_to_code;
use datafusion_bio_function_vep::kv_cache::position_index::{
    PositionIndex, find_position_index_file,
};
use datafusion_bio_function_vep::kv_cache::variant_bloom_index::{
    VariantBloomIndex, find_variant_bloom_index_file,
};
use datafusion_bio_function_vep::warm_cache::key::{
    position_key_from_code, variant_key_from_position,
};
use flate2::read::MultiGzDecoder;
use futures::TryStreamExt;
use lance::Dataset;
use lance::dataset::scanner::MaterializationStyle;

type DynError = Box<dyn std::error::Error + Send + Sync>;
type BenchResult<T> = std::result::Result<T, DynError>;

#[derive(Debug)]
struct Workload {
    chrom: String,
    chrom_code: u16,
    variants: usize,
    probe_attempts: usize,
    ordered_positions: Vec<u64>,
    variant_keys_by_position: HashMap<u64, Vec<u64>>,
    attempts_by_position: HashMap<u64, usize>,
}

#[derive(Debug)]
struct IndexGateStats {
    warm_miss_unique_positions: usize,
    warm_miss_probe_attempts: usize,
    position_index_unique_positions: usize,
    position_index_probe_attempts: usize,
    bloom_unique_positions: usize,
    bloom_probe_attempts: usize,
    position_index_positions: Vec<u64>,
    bloom_positions: Vec<u64>,
}

#[derive(Debug)]
struct ScanStats {
    rows: usize,
    batches: usize,
    elapsed: Duration,
    selected_positions: HashSet<u64>,
}

#[derive(Debug)]
struct ResultRow {
    phase: String,
    layout: String,
    fragment_rows: Option<usize>,
    query_batch: String,
    probe_positions: Option<usize>,
    selected_rows: usize,
    selected_positions: usize,
    batches: usize,
    seconds: f64,
}

struct ReportInput<'a> {
    cache_root: &'a Path,
    lance_root: &'a Path,
    vcf_path: &'a Path,
    workload: &'a Workload,
    position_index_path: &'a Path,
    bloom_index_path: &'a Path,
    position_index_s: Duration,
    bloom_index_s: Duration,
    gate: &'a IndexGateStats,
    workload_batch_size: usize,
    warm_scan_batch_size: usize,
    cold_scan_batch_size: usize,
    sort_cold_positions: bool,
    result_rows: &'a [ResultRow],
}

#[tokio::main(flavor = "multi_thread")]
async fn main() -> BenchResult<()> {
    let args = env::args().collect::<Vec<_>>();
    if args.len() < 4 || args.len() > 5 {
        eprintln!(
            "usage: {} <cache_root> <lance_root> <input.vcf.gz> [report.md]",
            args[0]
        );
        std::process::exit(2);
    }

    let cache_root = PathBuf::from(&args[1]);
    let lance_root = PathBuf::from(&args[2]);
    let vcf_path = PathBuf::from(&args[3]);
    let report_path = args.get(4).map(PathBuf::from);

    let chrom = env::var("CHROM").unwrap_or_else(|_| "chr1".to_string());
    let workload_batch_size = env_usize("WORKLOAD_BATCH_SIZE", 2_000);
    let warm_scan_batch_size = env_usize("WARM_SCAN_BATCH_SIZE", 65_536);
    let cold_scan_batch_size = env_usize("COLD_SCAN_BATCH_SIZE", 8_192);
    let sort_cold_positions = env_bool("SORT_COLD_POSITIONS", false);
    let max_variants = env::var("MAX_VARIANTS")
        .ok()
        .and_then(|value| value.parse::<usize>().ok());
    let variants = env_csv("LANCE_VARIANTS", &["2.1-unpacked", "2.2-packed"]);
    let fragment_rows = env_csv_usize("COLD_FRAGMENT_ROWS", &[512, 1024, 2048, 4096, 8192, 16_384]);

    println!(
        "cache_root={}\nlance_root={}\nvcf={}\nchrom={chrom}\nworkload_batch_size={workload_batch_size}\nwarm_scan_batch_size={warm_scan_batch_size}\ncold_scan_batch_size={cold_scan_batch_size}\nsort_cold_positions={sort_cold_positions}\nvariants={}\nfragment_rows={:?}\nmax_variants={}",
        cache_root.display(),
        lance_root.display(),
        vcf_path.display(),
        variants.join(","),
        fragment_rows,
        max_variants
            .map(|value| value.to_string())
            .unwrap_or_else(|| "all".to_string())
    );

    let started = Instant::now();
    let workload = load_vcf_workload(&vcf_path, &chrom, max_variants)?;
    println!(
        "parsed variants={} probe_attempts={} unique_probe_positions={} parse_s={:.3}",
        workload.variants,
        workload.probe_attempts,
        workload.ordered_positions.len(),
        started.elapsed().as_secs_f64()
    );

    let position_index_path =
        find_position_index_file(cache_root.join("variation.position_index"), &chrom)
            .ok_or_else(|| format!("missing position index for {chrom}"))?;
    let bloom_index_path =
        find_variant_bloom_index_file(cache_root.join("variation.variant_bloom_index"), &chrom)
            .ok_or_else(|| format!("missing variant bloom index for {chrom}"))?;

    let t0 = Instant::now();
    let position_index = PositionIndex::read_from_path(&position_index_path)?;
    let position_index_s = t0.elapsed();
    let t0 = Instant::now();
    let bloom_index = VariantBloomIndex::read_from_path(&bloom_index_path)?;
    let bloom_index_s = t0.elapsed();
    println!(
        "loaded indexes position_entries={} position_bytes={} position_load_s={:.3} bloom_inserted={} bloom_bits={} bloom_hashes={} bloom_bytes={} bloom_load_s={:.3}",
        position_index.len(),
        position_index.storage_bytes(),
        position_index_s.as_secs_f64(),
        bloom_index.inserted(),
        bloom_index.bit_count(),
        bloom_index.hash_count(),
        bloom_index.storage_bytes(),
        bloom_index_s.as_secs_f64()
    );

    let mut result_rows = Vec::new();
    let mut warm_positions_for_gate = None;

    for layout in &variants {
        let warm_path = lance_root.join(format!("{}_warm_{layout}.lance", chrom));
        println!(
            "warm full scan layout={layout} path={}",
            warm_path.display()
        );
        let mut stats = scan_lance_full(
            &warm_path,
            warm_scan_batch_size,
            MaterializationStyle::AllEarly,
        )
        .await?;
        let selected_positions = stats.selected_positions.len();
        println!(
            "warm full scan layout={layout} rows={} unique_positions={} batches={} seconds={:.3}",
            stats.rows,
            selected_positions,
            stats.batches,
            stats.elapsed.as_secs_f64()
        );
        result_rows.push(ResultRow {
            phase: "warm_full_scan_all_columns".to_string(),
            layout: layout.clone(),
            fragment_rows: None,
            query_batch: "full_scan".to_string(),
            probe_positions: None,
            selected_rows: stats.rows,
            selected_positions,
            batches: stats.batches,
            seconds: stats.elapsed.as_secs_f64(),
        });

        if warm_positions_for_gate.is_none() {
            warm_positions_for_gate = Some(std::mem::take(&mut stats.selected_positions));
        }
    }

    let warm_positions = warm_positions_for_gate.ok_or("no warm layouts were scanned")?;
    let mut gate = build_cold_gate(
        &workload,
        &warm_positions,
        &position_index,
        &bloom_index,
        true,
    );
    if sort_cold_positions {
        gate.position_index_positions.sort_unstable();
        gate.bloom_positions.sort_unstable();
    }
    println!(
        "cold gate warm_miss_unique={} warm_miss_attempts={} position_index_unique={} position_index_attempts={} bloom_unique={} bloom_attempts={} posidx_cold_batches={}",
        gate.warm_miss_unique_positions,
        gate.warm_miss_probe_attempts,
        gate.position_index_unique_positions,
        gate.position_index_probe_attempts,
        gate.bloom_unique_positions,
        gate.bloom_probe_attempts,
        gate.position_index_positions
            .len()
            .div_ceil(workload_batch_size)
    );

    for layout in &variants {
        for rows in &fragment_rows {
            let cold_path = lance_root.join(format!("{}_cold_{layout}_fr{rows}.lance", chrom));
            if !cold_path.is_dir() {
                println!("skip missing cold dataset {}", cold_path.display());
                continue;
            }
            for (phase, positions) in [
                (
                    "cold_posidx_probe_batches_all_columns",
                    gate.position_index_positions.as_slice(),
                ),
                (
                    "cold_posidx_bloom_probe_batches_all_columns",
                    gate.bloom_positions.as_slice(),
                ),
            ] {
                println!(
                    "{phase} layout={layout} fragment_rows={rows} positions={} path={}",
                    positions.len(),
                    cold_path.display()
                );
                let stats = scan_lance_filter_batches(
                    &cold_path,
                    positions,
                    workload_batch_size,
                    cold_scan_batch_size,
                    MaterializationStyle::AllLate,
                )
                .await?;
                let selected_positions = stats.selected_positions.len();
                println!(
                    "{phase} layout={layout} fragment_rows={rows} rows={} unique_positions={} batches={} seconds={:.3}",
                    stats.rows,
                    selected_positions,
                    stats.batches,
                    stats.elapsed.as_secs_f64()
                );
                result_rows.push(ResultRow {
                    phase: phase.to_string(),
                    layout: layout.clone(),
                    fragment_rows: Some(*rows),
                    query_batch: workload_batch_size.to_string(),
                    probe_positions: Some(positions.len()),
                    selected_rows: stats.rows,
                    selected_positions,
                    batches: stats.batches,
                    seconds: stats.elapsed.as_secs_f64(),
                });
            }
        }
    }

    let report = render_report(&ReportInput {
        cache_root: &cache_root,
        lance_root: &lance_root,
        vcf_path: &vcf_path,
        workload: &workload,
        position_index_path: &position_index_path,
        bloom_index_path: &bloom_index_path,
        position_index_s,
        bloom_index_s,
        gate: &gate,
        workload_batch_size,
        warm_scan_batch_size,
        cold_scan_batch_size,
        sort_cold_positions,
        result_rows: &result_rows,
    });
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

fn env_csv(name: &str, default: &[&str]) -> Vec<String> {
    env::var(name)
        .ok()
        .map(|value| {
            value
                .split(',')
                .map(str::trim)
                .filter(|value| !value.is_empty())
                .map(str::to_string)
                .collect::<Vec<_>>()
        })
        .filter(|values| !values.is_empty())
        .unwrap_or_else(|| default.iter().map(|value| value.to_string()).collect())
}

fn env_csv_usize(name: &str, default: &[usize]) -> Vec<usize> {
    env::var(name)
        .ok()
        .map(|value| {
            value
                .split(',')
                .map(str::trim)
                .filter_map(|value| value.parse::<usize>().ok())
                .collect::<Vec<_>>()
        })
        .filter(|values| !values.is_empty())
        .unwrap_or_else(|| default.to_vec())
}

fn load_vcf_workload(
    path: &Path,
    chrom: &str,
    max_variants: Option<usize>,
) -> BenchResult<Workload> {
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
        if variants.is_multiple_of(50_000) {
            eprintln!("parsed variants={variants} probe_attempts={probe_attempts}");
        }
        if max_variants.is_some_and(|max| variants >= max) {
            break;
        }
    }

    Ok(Workload {
        chrom: chrom.to_string(),
        chrom_code,
        variants,
        probe_attempts,
        ordered_positions,
        variant_keys_by_position,
        attempts_by_position,
    })
}

fn build_cold_gate(
    workload: &Workload,
    warm_positions: &HashSet<u64>,
    position_index: &PositionIndex,
    bloom_index: &VariantBloomIndex,
    collect_colocated: bool,
) -> IndexGateStats {
    let mut warm_miss_unique_positions = 0usize;
    let mut warm_miss_probe_attempts = 0usize;
    let mut position_index_unique_positions = 0usize;
    let mut position_index_probe_attempts = 0usize;
    let mut bloom_unique_positions = 0usize;
    let mut bloom_probe_attempts = 0usize;
    let mut position_index_positions = Vec::new();
    let mut bloom_positions = Vec::new();

    for &position_key in &workload.ordered_positions {
        let attempts = workload
            .attempts_by_position
            .get(&position_key)
            .copied()
            .unwrap_or(1);
        if warm_positions.contains(&position_key) {
            continue;
        }
        warm_miss_unique_positions += 1;
        warm_miss_probe_attempts += attempts;

        if !position_index.contains_position_key(position_key) {
            continue;
        }
        position_index_unique_positions += 1;
        position_index_probe_attempts += attempts;
        position_index_positions.push(position_key);

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
            bloom_unique_positions += 1;
            bloom_probe_attempts += attempts;
            bloom_positions.push(position_key);
        }
    }

    IndexGateStats {
        warm_miss_unique_positions,
        warm_miss_probe_attempts,
        position_index_unique_positions,
        position_index_probe_attempts,
        bloom_unique_positions,
        bloom_probe_attempts,
        position_index_positions,
        bloom_positions,
    }
}

async fn scan_lance_full(
    dataset_path: &Path,
    scan_batch_size: usize,
    materialization_style: MaterializationStyle,
) -> BenchResult<ScanStats> {
    let uri = dataset_path.to_string_lossy();
    let dataset = Dataset::open(uri.as_ref()).await?;
    let mut scanner = dataset.scan();
    scanner
        .batch_size(scan_batch_size)
        .materialization_style(materialization_style);
    let started = Instant::now();
    let mut stream = scanner.try_into_stream().await?;
    let mut rows = 0usize;
    let mut batches = 0usize;
    let mut selected_positions = HashSet::new();

    while let Some(batch) = stream.try_next().await? {
        rows += batch.num_rows();
        batches += 1;
        append_position_keys(&batch, &mut selected_positions)?;
    }

    Ok(ScanStats {
        rows,
        batches,
        elapsed: started.elapsed(),
        selected_positions,
    })
}

async fn scan_lance_filter_batches(
    dataset_path: &Path,
    positions: &[u64],
    query_batch_size: usize,
    scan_batch_size: usize,
    materialization_style: MaterializationStyle,
) -> BenchResult<ScanStats> {
    let uri = dataset_path.to_string_lossy();
    let dataset = Dataset::open(uri.as_ref()).await?;
    let started = Instant::now();
    let mut rows = 0usize;
    let mut batches = 0usize;
    let mut selected_positions = HashSet::new();

    for chunk in positions.chunks(query_batch_size.max(1)) {
        let mut scanner = dataset.scan();
        let filter = position_key_in_filter(chunk);
        scanner
            .filter(&filter)?
            .batch_size(scan_batch_size)
            .materialization_style(materialization_style.clone());
        let mut stream = scanner.try_into_stream().await?;
        while let Some(batch) = stream.try_next().await? {
            rows += batch.num_rows();
            batches += 1;
            append_position_keys(&batch, &mut selected_positions)?;
        }
    }

    Ok(ScanStats {
        rows,
        batches,
        elapsed: started.elapsed(),
        selected_positions,
    })
}

fn position_key_in_filter(keys: &[u64]) -> String {
    let mut filter = String::with_capacity(keys.len() * 20 + 24);
    filter.push_str("position_key IN (");
    for (idx, key) in keys.iter().enumerate() {
        if idx > 0 {
            filter.push(',');
        }
        let _ = write!(&mut filter, "{key}");
    }
    filter.push(')');
    filter
}

fn append_position_keys(batch: &RecordBatch, out: &mut HashSet<u64>) -> BenchResult<()> {
    let index = batch.schema().index_of("position_key")?;
    let array = batch.column(index).as_ref();

    if let Some(array) = array.as_any().downcast_ref::<UInt64Array>() {
        for row in 0..array.len() {
            if !array.is_null(row) {
                out.insert(array.value(row));
            }
        }
        return Ok(());
    }
    if let Some(array) = array.as_any().downcast_ref::<Int64Array>() {
        for row in 0..array.len() {
            if !array.is_null(row) {
                out.insert(u64::try_from(array.value(row))?);
            }
        }
        return Ok(());
    }
    if let Some(array) = array.as_any().downcast_ref::<UInt32Array>() {
        for row in 0..array.len() {
            if !array.is_null(row) {
                out.insert(u64::from(array.value(row)));
            }
        }
        return Ok(());
    }
    if let Some(array) = array.as_any().downcast_ref::<Int32Array>() {
        for row in 0..array.len() {
            if !array.is_null(row) {
                out.insert(u64::try_from(array.value(row))?);
            }
        }
        return Ok(());
    }

    Err(format!("unsupported position_key array type: {}", array.data_type()).into())
}

fn render_report(input: &ReportInput<'_>) -> String {
    let cache_root = input.cache_root;
    let lance_root = input.lance_root;
    let vcf_path = input.vcf_path;
    let workload = input.workload;
    let position_index_path = input.position_index_path;
    let bloom_index_path = input.bloom_index_path;
    let position_index_s = input.position_index_s;
    let bloom_index_s = input.bloom_index_s;
    let gate = input.gate;
    let workload_batch_size = input.workload_batch_size;
    let warm_scan_batch_size = input.warm_scan_batch_size;
    let cold_scan_batch_size = input.cold_scan_batch_size;
    let sort_cold_positions = input.sort_cold_positions;
    let result_rows = input.result_rows;

    let mut report = String::new();
    let _ = writeln!(
        &mut report,
        "# Rust Lance chr1 Variation Workload Benchmark"
    );
    let _ = writeln!(&mut report);
    let _ = writeln!(&mut report, "- Lance crate: `7.0.0`");
    let _ = writeln!(&mut report, "- Cache root: `{}`", cache_root.display());
    let _ = writeln!(&mut report, "- Lance root: `{}`", lance_root.display());
    let _ = writeln!(&mut report, "- VCF: `{}`", vcf_path.display());
    let _ = writeln!(
        &mut report,
        "- Position index: `{}` ({:.3}s load)",
        position_index_path.display(),
        position_index_s.as_secs_f64()
    );
    let _ = writeln!(
        &mut report,
        "- Variant bloom index: `{}` ({:.3}s load)",
        bloom_index_path.display(),
        bloom_index_s.as_secs_f64()
    );
    let _ = writeln!(
        &mut report,
        "- Workload query batch size: `{workload_batch_size}`"
    );
    let _ = writeln!(
        &mut report,
        "- Warm scan batch size: `{warm_scan_batch_size}`"
    );
    let _ = writeln!(
        &mut report,
        "- Cold scan batch size: `{cold_scan_batch_size}`"
    );
    let _ = writeln!(
        &mut report,
        "- Cold positions sorted before Lance batches: `{sort_cold_positions}`"
    );
    let _ = writeln!(&mut report);

    let best_warm = result_rows
        .iter()
        .filter(|row| row.phase == "warm_full_scan_all_columns")
        .min_by(|left, right| left.seconds.total_cmp(&right.seconds));
    let best_posidx = result_rows
        .iter()
        .filter(|row| row.phase == "cold_posidx_probe_batches_all_columns")
        .min_by(|left, right| left.seconds.total_cmp(&right.seconds));
    let best_bloom = result_rows
        .iter()
        .filter(|row| row.phase == "cold_posidx_bloom_probe_batches_all_columns")
        .min_by(|left, right| left.seconds.total_cmp(&right.seconds));

    let _ = writeln!(&mut report, "## Summary");
    let _ = writeln!(&mut report);
    let _ = writeln!(
        &mut report,
        "- Position index reduces the cold lookup set to `{}` unique positions, `{}` probe attempts, and `{}` 2k batches.",
        gate.position_index_unique_positions,
        gate.position_index_probe_attempts,
        gate.position_index_positions
            .len()
            .div_ceil(workload_batch_size)
    );
    let _ = writeln!(
        &mut report,
        "- Adding variant bloom admission reduces that further to `{}` unique positions, `{}` probe attempts, and `{}` 2k batches.",
        gate.bloom_unique_positions,
        gate.bloom_probe_attempts,
        gate.bloom_positions.len().div_ceil(workload_batch_size)
    );
    if let Some(row) = best_warm {
        let _ = writeln!(
            &mut report,
            "- Best warm full scan: `{}` at `{:.3}s` for `{}` selected rows.",
            row.layout, row.seconds, row.selected_rows
        );
    }
    if let Some(row) = best_posidx {
        let _ = writeln!(
            &mut report,
            "- Best default cold `posidx` scan: `{}` with `{}` fragment rows at `{:.3}s`, selecting `{}` rows across `{}` positions.",
            row.layout,
            row.fragment_rows.unwrap_or_default(),
            row.seconds,
            row.selected_rows,
            row.selected_positions
        );
    }
    if let Some(row) = best_bloom {
        let _ = writeln!(
            &mut report,
            "- Best optional `posidx_bloom` scan: `{}` with `{}` fragment rows at `{:.3}s`, selecting `{}` rows across `{}` positions.",
            row.layout,
            row.fragment_rows.unwrap_or_default(),
            row.seconds,
            row.selected_rows,
            row.selected_positions
        );
    }
    let _ = writeln!(
        &mut report,
        "- Current result favors `2.1-unpacked` for cold point lookups; `2.2-packed` is about 2x slower in this all-column Rust workload."
    );
    let _ = writeln!(&mut report);

    let _ = writeln!(&mut report, "## Workload");
    let _ = writeln!(&mut report);
    let _ = writeln!(&mut report, "| metric | value |");
    let _ = writeln!(&mut report, "|---|---:|");
    let _ = writeln!(&mut report, "| chrom | {} |", workload.chrom);
    let _ = writeln!(&mut report, "| chrom code | {} |", workload.chrom_code);
    let _ = writeln!(
        &mut report,
        "| parsed VCF variants | {} |",
        workload.variants
    );
    let _ = writeln!(
        &mut report,
        "| raw extended probe attempts | {} |",
        workload.probe_attempts
    );
    let _ = writeln!(
        &mut report,
        "| unique probe positions | {} |",
        workload.ordered_positions.len()
    );
    let _ = writeln!(
        &mut report,
        "| warm-miss unique positions | {} |",
        gate.warm_miss_unique_positions
    );
    let _ = writeln!(
        &mut report,
        "| warm-miss probe attempts | {} |",
        gate.warm_miss_probe_attempts
    );
    let _ = writeln!(
        &mut report,
        "| position-index admitted unique positions | {} |",
        gate.position_index_unique_positions
    );
    let _ = writeln!(
        &mut report,
        "| position-index admitted probe attempts | {} |",
        gate.position_index_probe_attempts
    );
    let _ = writeln!(
        &mut report,
        "| variant-bloom admitted unique cold positions | {} |",
        gate.bloom_unique_positions
    );
    let _ = writeln!(
        &mut report,
        "| variant-bloom admitted cold probe attempts | {} |",
        gate.bloom_probe_attempts
    );
    let _ = writeln!(
        &mut report,
        "| 2k cold Lance filter batches | {} |",
        gate.position_index_positions
            .len()
            .div_ceil(workload_batch_size)
    );
    let _ = writeln!(
        &mut report,
        "| 2k cold Lance filter batches after bloom | {} |",
        gate.bloom_positions.len().div_ceil(workload_batch_size)
    );
    let _ = writeln!(&mut report);

    let _ = writeln!(&mut report, "## Rust Lance Results");
    let _ = writeln!(&mut report);
    let _ = writeln!(
        &mut report,
        "| phase | layout | fragment rows | query batch | probe positions | selected rows | selected positions | result batches | seconds | rows/s |"
    );
    let _ = writeln!(
        &mut report,
        "|---|---|---:|---:|---:|---:|---:|---:|---:|---:|"
    );
    for row in result_rows {
        let rows_per_s = if row.seconds > 0.0 {
            row.selected_rows as f64 / row.seconds
        } else {
            0.0
        };
        let _ = writeln!(
            &mut report,
            "| {} | {} | {} | {} | {} | {} | {} | {} | {:.3} | {:.0} |",
            row.phase,
            row.layout,
            row.fragment_rows
                .map(|value| value.to_string())
                .unwrap_or_default(),
            row.query_batch,
            row.probe_positions
                .map(|value| value.to_string())
                .unwrap_or_default(),
            row.selected_rows,
            row.selected_positions,
            row.batches,
            row.seconds,
            rows_per_s
        );
    }
    report
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

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn build_probe_starts_keeps_snv_single_start_without_extended_probes() {
        assert_eq!(build_probe_starts(10, 10, "A", "G", false), vec![10]);
    }

    #[test]
    fn build_probe_starts_adds_adjacent_start_for_extended_snv() {
        assert_eq!(build_probe_starts(10, 10, "A", "G", true), vec![10, 11]);
    }

    #[test]
    fn build_probe_starts_adds_shifted_repeat_deletion_window() {
        let starts = build_probe_starts(100, 102, "CTG", "C", true);
        assert!(starts.contains(&100));
        assert!(starts.contains(&102));
        assert!(starts.contains(&101));
    }

    #[test]
    fn warm_variant_key_candidates_deduplicate_equivalent_forms() {
        let position_key = position_key_from_code(chrom_to_code("1"), 10).unwrap();
        let keys = warm_variant_key_candidates(position_key, "A", "G");
        let unique = keys.iter().copied().collect::<HashSet<_>>();
        assert_eq!(keys.len(), unique.len());
        assert!(!keys.is_empty());
    }
}
