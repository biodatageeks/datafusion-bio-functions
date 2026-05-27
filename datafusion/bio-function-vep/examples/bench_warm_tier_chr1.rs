//! Benchmark the warm variation tier against direct Fjall point lookups on chr1.
//!
//! Usage:
//!   RUSTFLAGS="-C target-cpu=native" cargo run --release \
//!     -p datafusion-bio-function-vep --features kv-cache \
//!     --example bench_warm_tier_chr1 -- \
//!     /path/to/variation.fjall /tmp/vepyr-warm-prototype/variation/chr1_warm.parquet \
//!     /path/to/input.vcf.gz
//!
//! Optional env:
//!   MAX_VARIANTS=200000     cap parsed chr1 variants
//!   CACHE_MB=1024           Fjall block cache size
//!   WARM_BATCH_SIZE=65536   Arrow batch size when reading warm row groups

use std::env;
use std::fs::File;
use std::io::{BufRead, BufReader};
use std::path::{Path, PathBuf};
use std::time::{Duration, Instant};

use datafusion::arrow::array::{Array, StringArray};
use datafusion::common::{DataFusionError, Result};
use datafusion_bio_function_vep::allele::{allele_matches, vcf_to_vep_input_allele};
use datafusion_bio_function_vep::kv_cache::key_encoding::{chrom_to_code, encode_position_key_buf};
use datafusion_bio_function_vep::kv_cache::kv_store::VepKvStore;
use datafusion_bio_function_vep::kv_cache::position_entry::PositionEntryReader;
use datafusion_bio_function_vep::warm_cache::chunk::{WarmChunkContext, WarmProbeResult};
use datafusion_bio_function_vep::warm_cache::key::{position_key, variant_key};
use datafusion_bio_function_vep::warm_cache::reader::load_warm_chunk_row_group;
use flate2::read::MultiGzDecoder;
use parquet::arrow::arrow_reader::ParquetRecordBatchReaderBuilder;

type DynError = Box<dyn std::error::Error + Send + Sync>;
type BenchResult<T> = std::result::Result<T, DynError>;

#[derive(Debug)]
struct ProbeSet {
    chrom: String,
    chrom_code: u16,
    variants: usize,
    probes: Vec<Probe>,
}

#[derive(Debug)]
struct Probe {
    start: i64,
    ref_allele: String,
    alt_allele: String,
    candidate_keys: Vec<u64>,
}

#[derive(Debug, Default)]
struct LookupStats {
    probes: usize,
    found_entries: usize,
    missing_entries: usize,
    matched_alleles: usize,
    bytes: u64,
    get_s: Duration,
    decompress_s: Duration,
    reader_s: Duration,
    match_s: Duration,
    total_s: Duration,
}

#[derive(Debug, Default)]
struct WarmStats {
    probes: usize,
    warm_hits: usize,
    warm_unverified_hits: usize,
    warm_definitive_misses: usize,
    warm_not_covered: usize,
    fjall_fallbacks: usize,
    fjall_saved: usize,
    chunks_loaded: usize,
    chunk_rows: usize,
    chunk_load_s: Duration,
    warm_probe_s: Duration,
    fallback: LookupStats,
    total_s: Duration,
}

struct WarmWindow {
    path: PathBuf,
    batch_size: usize,
    total_row_groups: usize,
    next_row_group: usize,
    previous: Option<WarmChunkContext>,
    current: Option<WarmChunkContext>,
    chunks_loaded: usize,
    chunk_rows: usize,
    chunk_load_s: Duration,
}

enum WarmProbeOutcome {
    Hit,
    DefinitiveMiss,
    NotCovered,
    UnverifiedHit,
}

fn main() -> BenchResult<()> {
    let args = env::args().collect::<Vec<_>>();
    if args.len() != 4 {
        eprintln!(
            "usage: {} <variation.fjall> <chr1_warm.parquet> <input.vcf.gz>",
            args[0]
        );
        std::process::exit(2);
    }

    let fjall_path = PathBuf::from(&args[1]);
    let warm_path = PathBuf::from(&args[2]);
    let vcf_path = PathBuf::from(&args[3]);
    let max_variants = env::var("MAX_VARIANTS")
        .ok()
        .and_then(|v| v.parse::<usize>().ok());
    let cache_mb = env::var("CACHE_MB")
        .ok()
        .and_then(|v| v.parse::<u64>().ok())
        .unwrap_or(1024);
    let warm_batch_size = env::var("WARM_BATCH_SIZE")
        .ok()
        .and_then(|v| v.parse::<usize>().ok())
        .unwrap_or(65_536);

    println!(
        "fjall={}\nwarm={}\nvcf={}\ncache_mb={cache_mb}\nwarm_batch_size={warm_batch_size}\nmax_variants={}",
        fjall_path.display(),
        warm_path.display(),
        vcf_path.display(),
        max_variants
            .map(|v| v.to_string())
            .unwrap_or_else(|| "all chr1".to_string())
    );

    let probes = parse_chr1_probe_stream(&vcf_path, max_variants)?;
    println!(
        "parsed chr1 variants={} probes={} avg_probes_per_variant={:.2}",
        probes.variants,
        probes.probes.len(),
        probes.probes.len() as f64 / probes.variants.max(1) as f64
    );

    let store = VepKvStore::open_with_cache_size(&fjall_path, cache_mb * 1024 * 1024)?;

    println!("\n=== fjall baseline ===");
    let baseline = run_fjall_lookup(&store, probes.chrom_code, &probes.probes)?;
    print_lookup_stats("fjall", &baseline);

    println!("\n=== warm tier + fjall fallback ===");
    let warm = run_warm_tier_lookup(&store, &warm_path, warm_batch_size, &probes)?;
    print_warm_stats(&warm);

    Ok(())
}

fn parse_chr1_probe_stream(path: &Path, max_variants: Option<usize>) -> BenchResult<ProbeSet> {
    let file = File::open(path)?;
    let decoder = MultiGzDecoder::new(file);
    let reader = BufReader::new(decoder);
    let chrom = "1".to_string();
    let chrom_code = chrom_to_code(&chrom);
    let mut variants = 0usize;
    let mut probes = Vec::new();
    let mut seen_chr1 = false;

    for line in reader.lines() {
        let line = line?;
        if line.starts_with('#') {
            continue;
        }

        let mut fields = line.split('\t');
        let Some(row_chrom) = fields.next() else {
            continue;
        };
        let is_chr1 = row_chrom == "1" || row_chrom == "chr1";
        if !is_chr1 {
            if seen_chr1 {
                break;
            }
            continue;
        }
        seen_chr1 = true;

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
        for start in build_probe_starts(pos, end, ref_allele, alt_allele, true) {
            probes.push(Probe {
                start,
                ref_allele: ref_allele.to_string(),
                alt_allele: alt_allele.to_string(),
                candidate_keys: candidate_variant_keys(&chrom, pos, start, ref_allele, alt_allele),
            });
        }

        variants += 1;
        if variants.is_multiple_of(50_000) {
            eprintln!("parsed chr1 variants={variants} probes={}", probes.len());
        }
        if max_variants.is_some_and(|max| variants >= max) {
            break;
        }
    }

    Ok(ProbeSet {
        chrom,
        chrom_code,
        variants,
        probes,
    })
}

fn candidate_variant_keys(
    chrom: &str,
    original_pos: i64,
    probe_start: i64,
    ref_allele: &str,
    alt_allele: &str,
) -> Vec<u64> {
    let mut keys = Vec::new();
    for alt in alt_allele.split(['|', ',']).filter(|a| !a.is_empty()) {
        push_variant_key(&mut keys, chrom, probe_start, ref_allele, alt);

        let (input_ref, input_alt, input_start) =
            vcf_to_vep_input_allele(original_pos, ref_allele, alt);
        if input_start == probe_start {
            push_variant_key(&mut keys, chrom, input_start, &input_ref, &input_alt);
        }
    }
    keys.sort_unstable();
    keys.dedup();
    keys
}

fn push_variant_key(
    keys: &mut Vec<u64>,
    chrom: &str,
    start: i64,
    reference: &str,
    alternate: &str,
) {
    if let Ok(key) = variant_key(chrom, start, reference, alternate) {
        keys.push(key);
    }
}

fn run_fjall_lookup(store: &VepKvStore, chrom_code: u16, probes: &[Probe]) -> Result<LookupStats> {
    let total_started = Instant::now();
    let mut stats = LookupStats {
        probes: probes.len(),
        ..Default::default()
    };
    let mut key_buf = Vec::with_capacity(10);
    let mut decompress_buf = Vec::new();
    let mut decompressor = store.create_decompressor()?;

    for probe in probes {
        run_one_fjall_probe(
            store,
            chrom_code,
            probe,
            &mut key_buf,
            decompressor.as_mut(),
            &mut decompress_buf,
            &mut stats,
        )?;
    }

    stats.total_s = total_started.elapsed();
    Ok(stats)
}

fn run_warm_tier_lookup(
    store: &VepKvStore,
    warm_path: &Path,
    warm_batch_size: usize,
    probes: &ProbeSet,
) -> Result<WarmStats> {
    let total_started = Instant::now();
    let mut warm = WarmWindow::open(warm_path, warm_batch_size)?;
    let mut stats = WarmStats {
        probes: probes.probes.len(),
        ..Default::default()
    };
    let mut key_buf = Vec::with_capacity(10);
    let mut decompress_buf = Vec::new();
    let mut decompressor = store.create_decompressor()?;

    for probe in &probes.probes {
        let probe_started = Instant::now();
        let outcome = warm.probe(&probes.chrom, probe)?;
        stats.warm_probe_s += probe_started.elapsed();

        match outcome {
            WarmProbeOutcome::Hit => {
                stats.warm_hits += 1;
                stats.fjall_saved += 1;
            }
            WarmProbeOutcome::UnverifiedHit => {
                stats.warm_unverified_hits += 1;
                stats.warm_not_covered += 1;
                stats.fjall_fallbacks += 1;
                run_one_fjall_probe(
                    store,
                    probes.chrom_code,
                    probe,
                    &mut key_buf,
                    decompressor.as_mut(),
                    &mut decompress_buf,
                    &mut stats.fallback,
                )?;
            }
            WarmProbeOutcome::DefinitiveMiss => {
                stats.warm_definitive_misses += 1;
                stats.fjall_saved += 1;
            }
            WarmProbeOutcome::NotCovered => {
                stats.warm_not_covered += 1;
                stats.fjall_fallbacks += 1;
                run_one_fjall_probe(
                    store,
                    probes.chrom_code,
                    probe,
                    &mut key_buf,
                    decompressor.as_mut(),
                    &mut decompress_buf,
                    &mut stats.fallback,
                )?;
            }
        }
    }

    stats.chunks_loaded = warm.chunks_loaded;
    stats.chunk_rows = warm.chunk_rows;
    stats.chunk_load_s = warm.chunk_load_s;
    stats.fallback.probes = stats.fjall_fallbacks;
    stats.fallback.total_s = stats
        .fallback
        .get_s
        .saturating_add(stats.fallback.decompress_s)
        .saturating_add(stats.fallback.reader_s)
        .saturating_add(stats.fallback.match_s);
    stats.total_s = total_started.elapsed();
    Ok(stats)
}

fn run_one_fjall_probe(
    store: &VepKvStore,
    chrom_code: u16,
    probe: &Probe,
    key_buf: &mut Vec<u8>,
    decompressor: Option<&mut zstd::bulk::Decompressor<'_>>,
    decompress_buf: &mut Vec<u8>,
    stats: &mut LookupStats,
) -> Result<()> {
    encode_position_key_buf(chrom_code, probe.start, key_buf);

    let get_started = Instant::now();
    let raw = store
        .data_partition()
        .get(key_buf.as_slice())
        .map_err(|e| DataFusionError::External(Box::new(e)))?;
    stats.get_s += get_started.elapsed();

    let Some(raw) = raw else {
        stats.missing_entries += 1;
        return Ok(());
    };

    stats.found_entries += 1;
    stats.bytes += raw.len() as u64;

    let decompress_started = Instant::now();
    store.decode_position_entry_value(&raw, decompressor, decompress_buf)?;
    stats.decompress_s += decompress_started.elapsed();

    let reader_started = Instant::now();
    let reader = PositionEntryReader::new(decompress_buf)?;
    stats.reader_s += reader_started.elapsed();

    let match_started = Instant::now();
    for allele_idx in 0..reader.num_alleles() {
        if allele_matches(
            &probe.ref_allele,
            &probe.alt_allele,
            reader.allele_string(allele_idx),
        ) {
            stats.matched_alleles += 1;
        }
    }
    stats.match_s += match_started.elapsed();

    Ok(())
}

impl WarmWindow {
    fn open(path: impl AsRef<Path>, batch_size: usize) -> Result<Self> {
        let path = path.as_ref().to_path_buf();
        let file = File::open(&path).map_err(|e| {
            DataFusionError::Execution(format!("failed to open warm parquet file: {e}"))
        })?;
        let builder = ParquetRecordBatchReaderBuilder::try_new(file)?;
        let total_row_groups = builder.metadata().num_row_groups();
        let mut window = Self {
            path,
            batch_size,
            total_row_groups,
            next_row_group: 0,
            previous: None,
            current: None,
            chunks_loaded: 0,
            chunk_rows: 0,
            chunk_load_s: Duration::ZERO,
        };
        window.load_next()?;
        Ok(window)
    }

    fn load_next(&mut self) -> Result<bool> {
        if self.next_row_group >= self.total_row_groups {
            return Ok(false);
        }
        let row_group_id = self.next_row_group;
        self.next_row_group += 1;

        let started = Instant::now();
        let chunk = load_warm_chunk_row_group(&self.path, row_group_id, self.batch_size)?;
        self.chunk_load_s += started.elapsed();
        self.chunk_rows += chunk.batch.num_rows();
        self.chunks_loaded += 1;
        self.previous = self.current.take();
        self.current = Some(chunk);
        Ok(true)
    }

    fn probe(&mut self, chrom: &str, probe: &Probe) -> Result<WarmProbeOutcome> {
        let position = position_key(chrom, probe.start).map_err(|e| {
            DataFusionError::Execution(format!("failed to build warm position key: {e}"))
        })?;

        while self
            .current
            .as_ref()
            .is_some_and(|chunk| position > chunk.max_position_key)
            && self.load_next()?
        {}

        let mut saw_unverified_hit = false;
        let mut covered = false;

        if let Some(outcome) = self.probe_chunk(self.current.as_ref(), position, probe)? {
            match outcome {
                WarmProbeOutcome::Hit => return Ok(WarmProbeOutcome::Hit),
                WarmProbeOutcome::UnverifiedHit => saw_unverified_hit = true,
                WarmProbeOutcome::DefinitiveMiss => covered = true,
                WarmProbeOutcome::NotCovered => {}
            }
        }

        if let Some(outcome) = self.probe_chunk(self.previous.as_ref(), position, probe)? {
            match outcome {
                WarmProbeOutcome::Hit => return Ok(WarmProbeOutcome::Hit),
                WarmProbeOutcome::UnverifiedHit => saw_unverified_hit = true,
                WarmProbeOutcome::DefinitiveMiss => covered = true,
                WarmProbeOutcome::NotCovered => {}
            }
        }

        if saw_unverified_hit {
            Ok(WarmProbeOutcome::UnverifiedHit)
        } else if covered {
            Ok(WarmProbeOutcome::DefinitiveMiss)
        } else {
            Ok(WarmProbeOutcome::NotCovered)
        }
    }

    fn probe_chunk(
        &self,
        chunk: Option<&WarmChunkContext>,
        position: u64,
        probe: &Probe,
    ) -> Result<Option<WarmProbeOutcome>> {
        let Some(chunk) = chunk else {
            return Ok(None);
        };
        let mut covered = false;

        for &candidate_key in &probe.candidate_keys {
            match chunk.probe(position, candidate_key) {
                WarmProbeResult::Hit(rows) => {
                    if warm_rows_match(chunk, rows.as_slice(), probe)? {
                        return Ok(Some(WarmProbeOutcome::Hit));
                    }
                    covered = true;
                }
                WarmProbeResult::DefinitiveMiss => covered = true,
                WarmProbeResult::NotCovered => {}
            }
        }

        if !covered && probe.candidate_keys.is_empty() && chunk.contains_position(position) {
            covered = true;
        }

        if covered {
            Ok(Some(WarmProbeOutcome::DefinitiveMiss))
        } else {
            Ok(Some(WarmProbeOutcome::NotCovered))
        }
    }
}

fn warm_rows_match(chunk: &WarmChunkContext, rows: &[u32], probe: &Probe) -> Result<bool> {
    let allele_idx = chunk
        .batch
        .schema()
        .index_of("allele_string")
        .map_err(|_| {
            DataFusionError::Execution("warm chunk is missing allele_string column".to_string())
        })?;
    let alleles = chunk
        .batch
        .column(allele_idx)
        .as_any()
        .downcast_ref::<StringArray>()
        .ok_or_else(|| {
            DataFusionError::Execution("warm allele_string column must be Utf8".to_string())
        })?;

    for &row in rows {
        let row = row as usize;
        if row < alleles.len()
            && !alleles.is_null(row)
            && allele_matches(&probe.ref_allele, &probe.alt_allele, alleles.value(row))
        {
            return Ok(true);
        }
    }

    Ok(false)
}

fn print_lookup_stats(label: &str, stats: &LookupStats) {
    let probes_s = stats.probes as f64 / stats.total_s.as_secs_f64().max(f64::EPSILON);
    println!(
        "{label}: total={:.3}s probes={} probes/s={probes_s:.0} found={} missing={} matched_alleles={} bytes={:.2}GiB",
        stats.total_s.as_secs_f64(),
        stats.probes,
        stats.found_entries,
        stats.missing_entries,
        stats.matched_alleles,
        stats.bytes as f64 / 1_073_741_824.0,
    );
    println!(
        "{label}-detail: get={:.3}s decompress={:.3}s reader={:.3}s match={:.3}s",
        stats.get_s.as_secs_f64(),
        stats.decompress_s.as_secs_f64(),
        stats.reader_s.as_secs_f64(),
        stats.match_s.as_secs_f64(),
    );
}

fn print_warm_stats(stats: &WarmStats) {
    let probes_s = stats.probes as f64 / stats.total_s.as_secs_f64().max(f64::EPSILON);
    let saved_pct = stats.fjall_saved as f64 * 100.0 / stats.probes.max(1) as f64;
    println!(
        "warm: total={:.3}s probes={} probes/s={probes_s:.0} saved={} ({saved_pct:.1}%) fallbacks={} hits={} unverified_hits={} definitive_misses={} not_covered={}",
        stats.total_s.as_secs_f64(),
        stats.probes,
        stats.fjall_saved,
        stats.fjall_fallbacks,
        stats.warm_hits,
        stats.warm_unverified_hits,
        stats.warm_definitive_misses,
        stats.warm_not_covered,
    );
    println!(
        "warm-detail: chunks_loaded={} chunk_rows={} chunk_load_index={:.3}s warm_probe={:.3}s",
        stats.chunks_loaded,
        stats.chunk_rows,
        stats.chunk_load_s.as_secs_f64(),
        stats.warm_probe_s.as_secs_f64(),
    );
    print_lookup_stats("warm-fallback", &stats.fallback);
}

#[inline]
fn push_unique_probe_start(probe_starts: &mut Vec<i64>, start: i64) {
    if !probe_starts.contains(&start) {
        probe_starts.push(start);
    }
}

fn build_probe_starts(
    norm_start_i64: i64,
    norm_end_i64: i64,
    vcf_ref: &str,
    vcf_alt: &str,
    extended_probes: bool,
) -> Vec<i64> {
    let mut probe_starts: Vec<i64> = Vec::with_capacity(6);
    push_unique_probe_start(&mut probe_starts, norm_start_i64);

    if !extended_probes {
        return probe_starts;
    }

    if norm_start_i64 == norm_end_i64 {
        push_unique_probe_start(&mut probe_starts, norm_start_i64.saturating_add(1));
    } else {
        push_unique_probe_start(&mut probe_starts, norm_end_i64);
    }

    for alt in vcf_alt.split(['|', ',']).filter(|a| !a.is_empty()) {
        let (_, _, input_start) = vcf_to_vep_input_allele(norm_start_i64, vcf_ref, alt);
        push_unique_probe_start(&mut probe_starts, input_start);

        let shift_usize = common_prefix_len(vcf_ref, alt);
        if shift_usize == 0 {
            continue;
        }
        let shift = shift_usize as i64;
        if let Some(shifted_start) = norm_start_i64.checked_add(shift) {
            push_unique_probe_start(&mut probe_starts, shifted_start);
        }
    }

    for alt in vcf_alt.split(['|', ',']).filter(|a| !a.is_empty()) {
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
                push_unique_probe_start(&mut probe_starts, candidate_start);
            }
        }
    }

    probe_starts
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
