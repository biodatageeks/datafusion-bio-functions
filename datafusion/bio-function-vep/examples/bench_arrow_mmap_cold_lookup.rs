//! Compare Lance cold lookups against a benchmark-only Arrow IPC mmap store
//! indexed by a shared packed ZSTD position-range sidecar.
//!
//! Usage:
//!   RUSTFLAGS="-C target-cpu=native" cargo run --release \
//!     -p datafusion-bio-function-vep --features lance-cache \
//!     --example bench_arrow_mmap_cold_lookup -- \
//!     /path/to/115_GRCh38_merged \
//!     /path/to/input_chr1.vcf.gz \
//!     /path/to/variation.lance/chr1.lance \
//!     /path/to/arrow_mmap_chr1_bench \
//!     research/reports/arrow_mmap_cold_lookup_20260611/chr1_arrow_mmap_cold_lookup.md

use std::collections::{HashMap, HashSet};
use std::env;
use std::fmt::Write as _;
use std::fs::{self, File};
use std::io::{BufRead, BufReader, Cursor, Read, Write};
use std::path::{Path, PathBuf};
use std::sync::{Arc, Mutex};
use std::time::{Duration, Instant};

use arrow_ipc::reader::FileReader as IpcFileReader;
use arrow_ipc::writer::{FileWriter as IpcFileWriter, IpcWriteOptions};
use arrow_ipc::{CompressionType, MetadataVersion};
use datafusion::arrow::array::{
    Array, BinaryArray, BooleanArray, DictionaryArray, FixedSizeBinaryArray, Float32Array,
    Float64Array, GenericListArray, Int8Array, Int16Array, Int32Array, Int64Array,
    LargeBinaryArray, LargeListArray, LargeStringArray, ListArray, NullArray, RecordBatch,
    StringArray, StructArray, UInt8Array, UInt16Array, UInt32Array, UInt64Array, make_array,
};
use datafusion::arrow::compute::concat_batches;
use datafusion::arrow::datatypes::{ArrowNativeType, DataType, Field, Fields, Schema, SchemaRef};
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
use memmap2::Mmap;
use parquet::arrow::ProjectionMask;
use parquet::arrow::arrow_reader::ParquetRecordBatchReaderBuilder;
use serde::{Deserialize, Serialize};

type DynError = Box<dyn std::error::Error + Send + Sync>;
type BenchResult<T> = std::result::Result<T, DynError>;

const FORMAT_NAME: &str = "vep-arrow-mmap-cold";
const FORMAT_VERSION: u32 = 1;
const INDEX_KIND_PACKED_ZSTD: &str = "packed_zstd";
const PACKED_INDEX_MAGIC: &[u8; 8] = b"VPRNGZ1\0";
const DEFAULT_INDEX_BLOCK_ENTRIES: usize = 4096;
const DEFAULT_INDEX_ZSTD_LEVEL: i32 = 3;

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
struct LancePlanResult {
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

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
struct PositionRange {
    batch_id: u32,
    start_row: u32,
    len: u32,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
enum ArrowStoreCompression {
    Uncompressed,
    Lz4,
    Zstd,
}

impl ArrowStoreCompression {
    fn all() -> [Self; 3] {
        [Self::Uncompressed, Self::Lz4, Self::Zstd]
    }

    fn as_str(self) -> &'static str {
        match self {
            Self::Uncompressed => "uncompressed",
            Self::Lz4 => "lz4",
            Self::Zstd => "zstd",
        }
    }

    fn dirname(self, chrom: &str) -> String {
        format!("{chrom}.arrow_mmap_{}", self.as_str())
    }

    fn ipc_compression(self) -> Option<CompressionType> {
        match self {
            Self::Uncompressed => None,
            Self::Lz4 => Some(CompressionType::LZ4_FRAME),
            Self::Zstd => Some(CompressionType::ZSTD),
        }
    }
}

#[derive(Debug, Serialize, Deserialize, Clone)]
struct ArrowMmapManifest {
    format: String,
    version: u32,
    chrom: String,
    tier: String,
    compression: String,
    projected_columns: Vec<String>,
    target_batch_rows: usize,
    batch_count: usize,
    total_rows: usize,
    unique_positions: usize,
    max_position_group_len: usize,
    data_file_size: u64,
    index_kind: String,
    index_path: String,
    index_block_entries: usize,
    index_size: u64,
    total_size: u64,
}

#[derive(Debug)]
struct ArrowBuildStats {
    compression: ArrowStoreCompression,
    elapsed: Duration,
    manifest: ArrowMmapManifest,
}

#[derive(Default, Debug, Clone)]
struct TouchAccumulator {
    checksum: u64,
    estimated_bytes_touched: u64,
    rows: usize,
}

impl TouchAccumulator {
    fn mix_u64(&mut self, value: u64) {
        self.checksum = self
            .checksum
            .wrapping_mul(0x9E37_79B1_85EB_CA87)
            .wrapping_add(value.rotate_left(17));
    }

    fn touch_bytes(&mut self, bytes: &[u8]) {
        self.estimated_bytes_touched += bytes.len() as u64;
        self.mix_u64(bytes.len() as u64);
        for chunk in bytes.chunks(8) {
            let mut padded = [0_u8; 8];
            padded[..chunk.len()].copy_from_slice(chunk);
            self.mix_u64(u64::from_le_bytes(padded));
        }
    }
}

#[derive(Debug)]
struct ArrowLookupResult {
    compression: ArrowStoreCompression,
    total_size: u64,
    data_size: u64,
    index_size: u64,
    open_elapsed: Duration,
    lookup_elapsed: Duration,
    positions_requested: usize,
    positions_found: usize,
    selected_rows: usize,
    projected_columns: usize,
    estimated_bytes_touched: u64,
    batch_files_opened: usize,
    index_gets: usize,
    missing_positions: usize,
    checksum: u64,
}

struct CachedBatch {
    _mmap: Mmap,
    batch: RecordBatch,
}

struct ArrowMmapStore {
    root: PathBuf,
    manifest: ArrowMmapManifest,
    index: PackedRangeIndexReader,
    batch_cache: HashMap<u32, CachedBatch>,
    batch_files_opened: usize,
    index_gets: usize,
    missing_positions: usize,
}

#[tokio::main(flavor = "multi_thread")]
async fn main() -> BenchResult<()> {
    let args = env::args().collect::<Vec<_>>();
    if args.len() < 5 || args.len() > 6 {
        eprintln!(
            "usage: {} <cache_root> <input.vcf.gz> <merged_lance_chr_dataset> <arrow_store_root> [report.md]",
            args[0]
        );
        std::process::exit(2);
    }

    let cache_root = PathBuf::from(&args[1]);
    let vcf_path = PathBuf::from(&args[2]);
    let dataset_path = PathBuf::from(&args[3]);
    let arrow_store_root = PathBuf::from(&args[4]);
    let report_path = args.get(5).map(PathBuf::from);
    let chrom = env::var("CHROM").unwrap_or_else(|_| "chr1".to_string());
    let batch_rows = env_usize("ARROW_MMAP_BATCH_ROWS", 16_384);
    let cold_scan_batch_size = env_usize("COLD_SCAN_BATCH_SIZE", 2_000);
    let small_batch_size = env_usize("SMALL_QUERY_BATCH_SIZE", 238);
    let large_batch_size = env_usize("LARGE_QUERY_BATCH_SIZE", 2_000);
    let index_block_entries = env_usize(
        "ARROW_MMAP_INDEX_BLOCK_ENTRIES",
        DEFAULT_INDEX_BLOCK_ENTRIES,
    );
    let sort_positions = env_bool("SORT_COLD_POSITIONS", true);
    let rebuild = env_bool("REBUILD_ARROW_MMAP", false);

    println!(
        "cache_root={}\nvcf={}\nlance_dataset={}\narrow_store_root={}\nchrom={chrom}\narrow_mmap_batch_rows={batch_rows}\ncold_scan_batch_size={cold_scan_batch_size}\nsmall_query_batch_size={small_batch_size}\nlarge_query_batch_size={large_batch_size}\narrow_mmap_index_kind={INDEX_KIND_PACKED_ZSTD}\narrow_mmap_index_block_entries={index_block_entries}\nsort_positions={sort_positions}\nrebuild_arrow_mmap={rebuild}",
        cache_root.display(),
        vcf_path.display(),
        dataset_path.display(),
        arrow_store_root.display(),
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

    let (projection, omitted_projection) =
        projected_cold_everything_columns(&dataset, &cache_root.join("variation"), &chrom)?;
    println!(
        "projection_cols={} omitted_cold_projection_cols={} projection={} omitted={}",
        projection.len(),
        omitted_projection.len(),
        projection.join(","),
        omitted_projection.join(",")
    );

    let mut build_stats = Vec::new();
    fs::create_dir_all(&arrow_store_root)?;
    let shared_index_path = arrow_store_root.join(format!("{chrom}.range_index.pzstd"));
    if rebuild && shared_index_path.exists() {
        fs::remove_file(&shared_index_path)?;
    }
    for compression in ArrowStoreCompression::all() {
        let store_path = arrow_store_root.join(compression.dirname(&chrom));
        if rebuild || !store_path.join("manifest.json").exists() {
            let build_index = !shared_index_path.exists();
            let stats = build_arrow_store(
                &cache_root,
                &chrom,
                &projection,
                &store_path,
                compression,
                batch_rows,
                &shared_index_path,
                build_index,
                index_block_entries,
            )?;
            println!(
                "built_arrow_store compression={} rows={} positions={} batches={} data_size={} index_size={} total_size={} build_s={:.3}",
                compression.as_str(),
                stats.manifest.total_rows,
                stats.manifest.unique_positions,
                stats.manifest.batch_count,
                stats.manifest.data_file_size,
                stats.manifest.index_size,
                stats.manifest.total_size,
                stats.elapsed.as_secs_f64()
            );
            build_stats.push(stats);
        } else {
            println!(
                "reuse_arrow_store compression={} path={}",
                compression.as_str(),
                store_path.display()
            );
            let manifest = read_manifest(&store_path)?;
            build_stats.push(ArrowBuildStats {
                compression,
                elapsed: Duration::ZERO,
                manifest,
            });
        }
    }

    let lance_results = vec![
        scan_lance_plan(
            &dataset,
            "small_chunks",
            &gate.bloom_positions,
            small_batch_size,
            cold_scan_batch_size,
            &projection,
            MaterializationStyle::AllLate,
        )
        .await?,
        scan_lance_plan(
            &dataset,
            "large_chunks",
            &gate.bloom_positions,
            large_batch_size,
            cold_scan_batch_size,
            &projection,
            MaterializationStyle::AllLate,
        )
        .await?,
    ];

    let mut arrow_results = Vec::new();
    for compression in ArrowStoreCompression::all() {
        let store_path = arrow_store_root.join(compression.dirname(&chrom));
        let open_started = Instant::now();
        let mut store = ArrowMmapStore::open(&store_path)?;
        let open_elapsed = open_started.elapsed();
        let result = store.run_lookup(
            compression,
            &gate.bloom_positions,
            projection.len(),
            open_elapsed,
        )?;
        println!(
            "arrow_lookup compression={} rows={} positions_found={} lookup_s={:.3} checksum={} estimated_bytes_touched={} batch_files_opened={} index_gets={} misses={}",
            compression.as_str(),
            result.selected_rows,
            result.positions_found,
            result.lookup_elapsed.as_secs_f64(),
            result.checksum,
            result.estimated_bytes_touched,
            result.batch_files_opened,
            result.index_gets,
            result.missing_positions
        );
        arrow_results.push(result);
    }

    let report = render_report(
        &cache_root,
        &vcf_path,
        &dataset_path,
        &arrow_store_root,
        &workload,
        &gate,
        &projection,
        &omitted_projection,
        &build_stats,
        &lance_results,
        &arrow_results,
    );
    print!("{report}");
    if let Some(path) = report_path {
        if let Some(parent) = path.parent() {
            fs::create_dir_all(parent)?;
        }
        fs::write(&path, report)?;
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
        append_position_keys(&batch, &mut positions)?;
    }
    Ok(positions)
}

#[allow(clippy::too_many_arguments)]
fn build_arrow_store(
    cache_root: &Path,
    chrom: &str,
    projection: &[String],
    store_path: &Path,
    compression: ArrowStoreCompression,
    target_batch_rows: usize,
    index_path: &Path,
    build_index: bool,
    index_block_entries: usize,
) -> BenchResult<ArrowBuildStats> {
    let started = Instant::now();
    if store_path.exists() {
        fs::remove_dir_all(store_path)?;
    }
    fs::create_dir_all(store_path.join("batches"))?;
    if build_index {
        if let Some(parent) = index_path.parent() {
            fs::create_dir_all(parent)?;
        }
        if index_path.exists() {
            fs::remove_file(index_path)?;
        }
    } else if !index_path.exists() {
        return Err(format!(
            "shared packed range index does not exist: {}",
            index_path.display()
        )
        .into());
    }

    let cold_files = cold_variation_files_for_chrom(cache_root.join("variation"), chrom);
    if cold_files.is_empty() {
        return Err(format!(
            "no cold variation parquet files found for {chrom} under {}",
            cache_root.join("variation").display()
        )
        .into());
    }

    let index_builder = if build_index {
        Some(PackedRangeIndexBuilder::create(
            index_path,
            index_block_entries,
            DEFAULT_INDEX_ZSTD_LEVEL,
        )?)
    } else {
        None
    };
    let mut builder = ArrowStoreBuilder::new(
        store_path.to_path_buf(),
        index_path.to_path_buf(),
        projection.to_vec(),
        compression,
        target_batch_rows,
        index_block_entries,
        index_builder,
    );

    for path in cold_files {
        let file = File::open(&path)?;
        let parquet_builder = ParquetRecordBatchReaderBuilder::try_new(file)?;
        let mask = projection_mask_for_columns(
            parquet_builder.schema().as_ref(),
            parquet_builder.parquet_schema(),
            projection,
        )?;
        let reader = parquet_builder
            .with_projection(mask)
            .with_batch_size(target_batch_rows.max(1))
            .build()?;
        for batch in reader {
            builder.append_batch(batch?)?;
        }
    }
    let mut manifest = builder.finish(chrom)?;

    manifest.index_size = dir_size(index_path)?;
    manifest.total_size = dir_size(store_path)? + manifest.index_size;
    fs::write(
        store_path.join("manifest.json"),
        serde_json::to_string_pretty(&manifest)?,
    )?;

    Ok(ArrowBuildStats {
        compression,
        elapsed: started.elapsed(),
        manifest,
    })
}

struct ArrowStoreBuilder {
    store_path: PathBuf,
    index_path: PathBuf,
    projection: Vec<String>,
    compression: ArrowStoreCompression,
    target_batch_rows: usize,
    index_block_entries: usize,
    index_builder: Option<PackedRangeIndexBuilder>,
    current_key: Option<u64>,
    current_segments: Vec<RecordBatch>,
    current_len: usize,
    pending_segments: Vec<RecordBatch>,
    pending_rows: usize,
    batch_id: u32,
    total_rows: usize,
    unique_positions: usize,
    max_position_group_len: usize,
    data_file_size: u64,
    schema: Option<SchemaRef>,
    previous_key: Option<u64>,
}

impl ArrowStoreBuilder {
    fn new(
        store_path: PathBuf,
        index_path: PathBuf,
        projection: Vec<String>,
        compression: ArrowStoreCompression,
        target_batch_rows: usize,
        index_block_entries: usize,
        index_builder: Option<PackedRangeIndexBuilder>,
    ) -> Self {
        Self {
            store_path,
            index_path,
            projection,
            compression,
            target_batch_rows: target_batch_rows.max(1),
            index_block_entries,
            index_builder,
            current_key: None,
            current_segments: Vec::new(),
            current_len: 0,
            pending_segments: Vec::new(),
            pending_rows: 0,
            batch_id: 0,
            total_rows: 0,
            unique_positions: 0,
            max_position_group_len: 0,
            data_file_size: 0,
            schema: None,
            previous_key: None,
        }
    }

    fn append_batch(&mut self, batch: RecordBatch) -> BenchResult<()> {
        if batch.num_rows() == 0 {
            return Ok(());
        }
        let schema = batch.schema();
        match self.schema.as_ref() {
            Some(existing) if existing.as_ref() != schema.as_ref() => {
                return Err("projected parquet schemas do not match".into());
            }
            None => self.schema = Some(schema.clone()),
            _ => {}
        }

        let position_idx = schema.index_of("position_key")?;
        let positions = batch
            .column(position_idx)
            .as_any()
            .downcast_ref::<UInt64Array>()
            .ok_or("position_key must be UInt64")?;
        let mut row = 0usize;
        while row < batch.num_rows() {
            if positions.is_null(row) {
                return Err("position_key must be non-null".into());
            }
            let key = positions.value(row);
            if self.previous_key.is_some_and(|previous| key < previous) {
                return Err(format!(
                    "rows out of position_key order: previous={} current={}",
                    self.previous_key.unwrap(),
                    key
                )
                .into());
            }
            let mut end = row + 1;
            while end < batch.num_rows() && !positions.is_null(end) && positions.value(end) == key {
                end += 1;
            }
            self.append_run(key, batch.slice(row, end - row))?;
            self.previous_key = Some(key);
            row = end;
        }
        Ok(())
    }

    fn append_run(&mut self, key: u64, segment: RecordBatch) -> BenchResult<()> {
        match self.current_key {
            Some(current) if current == key => {
                self.current_len += segment.num_rows();
                self.current_segments.push(segment);
            }
            Some(_) => {
                self.finish_current_group()?;
                self.current_key = Some(key);
                self.current_len = segment.num_rows();
                self.current_segments.push(segment);
            }
            None => {
                self.current_key = Some(key);
                self.current_len = segment.num_rows();
                self.current_segments.push(segment);
            }
        }
        Ok(())
    }

    fn finish_current_group(&mut self) -> BenchResult<()> {
        let Some(key) = self.current_key.take() else {
            return Ok(());
        };
        let group_len = self.current_len;
        if self.pending_rows > 0 && self.pending_rows + group_len > self.target_batch_rows {
            self.flush_pending()?;
        }
        let start_row = self.pending_rows;
        if let Some(index_builder) = self.index_builder.as_mut() {
            index_builder.append(
                key,
                self.batch_id,
                u32::try_from(start_row)?,
                u16::try_from(group_len)?,
            )?;
        }
        self.pending_rows += group_len;
        self.pending_segments.append(&mut self.current_segments);
        self.total_rows += group_len;
        self.unique_positions += 1;
        self.max_position_group_len = self.max_position_group_len.max(group_len);
        self.current_len = 0;
        if self.pending_rows >= self.target_batch_rows {
            self.flush_pending()?;
        }
        Ok(())
    }

    fn flush_pending(&mut self) -> BenchResult<()> {
        if self.pending_segments.is_empty() {
            return Ok(());
        }
        let schema = self
            .schema
            .as_ref()
            .ok_or("cannot flush Arrow store without schema")?
            .clone();
        let batch = concat_batches(&schema, &self.pending_segments)?;
        let path = self
            .store_path
            .join("batches")
            .join(format!("batch_{:06}.arrow", self.batch_id));
        self.data_file_size += write_arrow_batch_file(&path, &batch, self.compression)?;
        self.batch_id += 1;
        self.pending_segments.clear();
        self.pending_rows = 0;
        Ok(())
    }

    fn finish(mut self, chrom: &str) -> BenchResult<ArrowMmapManifest> {
        self.finish_current_group()?;
        self.flush_pending()?;
        if let Some(index_builder) = self.index_builder.take() {
            index_builder.finish()?;
        }
        Ok(ArrowMmapManifest {
            format: FORMAT_NAME.to_string(),
            version: FORMAT_VERSION,
            chrom: chrom.to_string(),
            tier: "cold".to_string(),
            compression: self.compression.as_str().to_string(),
            projected_columns: self.projection,
            target_batch_rows: self.target_batch_rows,
            batch_count: self.batch_id as usize,
            total_rows: self.total_rows,
            unique_positions: self.unique_positions,
            max_position_group_len: self.max_position_group_len,
            data_file_size: self.data_file_size,
            index_kind: INDEX_KIND_PACKED_ZSTD.to_string(),
            index_path: self.index_path.display().to_string(),
            index_block_entries: self.index_block_entries,
            index_size: dir_size(&self.index_path)?,
            total_size: 0,
        })
    }
}

fn write_arrow_batch_file(
    path: &Path,
    batch: &RecordBatch,
    compression: ArrowStoreCompression,
) -> BenchResult<u64> {
    let mut file = File::create(path)?;
    let mut options = IpcWriteOptions::try_new(8, false, MetadataVersion::V5)?;
    options = options.try_with_compression(compression.ipc_compression())?;
    {
        let mut writer =
            IpcFileWriter::try_new_with_options(&mut file, batch.schema().as_ref(), options)?;
        writer.write(batch)?;
        writer.finish()?;
    }
    Ok(file.metadata()?.len())
}

#[derive(Debug, Clone, Copy)]
struct PackedRangeEntry {
    position_key: u64,
    batch_id: u32,
    start_row: u32,
    len: u16,
}

#[derive(Debug, Clone)]
struct PackedBlockMeta {
    first_key: u64,
    last_key: u64,
    offset: u64,
    compressed_len: u32,
    uncompressed_len: u32,
    count: u32,
}

struct PackedRangeIndexBuilder {
    path: PathBuf,
    temp_path: PathBuf,
    temp_file: File,
    block_entries: usize,
    zstd_level: i32,
    entries: Vec<PackedRangeEntry>,
    blocks: Vec<PackedBlockMeta>,
    entry_count: u64,
    data_len: u64,
    previous_key: Option<u64>,
}

impl PackedRangeIndexBuilder {
    fn create(path: &Path, block_entries: usize, zstd_level: i32) -> BenchResult<Self> {
        let block_entries = block_entries.max(1);
        let temp_path = path.with_extension("pzstd.tmp");
        if temp_path.exists() {
            fs::remove_file(&temp_path)?;
        }
        let temp_file = File::create(&temp_path)?;
        Ok(Self {
            path: path.to_path_buf(),
            temp_path,
            temp_file,
            block_entries,
            zstd_level,
            entries: Vec::with_capacity(block_entries),
            blocks: Vec::new(),
            entry_count: 0,
            data_len: 0,
            previous_key: None,
        })
    }

    fn append(
        &mut self,
        position_key: u64,
        batch_id: u32,
        start_row: u32,
        len: u16,
    ) -> BenchResult<()> {
        if let Some(previous_key) = self.previous_key
            && position_key <= previous_key
        {
            return Err(format!(
                "packed range index keys must be strictly increasing: previous={} current={}",
                previous_key, position_key
            )
            .into());
        }
        self.entries.push(PackedRangeEntry {
            position_key,
            batch_id,
            start_row,
            len,
        });
        self.previous_key = Some(position_key);
        self.entry_count += 1;
        if self.entries.len() >= self.block_entries {
            self.flush_block()?;
        }
        Ok(())
    }

    fn flush_block(&mut self) -> BenchResult<()> {
        if self.entries.is_empty() {
            return Ok(());
        }
        let raw = encode_packed_block(&self.entries)?;
        let compressed = zstd::bulk::compress(&raw, self.zstd_level)?;
        let first_key = self
            .entries
            .first()
            .expect("checked non-empty")
            .position_key;
        let last_key = self.entries.last().expect("checked non-empty").position_key;
        let offset = self.data_len;
        self.temp_file.write_all(&compressed)?;
        self.data_len += u64::try_from(compressed.len())?;
        self.blocks.push(PackedBlockMeta {
            first_key,
            last_key,
            offset,
            compressed_len: u32::try_from(compressed.len())?,
            uncompressed_len: u32::try_from(raw.len())?,
            count: u32::try_from(self.entries.len())?,
        });
        self.entries.clear();
        Ok(())
    }

    fn finish(mut self) -> BenchResult<()> {
        self.flush_block()?;
        self.temp_file.flush()?;
        drop(self.temp_file);

        let mut out = File::create(&self.path)?;
        out.write_all(PACKED_INDEX_MAGIC)?;
        write_u32(&mut out, u32::try_from(self.block_entries)?)?;
        write_i32(&mut out, self.zstd_level)?;
        write_u64(&mut out, self.entry_count)?;
        write_u64(&mut out, u64::try_from(self.blocks.len())?)?;
        let metadata_len = u64::try_from(self.blocks.len())? * 36;
        write_u64(&mut out, metadata_len)?;
        write_u64(&mut out, self.data_len)?;
        for block in &self.blocks {
            write_u64(&mut out, block.first_key)?;
            write_u64(&mut out, block.last_key)?;
            write_u64(&mut out, block.offset)?;
            write_u32(&mut out, block.compressed_len)?;
            write_u32(&mut out, block.uncompressed_len)?;
            write_u32(&mut out, block.count)?;
        }
        let mut temp = File::open(&self.temp_path)?;
        std::io::copy(&mut temp, &mut out)?;
        out.flush()?;
        fs::remove_file(&self.temp_path)?;
        Ok(())
    }
}

struct DecodedPackedBlock {
    base_key: u64,
    key_deltas: Vec<u32>,
    batch_ids: Vec<u32>,
    start_rows: Vec<u32>,
    lens: Vec<u16>,
}

struct PackedRangeIndexReader {
    _mmap: Mmap,
    blocks: Vec<PackedBlockMeta>,
    data_start: usize,
    decoded_block: Option<(usize, DecodedPackedBlock)>,
}

impl PackedRangeIndexReader {
    fn open(path: &Path, _target_batch_rows: usize) -> BenchResult<Self> {
        let file = File::open(path)?;
        let mmap = unsafe { Mmap::map(&file) }?;
        let mut cursor = Cursor::new(&mmap[..]);
        let mut magic = [0_u8; 8];
        cursor.read_exact(&mut magic)?;
        if &magic != PACKED_INDEX_MAGIC {
            return Err(format!("invalid packed range index magic in {}", path.display()).into());
        }
        let _block_entries = read_u32(&mut cursor)?;
        let _zstd_level = read_i32(&mut cursor)?;
        let _entry_count = read_u64(&mut cursor)?;
        let block_count = usize::try_from(read_u64(&mut cursor)?)?;
        let metadata_len = usize::try_from(read_u64(&mut cursor)?)?;
        let _data_len = usize::try_from(read_u64(&mut cursor)?)?;
        let expected_metadata_len = block_count
            .checked_mul(36)
            .ok_or("packed range index metadata length overflow")?;
        if metadata_len != expected_metadata_len {
            return Err(format!(
                "packed range index metadata length mismatch: expected={} actual={}",
                expected_metadata_len, metadata_len
            )
            .into());
        }
        let mut blocks = Vec::with_capacity(block_count);
        for _ in 0..block_count {
            blocks.push(PackedBlockMeta {
                first_key: read_u64(&mut cursor)?,
                last_key: read_u64(&mut cursor)?,
                offset: read_u64(&mut cursor)?,
                compressed_len: read_u32(&mut cursor)?,
                uncompressed_len: read_u32(&mut cursor)?,
                count: read_u32(&mut cursor)?,
            });
        }
        let data_start = usize::try_from(cursor.position())?;
        if data_start + metadata_len.saturating_sub(metadata_len) > mmap.len() {
            return Err(format!("packed range index {} is truncated", path.display()).into());
        }
        Ok(Self {
            _mmap: mmap,
            blocks,
            data_start,
            decoded_block: None,
        })
    }

    fn lookup_position(&mut self, position_key: u64) -> BenchResult<Option<PositionRange>> {
        let Some(block_id) = self.find_block(position_key) else {
            return Ok(None);
        };
        if self
            .decoded_block
            .as_ref()
            .is_none_or(|(decoded_id, _)| *decoded_id != block_id)
        {
            let decoded = self.decode_block(block_id)?;
            self.decoded_block = Some((block_id, decoded));
        }
        let decoded = &self
            .decoded_block
            .as_ref()
            .expect("decoded block inserted above")
            .1;
        let Some(delta) = position_key.checked_sub(decoded.base_key) else {
            return Ok(None);
        };
        let Ok(delta) = u32::try_from(delta) else {
            return Ok(None);
        };
        let Ok(offset) = decoded.key_deltas.binary_search(&delta) else {
            return Ok(None);
        };
        Ok(Some(PositionRange {
            batch_id: decoded.batch_ids[offset],
            start_row: decoded.start_rows[offset],
            len: u32::from(decoded.lens[offset]),
        }))
    }

    fn find_block(&self, position_key: u64) -> Option<usize> {
        let block_id = self
            .blocks
            .partition_point(|block| block.last_key < position_key);
        let block = self.blocks.get(block_id)?;
        (position_key >= block.first_key).then_some(block_id)
    }

    fn decode_block(&self, block_id: usize) -> BenchResult<DecodedPackedBlock> {
        let block = self
            .blocks
            .get(block_id)
            .ok_or_else(|| format!("packed range block out of bounds: {block_id}"))?;
        let start = self.data_start + usize::try_from(block.offset)?;
        let end = start + usize::try_from(block.compressed_len)?;
        let compressed = self
            ._mmap
            .get(start..end)
            .ok_or("packed range index block extends past file end")?;
        let raw = zstd::bulk::decompress(compressed, usize::try_from(block.uncompressed_len)?)?;
        decode_packed_block(&raw, usize::try_from(block.count)?)
    }
}

fn encode_packed_block(entries: &[PackedRangeEntry]) -> BenchResult<Vec<u8>> {
    let first_key = entries
        .first()
        .ok_or("cannot encode empty packed range block")?
        .position_key;
    let mut raw = Vec::with_capacity(12 + entries.len() * 14);
    raw.extend_from_slice(&u32::try_from(entries.len())?.to_le_bytes());
    raw.extend_from_slice(&first_key.to_le_bytes());
    for entry in entries {
        let delta = entry
            .position_key
            .checked_sub(first_key)
            .ok_or("packed range block keys are not sorted")?;
        raw.extend_from_slice(&u32::try_from(delta)?.to_le_bytes());
    }
    for entry in entries {
        raw.extend_from_slice(&entry.batch_id.to_le_bytes());
    }
    for entry in entries {
        raw.extend_from_slice(&entry.start_row.to_le_bytes());
    }
    for entry in entries {
        raw.extend_from_slice(&entry.len.to_le_bytes());
    }
    Ok(raw)
}

fn decode_packed_block(raw: &[u8], expected_count: usize) -> BenchResult<DecodedPackedBlock> {
    let mut cursor = Cursor::new(raw);
    let count = usize::try_from(read_u32(&mut cursor)?)?;
    if count != expected_count {
        return Err(format!(
            "packed range block count mismatch: expected={} actual={}",
            expected_count, count
        )
        .into());
    }
    let base_key = read_u64(&mut cursor)?;
    let mut key_deltas = Vec::with_capacity(count);
    for _ in 0..count {
        key_deltas.push(read_u32(&mut cursor)?);
    }
    let mut batch_ids = Vec::with_capacity(count);
    for _ in 0..count {
        batch_ids.push(read_u32(&mut cursor)?);
    }
    let mut start_rows = Vec::with_capacity(count);
    for _ in 0..count {
        start_rows.push(read_u32(&mut cursor)?);
    }
    let mut lens = Vec::with_capacity(count);
    for _ in 0..count {
        lens.push(read_u16(&mut cursor)?);
    }
    if usize::try_from(cursor.position())? != raw.len() {
        return Err("packed range block has trailing bytes".into());
    }
    Ok(DecodedPackedBlock {
        base_key,
        key_deltas,
        batch_ids,
        start_rows,
        lens,
    })
}

fn write_u32(writer: &mut impl Write, value: u32) -> BenchResult<()> {
    Ok(writer.write_all(&value.to_le_bytes())?)
}

fn write_i32(writer: &mut impl Write, value: i32) -> BenchResult<()> {
    Ok(writer.write_all(&value.to_le_bytes())?)
}

fn write_u64(writer: &mut impl Write, value: u64) -> BenchResult<()> {
    Ok(writer.write_all(&value.to_le_bytes())?)
}

fn read_u16(reader: &mut impl Read) -> BenchResult<u16> {
    let mut bytes = [0_u8; 2];
    reader.read_exact(&mut bytes)?;
    Ok(u16::from_le_bytes(bytes))
}

fn read_u32(reader: &mut impl Read) -> BenchResult<u32> {
    let mut bytes = [0_u8; 4];
    reader.read_exact(&mut bytes)?;
    Ok(u32::from_le_bytes(bytes))
}

fn read_i32(reader: &mut impl Read) -> BenchResult<i32> {
    let mut bytes = [0_u8; 4];
    reader.read_exact(&mut bytes)?;
    Ok(i32::from_le_bytes(bytes))
}

fn read_u64(reader: &mut impl Read) -> BenchResult<u64> {
    let mut bytes = [0_u8; 8];
    reader.read_exact(&mut bytes)?;
    Ok(u64::from_le_bytes(bytes))
}

impl ArrowMmapStore {
    fn open(root: &Path) -> BenchResult<Self> {
        let manifest = read_manifest(root)?;
        if manifest.format != FORMAT_NAME || manifest.version != FORMAT_VERSION {
            return Err(format!(
                "unsupported Arrow mmap store format={} version={}",
                manifest.format, manifest.version
            )
            .into());
        }
        if manifest.index_kind != INDEX_KIND_PACKED_ZSTD {
            return Err(format!("unsupported range index kind {}", manifest.index_kind).into());
        }
        let index_path = PathBuf::from(&manifest.index_path);
        let index = PackedRangeIndexReader::open(&index_path, manifest.target_batch_rows)?;
        Ok(Self {
            root: root.to_path_buf(),
            manifest,
            index,
            batch_cache: HashMap::new(),
            batch_files_opened: 0,
            index_gets: 0,
            missing_positions: 0,
        })
    }

    fn run_lookup(
        &mut self,
        compression: ArrowStoreCompression,
        positions: &[u64],
        projected_columns: usize,
        open_elapsed: Duration,
    ) -> BenchResult<ArrowLookupResult> {
        let started = Instant::now();
        let mut ranges = Vec::new();
        let mut positions_found = 0usize;
        for &position_key in positions {
            self.index_gets += 1;
            let Some(range) = self
                .index
                .lookup_position(position_key)
                .map_err(|err| format!("failed to lookup position_key {position_key}: {err}"))?
            else {
                self.missing_positions += 1;
                continue;
            };
            positions_found += 1;
            ranges.push((position_key, range));
        }

        let mut acc = TouchAccumulator::default();
        for (position_key, range) in ranges {
            let batch = self.load_batch(range.batch_id).map_err(|err| {
                format!(
                    "failed to load batch {} for position_key {}: {err}",
                    range.batch_id, position_key
                )
            })?;
            let start = usize::try_from(range.start_row)?;
            let len = usize::try_from(range.len)?;
            if start + len > batch.num_rows() {
                return Err(format!(
                    "range out of bounds batch={} start={} len={} rows={}",
                    range.batch_id,
                    start,
                    len,
                    batch.num_rows()
                )
                .into());
            }
            for row in start..start + len {
                for (field, column) in batch.schema().fields().iter().zip(batch.columns()) {
                    touch_array_value(column.as_ref(), field.data_type(), row, &mut acc)?;
                }
                acc.rows += 1;
            }
        }
        Ok(ArrowLookupResult {
            compression,
            total_size: self.manifest.total_size,
            data_size: self.manifest.data_file_size,
            index_size: self.manifest.index_size,
            open_elapsed,
            lookup_elapsed: started.elapsed(),
            positions_requested: positions.len(),
            positions_found,
            selected_rows: acc.rows,
            projected_columns,
            estimated_bytes_touched: acc.estimated_bytes_touched,
            batch_files_opened: self.batch_files_opened,
            index_gets: self.index_gets,
            missing_positions: self.missing_positions,
            checksum: acc.checksum,
        })
    }

    fn load_batch(&mut self, batch_id: u32) -> BenchResult<&RecordBatch> {
        if !self.batch_cache.contains_key(&batch_id) {
            let path = self
                .root
                .join("batches")
                .join(format!("batch_{batch_id:06}.arrow"));
            let file = File::open(&path)?;
            let file_len = file.metadata()?.len();
            if file_len == 0 {
                return Err(format!("Arrow IPC batch file is empty: {}", path.display()).into());
            }
            let mmap = unsafe { Mmap::map(&file) }
                .map_err(|err| format!("failed to mmap {}: {err}", path.display()))?;
            let cursor = Cursor::new(&mmap[..]);
            let mut reader = IpcFileReader::try_new(cursor, None)
                .map_err(|err| format!("failed to open Arrow IPC {}: {err}", path.display()))?;
            let batch = reader.next().ok_or_else(|| {
                format!(
                    "Arrow IPC batch file {} did not contain a RecordBatch",
                    path.display()
                )
            })??;
            self.batch_files_opened += 1;
            self.batch_cache
                .insert(batch_id, CachedBatch { _mmap: mmap, batch });
        }
        Ok(&self
            .batch_cache
            .get(&batch_id)
            .expect("batch inserted above")
            .batch)
    }
}

fn touch_array_value(
    array: &dyn Array,
    data_type: &DataType,
    row: usize,
    acc: &mut TouchAccumulator,
) -> BenchResult<()> {
    if array.is_null(row) {
        acc.mix_u64(0xFFFF_FFFF_FFFF_FFFF);
        acc.estimated_bytes_touched += 1;
        return Ok(());
    }
    match data_type {
        DataType::Null => {
            let _ = array.as_any().downcast_ref::<NullArray>();
            acc.mix_u64(0);
        }
        DataType::Boolean => {
            let values = downcast::<BooleanArray>(array, data_type)?;
            acc.mix_u64(u64::from(values.value(row)));
            acc.estimated_bytes_touched += 1;
        }
        DataType::Int8 => {
            touch_primitive::<Int8Array, _>(array, data_type, row, acc, |v| v as i64 as u64)?
        }
        DataType::Int16 => {
            touch_primitive::<Int16Array, _>(array, data_type, row, acc, |v| v as i64 as u64)?
        }
        DataType::Int32 => {
            touch_primitive::<Int32Array, _>(array, data_type, row, acc, |v| v as i64 as u64)?
        }
        DataType::Int64 => {
            touch_primitive::<Int64Array, _>(array, data_type, row, acc, |v| v as u64)?
        }
        DataType::UInt8 => touch_primitive::<UInt8Array, _>(array, data_type, row, acc, u64::from)?,
        DataType::UInt16 => {
            touch_primitive::<UInt16Array, _>(array, data_type, row, acc, u64::from)?
        }
        DataType::UInt32 => {
            touch_primitive::<UInt32Array, _>(array, data_type, row, acc, u64::from)?
        }
        DataType::UInt64 => touch_primitive::<UInt64Array, _>(array, data_type, row, acc, |v| v)?,
        DataType::Float32 => {
            touch_primitive::<Float32Array, _>(array, data_type, row, acc, |v| v.to_bits() as u64)?
        }
        DataType::Float64 => {
            touch_primitive::<Float64Array, _>(array, data_type, row, acc, |v| v.to_bits())?
        }
        DataType::Utf8 => {
            let values = downcast::<StringArray>(array, data_type)?;
            acc.touch_bytes(values.value(row).as_bytes());
        }
        DataType::LargeUtf8 => {
            let values = downcast::<LargeStringArray>(array, data_type)?;
            acc.touch_bytes(values.value(row).as_bytes());
        }
        DataType::Binary => {
            let values = downcast::<BinaryArray>(array, data_type)?;
            acc.touch_bytes(values.value(row));
        }
        DataType::LargeBinary => {
            let values = downcast::<LargeBinaryArray>(array, data_type)?;
            acc.touch_bytes(values.value(row));
        }
        DataType::FixedSizeBinary(_) => {
            let values = downcast::<FixedSizeBinaryArray>(array, data_type)?;
            acc.touch_bytes(values.value(row));
        }
        DataType::List(field) => {
            let values = downcast::<ListArray>(array, data_type)?;
            touch_list_value::<i32>(values, field, row, acc)?;
        }
        DataType::LargeList(field) => {
            let values = downcast::<LargeListArray>(array, data_type)?;
            touch_list_value::<i64>(values, field, row, acc)?;
        }
        DataType::Struct(fields) => {
            let values = downcast::<StructArray>(array, data_type)?;
            touch_struct_value(values, fields, row, acc)?;
        }
        DataType::Dictionary(key_type, value_type) => {
            touch_dictionary_value(array, key_type.as_ref(), value_type.as_ref(), row, acc)?;
        }
        other => {
            let sliced = array.slice(row, 1);
            let data = sliced.to_data();
            acc.mix_u64(row as u64);
            acc.mix_u64(data.len() as u64);
            for buffer in data.buffers() {
                acc.touch_bytes(buffer.as_slice());
            }
            for child in data.child_data() {
                let child_array = make_array(child.clone());
                touch_array_value(child_array.as_ref(), other, 0, acc)?;
            }
        }
    }
    Ok(())
}

fn touch_primitive<A, F>(
    array: &dyn Array,
    data_type: &DataType,
    row: usize,
    acc: &mut TouchAccumulator,
    convert: F,
) -> BenchResult<()>
where
    A: Array + ArrowPrimitiveValue + 'static,
    F: Fn(<A as ArrowPrimitiveValue>::Native) -> u64,
{
    let values = downcast::<A>(array, data_type)?;
    let value = primitive_value(values, row)?;
    acc.mix_u64(convert(value));
    acc.estimated_bytes_touched += std::mem::size_of_val(&value) as u64;
    Ok(())
}

trait ArrowPrimitiveValue: Array {
    type Native: Copy;
    fn primitive_value(&self, row: usize) -> Self::Native;
}

macro_rules! impl_primitive_value {
    ($array:ty, $native:ty) => {
        impl ArrowPrimitiveValue for $array {
            type Native = $native;
            fn primitive_value(&self, row: usize) -> Self::Native {
                self.value(row)
            }
        }
    };
}

impl_primitive_value!(Int8Array, i8);
impl_primitive_value!(Int16Array, i16);
impl_primitive_value!(Int32Array, i32);
impl_primitive_value!(Int64Array, i64);
impl_primitive_value!(UInt8Array, u8);
impl_primitive_value!(UInt16Array, u16);
impl_primitive_value!(UInt32Array, u32);
impl_primitive_value!(UInt64Array, u64);
impl_primitive_value!(Float32Array, f32);
impl_primitive_value!(Float64Array, f64);

fn primitive_value<A: ArrowPrimitiveValue>(array: &A, row: usize) -> BenchResult<A::Native> {
    Ok(array.primitive_value(row))
}

fn touch_list_value<O>(
    array: &GenericListArray<O>,
    field: &Arc<Field>,
    row: usize,
    acc: &mut TouchAccumulator,
) -> BenchResult<()>
where
    O: datafusion::arrow::array::OffsetSizeTrait,
{
    let offsets = array.value_offsets();
    let start = offsets[row]
        .to_usize()
        .ok_or("list offset out of usize range")?;
    let end = offsets[row + 1]
        .to_usize()
        .ok_or("list offset out of usize range")?;
    acc.mix_u64((end - start) as u64);
    acc.estimated_bytes_touched += std::mem::size_of_val(&offsets[row]) as u64 * 2;
    let values = array.values();
    for child_row in start..end {
        touch_array_value(values.as_ref(), field.data_type(), child_row, acc)?;
    }
    Ok(())
}

fn touch_struct_value(
    array: &StructArray,
    fields: &Fields,
    row: usize,
    acc: &mut TouchAccumulator,
) -> BenchResult<()> {
    for (field, column) in fields.iter().zip(array.columns()) {
        touch_array_value(column.as_ref(), field.data_type(), row, acc)?;
    }
    Ok(())
}

fn touch_dictionary_value(
    array: &dyn Array,
    key_type: &DataType,
    value_type: &DataType,
    row: usize,
    acc: &mut TouchAccumulator,
) -> BenchResult<()> {
    macro_rules! touch_dict {
        ($key_ty:ty) => {{
            let dict = downcast::<DictionaryArray<$key_ty>>(
                array,
                &DataType::Dictionary(Box::new(key_type.clone()), Box::new(value_type.clone())),
            )?;
            let key = dict
                .keys()
                .value(row)
                .to_usize()
                .ok_or("dictionary key out of usize range")?;
            acc.mix_u64(key as u64);
            touch_array_value(dict.values().as_ref(), value_type, key, acc)
        }};
    }
    match key_type {
        DataType::Int8 => touch_dict!(datafusion::arrow::datatypes::Int8Type),
        DataType::Int16 => touch_dict!(datafusion::arrow::datatypes::Int16Type),
        DataType::Int32 => touch_dict!(datafusion::arrow::datatypes::Int32Type),
        DataType::Int64 => touch_dict!(datafusion::arrow::datatypes::Int64Type),
        DataType::UInt8 => touch_dict!(datafusion::arrow::datatypes::UInt8Type),
        DataType::UInt16 => touch_dict!(datafusion::arrow::datatypes::UInt16Type),
        DataType::UInt32 => touch_dict!(datafusion::arrow::datatypes::UInt32Type),
        DataType::UInt64 => touch_dict!(datafusion::arrow::datatypes::UInt64Type),
        other => Err(format!("unsupported dictionary key type: {other:?}").into()),
    }
}

fn downcast<'a, A: Array + 'static>(
    array: &'a dyn Array,
    data_type: &DataType,
) -> BenchResult<&'a A> {
    array
        .as_any()
        .downcast_ref::<A>()
        .ok_or_else(|| format!("array downcast failed for {data_type:?}").into())
}

async fn scan_lance_plan(
    dataset: &Dataset,
    name: &str,
    positions: &[u64],
    query_batch_size: usize,
    scan_batch_size: usize,
    projection: &[String],
    materialization_style: MaterializationStyle,
) -> BenchResult<LancePlanResult> {
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

    Ok(LancePlanResult {
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

fn projected_cold_everything_columns(
    dataset: &Dataset,
    variation_root: &Path,
    chrom: &str,
) -> BenchResult<(Vec<String>, Vec<String>)> {
    let lance_available = dataset
        .schema()
        .fields
        .iter()
        .map(|field| field.name.clone())
        .collect::<HashSet<_>>();
    let cold_files = cold_variation_files_for_chrom(variation_root, chrom);
    let first_cold_file = cold_files
        .first()
        .ok_or_else(|| format!("no cold variation parquet files found for {chrom}"))?;
    let cold_file = File::open(first_cold_file)?;
    let cold_schema = ParquetRecordBatchReaderBuilder::try_new(cold_file)?
        .schema()
        .clone();
    let cold_available = cold_schema
        .fields()
        .iter()
        .map(|field| field.name().to_string())
        .collect::<HashSet<_>>();

    let mut preferred = vec!["consequence_types", "most_severe_consequence"];
    for column in cache_lookup_column_names() {
        if !preferred.contains(&column) {
            preferred.push(column);
        }
    }
    let cache_columns = preferred
        .into_iter()
        .filter(|column| lance_available.contains(*column) && cold_available.contains(*column))
        .map(str::to_string)
        .collect::<Vec<_>>();
    let mut projection = Vec::new();
    let mut omitted = Vec::new();
    for column in projection_columns_for_cache(&cache_columns, true) {
        if lance_available.contains(&column) && cold_available.contains(&column) {
            projection.push(column);
        } else if lance_available.contains(&column) && !cold_available.contains(&column) {
            push_unique(&mut omitted, column);
        }
    }
    Ok((projection, omitted))
}

#[allow(clippy::too_many_arguments)]
fn render_report(
    cache_root: &Path,
    vcf_path: &Path,
    lance_dataset_path: &Path,
    arrow_store_root: &Path,
    workload: &Workload,
    gate: &GateStats,
    projection: &[String],
    omitted_projection: &[String],
    build_stats: &[ArrowBuildStats],
    lance_results: &[LancePlanResult],
    arrow_results: &[ArrowLookupResult],
) -> String {
    let mut out = String::new();
    let _ = writeln!(out, "# Arrow Mmap Cold Lookup Benchmark");
    let _ = writeln!(out);
    let _ = writeln!(out, "## Workload");
    let _ = writeln!(out);
    let _ = writeln!(out, "- Cache root: `{}`", cache_root.display());
    let _ = writeln!(out, "- VCF: `{}`", vcf_path.display());
    let _ = writeln!(out, "- Lance dataset: `{}`", lance_dataset_path.display());
    let _ = writeln!(out, "- Arrow store root: `{}`", arrow_store_root.display());
    let _ = writeln!(out, "- Chrom: `{}`", workload.chrom);
    let _ = writeln!(out, "- Parsed variants: `{}`", workload.variants);
    let _ = writeln!(out, "- Raw probe attempts: `{}`", workload.probe_attempts);
    let _ = writeln!(
        out,
        "- Warm-miss unique positions: `{}`",
        gate.warm_miss_unique
    );
    let _ = writeln!(
        out,
        "- Position-index admitted unique positions: `{}`",
        gate.position_index_unique
    );
    let _ = writeln!(
        out,
        "- Bloom-admitted unique positions: `{}`",
        gate.bloom_unique
    );
    let _ = writeln!(out, "- Bloom-admitted attempts: `{}`", gate.bloom_attempts);
    let _ = writeln!(out, "- Projected columns: `{}`", projection.len());
    if let Some(stats) = build_stats.first() {
        let _ = writeln!(
            out,
            "- Row-location index kind: `{}`",
            stats.manifest.index_kind
        );
        let _ = writeln!(
            out,
            "- Row-location index path: `{}`",
            stats.manifest.index_path
        );
        let _ = writeln!(
            out,
            "- Row-location index block entries: `{}`",
            stats.manifest.index_block_entries
        );
    }
    let _ = writeln!(
        out,
        "- Omitted requested columns absent from cold parquet: `{}`",
        if omitted_projection.is_empty() {
            "none".to_string()
        } else {
            omitted_projection.join(",")
        }
    );
    let _ = writeln!(out);
    let _ = writeln!(out, "## Store Sizes");
    let _ = writeln!(out);
    let _ = writeln!(
        out,
        "Total size is data-store directory size plus the shared packed range index."
    );
    let _ = writeln!(out);
    let _ = writeln!(
        out,
        "| compression | build seconds | rows | unique positions | batches | data size | index size | total size |"
    );
    let _ = writeln!(out, "|---|---:|---:|---:|---:|---:|---:|---:|");
    for stats in build_stats {
        let _ = writeln!(
            out,
            "| {} | {:.3} | {} | {} | {} | {} | {} | {} |",
            stats.compression.as_str(),
            stats.elapsed.as_secs_f64(),
            stats.manifest.total_rows,
            stats.manifest.unique_positions,
            stats.manifest.batch_count,
            stats.manifest.data_file_size,
            stats.manifest.index_size,
            stats.manifest.total_size
        );
    }
    let _ = writeln!(out);
    let _ = writeln!(out, "## Lance Cold Lookup");
    let _ = writeln!(out);
    let _ = writeln!(
        out,
        "| plan | query batch | scans | filter keys | rows | result batches | selected positions | seconds | bytes read | requests | ranges | fragments | parts loaded | index comparisons |"
    );
    let _ = writeln!(
        out,
        "|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|"
    );
    for result in lance_results {
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
    let _ = writeln!(out);
    let _ = writeln!(out, "## Arrow Mmap Lookup");
    let _ = writeln!(out);
    let _ = writeln!(
        out,
        "| compression | mmap direct | total size | data size | index size | open seconds | lookup seconds | positions requested | positions found | selected rows | projected columns | estimated bytes touched | batch files opened | index gets | missing positions | checksum |"
    );
    let _ = writeln!(
        out,
        "|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|"
    );
    for result in arrow_results {
        let mmap_direct = result.compression == ArrowStoreCompression::Uncompressed;
        let _ = writeln!(
            out,
            "| {} | {} | {} | {} | {} | {:.3} | {:.3} | {} | {} | {} | {} | {} | {} | {} | {} | {} |",
            result.compression.as_str(),
            mmap_direct,
            result.total_size,
            result.data_size,
            result.index_size,
            result.open_elapsed.as_secs_f64(),
            result.lookup_elapsed.as_secs_f64(),
            result.positions_requested,
            result.positions_found,
            result.selected_rows,
            result.projected_columns,
            result.estimated_bytes_touched,
            result.batch_files_opened,
            result.index_gets,
            result.missing_positions,
            result.checksum
        );
    }
    let _ = writeln!(out);
    let _ = writeln!(out, "## Verification");
    let _ = writeln!(out);
    let arrow_rows_match = arrow_results.first().is_none_or(|first| {
        arrow_results
            .iter()
            .all(|row| row.selected_rows == first.selected_rows)
    });
    let arrow_positions_match = arrow_results.first().is_none_or(|first| {
        arrow_results
            .iter()
            .all(|row| row.positions_found == first.positions_found)
    });
    let arrow_checksums_match = arrow_results.first().is_none_or(|first| {
        arrow_results
            .iter()
            .all(|row| row.checksum == first.checksum)
    });
    let lance_rows = lance_results
        .iter()
        .find(|row| row.name == "large_chunks")
        .map(|row| row.rows);
    let arrow_rows = arrow_results.first().map(|row| row.selected_rows);
    let lance_positions = lance_results
        .iter()
        .find(|row| row.name == "large_chunks")
        .map(|row| row.selected_positions);
    let arrow_positions = arrow_results.first().map(|row| row.positions_found);
    let _ = writeln!(out, "- Arrow row counts match: `{arrow_rows_match}`");
    let _ = writeln!(
        out,
        "- Arrow found-position counts match: `{arrow_positions_match}`"
    );
    let _ = writeln!(out, "- Arrow checksums match: `{arrow_checksums_match}`");
    let _ = writeln!(
        out,
        "- Lance large-chunks rows vs Arrow rows: `{:?}` vs `{:?}`",
        lance_rows, arrow_rows
    );
    let _ = writeln!(
        out,
        "- Lance large-chunks positions vs Arrow found positions: `{:?}` vs `{:?}`",
        lance_positions, arrow_positions
    );
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

fn projection_mask_for_columns(
    arrow_schema: &Schema,
    parquet_schema: &parquet::schema::types::SchemaDescriptor,
    columns: &[String],
) -> BenchResult<ProjectionMask> {
    let mut roots = Vec::with_capacity(columns.len());
    for column in columns {
        roots.push(arrow_schema.index_of(column)?);
    }
    Ok(ProjectionMask::roots(parquet_schema, roots))
}

fn cold_variation_files_for_chrom(variation_dir: impl AsRef<Path>, chrom: &str) -> Vec<PathBuf> {
    let dir = variation_dir.as_ref();
    let stripped = chrom.strip_prefix("chr").unwrap_or(chrom);
    let mut candidates = vec![
        dir.join(format!("{chrom}_cold.parquet")),
        dir.join(format!("{stripped}_cold.parquet")),
        dir.join(format!("chr{stripped}_cold.parquet")),
    ];
    let mut files = Vec::new();
    for candidate in candidates.drain(..) {
        if candidate.is_file() {
            files.push(candidate);
        }
    }
    if let Ok(entries) = fs::read_dir(dir) {
        for entry in entries.flatten() {
            let path = entry.path();
            let Some(name) = path.file_name().and_then(|name| name.to_str()) else {
                continue;
            };
            let is_match = name.starts_with(&format!("{chrom}_cold_"))
                || name.starts_with(&format!("{stripped}_cold_"))
                || name.starts_with(&format!("chr{stripped}_cold_"));
            if is_match && name.ends_with(".parquet") && path.is_file() {
                files.push(path);
            }
        }
    }
    files.sort();
    files.dedup();
    files
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

fn append_position_keys(batch: &RecordBatch, out: &mut HashSet<u64>) -> BenchResult<()> {
    let index = batch.schema().index_of("position_key")?;
    let array = batch
        .column(index)
        .as_any()
        .downcast_ref::<UInt64Array>()
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
    start: i64,
    end: i64,
    reference: &str,
    alternate: &str,
    extended: bool,
) -> Vec<i64> {
    let mut starts = Vec::with_capacity(4);
    push_unique(&mut starts, start);
    if extended {
        push_unique(&mut starts, end);
        let (input_ref, input_alt, input_start) =
            vcf_to_vep_input_allele(start, reference, alternate);
        push_unique(&mut starts, input_start);
        let vep_start = vep_norm_start(start, &input_ref, &input_alt);
        push_unique(&mut starts, vep_start);
    }
    starts
}

fn vep_norm_start(start: i64, reference: &str, _alternate: &str) -> i64 {
    if reference == "-" || reference.is_empty() {
        start + 1
    } else {
        start
    }
}

fn push_unique<T: PartialEq>(values: &mut Vec<T>, value: T) {
    if !values.iter().any(|existing| existing == &value) {
        values.push(value);
    }
}

fn read_manifest(path: &Path) -> BenchResult<ArrowMmapManifest> {
    Ok(serde_json::from_slice(&fs::read(
        path.join("manifest.json"),
    )?)?)
}

fn dir_size(path: impl AsRef<Path>) -> BenchResult<u64> {
    let path = path.as_ref();
    if path.is_file() {
        return Ok(path.metadata()?.len());
    }
    let mut total = 0_u64;
    if !path.exists() {
        return Ok(0);
    }
    for entry in fs::read_dir(path)? {
        let entry = entry?;
        let metadata = entry.metadata()?;
        if metadata.is_dir() {
            total += dir_size(entry.path())?;
        } else {
            total += metadata.len();
        }
    }
    Ok(total)
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
        .map(|value| value != "0" && !value.eq_ignore_ascii_case("false"))
        .unwrap_or(default)
}

#[cfg(test)]
mod tests {
    use super::*;
    use datafusion::arrow::array::{ArrayRef, Int64Array, StringArray};
    use datafusion::arrow::datatypes::{Field, Schema};
    use tempfile::TempDir;

    #[test]
    fn builder_keeps_position_group_inside_one_batch_file() {
        let tmp = TempDir::new().unwrap();
        let schema = Arc::new(Schema::new(vec![
            Field::new("position_key", DataType::UInt64, false),
            Field::new("end", DataType::Int64, false),
            Field::new("allele_string", DataType::Utf8, false),
        ]));
        let batch = RecordBatch::try_new(
            schema,
            vec![
                Arc::new(UInt64Array::from(vec![10, 10, 20, 20, 20])) as ArrayRef,
                Arc::new(Int64Array::from(vec![1, 2, 3, 4, 5])) as ArrayRef,
                Arc::new(StringArray::from(vec!["A/G", "A/T", "C/G", "C/T", "C/A"])) as ArrayRef,
            ],
        )
        .unwrap();
        let store = tmp.path().join("store");
        fs::create_dir_all(store.join("batches")).unwrap();
        let index_path = tmp.path().join("chr1.range_index.pzstd");
        let mut builder = ArrowStoreBuilder::new(
            store.clone(),
            index_path.clone(),
            vec![
                "position_key".to_string(),
                "end".to_string(),
                "allele_string".to_string(),
            ],
            ArrowStoreCompression::Uncompressed,
            3,
            2,
            Some(
                PackedRangeIndexBuilder::create(&index_path, 2, DEFAULT_INDEX_ZSTD_LEVEL).unwrap(),
            ),
        );
        builder.append_batch(batch).unwrap();
        let manifest = builder.finish("chr1").unwrap();
        assert_eq!(manifest.batch_count, 2);
        let mut index = PackedRangeIndexReader::open(&index_path, 3).unwrap();
        let range10 = index.lookup_position(10).unwrap().unwrap();
        let range20 = index.lookup_position(20).unwrap().unwrap();
        assert_eq!(range10.batch_id, 0);
        assert_eq!(range10.start_row, 0);
        assert_eq!(range10.len, 2);
        assert_eq!(range20.batch_id, 1);
        assert_eq!(range20.start_row, 0);
        assert_eq!(range20.len, 3);
    }

    #[test]
    fn builder_rejects_out_of_order_position_key() {
        let tmp = TempDir::new().unwrap();
        let schema = Arc::new(Schema::new(vec![Field::new(
            "position_key",
            DataType::UInt64,
            false,
        )]));
        let batch = RecordBatch::try_new(
            schema,
            vec![Arc::new(UInt64Array::from(vec![20, 10])) as ArrayRef],
        )
        .unwrap();
        let store = tmp.path().join("store");
        fs::create_dir_all(store.join("batches")).unwrap();
        let index_path = tmp.path().join("chr1.range_index.pzstd");
        let mut builder = ArrowStoreBuilder::new(
            store,
            index_path.clone(),
            vec!["position_key".to_string()],
            ArrowStoreCompression::Uncompressed,
            3,
            2,
            Some(
                PackedRangeIndexBuilder::create(&index_path, 2, DEFAULT_INDEX_ZSTD_LEVEL).unwrap(),
            ),
        );
        let err = builder.append_batch(batch).unwrap_err();
        assert!(err.to_string().contains("rows out of position_key order"));
    }

    #[test]
    fn packed_block_round_trips_ranges() {
        let entries = vec![
            PackedRangeEntry {
                position_key: 10,
                batch_id: 0,
                start_row: 0,
                len: 2,
            },
            PackedRangeEntry {
                position_key: 20,
                batch_id: 1,
                start_row: 0,
                len: 3,
            },
        ];
        let raw = encode_packed_block(&entries).unwrap();
        let decoded = decode_packed_block(&raw, entries.len()).unwrap();
        assert_eq!(decoded.base_key, 10);
        assert_eq!(decoded.key_deltas, vec![0, 10]);
        assert_eq!(decoded.batch_ids, vec![0, 1]);
        assert_eq!(decoded.start_rows, vec![0, 0]);
        assert_eq!(decoded.lens, vec![2, 3]);
    }

    #[test]
    fn arrow_variants_return_same_checksum() {
        let tmp = TempDir::new().unwrap();
        let batch = test_batch();
        let mut checksums = Vec::new();
        let mut rows = Vec::new();
        for compression in ArrowStoreCompression::all() {
            let store_path = build_test_store(tmp.path(), compression, batch.clone(), 3);
            let mut store = ArrowMmapStore::open(&store_path).unwrap();
            let result = store
                .run_lookup(compression, &[10, 20], 3, Duration::ZERO)
                .unwrap();
            checksums.push(result.checksum);
            rows.push(result.selected_rows);
        }
        assert!(checksums.iter().all(|value| *value == checksums[0]));
        assert_eq!(rows, vec![5, 5, 5]);
    }

    #[test]
    fn missing_position_is_counted_as_miss() {
        let tmp = TempDir::new().unwrap();
        let store_path = build_test_store(
            tmp.path(),
            ArrowStoreCompression::Uncompressed,
            test_batch(),
            3,
        );
        let mut store = ArrowMmapStore::open(&store_path).unwrap();
        let result = store
            .run_lookup(
                ArrowStoreCompression::Uncompressed,
                &[10, 999],
                3,
                Duration::ZERO,
            )
            .unwrap();
        assert_eq!(result.positions_found, 1);
        assert_eq!(result.missing_positions, 1);
        assert_eq!(result.index_gets, 2);
        assert_eq!(result.selected_rows, 2);
    }

    fn test_batch() -> RecordBatch {
        let schema = Arc::new(Schema::new(vec![
            Field::new("position_key", DataType::UInt64, false),
            Field::new("end", DataType::Int64, false),
            Field::new("allele_string", DataType::Utf8, false),
        ]));
        RecordBatch::try_new(
            schema,
            vec![
                Arc::new(UInt64Array::from(vec![10, 10, 20, 20, 20])) as ArrayRef,
                Arc::new(Int64Array::from(vec![1, 2, 3, 4, 5])) as ArrayRef,
                Arc::new(StringArray::from(vec!["A/G", "A/T", "C/G", "C/T", "C/A"])) as ArrayRef,
            ],
        )
        .unwrap()
    }

    fn build_test_store(
        root: &Path,
        compression: ArrowStoreCompression,
        batch: RecordBatch,
        target_rows: usize,
    ) -> PathBuf {
        let store = root.join(compression.dirname("chr1"));
        fs::create_dir_all(store.join("batches")).unwrap();
        let index_path = root.join(format!("{}.range_index.pzstd", compression.as_str()));
        let mut builder = ArrowStoreBuilder::new(
            store.clone(),
            index_path.clone(),
            vec![
                "position_key".to_string(),
                "end".to_string(),
                "allele_string".to_string(),
            ],
            compression,
            target_rows,
            2,
            Some(
                PackedRangeIndexBuilder::create(&index_path, 2, DEFAULT_INDEX_ZSTD_LEVEL).unwrap(),
            ),
        );
        builder.append_batch(batch).unwrap();
        let mut manifest = builder.finish("chr1").unwrap();
        manifest.index_size = dir_size(&index_path).unwrap();
        manifest.total_size = dir_size(&store).unwrap() + manifest.index_size;
        fs::write(
            store.join("manifest.json"),
            serde_json::to_string_pretty(&manifest).unwrap(),
        )
        .unwrap();
        store
    }
}
