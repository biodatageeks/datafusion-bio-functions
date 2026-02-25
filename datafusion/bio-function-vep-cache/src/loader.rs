//! Parquet -> fjall cache ingestion pipeline (streaming).
//!
//! Loads a VEP variation cache from Parquet into a fjall KV store,
//! grouping variants into windows and parallelizing across DataFusion
//! physical partitions.
//!
//! Streaming approach:
//! - Each partition processes its stream batch-by-batch, splitting rows into
//!   per-window buffers local to that partition.
//! - Windows are flushed to fjall as soon as the stream advances past them
//!   (position-order tracking), keeping memory bounded to a few active windows.
//! - Cross-partition window overlap (rare — only at partition boundaries) is
//!   handled via atomic read-modify-write under sharded window locks.

use std::collections::HashMap;
use std::collections::hash_map::DefaultHasher;
use std::hash::{Hash, Hasher};
use std::path::Path;
use std::sync::Arc;
use std::time::Instant;

use datafusion::arrow::array::{
    Array, AsArray, RecordBatch, StringArray, StringViewArray, UInt32Array,
};
use datafusion::arrow::datatypes::{Int64Type, SchemaRef};
use datafusion::common::{DataFusionError, Result};
use datafusion::physical_plan::{ExecutionPlan, ExecutionPlanProperties};
use datafusion::prelude::SessionContext;
use futures::StreamExt;
use log::{debug, info};

use crate::key_encoding::{
    DEFAULT_WINDOW_SIZE, EntryType, encode_entry_key, encode_position_key, window_id_for_position,
};
use crate::kv_store::{
    FORMAT_V1, FORMAT_V2, FORMAT_V3, FORMAT_V4, FORMAT_V5, VepKvStore, serialize_ipc_pub,
    to_v1_column_index, validate_v1_schema_width,
};
use crate::mmap_block_store::WindowBlockCodec;
use crate::position_entry::serialize_position_entry;
use crate::position_index::PositionIndex;

/// Statistics returned after loading.
#[derive(Debug, Clone)]
pub struct LoadStats {
    pub total_variants: u64,
    pub total_windows: u64,
    pub total_bytes: u64,
    pub elapsed_secs: f64,
}

/// Per-partition stats accumulated during streaming.
struct PartitionStats {
    total_variants: u64,
    total_windows: u64,
    total_bytes: u64,
}

const FLUSH_LOCK_SHARDS: usize = 256;

/// Striped lock set for window-level read-modify-write coordination.
struct FlushShards {
    locks: Vec<std::sync::Mutex<()>>,
}

impl FlushShards {
    fn new(shards: usize) -> Self {
        let shard_count = shards.max(1);
        let mut locks = Vec::with_capacity(shard_count);
        for _ in 0..shard_count {
            locks.push(std::sync::Mutex::new(()));
        }
        Self { locks }
    }

    fn shard_index(&self, chrom: &str, window_id: u64) -> usize {
        let mut hasher = DefaultHasher::new();
        chrom.hash(&mut hasher);
        window_id.hash(&mut hasher);
        (hasher.finish() as usize) % self.locks.len()
    }
}

/// Per-window buffer drained from a partition's pending map and flushed.
struct PendingWindow {
    chrom: String,
    window_id: u64,
    batches: Vec<RecordBatch>,
}

/// Shared context passed to each partition's streaming task.
struct StreamContext {
    store: Arc<VepKvStore>,
    schema: SchemaRef,
    flush_locks: Arc<FlushShards>,
    window_size: u64,
    v5_zstd_level: i32,
    v5_dict_size_kb: u32,
}

/// Builder for loading a Parquet-based VEP cache into fjall.
pub struct CacheLoader {
    source_table: String,
    target_path: String,
    window_size: u64,
    parallelism: Option<usize>,
    format_version: u8,
    v2_block_codec: WindowBlockCodec,
    v5_zstd_level: i32,
    v5_dict_size_kb: u32,
}

impl CacheLoader {
    /// Create a new loader.
    ///
    /// `source_table` is the name of a table already registered in the session.
    /// `target_path` is the filesystem path for the fjall database.
    pub fn new(source_table: impl Into<String>, target_path: impl Into<String>) -> Self {
        Self {
            source_table: source_table.into(),
            target_path: target_path.into(),
            window_size: DEFAULT_WINDOW_SIZE,
            parallelism: None,
            format_version: FORMAT_V1,
            v2_block_codec: WindowBlockCodec::None,
            v5_zstd_level: 3,
            v5_dict_size_kb: 112,
        }
    }

    /// Set the window size (default: 10 Kb).
    pub fn with_window_size(mut self, size: u64) -> Self {
        self.window_size = size;
        self
    }

    /// Set an optional concurrency cap on partition-level parallelism.
    /// Default is `None` — all partitions run concurrently.
    pub fn with_parallelism(mut self, n: usize) -> Self {
        self.parallelism = Some(n);
        self
    }

    /// Select cache format version.
    ///
    /// Supported:
    /// - `FORMAT_V1`: position index + per-column Arrow IPC in Fjall
    /// - `FORMAT_V2`: position index in Fjall + mmap per-column window payload blocks
    /// - `FORMAT_V3`: position index in Fjall + mmap whole-batch Arrow IPC blocks (fastest)
    pub fn with_format_version(mut self, format_version: u8) -> Self {
        self.format_version = format_version;
        self
    }

    /// Set v2 window payload IPC compression codec.
    ///
    /// Ignored when format_version is FORMAT_V1.
    pub fn with_v2_block_codec(mut self, codec: WindowBlockCodec) -> Self {
        self.v2_block_codec = codec;
        self
    }

    /// Set the zstd compression level for V5 caches (default: 3).
    ///
    /// Higher levels produce smaller caches at the cost of slower writes.
    /// Decompression speed is constant regardless of level.
    pub fn with_v5_zstd_level(mut self, level: i32) -> Self {
        self.v5_zstd_level = level;
        self
    }

    /// Set the zstd dictionary size in KB for V5 caches (default: 112).
    pub fn with_v5_dict_size_kb(mut self, size_kb: u32) -> Self {
        self.v5_dict_size_kb = size_kb;
        self
    }

    /// Load the cache from the registered source table into fjall.
    ///
    /// Spawns one streaming task per physical partition. Each task reads
    /// batches, splits them into windows, and flushes completed windows
    /// to fjall incrementally — memory stays bounded to a few active
    /// windows per partition rather than the entire dataset.
    pub async fn load(&self, session: &SessionContext) -> Result<LoadStats> {
        let start_time = Instant::now();

        let source_df = session.table(&self.source_table).await?;
        let schema: SchemaRef = Arc::new(source_df.schema().as_arrow().clone());

        if self.format_version == FORMAT_V1 {
            validate_v1_schema_width(schema.fields().len())?;
        }

        let store = Arc::new(match self.format_version {
            FORMAT_V1 => VepKvStore::create(
                Path::new(&self.target_path),
                schema.clone(),
                self.window_size,
            )?,
            FORMAT_V2 => VepKvStore::create_v2_with_codec(
                Path::new(&self.target_path),
                schema.clone(),
                self.window_size,
                self.v2_block_codec,
            )?,
            FORMAT_V3 => VepKvStore::create_v3_with_codec(
                Path::new(&self.target_path),
                schema.clone(),
                self.window_size,
                self.v2_block_codec,
            )?,
            FORMAT_V4 => VepKvStore::create_v4_with_codec(
                Path::new(&self.target_path),
                schema.clone(),
                self.window_size,
                self.v2_block_codec,
            )?,
            FORMAT_V5 => VepKvStore::create_v5(Path::new(&self.target_path), schema.clone())?,
            other => {
                return Err(DataFusionError::Execution(format!(
                    "unsupported cache format version in loader: {other} (supported: {FORMAT_V1}, {FORMAT_V2}, {FORMAT_V3}, {FORMAT_V4}, {FORMAT_V5})"
                )));
            }
        });

        let plan = source_df.create_physical_plan().await?;
        let partition_count = plan.output_partitioning().partition_count();
        let task_ctx = session.task_ctx();

        info!(
            "Loading into KV cache: {} partitions, window_size={}, format_version={}, v2_block_codec={}",
            partition_count,
            self.window_size,
            self.format_version,
            self.v2_block_codec.as_str()
        );
        debug!(
            "Physical plan:\n{}",
            datafusion::physical_plan::displayable(plan.as_ref()).indent(true)
        );

        // Sharded locks for cross-partition read-modify-write safety.
        let flush_locks = Arc::new(FlushShards::new(FLUSH_LOCK_SHARDS));
        let semaphore = self
            .parallelism
            .map(|n| Arc::new(tokio::sync::Semaphore::new(n.min(partition_count))));

        let mut handles = Vec::with_capacity(partition_count);
        let window_size = self.window_size;

        for partition_id in 0..partition_count {
            let plan = Arc::clone(&plan);
            let task_ctx = Arc::clone(&task_ctx);
            let sem = semaphore.clone();
            let ctx = StreamContext {
                store: Arc::clone(&store),
                schema: schema.clone(),
                flush_locks: Arc::clone(&flush_locks),
                window_size,
                v5_zstd_level: self.v5_zstd_level,
                v5_dict_size_kb: self.v5_dict_size_kb,
            };

            handles.push(tokio::spawn(async move {
                if let Some(ref sem) = sem {
                    let _permit = sem
                        .acquire()
                        .await
                        .map_err(|e| DataFusionError::External(Box::new(e)))?;
                }

                stream_partition(plan, task_ctx, partition_id, ctx).await
            }));
        }

        // Collect stats from all partitions.
        let mut total_variants = 0u64;
        let mut total_windows = 0u64;
        let mut total_bytes = 0u64;

        for handle in handles {
            let ps = handle
                .await
                .map_err(|e| DataFusionError::External(Box::new(e)))??;
            total_variants += ps.total_variants;
            total_windows += ps.total_windows;
            total_bytes += ps.total_bytes;
        }

        store.persist()?;

        let elapsed = start_time.elapsed().as_secs_f64();
        info!(
            "Loaded {total_variants} variants into {total_windows} windows ({total_bytes} bytes) in {elapsed:.1}s"
        );

        Ok(LoadStats {
            total_variants,
            total_windows,
            total_bytes,
            elapsed_secs: elapsed,
        })
    }
}

/// Stream a single partition: read batches, split into windows, flush to fjall incrementally.
///
/// Uses position-order tracking to flush completed windows eagerly: once the
/// stream advances past a window (sorted input), that window is sealed and
/// written. At end-of-stream, all remaining windows are flushed.
///
/// Cross-partition overlap (rare — only at partition boundaries) is handled
/// by read-modify-write under sharded window locks.
async fn stream_partition(
    plan: Arc<dyn ExecutionPlan>,
    task_ctx: Arc<datafusion::execution::TaskContext>,
    partition_id: usize,
    ctx: StreamContext,
) -> Result<PartitionStats> {
    let part_start = Instant::now();
    let mut stream = plan.execute(partition_id, task_ctx)?;

    let chrom_col_idx = ctx.schema.index_of("chrom")?;
    let start_col_idx = ctx.schema.index_of("start")?;

    // Local per-window buffers (NOT shared across partitions).
    let mut pending: HashMap<(String, u64), Vec<RecordBatch>> = HashMap::new();

    // Track the maximum primary window ID seen per chromosome (frontier).
    // Windows behind the frontier are sealed — no more data will arrive.
    let mut frontier: HashMap<String, u64> = HashMap::new();

    let mut stats = PartitionStats {
        total_variants: 0,
        total_windows: 0,
        total_bytes: 0,
    };

    let is_v5 = ctx.store.format_version() == FORMAT_V5;

    // V5 zstd dictionary compressor — built from first batch, reused for all subsequent.
    // None until the first batch trains the dictionary. If training is skipped (too few
    // samples), this stays None and all entries are stored uncompressed.
    let mut v5_compressor: Option<zstd::bulk::Compressor<'static>> = None;
    let mut v5_first_batch_done = false;

    while let Some(batch_result) = stream.next().await {
        let batch = batch_result?;
        if batch.num_rows() == 0 {
            continue;
        }

        stats.total_variants += batch.num_rows() as u64;

        if is_v5 {
            let store = Arc::clone(&ctx.store);
            let schema = ctx.schema.clone();
            let batch_clone = batch.clone();

            if !v5_first_batch_done {
                v5_first_batch_done = true;
                let zstd_level = ctx.v5_zstd_level;
                let dict_size_kb = ctx.v5_dict_size_kb;
                // First batch: try to train dictionary, store it, build compressor.
                // If training fails (too few samples), entries are inserted uncompressed.
                let (compressor, positions, bytes) = tokio::task::spawn_blocking(move || {
                    train_and_flush_first_v5_batch(
                        &store,
                        &schema,
                        &batch_clone,
                        zstd_level,
                        dict_size_kb,
                    )
                })
                .await
                .map_err(|e| DataFusionError::External(Box::new(e)))??;
                v5_compressor = compressor;
                stats.total_windows += positions;
                stats.total_bytes += bytes;
            } else if let Some(mut compressor) = v5_compressor.take() {
                // Subsequent batches with dictionary: compress and insert.
                let (comp_back, positions, bytes) = tokio::task::spawn_blocking(move || {
                    let result = flush_positions_v5_compressed(
                        &store,
                        &schema,
                        &batch_clone,
                        &mut compressor,
                    );
                    result.map(|(p, b)| (compressor, p, b))
                })
                .await
                .map_err(|e| DataFusionError::External(Box::new(e)))??;
                v5_compressor = Some(comp_back);
                stats.total_windows += positions;
                stats.total_bytes += bytes;
            } else {
                // Subsequent batches without dictionary: insert uncompressed.
                let (positions, bytes) = tokio::task::spawn_blocking(move || {
                    flush_positions_v5_uncompressed(&store, &schema, &batch_clone)
                })
                .await
                .map_err(|e| DataFusionError::External(Box::new(e)))??;
                stats.total_windows += positions;
                stats.total_bytes += bytes;
            }
            continue;
        }

        // V1-V4: window-based path.
        // Update frontier with max primary window per chrom in this batch.
        let starts = batch.column(start_col_idx).as_primitive::<Int64Type>();
        let chroms = batch.column(chrom_col_idx);
        for row in 0..batch.num_rows() {
            let chrom = string_value(chroms.as_ref(), row);
            let wid = window_id_for_position(starts.value(row), ctx.window_size);
            frontier
                .entry(chrom.to_string())
                .and_modify(|f| *f = (*f).max(wid))
                .or_insert(wid);
        }

        // Split batch into per-window sub-batches and accumulate.
        let window_batches = split_batch_into_windows(&batch, ctx.window_size)?;
        for ((chrom, wid), sub_batch) in window_batches {
            pending.entry((chrom, wid)).or_default().push(sub_batch);
        }

        // Flush windows that are behind the frontier (sealed — no more data coming).
        let flush_keys: Vec<_> = pending
            .keys()
            .filter(|(chrom, wid)| frontier.get(chrom).is_some_and(|f| *wid < *f))
            .cloned()
            .collect();

        if !flush_keys.is_empty() {
            let (new_windows, bytes) = flush_windows(&ctx, &mut pending, &flush_keys).await?;
            stats.total_windows += new_windows;
            stats.total_bytes += bytes;
        }
    }

    // Flush all remaining windows (the frontier windows themselves).
    // V5 flushes inline, so pending is always empty.
    let remaining: Vec<_> = pending.keys().cloned().collect();
    if !remaining.is_empty() {
        let (new_windows, bytes) = flush_windows(&ctx, &mut pending, &remaining).await?;
        stats.total_windows += new_windows;
        stats.total_bytes += bytes;
    }

    info!(
        "Partition {partition_id}: {} variants, {} windows in {:.1}s",
        stats.total_variants,
        stats.total_windows,
        part_start.elapsed().as_secs_f64()
    );

    Ok(stats)
}

/// Read an existing window's full data from fjall, if present.
///
/// Used for read-modify-write when cross-partition overlap occurs.
fn read_existing_window(
    store: &VepKvStore,
    schema: &SchemaRef,
    chrom: &str,
    window_id: u64,
) -> Result<Option<RecordBatch>> {
    match store.format_version() {
        FORMAT_V1 => {
            if store.get_position_index(chrom, window_id)?.is_none() {
                return Ok(None);
            }
            let col_indices: Vec<u8> = (0..schema.fields().len())
                .map(to_v1_column_index)
                .collect::<Result<Vec<_>>>()?;
            match store.get_columns(chrom, window_id, &col_indices)? {
                Some(col_data) => {
                    let arrays: Vec<_> = col_data.into_iter().map(|(arr, _)| arr).collect();
                    let batch = RecordBatch::try_new(schema.clone(), arrays)
                        .map_err(|e| DataFusionError::ArrowError(Box::new(e), None))?;
                    Ok(Some(batch))
                }
                None => Ok(None),
            }
        }
        FORMAT_V2 => store.get_window_block_v2(chrom, window_id),
        FORMAT_V3 => store.get_window_batch_v3(chrom, window_id),
        FORMAT_V4 => store.get_window_batch_v4(chrom, window_id),
        other => Err(DataFusionError::Execution(format!(
            "unsupported cache format version in read_existing_window: {other}"
        ))),
    }
}

/// Train a zstd dictionary from the first V5 batch, store it, then compress + insert all entries.
///
/// Returns `(Option<Compressor>, num_positions, total_bytes)`. The compressor is `None` if
/// there weren't enough samples to train a meaningful dictionary (entries stored uncompressed).
fn train_and_flush_first_v5_batch(
    store: &VepKvStore,
    schema: &SchemaRef,
    batch: &RecordBatch,
    zstd_level: i32,
    dict_size_kb: u32,
) -> Result<(Option<zstd::bulk::Compressor<'static>>, u64, u64)> {
    let chrom_col_idx = schema.index_of("chrom")?;
    let start_col_idx = schema.index_of("start")?;
    let end_col_idx = schema.index_of("end")?;
    let allele_col_idx = schema.index_of("allele_string")?;

    let chrom_col = batch.column(chrom_col_idx);
    let starts = batch.column(start_col_idx).as_primitive::<Int64Type>();
    let ends = batch.column(end_col_idx).as_primitive::<Int64Type>();

    let col_indices: Vec<usize> = (0..schema.fields().len())
        .filter(|&i| i != chrom_col_idx && i != start_col_idx && i != end_col_idx)
        .collect();

    // Group row indices by (chrom, start, end).
    let mut groups: HashMap<(String, i64, i64), Vec<usize>> = HashMap::new();
    for row in 0..batch.num_rows() {
        let chrom = string_value(chrom_col.as_ref(), row);
        let start = starts.value(row);
        let end = ends.value(row);
        groups
            .entry((chrom.to_string(), start, end))
            .or_default()
            .push(row);
    }

    // Phase 1: serialize all entries (uncompressed) to collect training samples.
    let mut serialized: Vec<(Vec<u8>, Vec<u8>)> = Vec::with_capacity(groups.len());
    for ((chrom, start, end), rows) in &groups {
        let value = serialize_position_entry(rows, batch, &col_indices, allele_col_idx)?;
        let key = encode_position_key(chrom, *start, *end);
        serialized.push((key, value));
    }

    // Phase 2: train dictionary from samples (up to 10,000).
    // Require a minimum of 100 samples for meaningful dictionary training.
    let min_training_samples = 100;
    let dict_bytes = if serialized.len() >= min_training_samples {
        let max_samples = 10_000.min(serialized.len());
        let samples: Vec<&[u8]> = serialized[..max_samples]
            .iter()
            .map(|(_, v)| v.as_slice())
            .collect();
        match train_v5_dict(&samples, dict_size_kb) {
            Ok(dict) => {
                info!(
                    "V5 zstd dict trained: {} samples, dict {} bytes",
                    max_samples,
                    dict.len()
                );
                Some(dict)
            }
            Err(e) => {
                info!("V5 zstd dict training failed (falling back to uncompressed): {e}");
                None
            }
        }
    } else {
        info!(
            "V5: only {} entries, skipping dictionary training (need >= {})",
            serialized.len(),
            min_training_samples
        );
        None
    };

    if let Some(ref dict) = dict_bytes {
        // Store dictionary and compression level in meta, then compress entries.
        store.store_v5_dict(dict)?;
        store.store_v5_zstd_level(zstd_level)?;

        let mut compressor =
            zstd::bulk::Compressor::with_dictionary(zstd_level, dict).map_err(|e| {
                DataFusionError::Execution(format!("failed to create zstd compressor: {e}"))
            })?;

        let mut entries: Vec<(Vec<u8>, Vec<u8>)> = Vec::with_capacity(serialized.len());
        let mut total_bytes = 0u64;
        for (key, raw_value) in serialized {
            let compressed = compressor
                .compress(&raw_value)
                .map_err(|e| DataFusionError::Execution(format!("zstd compression failed: {e}")))?;
            total_bytes += compressed.len() as u64;
            entries.push((key, compressed));
        }

        let num_positions = entries.len() as u64;
        store.batch_insert_position_entries(&entries)?;

        Ok((Some(compressor), num_positions, total_bytes))
    } else {
        // Not enough samples — insert uncompressed, no compressor for subsequent batches.
        let mut total_bytes = 0u64;
        for (_, v) in &serialized {
            total_bytes += v.len() as u64;
        }
        let num_positions = serialized.len() as u64;
        store.batch_insert_position_entries(&serialized)?;

        Ok((None, num_positions, total_bytes))
    }
}

/// Flush V5 position entries without compression (fallback when dictionary training was skipped).
fn flush_positions_v5_uncompressed(
    store: &VepKvStore,
    schema: &SchemaRef,
    batch: &RecordBatch,
) -> Result<(u64, u64)> {
    let chrom_col_idx = schema.index_of("chrom")?;
    let start_col_idx = schema.index_of("start")?;
    let end_col_idx = schema.index_of("end")?;
    let allele_col_idx = schema.index_of("allele_string")?;

    let chrom_col = batch.column(chrom_col_idx);
    let starts = batch.column(start_col_idx).as_primitive::<Int64Type>();
    let ends = batch.column(end_col_idx).as_primitive::<Int64Type>();

    let col_indices: Vec<usize> = (0..schema.fields().len())
        .filter(|&i| i != chrom_col_idx && i != start_col_idx && i != end_col_idx)
        .collect();

    let mut groups: HashMap<(String, i64, i64), Vec<usize>> = HashMap::new();
    for row in 0..batch.num_rows() {
        let chrom = string_value(chrom_col.as_ref(), row);
        let start = starts.value(row);
        let end = ends.value(row);
        groups
            .entry((chrom.to_string(), start, end))
            .or_default()
            .push(row);
    }

    let mut entries: Vec<(Vec<u8>, Vec<u8>)> = Vec::with_capacity(groups.len());
    let mut total_bytes = 0u64;

    for ((chrom, start, end), rows) in &groups {
        let value = serialize_position_entry(rows, batch, &col_indices, allele_col_idx)?;
        let key = encode_position_key(chrom, *start, *end);
        total_bytes += value.len() as u64;
        entries.push((key, value));
    }

    let num_positions = entries.len() as u64;
    store.batch_insert_position_entries(&entries)?;

    Ok((num_positions, total_bytes))
}

/// Compress and flush V5 position entries using a pre-trained zstd dictionary compressor.
fn flush_positions_v5_compressed(
    store: &VepKvStore,
    schema: &SchemaRef,
    batch: &RecordBatch,
    compressor: &mut zstd::bulk::Compressor<'_>,
) -> Result<(u64, u64)> {
    let chrom_col_idx = schema.index_of("chrom")?;
    let start_col_idx = schema.index_of("start")?;
    let end_col_idx = schema.index_of("end")?;
    let allele_col_idx = schema.index_of("allele_string")?;

    let chrom_col = batch.column(chrom_col_idx);
    let starts = batch.column(start_col_idx).as_primitive::<Int64Type>();
    let ends = batch.column(end_col_idx).as_primitive::<Int64Type>();

    let col_indices: Vec<usize> = (0..schema.fields().len())
        .filter(|&i| i != chrom_col_idx && i != start_col_idx && i != end_col_idx)
        .collect();

    let mut groups: HashMap<(String, i64, i64), Vec<usize>> = HashMap::new();
    for row in 0..batch.num_rows() {
        let chrom = string_value(chrom_col.as_ref(), row);
        let start = starts.value(row);
        let end = ends.value(row);
        groups
            .entry((chrom.to_string(), start, end))
            .or_default()
            .push(row);
    }

    let mut entries: Vec<(Vec<u8>, Vec<u8>)> = Vec::with_capacity(groups.len());
    let mut total_bytes = 0u64;

    for ((chrom, start, end), rows) in &groups {
        let raw_value = serialize_position_entry(rows, batch, &col_indices, allele_col_idx)?;
        let compressed = compressor
            .compress(&raw_value)
            .map_err(|e| DataFusionError::Execution(format!("zstd compression failed: {e}")))?;
        let key = encode_position_key(chrom, *start, *end);
        total_bytes += compressed.len() as u64;
        entries.push((key, compressed));
    }

    let num_positions = entries.len() as u64;
    store.batch_insert_position_entries(&entries)?;

    Ok((num_positions, total_bytes))
}

/// Train a zstd dictionary from serialized position entry samples.
fn train_v5_dict(samples: &[&[u8]], dict_size_kb: u32) -> Result<Vec<u8>> {
    // Flatten samples into continuous buffer with sizes array for zstd training API.
    let mut continuous = Vec::new();
    let mut sizes = Vec::with_capacity(samples.len());
    for sample in samples {
        continuous.extend_from_slice(sample);
        sizes.push(sample.len());
    }

    let dict_size = dict_size_kb as usize * 1024;
    let dict = zstd::dict::from_continuous(&continuous, &sizes, dict_size)
        .map_err(|e| DataFusionError::Execution(format!("zstd dictionary training failed: {e}")))?;
    Ok(dict)
}

/// Flush pending window data to fjall.
///
/// Drains the requested keys from `pending`, then offloads read/merge/serialize/write
/// work to `spawn_blocking` so async stream tasks stay responsive.
async fn flush_windows(
    ctx: &StreamContext,
    pending: &mut HashMap<(String, u64), Vec<RecordBatch>>,
    keys: &[(String, u64)],
) -> Result<(u64, u64)> {
    let mut windows = Vec::with_capacity(keys.len());

    for (chrom, wid) in keys {
        if let Some(batches) = pending.remove(&(chrom.clone(), *wid)) {
            windows.push(PendingWindow {
                chrom: chrom.clone(),
                window_id: *wid,
                batches,
            });
        }
    }

    if windows.is_empty() {
        return Ok((0, 0));
    }

    let store = Arc::clone(&ctx.store);
    let schema = ctx.schema.clone();
    let flush_locks = Arc::clone(&ctx.flush_locks);

    tokio::task::spawn_blocking(move || flush_windows_blocking(store, schema, windows, flush_locks))
        .await
        .map_err(|e| DataFusionError::External(Box::new(e)))?
}

/// Blocking read-modify-write flush path.
///
/// Windows are grouped by lock shard to avoid a single global serialization point
/// while preserving atomicity for each window's merge+write sequence.
fn flush_windows_blocking(
    store: Arc<VepKvStore>,
    schema: SchemaRef,
    windows: Vec<PendingWindow>,
    flush_locks: Arc<FlushShards>,
) -> Result<(u64, u64)> {
    let mut windows_by_shard: HashMap<usize, Vec<PendingWindow>> = HashMap::new();
    for window in windows {
        let shard_idx = flush_locks.shard_index(&window.chrom, window.window_id);
        windows_by_shard.entry(shard_idx).or_default().push(window);
    }

    let mut total_windows = 0u64;
    let mut total_bytes = 0u64;

    for (shard_idx, shard_windows) in windows_by_shard {
        let _guard = flush_locks.locks[shard_idx].lock().unwrap();
        match store.format_version() {
            FORMAT_V1 => {
                let mut all_entries: Vec<(Vec<u8>, Vec<u8>)> = Vec::new();
                for mut window in shard_windows {
                    // Read-modify-write: prepend existing data if another partition
                    // already wrote this window.
                    let existing =
                        read_existing_window(&store, &schema, &window.chrom, window.window_id)?;
                    let is_new = existing.is_none();
                    if let Some(existing_batch) = existing {
                        window.batches.insert(0, existing_batch);
                    }

                    let prepared = prepare_window_entries(
                        &schema,
                        &window.chrom,
                        window.window_id,
                        window.batches,
                    )?;
                    total_bytes += prepared.total_bytes;
                    if is_new {
                        total_windows += 1;
                    }
                    all_entries.extend(prepared.entries);
                }

                if !all_entries.is_empty() {
                    store.batch_insert_raw(&all_entries)?;
                }
            }
            FORMAT_V2 => {
                for mut window in shard_windows {
                    let existing =
                        read_existing_window(&store, &schema, &window.chrom, window.window_id)?;
                    let is_new = existing.is_none();
                    if let Some(existing_batch) = existing {
                        window.batches.insert(0, existing_batch);
                    }

                    let Some(merged) = merge_window_batches(&schema, window.batches)? else {
                        continue;
                    };
                    total_bytes += merged.get_array_memory_size() as u64;
                    if is_new {
                        total_windows += 1;
                    }

                    let pos_index = PositionIndex::from_batch(&merged)?;
                    store.put_position_index(&window.chrom, window.window_id, &pos_index)?;
                    store.put_window_block_v2(&window.chrom, window.window_id, &merged)?;
                }
            }
            FORMAT_V3 => {
                for mut window in shard_windows {
                    let existing =
                        read_existing_window(&store, &schema, &window.chrom, window.window_id)?;
                    let is_new = existing.is_none();
                    if let Some(existing_batch) = existing {
                        window.batches.insert(0, existing_batch);
                    }

                    let Some(merged) = merge_window_batches(&schema, window.batches)? else {
                        continue;
                    };
                    total_bytes += merged.get_array_memory_size() as u64;
                    if is_new {
                        total_windows += 1;
                    }

                    let pos_index = PositionIndex::from_batch(&merged)?;
                    store.put_position_index(&window.chrom, window.window_id, &pos_index)?;
                    store.put_window_block_v3(&window.chrom, window.window_id, &merged)?;
                }
            }
            FORMAT_V4 => {
                for mut window in shard_windows {
                    let existing =
                        read_existing_window(&store, &schema, &window.chrom, window.window_id)?;
                    let is_new = existing.is_none();
                    if let Some(existing_batch) = existing {
                        window.batches.insert(0, existing_batch);
                    }

                    let Some(merged) = merge_window_batches(&schema, window.batches)? else {
                        continue;
                    };
                    total_bytes += merged.get_array_memory_size() as u64;
                    if is_new {
                        total_windows += 1;
                    }

                    let pos_index = PositionIndex::from_batch(&merged)?;
                    store.put_position_index(&window.chrom, window.window_id, &pos_index)?;
                    store.put_window_block_v4(&window.chrom, window.window_id, &merged)?;
                }
            }
            other => {
                return Err(DataFusionError::Execution(format!(
                    "unsupported cache format version in flush_windows_blocking: {other}"
                )));
            }
        }
    }

    Ok((total_windows, total_bytes))
}

/// Split a RecordBatch into sub-batches keyed by (chrom, window_id).
///
/// Handles cross-window indels by duplicating rows into extra windows.
fn split_batch_into_windows(
    batch: &RecordBatch,
    window_size: u64,
) -> Result<HashMap<(String, u64), RecordBatch>> {
    if batch.num_rows() == 0 {
        return Ok(HashMap::new());
    }

    let schema = batch.schema();
    let chrom_col_idx = schema.index_of("chrom")?;
    let start_col_idx = schema.index_of("start")?;
    let end_col_idx = schema.index_of("end")?;

    let chrom_col = batch.column(chrom_col_idx);
    let starts = batch.column(start_col_idx).as_primitive::<Int64Type>();
    let ends = batch.column(end_col_idx).as_primitive::<Int64Type>();

    // Group row indices by (chrom, window_id).
    let mut groups: HashMap<(String, u64), Vec<u32>> = HashMap::new();

    for row in 0..batch.num_rows() {
        let chrom = string_value(chrom_col.as_ref(), row);
        let start = starts.value(row);
        let end = ends.value(row);

        let wid = window_id_for_position(start, window_size);

        groups
            .entry((chrom.to_string(), wid))
            .or_default()
            .push(row as u32);

        // Cross-window indels: duplicate into extra windows.
        let end_wid = window_id_for_position(end, window_size);
        for extra_wid in (wid + 1)..=end_wid {
            groups
                .entry((chrom.to_string(), extra_wid))
                .or_default()
                .push(row as u32);
        }
    }

    // Build sub-batches using arrow::compute::take.
    let mut result = HashMap::with_capacity(groups.len());
    for (key, indices) in groups {
        let indices_arr = UInt32Array::from(indices);
        let columns: Vec<_> = batch
            .columns()
            .iter()
            .map(|col| datafusion::arrow::compute::take(col.as_ref(), &indices_arr, None))
            .collect::<std::result::Result<Vec<_>, _>>()
            .map_err(|e| DataFusionError::ArrowError(Box::new(e), None))?;
        let sub_batch = RecordBatch::try_new(schema.clone(), columns)
            .map_err(|e| DataFusionError::ArrowError(Box::new(e), None))?;
        result.insert(key, sub_batch);
    }

    Ok(result)
}

/// Extract a string value from either Utf8 or Utf8View array.
fn string_value(col: &dyn Array, row: usize) -> &str {
    if let Some(arr) = col.as_any().downcast_ref::<StringArray>() {
        arr.value(row)
    } else if let Some(arr) = col.as_any().downcast_ref::<StringViewArray>() {
        arr.value(row)
    } else {
        panic!("Expected Utf8 or Utf8View column for chrom")
    }
}

/// Pre-serialized window data ready for batch insertion into fjall.
struct PreparedWindow {
    /// (key, value) pairs to insert into the data partition.
    entries: Vec<(Vec<u8>, Vec<u8>)>,
    /// Approximate in-memory size of the source data.
    total_bytes: u64,
}

/// Serialize a window's data into (key, value) pairs for batch insertion.
///
/// This is CPU-bound (Arrow IPC + LZ4 compression) and designed to run on
/// `tokio::task::spawn_blocking`.
fn prepare_window_entries(
    schema: &SchemaRef,
    chrom: &str,
    window_id: u64,
    batches: Vec<RecordBatch>,
) -> Result<PreparedWindow> {
    let Some(merged) = merge_window_batches(schema, batches)? else {
        return Ok(PreparedWindow {
            entries: vec![],
            total_bytes: 0,
        });
    };

    let total_bytes = merged.get_array_memory_size() as u64;

    prepare_v1_entries(chrom, window_id, &merged, total_bytes)
}

/// Merge all sub-batches for a single `(chrom, window)` buffer.
fn merge_window_batches(
    schema: &SchemaRef,
    batches: Vec<RecordBatch>,
) -> Result<Option<RecordBatch>> {
    if batches.is_empty() {
        return Ok(None);
    }

    let merged = if batches.len() == 1 {
        batches.into_iter().next().unwrap()
    } else {
        datafusion::arrow::compute::concat_batches(schema, &batches)
            .map_err(|e| DataFusionError::ArrowError(Box::new(e), None))?
    };
    Ok(Some(merged))
}

/// v1 (columnar): serialize position index + per-column Arrow IPC entries.
fn prepare_v1_entries(
    chrom: &str,
    window_id: u64,
    merged: &RecordBatch,
    total_bytes: u64,
) -> Result<PreparedWindow> {
    let num_cols = merged.schema().fields().len();
    let mut entries = Vec::with_capacity(num_cols + 1);

    // Position index entry.
    let pos_index = PositionIndex::from_batch(merged)?;
    let pos_key = encode_entry_key(chrom, window_id, EntryType::PositionIndex);
    entries.push((pos_key, pos_index.to_bytes()));

    // Per-column entries.
    for (col_idx, field) in merged.schema().fields().iter().enumerate() {
        let col_key = encode_entry_key(
            chrom,
            window_id,
            EntryType::Column(to_v1_column_index(col_idx)?),
        );
        let col_batch = RecordBatch::try_new(
            Arc::new(datafusion::arrow::datatypes::Schema::new(vec![
                field.as_ref().clone(),
            ])),
            vec![merged.column(col_idx).clone()],
        )
        .map_err(|e| DataFusionError::ArrowError(Box::new(e), None))?;
        let col_value = serialize_ipc_pub(&col_batch)?;
        entries.push((col_key, col_value));
    }

    Ok(PreparedWindow {
        entries,
        total_bytes,
    })
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::key_encoding::EntryType;
    use datafusion::arrow::array::{ArrayRef, Int64Array, StringArray};
    use datafusion::arrow::datatypes::{DataType, Field, Schema};
    use datafusion::datasource::MemTable;
    use std::sync::Arc;

    fn make_test_table() -> (SchemaRef, MemTable) {
        let schema = Arc::new(Schema::new(vec![
            Field::new("chrom", DataType::Utf8, false),
            Field::new("start", DataType::Int64, false),
            Field::new("end", DataType::Int64, false),
            Field::new("variation_name", DataType::Utf8, true),
            Field::new("allele_string", DataType::Utf8, false),
        ]));
        let batch = RecordBatch::try_new(
            schema.clone(),
            vec![
                Arc::new(StringArray::from(vec!["1", "1", "1", "2"])),
                Arc::new(Int64Array::from(vec![100, 500_000, 1_500_000, 100])),
                Arc::new(Int64Array::from(vec![100, 500_000, 1_500_000, 100])),
                Arc::new(StringArray::from(vec!["rs1", "rs2", "rs3", "rs4"])),
                Arc::new(StringArray::from(vec!["A/G", "C/T", "G/A", "T/C"])),
            ],
        )
        .unwrap();
        let table = MemTable::try_new(schema.clone(), vec![vec![batch]]).unwrap();
        (schema, table)
    }

    fn make_wide_test_table(total_cols: usize) -> (SchemaRef, MemTable) {
        assert!(total_cols >= 5);

        let mut fields = vec![
            Field::new("chrom", DataType::Utf8, false),
            Field::new("start", DataType::Int64, false),
            Field::new("end", DataType::Int64, false),
            Field::new("variation_name", DataType::Utf8, true),
            Field::new("allele_string", DataType::Utf8, false),
        ];
        for idx in 5..total_cols {
            fields.push(Field::new(format!("extra_{idx}"), DataType::Utf8, true));
        }

        let schema = Arc::new(Schema::new(fields));
        let mut columns: Vec<ArrayRef> = vec![
            Arc::new(StringArray::from(vec!["1"])),
            Arc::new(Int64Array::from(vec![100])),
            Arc::new(Int64Array::from(vec![100])),
            Arc::new(StringArray::from(vec![Some("rs_test")])),
            Arc::new(StringArray::from(vec!["A/G"])),
        ];
        for idx in 5..total_cols {
            columns.push(Arc::new(StringArray::from(vec![format!("v{idx}")])));
        }

        let batch = RecordBatch::try_new(schema.clone(), columns).unwrap();
        let table = MemTable::try_new(schema.clone(), vec![vec![batch]]).unwrap();
        (schema, table)
    }

    #[tokio::test(flavor = "multi_thread")]
    async fn test_loader_basic_v1() {
        let dir = tempfile::tempdir().unwrap();
        let ctx = SessionContext::new();

        let (schema, table) = make_test_table();
        ctx.register_table("source", Arc::new(table)).unwrap();

        let loader =
            CacheLoader::new("source", dir.path().to_str().unwrap()).with_window_size(1_000_000);
        let stats = loader.load(&ctx).await.unwrap();

        assert_eq!(stats.total_windows, 3);

        let store = VepKvStore::open(dir.path()).unwrap();
        assert_eq!(store.schema().as_ref(), schema.as_ref());

        // v1: position index should be available
        let idx = store.get_position_index("1", 0).unwrap().unwrap();
        assert_eq!(idx.num_rows(), 2);

        let idx = store.get_position_index("1", 1).unwrap().unwrap();
        assert_eq!(idx.num_rows(), 1);

        let idx = store.get_position_index("2", 0).unwrap().unwrap();
        assert_eq!(idx.num_rows(), 1);

        // v1: can read individual columns
        let cols = store.get_columns("1", 0, &[0, 4]).unwrap().unwrap();
        assert_eq!(cols.len(), 2);
        assert_eq!(cols[0].0.len(), 2); // 2 rows in window 0
    }

    #[tokio::test(flavor = "multi_thread")]
    async fn test_loader_multi_partition() {
        let dir = tempfile::tempdir().unwrap();
        let ctx = SessionContext::new();

        let schema = Arc::new(Schema::new(vec![
            Field::new("chrom", DataType::Utf8, false),
            Field::new("start", DataType::Int64, false),
            Field::new("end", DataType::Int64, false),
            Field::new("variation_name", DataType::Utf8, true),
            Field::new("allele_string", DataType::Utf8, false),
        ]));

        // Partition 1: chrom 1 variants
        let batch1 = RecordBatch::try_new(
            schema.clone(),
            vec![
                Arc::new(StringArray::from(vec!["1", "1"])),
                Arc::new(Int64Array::from(vec![100, 500_000])),
                Arc::new(Int64Array::from(vec![100, 500_000])),
                Arc::new(StringArray::from(vec!["rs1", "rs2"])),
                Arc::new(StringArray::from(vec!["A/G", "C/T"])),
            ],
        )
        .unwrap();

        // Partition 2: mix of chrom 1 and chrom 2
        let batch2 = RecordBatch::try_new(
            schema.clone(),
            vec![
                Arc::new(StringArray::from(vec!["1", "2"])),
                Arc::new(Int64Array::from(vec![1_500_000, 100])),
                Arc::new(Int64Array::from(vec![1_500_000, 100])),
                Arc::new(StringArray::from(vec!["rs3", "rs4"])),
                Arc::new(StringArray::from(vec!["G/A", "T/C"])),
            ],
        )
        .unwrap();

        // Create MemTable with 2 partitions
        let table = MemTable::try_new(schema.clone(), vec![vec![batch1], vec![batch2]]).unwrap();
        ctx.register_table("source", Arc::new(table)).unwrap();

        let loader =
            CacheLoader::new("source", dir.path().to_str().unwrap()).with_window_size(1_000_000);
        let stats = loader.load(&ctx).await.unwrap();

        assert_eq!(stats.total_windows, 3); // (1,0), (1,1), (2,0)
        assert_eq!(stats.total_variants, 4);

        let store = VepKvStore::open(dir.path()).unwrap();

        // Window (1,0) should have 2 rows from partition 1
        let idx = store.get_position_index("1", 0).unwrap().unwrap();
        assert_eq!(idx.num_rows(), 2);

        // Window (1,1) should have 1 row from partition 2
        let idx = store.get_position_index("1", 1).unwrap().unwrap();
        assert_eq!(idx.num_rows(), 1);

        // Window (2,0) should have 1 row from partition 2
        let idx = store.get_position_index("2", 0).unwrap().unwrap();
        assert_eq!(idx.num_rows(), 1);
    }

    #[tokio::test(flavor = "multi_thread")]
    async fn test_loader_basic_v2() {
        let dir = tempfile::tempdir().unwrap();
        let ctx = SessionContext::new();

        let (schema, table) = make_test_table();
        ctx.register_table("source", Arc::new(table)).unwrap();

        let loader = CacheLoader::new("source", dir.path().to_str().unwrap())
            .with_window_size(1_000_000)
            .with_format_version(FORMAT_V2);
        let stats = loader.load(&ctx).await.unwrap();

        assert_eq!(stats.total_windows, 3);

        let store = VepKvStore::open(dir.path()).unwrap();
        assert_eq!(store.format_version(), FORMAT_V2);
        assert_eq!(store.schema().as_ref(), schema.as_ref());

        let idx = store.get_position_index("1", 0).unwrap().unwrap();
        assert_eq!(idx.num_rows(), 2);

        let batch = store.get_window_block_v2("1", 0).unwrap().unwrap();
        assert_eq!(batch.num_rows(), 2);
        assert_eq!(batch.num_columns(), schema.fields().len());
    }

    #[tokio::test(flavor = "multi_thread")]
    async fn test_loader_cross_window_indel_concurrent() {
        let dir = tempfile::tempdir().unwrap();
        let ctx = SessionContext::new();

        let schema = Arc::new(Schema::new(vec![
            Field::new("chrom", DataType::Utf8, false),
            Field::new("start", DataType::Int64, false),
            Field::new("end", DataType::Int64, false),
            Field::new("variation_name", DataType::Utf8, true),
            Field::new("allele_string", DataType::Utf8, false),
        ]));

        // Partition 1: has a cross-window indel on chrom 1 (window 0 -> 1)
        let batch1 = RecordBatch::try_new(
            schema.clone(),
            vec![
                Arc::new(StringArray::from(vec!["1", "1"])),
                Arc::new(Int64Array::from(vec![100, 999_990])),
                Arc::new(Int64Array::from(vec![100, 1_000_010])), // crosses into window 1
                Arc::new(StringArray::from(vec!["rs1", "indel1"])),
                Arc::new(StringArray::from(vec!["A/G", "ATCGATCG/-"])),
            ],
        )
        .unwrap();

        // Partition 2: also has variants in window 1 of chrom 1
        let batch2 = RecordBatch::try_new(
            schema.clone(),
            vec![
                Arc::new(StringArray::from(vec!["1"])),
                Arc::new(Int64Array::from(vec![1_000_050])),
                Arc::new(Int64Array::from(vec![1_000_050])),
                Arc::new(StringArray::from(vec!["rs5"])),
                Arc::new(StringArray::from(vec!["G/C"])),
            ],
        )
        .unwrap();

        let table = MemTable::try_new(schema.clone(), vec![vec![batch1], vec![batch2]]).unwrap();
        ctx.register_table("source", Arc::new(table)).unwrap();

        let loader =
            CacheLoader::new("source", dir.path().to_str().unwrap()).with_window_size(1_000_000);
        let stats = loader.load(&ctx).await.unwrap();

        let store = VepKvStore::open(dir.path()).unwrap();

        // Window (1,0) should have rs1 + indel1 = 2 rows
        let idx = store.get_position_index("1", 0).unwrap().unwrap();
        assert_eq!(idx.num_rows(), 2);

        // Window (1,1) should have indel1 (cross-window) + rs5 = 2 rows
        let idx = store.get_position_index("1", 1).unwrap().unwrap();
        assert_eq!(idx.num_rows(), 2);

        // Total: 2 windows
        assert_eq!(stats.total_windows, 2);
    }

    #[tokio::test(flavor = "multi_thread")]
    async fn test_loader_v1_accepts_max_supported_columns() {
        let dir = tempfile::tempdir().unwrap();
        let ctx = SessionContext::new();

        let max_cols = EntryType::MAX_COLUMN_INDEX as usize + 1;
        let (schema, table) = make_wide_test_table(max_cols);
        ctx.register_table("source", Arc::new(table)).unwrap();

        let loader =
            CacheLoader::new("source", dir.path().to_str().unwrap()).with_window_size(1_000_000);
        let stats = loader.load(&ctx).await.unwrap();
        assert_eq!(stats.total_variants, 1);

        let store = VepKvStore::open(dir.path()).unwrap();
        assert_eq!(store.schema().fields().len(), schema.fields().len());

        let cols = store
            .get_columns("1", 0, &[0, EntryType::MAX_COLUMN_INDEX])
            .unwrap()
            .unwrap();
        assert_eq!(cols.len(), 2);
    }

    #[tokio::test(flavor = "multi_thread")]
    async fn test_loader_v1_rejects_too_wide_schema() {
        let dir = tempfile::tempdir().unwrap();
        let ctx = SessionContext::new();

        let too_wide_cols = EntryType::MAX_COLUMN_INDEX as usize + 2;
        let (_, table) = make_wide_test_table(too_wide_cols);
        ctx.register_table("source", Arc::new(table)).unwrap();

        let loader =
            CacheLoader::new("source", dir.path().to_str().unwrap()).with_window_size(1_000_000);
        let err = loader.load(&ctx).await.unwrap_err().to_string();
        assert!(err.contains("supports at most"));
    }

    #[tokio::test(flavor = "multi_thread")]
    async fn test_loader_basic_v3() {
        let dir = tempfile::tempdir().unwrap();
        let ctx = SessionContext::new();

        let (schema, table) = make_test_table();
        ctx.register_table("source", Arc::new(table)).unwrap();

        let loader = CacheLoader::new("source", dir.path().to_str().unwrap())
            .with_window_size(1_000_000)
            .with_format_version(FORMAT_V3);
        let stats = loader.load(&ctx).await.unwrap();

        assert_eq!(stats.total_windows, 3);

        let store = VepKvStore::open(dir.path()).unwrap();
        assert_eq!(store.format_version(), FORMAT_V3);
        assert_eq!(store.schema().as_ref(), schema.as_ref());

        let idx = store.get_position_index("1", 0).unwrap().unwrap();
        assert_eq!(idx.num_rows(), 2);

        let batch = store.get_window_batch_v3("1", 0).unwrap().unwrap();
        assert_eq!(batch.num_rows(), 2);
        assert_eq!(batch.num_columns(), schema.fields().len());
    }

    #[tokio::test(flavor = "multi_thread")]
    async fn test_loader_v3_with_lz4() {
        let dir = tempfile::tempdir().unwrap();
        let ctx = SessionContext::new();

        let (schema, table) = make_test_table();
        ctx.register_table("source", Arc::new(table)).unwrap();

        let loader = CacheLoader::new("source", dir.path().to_str().unwrap())
            .with_window_size(1_000_000)
            .with_format_version(FORMAT_V3)
            .with_v2_block_codec(WindowBlockCodec::Lz4);
        let stats = loader.load(&ctx).await.unwrap();

        assert_eq!(stats.total_windows, 3);

        let store = VepKvStore::open(dir.path()).unwrap();
        assert_eq!(store.format_version(), FORMAT_V3);

        let batch = store.get_window_batch_v3("1", 0).unwrap().unwrap();
        assert_eq!(batch.num_rows(), 2);
        assert_eq!(batch.num_columns(), schema.fields().len());

        // Verify column projection works
        let cols = store
            .get_window_columns_v3("1", 0, &[0, 3])
            .unwrap()
            .unwrap();
        assert_eq!(cols.len(), 2);
        assert_eq!(cols[0].len(), 2);
    }

    #[tokio::test(flavor = "multi_thread")]
    async fn test_loader_basic_v4() {
        let dir = tempfile::tempdir().unwrap();
        let ctx = SessionContext::new();

        let (schema, table) = make_test_table();
        ctx.register_table("source", Arc::new(table)).unwrap();

        let loader = CacheLoader::new("source", dir.path().to_str().unwrap())
            .with_window_size(1_000_000)
            .with_format_version(FORMAT_V4);
        let stats = loader.load(&ctx).await.unwrap();

        assert_eq!(stats.total_windows, 3);

        let store = VepKvStore::open(dir.path()).unwrap();
        assert_eq!(store.format_version(), FORMAT_V4);
        assert_eq!(store.schema().as_ref(), schema.as_ref());

        let idx = store.get_position_index("1", 0).unwrap().unwrap();
        assert_eq!(idx.num_rows(), 2);

        let batch = store.get_window_batch_v4("1", 0).unwrap().unwrap();
        assert_eq!(batch.num_rows(), 2);
        assert_eq!(batch.num_columns(), schema.fields().len());
    }

    #[tokio::test(flavor = "multi_thread")]
    async fn test_loader_v4_with_lz4() {
        let dir = tempfile::tempdir().unwrap();
        let ctx = SessionContext::new();

        let (schema, table) = make_test_table();
        ctx.register_table("source", Arc::new(table)).unwrap();

        let loader = CacheLoader::new("source", dir.path().to_str().unwrap())
            .with_window_size(1_000_000)
            .with_format_version(FORMAT_V4)
            .with_v2_block_codec(WindowBlockCodec::Lz4);
        let stats = loader.load(&ctx).await.unwrap();

        assert_eq!(stats.total_windows, 3);

        let store = VepKvStore::open(dir.path()).unwrap();
        assert_eq!(store.format_version(), FORMAT_V4);

        let batch = store.get_window_batch_v4("1", 0).unwrap().unwrap();
        assert_eq!(batch.num_rows(), 2);
        assert_eq!(batch.num_columns(), schema.fields().len());

        // Verify column projection works
        let cols = store
            .get_window_columns_v4("1", 0, &[0, 3])
            .unwrap()
            .unwrap();
        assert_eq!(cols.len(), 2);
        assert_eq!(cols[0].len(), 2);
    }

    #[tokio::test(flavor = "multi_thread")]
    async fn test_loader_basic_v5() {
        let dir = tempfile::tempdir().unwrap();
        let ctx = SessionContext::new();

        let (schema, table) = make_test_table();
        ctx.register_table("source", Arc::new(table)).unwrap();

        let loader =
            CacheLoader::new("source", dir.path().to_str().unwrap()).with_format_version(FORMAT_V5);
        let stats = loader.load(&ctx).await.unwrap();

        // 4 rows -> 3 unique positions: (1,100,100), (1,500000,500000), (1,1500000,1500000), (2,100,100)
        assert_eq!(stats.total_variants, 4);
        // "windows" for v5 = positions written
        assert!(stats.total_windows > 0);

        let store = VepKvStore::open(dir.path()).unwrap();
        assert_eq!(store.format_version(), FORMAT_V5);
        assert_eq!(store.schema().as_ref(), schema.as_ref());

        // Verify position entries exist (using decompressed read since loader now compresses).
        let chrom_code_1 = crate::key_encoding::chrom_to_code("1");
        let chrom_code_2 = crate::key_encoding::chrom_to_code("2");

        // Test data has too few entries for dictionary training; dict is absent.
        // With real data (100+ positions), has_v5_dict() would be true.

        let entry = store
            .get_position_entry_decompressed(chrom_code_1, 100, 100)
            .unwrap();
        assert!(entry.is_some());

        let entry = store
            .get_position_entry_decompressed(chrom_code_1, 500_000, 500_000)
            .unwrap();
        assert!(entry.is_some());

        let entry = store
            .get_position_entry_decompressed(chrom_code_2, 100, 100)
            .unwrap();
        assert!(entry.is_some());

        // Missing position.
        let entry = store
            .get_position_entry_decompressed(chrom_code_1, 999, 999)
            .unwrap();
        assert!(entry.is_none());

        // Verify a position entry can be deserialized after decompression.
        let raw = store
            .get_position_entry_decompressed(chrom_code_1, 100, 100)
            .unwrap()
            .unwrap();
        let reader = crate::position_entry::PositionEntryReader::new(&raw).unwrap();
        // Position (1, 100, 100) has 2 rows in test data with same (start, end).
        // Actually it has just 1 row at (1, 100, 100) in the test data.
        assert!(reader.num_alleles() >= 1);
        assert!(reader.num_cols() > 0);
    }
}
