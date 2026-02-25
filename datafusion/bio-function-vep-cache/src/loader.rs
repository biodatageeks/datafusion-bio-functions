//! Parquet -> fjall cache ingestion pipeline (streaming).
//!
//! Loads a VEP variation cache from Parquet into a fjall KV store,
//! storing one zstd-compressed entry per genomic position and
//! parallelizing across DataFusion physical partitions.
//!
//! Streaming approach:
//! - Each partition processes its stream batch-by-batch.
//! - The first batch trains a zstd dictionary from serialized position entries,
//!   then all entries (first + subsequent) are compressed with that dictionary.
//! - If training fails (too few samples), entries are stored uncompressed.

use std::collections::HashMap;
use std::path::Path;
use std::sync::Arc;
use std::time::Instant;

use datafusion::arrow::array::{Array, AsArray, RecordBatch, StringArray, StringViewArray};
use datafusion::arrow::datatypes::{Int64Type, SchemaRef};
use datafusion::common::{DataFusionError, Result};
use datafusion::physical_plan::{ExecutionPlan, ExecutionPlanProperties};
use datafusion::prelude::SessionContext;
use futures::StreamExt;
use log::{debug, info};

use crate::key_encoding::encode_position_key;
use crate::kv_store::VepKvStore;
use crate::position_entry::serialize_position_entry;

/// Statistics returned after loading.
#[derive(Debug, Clone)]
pub struct LoadStats {
    pub total_variants: u64,
    pub total_positions: u64,
    pub total_bytes: u64,
    pub elapsed_secs: f64,
}

/// Per-partition stats accumulated during streaming.
struct PartitionStats {
    total_variants: u64,
    total_positions: u64,
    total_bytes: u64,
}

/// Shared context passed to each partition's streaming task.
struct StreamContext {
    store: Arc<VepKvStore>,
    schema: SchemaRef,
    zstd_level: i32,
    dict_size_kb: u32,
}

/// Builder for loading a Parquet-based VEP cache into fjall.
pub struct CacheLoader {
    source_table: String,
    target_path: String,
    parallelism: Option<usize>,
    zstd_level: i32,
    dict_size_kb: u32,
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
            parallelism: None,
            zstd_level: 3,
            dict_size_kb: 112,
        }
    }

    /// Set an optional concurrency cap on partition-level parallelism.
    /// Default is `None` — all partitions run concurrently.
    pub fn with_parallelism(mut self, n: usize) -> Self {
        self.parallelism = Some(n);
        self
    }

    /// Set the zstd compression level (default: 3).
    ///
    /// Higher levels produce smaller caches at the cost of slower writes.
    /// Decompression speed is constant regardless of level.
    pub fn with_zstd_level(mut self, level: i32) -> Self {
        self.zstd_level = level;
        self
    }

    /// Set the zstd dictionary size in KB (default: 112).
    pub fn with_dict_size_kb(mut self, size_kb: u32) -> Self {
        self.dict_size_kb = size_kb;
        self
    }

    /// Load the cache from the registered source table into fjall.
    ///
    /// Spawns one streaming task per physical partition. Each task reads
    /// batches and flushes position entries inline — memory stays bounded
    /// to one batch at a time per partition.
    pub async fn load(&self, session: &SessionContext) -> Result<LoadStats> {
        let start_time = Instant::now();

        let source_df = session.table(&self.source_table).await?;
        let schema: SchemaRef = Arc::new(source_df.schema().as_arrow().clone());

        let store = Arc::new(VepKvStore::create(
            Path::new(&self.target_path),
            schema.clone(),
        )?);

        let plan = source_df.create_physical_plan().await?;
        let partition_count = plan.output_partitioning().partition_count();
        let task_ctx = session.task_ctx();

        info!(
            "Loading into KV cache: {} partitions, zstd_level={}, dict_size_kb={}",
            partition_count, self.zstd_level, self.dict_size_kb,
        );
        debug!(
            "Physical plan:\n{}",
            datafusion::physical_plan::displayable(plan.as_ref()).indent(true)
        );

        let semaphore = self
            .parallelism
            .map(|n| Arc::new(tokio::sync::Semaphore::new(n.min(partition_count))));

        let mut handles = Vec::with_capacity(partition_count);

        for partition_id in 0..partition_count {
            let plan = Arc::clone(&plan);
            let task_ctx = Arc::clone(&task_ctx);
            let sem = semaphore.clone();
            let ctx = StreamContext {
                store: Arc::clone(&store),
                schema: schema.clone(),
                zstd_level: self.zstd_level,
                dict_size_kb: self.dict_size_kb,
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
        let mut total_positions = 0u64;
        let mut total_bytes = 0u64;

        for handle in handles {
            let ps = handle
                .await
                .map_err(|e| DataFusionError::External(Box::new(e)))??;
            total_variants += ps.total_variants;
            total_positions += ps.total_positions;
            total_bytes += ps.total_bytes;
        }

        store.persist()?;

        let elapsed = start_time.elapsed().as_secs_f64();
        info!(
            "Loaded {total_variants} variants into {total_positions} positions ({total_bytes} bytes) in {elapsed:.1}s"
        );

        Ok(LoadStats {
            total_variants,
            total_positions,
            total_bytes,
            elapsed_secs: elapsed,
        })
    }
}

/// Stream a single partition: read batches and flush position entries inline.
async fn stream_partition(
    plan: Arc<dyn ExecutionPlan>,
    task_ctx: Arc<datafusion::execution::TaskContext>,
    partition_id: usize,
    ctx: StreamContext,
) -> Result<PartitionStats> {
    let part_start = Instant::now();
    let mut stream = plan.execute(partition_id, task_ctx)?;

    let mut stats = PartitionStats {
        total_variants: 0,
        total_positions: 0,
        total_bytes: 0,
    };

    // Zstd dictionary compressor — built from first batch, reused for all subsequent.
    // None until the first batch trains the dictionary. If training is skipped (too few
    // samples), this stays None and all entries are stored uncompressed.
    let mut compressor: Option<zstd::bulk::Compressor<'static>> = None;
    let mut first_batch_done = false;

    while let Some(batch_result) = stream.next().await {
        let batch = batch_result?;
        if batch.num_rows() == 0 {
            continue;
        }

        stats.total_variants += batch.num_rows() as u64;

        let store = Arc::clone(&ctx.store);
        let schema = ctx.schema.clone();
        let batch_clone = batch.clone();

        if !first_batch_done {
            first_batch_done = true;
            let zstd_level = ctx.zstd_level;
            let dict_size_kb = ctx.dict_size_kb;
            // First batch: try to train dictionary, store it, build compressor.
            let (comp, positions, bytes) = tokio::task::spawn_blocking(move || {
                train_and_flush_first_batch(&store, &schema, &batch_clone, zstd_level, dict_size_kb)
            })
            .await
            .map_err(|e| DataFusionError::External(Box::new(e)))??;
            compressor = comp;
            stats.total_positions += positions;
            stats.total_bytes += bytes;
        } else if let Some(mut comp) = compressor.take() {
            // Subsequent batches with dictionary: compress and insert.
            let (comp_back, positions, bytes) = tokio::task::spawn_blocking(move || {
                let result = flush_positions_compressed(&store, &schema, &batch_clone, &mut comp);
                result.map(|(p, b)| (comp, p, b))
            })
            .await
            .map_err(|e| DataFusionError::External(Box::new(e)))??;
            compressor = Some(comp_back);
            stats.total_positions += positions;
            stats.total_bytes += bytes;
        } else {
            // Subsequent batches without dictionary: insert uncompressed.
            let (positions, bytes) = tokio::task::spawn_blocking(move || {
                flush_positions_uncompressed(&store, &schema, &batch_clone)
            })
            .await
            .map_err(|e| DataFusionError::External(Box::new(e)))??;
            stats.total_positions += positions;
            stats.total_bytes += bytes;
        }
    }

    info!(
        "Partition {partition_id}: {} variants, {} positions in {:.1}s",
        stats.total_variants,
        stats.total_positions,
        part_start.elapsed().as_secs_f64()
    );

    Ok(stats)
}

/// Train a zstd dictionary from the first batch, store it, then compress + insert all entries.
///
/// Returns `(Option<Compressor>, num_positions, total_bytes)`. The compressor is `None` if
/// there weren't enough samples to train a meaningful dictionary (entries stored uncompressed).
fn train_and_flush_first_batch(
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
        match train_dict(&samples, dict_size_kb) {
            Ok(dict) => {
                info!(
                    "zstd dict trained: {} samples, dict {} bytes",
                    max_samples,
                    dict.len()
                );
                Some(dict)
            }
            Err(e) => {
                info!("zstd dict training failed (falling back to uncompressed): {e}");
                None
            }
        }
    } else {
        info!(
            "only {} entries, skipping dictionary training (need >= {})",
            serialized.len(),
            min_training_samples
        );
        None
    };

    if let Some(ref dict) = dict_bytes {
        // Store dictionary and compression level in meta, then compress entries.
        store.store_dict(dict)?;
        store.store_zstd_level(zstd_level)?;

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

/// Flush position entries without compression (fallback when dictionary training was skipped).
fn flush_positions_uncompressed(
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

/// Compress and flush position entries using a pre-trained zstd dictionary compressor.
fn flush_positions_compressed(
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
fn train_dict(samples: &[&[u8]], dict_size_kb: u32) -> Result<Vec<u8>> {
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

#[cfg(test)]
mod tests {
    use super::*;
    use crate::kv_store::FORMAT_V0;
    use datafusion::arrow::array::{Int64Array, StringArray};
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

    #[tokio::test(flavor = "multi_thread")]
    async fn test_loader_basic() {
        let dir = tempfile::tempdir().unwrap();
        let ctx = SessionContext::new();

        let (schema, table) = make_test_table();
        ctx.register_table("source", Arc::new(table)).unwrap();

        let loader = CacheLoader::new("source", dir.path().to_str().unwrap());
        let stats = loader.load(&ctx).await.unwrap();

        // 4 rows -> 4 unique positions: (1,100,100), (1,500000,500000), (1,1500000,1500000), (2,100,100)
        assert_eq!(stats.total_variants, 4);
        assert!(stats.total_positions > 0);

        let store = VepKvStore::open(dir.path()).unwrap();
        assert_eq!(store.format_version(), FORMAT_V0);
        assert_eq!(store.schema().as_ref(), schema.as_ref());

        // Verify position entries exist.
        let chrom_code_1 = crate::key_encoding::chrom_to_code("1");
        let chrom_code_2 = crate::key_encoding::chrom_to_code("2");

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
        assert!(reader.num_alleles() >= 1);
        assert!(reader.num_cols() > 0);
    }
}
