//! Parquet -> fjall cache ingestion pipeline (streaming).
//!
//! Loads a VEP variation cache from Parquet into a fjall KV store,
//! storing one zstd-compressed entry per `(chrom, start)` position and
//! parallelizing across DataFusion physical partitions.
//! All variants at the same `(chrom, start)` — regardless of `end` —
//! are merged into a single entry with `end` stored as a regular column.
//!
//! Streaming approach:
//! 1. A sample of up to 10,000 rows is read from the source table to train
//!    a single zstd dictionary. The dictionary is stored in KV meta and shared
//!    across all partitions — this avoids the race condition where parallel
//!    partitions would each train and overwrite different dictionaries.
//! 2. Each partition processes its stream batch-by-batch, compressing entries
//!    with the shared dictionary. Memory stays bounded to one batch per partition.
//! 3. If training fails (too few samples), entries are stored uncompressed.

use std::collections::HashMap;
use std::path::Path;
use std::sync::Arc;
use std::time::Instant;

use datafusion::arrow::array::{
    Array, AsArray, LargeStringArray, RecordBatch, StringArray, StringViewArray,
};
use datafusion::arrow::datatypes::{DataType, Int64Type, Schema, SchemaRef};
use datafusion::common::{DataFusionError, Result};
use datafusion::physical_plan::{ExecutionPlan, ExecutionPlanProperties};
use datafusion::prelude::SessionContext;
use futures::StreamExt;
use log::{debug, info};

use super::key_encoding::{chrom_to_code, encode_position_key};
use super::kv_store::{VepKvStore, decompress_into_buffer_with_retry};
use super::position_entry::{PositionEntryReader, make_builder, serialize_position_entry};

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
    /// Pre-trained dictionary bytes shared by all partitions. `None` if training
    /// was skipped (too few samples) — partitions store entries uncompressed.
    dict: Option<Arc<Vec<u8>>>,
    zstd_level: i32,
    /// Serializes read-merge-write operations across partitions to prevent
    /// concurrent merges from silently dropping alleles when the same
    /// `(chrom, start)` key appears in multiple partitions (e.g., at
    /// partition boundaries for sorted input).
    merge_lock: Arc<std::sync::Mutex<()>>,
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
    /// First trains a zstd dictionary from a sample of the source data, then
    /// spawns one streaming task per physical partition. Each task reads batches
    /// and flushes position entries inline — memory stays bounded to one batch
    /// at a time per partition.
    pub async fn load(&self, session: &SessionContext) -> Result<LoadStats> {
        let start_time = Instant::now();

        let source_df = session.table(&self.source_table).await?;
        let schema: SchemaRef = Arc::new(source_df.schema().as_arrow().clone());

        let store = Arc::new(VepKvStore::create(
            Path::new(&self.target_path),
            schema.clone(),
        )?);

        // Phase 1: Train zstd dictionary from a sample of the source data.
        // This happens once before any parallel work, ensuring all partitions
        // share the same dictionary.
        let dict = self.train_dictionary(session, &store, &schema).await?;

        // Phase 2: Stream all partitions in parallel with the shared dictionary.
        let source_df = session.table(&self.source_table).await?;
        let plan = source_df.create_physical_plan().await?;
        let partition_count = plan.output_partitioning().partition_count();
        let task_ctx = session.task_ctx();

        info!(
            "Loading into KV cache: {} partitions, zstd_level={}, dict={}",
            partition_count,
            self.zstd_level,
            if dict.is_some() { "trained" } else { "none" },
        );
        debug!(
            "Physical plan:\n{}",
            datafusion::physical_plan::displayable(plan.as_ref()).indent(true)
        );

        let semaphore = self
            .parallelism
            .map(|n| Arc::new(tokio::sync::Semaphore::new(n.min(partition_count))));
        let merge_lock = Arc::new(std::sync::Mutex::new(()));

        let mut handles = Vec::with_capacity(partition_count);

        for partition_id in 0..partition_count {
            let plan = Arc::clone(&plan);
            let task_ctx = Arc::clone(&task_ctx);
            let sem = semaphore.clone();
            let ctx = StreamContext {
                store: Arc::clone(&store),
                schema: schema.clone(),
                dict: dict.clone(),
                zstd_level: self.zstd_level,
                merge_lock: Arc::clone(&merge_lock),
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

        // Major-compact the data keyspace so the LSM tree is fully optimized
        // for reads (bloom filters built, levels merged) before first use.
        info!("Running major compaction on data keyspace...");
        let compact_start = Instant::now();
        store
            .data_partition()
            .major_compact()
            .map_err(|e| DataFusionError::External(Box::new(e)))?;
        info!(
            "Major compaction completed in {:.1}s",
            compact_start.elapsed().as_secs_f64()
        );

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

    /// Train a zstd dictionary from a sample of the source table.
    ///
    /// Reads up to 10,000 rows via `LIMIT`, serializes them as position entries,
    /// and trains a dictionary. Stores the dictionary and compression level in
    /// the KV store meta partition. Returns `None` if too few samples.
    async fn train_dictionary(
        &self,
        session: &SessionContext,
        store: &VepKvStore,
        schema: &SchemaRef,
    ) -> Result<Option<Arc<Vec<u8>>>> {
        let sample_df = session
            .sql(&format!(
                "SELECT * FROM `{}` LIMIT 10000",
                self.source_table
            ))
            .await?;
        let batches = sample_df.collect().await?;

        let chrom_col_idx = schema.index_of("chrom")?;
        let start_col_idx = schema.index_of("start")?;
        let allele_col_idx = schema.index_of("allele_string")?;
        // `end` is stored as a regular column inside the entry (not part of the key).
        let col_indices: Vec<usize> = (0..schema.fields().len())
            .filter(|&i| i != chrom_col_idx && i != start_col_idx)
            .collect();

        // Serialize position entries from sample batches.
        let mut samples: Vec<Vec<u8>> = Vec::new();
        for batch in &batches {
            if batch.num_rows() == 0 {
                continue;
            }
            let chrom_col = batch.column(chrom_col_idx);
            let starts = batch.column(start_col_idx).as_primitive::<Int64Type>();

            let mut groups: HashMap<(String, i64), Vec<usize>> = HashMap::new();
            for row in 0..batch.num_rows() {
                let chrom = string_value(chrom_col.as_ref(), row);
                let start = starts.value(row);
                groups
                    .entry((chrom.to_string(), start))
                    .or_default()
                    .push(row);
            }

            for rows in groups.values() {
                let value = serialize_position_entry(rows, batch, &col_indices, allele_col_idx)?;
                samples.push(value);
            }
        }

        let min_training_samples = 100;
        if samples.len() < min_training_samples {
            info!(
                "only {} sample entries, skipping dictionary training (need >= {})",
                samples.len(),
                min_training_samples
            );
            return Ok(None);
        }

        let max_samples = 10_000.min(samples.len());
        let sample_refs: Vec<&[u8]> = samples[..max_samples]
            .iter()
            .map(|v| v.as_slice())
            .collect();

        let dict_size_kb = self.dict_size_kb;
        let zstd_level = self.zstd_level;

        match train_dict(&sample_refs, dict_size_kb) {
            Ok(dict) => {
                info!(
                    "zstd dict trained: {} samples, dict {} bytes",
                    max_samples,
                    dict.len()
                );
                store.store_dict(&dict)?;
                store.store_zstd_level(zstd_level)?;
                Ok(Some(Arc::new(dict)))
            }
            Err(e) => {
                info!("zstd dict training failed (falling back to uncompressed): {e}");
                Ok(None)
            }
        }
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

    // Create compressor from shared dictionary (if available).
    // Each partition gets its own compressor instance but all share the same dict.
    let mut compressor = if let Some(ref dict) = ctx.dict {
        Some(
            zstd::bulk::Compressor::with_dictionary(ctx.zstd_level, dict).map_err(|e| {
                DataFusionError::Execution(format!("failed to create zstd compressor: {e}"))
            })?,
        )
    } else {
        None
    };
    let mut decompressor = if let Some(ref dict) = ctx.dict {
        Some(
            zstd::bulk::Decompressor::with_dictionary(dict).map_err(|e| {
                DataFusionError::Execution(format!("failed to create zstd decompressor: {e}"))
            })?,
        )
    } else {
        None
    };

    while let Some(batch_result) = stream.next().await {
        let batch = batch_result?;
        if batch.num_rows() == 0 {
            continue;
        }

        stats.total_variants += batch.num_rows() as u64;

        let store = Arc::clone(&ctx.store);
        let schema = ctx.schema.clone();
        let batch_clone = batch.clone();

        if let Some(mut comp) = compressor.take() {
            let mut dec = decompressor.take().ok_or_else(|| {
                DataFusionError::Execution(
                    "missing zstd decompressor while compression is enabled".to_string(),
                )
            })?;
            let lock = Arc::clone(&ctx.merge_lock);
            let (comp_back, dec_back, positions, bytes) = tokio::task::spawn_blocking(move || {
                let result = flush_positions_compressed(
                    &store,
                    &schema,
                    &batch_clone,
                    &mut comp,
                    &mut dec,
                    &lock,
                );
                result.map(|(p, b)| (comp, dec, p, b))
            })
            .await
            .map_err(|e| DataFusionError::External(Box::new(e)))??;
            compressor = Some(comp_back);
            decompressor = Some(dec_back);
            stats.total_positions += positions;
            stats.total_bytes += bytes;
        } else {
            let lock = Arc::clone(&ctx.merge_lock);
            let (positions, bytes) = tokio::task::spawn_blocking(move || {
                flush_positions_uncompressed(&store, &schema, &batch_clone, &lock)
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

/// Flush position entries without compression (fallback when dictionary training was skipped).
fn flush_positions_uncompressed(
    store: &VepKvStore,
    schema: &SchemaRef,
    batch: &RecordBatch,
    merge_lock: &std::sync::Mutex<()>,
) -> Result<(u64, u64)> {
    let chrom_col_idx = schema.index_of("chrom")?;
    let start_col_idx = schema.index_of("start")?;
    let allele_col_idx = schema.index_of("allele_string")?;

    let chrom_col = batch.column(chrom_col_idx);
    let starts = batch.column(start_col_idx).as_primitive::<Int64Type>();

    // `end` is stored as a regular column inside the entry.
    let col_indices: Vec<usize> = (0..schema.fields().len())
        .filter(|&i| i != chrom_col_idx && i != start_col_idx)
        .collect();

    let mut groups: HashMap<(String, i64), Vec<usize>> = HashMap::new();
    for row in 0..batch.num_rows() {
        let chrom = string_value(chrom_col.as_ref(), row);
        let start = starts.value(row);
        groups
            .entry((chrom.to_string(), start))
            .or_default()
            .push(row);
    }

    let mut total_bytes = 0u64;
    let mut num_positions = 0u64;

    for ((chrom, start), rows) in &groups {
        let mut value = serialize_position_entry(rows, batch, &col_indices, allele_col_idx)?;
        // Read-merge-write must be atomic: hold the lock through the write
        // to prevent concurrent partitions from overwriting each other's
        // alleles at shared boundary positions.
        let _guard = merge_lock.lock().unwrap();
        if let Some(existing) =
            store.get_position_entry_decompressed(chrom_to_code(chrom), *start)?
        {
            value = merge_position_entries(&existing, &value, schema)?;
        }
        store.put_position_entry(chrom, *start, &value)?;
        drop(_guard);
        total_bytes += value.len() as u64;
        num_positions += 1;
    }

    Ok((num_positions, total_bytes))
}

/// Compress and flush position entries using a pre-trained zstd dictionary compressor.
fn flush_positions_compressed(
    store: &VepKvStore,
    schema: &SchemaRef,
    batch: &RecordBatch,
    compressor: &mut zstd::bulk::Compressor<'_>,
    decompressor: &mut zstd::bulk::Decompressor<'_>,
    merge_lock: &std::sync::Mutex<()>,
) -> Result<(u64, u64)> {
    let chrom_col_idx = schema.index_of("chrom")?;
    let start_col_idx = schema.index_of("start")?;
    let allele_col_idx = schema.index_of("allele_string")?;

    let chrom_col = batch.column(chrom_col_idx);
    let starts = batch.column(start_col_idx).as_primitive::<Int64Type>();

    // `end` is stored as a regular column inside the entry.
    let col_indices: Vec<usize> = (0..schema.fields().len())
        .filter(|&i| i != chrom_col_idx && i != start_col_idx)
        .collect();

    let mut groups: HashMap<(String, i64), Vec<usize>> = HashMap::new();
    for row in 0..batch.num_rows() {
        let chrom = string_value(chrom_col.as_ref(), row);
        let start = starts.value(row);
        groups
            .entry((chrom.to_string(), start))
            .or_default()
            .push(row);
    }

    let mut total_bytes = 0u64;
    let mut num_positions = 0u64;

    for ((chrom, start), rows) in &groups {
        let chrom_code = chrom_to_code(chrom);
        let mut raw_value = serialize_position_entry(rows, batch, &col_indices, allele_col_idx)?;
        // Read-merge-write must be atomic: hold the lock through the write
        // to prevent concurrent partitions from overwriting each other's
        // alleles at shared boundary positions.
        let _guard = merge_lock.lock().unwrap();
        if let Some(existing_compressed) = store.get_position_entry(chrom_code, *start)? {
            let mut existing_raw = Vec::new();
            decompress_into_buffer_with_retry(
                decompressor,
                &existing_compressed,
                &mut existing_raw,
                "zstd decompression failed while merging existing position entry",
            )?;
            raw_value = merge_position_entries(&existing_raw, &raw_value, schema)?;
        }
        let compressed = compressor
            .compress(&raw_value)
            .map_err(|e| DataFusionError::Execution(format!("zstd compression failed: {e}")))?;
        store.put_position_entry(chrom, *start, &compressed)?;
        drop(_guard);
        total_bytes += compressed.len() as u64;
        num_positions += 1;
    }

    Ok((num_positions, total_bytes))
}

/// Merge two serialized entries for the same genomic position.
///
/// This preserves all allele rows when the same `(chrom,start)` appears in
/// multiple input batches/partitions. Entries with different `end` values
/// are merged together since `end` is stored as a regular column.
fn merge_position_entries(existing: &[u8], incoming: &[u8], schema: &SchemaRef) -> Result<Vec<u8>> {
    let existing_reader = PositionEntryReader::new(existing)?;
    let incoming_reader = PositionEntryReader::new(incoming)?;

    if existing_reader.num_cols() != incoming_reader.num_cols() {
        return Err(DataFusionError::Execution(format!(
            "cannot merge position entries with different column counts: existing={} incoming={}",
            existing_reader.num_cols(),
            incoming_reader.num_cols()
        )));
    }

    let chrom_col_idx = schema.index_of("chrom")?;
    let start_col_idx = schema.index_of("start")?;
    // `end` is stored inside the entry as a regular column, not excluded.
    let stored_col_indices: Vec<usize> = (0..schema.fields().len())
        .filter(|&i| i != chrom_col_idx && i != start_col_idx)
        .collect();

    if stored_col_indices.len() != existing_reader.num_cols() {
        return Err(DataFusionError::Execution(format!(
            "schema/entry mismatch while merging position entries: schema stored cols={} entry cols={}",
            stored_col_indices.len(),
            existing_reader.num_cols()
        )));
    }

    let total_rows = existing_reader.num_alleles() + incoming_reader.num_alleles();
    let existing_rows: Vec<usize> = (0..existing_reader.num_alleles()).collect();
    let incoming_rows: Vec<usize> = (0..incoming_reader.num_alleles()).collect();

    let mut merged_columns = Vec::with_capacity(stored_col_indices.len());
    let mut merged_fields = Vec::with_capacity(stored_col_indices.len());
    for (entry_col_idx, &schema_col_idx) in stored_col_indices.iter().enumerate() {
        let src_field = schema.field(schema_col_idx);
        let normalized_dt = match src_field.data_type() {
            DataType::Utf8View | DataType::LargeUtf8 => DataType::Utf8,
            other => other.clone(),
        };
        let field = datafusion::arrow::datatypes::Field::new(
            src_field.name(),
            normalized_dt,
            src_field.is_nullable(),
        );
        let mut builder = make_builder(field.data_type(), total_rows)?;
        existing_reader.append_column_values(entry_col_idx, &existing_rows, builder.as_mut())?;
        incoming_reader.append_column_values(entry_col_idx, &incoming_rows, builder.as_mut())?;
        merged_columns.push(builder.finish());
        merged_fields.push(field);
    }

    let merged_schema = Arc::new(Schema::new(merged_fields));
    let merged_batch = RecordBatch::try_new(merged_schema.clone(), merged_columns)
        .map_err(|e| DataFusionError::ArrowError(Box::new(e), None))?;

    let allele_col_idx = merged_schema.index_of("allele_string").map_err(|_| {
        DataFusionError::Execution(
            "cannot merge position entries: allele_string column missing".to_string(),
        )
    })?;
    let col_indices: Vec<usize> = (0..merged_schema.fields().len()).collect();
    let rows: Vec<usize> = (0..total_rows).collect();
    serialize_position_entry(&rows, &merged_batch, &col_indices, allele_col_idx)
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

/// Extract a string value from a Utf8, LargeUtf8, or Utf8View array.
fn string_value(col: &dyn Array, row: usize) -> &str {
    if let Some(arr) = col.as_any().downcast_ref::<StringArray>() {
        arr.value(row)
    } else if let Some(arr) = col.as_any().downcast_ref::<StringViewArray>() {
        arr.value(row)
    } else if let Some(arr) = col.as_any().downcast_ref::<LargeStringArray>() {
        arr.value(row)
    } else {
        panic!("Expected Utf8, LargeUtf8, or Utf8View column for chrom")
    }
}

/// Load translation features into the fjall database's translations keyspace.
///
/// Opens (or creates) the translations keyspace and inserts each feature
/// keyed by `transcript_id`.
pub fn load_translations_into_kv(
    db: &fjall::Database,
    translations: &[crate::transcript_consequence::TranslationFeature],
) -> Result<usize> {
    use super::context_store::TranslationKvStore;

    let store = TranslationKvStore::create(db)?;
    for t in translations {
        store.put(t)?;
    }
    info!(
        "Loaded {} translations into fjall translations keyspace",
        translations.len()
    );
    Ok(translations.len())
}

/// Load exon features into the fjall database's exons keyspace.
///
/// Groups exons by `transcript_id`, sorts each group by `exon_number`,
/// then stores one entry per transcript.
pub fn load_exons_into_kv(
    db: &fjall::Database,
    exons: &[crate::transcript_consequence::ExonFeature],
) -> Result<usize> {
    use super::context_store::ExonKvStore;

    let store = ExonKvStore::create(db)?;

    // Group by transcript_id.
    let mut by_tx: HashMap<&str, Vec<crate::transcript_consequence::ExonFeature>> = HashMap::new();
    for e in exons {
        by_tx
            .entry(e.transcript_id.as_str())
            .or_default()
            .push(e.clone());
    }
    let tx_count = by_tx.len();
    for (tx_id, mut group) in by_tx {
        store.put(tx_id, &mut group)?;
    }
    info!(
        "Loaded {} exons for {} transcripts into fjall exons keyspace",
        exons.len(),
        tx_count
    );
    Ok(exons.len())
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::kv_cache::kv_store::FORMAT_V0;
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

        // 4 rows -> 4 unique positions: (1,100), (1,500000), (1,1500000), (2,100)
        assert_eq!(stats.total_variants, 4);
        assert!(stats.total_positions > 0);

        let store = VepKvStore::open(dir.path()).unwrap();
        assert_eq!(store.format_version(), FORMAT_V0);
        assert_eq!(store.schema().as_ref(), schema.as_ref());

        // Verify position entries exist.
        let chrom_code_1 = crate::kv_cache::key_encoding::chrom_to_code("1");
        let chrom_code_2 = crate::kv_cache::key_encoding::chrom_to_code("2");

        let entry = store
            .get_position_entry_decompressed(chrom_code_1, 100)
            .unwrap();
        assert!(entry.is_some());

        let entry = store
            .get_position_entry_decompressed(chrom_code_1, 500_000)
            .unwrap();
        assert!(entry.is_some());

        let entry = store
            .get_position_entry_decompressed(chrom_code_2, 100)
            .unwrap();
        assert!(entry.is_some());

        // Missing position.
        let entry = store
            .get_position_entry_decompressed(chrom_code_1, 999)
            .unwrap();
        assert!(entry.is_none());

        // Verify a position entry can be deserialized after decompression.
        let raw = store
            .get_position_entry_decompressed(chrom_code_1, 100)
            .unwrap()
            .unwrap();
        let reader = crate::kv_cache::position_entry::PositionEntryReader::new(&raw).unwrap();
        assert!(reader.num_alleles() >= 1);
        assert!(reader.num_cols() > 0);
    }
}
