//! Parquet -> fjall cache ingestion pipeline.
//!
//! Loads a VEP variation cache from Parquet into a fjall KV store,
//! grouping variants into windows and parallelizing across DataFusion
//! physical partitions.
//!
//! Two-phase approach:
//! 1. **Read phase** (parallel): each partition reads its stream and groups
//!    rows by (chrom, window_id) into a shared in-memory buffer.
//! 2. **Write phase** (sequential): each window is flushed to fjall exactly
//!    once — no read-modify-write overhead.

use std::collections::HashMap;
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
use log::info;

use crate::key_encoding::{DEFAULT_WINDOW_SIZE, window_id_for_position};
use crate::kv_store::{FORMAT_V1, VepKvStore, to_v1_column_index, validate_v1_schema_width};
use crate::position_index::PositionIndex;

/// Statistics returned after loading.
#[derive(Debug, Clone)]
pub struct LoadStats {
    pub total_variants: u64,
    pub total_windows: u64,
    pub total_bytes: u64,
    pub elapsed_secs: f64,
}

/// Thread-safe accumulator for window data from multiple partitions.
struct SharedWindowBuffers {
    buffers: std::sync::Mutex<HashMap<(String, u64), Vec<RecordBatch>>>,
}

impl SharedWindowBuffers {
    fn new() -> Self {
        Self {
            buffers: std::sync::Mutex::new(HashMap::new()),
        }
    }

    /// Push a sub-batch for a given (chrom, window_id).
    fn push(&self, chrom: String, wid: u64, batch: RecordBatch) {
        let mut map = self.buffers.lock().unwrap();
        map.entry((chrom, wid)).or_default().push(batch);
    }

    /// Drain all accumulated buffers.
    fn drain(&self) -> HashMap<(String, u64), Vec<RecordBatch>> {
        let mut map = self.buffers.lock().unwrap();
        std::mem::take(&mut *map)
    }
}

/// Builder for loading a Parquet-based VEP cache into fjall.
pub struct CacheLoader {
    source_table: String,
    target_path: String,
    window_size: u64,
    parallelism: Option<usize>,
    format_version: u8,
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

    /// Set the format version (default: FORMAT_V1 = columnar).
    pub fn with_format_version(mut self, version: u8) -> Self {
        self.format_version = version;
        self
    }

    /// Load the cache from the registered source table into fjall.
    ///
    /// Uses DataFusion's physical partitions for parallelism rather than
    /// per-chromosome queries. The number of partitions is controlled by
    /// the session's `target_partitions` config.
    pub async fn load(&self, session: &SessionContext) -> Result<LoadStats> {
        let start_time = Instant::now();

        // Get the source schema.
        let source_df = session.table(&self.source_table).await?;
        let schema: SchemaRef = Arc::new(source_df.schema().as_arrow().clone());

        if self.format_version >= FORMAT_V1 {
            validate_v1_schema_width(schema.fields().len())?;
        }

        // Create the KV store with the chosen format version.
        let store = Arc::new(VepKvStore::create_with_version(
            Path::new(&self.target_path),
            schema.clone(),
            self.window_size,
            self.format_version,
        )?);

        // Create physical plan — no WHERE, no ORDER BY.
        let plan = source_df.create_physical_plan().await?;
        let partition_count = plan.output_partitioning().partition_count();
        let task_ctx = session.task_ctx();

        info!(
            "Loading into KV cache: {} partitions, window_size={}",
            partition_count, self.window_size
        );

        // Phase 1: Read all partitions in parallel, accumulate in shared buffer.
        let shared_buffers = Arc::new(SharedWindowBuffers::new());
        let semaphore = self
            .parallelism
            .map(|n| Arc::new(tokio::sync::Semaphore::new(n.min(partition_count))));

        let mut handles = Vec::with_capacity(partition_count);
        let window_size = self.window_size;

        for partition_id in 0..partition_count {
            let plan = Arc::clone(&plan);
            let task_ctx = Arc::clone(&task_ctx);
            let buffers = Arc::clone(&shared_buffers);
            let sem = semaphore.clone();

            let handle = tokio::spawn(async move {
                if let Some(ref sem) = sem {
                    let _permit = sem
                        .acquire()
                        .await
                        .map_err(|e| DataFusionError::External(Box::new(e)))?;
                }

                read_partition(plan, task_ctx, partition_id, window_size, buffers).await
            });
            handles.push(handle);
        }

        // Collect variant counts from all partitions.
        let mut total_variants = 0u64;
        for handle in handles {
            total_variants += handle
                .await
                .map_err(|e| DataFusionError::External(Box::new(e)))??;
        }

        // Phase 2: Flush all accumulated windows to fjall (sequential).
        // Each window is written exactly once — no read-modify-write.
        let all_windows = shared_buffers.drain();
        let mut total_windows = 0u64;
        let mut total_bytes = 0u64;
        let format_version = self.format_version;

        info!("Flushing {} windows to fjall", all_windows.len());

        for ((chrom, wid), batches) in &all_windows {
            flush_to_store(
                &store,
                chrom,
                *wid,
                batches,
                &mut total_windows,
                &mut total_bytes,
                format_version,
            )?;
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

/// Read all data from a single physical partition and accumulate into shared buffers.
async fn read_partition(
    plan: Arc<dyn ExecutionPlan>,
    task_ctx: Arc<datafusion::execution::TaskContext>,
    partition_id: usize,
    window_size: u64,
    buffers: Arc<SharedWindowBuffers>,
) -> Result<u64> {
    let mut stream = plan.execute(partition_id, task_ctx)?;
    let mut total_variants = 0u64;

    while let Some(batch_result) = stream.next().await {
        let batch = batch_result?;
        if batch.num_rows() == 0 {
            continue;
        }

        total_variants += batch.num_rows() as u64;

        let window_batches = split_batch_into_windows(&batch, window_size)?;
        for ((chrom, wid), sub_batch) in window_batches {
            buffers.push(chrom, wid, sub_batch);
        }
    }

    Ok(total_variants)
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

/// Merge batches for a window and write to the KV store.
fn flush_to_store(
    store: &VepKvStore,
    chrom: &str,
    window_id: u64,
    batches: &[RecordBatch],
    total_windows: &mut u64,
    total_bytes: &mut u64,
    format_version: u8,
) -> Result<()> {
    if batches.is_empty() {
        return Ok(());
    }

    let merged = if batches.len() == 1 {
        batches[0].clone()
    } else {
        let schema = batches[0].schema();
        datafusion::arrow::compute::concat_batches(&schema, batches)
            .map_err(|e| DataFusionError::ArrowError(Box::new(e), None))?
    };

    if format_version >= FORMAT_V1 {
        flush_to_store_v1(store, chrom, window_id, merged, total_windows, total_bytes)?;
    } else {
        flush_to_store_v0(store, chrom, window_id, merged, total_windows, total_bytes)?;
    }

    Ok(())
}

/// v0 (legacy monolithic): write all columns as a single Arrow IPC blob.
fn flush_to_store_v0(
    store: &VepKvStore,
    chrom: &str,
    window_id: u64,
    merged: RecordBatch,
    total_windows: &mut u64,
    total_bytes: &mut u64,
) -> Result<()> {
    let mem_size = merged.get_array_memory_size() as u64;

    // Window should not already exist (all data accumulated before writing).
    // Keep the merge path for safety in case of cross-window indel edge cases.
    if let Some(existing) = store.get_window(chrom, window_id)? {
        let schema = existing.schema();
        let combined = datafusion::arrow::compute::concat_batches(&schema, &[existing, merged])
            .map_err(|e| DataFusionError::ArrowError(Box::new(e), None))?;
        store.put_window(chrom, window_id, &combined)?;
    } else {
        store.put_window(chrom, window_id, &merged)?;
        *total_windows += 1;
    }

    *total_bytes += mem_size;
    Ok(())
}

/// v1 (columnar): write position index + per-column Arrow IPC entries.
fn flush_to_store_v1(
    store: &VepKvStore,
    chrom: &str,
    window_id: u64,
    merged: RecordBatch,
    total_windows: &mut u64,
    total_bytes: &mut u64,
) -> Result<()> {
    // Window should not already exist (all data accumulated before writing).
    // Keep the merge path for safety in case of cross-window indel edge cases.
    let merged = if let Some(existing_idx) = store.get_position_index(chrom, window_id)? {
        let all_col_indices: Vec<u8> = (0..store.schema().fields().len())
            .map(to_v1_column_index)
            .collect::<Result<Vec<_>>>()?;
        if let Some(existing_cols) = store.get_columns(chrom, window_id, &all_col_indices)? {
            let existing_arrays: Vec<datafusion::arrow::array::ArrayRef> =
                existing_cols.into_iter().map(|(a, _)| a).collect();
            let existing_batch = RecordBatch::try_new(store.schema().clone(), existing_arrays)
                .map_err(|e| DataFusionError::ArrowError(Box::new(e), None))?;
            drop(existing_idx);
            datafusion::arrow::compute::concat_batches(
                &store.schema().clone(),
                &[existing_batch, merged],
            )
            .map_err(|e| DataFusionError::ArrowError(Box::new(e), None))?
        } else {
            merged
        }
    } else {
        *total_windows += 1;
        merged
    };

    // Build and write position index.
    let pos_index = PositionIndex::from_batch(&merged)?;
    store.put_position_index(chrom, window_id, &pos_index)?;

    // Write each column separately.
    for (col_idx, field) in merged.schema().fields().iter().enumerate() {
        store.put_column(
            chrom,
            window_id,
            to_v1_column_index(col_idx)?,
            merged.column(col_idx),
            field,
        )?;
    }

    *total_bytes += merged.get_array_memory_size() as u64;
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::key_encoding::EntryType;
    use crate::kv_store::FORMAT_V0;
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
    async fn test_loader_basic_v0() {
        let dir = tempfile::tempdir().unwrap();
        let ctx = SessionContext::new();

        let (schema, table) = make_test_table();
        ctx.register_table("source", Arc::new(table)).unwrap();

        let loader = CacheLoader::new("source", dir.path().to_str().unwrap())
            .with_window_size(1_000_000)
            .with_format_version(FORMAT_V0);
        let stats = loader.load(&ctx).await.unwrap();

        assert_eq!(stats.total_windows, 3);

        let store = VepKvStore::open(dir.path()).unwrap();
        assert_eq!(store.format_version(), FORMAT_V0);
        assert_eq!(store.schema().as_ref(), schema.as_ref());

        let batch = store.get_window("1", 0).unwrap().unwrap();
        assert_eq!(batch.num_rows(), 2);

        let batch = store.get_window("1", 1).unwrap().unwrap();
        assert_eq!(batch.num_rows(), 1);

        let batch = store.get_window("2", 0).unwrap().unwrap();
        assert_eq!(batch.num_rows(), 1);
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
        assert_eq!(store.format_version(), FORMAT_V1);
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
}
