//! Parquet -> fjall cache ingestion pipeline.
//!
//! Loads a VEP variation cache from Parquet into a fjall KV store,
//! grouping variants into windows and parallelizing by chromosome.

use std::path::Path;
use std::sync::Arc;
use std::time::Instant;

use datafusion::arrow::array::{Array, AsArray, RecordBatch};
use datafusion::arrow::datatypes::{Int64Type, SchemaRef};
use datafusion::common::{DataFusionError, Result};
use datafusion::prelude::SessionContext;
use futures::StreamExt;
use log::info;

use crate::key_encoding::{DEFAULT_WINDOW_SIZE, window_id_for_position};
use crate::kv_store::VepKvStore;

/// Statistics returned after loading.
#[derive(Debug, Clone)]
pub struct LoadStats {
    pub total_variants: u64,
    pub total_windows: u64,
    pub total_bytes: u64,
    pub elapsed_secs: f64,
}

/// Builder for loading a Parquet-based VEP cache into fjall.
pub struct CacheLoader {
    source_table: String,
    target_path: String,
    window_size: u64,
    parallelism: usize,
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
            parallelism: 25,
        }
    }

    /// Set the window size (default: 1 Mb).
    pub fn with_window_size(mut self, size: u64) -> Self {
        self.window_size = size;
        self
    }

    /// Set parallelism (default: 25, roughly one per chromosome).
    pub fn with_parallelism(mut self, n: usize) -> Self {
        self.parallelism = n;
        self
    }

    /// Load the cache from the registered Parquet table into fjall.
    pub async fn load(&self, session: &SessionContext) -> Result<LoadStats> {
        let start_time = Instant::now();

        // Get the source schema.
        let source_df = session.table(&self.source_table).await?;
        let schema: SchemaRef = Arc::new(source_df.schema().as_arrow().clone());

        // Create the KV store.
        let store = Arc::new(VepKvStore::create(
            Path::new(&self.target_path),
            schema.clone(),
            self.window_size,
        )?);

        // Discover distinct chromosomes.
        let chrom_df = session
            .sql(&format!(
                "SELECT DISTINCT chrom FROM `{}` ORDER BY chrom",
                self.source_table
            ))
            .await?;
        let chrom_batches = chrom_df.collect().await?;

        let mut chromosomes = Vec::new();
        for batch in &chrom_batches {
            let col = batch.column(0);
            if let Some(arr) = col
                .as_any()
                .downcast_ref::<datafusion::arrow::array::StringArray>()
            {
                for i in 0..arr.len() {
                    if !arr.is_null(i) {
                        chromosomes.push(arr.value(i).to_string());
                    }
                }
            } else if let Some(arr) = col
                .as_any()
                .downcast_ref::<datafusion::arrow::array::StringViewArray>()
            {
                for i in 0..arr.len() {
                    if !arr.is_null(i) {
                        chromosomes.push(arr.value(i).to_string());
                    }
                }
            }
        }

        info!(
            "Loading {} chromosomes into KV cache (window_size={})",
            chromosomes.len(),
            self.window_size
        );

        // Process chromosomes with bounded parallelism.
        let semaphore = Arc::new(tokio::sync::Semaphore::new(self.parallelism));
        let mut handles = Vec::new();
        let window_size = self.window_size;

        for chrom in chromosomes {
            let session = session.clone();
            let source_table = self.source_table.clone();
            let store = Arc::clone(&store);
            let sem = Arc::clone(&semaphore);

            let handle = tokio::spawn(async move {
                let _permit = sem
                    .acquire()
                    .await
                    .map_err(|e| DataFusionError::External(Box::new(e)))?;

                load_chromosome(&session, &source_table, &chrom, &store, window_size).await
            });
            handles.push(handle);
        }

        // Collect results.
        let mut total_variants = 0u64;
        let mut total_windows = 0u64;
        let mut total_bytes = 0u64;

        for handle in handles {
            let (variants, windows, bytes) = handle
                .await
                .map_err(|e| DataFusionError::External(Box::new(e)))??;
            total_variants += variants;
            total_windows += windows;
            total_bytes += bytes;
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

/// Load all variants for a single chromosome.
async fn load_chromosome(
    session: &SessionContext,
    source_table: &str,
    chrom: &str,
    store: &VepKvStore,
    window_size: u64,
) -> Result<(u64, u64, u64)> {
    let query = format!("SELECT * FROM `{source_table}` WHERE chrom = '{chrom}' ORDER BY start");
    let df = session.sql(&query).await?;
    let mut stream = df.execute_stream().await?;

    let mut total_variants = 0u64;
    let mut total_windows = 0u64;
    let mut total_bytes = 0u64;

    // Accumulate rows per window, flush when window changes.
    let mut current_window_id: Option<u64> = None;
    let mut window_batches: Vec<RecordBatch> = Vec::new();

    while let Some(batch_result) = stream.next().await {
        let batch = batch_result?;
        if batch.num_rows() == 0 {
            continue;
        }

        let schema = batch.schema();
        let start_col_idx = schema.index_of("start")?;
        let end_col_idx = schema.index_of("end")?;

        let starts = batch.column(start_col_idx).as_primitive::<Int64Type>();
        let ends = batch.column(end_col_idx).as_primitive::<Int64Type>();

        // Split this batch into sub-batches per window.
        let mut range_start = 0;
        for row in 0..batch.num_rows() {
            let start = starts.value(row);
            let end = ends.value(row);
            let wid = window_id_for_position(start, window_size);

            // Handle cross-window indels: store in all overlapping windows.
            let end_wid = window_id_for_position(end, window_size);
            if end_wid > wid {
                // This variant crosses a window boundary.
                // We'll store a copy in each extra window it touches.
                for extra_wid in (wid + 1)..=end_wid {
                    let single = batch.slice(row, 1);
                    flush_to_store(
                        store,
                        chrom,
                        extra_wid,
                        &[single],
                        &mut total_windows,
                        &mut total_bytes,
                    )?;
                }
            }

            if let Some(cur_wid) = current_window_id {
                if wid != cur_wid {
                    // Flush accumulated rows to the current window.
                    if range_start < row {
                        window_batches.push(batch.slice(range_start, row - range_start));
                    }
                    flush_to_store(
                        store,
                        chrom,
                        cur_wid,
                        &window_batches,
                        &mut total_windows,
                        &mut total_bytes,
                    )?;
                    window_batches.clear();
                    current_window_id = Some(wid);
                    range_start = row;
                }
            } else {
                current_window_id = Some(wid);
            }
        }

        // Remaining rows in this batch belong to current window.
        if range_start < batch.num_rows() {
            window_batches.push(batch.slice(range_start, batch.num_rows() - range_start));
            total_variants += (batch.num_rows() - range_start) as u64;
        }
        // Count all rows processed (including cross-window dupes)
        total_variants += range_start as u64;
    }

    // Flush final window.
    if let Some(wid) = current_window_id {
        if !window_batches.is_empty() {
            flush_to_store(
                store,
                chrom,
                wid,
                &window_batches,
                &mut total_windows,
                &mut total_bytes,
            )?;
        }
    }

    info!("Chromosome {chrom}: {total_variants} variants, {total_windows} windows");

    Ok((total_variants, total_windows, total_bytes))
}

/// Merge batches for a window and write to the KV store.
fn flush_to_store(
    store: &VepKvStore,
    chrom: &str,
    window_id: u64,
    batches: &[RecordBatch],
    total_windows: &mut u64,
    total_bytes: &mut u64,
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

    // If window already exists (from cross-window indel), merge with existing.
    if let Some(existing) = store.get_window(chrom, window_id)? {
        let schema = existing.schema();
        let combined = datafusion::arrow::compute::concat_batches(&schema, &[existing, merged])
            .map_err(|e| DataFusionError::ArrowError(Box::new(e), None))?;
        store.put_window(chrom, window_id, &combined)?;
    } else {
        store.put_window(chrom, window_id, &merged)?;
        *total_windows += 1;
    }

    *total_bytes += batches
        .iter()
        .map(|b| b.get_array_memory_size() as u64)
        .sum::<u64>();
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;
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

        assert_eq!(stats.total_windows, 3); // chr1:window0, chr1:window1, chr2:window0

        // Verify data
        let store = VepKvStore::open(dir.path()).unwrap();
        assert_eq!(store.schema().as_ref(), schema.as_ref());

        // Chr1 window 0 has 2 variants (pos 100, 500_000)
        let batch = store.get_window("1", 0).unwrap().unwrap();
        assert_eq!(batch.num_rows(), 2);

        // Chr1 window 1 has 1 variant (pos 1_500_000)
        let batch = store.get_window("1", 1).unwrap().unwrap();
        assert_eq!(batch.num_rows(), 1);

        // Chr2 window 0 has 1 variant
        let batch = store.get_window("2", 0).unwrap().unwrap();
        assert_eq!(batch.num_rows(), 1);
    }
}
