#![allow(
    clippy::collapsible_if,
    clippy::too_many_arguments,
    unused_variables,
    clippy::manual_saturating_arithmetic,
    clippy::implicit_saturating_sub
)]
//! Benchmark: individual SIFT/PolyPhen point queries on translation parquet.
//!
//! Simulates the lazy per-transcript loading approach where each transcript's
//! predictions are fetched on demand with a single-row SQL query.
//!
//! Usage:
//!   cargo run --release -p datafusion-bio-function-vep --example bench_sift_queries \
//!     /path/to/translation.parquet

use std::time::Instant;

use datafusion::arrow::array::Array;
use datafusion::common::Result;
use datafusion::prelude::{ParquetReadOptions, SessionConfig, SessionContext};

#[tokio::main]
async fn main() -> Result<()> {
    let args: Vec<String> = std::env::args().collect();
    let parquet_path = args
        .get(1)
        .expect("usage: bench_sift_queries <translation.parquet>");

    let config = SessionConfig::new().with_target_partitions(1);
    let ctx = SessionContext::new_with_config(config);
    ctx.register_parquet("translations", parquet_path, ParquetReadOptions::default())
        .await?;

    // Collect all transcript_ids from the table.
    let all_ids: Vec<String> = {
        let batches = ctx
            .sql("SELECT DISTINCT transcript_id FROM translations")
            .await?
            .collect()
            .await?;
        let mut ids = Vec::new();
        for batch in &batches {
            let col = batch
                .column_by_name("transcript_id")
                .expect("missing transcript_id");
            for row in 0..batch.num_rows() {
                if let Some(arr) = col
                    .as_any()
                    .downcast_ref::<datafusion::arrow::array::StringArray>()
                {
                    if !arr.is_null(row) {
                        ids.push(arr.value(row).to_string());
                    }
                } else if let Some(arr) = col
                    .as_any()
                    .downcast_ref::<datafusion::arrow::array::StringViewArray>()
                {
                    if !arr.is_null(row) {
                        ids.push(arr.value(row).to_string());
                    }
                }
            }
        }
        ids
    };
    println!("total transcripts: {}", all_ids.len());

    // Deterministic sample: take every Nth element.
    for &count in &[500usize, 1000, 5000] {
        let count = count.min(all_ids.len());
        let step = all_ids.len() / count;
        let sample: Vec<&str> = all_ids
            .iter()
            .step_by(step)
            .take(count)
            .map(|s| s.as_str())
            .collect();
        println!("\n=== {} individual point queries ===", sample.len());

        let start = Instant::now();
        let mut total_rows = 0u64;
        let mut peak_batch_bytes = 0usize;
        for tx_id in &sample {
            let sql = format!(
                "SELECT sift_predictions, polyphen_predictions \
                 FROM translations WHERE transcript_id = '{tx_id}'"
            );
            let batches = ctx.sql(&sql).await?.collect().await?;
            for b in &batches {
                total_rows += b.num_rows() as u64;
                let batch_bytes: usize =
                    b.columns().iter().map(|c| c.get_buffer_memory_size()).sum();
                if batch_bytes > peak_batch_bytes {
                    peak_batch_bytes = batch_bytes;
                }
            }
        }
        let elapsed = start.elapsed();
        println!(
            "  total: {:.2}s, {:.1}ms/query, {} rows returned, peak batch ~{:.1}MB",
            elapsed.as_secs_f64(),
            elapsed.as_secs_f64() / sample.len() as f64 * 1000.0,
            total_rows,
            peak_batch_bytes as f64 / 1e6,
        );
    }

    Ok(())
}
