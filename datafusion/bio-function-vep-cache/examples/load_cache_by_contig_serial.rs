use std::path::Path;
use std::time::Instant;

use datafusion::arrow::array::{Array, Int64Array, StringArray, StringViewArray};
use datafusion::common::{DataFusionError, Result};
use datafusion::prelude::{SessionConfig, SessionContext};
use datafusion_bio_function_vep_cache::CacheLoader;

#[derive(Debug)]
struct ContigPlan {
    chrom: String,
    variants: i64,
}

#[tokio::main]
async fn main() -> Result<()> {
    env_logger::Builder::from_env(env_logger::Env::default().default_filter_or("info")).init();

    let args: Vec<String> = std::env::args().collect();
    if args.len() < 3 {
        eprintln!(
            "Usage: {} <parquet_path> <fjall_output_path> [window_size_kb]",
            args[0]
        );
        std::process::exit(1);
    }

    let parquet_path = &args[1];
    let output_path = &args[2];
    let window_size_kb: u64 = args.get(3).and_then(|s| s.parse().ok()).unwrap_or(1000);
    let window_size = window_size_kb * 1000;

    let config = SessionConfig::new().with_target_partitions(1);
    let ctx = SessionContext::new_with_config(config);
    ctx.register_parquet("vep_source", parquet_path, Default::default())
        .await?;

    let contigs = load_contigs(&ctx).await?;
    if contigs.is_empty() {
        return Err(DataFusionError::Execution(
            "No contigs found in source parquet".to_string(),
        ));
    }

    if Path::new(output_path).exists() {
        std::fs::remove_dir_all(output_path)
            .map_err(|e| DataFusionError::Execution(format!("failed to clean output path: {e}")))?;
    }

    println!(
        "contigs={} window_kb={} target_partitions=1 loader_parallelism=1",
        contigs.len(),
        window_size_kb
    );
    println!(
        "{:<8} {:>14} {:>10} {:>10} {:>9} {:>13} {:>13}",
        "chrom", "variants", "windows", "bytes", "elapsed_s", "delta_gb", "total_gb"
    );

    let mut prev_total_bytes = 0u64;
    let run_start = Instant::now();
    for (idx, plan) in contigs.iter().enumerate() {
        let escaped = plan.chrom.replace('\'', "''");
        let view_sql = format!(
            "CREATE OR REPLACE VIEW vep_contig AS SELECT * FROM vep_source WHERE chrom = '{escaped}'"
        );
        ctx.sql(&view_sql).await?;

        let start = Instant::now();
        let stats = CacheLoader::new("vep_contig", output_path)
            .with_window_size(window_size)
            .with_parallelism(1)
            .load(&ctx)
            .await?;
        let elapsed = start.elapsed().as_secs_f64();

        let total_bytes = dir_size_bytes(Path::new(output_path));
        let delta_bytes = total_bytes.saturating_sub(prev_total_bytes);
        prev_total_bytes = total_bytes;

        println!(
            "{:<8} {:>14} {:>10} {:>10} {:>9.1} {:>13.3} {:>13.3}",
            plan.chrom,
            plan.variants,
            stats.total_windows,
            stats.total_bytes,
            elapsed,
            bytes_to_gb(delta_bytes),
            bytes_to_gb(total_bytes)
        );

        if idx + 1 == contigs.len() {
            println!(
                "done elapsed_total_s={:.1}",
                run_start.elapsed().as_secs_f64()
            );
        }
    }

    Ok(())
}

async fn load_contigs(ctx: &SessionContext) -> Result<Vec<ContigPlan>> {
    let df = ctx
        .sql(
            "SELECT chrom, COUNT(*)::BIGINT AS variants \
             FROM vep_source \
             WHERE chrom IS NOT NULL \
             GROUP BY chrom \
             ORDER BY chrom",
        )
        .await?;
    let batches = df.collect().await?;

    let mut out = Vec::new();
    for batch in batches {
        let chrom_col = batch.column(0);
        let variants = batch
            .column(1)
            .as_any()
            .downcast_ref::<Int64Array>()
            .ok_or_else(|| DataFusionError::Execution("variants must be Int64".to_string()))?;

        for row in 0..batch.num_rows() {
            let chrom = string_value(chrom_col.as_ref(), row).to_string();
            out.push(ContigPlan {
                chrom,
                variants: variants.value(row),
            });
        }
    }
    Ok(out)
}

fn string_value(col: &dyn Array, row: usize) -> &str {
    if let Some(arr) = col.as_any().downcast_ref::<StringArray>() {
        arr.value(row)
    } else if let Some(arr) = col.as_any().downcast_ref::<StringViewArray>() {
        arr.value(row)
    } else {
        panic!("Expected Utf8 or Utf8View column for chrom")
    }
}

fn dir_size_bytes(path: &Path) -> u64 {
    let mut total = 0u64;
    if let Ok(entries) = std::fs::read_dir(path) {
        for entry in entries.flatten() {
            if let Ok(meta) = entry.metadata() {
                if meta.is_file() {
                    total = total.saturating_add(meta.len());
                } else if meta.is_dir() {
                    total = total.saturating_add(dir_size_bytes(&entry.path()));
                }
            }
        }
    }
    total
}

fn bytes_to_gb(bytes: u64) -> f64 {
    bytes as f64 / (1024.0 * 1024.0 * 1024.0)
}
