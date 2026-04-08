use std::path::Path;
use std::time::Instant;

use datafusion::prelude::{SessionConfig, SessionContext};
use datafusion_bio_function_vep::kv_cache::CacheLoader;

#[tokio::main]
async fn main() -> datafusion::common::Result<()> {
    env_logger::Builder::from_env(env_logger::Env::default().default_filter_or("info")).init();

    let args: Vec<String> = std::env::args().collect();
    if args.len() < 4 {
        eprintln!(
            "Usage: {} <parquet_path> <fjall_output_path> <threads> [chromosomes_csv]",
            args[0]
        );
        std::process::exit(1);
    }

    let parquet_path = &args[1];
    let output_path = &args[2];
    let threads: usize = args[3].parse().map_err(|e| {
        datafusion::common::DataFusionError::Execution(format!("invalid threads: {e}"))
    })?;
    let chromosomes: Vec<String> = args
        .get(4)
        .map(|value| {
            value
                .split(',')
                .map(str::trim)
                .filter(|chrom| !chrom.is_empty())
                .map(ToOwned::to_owned)
                .collect()
        })
        .unwrap_or_default();

    let config = SessionConfig::new().with_target_partitions(threads);
    let ctx = SessionContext::new_with_config(config);
    ctx.register_parquet("vep_source", parquet_path, Default::default())
        .await?;

    let source_table = if chromosomes.is_empty() {
        "vep_source"
    } else {
        let filter = chromosomes
            .iter()
            .map(|chrom| format!("'{chrom}'"))
            .collect::<Vec<_>>()
            .join(", ");
        ctx.sql(&format!(
            "CREATE VIEW vep_filtered AS SELECT * FROM vep_source WHERE chrom IN ({filter})"
        ))
        .await?;
        "vep_filtered"
    };

    if Path::new(output_path).exists() {
        std::fs::remove_dir_all(output_path).ok();
    }

    let start = Instant::now();
    let loader = CacheLoader::new(source_table, output_path).with_parallelism(threads);
    let stats = loader.load(&ctx).await?;

    println!(
        "Loaded: variants={} positions={} bytes={} elapsed={:.1}s threads={} chromosomes={}",
        stats.total_variants,
        stats.total_positions,
        stats.total_bytes,
        start.elapsed().as_secs_f64(),
        threads,
        if chromosomes.is_empty() {
            "all".to_string()
        } else {
            chromosomes.join(",")
        }
    );

    Ok(())
}
