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
            "Usage: {} <parquet_path> <fjall_output_path> <threads>",
            args[0]
        );
        std::process::exit(1);
    }

    let parquet_path = &args[1];
    let output_path = &args[2];
    let threads: usize = args[3].parse().map_err(|e| {
        datafusion::common::DataFusionError::Execution(format!("invalid threads: {e}"))
    })?;

    let config = SessionConfig::new().with_target_partitions(threads);
    let ctx = SessionContext::new_with_config(config);
    ctx.register_parquet("vep_source", parquet_path, Default::default())
        .await?;

    if Path::new(output_path).exists() {
        std::fs::remove_dir_all(output_path).ok();
    }

    let start = Instant::now();
    let loader = CacheLoader::new("vep_source", output_path).with_parallelism(threads);
    let stats = loader.load(&ctx).await?;

    println!(
        "Loaded: variants={} positions={} bytes={} elapsed={:.1}s threads={}",
        stats.total_variants,
        stats.total_positions,
        stats.total_bytes,
        start.elapsed().as_secs_f64(),
        threads
    );

    Ok(())
}
