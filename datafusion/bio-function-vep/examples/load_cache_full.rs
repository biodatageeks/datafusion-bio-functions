use std::path::Path;
use std::time::Instant;

use datafusion::prelude::{SessionConfig, SessionContext};
use datafusion_bio_function_vep::kv_cache::{CacheBackend, CacheLoader};

#[tokio::main]
async fn main() -> datafusion::common::Result<()> {
    env_logger::Builder::from_env(env_logger::Env::default().default_filter_or("info")).init();

    let args: Vec<String> = std::env::args().collect();
    if args.len() < 4 {
        eprintln!(
            "Usage: {} <parquet_path> <cache_output_path> <threads> [fjall|redb]",
            args[0]
        );
        std::process::exit(1);
    }

    let parquet_path = &args[1];
    let output_path = &args[2];
    let threads: usize = args[3].parse().map_err(|e| {
        datafusion::common::DataFusionError::Execution(format!("invalid threads: {e}"))
    })?;
    let backend = args
        .get(4)
        .map(|value| CacheBackend::parse(value))
        .transpose()?
        .unwrap_or(CacheBackend::Fjall);

    let config = SessionConfig::new().with_target_partitions(threads);
    let ctx = SessionContext::new_with_config(config);
    ctx.register_parquet("vep_source", parquet_path, Default::default())
        .await?;

    if Path::new(output_path).exists() {
        remove_existing(output_path);
    }

    let start = Instant::now();
    let loader = CacheLoader::new("vep_source", output_path)
        .with_parallelism(threads)
        .with_backend(backend);
    let stats = loader.load(&ctx).await?;

    println!(
        "Loaded: backend={} variants={} positions={} bytes={} elapsed={:.1}s threads={}",
        backend.as_str(),
        stats.total_variants,
        stats.total_positions,
        stats.total_bytes,
        start.elapsed().as_secs_f64(),
        threads
    );

    Ok(())
}

fn remove_existing(path: &str) {
    let path = Path::new(path);
    if path.is_dir() {
        std::fs::remove_dir_all(path).ok();
    } else {
        std::fs::remove_file(path).ok();
    }
}
