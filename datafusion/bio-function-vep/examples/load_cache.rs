use std::path::Path;
use std::time::Instant;

use datafusion::prelude::{SessionConfig, SessionContext};
use datafusion_bio_function_vep::kv_cache::{CacheLoader, VepKvStore};

#[tokio::main]
async fn main() -> datafusion::common::Result<()> {
    env_logger::Builder::from_env(env_logger::Env::default().default_filter_or("warn")).init();

    let args: Vec<String> = std::env::args().collect();
    if args.len() < 3 {
        eprintln!(
            "Usage: {} <parquet_path> <fjall_output_path> [chrom_filter] [partitions] [zstd_level] [dict_size_kb]",
            args[0]
        );
        eprintln!(
            "Example: {} /data/vep/115_GRCh38.parquet /tmp/vep_cache 22 8 9 256",
            args[0]
        );
        std::process::exit(1);
    }

    let parquet_path = &args[1];
    let output_path = &args[2];
    let chrom_filter = args.get(3).map(|s| s.as_str());
    let target_partitions: Option<usize> = args.get(4).and_then(|s| s.parse().ok());
    let zstd_level: i32 = args.get(5).and_then(|s| s.parse().ok()).unwrap_or(3);
    let dict_size_kb: u32 = args.get(6).and_then(|s| s.parse().ok()).unwrap_or(112);

    let ctx = if let Some(partitions) = target_partitions {
        let config = SessionConfig::new().with_target_partitions(partitions);
        SessionContext::new_with_config(config)
    } else {
        SessionContext::new()
    };

    // Register parquet source
    ctx.register_parquet("vep_source", parquet_path, Default::default())
        .await?;

    // If a chromosome filter is provided, create a filtered view
    let source_table = if let Some(chrom) = chrom_filter {
        ctx.sql(&format!(
            "CREATE VIEW vep_filtered AS SELECT * FROM vep_source WHERE chrom = '{chrom}'"
        ))
        .await?;
        "vep_filtered"
    } else {
        "vep_source"
    };

    // Count rows
    let count_df = ctx
        .sql(&format!("SELECT COUNT(*) as cnt FROM {source_table}"))
        .await?;
    let count_batches = count_df.collect().await?;
    let row_count = count_batches[0]
        .column(0)
        .as_any()
        .downcast_ref::<datafusion::arrow::array::Int64Array>()
        .unwrap()
        .value(0);

    // Clean output dir if exists
    if Path::new(output_path).exists() {
        std::fs::remove_dir_all(output_path).ok();
    }

    // Load into fjall
    let load_start = Instant::now();
    let loader = CacheLoader::new(source_table, output_path)
        .with_zstd_level(zstd_level)
        .with_dict_size_kb(dict_size_kb);
    let stats = loader.load(&ctx).await?;
    let load_elapsed = load_start.elapsed().as_secs_f64();

    // Verify + measure disk
    let store = VepKvStore::open(output_path)?;
    let dir_size = walkdir(output_path);
    let avg_variants_per_pos = if stats.total_positions > 0 {
        row_count as f64 / stats.total_positions as f64
    } else {
        0.0
    };
    let avg_pos_disk_kb = if stats.total_positions > 0 {
        dir_size as f64 / stats.total_positions as f64 / 1024.0
    } else {
        0.0
    };

    // Measure a few position lookups
    let sample_chrom = chrom_filter.unwrap_or("22");
    let chrom_code =
        datafusion_bio_function_vep::kv_cache::key_encoding::chrom_to_code(sample_chrom);
    let read_start = Instant::now();
    let mut read_iters = 0u32;
    for pos in 0..20u64 {
        if store
            .get_position_entry_decompressed(chrom_code, pos as i64)?
            .is_some()
        {
            read_iters += 1;
        }
    }
    let avg_read_ms = if read_iters > 0 {
        read_start.elapsed().as_secs_f64() * 1000.0 / read_iters as f64
    } else {
        0.0
    };

    // Single lookup from mid-range
    let mid_pos = stats.total_positions / 2;
    let single_start = Instant::now();
    let _ = store.get_position_entry_decompressed(chrom_code, mid_pos as i64)?;
    let single_read_ms = single_start.elapsed().as_secs_f64() * 1000.0;

    drop(store);

    println!(
        "zstd_level={} dict_kb={} | variants={:>10} | positions={:>6} | avg/pos={:>8.0} | load={:>5.1}s | disk={:>7.1}MB | avg_pos_disk={:>6.0}KB | read={:>6.1}ms avg | single={:>6.1}ms",
        zstd_level,
        dict_size_kb,
        row_count,
        stats.total_positions,
        avg_variants_per_pos,
        load_elapsed,
        dir_size as f64 / 1_048_576.0,
        avg_pos_disk_kb,
        avg_read_ms,
        single_read_ms,
    );

    Ok(())
}

fn walkdir(path: &str) -> u64 {
    let mut total = 0u64;
    if let Ok(entries) = std::fs::read_dir(path) {
        for entry in entries.flatten() {
            let meta = entry.metadata();
            if let Ok(m) = meta {
                if m.is_file() {
                    total += m.len();
                } else if m.is_dir() {
                    total += walkdir(&entry.path().to_string_lossy());
                }
            }
        }
    }
    total
}
