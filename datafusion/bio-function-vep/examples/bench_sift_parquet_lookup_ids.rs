//! Benchmark compact parquet SIFT/PolyPhen transcript-id lookups.
//!
//! Usage:
//!   ITERS=5 cargo run --release -p datafusion-bio-function-vep --example bench_sift_parquet_lookup_ids -- \
//!     /path/to/translation_sift_lookup <ids.txt>

use std::fs;
use std::path::PathBuf;
use std::time::Instant;

use datafusion::common::{DataFusionError, Result};
use datafusion_bio_function_vep::kv_cache::{SiftParquetStore, SiftPredictionStore};

fn main() -> Result<()> {
    let mut args = std::env::args().skip(1);
    let parquet_path = args.next().map(PathBuf::from).ok_or_else(usage)?;
    let ids_path = args.next().map(PathBuf::from).ok_or_else(usage)?;
    let ids_text = fs::read_to_string(&ids_path).map_err(|error| {
        DataFusionError::Execution(format!("failed to read '{}': {error}", ids_path.display()))
    })?;
    let ids = ids_text
        .lines()
        .map(str::trim)
        .filter(|line| !line.is_empty())
        .map(ToString::to_string)
        .collect::<Vec<_>>();
    let iters = std::env::var("ITERS")
        .ok()
        .and_then(|value| value.parse::<usize>().ok())
        .unwrap_or(3);

    let open_started = Instant::now();
    let store = SiftParquetStore::open_dir(&parquet_path)?.ok_or_else(|| {
        DataFusionError::Execution(format!(
            "no sift lookup parquet files found in '{}'",
            parquet_path.display()
        ))
    })?;
    let open_s = open_started.elapsed().as_secs_f64();

    println!("ids={} iters={} open_s={:.6}", ids.len(), iters, open_s);
    println!("iter\telapsed_s\tgets_s\tfound\tmissing\tsift_entries\tpolyphen_entries");
    for iter in 0..iters {
        let started = Instant::now();
        let found = store.get_many(&ids)?;
        let elapsed_s = started.elapsed().as_secs_f64();
        let sift_entries = found.values().map(|preds| preds.sift.len()).sum::<usize>();
        let polyphen_entries = found
            .values()
            .map(|preds| preds.polyphen.len())
            .sum::<usize>();
        let missing = ids.len().saturating_sub(found.len());
        println!(
            "{}\t{:.6}\t{:.0}\t{}\t{}\t{}\t{}",
            iter,
            elapsed_s,
            ids.len() as f64 / elapsed_s,
            found.len(),
            missing,
            sift_entries,
            polyphen_entries
        );
    }

    Ok(())
}

fn usage() -> DataFusionError {
    DataFusionError::Execution(
        "usage: bench_sift_parquet_lookup_ids <lookup-parquet-or-dir> <ids.txt>".to_string(),
    )
}
