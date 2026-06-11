//! Build a merged Lance variation cache from existing warm/cold Parquet tiers.
//!
//! Usage:
//!   cargo run --release -p datafusion-bio-function-vep --features lance-cache \
//!     --example build_lance_variation_cache -- \
//!     --cache /path/to/cache --chrom chr1 \
//!     [--warm-fragment-rows 500000] [--cold-fragment-rows 1024] [--batch-size 65536] \
//!     [--overwrite]

use std::fs::File;
use std::path::{Path, PathBuf};
use std::time::Instant;

use datafusion::arrow::record_batch::RecordBatch;
use datafusion::common::{DataFusionError, Result};
use datafusion_bio_function_vep::kv_cache::position_index::cold_variation_files_for_chrom;
use datafusion_bio_function_vep::warm_cache::lance_variation::{
    DEFAULT_WARM_LANCE_ROWS_PER_FILE, lance_variation_dataset_path,
    write_merged_lance_variation_dataset,
};
use datafusion_bio_function_vep::warm_cache::reader::{
    WARM_RUNTIME_COLUMNS, projection_for_existing_roots,
};
use parquet::arrow::arrow_reader::{ArrowReaderMetadata, ParquetRecordBatchReaderBuilder};

#[derive(Debug)]
struct Args {
    cache: PathBuf,
    chrom: String,
    warm_fragment_rows: usize,
    cold_fragment_rows: usize,
    batch_size: usize,
    overwrite: bool,
}

#[tokio::main]
async fn main() -> Result<()> {
    let args = parse_args()?;
    let started = Instant::now();
    let variation_dir = args.cache.join("variation");
    let output = lance_variation_dataset_path(&args.cache, &args.chrom);

    if output.exists() && !args.overwrite {
        return Err(DataFusionError::Execution(format!(
            "Lance variation dataset '{}' already exists; pass --overwrite",
            output.display()
        )));
    }

    let warm_path = warm_file_for_chrom(&variation_dir, &args.chrom).ok_or_else(|| {
        DataFusionError::Execution(format!(
            "missing warm variation parquet for {} under '{}'",
            args.chrom,
            variation_dir.display()
        ))
    })?;
    let cold_paths = cold_variation_files_for_chrom(&variation_dir, &args.chrom);
    if cold_paths.is_empty() {
        return Err(DataFusionError::Execution(format!(
            "missing cold variation parquet for {} under '{}'",
            args.chrom,
            variation_dir.display()
        )));
    }

    eprintln!("warm={}", warm_path.display());
    for path in &cold_paths {
        eprintln!("cold={}", path.display());
    }
    eprintln!("output={}", output.display());
    eprintln!("warm_fragment_rows={}", args.warm_fragment_rows);
    eprintln!("cold_fragment_rows={}", args.cold_fragment_rows);
    eprintln!("batch_size={}", args.batch_size);

    let projection_columns = runtime_projection_columns();
    eprintln!("projected_columns={}", projection_columns.len());

    let warm_batches = read_parquet_batches(&warm_path, args.batch_size, &projection_columns)?;
    let mut cold_batches = Vec::new();
    for path in &cold_paths {
        cold_batches.extend(read_parquet_batches(
            path,
            args.batch_size,
            &projection_columns,
        )?);
    }
    let warm_rows = row_count(&warm_batches);
    let cold_rows = row_count(&cold_batches);

    write_merged_lance_variation_dataset(
        &output,
        warm_batches,
        cold_batches,
        args.warm_fragment_rows,
        args.cold_fragment_rows,
    )
    .await?;

    eprintln!(
        "wrote warm_rows={} cold_rows={} elapsed={:.2}s",
        warm_rows,
        cold_rows,
        started.elapsed().as_secs_f64()
    );
    Ok(())
}

fn parse_args() -> Result<Args> {
    let mut cache = None;
    let mut chrom = None;
    let mut warm_fragment_rows = DEFAULT_WARM_LANCE_ROWS_PER_FILE;
    let mut cold_fragment_rows = 1_024usize;
    let mut batch_size = 65_536usize;
    let mut overwrite = false;

    let mut args = std::env::args().skip(1);
    while let Some(arg) = args.next() {
        match arg.as_str() {
            "--cache" => cache = Some(PathBuf::from(require_value(&mut args, "--cache")?)),
            "--chrom" => chrom = Some(require_value(&mut args, "--chrom")?),
            "--warm-fragment-rows" => {
                warm_fragment_rows = parse_usize(require_value(&mut args, "--warm-fragment-rows")?)?
            }
            "--cold-fragment-rows" => {
                cold_fragment_rows = parse_usize(require_value(&mut args, "--cold-fragment-rows")?)?
            }
            "--batch-size" => batch_size = parse_usize(require_value(&mut args, "--batch-size")?)?,
            "--overwrite" => overwrite = true,
            "--help" | "-h" => {
                print_usage();
                std::process::exit(0);
            }
            other => {
                return Err(DataFusionError::Execution(format!(
                    "unknown argument: {other}"
                )));
            }
        }
    }

    Ok(Args {
        cache: cache.ok_or_else(|| DataFusionError::Execution("--cache is required".into()))?,
        chrom: chrom.ok_or_else(|| DataFusionError::Execution("--chrom is required".into()))?,
        warm_fragment_rows: warm_fragment_rows.max(1),
        cold_fragment_rows: cold_fragment_rows.max(1),
        batch_size: batch_size.max(1),
        overwrite,
    })
}

fn require_value(args: &mut impl Iterator<Item = String>, flag: &str) -> Result<String> {
    args.next()
        .ok_or_else(|| DataFusionError::Execution(format!("{flag} requires a value")))
}

fn parse_usize(value: String) -> Result<usize> {
    value
        .parse::<usize>()
        .map_err(|error| DataFusionError::Execution(format!("invalid integer '{value}': {error}")))
}

fn print_usage() {
    eprintln!(
        "Usage: build_lance_variation_cache --cache /path/to/cache --chrom chr1 \
         [--warm-fragment-rows 500000] [--cold-fragment-rows 1024] [--batch-size 65536] \
         [--overwrite]"
    );
}

fn warm_file_for_chrom(variation_dir: &Path, chrom: &str) -> Option<PathBuf> {
    variation_split_file_for_chrom(variation_dir, chrom, "warm")
}

fn variation_split_file_for_chrom(
    variation_dir: &Path,
    chrom: &str,
    tier: &str,
) -> Option<PathBuf> {
    let direct = variation_dir.join(format!("{chrom}_{tier}.parquet"));
    if direct.is_file() {
        return Some(direct);
    }

    if let Some(stripped) = chrom.strip_prefix("chr") {
        let stripped = variation_dir.join(format!("{stripped}_{tier}.parquet"));
        if stripped.is_file() {
            return Some(stripped);
        }
    } else {
        let prefixed = variation_dir.join(format!("chr{chrom}_{tier}.parquet"));
        if prefixed.is_file() {
            return Some(prefixed);
        }
    }

    None
}

fn runtime_projection_columns() -> Vec<String> {
    WARM_RUNTIME_COLUMNS
        .iter()
        .map(|name| (*name).to_string())
        .collect()
}

fn read_parquet_batches(
    path: &Path,
    batch_size: usize,
    projection_columns: &[String],
) -> Result<Vec<RecordBatch>> {
    let file = File::open(path).map_err(|error| {
        DataFusionError::Execution(format!("failed to open '{}': {error}", path.display()))
    })?;
    let metadata = ArrowReaderMetadata::load(&file, Default::default())?;
    let mask = projection_for_existing_roots(
        metadata.schema(),
        metadata.parquet_schema(),
        projection_columns,
    );
    let reader = ParquetRecordBatchReaderBuilder::new_with_metadata(file, metadata)
        .with_projection(mask)
        .with_batch_size(batch_size.max(1))
        .build()?;

    let mut batches = Vec::new();
    for batch in reader {
        batches.push(batch?);
    }
    Ok(batches)
}

fn row_count(batches: &[RecordBatch]) -> usize {
    batches.iter().map(RecordBatch::num_rows).sum()
}
