//! Rebuild the port-test Parquet fixture cache (`vep-benchmark/data/port/_cache115`)
//! from a raw Ensembl VEP cache.
//!
//! The port tests annotate against a committed (Git-LFS) Parquet cache. That cache
//! must be produced by the **current** builder: the read path needs a per-entity
//! `chrom_manifest.json` *and* a Parquet page index (`PageDir`), neither of which
//! pre-v0.14.0 caches carry. Regenerating it by hand is not possible — the variation
//! shard is written as two `start`-ordered runs (warm tier then cold tier) that the
//! read-side binary search relies on.
//!
//! To keep the fixture small enough for Git-LFS and CI, point `--cache-root` at a
//! *trimmed* native cache: the 1 Mb chunk files covering the loci the port tests
//! query, plus an `all_vars.gz` sliced to the same windows with `tabix`. See
//! `porting-tests/` for the trimming recipe.
//!
//! Usage:
//!   cargo run --release -p datafusion-bio-function-vep \
//!     --features cache-builder \
//!     --example rebuild_port_fixture -- \
//!     --cache-root /path/to/trimmed/homo_sapiens/115_GRCh38 \
//!     --output-dir vep-benchmark/data/port/_cache115/parquet/115_GRCh38_vep \
//!     --chroms 21,MT --overwrite

use std::path::PathBuf;
use std::time::Instant;

use datafusion::common::{DataFusionError, Result};
use datafusion_bio_function_vep::cache_builder::CacheBuilder;

#[tokio::main]
async fn main() -> Result<()> {
    let _ = env_logger::try_init();
    let args = parse_args()?;
    let started = Instant::now();

    eprintln!("cache_root={}", args.cache_root.display());
    eprintln!("output_dir={}", args.output_dir.display());
    eprintln!("chroms={:?}", args.chroms);

    let stats = CacheBuilder::new(
        args.cache_root.to_string_lossy().to_string(),
        args.output_dir.to_string_lossy().to_string(),
    )
    .with_chrom_filter(args.chroms.clone())
    .with_partitions(args.partitions)
    .with_overwrite(args.overwrite)
    .build_all()
    .await?;

    let mut total = 0usize;
    for entity in &stats {
        let rows: usize = entity.parquet_files.iter().map(|(_, rows)| rows).sum();
        total += rows;
        eprintln!("  {:<18} {rows:>9} rows", entity.entity);
    }
    eprintln!(
        "rebuilt {} entities, {total} rows, elapsed={:.1}s",
        stats.len(),
        started.elapsed().as_secs_f64()
    );
    Ok(())
}

#[derive(Debug)]
struct Args {
    cache_root: PathBuf,
    output_dir: PathBuf,
    chroms: Vec<String>,
    partitions: usize,
    overwrite: bool,
}

fn parse_args() -> Result<Args> {
    let mut cache_root = None;
    let mut output_dir = None;
    let mut chroms: Vec<String> = Vec::new();
    let mut partitions = 4usize;
    let mut overwrite = false;

    let mut args = std::env::args().skip(1);
    while let Some(arg) = args.next() {
        match arg.as_str() {
            "--cache-root" => cache_root = Some(PathBuf::from(require_value(&mut args, &arg)?)),
            "--output-dir" => output_dir = Some(PathBuf::from(require_value(&mut args, &arg)?)),
            "--chroms" => {
                chroms = require_value(&mut args, &arg)?
                    .split(',')
                    .map(|chrom| chrom.trim().to_string())
                    .filter(|chrom| !chrom.is_empty())
                    .collect();
            }
            "--partitions" => {
                partitions = require_value(&mut args, &arg)?.parse::<usize>().map_err(|err| {
                    DataFusionError::Execution(format!("invalid --partitions: {err}"))
                })?;
            }
            "--overwrite" => overwrite = true,
            other => {
                return Err(DataFusionError::Execution(format!(
                    "unknown argument: {other}"
                )));
            }
        }
    }

    Ok(Args {
        cache_root: cache_root
            .ok_or_else(|| DataFusionError::Execution("--cache-root is required".into()))?,
        output_dir: output_dir
            .ok_or_else(|| DataFusionError::Execution("--output-dir is required".into()))?,
        chroms,
        partitions: partitions.max(1),
        overwrite,
    })
}

fn require_value(args: &mut impl Iterator<Item = String>, flag: &str) -> Result<String> {
    args.next()
        .ok_or_else(|| DataFusionError::Execution(format!("{flag} requires a value")))
}
