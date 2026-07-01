//! Build the FULL single-path Lance cache (all entities, all contigs) from a raw
//! Ensembl/merged/refseq VEP cache. Wraps `CacheBuilder::build_all`.
//!
//! Usage:
//!   cargo run --release -p datafusion-bio-function-vep \
//!     --features lance-cache,cache-builder \
//!     --example build_lance_full_cache -- \
//!     --cache-root /path/to/homo_sapiens_merged/115_GRCh38 \
//!     --output-dir /path/to/115_GRCh38_merged \
//!     --cache-source-type merged --partitions 8 --overwrite

use std::str::FromStr;
use std::time::Instant;

use datafusion::common::{DataFusionError, Result};
use datafusion_bio_format_ensembl_cache::CacheSourceType;
use datafusion_bio_function_vep::cache_builder::{CacheBuilder, CacheFormat};

#[tokio::main]
async fn main() -> Result<()> {
    let mut cache_root = None;
    let mut output_dir = None;
    let mut cache_source_type = CacheSourceType::Ensembl;
    let mut partitions = 8usize;
    let mut overwrite = false;

    let mut args = std::env::args().skip(1);
    while let Some(arg) = args.next() {
        let mut value = || {
            args.next()
                .ok_or_else(|| DataFusionError::Execution(format!("{arg} requires a value")))
        };
        match arg.as_str() {
            "--cache-root" => cache_root = Some(value()?),
            "--output-dir" => output_dir = Some(value()?),
            "--cache-source-type" => {
                let v = value()?;
                cache_source_type = CacheSourceType::from_str(&v).map_err(|_| {
                    DataFusionError::Execution(format!(
                        "invalid --cache-source-type '{v}', expected ensembl, merged, or refseq"
                    ))
                })?;
            }
            "--partitions" => {
                partitions = value()?
                    .parse()
                    .map_err(|_| DataFusionError::Execution("invalid --partitions".into()))?;
            }
            "--overwrite" => overwrite = true,
            other => {
                return Err(DataFusionError::Execution(format!(
                    "unknown argument: {other}"
                )));
            }
        }
    }

    let cache_root =
        cache_root.ok_or_else(|| DataFusionError::Execution("--cache-root is required".into()))?;
    let output_dir =
        output_dir.ok_or_else(|| DataFusionError::Execution("--output-dir is required".into()))?;

    eprintln!("cache_root={cache_root}");
    eprintln!("output_dir={output_dir}");
    eprintln!("cache_source_type={cache_source_type}");
    eprintln!("partitions={partitions} overwrite={overwrite}");

    let started = Instant::now();
    let builder = CacheBuilder::new(cache_root, output_dir)
        .with_cache_format(CacheFormat::Lance)
        .with_cache_source_type(cache_source_type)
        .with_partitions(partitions.max(1))
        .with_overwrite(overwrite);

    let stats = builder.build_all().await?;
    for stat in &stats {
        let rows: usize = stat.parquet_files.iter().map(|(_, rows)| *rows).sum();
        eprintln!(
            "built entity={} files={} rows={}",
            stat.entity,
            stat.parquet_files.len(),
            rows
        );
    }
    eprintln!(
        "FULL BUILD DONE entities={} elapsed={:.1}s",
        stats.len(),
        started.elapsed().as_secs_f64()
    );
    Ok(())
}
