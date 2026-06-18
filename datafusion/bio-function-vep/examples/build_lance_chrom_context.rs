//! Rebuild a single chromosome's context entities (transcript, translation_core,
//! translation_sift) of the single-path Lance cache from the Ensembl-cache table
//! provider. Scoped to one chromosome via the build's `chrom_filter`, so it is
//! cheap enough for profiling / parity iteration on the position-sliced SIFT
//! layout without rebuilding the whole genome.
//!
//! Usage:
//!   cargo run --release -p datafusion-bio-function-vep \
//!     --features lance-cache,cache-builder \
//!     --example build_lance_chrom_context -- \
//!     --cache-root /path/to/homo_sapiens_merged/115_GRCh38 \
//!     --output-dir /path/to/scratch_chr1_cache \
//!     --chrom chr1 --cache-source-type merged --partitions 8 --overwrite \
//!     [--entities transcript,translation]

use std::str::FromStr;
use std::time::Instant;

use datafusion::common::{DataFusionError, Result};
use datafusion_bio_format_ensembl_cache::CacheSourceType;
use datafusion_bio_function_vep::cache_builder::{CacheBuilder, CacheFormat};

#[derive(Debug)]
struct Args {
    cache_root: String,
    output_dir: String,
    chrom: String,
    cache_source_type: CacheSourceType,
    partitions: usize,
    overwrite: bool,
    entities: Vec<String>,
}

#[tokio::main]
async fn main() -> Result<()> {
    let args = parse_args()?;
    let started = Instant::now();

    eprintln!("cache_root={}", args.cache_root);
    eprintln!("output_dir={}", args.output_dir);
    eprintln!("chrom={}", args.chrom);
    eprintln!("cache_source_type={}", args.cache_source_type);
    eprintln!("partitions={}", args.partitions);
    eprintln!("overwrite={}", args.overwrite);
    eprintln!("entities={}", args.entities.join(","));

    let builder = CacheBuilder::new(args.cache_root.clone(), args.output_dir.clone())
        .with_cache_format(CacheFormat::Lance)
        .with_cache_source_type(args.cache_source_type)
        .with_partitions(args.partitions)
        .with_overwrite(args.overwrite)
        .with_chrom_filter([args.chrom.clone()]);

    for entity in &args.entities {
        let entity_started = Instant::now();
        let stats = builder.build_entity(entity).await?;
        for stat in &stats {
            let rows: usize = stat.parquet_files.iter().map(|(_, rows)| *rows).sum();
            eprintln!(
                "built entity={} dataset={} files={} rows={} elapsed={:.2}s",
                entity,
                stat.entity,
                stat.parquet_files.len(),
                rows,
                entity_started.elapsed().as_secs_f64()
            );
        }
    }

    eprintln!(
        "done chrom={} entities={} elapsed={:.2}s",
        args.chrom,
        args.entities.join(","),
        started.elapsed().as_secs_f64()
    );
    Ok(())
}

fn parse_args() -> Result<Args> {
    let mut cache_root = None;
    let mut output_dir = None;
    let mut chrom = None;
    let mut cache_source_type = CacheSourceType::Ensembl;
    let mut partitions = 8usize;
    let mut overwrite = false;
    let mut entities: Vec<String> = vec!["transcript".to_string(), "translation".to_string()];

    let mut args = std::env::args().skip(1);
    while let Some(arg) = args.next() {
        match arg.as_str() {
            "--cache-root" => cache_root = Some(require_value(&mut args, &arg)?),
            "--output-dir" => output_dir = Some(require_value(&mut args, &arg)?),
            "--chrom" => chrom = Some(require_value(&mut args, &arg)?),
            "--cache-source-type" => {
                let value = require_value(&mut args, &arg)?;
                cache_source_type = CacheSourceType::from_str(&value).map_err(|_| {
                    DataFusionError::Execution(format!(
                        "invalid --cache-source-type '{value}', expected ensembl, merged, or refseq"
                    ))
                })?;
            }
            "--partitions" => partitions = parse_usize(require_value(&mut args, &arg)?)?,
            "--overwrite" => overwrite = true,
            "--entities" => {
                let value = require_value(&mut args, &arg)?;
                entities = value
                    .split(',')
                    .map(|s| s.trim().to_string())
                    .filter(|s| !s.is_empty())
                    .collect();
            }
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
        cache_root: cache_root
            .ok_or_else(|| DataFusionError::Execution("--cache-root is required".into()))?,
        output_dir: output_dir
            .ok_or_else(|| DataFusionError::Execution("--output-dir is required".into()))?,
        chrom: chrom.ok_or_else(|| DataFusionError::Execution("--chrom is required".into()))?,
        cache_source_type,
        partitions: partitions.max(1),
        overwrite,
        entities,
    })
}

fn require_value(args: &mut impl Iterator<Item = String>, flag: &str) -> Result<String> {
    args.next()
        .ok_or_else(|| DataFusionError::Execution(format!("{flag} requires a value")))
}

fn parse_usize(value: String) -> Result<usize> {
    value
        .parse::<usize>()
        .map_err(|_| DataFusionError::Execution(format!("invalid integer: {value}")))
}

fn print_usage() {
    eprintln!(
        "build_lance_chrom_context --cache-root <dir> --output-dir <dir> --chrom <chrom> \
         [--cache-source-type ensembl|merged|refseq] [--partitions N] [--overwrite] \
         [--entities transcript,translation]"
    );
}
