//! Rebuild one chromosome of the partitioned Lance variation cache from the
//! Ensembl-cache table provider.
//!
//! Usage:
//!   cargo run --release -p datafusion-bio-function-vep \
//!     --features lance-cache,cache-builder \
//!     --example build_lance_variation_chrom -- \
//!     --cache-root /path/to/homo_sapiens_merged/115_GRCh38 \
//!     --output-dir /path/to/115_GRCh38_merged \
//!     --chrom chr1 --cache-source-type merged --partitions 8 --overwrite

use std::path::{Path, PathBuf};
use std::str::FromStr;
use std::time::Instant;

use datafusion::common::{DataFusionError, Result};
use datafusion_bio_format_ensembl_cache::CacheSourceType;
use datafusion_bio_function_vep::lance_cache::build::{
    LanceCacheBuildOptions, build_lance_variation_chrom, lance_entity_dir_name,
};
use datafusion_bio_function_vep::lance_cache::manifest::{
    CHROM_MANIFEST_FILE, ChromDatasetEntry, ChromManifest,
};

#[derive(Debug)]
struct Args {
    cache_root: PathBuf,
    output_dir: PathBuf,
    chrom: String,
    cache_source_type: CacheSourceType,
    partitions: usize,
    overwrite: bool,
}

#[tokio::main]
async fn main() -> Result<()> {
    let args = parse_args()?;
    let started = Instant::now();
    let options = LanceCacheBuildOptions {
        cache_root: args.cache_root.to_string_lossy().to_string(),
        output_dir: args.output_dir.to_string_lossy().to_string(),
        partitions: args.partitions,
        cache_source_type: args.cache_source_type,
        overwrite: args.overwrite,
        chrom_filter: None,
    };

    eprintln!("cache_root={}", args.cache_root.display());
    eprintln!("output_dir={}", args.output_dir.display());
    eprintln!("chrom={}", args.chrom);
    eprintln!("cache_source_type={}", args.cache_source_type);
    eprintln!("partitions={}", args.partitions);
    eprintln!("overwrite={}", args.overwrite);

    let entry = build_lance_variation_chrom(&options, &args.chrom).await?;
    if entry.rows == 0 {
        return Err(DataFusionError::Execution(format!(
            "variation Lance build for {} wrote 0 rows",
            args.chrom
        )));
    }

    let entity_dir = args.output_dir.join(lance_entity_dir_name("variation"));
    upsert_manifest_entry(&entity_dir, entry.clone())?;

    eprintln!(
        "wrote chrom={} dataset={} rows={} elapsed={:.2}s",
        entry.chrom,
        entity_dir.join(&entry.dataset).display(),
        entry.rows,
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

    let mut args = std::env::args().skip(1);
    while let Some(arg) = args.next() {
        match arg.as_str() {
            "--cache-root" => cache_root = Some(PathBuf::from(require_value(&mut args, &arg)?)),
            "--output-dir" => output_dir = Some(PathBuf::from(require_value(&mut args, &arg)?)),
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
    })
}

fn upsert_manifest_entry(entity_dir: &Path, entry: ChromDatasetEntry) -> Result<()> {
    let manifest_path = entity_dir.join(CHROM_MANIFEST_FILE);
    let mut entries = if manifest_path.exists() {
        ChromManifest::read_from_entity_dir(entity_dir)?.entries
    } else {
        Vec::new()
    };
    if let Some(existing) = entries
        .iter_mut()
        .find(|existing| existing.chrom == entry.chrom)
    {
        *existing = entry;
    } else {
        entries.push(entry);
    }
    ChromManifest::new(entries).write_to_entity_dir(entity_dir)
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
        "Usage: build_lance_variation_chrom --cache-root /path/to/homo_sapiens_merged/115_GRCh38 \
         --output-dir /path/to/115_GRCh38_merged --chrom chr1 \
         [--cache-source-type merged] [--partitions 8] [--overwrite]"
    );
}
