//! Rebuild one chromosome of a context (scan) cache entity as dictionary-enabled
//! Parquet: transcript, exon, regulatory, motif, or translation_core.
//!
//! Writes `parquet.<entity>/<chrom>.parquet` + upserts `chrom_manifest.json`.
//!
//! Usage:
//!   cargo run --release -p datafusion-bio-function-vep \
//!     --features lance-cache,cache-builder \
//!     --example build_parquet_context_chrom -- \
//!     --entity transcript \
//!     --cache-root /path/to/homo_sapiens_merged/115_GRCh38 \
//!     --output-dir /path/to/115_GRCh38_merged \
//!     --chrom chr1 --cache-source-type merged --cache-version 115 \
//!     --partitions 8 --overwrite

use std::path::{Path, PathBuf};
use std::str::FromStr;
use std::time::Instant;

use datafusion::common::{DataFusionError, Result};
use datafusion_bio_format_ensembl_cache::{CacheSourceType, EnsemblEntityKind};
use datafusion_bio_function_vep::cache::build::{
    CacheBuildOptions, build_parquet_context_entity_chrom, build_parquet_translation_core_chrom,
};
use datafusion_bio_function_vep::cache::manifest::{
    CHROM_MANIFEST_FILE, ChromDatasetEntry, ChromManifest,
};
use datafusion_bio_function_vep::vep_semantics::target_for_cache_version;

#[derive(Debug)]
struct Args {
    entity: String,
    cache_root: PathBuf,
    output_dir: PathBuf,
    chrom: String,
    cache_source_type: CacheSourceType,
    cache_version: String,
    partitions: usize,
    overwrite: bool,
}

/// Map an entity CLI name to its `EnsemblEntityKind` and `parquet.<dir>` name.
/// `translation_core` is special-cased by the caller.
fn entity_kind(entity: &str) -> Result<(EnsemblEntityKind, &'static str)> {
    match entity {
        "transcript" => Ok((EnsemblEntityKind::Transcript, "transcript")),
        "exon" => Ok((EnsemblEntityKind::Exon, "exon")),
        "regulatory" => Ok((EnsemblEntityKind::RegulatoryFeature, "regulatory")),
        "motif" => Ok((EnsemblEntityKind::MotifFeature, "motif")),
        other => Err(DataFusionError::Execution(format!(
            "unknown --entity '{other}' (expected transcript, exon, regulatory, motif, translation_core)"
        ))),
    }
}

#[tokio::main]
async fn main() -> Result<()> {
    let _ = env_logger::try_init();
    let args = parse_args()?;
    let started = Instant::now();
    let options = CacheBuildOptions {
        cache_root: args.cache_root.to_string_lossy().to_string(),
        output_dir: args.output_dir.to_string_lossy().to_string(),
        partitions: args.partitions,
        cache_source_type: args.cache_source_type,
        cache_version: args.cache_version.clone(),
        overwrite: args.overwrite,
        chrom_filter: None,
    };

    eprintln!("entity={}", args.entity);
    eprintln!("cache_root={}", args.cache_root.display());
    eprintln!("output_dir={}", args.output_dir.display());
    eprintln!("chrom={}", args.chrom);
    eprintln!("cache_source_type={}", args.cache_source_type);
    eprintln!("cache_version={}", args.cache_version);

    let (entry, dir_name) = if args.entity == "translation_core" {
        (
            build_parquet_translation_core_chrom(&options, &args.chrom).await?,
            "translation_core",
        )
    } else {
        let (kind, dir_name) = entity_kind(&args.entity)?;
        (
            build_parquet_context_entity_chrom(&options, kind, &args.chrom).await?,
            dir_name,
        )
    };

    if entry.rows == 0 {
        eprintln!(
            "WARNING: {} {} wrote 0 rows (no shard)",
            args.entity, args.chrom
        );
        return Ok(());
    }

    let entity_dir = args.output_dir.join(dir_name);
    upsert_manifest_entry(&entity_dir, entry.clone())?;

    eprintln!(
        "wrote entity={} chrom={} shard={} rows={} elapsed={:.2}s",
        args.entity,
        entry.chrom,
        entity_dir.join(&entry.dataset).display(),
        entry.rows,
        started.elapsed().as_secs_f64()
    );
    Ok(())
}

fn parse_args() -> Result<Args> {
    let mut entity = None;
    let mut cache_root = None;
    let mut output_dir = None;
    let mut chrom = None;
    let mut cache_source_type = CacheSourceType::Ensembl;
    let mut cache_version = None;
    let mut partitions = 8usize;
    let mut overwrite = false;

    let mut args = std::env::args().skip(1);
    while let Some(arg) = args.next() {
        match arg.as_str() {
            "--entity" => entity = Some(require_value(&mut args, &arg)?),
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
            "--cache-version" => cache_version = Some(require_value(&mut args, &arg)?),
            "--partitions" => partitions = parse_usize(require_value(&mut args, &arg)?)?,
            "--overwrite" => overwrite = true,
            "--help" | "-h" => {
                eprintln!(
                    "Usage: build_parquet_context_chrom --entity <transcript|exon|regulatory|motif|translation_core> \
                     --cache-root <dir> --output-dir <dir> --chrom chr1 --cache-version 115 \
                     [--cache-source-type merged] [--partitions 8] [--overwrite]"
                );
                std::process::exit(0);
            }
            other => {
                return Err(DataFusionError::Execution(format!(
                    "unknown argument: {other}"
                )));
            }
        }
    }

    let cache_version = cache_version
        .ok_or_else(|| DataFusionError::Execution("--cache-version is required".into()))?;
    target_for_cache_version(&cache_version)?;

    Ok(Args {
        entity: entity.ok_or_else(|| DataFusionError::Execution("--entity is required".into()))?,
        cache_root: cache_root
            .ok_or_else(|| DataFusionError::Execution("--cache-root is required".into()))?,
        output_dir: output_dir
            .ok_or_else(|| DataFusionError::Execution("--output-dir is required".into()))?,
        chrom: chrom.ok_or_else(|| DataFusionError::Execution("--chrom is required".into()))?,
        cache_source_type,
        cache_version,
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
