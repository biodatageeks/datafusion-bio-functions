//! Build one chromosome of a plugin cache from a source manifest.
//!
//! ```text
//! cargo run -p datafusion-bio-function-vep --example build_plugin -- \
//!   --manifest <vepyr-plugins>/plugins/alphamissense/alphamissense.source.toml \
//!   --source-path /tmp/am_chr22.tsv.gz \
//!   --variation-shard <cache>/variation/chr22.parquet \
//!   --out /tmp/plugin_cache --chrom 22 \
//!   --af-max-sql "coalesce(minor_allele_freq, 0.0)"
//! ```

use std::path::PathBuf;

use datafusion::common::{DataFusionError, Result};
use datafusion_bio_function_vep::plugin_cache::build::build_plugin_chrom;
use datafusion_bio_function_vep::plugin_cache::cache_manifest::CacheManifest;
use datafusion_bio_function_vep::plugin_cache::source_manifest::SourceManifest;

fn arg(args: &[String], key: &str) -> Option<String> {
    args.iter()
        .position(|a| a == key)
        .and_then(|i| args.get(i + 1))
        .cloned()
}

#[tokio::main(flavor = "multi_thread")]
async fn main() -> Result<()> {
    let args: Vec<String> = std::env::args().collect();
    let manifest_path = arg(&args, "--manifest")
        .ok_or_else(|| DataFusionError::Execution("--manifest required".into()))?;
    let variation_shard = PathBuf::from(
        arg(&args, "--variation-shard")
            .ok_or_else(|| DataFusionError::Execution("--variation-shard required".into()))?,
    );
    let out = PathBuf::from(
        arg(&args, "--out").ok_or_else(|| DataFusionError::Execution("--out required".into()))?,
    );
    let chrom = arg(&args, "--chrom")
        .ok_or_else(|| DataFusionError::Execution("--chrom required".into()))?;
    let af_max_sql =
        arg(&args, "--af-max-sql").unwrap_or_else(|| "coalesce(minor_allele_freq, 0.0)".into());

    let mut manifest = SourceManifest::load(&PathBuf::from(&manifest_path))?;
    if let Some(source_path) = arg(&args, "--source-path")
        && let Some(first) = manifest.sources.first_mut()
    {
        first.path = source_path;
    }

    let manifest_file = PathBuf::from(&manifest_path)
        .file_name()
        .map(|s| s.to_string_lossy().into_owned())
        .unwrap_or_else(|| manifest_path.clone());

    println!(
        "Building plugin '{}' chrom {chrom} from '{}'",
        manifest.plugin_name, manifest.sources[0].path
    );
    let entry = build_plugin_chrom(
        &manifest,
        &manifest_file,
        &variation_shard,
        &af_max_sql,
        &out,
        &chrom,
    )
    .await?;
    println!(
        "  rows={} warm={} cold={} -> plugin/{}/{}",
        entry.rows, entry.warm, entry.cold, manifest.plugin_name, entry.file
    );

    // Write / update the cache manifest for this plugin.
    let plugin_dir = out.join("plugin").join(&manifest.plugin_name);
    let mut cache = CacheManifest::from_source(&manifest, &manifest_file);
    cache.chroms = vec![entry];
    cache.write(&plugin_dir)?;
    println!("  wrote {}", plugin_dir.join("manifest.json").display());
    Ok(())
}
