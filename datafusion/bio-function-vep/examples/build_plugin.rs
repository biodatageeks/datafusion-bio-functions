//! Build a plugin cache from a source manifest (all chroms, or a filtered set).
//!
//! ```text
//! cargo run -p datafusion-bio-function-vep --example build_plugin -- \
//!   --manifest <vepyr-plugins>/plugins/alphamissense/alphamissense.source.toml \
//!   --source-path /tmp/AlphaMissense_hg38.bgz.tsv.gz \
//!   --variation-cache-dir <cache root containing variation/> \
//!   --out /tmp/plugin_cache [--chrom 22] [--verify-source strict|warn|skip]
//!
//! Multi-source manifests use one `--source-part <name>=<path>` per part.
//! Each source is MD5-hashed once and checked against the manifest's
//! `md5` before anything is built (`--verify-source strict`, the
//! default). Pass `warn` for a deliberate build from a re-compressed or
//! derived file, or `skip` for a chromosome slice whose digest can never
//! match the whole file's.
//! Indexed TSV manifests query the requested chromosome natively through the
//! source's sibling `.tbi`; no external `tabix` process or pre-slice is needed.
//! ```

use std::path::PathBuf;

use datafusion::common::{DataFusionError, Result};
use datafusion_bio_function_vep::plugin_cache::builder::PluginCacheBuilder;
use datafusion_bio_function_vep::plugin_cache::source_manifest::SourceManifest;
use datafusion_bio_function_vep::plugin_cache::source_verify::SourceVerification;

fn arg(args: &[String], key: &str) -> Option<String> {
    args.iter()
        .position(|a| a == key)
        .and_then(|i| args.get(i + 1))
        .cloned()
}

fn arg_values(args: &[String], key: &str) -> Vec<String> {
    args.iter()
        .enumerate()
        .filter(|(_, value)| value.as_str() == key)
        .filter_map(|(index, _)| args.get(index + 1).cloned())
        .collect()
}

#[tokio::main(flavor = "multi_thread")]
async fn main() -> Result<()> {
    let _ = env_logger::try_init();
    let args: Vec<String> = std::env::args().collect();
    let manifest_path = arg(&args, "--manifest")
        .ok_or_else(|| DataFusionError::Execution("--manifest required".into()))?;
    let variation_cache_dir = PathBuf::from(
        arg(&args, "--variation-cache-dir")
            .ok_or_else(|| DataFusionError::Execution("--variation-cache-dir required".into()))?,
    );
    let out = PathBuf::from(
        arg(&args, "--out").ok_or_else(|| DataFusionError::Execution("--out required".into()))?,
    );

    let mut manifest = SourceManifest::load(&PathBuf::from(&manifest_path))?;
    if let Some(source_path) = arg(&args, "--source-path")
        && let Some(first) = manifest.sources.first_mut()
    {
        first.path = source_path;
    }
    for source_part in arg_values(&args, "--source-part") {
        let (part, path) = source_part.split_once('=').ok_or_else(|| {
            DataFusionError::Execution(format!(
                "--source-part expects <name>=<path>, got '{source_part}'"
            ))
        })?;
        let source = manifest
            .sources
            .iter_mut()
            .find(|source| source.part.as_deref() == Some(part))
            .ok_or_else(|| {
                DataFusionError::Execution(format!("manifest has no [[source]] with part='{part}'"))
            })?;
        source.path = path.to_string();
    }
    let manifest_file = PathBuf::from(&manifest_path)
        .file_name()
        .map(|s| s.to_string_lossy().into_owned())
        .unwrap_or_else(|| manifest_path.clone());

    println!(
        "Building plugin '{}' from {:?}",
        manifest.plugin_name,
        manifest
            .sources
            .iter()
            .map(|source| (&source.part, &source.path))
            .collect::<Vec<_>>()
    );
    let mut builder =
        PluginCacheBuilder::new(&manifest, &manifest_file, &variation_cache_dir, &out);
    if let Some(chrom) = arg(&args, "--chrom") {
        builder = builder.with_chrom_filter([chrom]);
    }
    if let Some(mode) = arg(&args, "--verify-source") {
        let verification: SourceVerification = mode
            .parse()
            .map_err(|e| DataFusionError::Execution(format!("--verify-source: {e}")))?;
        builder = builder.with_source_verification(verification);
    }
    let cache = builder.build_all().await?;
    for c in &cache.chroms {
        println!(
            "  {} rows={} warm={} cold={}",
            c.chrom, c.rows, c.warm, c.cold
        );
    }
    println!("  wrote plugin/{}/manifest.json", manifest.plugin_name);
    Ok(())
}
