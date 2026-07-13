//! Build a plugin cache from a source manifest (all chroms, or a filtered set).
//!
//! ```text
//! cargo run -p datafusion-bio-function-vep --example build_plugin -- \
//!   --manifest <vepyr-plugins>/plugins/alphamissense/alphamissense.source.toml \
//!   --source-path /tmp/AlphaMissense_hg38.tsv.gz \
//!   --variation-cache-dir <cache root containing variation/> \
//!   --out /tmp/plugin_cache [--chrom 22]
//! ```

use std::path::PathBuf;

use datafusion::common::{DataFusionError, Result};
use datafusion_bio_function_vep::plugin_cache::builder::PluginCacheBuilder;
use datafusion_bio_function_vep::plugin_cache::source_manifest::SourceManifest;

fn arg(args: &[String], key: &str) -> Option<String> {
    args.iter()
        .position(|a| a == key)
        .and_then(|i| args.get(i + 1))
        .cloned()
}

/// Every occurrence of `--key <value>` (the flag may repeat).
fn args_all(args: &[String], key: &str) -> Vec<String> {
    args.iter()
        .enumerate()
        .filter(|(_, a)| a.as_str() == key)
        .filter_map(|(i, _)| args.get(i + 1).cloned())
        .collect()
}

#[tokio::main(flavor = "multi_thread")]
async fn main() -> Result<()> {
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
    for spec in args_all(&args, "--source-path") {
        match spec.split_once('=') {
            // --source-path <part>=<path>
            Some((part, path)) => {
                let target = manifest
                    .sources
                    .iter_mut()
                    .find(|s| s.part.as_deref() == Some(part))
                    .ok_or_else(|| {
                        DataFusionError::Execution(format!(
                            "--source-path '{part}=...': no [[source]] with part = \"{part}\""
                        ))
                    })?;
                target.path = path.to_string();
            }
            // --source-path <path> — unambiguous only for a single-source manifest
            None => {
                if manifest.sources.len() != 1 {
                    return Err(DataFusionError::Execution(format!(
                        "--source-path <path> is ambiguous: the manifest has {} sources; \
                         use --source-path <part>=<path> (parts: {})",
                        manifest.sources.len(),
                        manifest
                            .sources
                            .iter()
                            .map(|s| s.part.as_deref().unwrap_or("<none>"))
                            .collect::<Vec<_>>()
                            .join(", ")
                    )));
                }
                manifest.sources[0].path = spec;
            }
        }
    }
    let manifest_file = PathBuf::from(&manifest_path)
        .file_name()
        .map(|s| s.to_string_lossy().into_owned())
        .unwrap_or_else(|| manifest_path.clone());

    println!(
        "Building plugin '{}' from '{}'",
        manifest.plugin_name, manifest.sources[0].path
    );
    let mut builder =
        PluginCacheBuilder::new(&manifest, &manifest_file, &variation_cache_dir, &out);
    let chroms = args_all(&args, "--chrom");
    if !chroms.is_empty() {
        builder = builder.with_chrom_filter(chroms);
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
