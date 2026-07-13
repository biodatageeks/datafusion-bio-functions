//! Build a plugin cache from a source manifest (all chroms, or a filtered set).
//!
//! ```text
//! cargo run -p datafusion-bio-function-vep --features parquet-cache --example build_plugin -- \
//!   --manifest <vepyr-plugins>/plugins/alphamissense/alphamissense.source.toml \
//!   --source-path /tmp/AlphaMissense_hg38.tsv.gz \
//!   --variation-cache-dir <cache root containing variation/> \
//!   --out /tmp/plugin_cache \
//!   [--chrom 21 --chrom 22]        # repeatable; omit to build every chrom in the cache
//!   [--overwrite]                  # clean rebuild instead of an UPSERT
//! ```
//!
//! `--source-path` takes a bare path only when the manifest declares a single `[[source]]`;
//! for multi-part manifests use `--source-path <part>=<path>`, repeated once per part.
//!
//! `--overwrite` starts from an empty chrom list (a clean rebuild). Without it, a filtered
//! build UPSERTs into the previous `manifest.json`, preserving chroms it did not rebuild.
//!
//! This example is deliberately thin: all argv handling lives in
//! `plugin_cache::cli::PluginBuildArgs`, where it is unit-tested.

use datafusion::common::Result;
use datafusion_bio_function_vep::plugin_cache::cli::PluginBuildArgs;

#[tokio::main(flavor = "multi_thread")]
async fn main() -> Result<()> {
    let argv: Vec<String> = std::env::args().skip(1).collect();
    let args = PluginBuildArgs::parse(&argv)?;
    let manifest = args.load_manifest()?;

    println!(
        "Building plugin '{}' from '{}'",
        manifest.plugin_name, manifest.sources[0].path
    );
    let cache = args.build(&manifest).await?;
    for c in &cache.chroms {
        println!(
            "  {} rows={} warm={} cold={}",
            c.chrom, c.rows, c.warm, c.cold
        );
    }
    println!("  wrote plugin/{}/manifest.json", manifest.plugin_name);
    Ok(())
}
