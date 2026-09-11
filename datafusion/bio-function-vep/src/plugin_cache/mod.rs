//! Custom VEP plugin caches: declarative build + tiered point-lookup.
//!
//! See `docs/superpowers/specs/2026-07-05-custom-vep-plugin-caches-design.md`.
//! A TOML source manifest declares a table provider (`csv`/`tsv`/`parquet`/
//! `vcf`/`bed`/`gff`) + `ingest_sql`; the builder registers the raw source,
//! normalizes contig/coordinates to the variation convention, LEFT-joins the
//! variation shard to inherit its warm/cold `tier`, and writes
//! `plugin/<name>/<chrom>.parquet` reusing the variation shard's
//! lookup-optimized writer properties. A manifest with `lookup = "interval"`
//! (gene or region tracks) skips the allele and the tier join: its shards are
//! probed at runtime through per-discriminator COITrees by overlap with the
//! variant's VEP-normalised span, first row in file order winning.

pub mod build;
pub mod builder;
pub mod cache_manifest;
pub mod csq;
pub mod dedup;
pub mod join;
mod join_strategy;
pub mod lookup;
pub mod mem_trace;
pub mod normalize;
pub mod provider;
pub mod registry;
pub mod source_manifest;
pub mod source_verify;
pub mod template;
pub mod write;
