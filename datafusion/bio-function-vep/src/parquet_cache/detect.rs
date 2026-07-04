//! Partitioned per-chromosome Parquet cache detection.
//!
//! Each entity is a `parquet.<entity>/` directory with a `chrom_manifest.json`
//! mapping chromosomes to per-chromosome `.parquet` files:
//!
//! ```text
//! 115_GRCh38_vep/
//!   parquet.variation/chrom_manifest.json
//!   parquet.variation/chr1.parquet
//!   parquet.transcript/chr1.parquet
//!   ...
//! ```
//!
//! The `ChromManifest` machinery is format-agnostic and reused verbatim.

use std::path::{Path, PathBuf};

use crate::lance_cache::manifest::ChromManifest;

/// Directory name for a Parquet cache entity (e.g. `parquet.variation`).
fn entity_dir_name(entity: &str) -> String {
    format!("parquet.{entity}")
}

/// A partitioned per-chromosome Parquet cache directory.
#[derive(Debug, Clone)]
pub struct PartitionedParquetCache {
    base_dir: PathBuf,
    variation_manifest: ChromManifest,
}

impl PartitionedParquetCache {
    /// Detect a Parquet cache layout at `cache_source`.
    ///
    /// Returns `Some` when `parquet.variation/chrom_manifest.json` can be read.
    pub fn detect(cache_source: &str) -> Option<Self> {
        let base_dir = PathBuf::from(cache_source);
        let variation_dir = base_dir.join(entity_dir_name("variation"));
        let variation_manifest = ChromManifest::read_from_entity_dir(&variation_dir).ok()?;
        Some(Self {
            base_dir,
            variation_manifest,
        })
    }

    pub fn base_dir(&self) -> &Path {
        &self.base_dir
    }

    pub fn available_chroms(&self) -> Vec<&str> {
        self.variation_manifest.available_chroms()
    }

    /// Path to the variation `.parquet` shard for `chrom`, if present.
    pub fn variation_path(&self, chrom: &str) -> Option<PathBuf> {
        self.variation_manifest
            .path_for_chrom(chrom)
            .map(|path| self.base_dir.join(entity_dir_name("variation")).join(path))
    }

    /// Path to a context entity's `.parquet` shard for `chrom`, if present.
    pub fn context_path(&self, context_type: &str, chrom: &str) -> Option<PathBuf> {
        let entity_dir = self.base_dir.join(entity_dir_name(context_type));
        let manifest = ChromManifest::read_from_entity_dir(&entity_dir).ok()?;
        manifest
            .path_for_chrom(chrom)
            .map(|path| entity_dir.join(path))
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::lance_cache::manifest::{ChromDatasetEntry, ChromManifest};

    fn write_entity(base: &Path, entity: &str, chrom: &str, file: &str) {
        let dir = base.join(format!("parquet.{entity}"));
        std::fs::create_dir_all(&dir).unwrap();
        std::fs::write(dir.join(file), b"parquet-shard-placeholder").unwrap();
        ChromManifest::new(vec![ChromDatasetEntry::new(chrom, file, 1)])
            .write_to_entity_dir(&dir)
            .unwrap();
    }

    #[test]
    fn detects_partitioned_parquet_cache_from_manifest() {
        let tmp = tempfile::tempdir().unwrap();
        write_entity(tmp.path(), "variation", "chr1", "chr1.parquet");

        let cache = PartitionedParquetCache::detect(tmp.path().to_str().unwrap()).unwrap();
        assert_eq!(cache.available_chroms(), ["chr1"]);
        assert_eq!(
            cache.variation_path("chr1").unwrap(),
            tmp.path().join("parquet.variation").join("chr1.parquet")
        );
        assert!(cache.variation_path("chr2").is_none());
    }

    #[test]
    fn detect_returns_none_without_variation_manifest() {
        let tmp = tempfile::tempdir().unwrap();
        assert!(PartitionedParquetCache::detect(tmp.path().to_str().unwrap()).is_none());
    }

    #[test]
    fn resolves_context_paths_from_entity_manifests() {
        let tmp = tempfile::tempdir().unwrap();
        write_entity(tmp.path(), "variation", "chr1", "chr1.parquet");
        write_entity(tmp.path(), "transcript", "chr1", "chr1.parquet");

        let cache = PartitionedParquetCache::detect(tmp.path().to_str().unwrap()).unwrap();
        assert_eq!(
            cache.context_path("transcript", "chr1").unwrap(),
            tmp.path().join("parquet.transcript").join("chr1.parquet")
        );
        assert!(cache.context_path("transcript", "chr2").is_none());
    }
}
