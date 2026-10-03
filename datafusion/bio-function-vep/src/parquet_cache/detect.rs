//! Partitioned per-chromosome Parquet cache detection.
//!
//! Each entity is a `parquet.<entity>/` directory with a `chrom_manifest.json`
//! mapping chromosomes to per-chromosome `.parquet` files:
//!
//! ```text
//! 115_GRCh38_vep/
//!   variation/chrom_manifest.json
//!   variation/chr1.parquet
//!   transcript/chr1.parquet
//!   ...
//! ```
//!
//! The `ChromManifest` machinery is format-agnostic and reused verbatim.

use std::collections::HashSet;
use std::path::{Path, PathBuf};
use std::sync::Arc;

use crate::cache::manifest::ChromManifest;
use crate::cache::synonyms::ChromosomeSynonyms;

/// Directory name for a Parquet cache entity (e.g. `variation`).
fn entity_dir_name(entity: &str) -> String {
    entity.to_string()
}

/// A partitioned per-chromosome Parquet cache directory.
#[derive(Debug, Clone)]
pub struct PartitionedParquetCache {
    base_dir: PathBuf,
    variation_manifest: ChromManifest,
    synonyms: Arc<ChromosomeSynonyms>,
}

impl PartitionedParquetCache {
    /// Detect a Parquet cache layout at `cache_source`.
    ///
    /// Returns `Some` when `variation/chrom_manifest.json` can be read.
    pub fn detect(cache_source: &str) -> Option<Self> {
        let base_dir = PathBuf::from(cache_source);
        let variation_dir = base_dir.join(entity_dir_name("variation"));
        let variation_manifest = ChromManifest::read_from_entity_dir(&variation_dir).ok()?;
        let synonyms = Arc::new(ChromosomeSynonyms::read(&base_dir).ok()?);
        Some(Self {
            base_dir,
            variation_manifest,
            synonyms,
        })
    }

    pub fn base_dir(&self) -> &Path {
        &self.base_dir
    }

    pub fn available_chroms(&self) -> Vec<&str> {
        self.variation_manifest.available_chroms()
    }

    /// Resolve only against this cache's actual contigs; accessions are never
    /// decoded into chromosome numbers or inferred from another assembly.
    pub fn resolve_chrom(&self, chrom: &str) -> Option<&str> {
        self.resolve_in_manifest(&self.variation_manifest, chrom)
    }

    fn resolve_in_manifest<'a>(&self, manifest: &'a ChromManifest, chrom: &str) -> Option<&'a str> {
        let ordinary = manifest.chrom_for_alias(chrom);
        if ordinary == Some(chrom) {
            return ordinary;
        }
        self.synonyms
            .aliases(chrom)
            .find_map(|alias| manifest.chrom_for_alias(alias))
            .or(ordinary)
    }

    pub(crate) fn accepted_chroms(&self) -> HashSet<String> {
        let mut names: HashSet<String> = self
            .available_chroms()
            .into_iter()
            .flat_map(crate::cache::manifest::contig_alias_set)
            .collect();
        names.extend(
            self.synonyms
                .names()
                .filter(|name| self.resolve_chrom(name).is_some())
                .cloned(),
        );
        names
    }

    /// Ordinary bare/chr/mitochondrial handling stays on its existing path.
    pub(crate) fn annotation_chrom_override(&self, chrom: &str) -> Option<String> {
        let resolved = self.resolve_chrom(chrom)?;
        (!crate::cache::manifest::contig_alias_set(resolved).contains(chrom))
            .then(|| resolved.to_owned())
    }

    /// Path to the variation `.parquet` shard for `chrom`, if present.
    pub fn variation_path(&self, chrom: &str) -> Option<PathBuf> {
        self.variation_manifest
            .path_for_chrom(self.resolve_chrom(chrom)?)
            .map(|path| self.base_dir.join(entity_dir_name("variation")).join(path))
    }

    /// Path to a context entity's `.parquet` shard for `chrom`, if present.
    pub fn context_path(&self, context_type: &str, chrom: &str) -> Option<PathBuf> {
        let entity_dir = self.base_dir.join(entity_dir_name(context_type));
        let manifest = ChromManifest::read_from_entity_dir(&entity_dir).ok()?;
        manifest
            .path_for_chrom(self.resolve_in_manifest(&manifest, chrom)?)
            .map(|path| entity_dir.join(path))
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::cache::manifest::{ChromDatasetEntry, ChromManifest};

    fn write_entity(base: &Path, entity: &str, chrom: &str, file: &str) {
        let dir = base.join(entity);
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
            tmp.path().join("variation").join("chr1.parquet")
        );
        assert!(cache.variation_path("chr2").is_none());
    }

    #[test]
    fn detect_returns_none_without_variation_manifest() {
        let tmp = tempfile::tempdir().unwrap();
        assert!(PartitionedParquetCache::detect(tmp.path().to_str().unwrap()).is_none());
    }

    #[test]
    fn resolves_explicit_cache_synonyms_for_variation_and_context() {
        let tmp = tempfile::tempdir().unwrap();
        write_entity(tmp.path(), "variation", "chr21", "chr21.parquet");
        write_entity(tmp.path(), "transcript", "chr21", "chr21.parquet");
        std::fs::write(
            tmp.path().join("chr_synonyms.txt"),
            "21 NC_000021.9 assembly_accession\nreverse_accession 21\nassembly_accession indirect_only\n",
        )
        .unwrap();
        let cache = PartitionedParquetCache::detect(tmp.path().to_str().unwrap()).unwrap();
        for alias in ["NC_000021.9", "assembly_accession", "reverse_accession"] {
            assert_eq!(
                cache.variation_path(alias),
                cache.variation_path("21"),
                "{alias}"
            );
            assert_eq!(
                cache.context_path("transcript", alias),
                cache.context_path("transcript", "21"),
                "{alias}"
            );
        }
        // VEP's file reader creates direct bidirectional pairs, not graph closure.
        assert!(cache.variation_path("indirect_only").is_none());
        assert!(cache.variation_path("unknown_accession").is_none());
    }

    #[test]
    fn exact_cache_contig_wins_over_a_synonym() {
        let tmp = tempfile::tempdir().unwrap();
        write_entity(tmp.path(), "variation", "primary", "primary.parquet");
        crate::cache::manifest::ChromManifest::new(vec![
            ChromDatasetEntry::new("primary", "primary.parquet", 1),
            ChromDatasetEntry::new("other", "other.parquet", 1),
        ])
        .write_to_entity_dir(&tmp.path().join("variation"))
        .unwrap();
        std::fs::write(tmp.path().join("chr_synonyms.txt"), "primary other\n").unwrap();
        let cache = PartitionedParquetCache::detect(tmp.path().to_str().unwrap()).unwrap();
        assert_eq!(
            cache.variation_path("primary").unwrap(),
            tmp.path().join("variation/primary.parquet")
        );
        assert_eq!(
            cache.variation_path("other").unwrap(),
            tmp.path().join("variation/other.parquet")
        );
    }

    #[test]
    fn resolves_context_paths_from_entity_manifests() {
        let tmp = tempfile::tempdir().unwrap();
        write_entity(tmp.path(), "variation", "chr1", "chr1.parquet");
        write_entity(tmp.path(), "transcript", "chr1", "chr1.parquet");

        let cache = PartitionedParquetCache::detect(tmp.path().to_str().unwrap()).unwrap();
        assert_eq!(
            cache.context_path("transcript", "chr1").unwrap(),
            tmp.path().join("transcript").join("chr1.parquet")
        );
        assert!(cache.context_path("transcript", "chr2").is_none());
    }
}
