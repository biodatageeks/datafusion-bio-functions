//! Partitioned per-chromosome Lance cache detection and registration.
//!
//! The partitioned Lance layout stores each context type in a `*.lance`
//! entity directory with a `chrom_manifest.json` mapping chromosomes to
//! per-chromosome Lance datasets:
//!
//! ```text
//! 115_GRCh38_vep/
//!   variation.lance/chrom_manifest.json
//!   variation.lance/chr1.lance/
//!   transcript.lance/chr1.lance/
//!   exon.lance/chr1.lance/
//!   translation_core.lance/chr1.lance/
//!   translation_sift.lance/chr1.lance/
//!   regulatory.lance/chr1.lance/
//!   motif.lance/chr1.lance/
//! ```

use std::path::{Path, PathBuf};

use datafusion::common::Result;
use datafusion::prelude::SessionContext;

/// Represents a partitioned per-chromosome Lance cache directory.
#[cfg(feature = "lance-cache")]
#[derive(Debug, Clone)]
pub struct PartitionedLanceCache {
    base_dir: PathBuf,
    variation_manifest: crate::lance_cache::manifest::ChromManifest,
}

#[cfg(feature = "lance-cache")]
impl PartitionedLanceCache {
    /// Detect a Lance cache layout at `cache_source`.
    ///
    /// Returns `Some` when `variation.lance/chrom_manifest.json` can be read.
    pub fn detect(cache_source: &str) -> Option<Self> {
        let base_dir = PathBuf::from(cache_source);
        let variation_dir = base_dir.join("variation.lance");
        let variation_manifest =
            crate::lance_cache::manifest::ChromManifest::read_from_entity_dir(&variation_dir)
                .ok()?;
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

    pub fn variation_path(&self, chrom: &str) -> Option<PathBuf> {
        self.variation_manifest
            .path_for_chrom(chrom)
            .map(|path| self.base_dir.join("variation.lance").join(path))
    }

    pub fn context_path(&self, context_type: &str, chrom: &str) -> Option<PathBuf> {
        let entity_dir = self.base_dir.join(format!("{context_type}.lance"));
        let manifest =
            crate::lance_cache::manifest::ChromManifest::read_from_entity_dir(&entity_dir).ok()?;
        manifest
            .path_for_chrom(chrom)
            .map(|path| entity_dir.join(path))
    }
}

/// Deregister an ephemeral table from the session.
pub async fn deregister_table(session: &SessionContext, name: &str) -> Result<()> {
    // deregister_table returns Option<Arc<dyn TableProvider>>; ignore it.
    let _ = session.deregister_table(name)?;
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;

    #[cfg(feature = "lance-cache")]
    #[test]
    fn detects_partitioned_lance_cache_from_manifest() {
        let tmp = tempfile::tempdir().unwrap();
        let variation = tmp.path().join("variation.lance");
        std::fs::create_dir_all(variation.join("chr1.lance")).unwrap();
        let manifest = crate::lance_cache::manifest::ChromManifest::new(vec![
            crate::lance_cache::manifest::ChromDatasetEntry::new("chr1", "chr1.lance", 1),
        ]);
        manifest.write_to_entity_dir(&variation).unwrap();

        let cache = PartitionedLanceCache::detect(tmp.path().to_str().unwrap()).unwrap();

        assert_eq!(cache.available_chroms(), ["chr1"]);
        assert_eq!(
            cache.variation_path("chr1").unwrap(),
            variation.join("chr1.lance")
        );
        assert!(cache.variation_path("chr2").is_none());
    }

    #[cfg(feature = "lance-cache")]
    #[test]
    fn resolves_partitioned_lance_context_paths_from_entity_manifests() {
        let tmp = tempfile::tempdir().unwrap();
        let variation = tmp.path().join("variation.lance");
        std::fs::create_dir_all(variation.join("chr1.lance")).unwrap();
        crate::lance_cache::manifest::ChromManifest::new(vec![
            crate::lance_cache::manifest::ChromDatasetEntry::new("chr1", "chr1.lance", 1),
        ])
        .write_to_entity_dir(&variation)
        .unwrap();

        let transcript = tmp.path().join("transcript.lance");
        std::fs::create_dir_all(transcript.join("chr1.lance")).unwrap();
        crate::lance_cache::manifest::ChromManifest::new(vec![
            crate::lance_cache::manifest::ChromDatasetEntry::new("chr1", "chr1.lance", 1),
        ])
        .write_to_entity_dir(&transcript)
        .unwrap();

        let cache = PartitionedLanceCache::detect(tmp.path().to_str().unwrap()).unwrap();
        assert_eq!(
            cache.context_path("transcript", "chr1").unwrap(),
            transcript.join("chr1.lance")
        );
        assert!(cache.context_path("transcript", "chr2").is_none());
    }
}
