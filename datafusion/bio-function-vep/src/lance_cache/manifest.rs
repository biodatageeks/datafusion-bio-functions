use std::collections::HashMap;
use std::fs;
use std::path::Path;

use datafusion::common::{DataFusionError, Result};
use serde::{Deserialize, Serialize};

pub const CHROM_MANIFEST_FILE: &str = "chrom_manifest.json";

#[derive(Debug, Clone, PartialEq, Eq, Serialize, Deserialize)]
pub struct ChromDatasetEntry {
    pub chrom: String,
    pub dataset: String,
    pub rows: usize,
}

impl ChromDatasetEntry {
    pub fn new(chrom: impl Into<String>, dataset: impl Into<String>, rows: usize) -> Self {
        Self {
            chrom: chrom.into(),
            dataset: dataset.into(),
            rows,
        }
    }
}

#[derive(Debug, Clone, PartialEq, Eq)]
pub struct ChromManifest {
    pub entries: Vec<ChromDatasetEntry>,
    by_chrom: HashMap<String, usize>,
}

impl ChromManifest {
    pub fn new(entries: Vec<ChromDatasetEntry>) -> Self {
        let by_chrom = entries
            .iter()
            .enumerate()
            .map(|(idx, entry)| (entry.chrom.clone(), idx))
            .collect();
        Self { entries, by_chrom }
    }

    pub fn read_from_entity_dir(entity_dir: &Path) -> Result<Self> {
        let path = entity_dir.join(CHROM_MANIFEST_FILE);
        let bytes = fs::read(&path).map_err(|err| {
            DataFusionError::Execution(format!("failed to read {}: {err}", path.display()))
        })?;
        let entries: Vec<ChromDatasetEntry> = serde_json::from_slice(&bytes).map_err(|err| {
            DataFusionError::Execution(format!("failed to parse {}: {err}", path.display()))
        })?;
        Ok(Self::new(entries))
    }

    pub fn write_to_entity_dir(&self, entity_dir: &Path) -> Result<()> {
        fs::create_dir_all(entity_dir).map_err(|err| {
            DataFusionError::Execution(format!("failed to create {}: {err}", entity_dir.display()))
        })?;
        let path = entity_dir.join(CHROM_MANIFEST_FILE);
        let json = serde_json::to_vec_pretty(&self.entries).map_err(|err| {
            DataFusionError::Execution(format!("failed to serialize chrom manifest: {err}"))
        })?;
        fs::write(&path, json).map_err(|err| {
            DataFusionError::Execution(format!("failed to write {}: {err}", path.display()))
        })
    }

    /// Resolve the per-chromosome dataset path, tolerating `chr`-prefix spelling
    /// differences between the query and the manifest (e.g. VCF contig `1` vs a
    /// manifest canonicalized to `chr1`, or vice versa). Exact match wins; the
    /// bare/`chr`-prefixed spellings are only tried as fallbacks. Mirrors the
    /// fallback previously applied only at the variation-lookup call site, so
    /// context entities (transcript/exon/translation/regulatory/motif/SIFT)
    /// resolve the same way instead of silently loading empty.
    pub fn path_for_chrom(&self, chrom: &str) -> Option<&str> {
        let idx = self.by_chrom.get(chrom).or_else(|| {
            chrom
                .strip_prefix("chr")
                .and_then(|bare| self.by_chrom.get(bare))
                .or_else(|| {
                    if chrom.starts_with("chr") {
                        None
                    } else {
                        self.by_chrom.get(&format!("chr{chrom}"))
                    }
                })
        })?;
        Some(self.entries[*idx].dataset.as_str())
    }

    pub fn available_chroms(&self) -> Vec<&str> {
        self.entries
            .iter()
            .map(|entry| entry.chrom.as_str())
            .collect()
    }
}

pub fn dataset_dir_name(chrom: &str) -> String {
    let label = canonical_chrom_label(chrom);
    let mut encoded = String::with_capacity(label.len() + ".lance".len());
    for byte in label.bytes() {
        let keep = byte.is_ascii_alphanumeric() || matches!(byte, b'_' | b'-' | b'.');
        if keep {
            encoded.push(byte as char);
        } else {
            encoded.push('_');
            encoded.push_str(&format!("{byte:02X}"));
        }
    }
    encoded.push_str(".lance");
    encoded
}

pub fn canonical_chrom_label(chrom: &str) -> String {
    let bare = chrom.strip_prefix("chr").unwrap_or(chrom);
    if let Ok(number) = bare.parse::<u8>()
        && (1..=22).contains(&number)
    {
        return format!("chr{number}");
    }
    match bare {
        "X" | "Y" => format!("chr{bare}"),
        "M" | "MT" => "chrMT".to_string(),
        _ => chrom.to_string(),
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use tempfile::TempDir;

    #[test]
    fn encodes_labels_without_losing_original_chrom() {
        assert_eq!(dataset_dir_name("chr1"), "chr1.lance");
        assert_eq!(dataset_dir_name("1"), "chr1.lance");
        assert_eq!(canonical_chrom_label("1"), "chr1");
        assert_eq!(canonical_chrom_label("chr1"), "chr1");
        assert_eq!(canonical_chrom_label("X"), "chrX");
        assert_eq!(canonical_chrom_label("MT"), "chrMT");
        assert_eq!(
            dataset_dir_name("HSCHR6_MHC_COX_CTG1"),
            "HSCHR6_MHC_COX_CTG1.lance"
        );
        assert_eq!(dataset_dir_name("foo/bar"), "foo_2Fbar.lance");
    }

    #[test]
    fn manifest_round_trips_chrom_paths() {
        let dir = TempDir::new().unwrap();
        let manifest = ChromManifest::new(vec![
            ChromDatasetEntry::new("chr1", "chr1.lance", 10),
            ChromDatasetEntry::new("HSCHR6_MHC_COX_CTG1", "HSCHR6_MHC_COX_CTG1.lance", 20),
        ]);

        manifest.write_to_entity_dir(dir.path()).unwrap();
        let loaded = ChromManifest::read_from_entity_dir(dir.path()).unwrap();

        assert_eq!(loaded.path_for_chrom("chr1").unwrap(), "chr1.lance");
        assert_eq!(
            loaded.path_for_chrom("HSCHR6_MHC_COX_CTG1").unwrap(),
            "HSCHR6_MHC_COX_CTG1.lance"
        );
        assert!(loaded.path_for_chrom("other").is_none());
        assert_eq!(loaded.available_chroms(), ["chr1", "HSCHR6_MHC_COX_CTG1"]);
    }

    #[test]
    fn path_for_chrom_tolerates_chr_prefix_spelling() {
        // Manifest canonicalized to `chr`-prefixed labels.
        let chr_prefixed =
            ChromManifest::new(vec![ChromDatasetEntry::new("chr1", "chr1.lance", 1)]);
        assert_eq!(chr_prefixed.path_for_chrom("chr1").unwrap(), "chr1.lance");
        // A bare VCF contig spelling still resolves.
        assert_eq!(chr_prefixed.path_for_chrom("1").unwrap(), "chr1.lance");
        assert!(chr_prefixed.path_for_chrom("2").is_none());

        // Manifest with bare labels resolves a `chr`-prefixed query too.
        let bare = ChromManifest::new(vec![ChromDatasetEntry::new("1", "1.lance", 1)]);
        assert_eq!(bare.path_for_chrom("1").unwrap(), "1.lance");
        assert_eq!(bare.path_for_chrom("chr1").unwrap(), "1.lance");
        assert!(bare.path_for_chrom("chrX").is_none());
    }
}
