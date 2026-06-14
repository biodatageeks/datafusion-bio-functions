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

    pub fn path_for_chrom(&self, chrom: &str) -> Option<&str> {
        self.by_chrom
            .get(chrom)
            .map(|idx| self.entries[*idx].dataset.as_str())
    }

    pub fn available_chroms(&self) -> Vec<&str> {
        self.entries
            .iter()
            .map(|entry| entry.chrom.as_str())
            .collect()
    }
}

pub fn dataset_dir_name(chrom: &str) -> String {
    let mut encoded = String::with_capacity(chrom.len() + ".lance".len());
    for byte in chrom.bytes() {
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

#[cfg(test)]
mod tests {
    use super::*;
    use tempfile::TempDir;

    #[test]
    fn encodes_labels_without_losing_original_chrom() {
        assert_eq!(dataset_dir_name("chr1"), "chr1.lance");
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
}
