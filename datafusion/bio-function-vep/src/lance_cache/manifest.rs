use std::collections::{HashMap, HashSet};
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

    /// Write these entries into `entity_dir`, MERGING with any existing manifest
    /// rather than replacing it: an entry whose `chrom` matches one already on
    /// disk overwrites it, other pre-existing entries are preserved. This keeps a
    /// filtered / single-chromosome rebuild (`--chrom-filter`, `--overwrite`)
    /// from dropping the other chromosomes' datasets — which remain on disk — out
    /// of the manifest and thus making their context unresolvable.
    pub fn merge_write_to_entity_dir(&self, entity_dir: &Path) -> Result<()> {
        let mut merged = if entity_dir.join(CHROM_MANIFEST_FILE).exists() {
            Self::read_from_entity_dir(entity_dir)?.entries
        } else {
            Vec::new()
        };
        for entry in &self.entries {
            match merged.iter_mut().find(|e| e.chrom == entry.chrom) {
                Some(slot) => *slot = entry.clone(),
                None => merged.push(entry.clone()),
            }
        }
        Self::new(merged).write_to_entity_dir(entity_dir)
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

    /// Resolve the per-chromosome dataset path, tolerating contig-spelling
    /// differences between the query and the manifest (e.g. VCF contig `1` vs a
    /// manifest canonicalized to `chr1`, or mitochondrial `M`/`chrM` vs `chrMT`).
    /// Exact match wins; the [`contig_alias_set`] spellings are only tried as
    /// fallbacks. Mirrors the fallback previously applied only at the
    /// variation-lookup call site, so context entities
    /// (transcript/exon/translation/regulatory/motif/SIFT) resolve the same way
    /// instead of silently loading empty.
    pub fn path_for_chrom(&self, chrom: &str) -> Option<&str> {
        if let Some(idx) = self.by_chrom.get(chrom) {
            return Some(self.entries[*idx].dataset.as_str());
        }
        contig_alias_set(chrom)
            .iter()
            .find_map(|alias| self.by_chrom.get(alias))
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

/// All contig spellings equivalent to `chrom` under our canonicalization: the
/// input itself, its `chr`-prefix add/strip counterpart, and — for the
/// mitochondrion — the full `M`/`MT`/`chrM`/`chrMT` set. A plain chr-prefix
/// add/strip cannot bridge `M` <-> `MT`, yet [`canonical_chrom_label`] folds all
/// four into `chrMT`, so a VCF using `M`/`chrM` must still match a manifest
/// canonicalized to `chrMT` (and vice versa). Used wherever manifest contigs are
/// expanded (contig filtering) or resolved (path lookup).
pub fn contig_alias_set(chrom: &str) -> HashSet<String> {
    let mut set = HashSet::new();
    set.insert(chrom.to_string());
    if let Some(bare) = chrom.strip_prefix("chr") {
        set.insert(bare.to_string());
    } else {
        set.insert(format!("chr{chrom}"));
    }
    if canonical_chrom_label(chrom) == "chrMT" {
        for alias in ["M", "MT", "chrM", "chrMT"] {
            set.insert(alias.to_string());
        }
    }
    set
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

    #[test]
    fn merge_write_preserves_other_chromosomes() {
        let dir = TempDir::new().unwrap();
        // Existing multi-chromosome manifest.
        ChromManifest::new(vec![
            ChromDatasetEntry::new("chr1", "chr1.lance", 10),
            ChromDatasetEntry::new("chr2", "chr2.lance", 20),
        ])
        .write_to_entity_dir(dir.path())
        .unwrap();

        // Filtered rebuild of just chr1 (new row count) must keep chr2.
        ChromManifest::new(vec![ChromDatasetEntry::new("chr1", "chr1.lance", 11)])
            .merge_write_to_entity_dir(dir.path())
            .unwrap();

        let merged = ChromManifest::read_from_entity_dir(dir.path()).unwrap();
        assert_eq!(merged.path_for_chrom("chr1").unwrap(), "chr1.lance");
        assert_eq!(merged.path_for_chrom("chr2").unwrap(), "chr2.lance");
        let chr1 = merged.entries.iter().find(|e| e.chrom == "chr1").unwrap();
        assert_eq!(chr1.rows, 11, "rebuilt chr1 entry overwrites the old one");
        assert_eq!(merged.entries.len(), 2);
    }

    #[test]
    fn merge_write_into_empty_dir_writes_subset() {
        let dir = TempDir::new().unwrap();
        ChromManifest::new(vec![ChromDatasetEntry::new("chr1", "chr1.lance", 1)])
            .merge_write_to_entity_dir(dir.path())
            .unwrap();
        let loaded = ChromManifest::read_from_entity_dir(dir.path()).unwrap();
        assert_eq!(loaded.available_chroms(), ["chr1"]);
    }

    #[test]
    fn path_for_chrom_resolves_mitochondrial_aliases() {
        // Builder canonicalizes the mitochondrion to `chrMT`; every VCF spelling
        // must still resolve to it.
        let mito = ChromManifest::new(vec![ChromDatasetEntry::new("chrMT", "chrMT.lance", 1)]);
        for query in ["chrMT", "MT", "chrM", "M"] {
            assert_eq!(
                mito.path_for_chrom(query).unwrap(),
                "chrMT.lance",
                "mito query {query} must resolve to chrMT"
            );
        }
        assert!(mito.path_for_chrom("chr1").is_none());
    }

    #[test]
    fn contig_alias_set_covers_mito_and_chr_prefix() {
        assert!(contig_alias_set("M").contains("chrMT"));
        assert!(contig_alias_set("chrMT").contains("M"));
        assert!(contig_alias_set("chrM").contains("MT"));
        assert!(contig_alias_set("1").contains("chr1"));
        assert!(contig_alias_set("chr1").contains("1"));
        // Non-mito contigs don't pick up mito aliases.
        assert!(!contig_alias_set("chr1").contains("chrMT"));
    }
}
