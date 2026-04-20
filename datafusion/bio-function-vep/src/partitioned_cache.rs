//! Partitioned per-chromosome parquet cache detection and registration.
//!
//! The partitioned layout stores each context type in a subdirectory with
//! per-chromosome parquet files:
//!
//! ```text
//! 115_GRCh38_vep/
//!   variation/chr1.parquet
//!   variation/chr2.parquet
//!   transcript/chr1.parquet
//!   exon/chr1.parquet
//!   translation_core/chr1.parquet
//!   translation_sift/chr1.parquet
//!   regulatory/chr1.parquet
//!   motif/chr1.parquet
//! ```

use std::path::{Path, PathBuf};

use datafusion::common::{DataFusionError, Result};
use datafusion::prelude::{ParquetReadOptions, SessionContext};

/// Known context subdirectories in the partitioned cache layout.
pub const CONTEXT_TYPES: &[&str] = &[
    "variation",
    "transcript",
    "exon",
    "translation_core",
    "translation_sift",
    "regulatory",
    "motif",
];

/// Represents a partitioned per-chromosome parquet cache directory.
#[derive(Clone)]
pub struct PartitionedParquetCache {
    base_dir: PathBuf,
    /// Chromosomes available in the variation subdirectory (from filenames).
    available_chroms: Vec<String>,
}

impl PartitionedParquetCache {
    /// Detect a partitioned cache layout at `path`.
    ///
    /// Returns `Some` if `path` is a directory containing a `variation/`
    /// subdirectory with at least one `*.parquet` file.
    pub fn detect(path: &str) -> Option<Self> {
        let base = Path::new(path);
        if !base.is_dir() {
            return None;
        }
        let variation_dir = base.join("variation");
        if !variation_dir.is_dir() {
            return None;
        }

        let mut chroms = Vec::new();
        let entries = std::fs::read_dir(&variation_dir).ok()?;
        for entry in entries.flatten() {
            let name = entry.file_name();
            let name_str = name.to_string_lossy();
            if let Some(stem) = name_str.strip_suffix(".parquet") {
                chroms.push(stem.to_string());
            }
        }

        if chroms.is_empty() {
            return None;
        }

        // Sort chromosomes naturally: chr1, chr2, ..., chr10, ..., chrX, chrY, chrMT
        chroms.sort_by_key(|a| natural_chrom_order(a));

        Some(Self {
            base_dir: base.to_path_buf(),
            available_chroms: chroms,
        })
    }

    /// All chromosomes available in the variation cache.
    pub fn available_chroms(&self) -> &[String] {
        &self.available_chroms
    }

    /// Path to a per-chromosome parquet file for a given context type.
    /// Returns `None` if the file does not exist.
    pub fn context_path(&self, context_type: &str, chrom: &str) -> Option<PathBuf> {
        let context_dir = self.base_dir.join(context_type);
        for candidate in chrom_path_candidates(chrom) {
            let path = context_dir.join(format!("{candidate}.parquet"));
            if path.exists() {
                return Some(path);
            }
        }
        None
    }

    /// Whether a per-chromosome parquet file exists for a given context type.
    pub fn has_chrom(&self, context_type: &str, chrom: &str) -> bool {
        self.context_path(context_type, chrom).is_some()
    }

    /// Base directory of the cache.
    pub fn base_dir(&self) -> &Path {
        &self.base_dir
    }
}

/// Register a per-chromosome parquet file as an ephemeral DataFusion table.
///
/// Returns the generated table name if registration succeeded, or `None` if
/// the parquet file does not exist for this chrom/context_type combination.
pub async fn register_chrom_parquet(
    session: &SessionContext,
    cache: &PartitionedParquetCache,
    context_type: &str,
    chrom: &str,
) -> Result<Option<String>> {
    let path = match cache.context_path(context_type, chrom) {
        Some(p) => p,
        None => return Ok(None),
    };

    let table_name = ephemeral_table_name(context_type, chrom);

    // Skip if already registered
    if session.table(&table_name).await.is_ok() {
        return Ok(Some(table_name));
    }

    session
        .register_parquet(
            &table_name,
            path.to_str().ok_or_else(|| {
                DataFusionError::Execution(format!(
                    "non-UTF8 path for {context_type}/{chrom}.parquet"
                ))
            })?,
            ParquetReadOptions::default(),
        )
        .await?;

    Ok(Some(table_name))
}

/// Deregister an ephemeral table from the session.
pub async fn deregister_table(session: &SessionContext, name: &str) -> Result<()> {
    // deregister_table returns Option<Arc<dyn TableProvider>>; ignore it.
    let _ = session.deregister_table(name)?;
    Ok(())
}

/// Generate an ephemeral table name for a context_type + chrom combination.
fn ephemeral_table_name(context_type: &str, chrom: &str) -> String {
    format!(
        "__vep_partitioned_{context_type}_{}",
        sanitize_identifier_component(chrom)
    )
}

fn sanitize_identifier_component(value: &str) -> String {
    let mut out = String::with_capacity(value.len());
    for ch in value.chars() {
        if ch.is_ascii_alphanumeric() {
            out.push(ch.to_ascii_lowercase());
        } else {
            out.push('_');
        }
    }
    out
}

fn chrom_path_candidates(chrom: &str) -> Vec<String> {
    let mut candidates = Vec::with_capacity(3);
    candidates.push(chrom.to_string());
    if let Some(bare) = chrom.strip_prefix("chr") {
        candidates.push(bare.to_string());
        candidates.push(format!("chr{bare}"));
    } else {
        candidates.push(format!("chr{chrom}"));
    }
    candidates.dedup();
    candidates
}

/// Natural chromosome ordering: numeric chroms first (sorted numerically),
/// then X, Y, MT/M, then anything else alphabetically.
fn natural_chrom_order(chrom: &str) -> (u8, u32, String) {
    let bare = chrom.strip_prefix("chr").unwrap_or(chrom);
    if let Ok(n) = bare.parse::<u32>() {
        (0, n, String::new())
    } else {
        let priority = match bare {
            "X" => 1,
            "Y" => 2,
            "MT" | "M" => 3,
            _ => 4,
        };
        (priority, 0, bare.to_string())
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_natural_chrom_order() {
        let mut chroms = vec![
            "chrX".to_string(),
            "chr2".to_string(),
            "chr10".to_string(),
            "chr1".to_string(),
            "chrY".to_string(),
            "chrMT".to_string(),
            "chr22".to_string(),
        ];
        chroms.sort_by_key(|a| natural_chrom_order(a));
        assert_eq!(
            chroms,
            vec!["chr1", "chr2", "chr10", "chr22", "chrX", "chrY", "chrMT"]
        );
    }

    #[test]
    fn test_detect_nonexistent_dir() {
        assert!(PartitionedParquetCache::detect("/nonexistent/path").is_none());
    }

    #[test]
    fn test_detect_file_not_dir() {
        // A regular file should not be detected as partitioned cache
        assert!(PartitionedParquetCache::detect("/dev/null").is_none());
    }

    #[test]
    fn test_ephemeral_table_name() {
        assert_eq!(
            ephemeral_table_name("variation", "chr1"),
            "__vep_partitioned_variation_chr1"
        );
        assert_eq!(
            ephemeral_table_name("transcript", "chrY"),
            "__vep_partitioned_transcript_chry"
        );
    }

    #[test]
    fn test_detect_dir_without_variation() {
        let tmp = tempfile::tempdir().unwrap();
        // Directory exists but no variation/ subdirectory
        assert!(PartitionedParquetCache::detect(tmp.path().to_str().unwrap()).is_none());
    }

    #[test]
    fn test_detect_variation_dir_empty() {
        let tmp = tempfile::tempdir().unwrap();
        std::fs::create_dir(tmp.path().join("variation")).unwrap();
        // variation/ exists but no parquet files
        assert!(PartitionedParquetCache::detect(tmp.path().to_str().unwrap()).is_none());
    }

    #[test]
    fn test_detect_valid_layout() {
        let tmp = tempfile::tempdir().unwrap();
        let var_dir = tmp.path().join("variation");
        std::fs::create_dir(&var_dir).unwrap();
        // Create dummy parquet files (content doesn't matter for detection)
        std::fs::write(var_dir.join("chr1.parquet"), b"dummy").unwrap();
        std::fs::write(var_dir.join("chr22.parquet"), b"dummy").unwrap();
        std::fs::write(var_dir.join("chrX.parquet"), b"dummy").unwrap();

        let cache = PartitionedParquetCache::detect(tmp.path().to_str().unwrap()).unwrap();
        assert_eq!(cache.available_chroms(), &["chr1", "chr22", "chrX"]);
    }

    #[test]
    fn test_context_path() {
        let tmp = tempfile::tempdir().unwrap();
        let var_dir = tmp.path().join("variation");
        std::fs::create_dir(&var_dir).unwrap();
        std::fs::write(var_dir.join("chr1.parquet"), b"dummy").unwrap();

        let tx_dir = tmp.path().join("transcript");
        std::fs::create_dir(&tx_dir).unwrap();
        std::fs::write(tx_dir.join("chr1.parquet"), b"dummy").unwrap();

        let cache = PartitionedParquetCache::detect(tmp.path().to_str().unwrap()).unwrap();

        assert!(cache.has_chrom("variation", "chr1"));
        assert!(cache.has_chrom("variation", "1"));
        assert!(cache.has_chrom("transcript", "chr1"));
        assert!(cache.has_chrom("transcript", "1"));
        assert!(!cache.has_chrom("transcript", "chr99"));
        assert!(!cache.has_chrom("nonexistent", "chr1"));
    }
}
