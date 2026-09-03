//! VEP cache source-mode metadata helpers.

#[cfg(feature = "parquet-cache")]
use std::path::{Path, PathBuf};
use std::str::FromStr;

use datafusion::arrow::datatypes::Schema;
use datafusion::common::{DataFusionError, Result};

pub(crate) const CACHE_SOURCE_METADATA_KEY: &str = "bio.vep.cache_source_type";

#[derive(Debug, Clone, Copy, PartialEq, Eq, Default)]
pub(crate) enum CacheSourceType {
    #[default]
    Ensembl,
    Merged,
    RefSeq,
}

impl CacheSourceType {
    pub(crate) fn as_str(self) -> &'static str {
        match self {
            Self::Ensembl => "ensembl",
            Self::Merged => "merged",
            Self::RefSeq => "refseq",
        }
    }

    pub(crate) fn from_schema(schema: &Schema) -> Result<Self> {
        let raw = schema.metadata().get(CACHE_SOURCE_METADATA_KEY).ok_or_else(|| {
            DataFusionError::Plan(format!(
                "annotate_vep(): cache table schema is missing {CACHE_SOURCE_METADATA_KEY}; expected one of ensembl, merged, refseq"
            ))
        })?;
        raw.parse()
    }

    /// Detect the cache source type from a partitioned Parquet cache directory by
    /// reading the `bio.vep.cache_source_type` metadata off the schema of the
    /// first `variation` shard that the manifest lists *and* that exists on disk.
    ///
    /// The metadata is uniform across a cache, so any shard answers. Selecting
    /// by presence rather than by manifest position keeps a per-contig download
    /// usable: it ships the published whole-cache manifest (463 contigs, `chr1`
    /// first) alongside only the shards the user fetched.
    #[cfg(feature = "parquet-cache")]
    pub(crate) fn from_partitioned_parquet_cache_source(cache_source: &str) -> Result<Self> {
        let shard = first_present_variation_parquet_shard(cache_source)?;
        let schema = read_parquet_shard_schema_sync(&shard)?;
        Self::from_schema(&schema)
    }
}

/// The first `variation` shard listed in the manifest that is present on disk.
///
/// Requires the manifest (a cache without one is not a partitioned cache) and
/// at least one listed entry. A manifest whose shards are all absent is a
/// cache-level defect and is reported as such, in terms of the manifest,
/// rather than as an `open()` failure on whichever contig happens to be listed
/// first. A contig that *is* requested but absent still fails at the per-contig
/// identity check, naming that contig.
#[cfg(feature = "parquet-cache")]
fn first_present_variation_parquet_shard(cache_source: &str) -> Result<PathBuf> {
    use crate::cache::manifest::{CHROM_MANIFEST_FILE, ChromManifest};

    let variation_dir = Path::new(cache_source).join("variation");
    let manifest = ChromManifest::read_from_entity_dir(&variation_dir)?;
    if manifest.entries.is_empty() {
        return Err(DataFusionError::Plan(format!(
            "annotate_vep(): Parquet cache source '{cache_source}' variation manifest contains no chromosomes"
        )));
    }
    manifest
        .entries
        .iter()
        .map(|entry| variation_dir.join(&entry.dataset))
        .find(|path| path.is_file())
        .ok_or_else(|| {
            DataFusionError::Execution(format!(
                "annotate_vep(): none of the {} variation shards listed in '{}' is present on disk; \
                 download at least one contig's shards or check the cache directory",
                manifest.entries.len(),
                variation_dir.join(CHROM_MANIFEST_FILE).display()
            ))
        })
}

/// Read the Arrow schema (including key-value metadata) of a `.parquet` shard
/// synchronously from its footer.
#[cfg(feature = "parquet-cache")]
fn read_parquet_shard_schema_sync(shard_path: &Path) -> Result<Schema> {
    use parquet::arrow::arrow_reader::ParquetRecordBatchReaderBuilder;
    let file = std::fs::File::open(shard_path).map_err(|err| {
        DataFusionError::Execution(format!(
            "annotate_vep(): failed to open Parquet cache shard '{}': {err}",
            shard_path.display()
        ))
    })?;
    let builder = ParquetRecordBatchReaderBuilder::try_new(file).map_err(|err| {
        DataFusionError::Execution(format!(
            "annotate_vep(): failed to read Parquet cache shard schema '{}': {err}",
            shard_path.display()
        ))
    })?;
    Ok(builder.schema().as_ref().clone())
}

impl FromStr for CacheSourceType {
    type Err = DataFusionError;

    fn from_str(raw: &str) -> Result<Self> {
        match raw {
            "ensembl" => Ok(Self::Ensembl),
            "merged" => Ok(Self::Merged),
            "refseq" => Ok(Self::RefSeq),
            other => Err(DataFusionError::Plan(format!(
                "annotate_vep(): unsupported {CACHE_SOURCE_METADATA_KEY} '{other}'; expected one of ensembl, merged, refseq"
            ))),
        }
    }
}

#[cfg(test)]
mod tests {
    use std::collections::HashMap;

    use datafusion::arrow::datatypes::{DataType, Field};

    use super::*;

    fn schema_with_source(value: &str) -> Schema {
        let mut metadata = HashMap::new();
        metadata.insert(CACHE_SOURCE_METADATA_KEY.to_string(), value.to_string());
        Schema::new(vec![Field::new("chrom", DataType::Utf8, false)]).with_metadata(metadata)
    }

    #[test]
    fn parses_supported_cache_source_metadata() {
        assert_eq!(
            CacheSourceType::from_schema(&schema_with_source("ensembl")).unwrap(),
            CacheSourceType::Ensembl
        );
        assert_eq!(
            CacheSourceType::from_schema(&schema_with_source("merged")).unwrap(),
            CacheSourceType::Merged
        );
        assert_eq!(
            CacheSourceType::from_schema(&schema_with_source("refseq")).unwrap(),
            CacheSourceType::RefSeq
        );
    }

    #[test]
    fn rejects_missing_invalid_and_case_variant_source_metadata() {
        let missing = Schema::new(vec![Field::new("chrom", DataType::Utf8, false)]);
        let err = CacheSourceType::from_schema(&missing)
            .unwrap_err()
            .to_string();
        assert!(err.contains("missing bio.vep.cache_source_type"));

        let err = CacheSourceType::from_schema(&schema_with_source("RefSeq"))
            .unwrap_err()
            .to_string();
        assert!(err.contains("unsupported bio.vep.cache_source_type 'RefSeq'"));

        let err = CacheSourceType::from_schema(&schema_with_source("other"))
            .unwrap_err()
            .to_string();
        assert!(err.contains("expected one of ensembl, merged, refseq"));
    }
}

#[cfg(all(test, feature = "parquet-cache"))]
mod partitioned_cache_tests {
    use std::collections::HashMap;
    use std::path::Path;
    use std::sync::Arc;

    use datafusion::arrow::array::{Int64Array, RecordBatch};
    use datafusion::arrow::datatypes::{DataType, Field};
    use parquet::arrow::ArrowWriter;

    use super::*;
    use crate::cache::manifest::{ChromDatasetEntry, ChromManifest};

    fn write_variation_manifest(root: &Path, chroms: &[&str]) {
        let entries = chroms
            .iter()
            .map(|chrom| ChromDatasetEntry::new(*chrom, format!("{chrom}.parquet"), 1))
            .collect();
        ChromManifest::new(entries)
            .write_to_entity_dir(&root.join("variation"))
            .unwrap();
    }

    fn write_variation_shard(root: &Path, chrom: &str, source_type: &str) {
        let mut metadata = HashMap::new();
        metadata.insert(
            CACHE_SOURCE_METADATA_KEY.to_string(),
            source_type.to_string(),
        );
        let schema = Arc::new(
            Schema::new(vec![Field::new("value", DataType::Int64, false)]).with_metadata(metadata),
        );
        let batch = RecordBatch::try_new(schema.clone(), vec![Arc::new(Int64Array::from(vec![1]))])
            .unwrap();
        let path = root.join("variation").join(format!("{chrom}.parquet"));
        let mut writer =
            ArrowWriter::try_new(std::fs::File::create(path).unwrap(), schema, None).unwrap();
        writer.write(&batch).unwrap();
        writer.close().unwrap();
    }

    /// A per-contig download keeps the published whole-cache manifest, so the
    /// manifest's first entry (`chr1`) is usually not on disk. The source-type
    /// probe must read a shard that exists, not the first one listed.
    #[test]
    fn partial_cache_reads_source_type_from_first_shard_present_on_disk() {
        let tmp = tempfile::tempdir().unwrap();
        write_variation_manifest(tmp.path(), &["chr1", "chr2", "chr21", "chr22"]);
        write_variation_shard(tmp.path(), "chr21", "refseq");
        write_variation_shard(tmp.path(), "chr22", "refseq");

        let source_type =
            CacheSourceType::from_partitioned_parquet_cache_source(tmp.path().to_str().unwrap())
                .unwrap();
        assert_eq!(source_type, CacheSourceType::RefSeq);
    }

    #[test]
    fn complete_cache_still_reads_the_first_listed_shard() {
        let tmp = tempfile::tempdir().unwrap();
        write_variation_manifest(tmp.path(), &["chr1", "chr2"]);
        write_variation_shard(tmp.path(), "chr1", "merged");
        write_variation_shard(tmp.path(), "chr2", "merged");

        let source_type =
            CacheSourceType::from_partitioned_parquet_cache_source(tmp.path().to_str().unwrap())
                .unwrap();
        assert_eq!(source_type, CacheSourceType::Merged);
    }

    /// No listed shard on disk is a cache-level defect: say so in terms of the
    /// manifest rather than reporting an `open()` failure on one arbitrary contig.
    #[test]
    fn cache_with_no_listed_shard_on_disk_reports_the_manifest_not_a_contig() {
        let tmp = tempfile::tempdir().unwrap();
        write_variation_manifest(tmp.path(), &["chr1", "chr2", "chr3"]);

        let err =
            CacheSourceType::from_partitioned_parquet_cache_source(tmp.path().to_str().unwrap())
                .unwrap_err()
                .to_string();
        assert!(
            err.contains("none of the 3 variation shards listed in"),
            "unexpected error: {err}"
        );
        assert!(
            err.contains("chrom_manifest.json"),
            "unexpected error: {err}"
        );
        assert!(
            !err.contains("No such file or directory"),
            "must not surface an open() error for an arbitrary contig: {err}"
        );
    }

    #[test]
    fn empty_manifest_is_still_an_error() {
        let tmp = tempfile::tempdir().unwrap();
        write_variation_manifest(tmp.path(), &[]);

        let err =
            CacheSourceType::from_partitioned_parquet_cache_source(tmp.path().to_str().unwrap())
                .unwrap_err()
                .to_string();
        assert!(
            err.contains("contains no chromosomes"),
            "unexpected error: {err}"
        );
    }

    /// The dataset card's broken recipe yields a shard and no manifest. That must
    /// keep failing on the manifest, unchanged by the shard-selection fix.
    #[test]
    fn missing_manifest_is_still_an_error() {
        let tmp = tempfile::tempdir().unwrap();
        std::fs::create_dir_all(tmp.path().join("variation")).unwrap();
        write_variation_shard(tmp.path(), "chr21", "ensembl");

        let err =
            CacheSourceType::from_partitioned_parquet_cache_source(tmp.path().to_str().unwrap())
                .unwrap_err()
                .to_string();
        assert!(err.contains("failed to read"), "unexpected error: {err}");
        assert!(
            err.contains("chrom_manifest.json"),
            "unexpected error: {err}"
        );
    }
}
