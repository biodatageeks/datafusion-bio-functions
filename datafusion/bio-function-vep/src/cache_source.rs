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
    /// reading the `bio.vep.cache_source_type` metadata off the first
    /// `parquet.variation` shard's schema.
    #[cfg(feature = "parquet-cache")]
    pub(crate) fn from_partitioned_parquet_cache_source(cache_source: &str) -> Result<Self> {
        let shard = first_variation_parquet_shard(cache_source)?;
        let schema = read_parquet_shard_schema_sync(&shard)?;
        Self::from_schema(&schema)
    }
}

#[cfg(feature = "parquet-cache")]
fn first_variation_parquet_shard(cache_source: &str) -> Result<PathBuf> {
    let variation_dir = Path::new(cache_source).join("parquet.variation");
    let manifest =
        crate::lance_cache::manifest::ChromManifest::read_from_entity_dir(&variation_dir)?;
    let first = manifest.entries.first().ok_or_else(|| {
        DataFusionError::Plan(format!(
            "annotate_vep(): Parquet cache source '{cache_source}' parquet.variation manifest contains no chromosomes"
        ))
    })?;
    Ok(variation_dir.join(&first.dataset))
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
