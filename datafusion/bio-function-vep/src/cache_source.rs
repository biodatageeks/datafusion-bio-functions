//! VEP cache source-mode metadata helpers.

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

    #[cfg(feature = "lance-cache")]
    pub(crate) fn from_partitioned_lance_cache_source(cache_source: &str) -> Result<Self> {
        let dataset = first_variation_lance_dataset(cache_source)?;
        Self::from_lance_dataset(&dataset)
    }

    #[cfg(feature = "lance-cache")]
    pub(crate) fn from_lance_dataset(dataset_path: &Path) -> Result<Self> {
        let schema = read_lance_dataset_schema_sync(dataset_path)?;
        Self::from_schema(&schema)
    }
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

#[cfg(feature = "lance-cache")]
fn first_variation_lance_dataset(cache_source: &str) -> Result<PathBuf> {
    let variation_dir = Path::new(cache_source).join("variation.lance");
    let manifest =
        crate::lance_cache::manifest::ChromManifest::read_from_entity_dir(&variation_dir)?;
    let first = manifest.entries.first().ok_or_else(|| {
        DataFusionError::Plan(format!(
            "annotate_vep(): Lance cache source '{}' variation.lance manifest contains no chromosomes",
            cache_source
        ))
    })?;
    Ok(variation_dir.join(&first.dataset))
}

#[cfg(feature = "lance-cache")]
fn read_lance_dataset_schema_sync(dataset_path: &Path) -> Result<Schema> {
    let path = dataset_path.to_path_buf();
    let open = async move {
        lance::Dataset::open(path.to_string_lossy().as_ref())
            .await
            .map_err(|err| {
                DataFusionError::Execution(format!(
                    "annotate_vep(): failed to open Lance cache dataset '{}': {err}",
                    path.display()
                ))
            })
            .map(|dataset| dataset.schema().into())
    };

    match tokio::runtime::Handle::try_current() {
        // Multi-thread runtime: `block_in_place` is supported and cheapest.
        Ok(handle) if handle.runtime_flavor() == tokio::runtime::RuntimeFlavor::MultiThread => {
            tokio::task::block_in_place(|| handle.block_on(open))
        }
        // Current-thread runtime (e.g. `#[tokio::test]`, embedded callers):
        // `block_in_place` would panic and we cannot nest a runtime on this
        // thread, so resolve on a dedicated OS thread with its own runtime.
        Ok(_handle) => std::thread::scope(|scope| {
            scope
                .spawn(|| {
                    let rt = tokio::runtime::Runtime::new()
                        .map_err(|err| DataFusionError::External(Box::new(err)))?;
                    rt.block_on(open)
                })
                .join()
                .unwrap_or_else(|_| {
                    Err(DataFusionError::Execution(
                        "Lance schema read worker thread panicked".to_string(),
                    ))
                })
        }),
        // No runtime in scope: create one here.
        Err(_) => {
            let rt = tokio::runtime::Runtime::new()
                .map_err(|err| DataFusionError::External(Box::new(err)))?;
            rt.block_on(open)
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
