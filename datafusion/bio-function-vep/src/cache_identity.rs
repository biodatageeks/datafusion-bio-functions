//! Lazy, shard-local cache identity validation.

use std::collections::HashMap;
use std::path::{Path, PathBuf};
use std::sync::Mutex;

use datafusion::arrow::datatypes::Schema;
use datafusion::common::{DataFusionError, Result};

use crate::cache::manifest::canonical_chrom_label;
use crate::cache_source::{CACHE_SOURCE_METADATA_KEY, CacheSourceType};
use crate::parquet_cache::detect::PartitionedParquetCache;
use crate::vep_semantics::{SupportedVepTarget, target_for_cache_version};

// Keep the release contract independent of whether the optional raw-cache
// builder dependency has already published the same metadata-key constant.
pub(crate) const CACHE_VERSION_METADATA_KEY: &str = "bio.vep.cache_version";

const CACHE_ENTITIES: &[&str] = &[
    "variation",
    "transcript",
    "exon",
    "translation_core",
    "translation_sift",
    "regulatory",
    "motif",
];

#[derive(Debug, Clone, PartialEq, Eq)]
pub(crate) struct CacheIdentity {
    pub(crate) source_type: CacheSourceType,
    pub(crate) cache_version: String,
    pub(crate) target: &'static SupportedVepTarget,
}

#[derive(Debug, Clone, PartialEq, Eq)]
pub struct ValidatedCacheIdentity {
    pub cache_source_type: String,
    pub cache_version: String,
    pub target: &'static SupportedVepTarget,
}

#[derive(Debug, Default)]
struct ValidationState {
    validated_contigs: HashMap<String, CacheIdentity>,
    invocation_identity: Option<CacheIdentity>,
}

/// Validate only the Parquet shards participating in one requested contig.
///
/// This is the reporting/diagnostic façade over the same validator used by the
/// annotation runtime. It intentionally reads manifests plus the selected
/// contig's footers and never probes an unrelated contig.
pub async fn validate_partitioned_cache_contig(
    cache_source: &str,
    chrom: &str,
    expected_cache_version: Option<String>,
) -> Result<ValidatedCacheIdentity> {
    let cache = PartitionedParquetCache::detect(cache_source).ok_or_else(|| {
        DataFusionError::Execution(format!(
            "no partitioned Parquet VEP cache found at '{cache_source}'"
        ))
    })?;
    let identity = LazyCacheIdentityValidator::new(expected_cache_version)?
        .validate_contig(&cache, chrom)
        .await?;
    Ok(ValidatedCacheIdentity {
        cache_source_type: identity.source_type.as_str().to_string(),
        cache_version: identity.cache_version,
        target: identity.target,
    })
}

#[derive(Debug)]
pub(crate) struct LazyCacheIdentityValidator {
    expected_cache_version: Option<String>,
    state: Mutex<ValidationState>,
}

impl LazyCacheIdentityValidator {
    pub(crate) fn new(expected_cache_version: Option<String>) -> Result<Self> {
        if let Some(expected) = expected_cache_version.as_deref() {
            target_for_cache_version(expected).map_err(|error| {
                DataFusionError::Plan(format!(
                    "annotate_vep(): invalid expected_cache_version assertion: {error}"
                ))
            })?;
        }
        Ok(Self {
            expected_cache_version,
            state: Mutex::new(ValidationState::default()),
        })
    }

    pub(crate) async fn validate_contig(
        &self,
        cache: &PartitionedParquetCache,
        chrom: &str,
    ) -> Result<CacheIdentity> {
        let canonical_chrom = canonical_chrom_label(chrom);
        if let Some(identity) = self
            .state
            .lock()
            .map_err(lock_error)?
            .validated_contigs
            .get(&canonical_chrom)
            .cloned()
        {
            return Ok(identity);
        }

        let shards = participating_shards(cache, chrom)?;
        let mut first: Option<(String, PathBuf, CacheSourceType, String)> = None;
        for (entity, path) in shards {
            let schema = read_parquet_shard_schema(&path).await?;
            let source_type = CacheSourceType::from_schema(&schema).map_err(|error| {
                DataFusionError::Execution(format!(
                    "annotate_vep(): cache identity validation failed for contig \
                     '{chrom}', entity '{entity}', shard '{}': {error}",
                    path.display()
                ))
            })?;
            let cache_version = schema
                .metadata()
                .get(CACHE_VERSION_METADATA_KEY)
                .ok_or_else(|| {
                    DataFusionError::Execution(format!(
                        "annotate_vep(): cache identity validation failed for contig \
                         '{chrom}', entity '{entity}', shard '{}': missing \
                         {CACHE_VERSION_METADATA_KEY}; rebuild this metadata-less cache",
                        path.display()
                    ))
                })?
                .to_string();
            target_for_cache_version(&cache_version).map_err(|error| {
                DataFusionError::Execution(format!(
                    "annotate_vep(): cache identity validation failed for contig \
                     '{chrom}', entity '{entity}', shard '{}': {error}",
                    path.display()
                ))
            })?;

            if let Some((first_entity, first_path, first_source, first_version)) = &first {
                if source_type != *first_source || cache_version != *first_version {
                    return Err(DataFusionError::Execution(format!(
                        "annotate_vep(): mixed cache identity for contig '{chrom}': \
                         {first_entity} shard '{}' declares source={} version={}, but \
                         {entity} shard '{}' declares source={} version={}",
                        first_path.display(),
                        first_source.as_str(),
                        first_version,
                        path.display(),
                        source_type.as_str(),
                        cache_version
                    )));
                }
            } else {
                first = Some((entity, path, source_type, cache_version));
            }
        }

        let (_, _, source_type, cache_version) = first.ok_or_else(|| {
            DataFusionError::Execution(format!(
                "annotate_vep(): no participating Parquet cache shards found for contig '{chrom}'"
            ))
        })?;
        if let Some(expected) = self.expected_cache_version.as_deref()
            && cache_version != expected
        {
            return Err(DataFusionError::Execution(format!(
                "annotate_vep(): cache version assertion failed for contig '{chrom}': \
                 embedded {CACHE_VERSION_METADATA_KEY}={cache_version}, expected {expected}"
            )));
        }

        let identity = CacheIdentity {
            source_type,
            target: target_for_cache_version(&cache_version)?,
            cache_version,
        };
        let mut state = self.state.lock().map_err(lock_error)?;
        if let Some(invocation) = &state.invocation_identity {
            if invocation.source_type != identity.source_type
                || invocation.cache_version != identity.cache_version
            {
                return Err(DataFusionError::Execution(format!(
                    "annotate_vep(): cache identity changed within one invocation at \
                     contig '{chrom}': first contig used source={} version={}, this \
                     contig uses source={} version={}",
                    invocation.source_type.as_str(),
                    invocation.cache_version,
                    identity.source_type.as_str(),
                    identity.cache_version
                )));
            }
        } else {
            state.invocation_identity = Some(identity.clone());
        }
        state
            .validated_contigs
            .insert(canonical_chrom, identity.clone());
        Ok(identity)
    }
}

fn lock_error<T>(error: std::sync::PoisonError<T>) -> DataFusionError {
    DataFusionError::Execution(format!("cache identity validator lock poisoned: {error}"))
}

fn participating_shards(
    cache: &PartitionedParquetCache,
    chrom: &str,
) -> Result<Vec<(String, PathBuf)>> {
    let variation = cache.variation_path(chrom).ok_or_else(|| {
        DataFusionError::Execution(format!(
            "annotate_vep(): Parquet cache has no variation shard for contig '{chrom}'"
        ))
    })?;
    let mut shards = vec![("variation".to_string(), variation)];
    for entity in CACHE_ENTITIES.iter().copied().skip(1) {
        if let Some(path) = cache.context_path(entity, chrom) {
            shards.push((entity.to_string(), path));
        }
    }
    Ok(shards)
}

async fn read_parquet_shard_schema(path: &Path) -> Result<Schema> {
    use parquet::arrow::async_reader::ParquetRecordBatchStreamBuilder;

    let file = tokio::fs::File::open(path).await.map_err(|error| {
        DataFusionError::Execution(format!(
            "annotate_vep(): failed to open Parquet cache shard footer '{}': {error}",
            path.display()
        ))
    })?;
    let builder = ParquetRecordBatchStreamBuilder::new(file)
        .await
        .map_err(|error| {
            DataFusionError::Execution(format!(
                "annotate_vep(): failed to read Parquet cache shard footer '{}': {error}",
                path.display()
            ))
        })?;
    Ok(builder.schema().as_ref().clone())
}

#[cfg(test)]
mod tests {
    use std::collections::HashMap;
    use std::sync::Arc;

    use datafusion::arrow::array::{Int64Array, RecordBatch};
    use datafusion::arrow::datatypes::{DataType, Field};
    use parquet::arrow::ArrowWriter;

    use super::*;
    use crate::cache::manifest::{ChromDatasetEntry, ChromManifest};

    fn write_shard(
        root: &Path,
        entity: &str,
        chrom: &str,
        source_type: Option<&str>,
        cache_version: Option<&str>,
    ) {
        let entity_dir = root.join(entity);
        std::fs::create_dir_all(&entity_dir).unwrap();
        let file_name = format!("{chrom}.parquet");
        let path = entity_dir.join(&file_name);
        let mut metadata = HashMap::new();
        if let Some(source_type) = source_type {
            metadata.insert(
                CACHE_SOURCE_METADATA_KEY.to_string(),
                source_type.to_string(),
            );
        }
        if let Some(cache_version) = cache_version {
            metadata.insert(
                CACHE_VERSION_METADATA_KEY.to_string(),
                cache_version.to_string(),
            );
        }
        let schema = Arc::new(
            Schema::new(vec![Field::new("value", DataType::Int64, false)]).with_metadata(metadata),
        );
        let batch = RecordBatch::try_new(schema.clone(), vec![Arc::new(Int64Array::from(vec![1]))])
            .unwrap();
        let mut writer =
            ArrowWriter::try_new(std::fs::File::create(&path).unwrap(), schema, None).unwrap();
        writer.write(&batch).unwrap();
        writer.close().unwrap();
        ChromManifest::new(vec![ChromDatasetEntry::new(chrom, file_name, 1)])
            .write_to_entity_dir(&entity_dir)
            .unwrap();
    }

    fn write_two_chrom_manifest(root: &Path, entity: &str) {
        ChromManifest::new(vec![
            ChromDatasetEntry::new("chr1", "chr1.parquet", 1),
            ChromDatasetEntry::new("chr2", "chr2.parquet", 1),
        ])
        .write_to_entity_dir(&root.join(entity))
        .unwrap();
    }

    #[tokio::test]
    async fn validates_supported_115_and_116_identities() {
        for version in ["115", "116"] {
            let tmp = tempfile::tempdir().unwrap();
            write_shard(
                tmp.path(),
                "variation",
                "chr1",
                Some("merged"),
                Some(version),
            );
            write_shard(
                tmp.path(),
                "transcript",
                "chr1",
                Some("merged"),
                Some(version),
            );
            let cache = PartitionedParquetCache::detect(tmp.path().to_str().unwrap()).unwrap();
            let identity = LazyCacheIdentityValidator::new(Some(version.to_string()))
                .unwrap()
                .validate_contig(&cache, "chr1")
                .await
                .unwrap();
            assert_eq!(identity.cache_version, version);
            assert_eq!(identity.source_type, CacheSourceType::Merged);
        }
    }

    #[tokio::test]
    async fn rejects_missing_and_mixed_metadata_with_shard_diagnostics() {
        let tmp = tempfile::tempdir().unwrap();
        write_shard(tmp.path(), "variation", "chr1", Some("merged"), Some("116"));
        write_shard(tmp.path(), "transcript", "chr1", Some("merged"), None);
        let cache = PartitionedParquetCache::detect(tmp.path().to_str().unwrap()).unwrap();
        let error = LazyCacheIdentityValidator::new(None)
            .unwrap()
            .validate_contig(&cache, "chr1")
            .await
            .unwrap_err()
            .to_string();
        assert!(error.contains("transcript"));
        assert!(error.contains(CACHE_VERSION_METADATA_KEY));

        write_shard(
            tmp.path(),
            "transcript",
            "chr1",
            Some("merged"),
            Some("115"),
        );
        let error = LazyCacheIdentityValidator::new(None)
            .unwrap()
            .validate_contig(&cache, "chr1")
            .await
            .unwrap_err()
            .to_string();
        assert!(error.contains("mixed cache identity"));
        assert!(error.contains("version=116"));
        assert!(error.contains("version=115"));
    }

    #[tokio::test]
    async fn chr1_validation_does_not_open_broken_chr2_footer() {
        let tmp = tempfile::tempdir().unwrap();
        write_shard(tmp.path(), "variation", "chr1", Some("merged"), Some("116"));
        std::fs::write(
            tmp.path().join("variation").join("chr2.parquet"),
            b"not parquet",
        )
        .unwrap();
        write_two_chrom_manifest(tmp.path(), "variation");
        let cache = PartitionedParquetCache::detect(tmp.path().to_str().unwrap()).unwrap();
        let validator = LazyCacheIdentityValidator::new(None).unwrap();

        validator.validate_contig(&cache, "chr1").await.unwrap();
        let error = validator
            .validate_contig(&cache, "chr2")
            .await
            .unwrap_err()
            .to_string();
        assert!(error.contains("chr2.parquet"));
    }

    #[tokio::test]
    async fn assertion_cannot_replace_missing_or_conflicting_identity() {
        let tmp = tempfile::tempdir().unwrap();
        write_shard(tmp.path(), "variation", "chr1", Some("ensembl"), None);
        let cache = PartitionedParquetCache::detect(tmp.path().to_str().unwrap()).unwrap();
        let error = LazyCacheIdentityValidator::new(Some("115".to_string()))
            .unwrap()
            .validate_contig(&cache, "chr1")
            .await
            .unwrap_err()
            .to_string();
        assert!(error.contains("missing bio.vep.cache_version"));

        write_shard(
            tmp.path(),
            "variation",
            "chr1",
            Some("ensembl"),
            Some("116"),
        );
        let error = LazyCacheIdentityValidator::new(Some("115".to_string()))
            .unwrap()
            .validate_contig(&cache, "chr1")
            .await
            .unwrap_err()
            .to_string();
        assert!(error.contains("embedded bio.vep.cache_version=116, expected 115"));
    }

    #[tokio::test]
    async fn later_contig_must_match_invocation_identity() {
        let tmp = tempfile::tempdir().unwrap();
        write_shard(tmp.path(), "variation", "chr1", Some("merged"), Some("115"));
        write_shard(tmp.path(), "variation", "chr2", Some("merged"), Some("116"));
        write_two_chrom_manifest(tmp.path(), "variation");
        let cache = PartitionedParquetCache::detect(tmp.path().to_str().unwrap()).unwrap();
        let validator = LazyCacheIdentityValidator::new(None).unwrap();
        validator.validate_contig(&cache, "chr1").await.unwrap();
        let error = validator
            .validate_contig(&cache, "chr2")
            .await
            .unwrap_err()
            .to_string();
        assert!(error.contains("changed within one invocation"));
    }
}
