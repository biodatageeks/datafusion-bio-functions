//! Cache-wide transcript-reference policy, qualified by the native cache identity.

use std::path::Path;

use datafusion::arrow::datatypes::Schema;
use datafusion::common::{DataFusionError, Result};
use serde::{Deserialize, Serialize};

use crate::cache_identity::CACHE_VERSION_METADATA_KEY;
use crate::cache_source::CacheSourceType;

pub(crate) const BAM_EDITED_METADATA_KEY: &str = "bio.vep.cache_bam_edited";
pub(crate) const REFERENCE_POLICY_FILE: &str = "reference_policy.json";

/// Missing policy in an older export remains UNKNOWN. Neither a directory name
/// nor a merged/RefSeq source declaration proves the native `info.txt` BAM value.
#[derive(Debug, Clone, PartialEq, Eq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub(crate) struct CacheReferencePolicy {
    schema_version: u32,
    cache_source_type: String,
    cache_version: String,
    bam_edited: bool,
}

pub(crate) fn bam_edited_from_schema(schema: &Schema) -> Result<Option<bool>> {
    schema
        .metadata()
        .get(BAM_EDITED_METADATA_KEY)
        .map(|value| match value.as_str() {
            "true" => Ok(true),
            "false" => Ok(false),
            _ => Err(DataFusionError::Execution(format!(
                "invalid {BAM_EDITED_METADATA_KEY} '{value}': expected true or false"
            ))),
        })
        .transpose()
}

impl CacheReferencePolicy {
    pub(crate) fn read(root: &Path) -> Result<Option<Self>> {
        let path = root.join(REFERENCE_POLICY_FILE);
        let bytes = match std::fs::read(&path) {
            Ok(bytes) => bytes,
            Err(error) if error.kind() == std::io::ErrorKind::NotFound => return Ok(None),
            Err(error) => {
                return Err(DataFusionError::Execution(format!(
                    "failed to read cache reference policy '{}': {error}",
                    path.display()
                )));
            }
        };
        let policy: Self = serde_json::from_slice(&bytes).map_err(|error| {
            DataFusionError::Execution(format!(
                "invalid cache reference policy '{}': {error}",
                path.display()
            ))
        })?;
        if policy.schema_version != 1
            || policy.cache_source_type.parse::<CacheSourceType>().is_err()
            || policy.cache_version.is_empty()
            || !policy
                .cache_version
                .bytes()
                .all(|byte| byte.is_ascii_digit())
        {
            return Err(DataFusionError::Execution(format!(
                "invalid cache reference policy '{}': expected schema_version=1, a cache source type and decimal cache_version",
                path.display()
            )));
        }
        Ok(Some(policy))
    }

    /// Root metadata can fill an absent optional shard key, but cannot override
    /// a contradictory source, release or explicit BAM declaration.
    pub(crate) fn reconcile(&self, schema: &Schema) -> Result<bool> {
        let source = CacheSourceType::from_schema(schema)?;
        let version = schema.metadata().get(CACHE_VERSION_METADATA_KEY);
        if source.as_str() != self.cache_source_type
            || version.map(String::as_str) != Some(self.cache_version.as_str())
        {
            return Err(DataFusionError::Execution(format!(
                "cache reference policy identity source={} version={} conflicts with shard source={} version={version:?}",
                self.cache_source_type,
                self.cache_version,
                source.as_str()
            )));
        }
        if let Some(shard_policy) = bam_edited_from_schema(schema)?
            && shard_policy != self.bam_edited
        {
            return Err(DataFusionError::Execution(format!(
                "cache reference policy bam_edited={} conflicts with shard {BAM_EDITED_METADATA_KEY}={shard_policy}",
                self.bam_edited
            )));
        }
        Ok(self.bam_edited)
    }

    #[cfg(feature = "cache-builder")]
    pub(crate) fn preserve(schema: &Schema, native_root: &Path, output: &Path) -> Result<()> {
        use std::io::Write;

        let policy = Self {
            schema_version: 1,
            cache_source_type: CacheSourceType::from_schema(schema)?.as_str().to_owned(),
            cache_version: schema
                .metadata()
                .get(CACHE_VERSION_METADATA_KEY)
                .ok_or_else(|| {
                    DataFusionError::Execution(
                        "native cache schema is missing cache version".to_string(),
                    )
                })?
                .clone(),
            bam_edited: bam_edited_from_schema(schema)?.ok_or_else(|| {
                DataFusionError::Execution("native cache schema is missing BAM policy".to_string())
            })?,
        };
        let needs_write = match Self::read(output)? {
            Some(existing) if existing == policy => false,
            Some(existing) => {
                return Err(DataFusionError::Execution(format!(
                    "cannot resume cache '{}': existing reference policy {existing:?} conflicts with native cache '{}' policy {policy:?}; use the original native cache or rebuild into a new destination",
                    output.display(),
                    native_root.display()
                )));
            }
            None => true,
        };
        // A resume must not relabel an older destination from a different raw
        // cache. Validate shards even when the root policy already agrees: a
        // copied or replaced shard may contradict that declaration.
        for entity in [
            "variation",
            "transcript",
            "exon",
            "translation_core",
            "translation_sift",
            "regulatory",
            "motif",
        ] {
            let entries = match std::fs::read_dir(output.join(entity)) {
                Ok(entries) => entries,
                Err(error) if error.kind() == std::io::ErrorKind::NotFound => continue,
                Err(error) => return Err(error.into()),
            };
            for entry in entries {
                let path = entry?.path();
                if path.extension().and_then(|ext| ext.to_str()) != Some("parquet") {
                    continue;
                }
                let shard = crate::cache_source::read_parquet_shard_schema_sync(&path)?;
                policy.reconcile(&shard).map_err(|error| {
                    DataFusionError::Execution(format!(
                        "cannot refresh '{}' from native cache '{}': existing shard '{}': {error}",
                        output.join(REFERENCE_POLICY_FILE).display(),
                        native_root.display(),
                        path.display()
                    ))
                })?;
            }
        }
        if !needs_write {
            return Ok(());
        }
        std::fs::create_dir_all(output)?;
        let destination = output.join(REFERENCE_POLICY_FILE);
        let mut temporary = tempfile::NamedTempFile::new_in(output)?;
        let mut bytes = serde_json::to_vec_pretty(&policy).map_err(|error| {
            DataFusionError::Execution(format!(
                "failed to serialize cache reference policy: {error}"
            ))
        })?;
        bytes.push(b'\n');
        temporary.write_all(&bytes)?;
        // Match native metadata permissions; tempfile's Unix default is 0600.
        temporary
            .as_file()
            .set_permissions(std::fs::metadata(native_root.join("info.txt"))?.permissions())?;
        temporary.persist(&destination).map_err(|error| {
            DataFusionError::Execution(format!(
                "failed to preserve cache reference policy '{}': {error}",
                destination.display()
            ))
        })?;
        Ok(())
    }
}

pub(crate) fn resolved_bam_policy(
    schema: &Schema,
    policy: Option<&CacheReferencePolicy>,
    shard: &Path,
) -> Result<Option<bool>> {
    match policy {
        Some(policy) => policy.reconcile(schema).map(Some),
        None => bam_edited_from_schema(schema),
    }
    .map_err(|error| {
        DataFusionError::Execution(format!(
            "invalid reference policy for cache shard '{}': {error}",
            shard.display()
        ))
    })
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::cache::manifest::{ChromDatasetEntry, ChromManifest};
    use crate::cache_source::{CACHE_SOURCE_METADATA_KEY, CacheMetadata};
    use datafusion::arrow::datatypes::{DataType, Field};
    use std::sync::Arc;

    fn schema(source: &str, version: &str, bam: Option<bool>) -> Schema {
        let mut metadata = [
            (CACHE_SOURCE_METADATA_KEY.to_string(), source.to_string()),
            (CACHE_VERSION_METADATA_KEY.to_string(), version.to_string()),
        ]
        .into_iter()
        .collect::<std::collections::HashMap<_, _>>();
        if let Some(bam) = bam {
            metadata.insert(BAM_EDITED_METADATA_KEY.to_string(), bam.to_string());
        }
        Schema::new(vec![Field::new("value", DataType::Int64, false)]).with_metadata(metadata)
    }

    fn policy(bam: bool) -> CacheReferencePolicy {
        CacheReferencePolicy {
            schema_version: 1,
            cache_source_type: "merged".into(),
            cache_version: "116".into(),
            bam_edited: bam,
        }
    }

    fn write_policy(root: &Path, bam: bool) {
        std::fs::write(
            root.join(REFERENCE_POLICY_FILE),
            serde_json::to_vec(&policy(bam)).unwrap(),
        )
        .unwrap();
    }

    fn write_shard(root: &Path, entity: &str, chrom: &str, schema: Schema) {
        let dir = root.join(entity);
        std::fs::create_dir_all(&dir).unwrap();
        let name = format!("{chrom}.parquet");
        let writer = parquet::arrow::ArrowWriter::try_new(
            std::fs::File::create(dir.join(&name)).unwrap(),
            Arc::new(schema),
            None,
        )
        .unwrap();
        writer.close().unwrap();
        ChromManifest::new(vec![ChromDatasetEntry::new(chrom, name, 0)])
            .merge_write_to_entity_dir(&dir)
            .unwrap();
    }

    #[test]
    fn eager_reference_policy_distinguishes_false_true_and_unknown() {
        let root = tempfile::tempdir().unwrap();
        for bam in [None, Some(false), Some(true)] {
            write_shard(
                root.path(),
                "variation",
                "chr21",
                schema("merged", "116", bam),
            );
            let metadata =
                CacheMetadata::from_partitioned_cache(root.path().to_str().unwrap()).unwrap();
            assert_eq!(metadata.bam_edited, bam);
            assert_eq!(metadata.source_type, CacheSourceType::Merged);
        }
        write_shard(
            root.path(),
            "variation",
            "chr21",
            schema("merged", "116", None),
        );
        for bam in [false, true] {
            write_policy(root.path(), bam);
            assert_eq!(
                CacheMetadata::from_partitioned_cache(root.path().to_str().unwrap())
                    .unwrap()
                    .bam_edited,
                Some(bam)
            );
        }
    }

    #[test]
    fn root_reference_policy_cannot_override_shard_identity_or_explicit_policy() {
        let root = tempfile::tempdir().unwrap();
        write_policy(root.path(), true);
        for shard in [
            schema("ensembl", "116", None),
            schema("merged", "115", None),
            schema("merged", "116", Some(false)),
        ] {
            write_shard(root.path(), "variation", "chr21", shard);
            let error = CacheMetadata::from_partitioned_cache(root.path().to_str().unwrap())
                .unwrap_err()
                .to_string();
            assert!(error.contains("chr21.parquet"), "{error}");
            assert!(error.contains("conflicts"), "{error}");
        }
    }

    #[test]
    fn malformed_present_reference_policy_is_not_legacy_unknown() {
        let root = tempfile::tempdir().unwrap();
        assert_eq!(CacheReferencePolicy::read(root.path()).unwrap(), None);
        for text in [
            "",
            "{}",
            "null",
            r#"{"schema_version":2,"cache_source_type":"merged","cache_version":"116","bam_edited":true}"#,
            r#"{"schema_version":1,"cache_source_type":"merged","cache_version":"v116","bam_edited":true}"#,
            r#"{"schema_version":1,"cache_source_type":"merged","cache_version":"116","bam_edited":"false"}"#,
        ] {
            std::fs::write(root.path().join(REFERENCE_POLICY_FILE), text).unwrap();
            let error = CacheReferencePolicy::read(root.path())
                .unwrap_err()
                .to_string();
            assert!(error.contains(REFERENCE_POLICY_FILE), "{error}");
        }
    }

    #[tokio::test]
    async fn all_participating_shards_and_later_contigs_must_agree_on_reference_policy() {
        use crate::cache_identity::LazyCacheIdentityValidator;
        use crate::parquet_cache::detect::PartitionedParquetCache;
        let root = tempfile::tempdir().unwrap();
        write_policy(root.path(), true);
        write_shard(
            root.path(),
            "variation",
            "chr1",
            schema("merged", "116", None),
        );
        write_shard(
            root.path(),
            "transcript",
            "chr1",
            schema("merged", "116", Some(false)),
        );
        let cache = PartitionedParquetCache::try_detect(root.path().to_str().unwrap())
            .unwrap()
            .unwrap();
        let error = LazyCacheIdentityValidator::new(None)
            .unwrap()
            .validate_contig(&cache, "chr1")
            .await
            .unwrap_err()
            .to_string();
        assert!(error.contains("transcript/chr1.parquet"), "{error}");
        write_shard(
            root.path(),
            "transcript",
            "chr1",
            schema("merged", "116", None),
        );
        let validator = LazyCacheIdentityValidator::new(None).unwrap();
        assert_eq!(
            validator
                .validate_contig(&cache, "chr1")
                .await
                .unwrap()
                .bam_edited,
            Some(true)
        );

        std::fs::remove_file(root.path().join(REFERENCE_POLICY_FILE)).unwrap();
        write_shard(
            root.path(),
            "variation",
            "chr2",
            schema("merged", "116", Some(true)),
        );
        let cache = PartitionedParquetCache::try_detect(root.path().to_str().unwrap())
            .unwrap()
            .unwrap();
        let validator = LazyCacheIdentityValidator::new(None).unwrap();
        assert_eq!(
            validator
                .validate_contig(&cache, "chr1")
                .await
                .unwrap()
                .bam_edited,
            None
        );
        let error = validator
            .validate_contig(&cache, "chr2")
            .await
            .unwrap_err()
            .to_string();
        assert!(error.contains("changed within one invocation"), "{error}");
    }

    #[cfg(feature = "cache-builder")]
    #[test]
    fn reference_policy_resume_rejects_conflicting_root_without_modifying_files() {
        let root = tempfile::tempdir().unwrap();
        let native = root.path().join("native");
        let output = root.path().join("output");
        std::fs::create_dir(&native).unwrap();
        std::fs::write(native.join("info.txt"), b"bam\t/path.bam\n").unwrap();
        // Legacy shards cannot identify the BAM policy on their own. The saved
        // root declaration is authoritative and must not be silently replaced.
        write_shard(
            &output,
            "transcript",
            "chr21",
            schema("merged", "116", None),
        );
        let shard = output.join("transcript/chr21.parquet");
        let original_shard = std::fs::read(&shard).unwrap();
        for incoming_bam in [false, true] {
            let native_schema = schema("merged", "116", Some(incoming_bam));
            for existing in [
                policy(!incoming_bam),
                CacheReferencePolicy {
                    cache_source_type: "ensembl".to_string(),
                    ..policy(incoming_bam)
                },
                CacheReferencePolicy {
                    cache_version: "115".to_string(),
                    ..policy(incoming_bam)
                },
            ] {
                let saved = serde_json::to_vec_pretty(&existing).unwrap();
                std::fs::write(output.join(REFERENCE_POLICY_FILE), &saved).unwrap();
                let error = CacheReferencePolicy::preserve(&native_schema, &native, &output)
                    .unwrap_err()
                    .to_string();
                assert!(error.contains("existing reference policy"), "{error}");
                assert!(error.contains("conflicts with native cache"), "{error}");
                assert!(error.contains("rebuild into a new destination"), "{error}");
                assert_eq!(
                    std::fs::read(output.join(REFERENCE_POLICY_FILE)).unwrap(),
                    saved
                );
                assert_eq!(std::fs::read(&shard).unwrap(), original_shard);
            }
        }
    }

    #[cfg(feature = "cache-builder")]
    #[test]
    fn reference_policy_refresh_validates_existing_shards_before_atomic_write() {
        let root = tempfile::tempdir().unwrap();
        let native = root.path().join("native");
        let output = root.path().join("output");
        std::fs::create_dir(&native).unwrap();
        std::fs::write(native.join("info.txt"), b"bam\t/path.bam\n").unwrap();
        write_shard(
            &output,
            "transcript",
            "chr21",
            schema("merged", "116", None),
        );
        let native_schema = schema("merged", "116", Some(true));
        CacheReferencePolicy::preserve(&native_schema, &native, &output).unwrap();
        assert_eq!(
            CacheReferencePolicy::read(&output).unwrap(),
            Some(policy(true))
        );
        let saved = std::fs::read(output.join(REFERENCE_POLICY_FILE)).unwrap();
        CacheReferencePolicy::preserve(&native_schema, &native, &output).unwrap();
        assert_eq!(
            std::fs::read(output.join(REFERENCE_POLICY_FILE)).unwrap(),
            saved
        );

        for conflicting in [
            schema("ensembl", "116", None),
            schema("merged", "115", None),
            schema("merged", "116", Some(false)),
        ] {
            std::fs::remove_file(output.join(REFERENCE_POLICY_FILE)).unwrap();
            write_shard(&output, "transcript", "chr21", conflicting);
            let error = CacheReferencePolicy::preserve(&native_schema, &native, &output)
                .unwrap_err()
                .to_string();
            assert!(error.contains("chr21.parquet"), "{error}");
            assert!(!output.join(REFERENCE_POLICY_FILE).exists());
            std::fs::write(output.join(REFERENCE_POLICY_FILE), &saved).unwrap();
        }

        // A matching root must not hide a conflicting shard on resume either.
        for conflicting in [
            schema("ensembl", "116", None),
            schema("merged", "115", None),
            schema("merged", "116", Some(false)),
        ] {
            write_shard(&output, "transcript", "chr21", conflicting);
            let shard = output.join("transcript/chr21.parquet");
            let shard_bytes = std::fs::read(&shard).unwrap();
            let error = CacheReferencePolicy::preserve(&native_schema, &native, &output)
                .unwrap_err()
                .to_string();
            assert!(error.contains("chr21.parquet"), "{error}");
            assert_eq!(std::fs::read(&shard).unwrap(), shard_bytes);
            assert_eq!(
                std::fs::read(output.join(REFERENCE_POLICY_FILE)).unwrap(),
                saved
            );
        }
    }
}
