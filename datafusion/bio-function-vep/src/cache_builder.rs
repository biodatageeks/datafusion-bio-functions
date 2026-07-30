//! Cache builder: converts a raw Ensembl VEP cache to the partitioned Parquet
//! cache. The actual build work is delegated to [`crate::cache::build`];
//! this module provides the public [`CacheBuilder`] entry point and entity
//! dispatch. Parquet is the only supported output format.

use datafusion::catalog::TableProvider;
use datafusion::common::{DataFusionError, Result};
use log::info;

use datafusion_bio_format_ensembl_cache::{
    CacheSourceType as BioFormatsCacheSourceType, EnsemblCacheOptions, EnsemblCacheTableProvider,
    EnsemblEntityKind,
};

use crate::cache_identity::CACHE_VERSION_METADATA_KEY as VEP_CACHE_VERSION_METADATA_KEY;
use crate::vep_semantics::target_for_cache_version;

#[derive(Debug, Clone, Copy, Default, PartialEq, Eq)]
pub enum CacheFormat {
    #[default]
    Parquet,
}

impl CacheFormat {
    pub fn parse(value: &str) -> Result<Self> {
        match value.to_ascii_lowercase().as_str() {
            // "cache" is accepted as a historical alias but always resolves to
            // the Parquet cache — the only supported output format.
            "lance" | "parquet" => Ok(Self::Parquet),
            other => Err(DataFusionError::Execution(format!(
                "cache_format must be 'parquet', got '{other}'"
            ))),
        }
    }

    pub fn as_str(self) -> &'static str {
        match self {
            Self::Parquet => "parquet",
        }
    }
}

/// Statistics returned after building one entity.
#[derive(Debug, Clone)]
pub struct EntityStats {
    pub entity: String,
    /// Files written for this entity, as `(path, row_count)`. Named
    /// `parquet_files` for backward compatibility with the Python wrapper;
    /// holds the entity's Parquet shard files.
    pub parquet_files: Vec<(String, usize)>,
}

/// Builder for converting a raw Ensembl VEP cache to the partitioned Parquet cache.
pub struct CacheBuilder {
    cache_root: String,
    output_dir: String,
    partitions: usize,
    cache_format: CacheFormat,
    overwrite: bool,
    cache_source_type: BioFormatsCacheSourceType,
    expected_cache_version: Option<String>,
    chrom_filter: Option<Vec<String>>,
}

impl CacheBuilder {
    pub fn new(cache_root: impl Into<String>, output_dir: impl Into<String>) -> Self {
        Self {
            cache_root: cache_root.into(),
            output_dir: output_dir.into(),
            partitions: 8,
            cache_format: CacheFormat::default(),
            overwrite: false,
            cache_source_type: BioFormatsCacheSourceType::Ensembl,
            expected_cache_version: None,
            chrom_filter: None,
        }
    }

    pub fn with_partitions(mut self, n: usize) -> Self {
        self.partitions = n;
        self
    }

    /// Restrict the build to the given chromosomes (matched after stripping any
    /// `chr` prefix). Useful for scoped rebuilds such as a single-chromosome
    /// profiling or parity build. An empty list is treated as "no filter".
    pub fn with_chrom_filter<I, S>(mut self, chroms: I) -> Self
    where
        I: IntoIterator<Item = S>,
        S: Into<String>,
    {
        let chroms: Vec<String> = chroms.into_iter().map(Into::into).collect();
        self.chrom_filter = (!chroms.is_empty()).then_some(chroms);
        self
    }

    pub fn with_cache_format(mut self, cache_format: CacheFormat) -> Self {
        self.cache_format = cache_format;
        self
    }

    pub fn with_overwrite(mut self, overwrite: bool) -> Self {
        self.overwrite = overwrite;
        self
    }

    pub fn with_cache_source_type(mut self, cache_source_type: BioFormatsCacheSourceType) -> Self {
        self.cache_source_type = cache_source_type;
        self
    }

    /// Assert the release detected independently from the raw cache metadata.
    ///
    /// The assertion cannot label an unknown raw cache. When omitted, the raw
    /// cache must still provide a supported release and that detected value is
    /// embedded in every generated Parquet shard.
    pub fn with_expected_cache_version(mut self, cache_version: impl Into<String>) -> Self {
        self.expected_cache_version = Some(cache_version.into());
        self
    }

    /// Build all entities into the Parquet cache.
    pub async fn build_all(&self) -> Result<Vec<EntityStats>> {
        let entities = [
            "variation",
            "transcript",
            "exon",
            "translation",
            "regulatory",
            "motif",
        ];
        let mut results = Vec::new();
        for entity in entities {
            match self.build_entity(entity).await {
                Ok(stats) => results.extend(stats),
                Err(e) => {
                    let msg = e.to_string();
                    if msg.contains("No source files discovered") || msg.contains("skipped") {
                        info!("{entity}: skipped (no source files)");
                    } else {
                        return Err(e);
                    }
                }
            }
        }
        Ok(results)
    }

    /// Build a single entity. Returns one or more EntityStats (translation splits into two).
    ///
    /// When `overwrite` is false (default), skips entities whose Parquet output
    /// already exists.
    pub async fn build_entity(&self, entity: &str) -> Result<Vec<EntityStats>> {
        let kind = parse_entity(entity)
            .ok_or_else(|| DataFusionError::Execution(format!("Unknown entity: {entity}")))?;

        // Parquet is the only supported cache format.
        let CacheFormat::Parquet = self.cache_format;

        #[cfg(feature = "parquet-cache")]
        {
            let cache_version = self.validate_raw_cache_identity(kind)?;
            let options = crate::cache::build::CacheBuildOptions {
                cache_root: self.cache_root.clone(),
                output_dir: self.output_dir.clone(),
                partitions: self.partitions,
                cache_source_type: self.cache_source_type,
                cache_version,
                overwrite: self.overwrite,
                chrom_filter: self.chrom_filter.clone(),
            };
            crate::cache::build::build_parquet_entity(&options, kind).await
        }

        #[cfg(not(feature = "parquet-cache"))]
        {
            let _ = kind;
            Err(DataFusionError::Execution(
                "cache builder requires the parquet-cache feature".to_string(),
            ))
        }
    }

    #[cfg(feature = "parquet-cache")]
    fn validate_raw_cache_identity(&self, kind: EnsemblEntityKind) -> Result<String> {
        let mut options = EnsemblCacheOptions::new(&self.cache_root)
            .with_cache_source_type(self.cache_source_type);
        if let Some(expected) = self.expected_cache_version.as_deref() {
            options = options.with_expected_cache_version(expected);
        }
        let provider = EnsemblCacheTableProvider::for_entity(kind, options)?;
        let schema = provider.schema();
        let detected = schema
            .metadata()
            .get(VEP_CACHE_VERSION_METADATA_KEY)
            .filter(|value| !value.is_empty() && value.as_str() != "unknown")
            .ok_or_else(|| {
                DataFusionError::Execution(format!(
                    "cache build refused before writing output: raw cache '{}' has no usable \
                     {VEP_CACHE_VERSION_METADATA_KEY}; add cache_version/version to info.txt or \
                     use a canonical raw VEP cache-root basename",
                    self.cache_root
                ))
            })?
            .to_string();
        target_for_cache_version(&detected).map_err(|error| {
            DataFusionError::Execution(format!(
                "cache build refused before writing output: raw cache '{}' declares {error}",
                self.cache_root
            ))
        })?;

        if let Some(expected) = self.expected_cache_version.as_deref() {
            target_for_cache_version(expected).map_err(|error| {
                DataFusionError::Execution(format!(
                    "cache build refused before writing output: invalid expected cache release: \
                     {error}"
                ))
            })?;
            if detected != expected {
                return Err(DataFusionError::Execution(format!(
                    "cache build refused before writing output: raw cache '{}' declares \
                     {VEP_CACHE_VERSION_METADATA_KEY}={detected}, requested release is {expected}",
                    self.cache_root
                )));
            }
        }
        Ok(detected)
    }
}

// ---------------------------------------------------------------------------
// Helper functions (ported from vepyr convert.rs)
// ---------------------------------------------------------------------------

fn parse_entity(name: &str) -> Option<EnsemblEntityKind> {
    match name {
        "variation" => Some(EnsemblEntityKind::Variation),
        "transcript" => Some(EnsemblEntityKind::Transcript),
        "exon" => Some(EnsemblEntityKind::Exon),
        "translation" => Some(EnsemblEntityKind::Translation),
        "regulatory" => Some(EnsemblEntityKind::RegulatoryFeature),
        "motif" => Some(EnsemblEntityKind::MotifFeature),
        _ => None,
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn parse_entity_known() {
        assert_eq!(
            parse_entity("variation"),
            Some(EnsemblEntityKind::Variation)
        );
        assert_eq!(
            parse_entity("transcript"),
            Some(EnsemblEntityKind::Transcript)
        );
        assert_eq!(parse_entity("exon"), Some(EnsemblEntityKind::Exon));
        assert_eq!(
            parse_entity("translation"),
            Some(EnsemblEntityKind::Translation)
        );
        assert_eq!(
            parse_entity("regulatory"),
            Some(EnsemblEntityKind::RegulatoryFeature)
        );
        assert_eq!(parse_entity("motif"), Some(EnsemblEntityKind::MotifFeature));
    }

    #[test]
    fn parse_entity_unknown() {
        assert_eq!(parse_entity("nonsense"), None);
    }

    #[test]
    fn cache_format_parser_accepts_parquet_and_lance_alias() {
        assert_eq!(CacheFormat::parse("parquet").unwrap(), CacheFormat::Parquet);
        assert_eq!(CacheFormat::parse("lance").unwrap(), CacheFormat::Parquet);
        assert_eq!(CacheFormat::parse("LANCE").unwrap(), CacheFormat::Parquet);
        assert!(CacheFormat::parse("indexed_parquet").is_err());
        assert!(CacheFormat::parse("legacy_fjall").is_err());
    }

    #[test]
    fn cache_format_default_is_parquet() {
        assert_eq!(CacheFormat::default(), CacheFormat::Parquet);
        assert_eq!(CacheFormat::Parquet.as_str(), "parquet");
    }

    #[test]
    fn cache_builder_defaults_to_parquet() {
        let builder = CacheBuilder::new("/cache", "/output");
        assert_eq!(builder.cache_format, CacheFormat::Parquet);
        assert_eq!(builder.partitions, 8);
        assert!(!builder.overwrite);
        assert!(builder.expected_cache_version.is_none());
        assert!(builder.chrom_filter.is_none());
    }

    #[test]
    fn cache_builder_with_overrides() {
        let builder = CacheBuilder::new("/cache", "/output")
            .with_partitions(4)
            .with_cache_format(CacheFormat::Parquet)
            .with_overwrite(true)
            .with_expected_cache_version("116")
            .with_chrom_filter(["chr1", "chr2"]);
        assert_eq!(builder.partitions, 4);
        assert_eq!(builder.cache_format, CacheFormat::Parquet);
        assert!(builder.overwrite);
        assert_eq!(builder.expected_cache_version.as_deref(), Some("116"));
        assert_eq!(
            builder.chrom_filter.as_deref(),
            Some(["chr1".to_string(), "chr2".to_string()].as_slice())
        );
    }

    #[tokio::test]
    async fn build_entity_rejects_unknown_entity() {
        let builder = CacheBuilder::new("/cache", "/output");
        let err = builder.build_entity("nonsense").await.unwrap_err();
        assert!(err.to_string().contains("Unknown entity"));
    }
}
