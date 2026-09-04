//! Cache manifest (build output) — emitted per plugin under `plugin/<name>/`,
//! discovered at runtime by scanning `plugin/*/manifest.json`. Carries the
//! value→CSQ mapping (single source of truth for the shard schema and emitted
//! fields) and per-chrom row/tier counts.

use std::path::Path;

use datafusion::common::{DataFusionError, Result};
use serde::{Deserialize, Serialize};

use crate::plugin_cache::source_manifest::{
    SourceManifest, SourceSpec, ValueType, validate_csq_field_name,
};

/// Per-chromosome build result.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct ChromEntry {
    pub chrom: String,
    pub file: String,
    pub rows: usize,
    pub warm: usize,
    pub cold: usize,
}

/// A value column and its CSQ field, as recorded in the built cache.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct ValueColumnRecord {
    pub column: String,
    pub csq_field: String,
    #[serde(rename = "type")]
    pub ty: String,
    /// Free text for this field's `##<CSQ_FIELD>=<description>` VCF header
    /// line. Ensembl plugins supply one via `get_header_info()`; absent here
    /// means no header line is written for the field.
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub description: Option<String>,
}

/// A per-transcript match discriminator binding, recorded so the runtime knows
/// which template builds each match column's discriminator value.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct MatchColumnRecord {
    pub column: String,
    pub template: String,
}

/// Provenance of one build input, as recorded in the built cache so a shard
/// can be traced to exact bytes without the source file.
///
/// `md5`/`url`/`path_md5` are copied from the source manifest; `verified_md5`
/// is the digest the build actually computed over the resolved file (absent
/// when verification was skipped or the manifest declared no digest), and
/// `file`/`size`/`mtime_ns`/`ino`/`ctime_ns` fingerprint that file so an
/// incremental per-chromosome build can trust an earlier verification instead
/// of re-hashing. `ino` and `ctime_ns` are the replacement-sensitive part: a
/// copy can preserve size and mtime, but a new inode and a fresh change time
/// come with any replacement or rewrite.
#[derive(Debug, Clone, PartialEq, Eq, Serialize, Deserialize)]
pub struct SourceRecord {
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub part: Option<String>,
    /// File name (not the full path) of the build input.
    pub file: String,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub url: Option<String>,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub md5: Option<String>,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub path_md5: Option<String>,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub verified_md5: Option<String>,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub size: Option<u64>,
    /// Modification time of the hashed file, nanoseconds since the Unix epoch.
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub mtime_ns: Option<i64>,
    /// Inode of the hashed file (Unix only).
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub ino: Option<u64>,
    /// Inode change time of the hashed file, nanoseconds since the Unix epoch
    /// (Unix); creation time on Windows. Absent where the platform offers
    /// neither, in which case the file is always re-hashed.
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub ctime_ns: Option<i64>,
    /// The source's index file, for tabix-indexed sources.
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub index: Option<IndexRecord>,
}

impl SourceRecord {
    /// Provenance as declared, before any verification.
    pub fn from_spec(spec: &SourceSpec) -> Self {
        let file = Path::new(&spec.path)
            .file_name()
            .map(|f| f.to_string_lossy().into_owned())
            .unwrap_or_else(|| spec.path.clone());
        let index = spec.index.map(|_| IndexRecord {
            file: format!("{file}.tbi"),
            md5: spec.index_md5.clone(),
            verified_md5: None,
            size: None,
            mtime_ns: None,
            ino: None,
            ctime_ns: None,
        });
        SourceRecord {
            part: spec.part.clone(),
            file,
            url: spec.url.clone(),
            md5: spec.md5.clone(),
            path_md5: spec.path_md5.clone(),
            verified_md5: None,
            size: None,
            mtime_ns: None,
            ino: None,
            ctime_ns: None,
            index,
        }
    }

    /// The digest this record declares for its build input: `path_md5`, else `md5`.
    pub fn expected_md5(&self) -> Option<&str> {
        self.path_md5.as_deref().or(self.md5.as_deref())
    }
}

/// Provenance of a source's index file (the tabix `.tbi`), which decides
/// which records each chromosome build reads and so is verified like the data.
#[derive(Debug, Clone, PartialEq, Eq, Serialize, Deserialize)]
pub struct IndexRecord {
    pub file: String,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub md5: Option<String>,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub verified_md5: Option<String>,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub size: Option<u64>,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub mtime_ns: Option<i64>,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub ino: Option<u64>,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub ctime_ns: Option<i64>,
}

/// The cache manifest written by the build and read at runtime.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct CacheManifest {
    pub plugin_name: String,
    pub source_manifest: String,
    pub key_columns: Vec<String>,
    #[serde(default)]
    pub match_columns: Vec<MatchColumnRecord>,
    pub value_columns: Vec<ValueColumnRecord>,
    pub chroms: Vec<ChromEntry>,
    /// Build-input provenance, one entry per `[[source]]`. Empty for caches
    /// built before provenance was recorded.
    #[serde(default, skip_serializing_if = "Vec::is_empty")]
    pub sources: Vec<SourceRecord>,
    pub cache_source_version: Option<String>,
    /// How Ensembl's own plugin matches a variant to a data row. Defaults to
    /// [`AlleleMatch::Exact`], which is what most plugins do.
    #[serde(default)]
    pub allele_match: AlleleMatch,
    /// Position of this plugin's field block in the CSQ string. Ensembl emits
    /// `--plugin` blocks in the order the flags were given and appends
    /// `--custom` blocks last, so the order is a convention rather than
    /// something derivable from the cache; it is pinned here. Plugins sharing a
    /// rank fall back to plugin-name order.
    #[serde(default = "default_csq_rank")]
    pub csq_rank: u32,
    /// Field order *within* this plugin's block. Ensembl's plugin framework
    /// emits a plugin's own fields sorted by name, while `--custom` fields keep
    /// the order the user listed them in.
    #[serde(default)]
    pub field_order: FieldOrder,
}

fn default_csq_rank() -> u32 {
    u32::MAX
}

/// Same default, reachable from the source manifest's serde attribute.
pub fn default_csq_rank_pub() -> u32 {
    default_csq_rank()
}

/// Field ordering within one plugin's CSQ block.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Default, Serialize, Deserialize)]
#[serde(rename_all = "lowercase")]
pub enum FieldOrder {
    /// Manifest declaration order — matches Ensembl `--custom` behaviour.
    #[default]
    Declared,
    /// Sorted by CSQ field name — matches Ensembl's `--plugin` behaviour.
    Alphabetical,
}

/// The allele-comparison rule a plugin's Ensembl implementation uses, which
/// decides whether an equal-length substitution may be reduced before probing.
///
/// `Minimised` plugins call `get_matched_variant_alleles()`, which runs
/// `trim_sequences()` over both the variant and the data row — so an untrimmed
/// MNV such as `TTGTGTGTGTGTG/GTGTGTGTGTGTG` matches CADD's `T/G` row.
/// `Exact` plugins compare `(start, ref, alt)` verbatim (SpliceAI, dbNSFP) and
/// must never be probed with a reduced allele, or they gain hits Ensembl does
/// not report.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Default, serde::Serialize, serde::Deserialize)]
#[serde(rename_all = "lowercase")]
pub enum AlleleMatch {
    #[default]
    Exact,
    Minimised,
}

fn ty_str(t: ValueType) -> &'static str {
    match t {
        ValueType::Utf8 => "Utf8",
        ValueType::Float32 => "Float32",
        ValueType::Int32 => "Int32",
    }
}

impl CacheManifest {
    fn validate(&self) -> Result<()> {
        for value in &self.value_columns {
            validate_csq_field_name(&self.plugin_name, &value.csq_field)?;
        }
        Ok(())
    }

    /// Seed a cache manifest from a source manifest (chroms filled in by the build).
    pub fn from_source(src: &SourceManifest, source_manifest_file: &str) -> Self {
        CacheManifest {
            plugin_name: src.plugin_name.clone(),
            source_manifest: source_manifest_file.to_string(),
            key_columns: vec![
                "chrom".into(),
                "start".into(),
                "end".into(),
                "allele_string".into(),
            ],
            match_columns: src
                .match_columns
                .iter()
                .map(|m| MatchColumnRecord {
                    column: m.column.clone(),
                    template: m.template.clone(),
                })
                .collect(),
            value_columns: src
                .value_columns
                .iter()
                .map(|v| ValueColumnRecord {
                    column: v.column.clone(),
                    csq_field: v.csq_field.clone(),
                    ty: ty_str(v.ty).into(),
                    description: v.description.clone(),
                })
                .collect(),
            chroms: vec![],
            // Filled in by the build, like `chroms`: `build_all` verifies the
            // sources and records what it found.
            sources: vec![],
            cache_source_version: None,
            allele_match: src.allele_match,
            csq_rank: src.csq_rank,
            field_order: src.field_order,
        }
    }

    /// Write `manifest.json` into `plugin_dir`.
    pub fn write(&self, plugin_dir: &Path) -> Result<()> {
        self.validate()?;
        let json = serde_json::to_string_pretty(self)
            .map_err(|e| DataFusionError::Execution(format!("serialize manifest: {e}")))?;
        std::fs::write(plugin_dir.join("manifest.json"), json)
            .map_err(|e| DataFusionError::Execution(format!("write manifest: {e}")))
    }
}

/// Discover built plugins under `<cache_root>/plugin/*/manifest.json`, ordered by
/// `csq_rank` then plugin name so the emitted CSQ field blocks are deterministic
/// and can be pinned to the order Ensembl VEP writes them in.
pub fn discover_plugins(cache_root: &Path) -> Result<Vec<CacheManifest>> {
    let plugin_root = cache_root.join("plugin");
    let mut out = Vec::new();
    if !plugin_root.exists() {
        return Ok(out);
    }
    for entry in std::fs::read_dir(&plugin_root)
        .map_err(|e| DataFusionError::Execution(format!("read {}: {e}", plugin_root.display())))?
    {
        let dir = entry
            .map_err(|e| DataFusionError::Execution(format!("dir entry: {e}")))?
            .path();
        let mf = dir.join("manifest.json");
        if mf.exists() {
            let text = std::fs::read_to_string(&mf)
                .map_err(|e| DataFusionError::Execution(format!("read {}: {e}", mf.display())))?;
            let manifest: CacheManifest = serde_json::from_str(&text)
                .map_err(|e| DataFusionError::Execution(format!("parse {}: {e}", mf.display())))?;
            manifest.validate()?;
            out.push(manifest);
        }
    }
    out.sort_by(|a: &CacheManifest, b: &CacheManifest| {
        (a.csq_rank, &a.plugin_name).cmp(&(b.csq_rank, &b.plugin_name))
    });
    Ok(out)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn writes_and_discovers_manifest() {
        let dir = tempfile::tempdir().unwrap();
        let plugin_dir = dir.path().join("plugin").join("demo");
        std::fs::create_dir_all(&plugin_dir).unwrap();
        let mut m = CacheManifest {
            plugin_name: "demo".into(),
            source_manifest: "demo.source.toml".into(),
            key_columns: vec![
                "chrom".into(),
                "start".into(),
                "end".into(),
                "allele_string".into(),
            ],
            match_columns: vec![],
            value_columns: vec![],
            chroms: vec![ChromEntry {
                chrom: "22".into(),
                file: "chr22.parquet".into(),
                rows: 3,
                warm: 1,
                cold: 2,
            }],
            sources: vec![],
            cache_source_version: None,
            allele_match: Default::default(),
            csq_rank: 0,
            field_order: Default::default(),
        };
        m.write(&plugin_dir).unwrap();
        let found = discover_plugins(dir.path()).unwrap();
        assert_eq!(found.len(), 1);
        assert_eq!(found[0].plugin_name, "demo");
        assert_eq!(found[0].chroms[0].rows, 3);

        m.value_columns = vec![ValueColumnRecord {
            column: "score".into(),
            csq_field: "fileformat".into(),
            ty: "Utf8".into(),
            description: Some("must not replace VCFv4.2".into()),
        }];
        let error = m.write(&plugin_dir).unwrap_err().to_string();
        assert!(
            error.contains("reserved VCF meta-information key"),
            "{error}"
        );

        // Runtime discovery also protects old or externally-created caches
        // that did not pass through this version's `write` method.
        std::fs::write(
            plugin_dir.join("manifest.json"),
            serde_json::to_string_pretty(&m).unwrap(),
        )
        .unwrap();
        let error = discover_plugins(dir.path()).unwrap_err().to_string();
        assert!(
            error.contains("reserved VCF meta-information key"),
            "{error}"
        );

        m.value_columns[0].csq_field = "SCORE\n##injected".into();
        std::fs::write(
            plugin_dir.join("manifest.json"),
            serde_json::to_string_pretty(&m).unwrap(),
        )
        .unwrap();
        let error = discover_plugins(dir.path()).unwrap_err().to_string();
        assert!(error.contains("invalid CSQ field name"), "{error}");
    }
}
