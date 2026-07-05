//! Cache manifest (build output) — emitted per plugin under `plugin/<name>/`,
//! discovered at runtime by scanning `plugin/*/manifest.json`. Carries the
//! value→CSQ mapping (single source of truth for the shard schema and emitted
//! fields) and per-chrom row/tier counts.

use std::path::Path;

use datafusion::common::{DataFusionError, Result};
use serde::{Deserialize, Serialize};

use crate::plugin_cache::source_manifest::{SourceManifest, ValueType};

/// Tier policy actually used at build time.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct TierRecord {
    pub threshold: f64,
    pub unmatched: String,
}

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
}

/// The cache manifest written by the build and read at runtime.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct CacheManifest {
    pub plugin_name: String,
    pub source_manifest: String,
    pub key_columns: Vec<String>,
    pub value_columns: Vec<ValueColumnRecord>,
    pub tier: TierRecord,
    pub chroms: Vec<ChromEntry>,
    pub cache_source_version: Option<String>,
}

fn ty_str(t: ValueType) -> &'static str {
    match t {
        ValueType::Utf8 => "Utf8",
        ValueType::Float32 => "Float32",
        ValueType::Int32 => "Int32",
    }
}

impl CacheManifest {
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
            value_columns: src
                .value_columns
                .iter()
                .map(|v| ValueColumnRecord {
                    column: v.column.clone(),
                    csq_field: v.csq_field.clone(),
                    ty: ty_str(v.ty).into(),
                })
                .collect(),
            tier: TierRecord {
                threshold: src.tier.threshold,
                unmatched: src.tier.unmatched.clone().unwrap_or_else(|| "cold".into()),
            },
            chroms: vec![],
            cache_source_version: None,
        }
    }

    /// Write `manifest.json` into `plugin_dir`.
    pub fn write(&self, plugin_dir: &Path) -> Result<()> {
        let json = serde_json::to_string_pretty(self)
            .map_err(|e| DataFusionError::Execution(format!("serialize manifest: {e}")))?;
        std::fs::write(plugin_dir.join("manifest.json"), json)
            .map_err(|e| DataFusionError::Execution(format!("write manifest: {e}")))
    }
}

/// Discover built plugins under `<cache_root>/plugin/*/manifest.json`, sorted by
/// plugin name for deterministic CSQ field ordering.
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
            out.push(
                serde_json::from_str(&text).map_err(|e| {
                    DataFusionError::Execution(format!("parse {}: {e}", mf.display()))
                })?,
            );
        }
    }
    out.sort_by(|a: &CacheManifest, b: &CacheManifest| a.plugin_name.cmp(&b.plugin_name));
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
        let m = CacheManifest {
            plugin_name: "demo".into(),
            source_manifest: "demo.source.toml".into(),
            key_columns: vec![
                "chrom".into(),
                "start".into(),
                "end".into(),
                "allele_string".into(),
            ],
            value_columns: vec![],
            tier: TierRecord {
                threshold: 0.01,
                unmatched: "cold".into(),
            },
            chroms: vec![ChromEntry {
                chrom: "22".into(),
                file: "chr22.parquet".into(),
                rows: 3,
                warm: 1,
                cold: 2,
            }],
            cache_source_version: None,
        };
        m.write(&plugin_dir).unwrap();
        let found = discover_plugins(dir.path()).unwrap();
        assert_eq!(found.len(), 1);
        assert_eq!(found[0].plugin_name, "demo");
        assert_eq!(found[0].chroms[0].rows, 3);
    }
}
