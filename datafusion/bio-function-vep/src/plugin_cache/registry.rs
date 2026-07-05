//! Runtime plugin registry: discover built plugins, open a per-chrom
//! [`PluginLookup`] for each, and expose the concatenated CSQ field list plus a
//! single `probe_all` aligned to it.

use std::path::Path;

use datafusion::common::Result;

use crate::cache::manifest::canonical_chrom_label;
use crate::plugin_cache::cache_manifest::discover_plugins;
use crate::plugin_cache::lookup::{PluginLookup, PluginScalar};

/// Per-transcript attributes the engine supplies at probe time, keyed by the
/// `engine_attr` names a plugin's match columns (§3.4) bind to. Extend as new
/// discriminators are added.
#[derive(Debug, Default, Clone)]
pub struct EngineAttrs {
    /// `{refAA}{protein_position}{altAA}` (e.g. `"V550A"`); `None` for a
    /// non-missense consequence (no amino-acid change) — the gate.
    pub amino_acid_change: Option<String>,
}

impl EngineAttrs {
    /// Resolve a named engine attribute; unknown names resolve to `None`.
    fn get(&self, engine_attr: &str) -> Option<String> {
        match engine_attr {
            "amino_acid_change" => self.amino_acid_change.clone(),
            _ => None,
        }
    }
}

/// One enabled plugin, with its per-chrom lookup (absent if this plugin has no
/// shard for the current chrom).
pub struct PluginEntry {
    pub name: String,
    pub csq_fields: Vec<String>,
    /// The `engine_attr` name for each match column, in order.
    match_engine_attrs: Vec<String>,
    lookup: Option<PluginLookup>,
}

/// All enabled plugins for one chromosome, in plugin-name order.
pub struct PluginRegistry {
    plugins: Vec<PluginEntry>,
}

impl PluginRegistry {
    /// Discover plugins under `cache_root` and open each one's shard for `chrom`.
    pub async fn open(cache_root: &Path, chrom: &str) -> Result<Self> {
        let want = canonical_chrom_label(chrom);
        let manifests = discover_plugins(cache_root)?;
        let mut plugins = Vec::with_capacity(manifests.len());
        for m in manifests {
            let csq_fields: Vec<String> = m
                .value_columns
                .iter()
                .map(|v| v.csq_field.clone())
                .collect();
            let value_columns: Vec<String> =
                m.value_columns.iter().map(|v| v.column.clone()).collect();
            let match_columns: Vec<String> =
                m.match_columns.iter().map(|mc| mc.column.clone()).collect();
            let match_engine_attrs: Vec<String> = m
                .match_columns
                .iter()
                .map(|mc| mc.engine_attr.clone())
                .collect();
            let chrom_entry = m.chroms.iter().find(|c| c.chrom == want);
            let lookup = match chrom_entry {
                Some(entry) => {
                    let shard = cache_root
                        .join("plugin")
                        .join(&m.plugin_name)
                        .join(&entry.file);
                    Some(PluginLookup::open(&shard, match_columns, value_columns).await?)
                }
                None => None,
            };
            plugins.push(PluginEntry {
                name: m.plugin_name,
                csq_fields,
                match_engine_attrs,
                lookup,
            });
        }
        Ok(Self { plugins })
    }

    /// True when no plugins are enabled (fast-path skip at the call site).
    pub fn is_empty(&self) -> bool {
        self.plugins.is_empty()
    }

    /// Concatenated CSQ field names across all plugins, in plugin-name order.
    pub fn csq_fields(&self) -> Vec<String> {
        self.plugins
            .iter()
            .flat_map(|p| p.csq_fields.iter().cloned())
            .collect()
    }

    /// Probe every plugin for `(start, allele_string, attrs)`, returning one
    /// [`PluginScalar`] per entry of [`Self::csq_fields`] (in the same order).
    /// Each plugin's match discriminator is built by resolving its match columns'
    /// `engine_attr`s against `attrs`. A plugin with no shard for this chrom, or a
    /// position/allele/discriminator miss, yields `PluginScalar::Null` per field
    /// (the per-transcript gate for match-column plugins).
    pub fn probe_all(
        &self,
        start: u32,
        allele_string: &str,
        attrs: &EngineAttrs,
    ) -> Result<Vec<PluginScalar>> {
        let mut out = Vec::new();
        for p in &self.plugins {
            let match_values: Vec<Option<String>> =
                p.match_engine_attrs.iter().map(|a| attrs.get(a)).collect();
            let hit = match &p.lookup {
                Some(lk) => lk.probe(start, allele_string, &match_values)?,
                None => None,
            };
            match hit {
                Some(row) => out.extend(row.values),
                None => out.extend(std::iter::repeat_n(PluginScalar::Null, p.csq_fields.len())),
            }
        }
        Ok(out)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::plugin_cache::cache_manifest::{
        CacheManifest, ChromEntry, TierRecord, ValueColumnRecord,
    };
    use crate::plugin_cache::source_manifest::{ValueColumn, ValueType};
    use crate::plugin_cache::write::{PluginShardWriter, plugin_output_schema};
    use datafusion::arrow::array::{Float32Array, Int8Array, StringArray, UInt32Array};
    use datafusion::arrow::record_batch::RecordBatch;
    use std::sync::Arc;

    #[tokio::test(flavor = "multi_thread")]
    async fn discovers_and_probes() {
        let dir = tempfile::tempdir().unwrap();
        let cache_root = dir.path();
        let plugin_dir = cache_root.join("plugin").join("alphamissense");
        std::fs::create_dir_all(&plugin_dir).unwrap();

        // Build the chr1 shard.
        let vals = vec![ValueColumn {
            column: "am_pathogenicity".into(),
            csq_field: "am_pathogenicity".into(),
            ty: ValueType::Float32,
        }];
        let schema = plugin_output_schema(&[], &vals);
        let batch = RecordBatch::try_new(
            schema.clone(),
            vec![
                Arc::new(StringArray::from(vec!["1"])),
                Arc::new(UInt32Array::from(vec![100u32])),
                Arc::new(UInt32Array::from(vec![100u32])),
                Arc::new(StringArray::from(vec!["A/G"])),
                Arc::new(Float32Array::from(vec![0.9f32])),
                Arc::new(Int8Array::from(vec![0i8])),
            ],
        )
        .unwrap();
        let mut w = PluginShardWriter::create(&plugin_dir.join("chr1.parquet"), schema).unwrap();
        w.write(&batch).unwrap();
        w.finish().unwrap();

        // Write its cache manifest.
        let manifest = CacheManifest {
            plugin_name: "alphamissense".into(),
            source_manifest: "alphamissense.source.toml".into(),
            key_columns: vec![
                "chrom".into(),
                "start".into(),
                "end".into(),
                "allele_string".into(),
            ],
            match_columns: vec![],
            value_columns: vec![ValueColumnRecord {
                column: "am_pathogenicity".into(),
                csq_field: "am_pathogenicity".into(),
                ty: "Float32".into(),
            }],
            tier: TierRecord {
                threshold: 0.01,
                unmatched: "cold".into(),
            },
            chroms: vec![ChromEntry {
                chrom: "chr1".into(),
                file: "chr1.parquet".into(),
                rows: 1,
                warm: 1,
                cold: 0,
            }],
            cache_source_version: None,
        };
        manifest.write(&plugin_dir).unwrap();

        let reg = PluginRegistry::open(cache_root, "1").await.unwrap();
        assert_eq!(reg.csq_fields(), vec!["am_pathogenicity".to_string()]);
        let attrs = EngineAttrs::default();
        let hit = reg.probe_all(100, "A/G", &attrs).unwrap();
        assert_eq!(hit.len(), 1);
        match &hit[0] {
            PluginScalar::F32(v) => assert!((*v - 0.9).abs() < 1e-6),
            other => panic!("{other:?}"),
        }
        // miss → Null aligned to the single field
        let miss = reg.probe_all(999, "A/G", &attrs).unwrap();
        assert_eq!(miss, vec![PluginScalar::Null]);
    }
}
