//! Runtime plugin registry: discover built plugins, open a per-chrom
//! [`PluginLookup`] for each, and drive the buffer-batched lookup — one
//! [`PluginLookup::take_buffer`] per plugin per buffer, then synchronous
//! per-transcript probes against the resulting [`PluginBufferSlice`]s.

use std::path::Path;

use datafusion::common::Result;

use crate::cache::manifest::canonical_chrom_label;
use crate::plugin_cache::cache_manifest::discover_plugins;
use crate::plugin_cache::lookup::{PluginBufferSlice, PluginLookup, PluginScalar};

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
struct PluginEntry {
    csq_fields: Vec<String>,
    /// The `engine_attr` name for each match column, in order.
    match_engine_attrs: Vec<String>,
    n_match: usize,
    n_values: usize,
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
            let n_match = match_columns.len();
            let n_values = value_columns.len();
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
                csq_fields,
                match_engine_attrs,
                n_match,
                n_values,
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

    /// Take the candidate rows for one buffer from every plugin shard — one
    /// page-scoped [`PluginLookup::take_buffer`] per plugin — into per-plugin
    /// [`PluginBufferSlice`]s. `sorted_unique_starts` must be sorted+deduped.
    pub async fn take_buffer_all(&self, sorted_unique_starts: &[u32]) -> Result<BufferSlices> {
        let mut entries = Vec::with_capacity(self.plugins.len());
        for p in &self.plugins {
            let slice = match &p.lookup {
                Some(lk) => {
                    let batch = lk.take_buffer(sorted_unique_starts).await?;
                    Some(PluginBufferSlice::from_batch(
                        &batch, p.n_match, p.n_values,
                    )?)
                }
                None => None,
            };
            entries.push(SliceEntry {
                csq_fields_len: p.csq_fields.len(),
                match_engine_attrs: p.match_engine_attrs.clone(),
                slice,
            });
        }
        Ok(BufferSlices { entries })
    }
}

/// Per-plugin buffer slice plus the metadata `probe_all` needs.
struct SliceEntry {
    csq_fields_len: usize,
    match_engine_attrs: Vec<String>,
    slice: Option<PluginBufferSlice>,
}

/// The per-buffer working set across all plugins. Probed synchronously per
/// transcript consequence.
pub struct BufferSlices {
    entries: Vec<SliceEntry>,
}

impl BufferSlices {
    /// Probe every plugin for `(start, allele_string, attrs)`, returning one
    /// [`PluginScalar`] per entry of [`PluginRegistry::csq_fields`] (same order).
    /// Each plugin's discriminator is built by resolving its match columns'
    /// `engine_attr`s against `attrs`; a plugin with no shard, or a
    /// position/allele/discriminator miss, yields `PluginScalar::Null` per field
    /// (the per-transcript gate for match-column plugins).
    pub fn probe_all(
        &self,
        start: u32,
        allele_string: &str,
        attrs: &EngineAttrs,
    ) -> Vec<PluginScalar> {
        let mut out = Vec::new();
        for e in &self.entries {
            let match_values: Vec<Option<String>> =
                e.match_engine_attrs.iter().map(|a| attrs.get(a)).collect();
            let hit = e
                .slice
                .as_ref()
                .and_then(|s| s.probe(start, allele_string, &match_values));
            match hit {
                Some(values) => out.extend(values),
                None => out.extend(std::iter::repeat_n(PluginScalar::Null, e.csq_fields_len)),
            }
        }
        out
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::plugin_cache::cache_manifest::{
        CacheManifest, ChromEntry, MatchColumnRecord, TierRecord, ValueColumnRecord,
    };
    use crate::plugin_cache::source_manifest::{MatchColumn, ValueColumn, ValueType};
    use crate::plugin_cache::write::{PluginShardWriter, plugin_output_schema};
    use datafusion::arrow::array::{Float32Array, Int8Array, StringArray, UInt32Array};
    use datafusion::arrow::record_batch::RecordBatch;
    use std::sync::Arc;

    #[tokio::test(flavor = "multi_thread")]
    async fn discovers_takes_and_probes() {
        let dir = tempfile::tempdir().unwrap();
        let cache_root = dir.path();
        let plugin_dir = cache_root.join("plugin").join("alphamissense");
        std::fs::create_dir_all(&plugin_dir).unwrap();

        // chr22 shard with a per-transcript discriminator; 100 is multi-isoform.
        let matches = vec![MatchColumn {
            column: "protein_variant".into(),
            engine_attr: "amino_acid_change".into(),
        }];
        let vals = vec![ValueColumn {
            column: "am_pathogenicity".into(),
            csq_field: "am_pathogenicity".into(),
            ty: ValueType::Float32,
        }];
        let schema = plugin_output_schema(&matches, &vals);
        let batch = RecordBatch::try_new(
            schema.clone(),
            vec![
                Arc::new(StringArray::from(vec!["22", "22"])),
                Arc::new(UInt32Array::from(vec![100u32, 100])),
                Arc::new(UInt32Array::from(vec![100u32, 100])),
                Arc::new(StringArray::from(vec!["A/G", "A/G"])),
                Arc::new(StringArray::from(vec!["R12G", "R78G"])),
                Arc::new(Float32Array::from(vec![0.0392f32, 0.0427])),
                Arc::new(Int8Array::from(vec![1i8, 1])),
            ],
        )
        .unwrap();
        let mut w = PluginShardWriter::create(&plugin_dir.join("chr22.parquet"), schema).unwrap();
        w.write(&batch).unwrap();
        w.finish().unwrap();

        let manifest = CacheManifest {
            plugin_name: "alphamissense".into(),
            source_manifest: "alphamissense.source.toml".into(),
            key_columns: vec![
                "chrom".into(),
                "start".into(),
                "end".into(),
                "allele_string".into(),
            ],
            match_columns: vec![MatchColumnRecord {
                column: "protein_variant".into(),
                engine_attr: "amino_acid_change".into(),
            }],
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
                chrom: "chr22".into(),
                file: "chr22.parquet".into(),
                rows: 2,
                warm: 0,
                cold: 2,
            }],
            cache_source_version: None,
        };
        manifest.write(&plugin_dir).unwrap();

        let reg = PluginRegistry::open(cache_root, "22").await.unwrap();
        assert_eq!(reg.csq_fields(), vec!["am_pathogenicity".to_string()]);

        let slices = reg.take_buffer_all(&[100]).await.unwrap();
        // isoform-specific hit
        let attrs = EngineAttrs {
            amino_acid_change: Some("R78G".into()),
        };
        let hit = slices.probe_all(100, "A/G", &attrs);
        match hit[0] {
            PluginScalar::F32(v) => assert!((v - 0.0427).abs() < 1e-6),
            ref other => panic!("{other:?}"),
        }
        // non-missense (no aa-change) → Null (gate)
        let none = slices.probe_all(100, "A/G", &EngineAttrs::default());
        assert_eq!(none, vec![PluginScalar::Null]);
    }
}
