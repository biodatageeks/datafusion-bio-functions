//! Runtime plugin registry: discover built plugins, open a per-chrom
//! [`PluginLookup`] for each, and drive the buffer-batched lookup — one
//! [`PluginLookup::take_buffer`] per plugin per buffer, then synchronous
//! per-transcript probes against the resulting [`PluginBufferSlice`]s.

use std::path::Path;

use datafusion::common::{DataFusionError, Result};

use crate::cache::manifest::canonical_chrom_label;
use crate::plugin_cache::cache_manifest::discover_plugins;
use crate::plugin_cache::lookup::{PluginBufferSlice, PluginLookup, PluginScalar};
use crate::plugin_cache::template::CompiledTemplate;

/// One enabled plugin, with its per-chrom lookup (absent if this plugin has no
/// shard for the current chrom).
struct PluginEntry {
    csq_fields: Vec<String>,
    /// The compiled discriminator template for each match column, in order.
    match_templates: Vec<CompiledTemplate>,
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
            let match_templates: Vec<CompiledTemplate> = m
                .match_columns
                .iter()
                .map(|mc| CompiledTemplate::compile(&mc.template))
                .collect::<Result<Vec<_>>>()?;
            let n_match = match_columns.len();
            let n_values = value_columns.len();
            let chrom_entry = m.chroms.iter().find(|c| c.chrom == want);
            // An empty chrom (rows == 0) legitimately has no shard on disk (build
            // removes any stale file) → `None` = empty plugin fields. But a chrom
            // the manifest says has rows MUST have its shard: a missing file then
            // means a partial/corrupt cache, and silently emitting nulls would
            // corrupt annotations while the header still lists the plugin fields —
            // so fail loudly instead.
            let lookup = match chrom_entry {
                Some(entry) if entry.rows > 0 => {
                    let shard = cache_root
                        .join("plugin")
                        .join(&m.plugin_name)
                        .join(&entry.file);
                    if !shard.exists() {
                        return Err(DataFusionError::Execution(format!(
                            "plugin '{}' manifest lists {} rows for chrom '{}' but its shard is \
                             missing (partial/corrupt cache): {}",
                            m.plugin_name,
                            entry.rows,
                            want,
                            shard.display()
                        )));
                    }
                    Some(PluginLookup::open(&shard, match_columns, value_columns).await?)
                }
                _ => None,
            };
            plugins.push(PluginEntry {
                csq_fields,
                match_templates,
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

    /// The concatenated CSQ field names for a plugin cache, read from manifests
    /// **without opening any shard** — cheap, contig-independent, for the VCF
    /// header. Best-effort: unreadable/absent → empty.
    pub fn field_names(cache_root: &Path) -> Vec<String> {
        discover_plugins(cache_root)
            .map(|manifests| {
                manifests
                    .iter()
                    .flat_map(|m| m.value_columns.iter().map(|v| v.csq_field.clone()))
                    .collect()
            })
            .unwrap_or_default()
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
                match_templates: p.match_templates.clone(),
                slice,
            });
        }
        Ok(BufferSlices { entries })
    }
}

/// Per-plugin buffer slice plus the metadata `probe_all` needs.
struct SliceEntry {
    csq_fields_len: usize,
    match_templates: Vec<CompiledTemplate>,
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
    /// Each plugin's discriminator is built by evaluating its match columns'
    /// templates against the engine-attribute namespace `attrs` (same order as
    /// [`crate::plugin_cache::template::ATTR_NAMES`]); a plugin with no shard, or
    /// a position/allele/discriminator miss, yields `PluginScalar::Null` per
    /// field (the per-transcript gate for match-column plugins).
    pub fn probe_all(
        &self,
        start: u32,
        allele_string: &str,
        attrs: &[Option<&str>],
    ) -> Vec<PluginScalar> {
        let mut out = Vec::new();
        for e in &self.entries {
            let match_values: Vec<Option<String>> =
                e.match_templates.iter().map(|t| t.eval(attrs)).collect();
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
        CacheManifest, ChromEntry, MatchColumnRecord, ValueColumnRecord,
    };
    use crate::plugin_cache::source_manifest::{MatchColumn, ValueColumn, ValueType};
    use crate::plugin_cache::template::build_attr_namespace;
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
            template: "{ref_aa}{Protein_position}{alt_aa}".into(),
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
                template: "{ref_aa}{Protein_position}{alt_aa}".into(),
            }],
            value_columns: vec![ValueColumnRecord {
                column: "am_pathogenicity".into(),
                csq_field: "am_pathogenicity".into(),
                ty: "Float32".into(),
            }],
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
        // isoform-specific hit: template {ref_aa}{Protein_position}{alt_aa} -> "R78G"
        let ns_hit = build_attr_namespace(
            "missense_variant",
            "",
            "",
            "",
            "",
            "",
            "",
            "",
            "",
            "78",
            "R/G",
            "",
            "A",
            "G",
        );
        let hit = slices.probe_all(100, "A/G", &ns_hit);
        match hit[0] {
            PluginScalar::F32(v) => assert!((v - 0.0427).abs() < 1e-6),
            ref other => panic!("{other:?}"),
        }
        // non-missense (no aa-change) → Null (gate)
        let ns_miss = build_attr_namespace(
            "synonymous_variant",
            "",
            "",
            "",
            "",
            "",
            "",
            "",
            "",
            "",
            "",
            "",
            "A",
            "G",
        );
        let none = slices.probe_all(100, "A/G", &ns_miss);
        assert_eq!(none, vec![PluginScalar::Null]);
    }

    // An empty chrom (rows: 0, no shard file) must open without error and probe
    // to empty fields — not fail on a missing/stale shard (PR #190 C3).
    #[tokio::test(flavor = "multi_thread")]
    async fn empty_chrom_entry_opens_without_shard() {
        let dir = tempfile::tempdir().unwrap();
        let cache_root = dir.path();
        let plugin_dir = cache_root.join("plugin").join("demo");
        std::fs::create_dir_all(&plugin_dir).unwrap();
        let manifest = CacheManifest {
            plugin_name: "demo".into(),
            source_manifest: "demo.source.toml".into(),
            key_columns: vec![
                "chrom".into(),
                "start".into(),
                "end".into(),
                "allele_string".into(),
            ],
            match_columns: vec![],
            value_columns: vec![ValueColumnRecord {
                column: "score".into(),
                csq_field: "SCORE".into(),
                ty: "Float32".into(),
            }],
            // chr22 present in the manifest with rows: 0 and NO file on disk.
            chroms: vec![ChromEntry {
                chrom: "chr22".into(),
                file: "chr22.parquet".into(),
                rows: 0,
                warm: 0,
                cold: 0,
            }],
            cache_source_version: None,
        };
        manifest.write(&plugin_dir).unwrap();

        let reg = PluginRegistry::open(cache_root, "22").await.unwrap();
        assert_eq!(reg.csq_fields(), vec!["SCORE".to_string()]);
        let slices = reg.take_buffer_all(&[100]).await.unwrap();
        // No shard → empty (Null) field, not an error.
        assert_eq!(slices.probe_all(100, "A/G", &[]), vec![PluginScalar::Null]);
    }

    // A manifest advertising rows > 0 with no shard on disk is a partial/corrupt
    // cache → open must error, not silently null the fields (PR #190 N2).
    #[tokio::test(flavor = "multi_thread")]
    async fn missing_shard_for_nonempty_chrom_errors() {
        let dir = tempfile::tempdir().unwrap();
        let cache_root = dir.path();
        let plugin_dir = cache_root.join("plugin").join("demo");
        std::fs::create_dir_all(&plugin_dir).unwrap();
        let manifest = CacheManifest {
            plugin_name: "demo".into(),
            source_manifest: "demo.source.toml".into(),
            key_columns: vec![
                "chrom".into(),
                "start".into(),
                "end".into(),
                "allele_string".into(),
            ],
            match_columns: vec![],
            value_columns: vec![ValueColumnRecord {
                column: "score".into(),
                csq_field: "SCORE".into(),
                ty: "Float32".into(),
            }],
            // rows: 5 but NO chr22.parquet on disk → corrupt cache.
            chroms: vec![ChromEntry {
                chrom: "chr22".into(),
                file: "chr22.parquet".into(),
                rows: 5,
                warm: 0,
                cold: 5,
            }],
            cache_source_version: None,
        };
        manifest.write(&plugin_dir).unwrap();

        match PluginRegistry::open(cache_root, "22").await {
            Err(e) => assert!(
                e.to_string().contains("shard is missing"),
                "error should name the missing shard, got: {e}"
            ),
            Ok(_) => panic!("expected error on missing non-empty shard"),
        }
    }
}
