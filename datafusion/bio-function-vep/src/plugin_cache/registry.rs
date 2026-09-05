//! Runtime plugin registry: discover built plugins, open a per-chrom
//! [`PluginLookup`] for each, and drive the buffer-batched lookup — one
//! [`PluginLookup::take_buffer`] per plugin per buffer, then synchronous
//! per-transcript probes against the resulting [`PluginBufferSlice`]s.

use std::collections::HashSet;
use std::path::Path;

use datafusion::common::{DataFusionError, Result};

use crate::cache::manifest::canonical_chrom_label;
use crate::plugin_cache::cache_manifest::{
    AlleleMatch, CacheManifest, FieldOrder, discover_plugins,
};
use crate::plugin_cache::lookup::{PluginBufferSlice, PluginLookup, PluginScalar};
use crate::plugin_cache::template::CompiledTemplate;

/// One enabled plugin, with its per-chrom lookup (absent if this plugin has no
/// shard for the current chrom).
struct PluginEntry {
    csq_fields: Vec<String>,
    /// Indices into the shard's value columns, in emitted-field order.
    emit_order: Vec<usize>,
    /// Ensembl's comparison rule for this plugin; gates allele reduction.
    allele_match: AlleleMatch,
    /// The compiled discriminator template for each match column, in order.
    match_templates: Vec<CompiledTemplate>,
    n_match: usize,
    n_values: usize,
    lookup: Option<PluginLookup>,
}

/// All enabled plugins for one chromosome, in requested order (or alphabetical
/// plugin-name order when no selection was supplied).
pub struct PluginRegistry {
    plugins: Vec<PluginEntry>,
}

/// Indices into a plugin's `value_columns`, in the order its CSQ fields are
/// emitted: declaration order, or sorted by field name when the plugin's
/// Ensembl counterpart sorts its own fields ([`FieldOrder::Alphabetical`]).
fn emit_order(m: &CacheManifest) -> Vec<usize> {
    let mut order: Vec<usize> = (0..m.value_columns.len()).collect();
    if m.field_order == FieldOrder::Alphabetical {
        order.sort_by(|&a, &b| {
            m.value_columns[a]
                .csq_field
                .cmp(&m.value_columns[b].csq_field)
        });
    }
    order
}

/// Return manifests in caller order when a selection is supplied, otherwise in
/// the deterministic alphabetical order returned by [`discover_plugins`].
fn select_manifests(
    cache_root: &Path,
    plugin_names: Option<&[String]>,
) -> Result<Vec<CacheManifest>> {
    let Some(plugin_names) = plugin_names else {
        // No selection means every discovered plugin is enabled. Keep this
        // path strict: silently skipping a corrupt enabled manifest would let
        // the runtime CSQ values drift from the advertised header layout.
        return discover_plugins(cache_root);
    };
    let mut seen = HashSet::with_capacity(plugin_names.len());
    for name in plugin_names {
        if !seen.insert(name.as_str()) {
            return Err(DataFusionError::Execution(format!(
                "duplicate plugin name in selection: '{name}'"
            )));
        }
    }
    if plugin_names.is_empty() {
        return Ok(Vec::new());
    }

    let plugin_root = cache_root.join("plugin");
    let mut available = Vec::new();
    if plugin_root.exists() {
        for entry in std::fs::read_dir(&plugin_root).map_err(|e| {
            DataFusionError::Execution(format!("read {}: {e}", plugin_root.display()))
        })? {
            let path = entry
                .map_err(|e| DataFusionError::Execution(format!("dir entry: {e}")))?
                .path();
            if path.join("manifest.json").exists()
                && let Some(name) = path.file_name().and_then(|name| name.to_str())
            {
                available.push(name.to_string());
            }
        }
    }
    available.sort();

    let available_set: HashSet<&str> = available.iter().map(String::as_str).collect();
    let mut selected = Vec::with_capacity(plugin_names.len());
    for name in plugin_names {
        if !available_set.contains(name.as_str()) {
            return Err(DataFusionError::Execution(format!(
                "plugin '{name}' was requested but is not available; available plugins: {}",
                if available.is_empty() {
                    "(none)".to_string()
                } else {
                    available.join(", ")
                }
            )));
        }
        let manifest_path = plugin_root.join(name).join("manifest.json");
        let manifest = CacheManifest::read(&manifest_path)?;
        if manifest.plugin_name != *name {
            return Err(DataFusionError::Execution(format!(
                "plugin directory '{}' contains a manifest for '{}'",
                name, manifest.plugin_name
            )));
        }
        selected.push(manifest);
    }
    Ok(selected)
}

impl PluginRegistry {
    /// Discover plugins under `cache_root` and open each one's shard for `chrom`.
    pub async fn open(
        cache_root: &Path,
        chrom: &str,
        plugin_names: Option<&[String]>,
    ) -> Result<Self> {
        let want = canonical_chrom_label(chrom);
        let manifests = select_manifests(cache_root, plugin_names)?;
        let mut plugins = Vec::with_capacity(manifests.len());
        for m in manifests {
            // The emitted field names are permuted, but the shard projection is
            // NOT: `ProjectionMask::leaves` yields columns in the file's physical
            // order whatever order the leaves are listed in, so permuting the
            // projection would move the names while leaving the values put, and
            // silently pair each field with another field's value. The values are
            // permuted after probing instead (see `probe_all`).
            let order = emit_order(&m);
            let csq_fields: Vec<String> = order
                .iter()
                .map(|&i| m.value_columns[i].csq_field.clone())
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
                emit_order: order,
                allele_match: m.allele_match,
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
    /// header.
    pub fn field_names(cache_root: &Path, plugin_names: Option<&[String]>) -> Result<Vec<String>> {
        Ok(select_manifests(cache_root, plugin_names)?
            .iter()
            .flat_map(|m| {
                emit_order(m)
                    .into_iter()
                    .map(|i| m.value_columns[i].csq_field.clone())
                    .collect::<Vec<_>>()
            })
            .collect())
    }

    /// Field names from every individually readable manifest, used only to
    /// remove stale plugin declarations from an input VCF header. A malformed
    /// disabled plugin must not prevent annotation with a selected subset.
    pub fn field_names_for_cleanup(cache_root: &Path) -> Vec<String> {
        let plugin_root = cache_root.join("plugin");
        let Ok(entries) = std::fs::read_dir(&plugin_root) else {
            return Vec::new();
        };
        let mut manifests = entries
            .filter_map(|entry| entry.ok())
            .filter_map(|entry| CacheManifest::read(&entry.path().join("manifest.json")).ok())
            .collect::<Vec<_>>();
        manifests.sort_by(|a, b| a.plugin_name.cmp(&b.plugin_name));
        manifests
            .iter()
            .flat_map(|manifest| {
                emit_order(manifest)
                    .into_iter()
                    .map(|index| manifest.value_columns[index].csq_field.clone())
                    .collect::<Vec<_>>()
            })
            .collect()
    }

    /// `(csq_field, description)` for every plugin field that declares one, in
    /// emitted order — the `##<FIELD>=<description>` header lines Ensembl VEP
    /// writes for its plugin fields. Same discovery path and ordering as
    /// [`Self::field_names`], so the header cannot drift from the CSQ layout.
    pub fn field_descriptions(
        cache_root: &Path,
        plugin_names: Option<&[String]>,
    ) -> Result<Vec<(String, String)>> {
        Ok(select_manifests(cache_root, plugin_names)?
            .iter()
            .flat_map(|m| {
                emit_order(m)
                    .into_iter()
                    .filter_map(|i| {
                        let v = &m.value_columns[i];
                        v.description
                            .as_ref()
                            .map(|d| (v.csq_field.clone(), d.clone()))
                    })
                    .collect::<Vec<_>>()
            })
            .collect())
    }

    /// Concatenated CSQ field names across all plugins, in emitted order.
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
                emit_order: p.emit_order.clone(),
                allele_match: p.allele_match,
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
    emit_order: Vec<usize>,
    allele_match: AlleleMatch,
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
        fallback_key: Option<(u32, &str)>,
        attrs: &[Option<&str>],
    ) -> Vec<PluginScalar> {
        let mut out = Vec::new();
        for e in &self.entries {
            let match_values: Vec<Option<String>> =
                e.match_templates.iter().map(|t| t.eval(attrs)).collect();
            let hit = e
                .slice
                .as_ref()
                .and_then(|s| s.probe(start, allele_string, &match_values))
                .or_else(|| {
                    // Sources differ in how they spell the same variant: most key
                    // the parser-level (anchor-trimmed) allele, but a per-base
                    // source such as CADD's whole-genome SNV file keys the fully
                    // minimal one, so an untrimmed MNV misses on the primary key.
                    // Only consulted on a miss, so no existing hit can change.
                    if e.allele_match != AlleleMatch::Minimised {
                        return None;
                    }
                    let (fb_start, fb_allele) = fallback_key?;
                    e.slice
                        .as_ref()
                        .and_then(|s| s.probe(fb_start, fb_allele, &match_values))
                });
            match hit {
                // `values` arrive in shard-column order; emit them in the same
                // order as this plugin's `csq_fields`.
                Some(values) => out.extend(e.emit_order.iter().map(|&i| values[i].clone())),
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

    fn write_empty_manifest(cache_root: &Path, plugin_name: &str, csq_field: &str) {
        let plugin_dir = cache_root.join("plugin").join(plugin_name);
        std::fs::create_dir_all(&plugin_dir).unwrap();
        CacheManifest {
            plugin_name: plugin_name.into(),
            source_manifest: format!("{plugin_name}.source.toml"),
            key_columns: vec![
                "chrom".into(),
                "start".into(),
                "end".into(),
                "allele_string".into(),
            ],
            match_columns: vec![],
            value_columns: vec![ValueColumnRecord {
                column: "score".into(),
                csq_field: csq_field.into(),
                ty: "Float32".into(),
                description: Some(format!("{plugin_name} description")),
            }],
            chroms: vec![],
            sources: vec![],
            cache_source_version: None,
            allele_match: Default::default(),
            field_order: Default::default(),
        }
        .write(&plugin_dir)
        .unwrap();
    }

    #[test]
    fn selection_preserves_caller_order_and_validates_names() {
        let dir = tempfile::tempdir().unwrap();
        write_empty_manifest(dir.path(), "zeta", "ZETA");
        write_empty_manifest(dir.path(), "alpha", "ALPHA");

        assert_eq!(
            PluginRegistry::field_names(dir.path(), None).unwrap(),
            vec!["ALPHA", "ZETA"]
        );

        let requested = vec!["zeta".to_string(), "alpha".to_string()];
        assert_eq!(
            PluginRegistry::field_names(dir.path(), Some(&requested)).unwrap(),
            vec!["ZETA", "ALPHA"]
        );
        assert_eq!(
            PluginRegistry::field_descriptions(dir.path(), Some(&requested)).unwrap(),
            vec![
                ("ZETA".to_string(), "zeta description".to_string()),
                ("ALPHA".to_string(), "alpha description".to_string()),
            ]
        );

        let duplicate = vec!["alpha".to_string(), "alpha".to_string()];
        let error = PluginRegistry::field_names(dir.path(), Some(&duplicate))
            .unwrap_err()
            .to_string();
        assert!(error.contains("duplicate plugin name"), "{error}");

        let missing = vec!["missing".to_string()];
        let error = PluginRegistry::field_names(dir.path(), Some(&missing))
            .unwrap_err()
            .to_string();
        assert!(error.contains("available plugins: alpha, zeta"), "{error}");
    }

    #[tokio::test]
    async fn selection_does_not_validate_disabled_plugin_manifests() {
        let dir = tempfile::tempdir().unwrap();
        write_empty_manifest(dir.path(), "enabled", "ENABLED");
        let broken_dir = dir.path().join("plugin").join("broken");
        std::fs::create_dir_all(&broken_dir).unwrap();
        std::fs::write(broken_dir.join("manifest.json"), "not valid JSON").unwrap();

        let selected = vec!["enabled".to_string()];
        assert_eq!(
            PluginRegistry::field_names(dir.path(), Some(&selected)).unwrap(),
            vec!["ENABLED"]
        );
        assert_eq!(
            PluginRegistry::field_descriptions(dir.path(), Some(&selected)).unwrap(),
            vec![("ENABLED".to_string(), "enabled description".to_string())]
        );
        assert!(
            !PluginRegistry::open(dir.path(), "22", Some(&selected))
                .await
                .unwrap()
                .is_empty()
        );

        let disabled = Vec::new();
        assert!(
            PluginRegistry::field_names(dir.path(), Some(&disabled))
                .unwrap()
                .is_empty()
        );
        assert!(
            PluginRegistry::open(dir.path(), "22", Some(&disabled))
                .await
                .unwrap()
                .is_empty()
        );
        assert_eq!(
            PluginRegistry::field_names_for_cleanup(dir.path()),
            vec!["ENABLED"]
        );
        // The unfiltered mode enables every plugin and therefore remains
        // intentionally strict about every manifest.
        assert!(PluginRegistry::field_names(dir.path(), None).is_err());
    }

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
            description: None,
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
                description: None,
            }],
            chroms: vec![ChromEntry {
                chrom: "chr22".into(),
                file: "chr22.parquet".into(),
                rows: 2,
                warm: 0,
                cold: 2,
            }],
            sources: vec![],
            cache_source_version: None,
            allele_match: Default::default(),
            field_order: Default::default(),
        };
        manifest.write(&plugin_dir).unwrap();

        let reg = PluginRegistry::open(cache_root, "22", None).await.unwrap();
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
            "",
            "78",
            "R/G",
            "",
            "A",
            "G",
        );
        let hit = slices.probe_all(100, "A/G", None, &ns_hit);
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
            "",
            "A",
            "G",
        );
        let none = slices.probe_all(100, "A/G", None, &ns_miss);
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
                description: None,
            }],
            // chr22 present in the manifest with rows: 0 and NO file on disk.
            chroms: vec![ChromEntry {
                chrom: "chr22".into(),
                file: "chr22.parquet".into(),
                rows: 0,
                warm: 0,
                cold: 0,
            }],
            sources: vec![],
            cache_source_version: None,
            allele_match: Default::default(),
            field_order: Default::default(),
        };
        manifest.write(&plugin_dir).unwrap();

        let reg = PluginRegistry::open(cache_root, "22", None).await.unwrap();
        assert_eq!(reg.csq_fields(), vec!["SCORE".to_string()]);
        let slices = reg.take_buffer_all(&[100]).await.unwrap();
        // No shard → empty (Null) field, not an error.
        assert_eq!(
            slices.probe_all(100, "A/G", None, &[]),
            vec![PluginScalar::Null]
        );
    }

    // A per-variant plugin (no match_column) is keyed only on (start,
    // allele_string), so an empty-namespace probe MUST hit — this is what the
    // no-transcript placeholder path relies on to emit per-variant values on
    // intergenic variants (Codex PR #190 per-variant placeholder finding).
    #[tokio::test(flavor = "multi_thread")]
    async fn per_variant_plugin_hits_with_empty_namespace() {
        let dir = tempfile::tempdir().unwrap();
        let cache_root = dir.path();
        let plugin_dir = cache_root.join("plugin").join("demo");
        std::fs::create_dir_all(&plugin_dir).unwrap();

        // Per-variant shard: no match columns, one value column.
        let vals = vec![ValueColumn {
            column: "score".into(),
            csq_field: "SCORE".into(),
            ty: ValueType::Float32,
            description: None,
        }];
        let schema = plugin_output_schema(&[], &vals);
        let batch = RecordBatch::try_new(
            schema.clone(),
            vec![
                Arc::new(StringArray::from(vec!["22"])),
                Arc::new(UInt32Array::from(vec![100u32])),
                Arc::new(UInt32Array::from(vec![100u32])),
                Arc::new(StringArray::from(vec!["A/G"])),
                Arc::new(Float32Array::from(vec![0.5f32])),
                Arc::new(Int8Array::from(vec![1i8])),
            ],
        )
        .unwrap();
        let mut w = PluginShardWriter::create(&plugin_dir.join("chr22.parquet"), schema).unwrap();
        w.write(&batch).unwrap();
        w.finish().unwrap();

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
                description: None,
            }],
            chroms: vec![ChromEntry {
                chrom: "chr22".into(),
                file: "chr22.parquet".into(),
                rows: 1,
                warm: 0,
                cold: 1,
            }],
            sources: vec![],
            cache_source_version: None,
            allele_match: Default::default(),
            field_order: Default::default(),
        };
        manifest.write(&plugin_dir).unwrap();

        let reg = PluginRegistry::open(cache_root, "22", None).await.unwrap();
        let slices = reg.take_buffer_all(&[100]).await.unwrap();
        // Empty namespace (no transcript) still hits the per-variant row.
        match slices.probe_all(100, "A/G", None, &[])[0] {
            PluginScalar::F32(v) => assert!((v - 0.5).abs() < 1e-6),
            ref other => panic!("{other:?}"),
        }
        // Wrong allele still misses.
        assert_eq!(
            slices.probe_all(100, "C/T", None, &[]),
            vec![PluginScalar::Null]
        );
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
                description: None,
            }],
            // rows: 5 but NO chr22.parquet on disk → corrupt cache.
            chroms: vec![ChromEntry {
                chrom: "chr22".into(),
                file: "chr22.parquet".into(),
                rows: 5,
                warm: 0,
                cold: 5,
            }],
            sources: vec![],
            cache_source_version: None,
            allele_match: Default::default(),
            field_order: Default::default(),
        };
        manifest.write(&plugin_dir).unwrap();

        match PluginRegistry::open(cache_root, "22", None).await {
            Err(e) => assert!(
                e.to_string().contains("shard is missing"),
                "error should name the missing shard, got: {e}"
            ),
            Ok(_) => panic!("expected error on missing non-empty shard"),
        }
    }
}
