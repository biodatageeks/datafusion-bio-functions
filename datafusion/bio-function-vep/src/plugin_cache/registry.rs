//! Runtime plugin registry: discover built plugins, open a per-chrom
//! [`PluginLookup`] for each, and drive the buffer-batched lookup — one
//! [`PluginLookup::take_buffer`] per plugin per buffer, then synchronous
//! per-transcript probes against the resulting [`PluginBufferSlice`]s.

use std::collections::HashSet;
use std::path::Path;
use std::sync::{Arc, OnceLock};
use std::time::Instant;

use datafusion::common::{DataFusionError, Result};

use crate::cache::manifest::canonical_chrom_label;
use crate::plugin_cache::cache_manifest::{
    AlleleMatch, CacheManifest, FieldOrder, LookupKind, discover_plugins,
};
use crate::plugin_cache::lookup::{
    IntervalLookup, PluginBufferSlice, PluginLookup, PluginScalar, TakeStats,
};
use crate::plugin_cache::template::CompiledTemplate;

/// The per-chrom shard handle of one plugin, by lookup kind.
enum LookupHandle {
    /// `Arc` so each buffer's take can run as its own spawned task.
    Point(Arc<PluginLookup>),
    Interval(Arc<IntervalLookup>),
}

/// One enabled plugin, with its per-chrom lookup (absent if this plugin has no
/// shard for the current chrom).
struct PluginEntry {
    /// Plugin name, for the `VEP_ENGINE_PROFILE` per-plugin take lines.
    name: String,
    csq_fields: Vec<String>,
    /// Indices into the shard's value columns, in emitted-field order.
    emit_order: Vec<usize>,
    /// Ensembl's comparison rule for this plugin; gates allele reduction.
    allele_match: AlleleMatch,
    /// The compiled discriminator template for each match column, in order.
    match_templates: Vec<CompiledTemplate>,
    n_match: usize,
    n_values: usize,
    lookup: Option<LookupHandle>,
}

/// All enabled plugins for one chromosome, in requested order (or alphabetical
/// plugin-name order when no selection was supplied).
pub struct PluginRegistry {
    plugins: Vec<PluginEntry>,
}

/// True when `VEP_ENGINE_PROFILE` is set; read once per process, since
/// [`PluginRegistry::take_buffer_all`] runs once per buffer.
fn plugin_profile_enabled() -> bool {
    static ENABLED: OnceLock<bool> = OnceLock::new();
    *ENABLED.get_or_init(|| std::env::var_os("VEP_ENGINE_PROFILE").is_some())
}

/// One `[VEP_PLUGIN_PROFILE]` line: a single plugin's take for one buffer.
/// `elapsed` is the take's own wall time; the split fields come from
/// [`TakeStats`].
fn plugin_profile_line(
    name: &str,
    probes: usize,
    elapsed: std::time::Duration,
    stats: &TakeStats,
) -> String {
    format!(
        "[VEP_PLUGIN_PROFILE] plugin={name} probes={probes} rows={} elapsed={:.6}s bytes={} ops={} open={:.6}s starts_read={:.6}s payload_read={:.6}s",
        stats.rows,
        elapsed.as_secs_f64(),
        stats.bytes,
        stats.ops,
        stats.open.as_secs_f64(),
        stats.starts_read.as_secs_f64(),
        stats.payload_read.as_secs_f64(),
    )
}

/// The process-wide pool the per-buffer plugin takes run on, created on first
/// use: a multi-thread runtime of `VEP_PLUGIN_TAKE_THREADS` workers (default
/// `min(4, available_parallelism)`), threads named `vep-plugin-take`.
///
/// Dedicated rather than the caller's runtime because the caller holds its
/// own worker while it waits (`lookup_exec::block_on` is `block_in_place` +
/// `block_on`): on a one-worker runtime the takes would still run one after
/// another. Tasks on this pool only read shards; they never call back into
/// the caller's runtime, so awaiting them from any context cannot deadlock.
fn plugin_take_pool() -> Result<&'static tokio::runtime::Handle> {
    static POOL: OnceLock<std::result::Result<tokio::runtime::Runtime, String>> = OnceLock::new();
    POOL.get_or_init(|| {
        let threads = std::env::var("VEP_PLUGIN_TAKE_THREADS")
            .ok()
            .and_then(|v| v.trim().parse::<usize>().ok())
            .filter(|&n| n > 0)
            .unwrap_or_else(|| {
                std::thread::available_parallelism()
                    .map(|n| n.get())
                    .unwrap_or(1)
                    .min(4)
            });
        tokio::runtime::Builder::new_multi_thread()
            .worker_threads(threads)
            .thread_name("vep-plugin-take")
            .enable_all()
            .build()
            .map_err(|e| e.to_string())
    })
    .as_ref()
    .map(|rt| rt.handle())
    .map_err(|e| DataFusionError::Execution(format!("build plugin take pool: {e}")))
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
    /// Validate the selected manifest set without opening any chromosome shard.
    /// Used during planning so invalid configuration is rejected even when the
    /// input has no data-bearing contigs.
    pub fn validate_selection(cache_root: &Path, plugin_names: Option<&[String]>) -> Result<()> {
        select_manifests(cache_root, plugin_names).map(|_| ())
    }

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
                    Some(match m.lookup {
                        LookupKind::Point => LookupHandle::Point(Arc::new(
                            PluginLookup::open(&shard, match_columns, value_columns).await?,
                        )),
                        LookupKind::Interval => LookupHandle::Interval(Arc::new(
                            IntervalLookup::open(&shard, match_columns, value_columns).await?,
                        )),
                    })
                }
                _ => None,
            };
            plugins.push(PluginEntry {
                name: m.plugin_name.clone(),
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
    ///
    /// The Point takes run concurrently, one task per plugin on the dedicated
    /// [`plugin_take_pool`], so their Parquet decode runs in parallel whatever
    /// the caller's runtime (a one-worker runtime included) instead of back
    /// to back on the caller's thread. Each task also builds its slice.
    /// Results are awaited in plugin order, so the entries, and on failure
    /// the error returned (the first failing plugin in that order), are the
    /// same as a serial loop's.
    pub async fn take_buffer_all(&self, sorted_unique_starts: &[u32]) -> Result<BufferSlices> {
        let starts: Arc<[u32]> = Arc::from(sorted_unique_starts);
        let pool = plugin_take_pool()?;
        let mut tasks = AbortOnDrop(
            self.plugins
                .iter()
                .map(|p| match &p.lookup {
                    Some(LookupHandle::Point(lk)) => {
                        let lk = Arc::clone(lk);
                        let starts = Arc::clone(&starts);
                        let (n_match, n_values) = (p.n_match, p.n_values);
                        Some(pool.spawn(async move {
                            let started = Instant::now();
                            let (batch, stats) = lk.take_buffer_with_stats(&starts).await?;
                            let elapsed = started.elapsed();
                            let slice = PluginBufferSlice::from_batch(&batch, n_match, n_values)?;
                            Ok::<_, DataFusionError>((slice, stats, elapsed))
                        }))
                    }
                    _ => None,
                })
                .collect(),
        );
        let mut entries = Vec::with_capacity(self.plugins.len());
        for (i, p) in self.plugins.iter().enumerate() {
            let (slice, interval) = match &p.lookup {
                Some(LookupHandle::Point(_)) => {
                    let handle = tasks.0[i]
                        .as_mut()
                        .expect("a take task is spawned for every Point plugin");
                    let (slice, stats, elapsed) = handle.await.map_err(|e| {
                        DataFusionError::Execution(format!(
                            "plugin '{}' take task failed: {e}",
                            p.name
                        ))
                    })??;
                    tasks.0[i] = None;
                    if plugin_profile_enabled() {
                        eprintln!(
                            "{}",
                            plugin_profile_line(&p.name, starts.len(), elapsed, &stats)
                        );
                    }
                    (Some(slice), None)
                }
                // Interval shards are resident for the whole contig; nothing
                // to take per buffer.
                Some(LookupHandle::Interval(il)) => (None, Some(Arc::clone(il))),
                None => (None, None),
            };
            entries.push(SliceEntry {
                csq_fields_len: p.csq_fields.len(),
                emit_order: p.emit_order.clone(),
                allele_match: p.allele_match,
                match_templates: p.match_templates.clone(),
                slice,
                interval,
            });
        }
        Ok(BufferSlices { entries })
    }
}

/// The per-plugin take tasks of one [`PluginRegistry::take_buffer_all`] call,
/// by plugin index. Aborts whatever is still pending when dropped, so an
/// early error return, or a caller dropping the future, does not leave shard
/// reads running detached.
struct AbortOnDrop<T>(Vec<Option<tokio::task::JoinHandle<T>>>);

impl<T> Drop for AbortOnDrop<T> {
    fn drop(&mut self) {
        for handle in self.0.iter().flatten() {
            handle.abort();
        }
    }
}

/// Per-plugin buffer slice plus the metadata `probe_all` needs.
struct SliceEntry {
    csq_fields_len: usize,
    emit_order: Vec<usize>,
    allele_match: AlleleMatch,
    match_templates: Vec<CompiledTemplate>,
    slice: Option<PluginBufferSlice>,
    interval: Option<Arc<IntervalLookup>>,
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
    /// field (the per-transcript gate for match-column plugins). `span` is the
    /// variant's VEP-normalised inclusive `[start, end]` with `span.0 <= span.1`
    /// (the caller swaps insertion coordinates); interval plugins probe by it,
    /// and `None` (a symbolic allele, which Ensembl never hands to a plugin)
    /// makes every interval plugin miss.
    pub fn probe_all(
        &self,
        start: u32,
        allele_string: &str,
        fallback_key: Option<(u32, &str)>,
        span: Option<(u32, u32)>,
        attrs: &[Option<&str>],
    ) -> Vec<PluginScalar> {
        let mut out = Vec::new();
        for e in &self.entries {
            let match_values: Vec<Option<String>> =
                e.match_templates.iter().map(|t| t.eval(attrs)).collect();
            if let Some(il) = &e.interval {
                match span.and_then(|(lo, hi)| il.probe(lo, hi, &match_values)) {
                    Some(values) => out.extend(e.emit_order.iter().map(|&i| values[i].clone())),
                    None => out.extend(std::iter::repeat_n(PluginScalar::Null, e.csq_fields_len)),
                }
                continue;
            }
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
            assume_unique: None,
            lookup: Default::default(),
        }
        .write(&plugin_dir)
        .unwrap();
    }

    #[tokio::test(flavor = "multi_thread")]
    async fn interval_plugin_probes_by_span_and_gene() {
        use crate::plugin_cache::cache_manifest::key_columns;
        use crate::plugin_cache::lookup::test_support::write_interval_shard;
        let dir = tempfile::tempdir().unwrap();
        let cache_root = dir.path();
        let plugin_dir = cache_root.join("plugin").join("po");
        std::fs::create_dir_all(&plugin_dir).unwrap();
        write_interval_shard(&plugin_dir.join("chr1.parquet"));
        CacheManifest {
            plugin_name: "po".into(),
            source_manifest: "po.source.toml".into(),
            key_columns: key_columns(LookupKind::Interval),
            match_columns: vec![MatchColumnRecord {
                column: "gene_id".into(),
                template: "{Gene}".into(),
            }],
            value_columns: vec![ValueColumnRecord {
                column: "rat".into(),
                csq_field: "PO_Rat".into(),
                ty: "Utf8".into(),
                description: None,
            }],
            chroms: vec![ChromEntry {
                chrom: "chr1".into(),
                file: "chr1.parquet".into(),
                rows: 4,
                warm: 0,
                cold: 4,
            }],
            sources: vec![],
            cache_source_version: None,
            allele_match: Default::default(),
            field_order: Default::default(),
            assume_unique: Some(true),
            lookup: LookupKind::Interval,
        }
        .write(&plugin_dir)
        .unwrap();

        let reg = PluginRegistry::open(cache_root, "1", None).await.unwrap();
        assert_eq!(reg.csq_fields(), vec!["PO_Rat"]);
        // Interval plugins ignore the buffer take; an empty start list is fine.
        let slices = reg.take_buffer_all(&[]).await.unwrap();
        let ns = build_attr_namespace(
            "intron_variant",
            "ENSG1",
            "SYM",
            "Transcript",
            "ENST1",
            "lncRNA",
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
        assert_eq!(
            slices.probe_all(160, "A/G", None, Some((160, 160)), &ns),
            vec![PluginScalar::Str("wide".into())]
        );
        assert_eq!(
            slices.probe_all(350, "A/G", None, Some((350, 350)), &ns),
            vec![PluginScalar::Str("late".into())]
        );
        assert_eq!(
            slices.probe_all(351, "A/G", None, Some((351, 351)), &ns),
            vec![PluginScalar::Null]
        );
        let other = build_attr_namespace(
            "intron_variant",
            "ENSG2",
            "SYM",
            "Transcript",
            "ENST2",
            "lncRNA",
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
        assert_eq!(
            slices.probe_all(160, "A/G", None, Some((160, 160)), &other),
            vec![PluginScalar::Null]
        );
        // Symbolic allele → no span → miss even with a matching gene.
        assert_eq!(
            slices.probe_all(160, "A/<DEL>", None, None, &ns),
            vec![PluginScalar::Null]
        );
        // No transcript → empty namespace → discriminator None → miss.
        assert_eq!(
            slices.probe_all(160, "A/G", None, Some((160, 160)), &[]),
            vec![PluginScalar::Null]
        );
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
        assert!(PluginRegistry::validate_selection(dir.path(), Some(&duplicate)).is_err());
        let error = PluginRegistry::field_names(dir.path(), Some(&duplicate))
            .unwrap_err()
            .to_string();
        assert!(error.contains("duplicate plugin name"), "{error}");

        let missing = vec!["missing".to_string()];
        assert!(PluginRegistry::validate_selection(dir.path(), Some(&missing)).is_err());
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
        PluginRegistry::validate_selection(dir.path(), Some(&selected)).unwrap();
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
        PluginRegistry::validate_selection(dir.path(), Some(&disabled)).unwrap();
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
        let broken = vec!["broken".to_string()];
        assert!(PluginRegistry::validate_selection(dir.path(), Some(&broken)).is_err());
        // The unfiltered mode enables every plugin and therefore remains
        // intentionally strict about every manifest.
        assert!(PluginRegistry::validate_selection(dir.path(), None).is_err());
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
        let schema = plugin_output_schema(LookupKind::Point, &matches, &vals);
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
            assume_unique: None,
            lookup: Default::default(),
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
        let hit = slices.probe_all(100, "A/G", None, Some((100, 100)), &ns_hit);
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
        let none = slices.probe_all(100, "A/G", None, Some((100, 100)), &ns_miss);
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
            assume_unique: None,
            lookup: Default::default(),
        };
        manifest.write(&plugin_dir).unwrap();

        let reg = PluginRegistry::open(cache_root, "22", None).await.unwrap();
        assert_eq!(reg.csq_fields(), vec!["SCORE".to_string()]);
        let slices = reg.take_buffer_all(&[100]).await.unwrap();
        // No shard → empty (Null) field, not an error.
        assert_eq!(
            slices.probe_all(100, "A/G", None, Some((100, 100)), &[]),
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
        let schema = plugin_output_schema(LookupKind::Point, &[], &vals);
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
            assume_unique: None,
            lookup: Default::default(),
        };
        manifest.write(&plugin_dir).unwrap();

        let reg = PluginRegistry::open(cache_root, "22", None).await.unwrap();
        let slices = reg.take_buffer_all(&[100]).await.unwrap();
        // Empty namespace (no transcript) still hits the per-variant row.
        match slices.probe_all(100, "A/G", None, Some((100, 100)), &[])[0] {
            PluginScalar::F32(v) => assert!((v - 0.5).abs() < 1e-6),
            ref other => panic!("{other:?}"),
        }
        // Wrong allele still misses.
        assert_eq!(
            slices.probe_all(100, "C/T", None, Some((100, 100)), &[]),
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
            assume_unique: None,
            lookup: Default::default(),
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

    /// A Point plugin shard of `n` rows at `base + i * stride`, two alleles
    /// per position, with an optional per-transcript match column and a value
    /// column of type `ty`. Rows are written sorted by `(tier, start)`.
    fn write_point_plugin(
        cache_root: &Path,
        name: &str,
        base: u32,
        stride: u32,
        n: usize,
        with_match: bool,
        ty: ValueType,
    ) {
        use datafusion::arrow::array::{ArrayRef, Int32Array};
        let plugin_dir = cache_root.join("plugin").join(name);
        std::fs::create_dir_all(&plugin_dir).unwrap();
        let matches: Vec<MatchColumn> = if with_match {
            vec![MatchColumn {
                column: "protein_variant".into(),
                template: "{ref_aa}{Protein_position}{alt_aa}".into(),
            }]
        } else {
            vec![]
        };
        let value_col = format!("{name}_score");
        let vals = vec![ValueColumn {
            column: value_col.clone(),
            csq_field: value_col.clone(),
            ty,
            description: None,
        }];
        let schema = plugin_output_schema(LookupKind::Point, &matches, &vals);
        let rows = 2 * n;
        let starts: Vec<u32> = (0..rows).map(|r| base + (r / 2) as u32 * stride).collect();
        let alleles: Vec<&str> = (0..rows)
            .map(|r| if r % 2 == 0 { "A/G" } else { "A/T" })
            .collect();
        let mut columns: Vec<ArrayRef> = vec![
            Arc::new(StringArray::from(vec!["22"; rows])),
            Arc::new(UInt32Array::from(starts.clone())),
            Arc::new(UInt32Array::from(starts)),
            Arc::new(StringArray::from(alleles)),
        ];
        if with_match {
            columns.push(Arc::new(StringArray::from(
                (0..rows).map(|r| format!("R{r}G")).collect::<Vec<_>>(),
            )));
        }
        columns.push(match ty {
            ValueType::Float32 => Arc::new(Float32Array::from(
                (0..rows).map(|r| r as f32 * 0.001).collect::<Vec<_>>(),
            )),
            ValueType::Int32 => Arc::new(Int32Array::from(
                (0..rows).map(|r| r as i32).collect::<Vec<_>>(),
            )),
            ValueType::Utf8 => Arc::new(StringArray::from(
                (0..rows).map(|r| format!("{name}-{r}")).collect::<Vec<_>>(),
            )),
        });
        columns.push(Arc::new(Int8Array::from(vec![1i8; rows])));
        let batch = RecordBatch::try_new(schema.clone(), columns).unwrap();
        let mut w = PluginShardWriter::create(&plugin_dir.join("chr22.parquet"), schema).unwrap();
        w.write(&batch).unwrap();
        w.finish().unwrap();

        CacheManifest {
            plugin_name: name.into(),
            source_manifest: format!("{name}.source.toml"),
            key_columns: vec![
                "chrom".into(),
                "start".into(),
                "end".into(),
                "allele_string".into(),
            ],
            match_columns: matches
                .iter()
                .map(|m| MatchColumnRecord {
                    column: m.column.clone(),
                    template: m.template.clone(),
                })
                .collect(),
            value_columns: vec![ValueColumnRecord {
                column: value_col.clone(),
                csq_field: value_col,
                ty: format!("{ty:?}"),
                description: None,
            }],
            chroms: vec![ChromEntry {
                chrom: "chr22".into(),
                file: "chr22.parquet".into(),
                rows,
                warm: 0,
                cold: rows,
            }],
            sources: vec![],
            cache_source_version: None,
            allele_match: Default::default(),
            field_order: Default::default(),
            assume_unique: None,
            lookup: Default::default(),
        }
        .write(&plugin_dir)
        .unwrap();
    }

    /// Four plugins over chr22 — three Point shards with different densities,
    /// match layouts and value types, plus one with no chr22 shard — and the
    /// buffers of sorted, deduped starts to take for them.
    fn multi_plugin_fixture() -> (tempfile::TempDir, Vec<String>, Vec<Vec<u32>>) {
        let dir = tempfile::tempdir().unwrap();
        let root = dir.path();
        write_point_plugin(root, "dense", 1_000, 1, 6_000, false, ValueType::Float32);
        write_point_plugin(root, "sparse", 1_000, 7, 3_000, true, ValueType::Utf8);
        write_point_plugin(root, "ints", 2_500, 3, 4_000, false, ValueType::Int32);
        write_empty_manifest(root, "absent", "absent_score");
        // Caller order, deliberately not alphabetical.
        let names: Vec<String> = ["sparse", "absent", "dense", "ints"]
            .iter()
            .map(|s| s.to_string())
            .collect();
        let mut buffers: Vec<Vec<u32>> = vec![
            vec![],
            (1_000..1_050).collect(),
            (0..400).map(|i| 1_000 + i * 13).collect(),
            vec![1, 999, 30_000, 40_000],
            (0..2_000).map(|i| 2_000 + i * 5).collect(),
            vec![6_999, 7_000, 14_503, 21_997, 22_000],
        ];
        for b in &mut buffers {
            b.sort_unstable();
            b.dedup();
        }
        (dir, names, buffers)
    }

    /// Assert that [`PluginRegistry::take_buffer_all`]'s entries equal, per
    /// plugin and in plugin order, a serial [`PluginLookup::take_buffer`] per
    /// plugin. Returns the total number of rows taken across Point plugins.
    async fn assert_all_matches_serial(
        reg: &PluginRegistry,
        slices: &BufferSlices,
        starts: &[u32],
    ) -> usize {
        assert_eq!(slices.entries.len(), reg.plugins.len());
        let mut total = 0;
        for (p, e) in reg.plugins.iter().zip(&slices.entries) {
            assert_eq!(e.csq_fields_len, p.csq_fields.len(), "plugin {}", p.name);
            assert_eq!(e.emit_order, p.emit_order, "plugin {}", p.name);
            assert_eq!(e.allele_match, p.allele_match, "plugin {}", p.name);
            assert_eq!(
                e.match_templates.len(),
                p.match_templates.len(),
                "plugin {}",
                p.name
            );
            match &p.lookup {
                Some(LookupHandle::Point(lk)) => {
                    let serial = PluginBufferSlice::from_batch(
                        &lk.take_buffer(starts).await.unwrap(),
                        p.n_match,
                        p.n_values,
                    )
                    .unwrap();
                    total += serial.len();
                    assert_eq!(e.slice.as_ref(), Some(&serial), "plugin {}", p.name);
                    assert!(e.interval.is_none());
                }
                Some(LookupHandle::Interval(_)) => panic!("fixture has no interval plugin"),
                None => {
                    assert!(e.slice.is_none(), "plugin {}", p.name);
                    assert!(e.interval.is_none(), "plugin {}", p.name);
                }
            }
        }
        total
    }

    async fn run_multi_plugin_buffers(via_block_on: bool) {
        let (dir, names, buffers) = multi_plugin_fixture();
        let reg = PluginRegistry::open(dir.path(), "22", Some(&names))
            .await
            .unwrap();
        assert_eq!(
            reg.plugins
                .iter()
                .map(|p| p.name.as_str())
                .collect::<Vec<_>>(),
            vec!["sparse", "absent", "dense", "ints"]
        );
        let mut total = 0;
        for starts in &buffers {
            let slices = if via_block_on {
                crate::cache::lookup_exec::block_on(reg.take_buffer_all(starts)).unwrap()
            } else {
                reg.take_buffer_all(starts).await.unwrap()
            };
            total += assert_all_matches_serial(&reg, &slices, starts).await;
        }
        // The buffers must actually hit rows, or the equality proves nothing.
        assert!(total > 1_000, "only {total} rows taken");
    }

    #[tokio::test]
    async fn take_buffer_all_matches_serial_current_thread() {
        run_multi_plugin_buffers(false).await;
        // Through the engine's `block_on`, which on a current-thread runtime
        // drives the take on a fresh runtime in a scoped thread.
        run_multi_plugin_buffers(true).await;
    }

    #[tokio::test(flavor = "multi_thread")]
    async fn take_buffer_all_matches_serial_multi_thread() {
        run_multi_plugin_buffers(false).await;
        // `block_on` here is `block_in_place` + `handle.block_on`, the path
        // the annotation provider takes.
        run_multi_plugin_buffers(true).await;
    }

    #[test]
    fn take_buffer_all_matches_serial_without_runtime() {
        let (dir, names, buffers) = multi_plugin_fixture();
        let rt = tokio::runtime::Runtime::new().unwrap();
        let reg = rt
            .block_on(PluginRegistry::open(dir.path(), "22", Some(&names)))
            .unwrap();
        drop(rt);
        // No runtime on this thread: `block_on` builds its own.
        let all: Vec<BufferSlices> = buffers
            .iter()
            .map(|starts| crate::cache::lookup_exec::block_on(reg.take_buffer_all(starts)).unwrap())
            .collect();
        let rt = tokio::runtime::Runtime::new().unwrap();
        let mut total = 0;
        for (starts, slices) in buffers.iter().zip(&all) {
            total += rt.block_on(assert_all_matches_serial(&reg, slices, starts));
        }
        assert!(total > 1_000, "only {total} rows taken");
    }

    /// With the takes running concurrently, the error returned is still the
    /// first failing plugin's in plugin order, as with a serial loop.
    #[tokio::test(flavor = "multi_thread")]
    async fn take_buffer_all_returns_first_error_in_plugin_order() {
        let (dir, names, _) = multi_plugin_fixture();
        let reg = PluginRegistry::open(dir.path(), "22", Some(&names))
            .await
            .unwrap();
        // Shards are reopened per take, so removing them after open makes
        // those plugins' takes fail. "dense" precedes "ints" in plugin order.
        for name in ["ints", "dense"] {
            std::fs::remove_file(dir.path().join("plugin").join(name).join("chr22.parquet"))
                .unwrap();
        }
        let starts: Vec<u32> = (1_000..1_050).collect();
        for _ in 0..20 {
            let err = match reg.take_buffer_all(&starts).await {
                Err(e) => e.to_string(),
                Ok(_) => panic!("expected the take to fail"),
            };
            assert!(err.contains("/dense/chr22.parquet"), "{err}");
        }
    }

    /// vepyr's `workers=1` setup: a one-worker multi-thread runtime whose only
    /// worker is held by the annotation task while it `block_on`s the take.
    fn one_worker_runtime() -> tokio::runtime::Runtime {
        tokio::runtime::Builder::new_multi_thread()
            .worker_threads(1)
            .enable_all()
            .build()
            .unwrap()
    }

    #[test]
    fn take_buffer_all_matches_serial_one_worker_runtime() {
        let rt = one_worker_runtime();
        rt.block_on(async {
            // Spawned, so the take is driven from the runtime's only worker,
            // as the annotation stream is.
            tokio::spawn(run_multi_plugin_buffers(true)).await.unwrap();
            tokio::spawn(run_multi_plugin_buffers(false)).await.unwrap();
        });
    }

    /// The pre-pool behaviour: every Point take awaited in turn on the
    /// calling task. Reference for the bench only.
    async fn take_all_serial(reg: &PluginRegistry, starts: &[u32]) -> usize {
        let mut rows = 0;
        for p in &reg.plugins {
            if let Some(LookupHandle::Point(lk)) = &p.lookup {
                let (batch, _) = lk.take_buffer_with_stats(starts).await.unwrap();
                rows += PluginBufferSlice::from_batch(&batch, p.n_match, p.n_values)
                    .unwrap()
                    .len();
            }
        }
        rows
    }

    /// Real-cache timing probe, not a correctness test. Set
    /// `VEP_PLUGIN_BENCH_ROOT` (a plugin cache root), `VEP_PLUGIN_BENCH_CHROM`
    /// and `VEP_PLUGIN_BENCH_STARTS` (a TSV of `buffer<TAB>start`), plus
    /// `VEP_ENGINE_PROFILE=1` for the per-plugin lines, then run with
    /// `--ignored --nocapture`, ideally in `--release`. Runs in vepyr's
    /// `workers=1` setup: a one-worker runtime whose worker is held by the
    /// task that `block_on`s the take. Compares the serial takes with
    /// `take_buffer_all`, `VEP_PLUGIN_BENCH_REPS` times each (default 3).
    #[test]
    #[ignore = "needs a real plugin cache; set VEP_PLUGIN_BENCH_*"]
    fn bench_take_buffer_all_real_cache() {
        let root = std::env::var("VEP_PLUGIN_BENCH_ROOT").expect("VEP_PLUGIN_BENCH_ROOT");
        let chrom = std::env::var("VEP_PLUGIN_BENCH_CHROM").unwrap_or_else(|_| "22".into());
        let tsv = std::env::var("VEP_PLUGIN_BENCH_STARTS").expect("VEP_PLUGIN_BENCH_STARTS");
        let reps: usize = std::env::var("VEP_PLUGIN_BENCH_REPS")
            .ok()
            .and_then(|v| v.parse().ok())
            .unwrap_or(3);
        let mut buffers: std::collections::BTreeMap<u32, Vec<u32>> = Default::default();
        for line in std::fs::read_to_string(tsv).unwrap().lines() {
            let (b, s) = line.split_once('\t').unwrap();
            buffers
                .entry(b.parse().unwrap())
                .or_default()
                .push(s.parse().unwrap());
        }
        let buffers: Vec<Vec<u32>> = buffers
            .into_values()
            .map(|mut b| {
                b.sort_unstable();
                b.dedup();
                b
            })
            .collect();
        let rt = one_worker_runtime();
        let reg = Arc::new(
            rt.block_on(PluginRegistry::open(Path::new(&root), &chrom, None))
                .unwrap(),
        );
        let n_buffers = buffers.len();
        let buffers = Arc::new(buffers);
        for rep in 0..reps {
            for mode in ["serial", "take_buffer_all"] {
                let (reg, buffers) = (Arc::clone(&reg), Arc::clone(&buffers));
                let (secs, rows) = rt.block_on(async move {
                    tokio::spawn(async move {
                        let started = Instant::now();
                        let mut rows = 0;
                        for starts in buffers.iter() {
                            rows += if mode == "serial" {
                                crate::cache::lookup_exec::block_on(async {
                                    Ok(take_all_serial(&reg, starts).await)
                                })
                                .unwrap()
                            } else {
                                let slices = crate::cache::lookup_exec::block_on(
                                    reg.take_buffer_all(starts),
                                )
                                .unwrap();
                                slices
                                    .entries
                                    .iter()
                                    .filter_map(|e| e.slice.as_ref())
                                    .map(PluginBufferSlice::len)
                                    .sum()
                            };
                        }
                        (started.elapsed().as_secs_f64(), rows)
                    })
                    .await
                    .unwrap()
                });
                eprintln!(
                    "[VEP_PLUGIN_BENCH] rep={rep} runtime=1-worker mode={mode} buffers={} rows={rows} total={secs:.6}s",
                    n_buffers
                );
            }
        }
    }
}
