//! Full plugin-cache build across chromosomes. Owns the per-chrom loop
//! (mirrors `cache_builder::CacheBuilder`): calls `build_plugin_chrom` per
//! chrom against `<variation_cache_dir>/variation/<chrom>.parquet`, accumulates
//! the per-chrom entries into one `CacheManifest`, and writes
//! `plugin/<name>/manifest.json`.

use std::path::{Path, PathBuf};

use datafusion::common::{DataFusionError, Result};
use log::info;

use crate::cache::manifest::canonical_chrom_label;
use crate::plugin_cache::build::build_plugin_chrom;
use crate::plugin_cache::cache_manifest::{CacheManifest, ChromEntry, SourceRecord};
use crate::plugin_cache::source_manifest::SourceManifest;
use crate::plugin_cache::source_verify::{SourceVerification, verify_sources};

/// Builds every requested chromosome's plugin shard against a variation cache.
pub struct PluginCacheBuilder<'a> {
    manifest: &'a SourceManifest,
    manifest_file: String,
    variation_cache_dir: PathBuf,
    out: PathBuf,
    chrom_filter: Option<Vec<String>>,
    overwrite: bool,
    verification: SourceVerification,
}

impl<'a> PluginCacheBuilder<'a> {
    pub fn new(
        manifest: &'a SourceManifest,
        manifest_file: impl Into<String>,
        variation_cache_dir: impl Into<PathBuf>,
        out: impl Into<PathBuf>,
    ) -> Self {
        Self {
            manifest,
            manifest_file: manifest_file.into(),
            variation_cache_dir: variation_cache_dir.into(),
            out: out.into(),
            chrom_filter: None,
            overwrite: false,
            verification: SourceVerification::default(),
        }
    }

    /// Restrict the build to the given chromosomes (empty = no filter → all).
    pub fn with_chrom_filter<I, S>(mut self, chroms: I) -> Self
    where
        I: IntoIterator<Item = S>,
        S: Into<String>,
    {
        let chroms: Vec<String> = chroms.into_iter().map(Into::into).collect();
        self.chrom_filter = (!chroms.is_empty()).then_some(chroms);
        self
    }

    pub fn with_overwrite(mut self, overwrite: bool) -> Self {
        self.overwrite = overwrite;
        self
    }

    /// How to treat each source's declared `md5`/`path_md5` before building.
    /// Defaults to [`SourceVerification::Strict`]; a manifest that declares no
    /// digest is never hashed regardless of mode.
    pub fn with_source_verification(mut self, verification: SourceVerification) -> Self {
        self.verification = verification;
        self
    }

    /// Per-chrom shard path: `<variation_cache_dir>/variation/<canonical>.parquet`.
    fn variation_shard(&self, chrom: &str) -> PathBuf {
        self.variation_cache_dir
            .join("variation")
            .join(format!("{}.parquet", canonical_chrom_label(chrom)))
    }

    /// Chroms to build: the explicit filter, else every `variation/chr*.parquet`.
    fn resolve_chroms(&self) -> Result<Vec<String>> {
        if let Some(f) = &self.chrom_filter {
            return Ok(f.clone());
        }
        let var_dir = self.variation_cache_dir.join("variation");
        let mut chroms = Vec::new();
        for e in std::fs::read_dir(&var_dir)
            .map_err(|e| DataFusionError::Execution(format!("read {}: {e}", var_dir.display())))?
        {
            let path = e
                .map_err(|e| DataFusionError::Execution(format!("dir entry: {e}")))?
                .path();
            if path.extension().and_then(|s| s.to_str()) == Some("parquet")
                && let Some(stem) = path.file_stem().and_then(|s| s.to_str())
            {
                chroms.push(stem.to_string());
            }
        }
        chroms.sort();
        Ok(chroms)
    }

    /// Build every resolved chromosome and write `plugin/<name>/manifest.json`.
    pub async fn build_all(&self) -> Result<CacheManifest> {
        // A manifest built programmatically may not have passed through
        // `SourceManifest::load`; the digest shape check has to run before
        // any digest is compared.
        self.manifest.validate()?;
        let plugin_dir = self.out.join("plugin").join(&self.manifest.plugin_name);
        let mut cache = CacheManifest::from_source(self.manifest, &self.manifest_file);
        let prior: Option<CacheManifest> =
            std::fs::read_to_string(plugin_dir.join("manifest.json"))
                .ok()
                .and_then(|t| serde_json::from_str::<CacheManifest>(&t).ok());
        // Verify the inputs before ingesting anything — once per source part,
        // not once per chromosome. An earlier build's records let an
        // incremental per-chrom build skip re-hashing an unchanged file.
        let prior_sources = prior.as_ref().map(|m| m.sources.as_slice()).unwrap_or(&[]);
        cache.sources = verify_sources(self.manifest, self.verification, prior_sources)?;
        // Preserve chromosomes from a prior build (their shards remain on disk),
        // so a filtered/incremental build UPSERTs the rebuilt chroms rather than
        // dropping the others from the manifest that runtime discovery relies on.
        // BUT only when the prior build is compatible with this one:
        // - same output schema — a filtered rebuild after a value/match-column
        //   change must NOT keep old shards under the new schema (that would
        //   misproject them / panic in `from_batch`);
        // - same declared source digests — a filtered rebuild against a new
        //   source release must NOT keep shards built from the old one under
        //   a manifest that claims the new digest.
        // In either case the stale entries are dropped.
        let mut chroms: Vec<ChromEntry> = prior
            .as_ref()
            .filter(|old| schema_matches(old, &cache))
            .filter(|old| {
                let same = provenance_matches(old, &cache.sources);
                if !same {
                    info!(
                        "plugin '{}': source digests differ from the earlier build; its \
                         other chromosomes are dropped from the manifest",
                        self.manifest.plugin_name
                    );
                }
                same
            })
            .map(|m| m.chroms.clone())
            .unwrap_or_default();
        let _ = self.overwrite; // build_plugin_chrom overwrites the shard file per chrom.
        for chrom in self.resolve_chroms()? {
            let shard = self.variation_shard(&chrom);
            if !shard.exists() {
                return Err(DataFusionError::Execution(format!(
                    "variation shard not found: {}",
                    shard.display()
                )));
            }
            let entry = build_plugin_chrom(
                self.manifest,
                &self.manifest_file,
                &shard,
                &self.out,
                &chrom,
            )
            .await?;
            chroms.retain(|c| c.chrom != entry.chrom);
            chroms.push(entry);
        }
        chroms.sort_by(|a, b| a.chrom.cmp(&b.chrom));
        cache.chroms = chroms;
        std::fs::create_dir_all(&plugin_dir).map_err(|e| {
            DataFusionError::Execution(format!("mkdir {}: {e}", plugin_dir.display()))
        })?;
        cache.write(&plugin_dir)?;
        Ok(cache)
    }
}

/// True when two cache manifests declare the same plugin output schema — same
/// value columns (name/csq_field/type) and match discriminators (column/template)
/// in the same order. Used to decide whether prior chrom shards can be preserved
/// across a filtered rebuild.
fn schema_matches(a: &CacheManifest, b: &CacheManifest) -> bool {
    a.value_columns.len() == b.value_columns.len()
        && a.value_columns
            .iter()
            .zip(&b.value_columns)
            .all(|(x, y)| x.column == y.column && x.csq_field == y.csq_field && x.ty == y.ty)
        && a.match_columns.len() == b.match_columns.len()
        && a.match_columns
            .iter()
            .zip(&b.match_columns)
            .all(|(x, y)| x.column == y.column && x.template == y.template)
}

/// True when the prior build's shards came from the same input as this
/// build's `sources`, part for part. Digests actually verified are compared
/// when both builds have them — two `Warn` builds can share a declared digest
/// and still hash different bytes — otherwise the declared digests are. A
/// prior manifest with no `sources` block (built before provenance was
/// recorded) cannot be compared and is trusted, as it was before. A part with
/// no digest on either side has nothing to compare and passes.
fn provenance_matches(prior: &CacheManifest, sources: &[SourceRecord]) -> bool {
    if prior.sources.is_empty() {
        return true;
    }
    prior.sources.len() == sources.len()
        && sources.iter().all(|new| {
            prior
                .sources
                .iter()
                .any(|old| old.part == new.part && same_source_digest(old, new))
        })
}

fn same_source_digest(old: &SourceRecord, new: &SourceRecord) -> bool {
    match (old.verified_md5.as_deref(), new.verified_md5.as_deref()) {
        (Some(a), Some(b)) => a == b,
        _ => match (old.expected_md5(), new.expected_md5()) {
            (Some(a), Some(b)) => a == b,
            _ => true,
        },
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use datafusion::arrow::array::{Int8Array, StringArray, UInt32Array};
    use datafusion::arrow::datatypes::{DataType, Field, Schema};
    use datafusion::arrow::record_batch::RecordBatch;
    use noodles_core_tabix::Position;
    use noodles_csi::{self as csi, binning_index::index::reference_sequence::bin::Chunk};
    use parquet::arrow::ArrowWriter;
    use std::io::Write;
    use std::sync::Arc;

    fn write_gz(path: &Path, body: &str) {
        let f = std::fs::File::create(path).unwrap();
        let mut enc = flate2::write::GzEncoder::new(f, flate2::Compression::default());
        enc.write_all(body.as_bytes()).unwrap();
        enc.finish().unwrap();
    }

    fn write_bgzf_tabix_bed(path: &Path, rows: &[(&str, usize, usize, &str, &str, &str)]) {
        let file = std::fs::File::create(path).unwrap();
        let mut writer = noodles_bgzf::io::Writer::new(file);
        let mut indexer = noodles_tabix::index::Indexer::default();
        indexer.set_header(csi::binning_index::index::header::Builder::bed().build());
        let mut chunk_start = writer.virtual_position();
        for &(chrom, start, end, reference, alternate, score) in rows {
            writeln!(
                writer,
                "{chrom}\t{start}\t{end}\t{reference}\t{alternate}\t{score}"
            )
            .unwrap();
            let chunk_end = writer.virtual_position();
            indexer
                .add_record(
                    chrom,
                    Position::try_from(start + 1).unwrap(),
                    Position::try_from(end).unwrap(),
                    Chunk::new(chunk_start, chunk_end),
                )
                .unwrap();
            chunk_start = chunk_end;
        }
        writer.finish().unwrap();
        let index_file = std::fs::File::create(format!("{}.tbi", path.display())).unwrap();
        noodles_tabix::io::Writer::new(index_file)
            .write_index(&indexer.build())
            .unwrap();
    }

    fn write_variation(path: &Path, rows: &[(&str, u32, &str, i8)]) {
        let schema = Arc::new(Schema::new(vec![
            Field::new("chrom", DataType::Utf8, false),
            Field::new("start", DataType::UInt32, false),
            Field::new("allele_string", DataType::Utf8, false),
            Field::new("tier", DataType::Int8, false),
        ]));
        let batch = RecordBatch::try_new(
            schema.clone(),
            vec![
                Arc::new(StringArray::from(
                    rows.iter().map(|r| r.0).collect::<Vec<_>>(),
                )),
                Arc::new(UInt32Array::from(
                    rows.iter().map(|r| r.1).collect::<Vec<_>>(),
                )),
                Arc::new(StringArray::from(
                    rows.iter().map(|r| r.2).collect::<Vec<_>>(),
                )),
                Arc::new(Int8Array::from(
                    rows.iter().map(|r| r.3).collect::<Vec<_>>(),
                )),
            ],
        )
        .unwrap();
        let file = std::fs::File::create(path).unwrap();
        let mut w = ArrowWriter::try_new(file, schema, None).unwrap();
        w.write(&batch).unwrap();
        w.close().unwrap();
    }

    #[tokio::test(flavor = "multi_thread")]
    async fn builds_all_chroms_and_writes_manifest() {
        let dir = tempfile::tempdir().unwrap();
        let src = dir.path().join("src.tsv.gz");
        write_bgzf_tabix_bed(
            &src,
            &[
                ("chr1", 99, 100, "A", "G", "0.9"),
                ("chr2", 199, 200, "C", "T", "0.5"),
            ],
        );

        let var_dir = dir.path().join("cache").join("variation");
        std::fs::create_dir_all(&var_dir).unwrap();
        write_variation(&var_dir.join("chr1.parquet"), &[("1", 100, "A/G", 0)]);
        write_variation(&var_dir.join("chr2.parquet"), &[("2", 200, "C/T", 1)]);

        let toml = format!(
            r##"
plugin_name = "demo"
coordinate_system = "0-based-half-open"
ingest_sql = """
SELECT chrom, CAST(start0 AS INT) AS start, CAST(end0 AS INT) AS end,
       concat(ref, '/', alt) AS allele_string, CAST(score AS FLOAT) AS demo_score
FROM plugin_demo_src
"""

[[source]]
provider = "csv"
path = "{}"
index = "tabix"
  [source.csv]
  delimiter = "\t"
  has_header = false
  compression = "gzip"
  schema = [
    {{ name = "chrom", type = "Utf8" }},
    {{ name = "start0", type = "Utf8" }},
    {{ name = "end0",   type = "Utf8" }},
    {{ name = "ref",   type = "Utf8" }},
    {{ name = "alt",   type = "Utf8" }},
    {{ name = "score", type = "Utf8" }},
  ]

[[value_columns]]
column = "demo_score"
csq_field = "DEMO"
type = "Float32"
"##,
            src.display()
        );
        let manifest: SourceManifest = toml::from_str(&toml).unwrap();
        let out = dir.path().join("out");
        let cache = PluginCacheBuilder::new(
            &manifest,
            "demo.source.toml",
            dir.path().join("cache"),
            &out,
        )
        .with_chrom_filter(["1", "2"])
        .build_all()
        .await
        .unwrap();
        assert_eq!(cache.chroms.len(), 2);
        assert!(
            out.join("plugin")
                .join("demo")
                .join("manifest.json")
                .exists()
        );
        assert!(
            out.join("plugin")
                .join("demo")
                .join("chr1.parquet")
                .exists()
        );
        assert!(
            out.join("plugin")
                .join("demo")
                .join("chr2.parquet")
                .exists()
        );
    }

    /// A tabix-indexed BED-like source plus a manifest that declares a digest
    /// for it; `md5` is substituted into the manifest text verbatim.
    fn verified_fixture(dir: &Path, md5: &str) -> (SourceManifest, PathBuf, PathBuf) {
        let src = dir.join("src.tsv.gz");
        write_bgzf_tabix_bed(&src, &[("chr1", 99, 100, "A", "G", "0.9")]);
        let var_dir = dir.join("cache").join("variation");
        std::fs::create_dir_all(&var_dir).unwrap();
        write_variation(&var_dir.join("chr1.parquet"), &[("1", 100, "A/G", 0)]);
        write_variation(&var_dir.join("chr2.parquet"), &[("2", 200, "C/T", 1)]);
        let toml = format!(
            r##"
plugin_name = "demo"
coordinate_system = "0-based-half-open"
ingest_sql = """
SELECT chrom, CAST(start0 AS INT) AS start, CAST(end0 AS INT) AS end,
       concat(ref, '/', alt) AS allele_string, CAST(score AS FLOAT) AS demo_score
FROM plugin_demo_src_scores
"""

[[source]]
part = "scores"
provider = "csv"
path = "{}"
url = "https://example.org/demo/scores.tsv.gz"
md5 = "{md5}"
index = "tabix"
  [source.csv]
  delimiter = "	"
  has_header = false
  compression = "gzip"
  schema = [
    {{ name = "chrom", type = "Utf8" }},
    {{ name = "start0", type = "Utf8" }},
    {{ name = "end0",   type = "Utf8" }},
    {{ name = "ref",   type = "Utf8" }},
    {{ name = "alt",   type = "Utf8" }},
    {{ name = "score", type = "Utf8" }},
  ]

[[value_columns]]
column = "demo_score"
csq_field = "DEMO"
type = "Float32"
"##,
            src.display()
        );
        let manifest: SourceManifest = toml::from_str(&toml).unwrap();
        (manifest, dir.join("cache"), dir.join("out"))
    }

    const WRONG_MD5: &str = "00000000000000000000000000000000";

    #[tokio::test(flavor = "multi_thread")]
    async fn strict_verification_rejects_a_source_whose_md5_differs() {
        let dir = tempfile::tempdir().unwrap();
        let (manifest, cache_dir, out) = verified_fixture(dir.path(), WRONG_MD5);
        let (actual, _) =
            crate::plugin_cache::source_verify::md5_file(Path::new(&manifest.sources[0].path))
                .unwrap();
        let error = PluginCacheBuilder::new(&manifest, "demo.source.toml", &cache_dir, &out)
            .with_chrom_filter(["1"])
            .build_all()
            .await
            .unwrap_err()
            .to_string();
        assert!(error.contains("MD5 mismatch"), "{error}");
        assert!(error.contains("part 'scores'"), "{error}");
        assert!(error.contains(WRONG_MD5), "{error}");
        assert!(error.contains(&actual), "{error}");
        assert!(
            error.contains("https://example.org/demo/scores.tsv.gz"),
            "{error}"
        );
        // Nothing was ingested or written: the check runs before the first chrom.
        assert!(!out.join("plugin").join("demo").exists());
    }

    #[tokio::test(flavor = "multi_thread")]
    async fn matching_md5_is_verified_and_recorded_in_the_manifest() {
        let dir = tempfile::tempdir().unwrap();
        let (probe, _, _) = verified_fixture(dir.path(), WRONG_MD5);
        let src = PathBuf::from(&probe.sources[0].path);
        let (actual, size) = crate::plugin_cache::source_verify::md5_file(&src).unwrap();
        let (manifest, cache_dir, out) = verified_fixture(dir.path(), &actual);
        let cache = PluginCacheBuilder::new(&manifest, "demo.source.toml", &cache_dir, &out)
            .with_chrom_filter(["1"])
            .build_all()
            .await
            .unwrap();
        assert_eq!(cache.sources.len(), 1);
        let record = &cache.sources[0];
        assert_eq!(record.part.as_deref(), Some("scores"));
        assert_eq!(record.file, "src.tsv.gz");
        assert_eq!(
            record.url.as_deref(),
            Some("https://example.org/demo/scores.tsv.gz")
        );
        assert_eq!(record.md5.as_deref(), Some(actual.as_str()));
        assert_eq!(record.path_md5, None);
        assert_eq!(record.verified_md5.as_deref(), Some(actual.as_str()));
        assert_eq!(record.size, Some(size));
        assert!(record.mtime_ns.is_some());

        // The record survives the round trip through manifest.json.
        let found = crate::plugin_cache::cache_manifest::discover_plugins(&out).unwrap();
        assert_eq!(found[0].sources, cache.sources);
    }

    #[tokio::test(flavor = "multi_thread")]
    async fn warn_mode_builds_and_records_the_digest_it_found() {
        let dir = tempfile::tempdir().unwrap();
        let (manifest, cache_dir, out) = verified_fixture(dir.path(), WRONG_MD5);
        let cache = PluginCacheBuilder::new(&manifest, "demo.source.toml", &cache_dir, &out)
            .with_chrom_filter(["1"])
            .with_source_verification(SourceVerification::Warn)
            .build_all()
            .await
            .unwrap();
        let record = &cache.sources[0];
        assert_eq!(record.md5.as_deref(), Some(WRONG_MD5));
        let verified = record.verified_md5.as_deref().unwrap();
        assert_ne!(verified, WRONG_MD5);
        assert_eq!(verified.len(), 32);
        assert_eq!(cache.chroms.len(), 1);
    }

    #[tokio::test(flavor = "multi_thread")]
    async fn skip_mode_keeps_provenance_but_records_no_digest() {
        let dir = tempfile::tempdir().unwrap();
        let (manifest, cache_dir, out) = verified_fixture(dir.path(), WRONG_MD5);
        let cache = PluginCacheBuilder::new(&manifest, "demo.source.toml", &cache_dir, &out)
            .with_chrom_filter(["1"])
            .with_source_verification(SourceVerification::Skip)
            .build_all()
            .await
            .unwrap();
        let record = &cache.sources[0];
        assert_eq!(
            record.url.as_deref(),
            Some("https://example.org/demo/scores.tsv.gz")
        );
        assert_eq!(record.md5.as_deref(), Some(WRONG_MD5));
        assert_eq!(record.verified_md5, None);
        assert_eq!(record.size, None);
    }

    #[tokio::test(flavor = "multi_thread")]
    async fn path_md5_is_checked_instead_of_the_upstream_md5() {
        let dir = tempfile::tempdir().unwrap();
        let (mut manifest, cache_dir, out) = verified_fixture(dir.path(), WRONG_MD5);
        let (actual, _) =
            crate::plugin_cache::source_verify::md5_file(Path::new(&manifest.sources[0].path))
                .unwrap();
        manifest.sources[0].path_md5 = Some(actual.clone());
        let cache = PluginCacheBuilder::new(&manifest, "demo.source.toml", &cache_dir, &out)
            .with_chrom_filter(["1"])
            .build_all()
            .await
            .unwrap();
        let record = &cache.sources[0];
        assert_eq!(record.md5.as_deref(), Some(WRONG_MD5));
        assert_eq!(record.path_md5.as_deref(), Some(actual.as_str()));
        assert_eq!(record.verified_md5.as_deref(), Some(actual.as_str()));
    }

    #[tokio::test(flavor = "multi_thread")]
    async fn a_manifest_without_a_digest_is_not_hashed() {
        let dir = tempfile::tempdir().unwrap();
        let (mut manifest, cache_dir, out) = verified_fixture(dir.path(), WRONG_MD5);
        manifest.sources[0].md5 = None;
        let cache = PluginCacheBuilder::new(&manifest, "demo.source.toml", &cache_dir, &out)
            .with_chrom_filter(["1"])
            .build_all()
            .await
            .unwrap();
        let record = &cache.sources[0];
        assert_eq!(record.md5, None);
        assert_eq!(record.verified_md5, None);
        assert!(record.url.is_some());
    }

    // A per-chromosome workflow calls the builder once per contig against the
    // same 87 GB input; it must hash that input once, not once per call.
    #[tokio::test(flavor = "multi_thread")]
    async fn incremental_build_trusts_an_earlier_verification_of_an_unchanged_file() {
        use crate::plugin_cache::cache_manifest::SourceRecord;

        let dir = tempfile::tempdir().unwrap();
        let (manifest, cache_dir, out) = verified_fixture(dir.path(), WRONG_MD5);
        let src = Path::new(&manifest.sources[0].path);
        let meta = std::fs::metadata(src).unwrap();
        let mtime_ns = meta
            .modified()
            .unwrap()
            .duration_since(std::time::UNIX_EPOCH)
            .unwrap()
            .as_nanos() as i64;

        // An earlier build's manifest claims to have verified the (wrong)
        // digest over this exact file name, size and mtime. A strict build now
        // succeeds only if it trusts that record instead of re-hashing.
        let plugin_dir = out.join("plugin").join("demo");
        std::fs::create_dir_all(&plugin_dir).unwrap();
        let mut earlier = CacheManifest::from_source(&manifest, "demo.source.toml");
        let trusted = SourceRecord {
            verified_md5: Some(WRONG_MD5.into()),
            size: Some(meta.len()),
            mtime_ns: Some(mtime_ns),
            ..SourceRecord::from_spec(&manifest.sources[0])
        };

        // The record is bound to the file name: a record for another file
        // with the same size and mtime is not trusted, and the re-hash fails.
        earlier.sources = vec![SourceRecord {
            file: "other.tsv.gz".into(),
            ..trusted.clone()
        }];
        earlier.write(&plugin_dir).unwrap();
        let error = PluginCacheBuilder::new(&manifest, "demo.source.toml", &cache_dir, &out)
            .with_chrom_filter(["1"])
            .build_all()
            .await
            .unwrap_err()
            .to_string();
        assert!(error.contains("MD5 mismatch"), "{error}");

        earlier.sources = vec![trusted];
        earlier.write(&plugin_dir).unwrap();
        let cache = PluginCacheBuilder::new(&manifest, "demo.source.toml", &cache_dir, &out)
            .with_chrom_filter(["1"])
            .build_all()
            .await
            .unwrap();
        assert_eq!(cache.sources[0].verified_md5.as_deref(), Some(WRONG_MD5));

        // Touching the file invalidates the fingerprint: the next strict build
        // re-hashes and reports the mismatch.
        let file = std::fs::File::options().write(true).open(src).unwrap();
        file.set_modified(
            std::time::SystemTime::UNIX_EPOCH + std::time::Duration::from_secs(1_700_000_000),
        )
        .unwrap();
        drop(file);
        let error = PluginCacheBuilder::new(&manifest, "demo.source.toml", &cache_dir, &out)
            .with_chrom_filter(["2"])
            .build_all()
            .await
            .unwrap_err()
            .to_string();
        assert!(error.contains("MD5 mismatch"), "{error}");
    }

    // A filtered rebuild from the SAME verified input keeps the other
    // chromosomes — and the second call trusts the first call's verification.
    #[tokio::test(flavor = "multi_thread")]
    async fn filtered_rebuild_from_the_same_source_preserves_other_chroms() {
        let dir = tempfile::tempdir().unwrap();
        let (probe, _, _) = verified_fixture(dir.path(), WRONG_MD5);
        let (actual, _) =
            crate::plugin_cache::source_verify::md5_file(Path::new(&probe.sources[0].path))
                .unwrap();
        let (manifest, cache_dir, out) = verified_fixture(dir.path(), &actual);
        PluginCacheBuilder::new(&manifest, "demo.source.toml", &cache_dir, &out)
            .with_chrom_filter(["1"])
            .build_all()
            .await
            .unwrap();
        let cache = PluginCacheBuilder::new(&manifest, "demo.source.toml", &cache_dir, &out)
            .with_chrom_filter(["2"])
            .build_all()
            .await
            .unwrap();
        let chroms: Vec<&str> = cache.chroms.iter().map(|c| c.chrom.as_str()).collect();
        assert_eq!(chroms, ["chr1", "chr2"]);
        assert_eq!(
            cache.sources[0].verified_md5.as_deref(),
            Some(actual.as_str())
        );
    }

    // A filtered rebuild against a NEW source release must not keep shards
    // built from the old one under a manifest that claims the new digest.
    #[tokio::test(flavor = "multi_thread")]
    async fn filtered_rebuild_after_a_source_change_drops_stale_chroms() {
        let dir = tempfile::tempdir().unwrap();
        let (probe, _, _) = verified_fixture(dir.path(), WRONG_MD5);
        let (actual, _) =
            crate::plugin_cache::source_verify::md5_file(Path::new(&probe.sources[0].path))
                .unwrap();
        let (manifest, cache_dir, out) = verified_fixture(dir.path(), &actual);
        PluginCacheBuilder::new(&manifest, "demo.source.toml", &cache_dir, &out)
            .with_chrom_filter(["1"])
            .build_all()
            .await
            .unwrap();

        // A new release lands at the same path and the manifest declares its
        // digest; build chr2 from it, strictly.
        let src = PathBuf::from(&manifest.sources[0].path);
        write_bgzf_tabix_bed(
            &src,
            &[
                ("chr1", 99, 100, "A", "G", "0.1"),
                ("chr2", 199, 200, "C", "T", "0.5"),
            ],
        );
        let (release2, _) = crate::plugin_cache::source_verify::md5_file(&src).unwrap();
        assert_ne!(release2, actual);
        let mut next = manifest.clone();
        next.sources[0].md5 = Some(release2.clone());
        let cache = PluginCacheBuilder::new(&next, "demo.source.toml", &cache_dir, &out)
            .with_chrom_filter(["2"])
            .build_all()
            .await
            .unwrap();
        let chroms: Vec<&str> = cache.chroms.iter().map(|c| c.chrom.as_str()).collect();
        assert_eq!(chroms, ["chr2"], "chr1 was built from another release");
        assert_eq!(
            cache.sources[0].verified_md5.as_deref(),
            Some(release2.as_str())
        );
    }

    // A cache built before provenance was recorded has no `sources` block;
    // an incremental build into it keeps its chromosomes as it always did.
    #[tokio::test(flavor = "multi_thread")]
    async fn filtered_rebuild_into_a_cache_without_provenance_preserves_chroms() {
        let dir = tempfile::tempdir().unwrap();
        let (probe, _, _) = verified_fixture(dir.path(), WRONG_MD5);
        let (actual, _) =
            crate::plugin_cache::source_verify::md5_file(Path::new(&probe.sources[0].path))
                .unwrap();
        let (manifest, cache_dir, out) = verified_fixture(dir.path(), &actual);
        let plugin_dir = out.join("plugin").join("demo");
        std::fs::create_dir_all(&plugin_dir).unwrap();
        let mut earlier = CacheManifest::from_source(&manifest, "demo.source.toml");
        earlier.chroms = vec![ChromEntry {
            chrom: "chr1".into(),
            file: "chr1.parquet".into(),
            rows: 1,
            warm: 1,
            cold: 0,
        }];
        earlier.write(&plugin_dir).unwrap();
        assert!(
            !std::fs::read_to_string(plugin_dir.join("manifest.json"))
                .unwrap()
                .contains("\"sources\"")
        );

        let cache = PluginCacheBuilder::new(&manifest, "demo.source.toml", &cache_dir, &out)
            .with_chrom_filter(["2"])
            .build_all()
            .await
            .unwrap();
        let chroms: Vec<&str> = cache.chroms.iter().map(|c| c.chrom.as_str()).collect();
        assert_eq!(chroms, ["chr1", "chr2"]);
    }

    #[test]
    fn provenance_comparison_is_per_part_and_tolerates_undeclared_digests() {
        let rec = |part: Option<&str>, md5: Option<&str>| SourceRecord {
            part: part.map(String::from),
            file: "x".into(),
            url: None,
            md5: md5.map(String::from),
            path_md5: None,
            verified_md5: None,
            size: None,
            mtime_ns: None,
        };
        let dir = tempfile::tempdir().unwrap();
        let (manifest, _, _) = verified_fixture(dir.path(), WRONG_MD5);
        let mut prior = CacheManifest::from_source(&manifest, "demo.source.toml");
        prior.sources = vec![rec(Some("snv"), Some("a")), rec(Some("indel"), Some("b"))];
        let same = [rec(Some("snv"), Some("a")), rec(Some("indel"), Some("b"))];
        let changed = [rec(Some("snv"), Some("a")), rec(Some("indel"), Some("c"))];
        let undeclared = [rec(Some("snv"), None), rec(Some("indel"), Some("b"))];
        let fewer = [rec(Some("snv"), Some("a"))];
        assert!(provenance_matches(&prior, &same));
        assert!(!provenance_matches(&prior, &changed));
        assert!(provenance_matches(&prior, &undeclared));
        assert!(!provenance_matches(&prior, &fewer));

        // Verified digests win over declared ones when both builds have them:
        // two warn-mode builds can share `md5 = "a"` and still hash different
        // bytes; a strict build after a warn-mode one likewise differs.
        let verified = |part: &str, md5: &str, actual: &str| SourceRecord {
            verified_md5: Some(actual.into()),
            ..rec(Some(part), Some(md5))
        };
        prior.sources = vec![verified("snv", "a", "a"), verified("indel", "b", "x")];
        assert!(provenance_matches(
            &prior,
            &[verified("snv", "a", "a"), verified("indel", "b", "x")]
        ));
        assert!(!provenance_matches(
            &prior,
            &[verified("snv", "a", "a"), verified("indel", "b", "y")]
        ));
        assert!(!provenance_matches(
            &prior,
            &[verified("snv", "a", "a"), verified("indel", "b", "b")]
        ));
        // A skipped verification on one side falls back to the declared digest.
        assert!(provenance_matches(&prior, &same));

        prior.sources.clear();
        assert!(provenance_matches(&prior, &changed));
    }

    // Two warn-mode builds against the same declared digest but different
    // bytes must not be merged into one cache.
    #[tokio::test(flavor = "multi_thread")]
    async fn warn_mode_rebuild_from_different_bytes_drops_stale_chroms() {
        let dir = tempfile::tempdir().unwrap();
        let (manifest, cache_dir, out) = verified_fixture(dir.path(), WRONG_MD5);
        let src = PathBuf::from(&manifest.sources[0].path);
        let first = PluginCacheBuilder::new(&manifest, "demo.source.toml", &cache_dir, &out)
            .with_chrom_filter(["1"])
            .with_source_verification(SourceVerification::Warn)
            .build_all()
            .await
            .unwrap();

        // Same path, same declared digest, different bytes.
        write_bgzf_tabix_bed(
            &src,
            &[
                ("chr1", 99, 100, "A", "G", "0.1"),
                ("chr2", 199, 200, "C", "T", "0.5"),
            ],
        );
        let second = PluginCacheBuilder::new(&manifest, "demo.source.toml", &cache_dir, &out)
            .with_chrom_filter(["2"])
            .with_source_verification(SourceVerification::Warn)
            .build_all()
            .await
            .unwrap();
        assert_ne!(
            first.sources[0].verified_md5,
            second.sources[0].verified_md5
        );
        let chroms: Vec<&str> = second.chroms.iter().map(|c| c.chrom.as_str()).collect();
        assert_eq!(chroms, ["chr2"], "chr1 was built from other bytes");
    }

    // A filtered/incremental rebuild must UPSERT into the existing manifest, not
    // drop previously built chromosomes (PR #190 N1).
    #[tokio::test(flavor = "multi_thread")]
    async fn filtered_rebuild_preserves_other_chroms() {
        let dir = tempfile::tempdir().unwrap();
        let src = dir.path().join("src.tsv.gz");
        write_gz(&src, "chr1\t100\tA\tG\t0.9\nchr2\t200\tC\tT\t0.5\n");
        let var_dir = dir.path().join("cache").join("variation");
        std::fs::create_dir_all(&var_dir).unwrap();
        write_variation(&var_dir.join("chr1.parquet"), &[("1", 100, "A/G", 0)]);
        write_variation(&var_dir.join("chr2.parquet"), &[("2", 200, "C/T", 1)]);
        let toml = format!(
            r##"
plugin_name = "demo"
coordinate_system = "1-based"
ingest_sql = """
SELECT chrom, CAST(pos AS INT) AS start, CAST(pos AS INT) AS end,
       concat(ref, '/', alt) AS allele_string, CAST(score AS FLOAT) AS demo_score
FROM plugin_demo_src
"""

[[source]]
provider = "csv"
path = "{}"
  [source.csv]
  delimiter = "\t"
  has_header = false
  compression = "gzip"
  schema = [
    {{ name = "chrom", type = "Utf8" }},
    {{ name = "pos",   type = "Utf8" }},
    {{ name = "ref",   type = "Utf8" }},
    {{ name = "alt",   type = "Utf8" }},
    {{ name = "score", type = "Utf8" }},
  ]

[[value_columns]]
column = "demo_score"
csq_field = "DEMO"
type = "Float32"
"##,
            src.display()
        );
        let manifest: SourceManifest = toml::from_str(&toml).unwrap();
        let cache_dir = dir.path().join("cache");
        let out = dir.path().join("out");

        // Build chr1 only, then rebuild chr2 only into the same output.
        PluginCacheBuilder::new(&manifest, "demo.source.toml", &cache_dir, &out)
            .with_chrom_filter(["1"])
            .build_all()
            .await
            .unwrap();
        let cache = PluginCacheBuilder::new(&manifest, "demo.source.toml", &cache_dir, &out)
            .with_chrom_filter(["2"])
            .build_all()
            .await
            .unwrap();

        let chroms: Vec<&str> = cache.chroms.iter().map(|c| c.chrom.as_str()).collect();
        assert!(
            chroms.contains(&"chr1"),
            "chr1 must be preserved: {chroms:?}"
        );
        assert!(chroms.contains(&"chr2"), "chr2 must be present: {chroms:?}");
        assert_eq!(cache.chroms.len(), 2);
    }

    // A filtered rebuild after a value/match-column change must NOT preserve old
    // chroms (their shards use the old schema) — PR #190 R2.
    #[test]
    fn schema_change_drops_stale_chroms() {
        use crate::plugin_cache::cache_manifest::ValueColumnRecord;
        let mk = |csq: &str| CacheManifest {
            plugin_name: "demo".into(),
            source_manifest: "demo.source.toml".into(),
            key_columns: vec![],
            match_columns: vec![],
            value_columns: vec![ValueColumnRecord {
                column: "score".into(),
                csq_field: csq.into(),
                ty: "Float32".into(),
                description: None,
            }],
            chroms: vec![],
            sources: vec![],
            cache_source_version: None,
            allele_match: Default::default(),
            csq_rank: 0,
            field_order: Default::default(),
        };
        assert!(schema_matches(&mk("DEMO"), &mk("DEMO")));
        assert!(!schema_matches(&mk("DEMO"), &mk("DEMO2")));
    }
}
