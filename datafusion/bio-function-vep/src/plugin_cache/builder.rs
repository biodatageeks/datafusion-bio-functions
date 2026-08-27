//! Full plugin-cache build across chromosomes. Owns the per-chrom loop
//! (mirrors `cache_builder::CacheBuilder`): calls `build_plugin_chrom` per
//! chrom against `<variation_cache_dir>/variation/<chrom>.parquet`, accumulates
//! the per-chrom entries into one `CacheManifest`, and writes
//! `plugin/<name>/manifest.json`.

use std::path::{Path, PathBuf};

use datafusion::common::{DataFusionError, Result};

use crate::cache::manifest::canonical_chrom_label;
use crate::plugin_cache::build::build_plugin_chrom;
use crate::plugin_cache::cache_manifest::{CacheManifest, ChromEntry};
use crate::plugin_cache::source_manifest::SourceManifest;

/// Builds every requested chromosome's plugin shard against a variation cache.
pub struct PluginCacheBuilder<'a> {
    manifest: &'a SourceManifest,
    manifest_file: String,
    variation_cache_dir: PathBuf,
    out: PathBuf,
    chrom_filter: Option<Vec<String>>,
    overwrite: bool,
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
        let plugin_dir = self.out.join("plugin").join(&self.manifest.plugin_name);
        let mut cache = CacheManifest::from_source(self.manifest, &self.manifest_file);
        // Preserve chromosomes from a prior build (their shards remain on disk),
        // so a filtered/incremental build UPSERTs the rebuilt chroms rather than
        // dropping the others from the manifest that runtime discovery relies on.
        // BUT only when the prior manifest's schema matches the new one — a
        // filtered rebuild after a value/match-column change must NOT keep old
        // shards under the new schema (that would misproject them / panic in
        // `from_batch`); in that case we drop the stale entries.
        let mut chroms: Vec<ChromEntry> = std::fs::read_to_string(plugin_dir.join("manifest.json"))
            .ok()
            .and_then(|t| serde_json::from_str::<CacheManifest>(&t).ok())
            .filter(|old| schema_matches(old, &cache))
            .map(|m| m.chroms)
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
            cache_source_version: None,
            allele_match: Default::default(),
            csq_rank: 0,
            field_order: Default::default(),
        };
        assert!(schema_matches(&mk("DEMO"), &mk("DEMO")));
        assert!(!schema_matches(&mk("DEMO"), &mk("DEMO2")));
    }
}
