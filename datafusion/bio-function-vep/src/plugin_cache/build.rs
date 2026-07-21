//! End-to-end per-chrom plugin cache build: register sources → ingest view →
//! normalization wrapper → variation-frequency join/tier → two-pass tiered
//! shard write → cache-manifest chrom entry.

use std::path::Path;
use std::sync::Arc;

use datafusion::arrow::array::{Array, BooleanArray, Int8Array, UInt32Array};
use datafusion::arrow::compute::{cast, filter_record_batch};
use datafusion::arrow::datatypes::SchemaRef;
use datafusion::arrow::record_batch::RecordBatch;
use datafusion::common::{DataFusionError, Result};
use datafusion::datasource::MemTable;
use datafusion::prelude::{SessionConfig, SessionContext};
use futures::{StreamExt, TryStreamExt};
use log::info;
use std::time::Instant;

use crate::cache::manifest::canonical_chrom_label;
use crate::plugin_cache::cache_manifest::ChromEntry;
use crate::plugin_cache::dedup::dedup_keep_first;
use crate::plugin_cache::join::tiered_stream;
use crate::plugin_cache::normalize::{
    canonical_contig_str, canonical_contig_udf, wrap_normalization,
};
use crate::plugin_cache::provider::register_sources;
use crate::plugin_cache::source_manifest::SourceManifest;
use crate::plugin_cache::write::{PluginShardWriter, plugin_output_schema};

/// Reproject `batch` to `schema`'s column order and types (casting where needed,
/// e.g. the normalized view's `Int64` `start`/`end` → the shard's `UInt32`).
fn reproject_cast(batch: &RecordBatch, schema: &SchemaRef) -> Result<RecordBatch> {
    let cols = schema
        .fields()
        .iter()
        .map(|f| {
            let idx = batch.schema().index_of(f.name())?;
            cast(batch.column(idx), f.data_type())
                .map_err(|e| DataFusionError::Execution(format!("cast {}: {e}", f.name())))
        })
        .collect::<Result<Vec<_>>>()?;
    RecordBatch::try_new(Arc::clone(schema), cols)
        .map_err(|e| DataFusionError::Execution(format!("reproject: {e}")))
}

/// Verify `start` is non-decreasing within `batch` and against `last_seen` (the
/// previous batch's final `start` for this tier), returning the batch's own
/// final `start` on success.
///
/// The streaming shard write below relies on the tiered join emitting rows in
/// position-ascending order per tier; that holds by construction for the
/// query shapes exercised so far (see the comment above the call site), but
/// `HashJoinExec`'s build/probe side can flip based on table-size statistics
/// for other plugin/variation size ratios, and a flipped build side does not
/// preserve row order. Rather than trust that silently, check it: a real
/// violation aborts the build immediately (no shard is written) instead of
/// producing a page directory whose binary-search lookup can silently miss
/// rows that are actually present (`PageDir::resolve_ranges` assumes a
/// genuinely sorted run).
fn assert_start_monotonic(
    batch: &RecordBatch,
    start_idx: usize,
    tier: i8,
    chrom: &str,
    last_seen: Option<u32>,
) -> Result<Option<u32>> {
    let starts = batch
        .column(start_idx)
        .as_any()
        .downcast_ref::<UInt32Array>()
        .ok_or_else(|| DataFusionError::Execution("start column not UInt32".into()))?;
    let mut prev = last_seen;
    for i in 0..starts.len() {
        let v = starts.value(i);
        if let Some(p) = prev
            && v < p
        {
            return Err(DataFusionError::Execution(format!(
                "plugin_cache[{chrom}]: tier {tier} shard write is not position-ascending \
                 (start {v} follows {p}) -- the tier join did not preserve row order for \
                 this plugin/variation size ratio; the on-disk point-lookup directory \
                 requires a sorted shard, so refusing to write a corrupt one. This needs \
                 a real re-sort of the tier's output, not just a bigger error message."
            )));
        }
        prev = Some(v);
    }
    Ok(prev)
}

/// Build one chromosome's plugin shard. Returns the cache-manifest chrom entry.
pub async fn build_plugin_chrom(
    src: &SourceManifest,
    _source_manifest_file: &str,
    variation_shard: &Path,
    output_cache_root: &Path,
    chrom: &str,
) -> Result<ChromEntry> {
    // Coarse-grained stage timing at `info` level — cheap (a handful of
    // `Instant::now()` calls) and the only thing that would have turned CADD
    // chr6's multi-hour "is it stuck?" investigation (no visibility into
    // which stage was running) into a 30-second log check.
    let t_start = Instant::now();
    // Read context: single partition ONLY for the source scan → normalize →
    // dedup leg, so the CSV scan yields rows in source-file order (a
    // multi-partition scan reads byte ranges concurrently and coalescing does not
    // restore file order). The dedup needs that order to keep VEP's first-in-file
    // record for a duplicate probe key (see `dedup::dedup_keep_first`). The
    // downstream tier join + write run in a separate, default-parallel context
    // (`build_ctx` below) since post-dedup there is one row per key and order no
    // longer matters — so only this cheap read leg is serialized, not the join
    // over the (multi-GB) variation shard.
    let read_ctx = SessionContext::new_with_config(SessionConfig::new().with_target_partitions(1));
    read_ctx.register_udf(canonical_contig_udf());
    // Held until this chrom's stream is fully materialized (dedup) below, then
    // dropped (deleting the decompressed temp) — so build_all keeps at most one temp.
    let src_temps = register_sources(&read_ctx, src).await?;

    // Ingest view (raw column mapping), then normalized view (contig + coords).
    read_ctx
        .sql(&format!(
            "CREATE OR REPLACE VIEW {} AS {}",
            src.ingest_view_name(),
            src.ingest_sql
        ))
        .await?;
    let value_cols: Vec<String> = src.value_columns.iter().map(|v| v.column.clone()).collect();
    let match_cols = src.match_column_names();
    let norm_sql = wrap_normalization(
        &src.ingest_view_name(),
        src.coordinate_system.clone(),
        &match_cols,
        &value_cols,
    );
    let norm_view = format!("plugin_{}_norm", src.plugin_name);
    // Filter on the same canonicalized contig the normalization applies to the
    // data (`canonical_contig_str` folds `M`/`chrM`/`chrMT` → `MT` and uppercases),
    // so `--chrom M`/`chrM`/`MT` all select the MT rows rather than silently
    // producing a 0-row shard.
    let source_chrom = canonical_contig_str(chrom);
    // Defense-in-depth: `chrom` is a trusted local arg, but reject anything
    // outside a safe contig charset before interpolating it into SQL.
    if source_chrom.is_empty()
        || !source_chrom
            .chars()
            .all(|c| c.is_ascii_alphanumeric() || matches!(c, '_' | '.' | '-'))
    {
        return Err(DataFusionError::Execution(format!(
            "invalid contig '{chrom}'"
        )));
    }
    read_ctx
        .sql(&format!(
            "CREATE OR REPLACE VIEW {norm_view} AS SELECT * FROM ({norm_sql}) WHERE chrom = '{source_chrom}'"
        ))
        .await?;

    // Collapse duplicate probe keys to their first source-file occurrence BEFORE
    // the tier join (which reorders rows and would destroy the file-order
    // tiebreak). Two overlapping genes can map the same genomic variant + aa-change
    // to different scores; VEP takes the first-in-file record, so we must too. This
    // reads the (single-partition, file-ordered) normalized stream and keeps the
    // first row per `(start, allele_string, <match cols>)`.
    let norm_stream = read_ctx
        .sql(&format!("SELECT * FROM {norm_view}"))
        .await?
        .execute_stream()
        .await?;
    let norm_schema = norm_stream.schema();
    // `assume_unique` sources are structurally guaranteed to never repeat a
    // probe key, so the keep-first pass (a HashSet<String> with one entry per
    // row — the dominant memory cost on the largest chromosomes) is skipped;
    // batches are collected as-is.
    let deduped = if src.assume_unique {
        norm_stream.try_collect::<Vec<_>>().await?
    } else {
        dedup_keep_first(norm_stream, &match_cols).await?
    };
    // The source scan is done — drop the decompressed temp before the join leg.
    drop(src_temps);
    let read_rows: usize = deduped.iter().map(|b| b.num_rows()).sum();
    info!(
        "plugin_cache[{}/{chrom}]: read+normalize+dedup done, {read_rows} rows, {:.1}s elapsed",
        src.plugin_name,
        t_start.elapsed().as_secs_f64()
    );

    // Build context: single partition, same reasoning as `read_ctx` above — the
    // streaming write below (see `tiered_stream` consumption) relies on the join
    // emitting rows in the same position-ascending order the dedup pass fed it,
    // so it never has to buffer more than one batch. With `target_partitions=1`
    // there is exactly one task for the whole join, so HashJoinExec has no
    // opportunity to reorder rows across partitions. The join's build side
    // (`plugin_variation_probe`, DISTINCT-projected and typically far smaller
    // than the plugin's own per-chrom row count for whole-genome plugins like
    // CADD/SpliceAI) is unaffected by this — it's still a single hash-table
    // build regardless of partition count. The deduped survivors are
    // re-registered here as a MemTable the join consumes in place of the raw
    // normalized view.
    let build_ctx = SessionContext::new_with_config(SessionConfig::new().with_target_partitions(1));
    let dedup_view = format!("plugin_{}_dedup", src.plugin_name);
    let mem = MemTable::try_new(norm_schema, vec![deduped])
        .map_err(|e| DataFusionError::Execution(format!("dedup memtable: {e}")))?;
    build_ctx.register_table(&dedup_view, Arc::new(mem))?;

    let out_schema = plugin_output_schema(&src.match_columns, &src.value_columns);
    let plugin_dir = output_cache_root.join("plugin").join(&src.plugin_name);
    std::fs::create_dir_all(&plugin_dir)
        .map_err(|e| DataFusionError::Execution(format!("mkdir {}: {e}", plugin_dir.display())))?;
    let file_name = format!("{}.parquet", canonical_chrom_label(chrom));
    let shard_path = plugin_dir.join(&file_name);

    // Stream the tiered join output straight into two per-tier temp shards
    // (warm run, cold run) instead of collecting the whole chromosome into one
    // Vec<RecordBatch>, concat_batches-ing it into a single allocation, and
    // sorting it in memory. That collect+concat+sort was the dominant remaining
    // memory cost after `assume_unique` removed the dedup HashSet, and still
    // OOM'd on the largest chromosomes (chr7, chr8, chrX for CADD). This keeps
    // peak memory at O(one batch) instead of O(whole chromosome).
    //
    // Correctness needs each tier's rows position-ascending on disk
    // (`PageDir::resolve_ranges`'s binary search assumes a sorted run). That
    // does NOT follow from `target_partitions=1` alone: `p LEFT JOIN v` plans
    // as `HashJoinExec`, and the query optimizer picks which side to hash
    // (`join_type`) based on relative table size -- for a whole-genome plugin
    // bigger than the variation probe (CADD, SpliceAI: the configs actually
    // built and shipped so far) it swaps to `join_type=Right` with the plugin
    // on the probe side, which DOES stream through unreordered; for a plugin
    // smaller than the probe (a sparse plugin, a small custom subset, a
    // low-coverage chrom) it can plan `join_type=Left` with the plugin on the
    // *build* side instead, whose row order a hash join does not preserve.
    // `assert_start_monotonic` below checks the actual invariant directly
    // rather than trusting which side the optimizer chose: a real violation
    // aborts the build loudly instead of silently writing a shard whose
    // point-lookup can miss rows that are actually present.
    let mut stream = tiered_stream(&build_ctx, &dedup_view, variation_shard).await?;

    let warm_tmp = plugin_dir.join(format!("{file_name}.warm.tmp"));
    let cold_tmp = plugin_dir.join(format!("{file_name}.cold.tmp"));
    let mut warm_writer = PluginShardWriter::create(&warm_tmp, Arc::clone(&out_schema))?;
    let mut cold_writer = PluginShardWriter::create(&cold_tmp, Arc::clone(&out_schema))?;
    let (mut warm, mut cold) = (0usize, 0usize);
    let start_idx = out_schema.index_of("start")?;
    let (mut warm_last_start, mut cold_last_start): (Option<u32>, Option<u32>) = (None, None);

    while let Some(b) = stream.next().await {
        let batch = b?;
        if batch.num_rows() == 0 {
            continue;
        }
        let reordered = reproject_cast(&batch, &out_schema)?;
        let tier_col = reordered
            .column(out_schema.index_of("tier")?)
            .as_any()
            .downcast_ref::<Int8Array>()
            .ok_or_else(|| DataFusionError::Execution("tier column not Int8".into()))?
            .clone();
        for keep in [0i8, 1i8] {
            let mask: BooleanArray = (0..reordered.num_rows())
                .map(|i| Some(tier_col.value(i) == keep))
                .collect();
            let filtered = filter_record_batch(&reordered, &mask)
                .map_err(|e| DataFusionError::Execution(format!("filter tier {keep}: {e}")))?;
            if filtered.num_rows() == 0 {
                continue;
            }
            if keep == 0 {
                warm_last_start =
                    assert_start_monotonic(&filtered, start_idx, keep, chrom, warm_last_start)?;
                warm += filtered.num_rows();
                warm_writer.write(&filtered)?;
            } else {
                cold_last_start =
                    assert_start_monotonic(&filtered, start_idx, keep, chrom, cold_last_start)?;
                cold += filtered.num_rows();
                cold_writer.write(&filtered)?;
            }
        }
    }
    warm_writer.finish()?;
    cold_writer.finish()?;
    info!(
        "plugin_cache[{}/{chrom}]: tier-join+streaming-write done, warm={warm} cold={cold}, {:.1}s elapsed",
        src.plugin_name,
        t_start.elapsed().as_secs_f64()
    );

    // Empty chrom → no shard (matches variation builder cleanup). Remove any
    // stale shard from a previous build so the manifest (rows: 0) matches disk
    // and the runtime never opens a leftover file for an empty chrom.
    if warm + cold == 0 {
        let _ = std::fs::remove_file(&warm_tmp);
        let _ = std::fs::remove_file(&cold_tmp);
        let _ = std::fs::remove_file(&shard_path);
        return Ok(ChromEntry {
            chrom: canonical_chrom_label(chrom),
            file: file_name,
            rows: 0,
            warm: 0,
            cold: 0,
        });
    }

    // Merge: final shard = warm run then cold run. Both temps are already
    // start-sorted by construction (see above), so this is a plain streaming
    // batch-by-batch copy into the final writer (which carries the point-lookup
    // page-index properties) — no sort, O(one batch) memory.
    let mut writer = PluginShardWriter::create(&shard_path, Arc::clone(&out_schema))?;
    for tmp in [&warm_tmp, &cold_tmp] {
        let file = std::fs::File::open(tmp)
            .map_err(|e| DataFusionError::Execution(format!("open temp {}: {e}", tmp.display())))?;
        let reader = parquet::arrow::arrow_reader::ParquetRecordBatchReaderBuilder::try_new(file)
            .map_err(|e| DataFusionError::Execution(format!("open temp reader: {e}")))?
            .build()
            .map_err(|e| DataFusionError::Execution(format!("build temp reader: {e}")))?;
        for batch in reader {
            let batch =
                batch.map_err(|e| DataFusionError::Execution(format!("read temp batch: {e}")))?;
            writer.write(&batch)?;
        }
        let _ = std::fs::remove_file(tmp);
    }
    let rows = writer.finish()?;
    if rows == 0 {
        let _ = std::fs::remove_file(&shard_path);
    }
    info!(
        "plugin_cache[{}/{chrom}]: done, rows={rows}, {:.1}s total",
        src.plugin_name,
        t_start.elapsed().as_secs_f64()
    );
    Ok(ChromEntry {
        chrom: canonical_chrom_label(chrom),
        file: file_name,
        rows,
        warm,
        cold,
    })
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::plugin_cache::source_manifest::SourceManifest;
    use datafusion::arrow::array::{Int8Array, StringArray, UInt32Array};
    use datafusion::arrow::datatypes::{DataType, Field, Schema};
    use parquet::arrow::ArrowWriter;
    use std::io::Write;
    use std::sync::Arc;

    fn write_gz(path: &std::path::Path, body: &str) {
        let f = std::fs::File::create(path).unwrap();
        let mut enc = flate2::write::GzEncoder::new(f, flate2::Compression::default());
        enc.write_all(body.as_bytes()).unwrap();
        enc.finish().unwrap();
    }

    /// Minimal variation-like shard: (chrom, start, allele_string, tier).
    fn write_synthetic_variation(path: &std::path::Path, rows: &[(&str, u32, &str, i8)]) {
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
    async fn builds_tiered_shard_with_counts() {
        let dir = tempfile::tempdir().unwrap();
        let tsv = dir.path().join("src.tsv.gz");
        // two rows on chr1: 100 (matches warm variation), 300 (miss)
        write_gz(&tsv, "chr1\t100\tA\tG\t0.9\nchr1\t300\tG\tA\t0.7\n");

        let var = dir.path().join("var.parquet");
        write_synthetic_variation(&var, &[("1", 100, "A/G", 0i8)]); // 100 warm (tier 0)

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
            tsv.display()
        );
        let manifest: SourceManifest = toml::from_str(&toml).unwrap();
        let out = dir.path().join("out");
        let entry = build_plugin_chrom(&manifest, "demo.source.toml", &var, &out, "1")
            .await
            .unwrap();
        assert_eq!(entry.rows, 2);
        assert_eq!(entry.warm, 1); // start 100 inherited tier 0
        assert_eq!(entry.cold, 1); // start 300 miss -> cold
        assert!(
            out.join("plugin")
                .join("demo")
                .join("chr1.parquet")
                .exists()
        );
    }

    // Two source rows sharing the runtime probe key
    // (start, allele_string, protein_variant) but different scores must collapse
    // to the FIRST in file order (VEP's first-in-file rule) — PR #190 dedup fix.
    #[tokio::test(flavor = "multi_thread")]
    async fn dedups_duplicate_aa_change_keeping_first_in_file() {
        use crate::plugin_cache::lookup::{PluginBufferSlice, PluginLookup, PluginScalar};

        let dir = tempfile::tempdir().unwrap();
        let tsv = dir.path().join("am.tsv.gz");
        // Same variant chr3:101 C>T, same aa-change H101Y, two UniProts, two
        // scores. VEP keeps the first (0.0431). A third distinct aa-change at the
        // same position (K55N) must survive.
        write_gz(
            &tsv,
            "chr3\t101\tC\tT\tH101Y\t0.0431\n\
             chr3\t101\tC\tT\tH101Y\t0.0898\n\
             chr3\t101\tC\tT\tK55N\t0.7000\n",
        );

        let var = dir.path().join("var.parquet");
        write_synthetic_variation(&var, &[("3", 101, "C/T", 0i8)]);

        let toml = format!(
            r##"
plugin_name = "am"
coordinate_system = "1-based"
ingest_sql = """
SELECT chrom, CAST(pos AS INT) AS start, CAST(pos AS INT) AS end,
       concat(ref, '/', alt) AS allele_string, pv AS protein_variant,
       CAST(score AS FLOAT) AS am_pathogenicity
FROM plugin_am_src
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
    {{ name = "pv",    type = "Utf8" }},
    {{ name = "score", type = "Utf8" }},
  ]

[[match_column]]
column = "protein_variant"
template = "{{ref_aa}}{{Protein_position}}{{alt_aa}}"

[[value_columns]]
column = "am_pathogenicity"
csq_field = "am_pathogenicity"
type = "Float32"
"##,
            tsv.display()
        );
        let manifest: SourceManifest = toml::from_str(&toml).unwrap();
        let out = dir.path().join("out");
        let entry = build_plugin_chrom(&manifest, "am.source.toml", &var, &out, "3")
            .await
            .unwrap();
        // H101Y duplicate collapsed → 2 rows survive (H101Y once + K55N).
        assert_eq!(entry.rows, 2, "duplicate H101Y row must be dropped");

        let shard = out.join("plugin").join("am").join("chr3.parquet");
        let lk = PluginLookup::open(
            &shard,
            vec!["protein_variant".into()],
            vec!["am_pathogenicity".into()],
        )
        .await
        .unwrap();
        let batch = lk.take_buffer(&[101]).await.unwrap();
        let slice = PluginBufferSlice::from_batch(&batch, 1, 1).unwrap();
        // The surviving H101Y row carries the FIRST score (0.0431), not 0.0898.
        match slice.probe(101, "C/T", &[Some("H101Y".into())]).unwrap()[0] {
            PluginScalar::F32(v) => assert!((v - 0.0431).abs() < 1e-6, "kept {v}"),
            ref other => panic!("{other:?}"),
        }
        // The distinct aa-change at the same position is untouched.
        match slice.probe(101, "C/T", &[Some("K55N".into())]).unwrap()[0] {
            PluginScalar::F32(v) => assert!((v - 0.7).abs() < 1e-6),
            ref other => panic!("{other:?}"),
        }
    }

    // C1 regression: a tiered join that does NOT preserve row order (e.g. the
    // planner puts the plugin on the hash-build side for a small-plugin/
    // large-probe query shape) must fail the build loudly rather than write a
    // shard whose on-disk (tier, start) sort claim is false.
    #[test]
    fn assert_start_monotonic_catches_within_batch_disorder() {
        let schema = Arc::new(Schema::new(vec![Field::new(
            "start",
            DataType::UInt32,
            false,
        )]));
        let batch = RecordBatch::try_new(
            schema,
            vec![Arc::new(UInt32Array::from(vec![10u32, 20, 15]))],
        )
        .unwrap();
        let err = assert_start_monotonic(&batch, 0, 1, "1", None).unwrap_err();
        assert!(err.to_string().contains("not position-ascending"));
    }

    #[test]
    fn assert_start_monotonic_catches_cross_batch_regression() {
        let schema = Arc::new(Schema::new(vec![Field::new(
            "start",
            DataType::UInt32,
            false,
        )]));
        let first = RecordBatch::try_new(
            schema.clone(),
            vec![Arc::new(UInt32Array::from(vec![100u32]))],
        )
        .unwrap();
        let last_seen = assert_start_monotonic(&first, 0, 0, "1", None).unwrap();
        assert_eq!(last_seen, Some(100));
        let second =
            RecordBatch::try_new(schema, vec![Arc::new(UInt32Array::from(vec![50u32]))]).unwrap();
        let err = assert_start_monotonic(&second, 0, 0, "1", last_seen).unwrap_err();
        assert!(err.to_string().contains("not position-ascending"));
    }

    #[test]
    fn assert_start_monotonic_accepts_sorted_input() {
        let schema = Arc::new(Schema::new(vec![Field::new(
            "start",
            DataType::UInt32,
            false,
        )]));
        let batch = RecordBatch::try_new(
            schema,
            vec![Arc::new(UInt32Array::from(vec![10u32, 10, 20, 30]))],
        )
        .unwrap();
        let last_seen = assert_start_monotonic(&batch, 0, 0, "1", None).unwrap();
        assert_eq!(last_seen, Some(30));
    }

    // --chrom M / chrM / MT must all select the MT rows (data is folded to "MT"
    // by canonical_contig), not silently produce a 0-row shard (PR #190 M1).
    #[tokio::test(flavor = "multi_thread")]
    async fn builds_mt_shard_from_any_mt_alias() {
        let dir = tempfile::tempdir().unwrap();
        let tsv = dir.path().join("mt.tsv.gz");
        write_gz(&tsv, "chrM\t100\tA\tG\t0.9\n");
        let var = dir.path().join("var.parquet");
        write_synthetic_variation(&var, &[("MT", 100, "A/G", 0i8)]);
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
            tsv.display()
        );
        let manifest: SourceManifest = toml::from_str(&toml).unwrap();
        for (i, alias) in ["M", "chrM", "MT"].iter().enumerate() {
            let out = dir.path().join(format!("out{i}"));
            let entry = build_plugin_chrom(&manifest, "demo.source.toml", &var, &out, alias)
                .await
                .unwrap();
            assert_eq!(entry.rows, 1, "alias '{alias}' should build the MT row");
        }
    }
}
