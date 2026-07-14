//! End-to-end per-chrom plugin cache build: register sources → ingest view →
//! normalization wrapper → variation-frequency join/tier → two-pass tiered
//! shard write → cache-manifest chrom entry.

use std::path::Path;
use std::sync::Arc;

use datafusion::arrow::array::{Array, AsArray, BooleanArray, Int8Array};
use datafusion::arrow::compute::{
    CastOptions, cast_with_options, concat_batches, filter_record_batch, sort_to_indices, take,
};
use datafusion::arrow::datatypes::{DataType, SchemaRef};
use datafusion::arrow::record_batch::RecordBatch;
use datafusion::common::{DataFusionError, Result};
use datafusion::datasource::MemTable;
use datafusion::prelude::{SessionConfig, SessionContext};
use futures::StreamExt;
use log::warn;

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
///
/// The cast is DELIBERATELY unsafe (`CastOptions { safe: false }`). Arrow's default —
/// what `arrow::compute::cast` gives you — is `safe: true`, which turns a value it
/// cannot cast into NULL rather than an error. Applied to a manifest's declared
/// `[[value_columns]] type`, that made a wrong type a *silent* failure, and the worst
/// one in the build: declaring `type = "Float32"` for a categorical string column
/// produced
///
/// ```text
/// build SUCCEEDED: rows=1 warm=1 cold=0    <-- even the counts look healthy
/// probe(100,"A/G") = Some([Null])          <-- every annotation empty
/// ```
///
/// and the real AlphaMissense manifest has `am_class` (Utf8) sitting directly above
/// `am_pathogenicity` (Float32), so swapping two adjacent `type =` lines was enough.
///
/// `safe: false` makes that a hard build error naming the column and the two types.
/// It does NOT reject legitimately-null source data: a null slot is not a value that
/// failed to cast, and arrow preserves it either way (pinned by
/// `correct_types_and_null_source_values_still_build`).
fn reproject_cast(batch: &RecordBatch, schema: &SchemaRef) -> Result<RecordBatch> {
    let strict = CastOptions {
        safe: false,
        ..Default::default()
    };
    let cols = schema
        .fields()
        .iter()
        .map(|f| {
            let idx = batch.schema().index_of(f.name())?;
            let col = batch.column(idx);
            cast_with_options(col, f.data_type(), &strict).map_err(|e| {
                DataFusionError::Execution(format!(
                    "plugin value column '{}' is declared type {} but the ingested data is {} \
                     and does not fit it: {e}. Fix the column's `type` in the manifest (or make \
                     `ingest_sql` produce the declared type) — a lenient cast would have written \
                     NULL for every offending row and the plugin would annotate nothing.",
                    f.name(),
                    f.data_type(),
                    col.data_type(),
                ))
            })
        })
        .collect::<Result<Vec<_>>>()?;
    RecordBatch::try_new(Arc::clone(schema), cols)
        .map_err(|e| DataFusionError::Execution(format!("reproject: {e}")))
}

/// Reject an `allele_string` that still carries the VCF reader's pipe-joined
/// multi-ALT form (`A/G|T`).
///
/// The bio-formats VCF reader collapses a record's ALT alleles into ONE `Utf8`
/// value joined by `|` (`chr1 100 . A G,T` → one row, `alt = "G|T"`), so the
/// obvious — and spec-canonical — mapping `concat(ref, '/', alt) AS allele_string`
/// stores `A/G|T`. The runtime probes one ALT at a time (`A/G`), and the variation
/// cache's `allele_string` is likewise per-allele, so such a row matches NOTHING:
/// the build succeeds, the shard is written, and the plugin annotates nothing,
/// forever, without a warning. Every VCF-backed plugin this provider unblocks
/// (Mastermind, gnomADMt, EVE, Geno2MP) is a multi-allelic source, so this is the
/// default outcome, not an edge case.
///
/// One predicate kills the whole class: no `allele_string` written to a plugin
/// shard may contain `|`. `ingest_sql` must split the ALTs into one row per
/// allele instead (see [`VcfParams`](crate::plugin_cache::source_manifest::VcfParams)).
///
/// Checked on the normalized, chrom-filtered rows already materialized in memory
/// (not by re-querying the ingest view), so the guard costs no extra pass over
/// what can be a tens-of-GB source.
///
/// # The guard must not fail OPEN
///
/// Which Arrow string encoding `allele_string` arrives in is DataFusion's choice, not
/// the manifest's: a parquet source yields `Utf8View` (`schema_force_view_types` is on
/// by default), CSV/TSV and the VCF reader yield `Utf8`, and `LargeUtf8` is reachable
/// through an `ingest_sql` cast. The first version of this guard downcast to
/// `StringArray` alone and `continue`d otherwise — so `provider = "parquet"` skipped the
/// check outright and built the very shard the guard exists to prevent:
///
/// ```text
/// PROBE TYPES: alt=Utf8View allele_string=Utf8View
/// PROBE RESULT: guard SILENTLY SKIPPED — built ChromEntry { rows: 1, warm: 0, cold: 1 }
///               with a pipe-joined allele_string
/// ```
///
/// Hence: every string encoding is checked, and anything that is NOT a string is an
/// ERROR naming the Arrow type found. "I do not recognize this column, so I will assume
/// it is fine" is precisely the silence this whole guard is here to break.
fn reject_pipe_joined_alleles(batches: &[RecordBatch], plugin: &str) -> Result<()> {
    /// First non-null value containing `|`, over any of Arrow's string iterators.
    fn first_pipe_joined<'a>(mut values: impl Iterator<Item = Option<&'a str>>) -> Option<&'a str> {
        values.find_map(|v| v.filter(|allele| allele.contains('|')))
    }

    for batch in batches {
        let idx = batch.schema().index_of("allele_string")?;
        let col = batch.column(idx);
        let offender = match col.data_type() {
            DataType::Utf8 => first_pipe_joined(col.as_string::<i32>().iter()),
            DataType::LargeUtf8 => first_pipe_joined(col.as_string::<i64>().iter()),
            DataType::Utf8View => first_pipe_joined(col.as_string_view().iter()),
            other => {
                return Err(DataFusionError::Execution(format!(
                    "plugin '{plugin}': ingest_sql produced an allele_string of Arrow type \
                     {other}, but allele_string must be a string (Utf8, Utf8View or LargeUtf8) — \
                     it is the key the runtime probes on and the variation cache joins on. Fix \
                     the column's expression in `ingest_sql` (canonically \
                     `concat(ref, '/', alt) AS allele_string`, one row per ALT)."
                )));
            }
        };
        let Some(offender) = offender else {
            continue;
        };
        return Err(DataFusionError::Execution(format!(
            "plugin '{plugin}': ingest_sql produced an allele_string containing '|' \
             (e.g. \"{offender}\"). The VCF reader joins a multi-allelic record's ALT \
             alleles into ONE `alt` value ('G|T'), so `concat(ref, '/', alt)` yields a key \
             that can never equal the runtime probe key ('{{ref}}/{{alt}}', one ALT at a time) \
             nor the variation cache's per-allele allele_string — the shard would build \
             cleanly and annotate nothing. Split the multi-allelic ALTs in `ingest_sql` so \
             each ALT is its own row, e.g.\n  \
             SELECT chrom, `start`, `end`, concat(`ref`, '/', alt_one) AS allele_string, …\n  \
             FROM (SELECT *, unnest(string_to_array(alt, '|')) AS alt_one FROM plugin_{plugin}_src)"
        )));
    }
    Ok(())
}

/// The warning text for a chrom that wrote rows of which NOT ONE joined the variation
/// cache — or `None` when there is nothing to say.
///
/// `warm == 0` with `rows > 0` is the shared signature of every silent failure this
/// build can produce: a wrong `coordinate_system` (every start off by one), a
/// `chr1`-vs-`1` contig mismatch, an `allele_string` in the wrong shape. The variation
/// LEFT JOIN is on `(chrom, start, allele_string)`, so if not one row of a whole
/// chromosome matched, the probe key almost certainly does not have the shape the
/// runtime will ask for — and the runtime will therefore find nothing either.
///
/// A WARNING, not an error: cold rows remain probeable, and a legitimately all-cold
/// plugin (every variant absent from the variation cache) is conceivable. But the
/// number was already computed and printed, and nothing acted on it. It must not be
/// silent.
///
/// Returned as a value rather than logged in place so the condition is unit-testable
/// without installing a process-global `log` sink.
fn no_warm_rows_warning(plugin: &str, chrom: &str, rows: usize, warm: usize) -> Option<String> {
    (rows > 0 && warm == 0).then(|| {
        format!(
            "plugin '{plugin}' chrom {chrom}: {rows} rows written but warm = 0 — NOT ONE row \
             joined the variation cache on (chrom, start, allele_string). The runtime probes on \
             that same key, so it will most likely find nothing either and this plugin will \
             annotate nothing. Usual causes: contig naming (`chr1` vs `1`), the manifest's \
             coordinate_system (every start off by one), or the allele_string format (it must be \
             per-allele `{{ref}}/{{alt}}`). If this source genuinely holds no variant present in \
             the variation cache, this warning is expected."
        )
    })
}

/// Sort a batch ascending by `start` (the per-tier run order the PageDir needs).
fn sort_by_start(batch: &RecordBatch) -> Result<RecordBatch> {
    let start = batch.column(batch.schema().index_of("start")?);
    let idx = sort_to_indices(start, None, None)
        .map_err(|e| DataFusionError::Execution(format!("sort_to_indices: {e}")))?;
    let cols = batch
        .columns()
        .iter()
        .map(|c| take(c, &idx, None))
        .collect::<std::result::Result<Vec<_>, _>>()
        .map_err(|e| DataFusionError::Execution(format!("take: {e}")))?;
    RecordBatch::try_new(batch.schema(), cols)
        .map_err(|e| DataFusionError::Execution(format!("rebuild sorted: {e}")))
}

/// Build one chromosome's plugin shard. Returns the cache-manifest chrom entry.
pub async fn build_plugin_chrom(
    src: &SourceManifest,
    _source_manifest_file: &str,
    variation_shard: &Path,
    output_cache_root: &Path,
    chrom: &str,
) -> Result<ChromEntry> {
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
    let deduped = dedup_keep_first(norm_stream, &match_cols).await?;
    // The source scan is done — drop the decompressed temp before the join leg.
    drop(src_temps);

    // Fail LOUDLY on a pipe-joined multi-allelic allele_string before spending the
    // (multi-GB) tier join on rows that could never match anything.
    reject_pipe_joined_alleles(&deduped, &src.plugin_name)?;

    // Build context: default parallelism for the tier join over the (multi-GB)
    // variation shard + write. The deduped survivors are re-registered here as a
    // MemTable the join consumes in place of the raw normalized view.
    let build_ctx = SessionContext::new();
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

    // Materialize tiered rows (tier inherited from the variation cache) from the
    // deduped rows.
    let mut stream = tiered_stream(&build_ctx, &dedup_view, variation_shard).await?;
    let mut batches = Vec::new();
    while let Some(b) = stream.next().await {
        batches.push(b?);
    }

    // Empty chrom → no shard (matches variation builder cleanup). Remove any
    // stale shard from a previous build so the manifest (rows: 0) matches disk
    // and the runtime never opens a leftover file for an empty chrom.
    let non_empty: Vec<RecordBatch> = batches.into_iter().filter(|b| b.num_rows() > 0).collect();
    if non_empty.is_empty() {
        let _ = std::fs::remove_file(&shard_path);
        return Ok(ChromEntry {
            chrom: canonical_chrom_label(chrom),
            file: file_name,
            rows: 0,
            warm: 0,
            cold: 0,
        });
    }
    let joined_schema = non_empty[0].schema();
    let full = concat_batches(&joined_schema, &non_empty)
        .map_err(|e| DataFusionError::Execution(format!("concat: {e}")))?;
    let reordered = reproject_cast(&full, &out_schema)?;

    // Split into warm (0) then cold (1) runs, each start-sorted, and write.
    let tier_col = reordered
        .column(out_schema.index_of("tier")?)
        .as_any()
        .downcast_ref::<Int8Array>()
        .ok_or_else(|| DataFusionError::Execution("tier column not Int8".into()))?
        .clone();
    let (mut warm, mut cold) = (0usize, 0usize);
    let mut writer = PluginShardWriter::create(&shard_path, Arc::clone(&out_schema))?;
    for keep in [0i8, 1i8] {
        let mask: BooleanArray = (0..reordered.num_rows())
            .map(|i| Some(tier_col.value(i) == keep))
            .collect();
        let filtered = filter_record_batch(&reordered, &mask)
            .map_err(|e| DataFusionError::Execution(format!("filter tier {keep}: {e}")))?;
        if filtered.num_rows() == 0 {
            continue;
        }
        let sorted = sort_by_start(&filtered)?;
        if keep == 0 {
            warm += sorted.num_rows();
        } else {
            cold += sorted.num_rows();
        }
        writer.write(&sorted)?;
    }
    let rows = writer.finish()?;
    if rows == 0 {
        let _ = std::fs::remove_file(&shard_path);
    }
    let label = canonical_chrom_label(chrom);
    // A shard whose every row missed the variation cache is the signature of a
    // mis-declared manifest — say so instead of printing `warm=0` and moving on.
    if let Some(msg) = no_warm_rows_warning(&src.plugin_name, &label, rows, warm) {
        warn!("{msg}");
    }
    Ok(ChromEntry {
        chrom: label,
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

    /// A `provider = "parquet"` plugin SOURCE (chrom, pos, allele_string, score).
    ///
    /// Written with `Utf8` string columns — but that is NOT what the build reads back:
    /// DataFusion's parquet reader materializes them as `Utf8View`
    /// (`schema_force_view_types`, on by default), which is exactly the encoding the
    /// pipe guard used to skip.
    fn write_parquet_source(path: &std::path::Path, rows: &[(&str, u32, &str, f32)]) {
        use datafusion::arrow::array::Float32Array;

        let schema = Arc::new(Schema::new(vec![
            Field::new("chrom", DataType::Utf8, false),
            Field::new("pos", DataType::UInt32, false),
            Field::new("allele_string", DataType::Utf8, false),
            Field::new("score", DataType::Float32, false),
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
                Arc::new(Float32Array::from(
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

    // The VCF reader pipe-joins a multi-allelic record's ALTs into ONE `alt`
    // value ("G|T"), so the canonical `concat(ref, '/', alt) AS allele_string`
    // stores "A/G|T" — a key the runtime probe ("A/G", one ALT at a time) can
    // never produce. Before the guard this built cleanly (rows=1) and annotated
    // nothing, forever, with no warning. It must now be a hard build error.
    #[tokio::test(flavor = "multi_thread")]
    async fn pipe_joined_allele_string_fails_the_build() {
        let dir = tempfile::tempdir().unwrap();
        let vcf = dir.path().join("multi.vcf");
        std::fs::write(
            &vcf,
            "##fileformat=VCFv4.2\n\
             ##INFO=<ID=SCORE,Number=1,Type=Float,Description=\"demo score\">\n\
             #CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n\
             chr1\t100\t.\tA\tG,T\t.\t.\tSCORE=0.9\n",
        )
        .unwrap();
        let var = dir.path().join("var.parquet");
        write_synthetic_variation(&var, &[("1", 100, "A/G", 0i8)]);

        // The naive mapping the spec's own canonical example uses.
        let toml = format!(
            r##"
plugin_name = "demo"
coordinate_system = "1-based"
ingest_sql = """
SELECT chrom, `start`, `end`,
       concat(`ref`, '/', alt) AS allele_string,
       CAST(`SCORE` AS FLOAT) AS demo_score
FROM plugin_demo_src
"""

[[source]]
provider = "vcf"
path = "{}"
  [source.vcf]
  info_fields = ["SCORE"]

[[value_columns]]
column = "demo_score"
csq_field = "DEMO"
type = "Float32"
"##,
            vcf.display()
        );
        let manifest: SourceManifest = toml::from_str(&toml).unwrap();
        manifest
            .validate()
            .expect("the manifest is structurally valid");

        let out = dir.path().join("out");
        let err = build_plugin_chrom(&manifest, "demo.source.toml", &var, &out, "1")
            .await
            .expect_err("a pipe-joined allele_string must fail the build, not build a dead shard");
        let msg = err.to_string();
        assert!(msg.contains("demo"), "must name the plugin: {msg}");
        assert!(
            msg.contains("A/G|T"),
            "must show an offending example value: {msg}"
        );
        assert!(
            msg.contains("ingest_sql"),
            "must tell the author where to fix it: {msg}"
        );
        assert!(
            msg.to_lowercase().contains("split"),
            "must tell the author to split multi-allelic ALTs: {msg}"
        );
    }

    /// The guard above, but reached through a `provider = "parquet"` source — the hole
    /// review found in it.
    ///
    /// DataFusion's parquet reader hands string columns back as `Utf8View`, not `Utf8`
    /// (`schema_force_view_types` is on by default), so the guard's
    /// `downcast_ref::<StringArray>()` returned `None` and its `else { continue }` SKIPPED
    /// the check entirely: a pipe-joined `allele_string` sailed through and built the exact
    /// dead shard the guard exists to prevent (`rows: 1, warm: 0, cold: 1`, annotating
    /// nothing forever). A guard whose whole thesis is "fail LOUDLY" must not fail OPEN on
    /// an encoding it did not anticipate.
    ///
    /// Pipe-joined `allele_string` values reach parquet sources for real: the VCF-derived
    /// dumps plugins ship (and any parquet materialized from one with the naive
    /// `concat(ref, '/', alt)`) carry `A/G|T` verbatim.
    #[tokio::test(flavor = "multi_thread")]
    async fn pipe_joined_allele_string_from_a_parquet_source_fails_the_build() {
        let dir = tempfile::tempdir().unwrap();
        let src = dir.path().join("src.parquet");
        write_parquet_source(&src, &[("chr1", 100, "A/G|T", 0.9)]);
        let var = dir.path().join("var.parquet");
        write_synthetic_variation(&var, &[("1", 100, "A/G", 0i8)]);

        let toml = format!(
            r##"
plugin_name = "demo"
coordinate_system = "1-based"
ingest_sql = """
SELECT chrom, CAST(pos AS INT) AS start, CAST(pos AS INT) AS end,
       allele_string, CAST(score AS FLOAT) AS demo_score
FROM plugin_demo_src
"""

[[source]]
provider = "parquet"
path = "{}"

[[value_columns]]
column = "demo_score"
csq_field = "DEMO"
type = "Float32"
"##,
            src.display()
        );
        let manifest: SourceManifest = toml::from_str(&toml).unwrap();
        manifest
            .validate()
            .expect("the manifest is structurally valid");

        let out = dir.path().join("out");
        let err = build_plugin_chrom(&manifest, "demo.source.toml", &var, &out, "1")
            .await
            .expect_err(
                "a pipe-joined allele_string must fail the build whatever string encoding \
                 DataFusion chose for it — parquet gives Utf8View, not Utf8",
            );
        let msg = err.to_string();
        assert!(msg.contains("demo"), "must name the plugin: {msg}");
        assert!(
            msg.contains("A/G|T"),
            "must show an offending example value: {msg}"
        );
        assert!(
            msg.contains("ingest_sql"),
            "must tell the author where to fix it: {msg}"
        );
    }

    /// The guard's fallback arm must be an ERROR, not a `continue`.
    ///
    /// An `allele_string` that is not a string at all is a broken manifest — and it is
    /// precisely the case the guard cannot check. Skipping it silently is how the
    /// `Utf8View` hole above stayed invisible; the build must instead say what it found.
    #[tokio::test(flavor = "multi_thread")]
    async fn non_string_allele_string_fails_the_build() {
        let dir = tempfile::tempdir().unwrap();
        let tsv = dir.path().join("src.tsv.gz");
        write_gz(&tsv, "chr1\t100\tA\tG\t0.9\n");
        let var = dir.path().join("var.parquet");
        write_synthetic_variation(&var, &[("1", 100, "A/G", 0i8)]);

        // `allele_string` is an INT here — a typo'd ingest_sql that names the wrong
        // column. The shard schema would happily stringify it ("100") and the plugin
        // would then annotate nothing.
        let toml = format!(
            r##"
plugin_name = "demo"
coordinate_system = "1-based"
ingest_sql = """
SELECT chrom, CAST(pos AS INT) AS start, CAST(pos AS INT) AS end,
       CAST(pos AS INT) AS allele_string, CAST(score AS FLOAT) AS demo_score
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
        let err = build_plugin_chrom(&manifest, "demo.source.toml", &var, &out, "1")
            .await
            .expect_err("a non-string allele_string is a broken manifest, not a check to skip");
        let msg = err.to_string();
        assert!(msg.contains("demo"), "must name the plugin: {msg}");
        assert!(msg.contains("allele_string"), "must name the column: {msg}");
        assert!(
            msg.contains("Int32"),
            "must name the Arrow type actually found: {msg}"
        );
    }

    // The mirror: an `ingest_sql` that DOES split the pipe-joined ALTs into one
    // row per ALT passes the guard and builds both alleles.
    #[tokio::test(flavor = "multi_thread")]
    async fn splitting_multiallelic_alts_passes_the_guard() {
        let dir = tempfile::tempdir().unwrap();
        let vcf = dir.path().join("multi.vcf");
        std::fs::write(
            &vcf,
            "##fileformat=VCFv4.2\n\
             ##INFO=<ID=SCORE,Number=1,Type=Float,Description=\"demo score\">\n\
             #CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n\
             chr1\t100\t.\tA\tG,T\t.\t.\tSCORE=0.9\n",
        )
        .unwrap();
        let var = dir.path().join("var.parquet");
        write_synthetic_variation(&var, &[("1", 100, "A/G", 0i8)]);

        let toml = format!(
            r##"
plugin_name = "demo"
coordinate_system = "1-based"
ingest_sql = """
SELECT chrom, `start`, `end`,
       concat(`ref`, '/', alt_one) AS allele_string,
       CAST(`SCORE` AS FLOAT) AS demo_score
FROM (
  SELECT *, unnest(string_to_array(alt, '|')) AS alt_one FROM plugin_demo_src
)
"""

[[source]]
provider = "vcf"
path = "{}"
  [source.vcf]
  info_fields = ["SCORE"]

[[value_columns]]
column = "demo_score"
csq_field = "DEMO"
type = "Float32"
"##,
            vcf.display()
        );
        let manifest: SourceManifest = toml::from_str(&toml).unwrap();
        let out = dir.path().join("out");
        let entry = build_plugin_chrom(&manifest, "demo.source.toml", &var, &out, "1")
            .await
            .expect("a per-ALT split must build");
        assert_eq!(entry.rows, 2, "one row per ALT (A/G and A/T)");
        assert_eq!(entry.warm, 1, "A/G matches the variation cache");
        assert_eq!(entry.cold, 1, "A/T does not");
    }

    /// The AlphaMissense manifest declares `am_class` (Utf8) directly above
    /// `am_pathogenicity` (Float32). Swap the two `type =` lines and `reproject_cast`'s
    /// arrow `cast` — whose default is `CastOptions { safe: true }` — turned every
    /// unparseable "likely_benign" into NULL instead of erroring:
    ///
    /// ```text
    /// build SUCCEEDED: rows=1 warm=1 cold=0    <-- even the counts look healthy
    /// probe(100,"A/G") = Some([Null])          <-- every annotation empty
    /// ```
    ///
    /// Worse than the multi-allelic trap, because there the row counts at least hinted
    /// at it. A declared type the data cannot hold must be a hard build error.
    #[tokio::test(flavor = "multi_thread")]
    async fn value_column_type_mismatch_fails_the_build() {
        let dir = tempfile::tempdir().unwrap();
        let tsv = dir.path().join("am.tsv.gz");
        // The AlphaMissense column layout: a categorical class and a float score.
        write_gz(&tsv, "chr1\t100\tA\tG\tlikely_benign\t0.0431\n");
        let var = dir.path().join("var.parquet");
        write_synthetic_variation(&var, &[("1", 100, "A/G", 0i8)]);

        // `am_class` is declared Float32 — the swap. "likely_benign" is not a float.
        let toml = format!(
            r##"
plugin_name = "am"
coordinate_system = "1-based"
ingest_sql = """
SELECT chrom, CAST(pos AS INT) AS start, CAST(pos AS INT) AS end,
       concat(ref, '/', alt) AS allele_string,
       am_class AS am_class,
       CAST(score AS FLOAT) AS am_pathogenicity
FROM plugin_am_src
"""

[[source]]
provider = "tsv"
path = "{}"
  [source.csv]
  delimiter = "\t"
  has_header = false
  compression = "gzip"
  schema = [
    {{ name = "chrom",    type = "Utf8" }},
    {{ name = "pos",      type = "Utf8" }},
    {{ name = "ref",      type = "Utf8" }},
    {{ name = "alt",      type = "Utf8" }},
    {{ name = "am_class", type = "Utf8" }},
    {{ name = "score",    type = "Utf8" }},
  ]

[[value_columns]]
column = "am_class"
csq_field = "am_class"
type = "Float32"

[[value_columns]]
column = "am_pathogenicity"
csq_field = "am_pathogenicity"
type = "Float32"
"##,
            tsv.display()
        );
        let manifest: SourceManifest = toml::from_str(&toml).unwrap();
        let out = dir.path().join("out");
        let err = build_plugin_chrom(&manifest, "am.source.toml", &var, &out, "1")
            .await
            .expect_err("a value column the declared type cannot hold must fail the build");
        let msg = err.to_string();
        assert!(msg.contains("am_class"), "must name the column: {msg}");
        assert!(
            msg.contains("Utf8") && msg.contains("Float32"),
            "must name both types: {msg}"
        );
    }

    /// The guard must not fire on the AlphaMissense manifest's REAL types (Utf8 class +
    /// Float32 score), nor on a legitimately NULL source value: `safe: false` rejects
    /// values that cannot be cast, and a null is not one of them.
    #[tokio::test(flavor = "multi_thread")]
    async fn correct_types_and_null_source_values_still_build() {
        let dir = tempfile::tempdir().unwrap();
        let tsv = dir.path().join("am.tsv.gz");
        // Row 2 has an EMPTY score field -> a NULL Float32 after the ingest CAST.
        write_gz(
            &tsv,
            "chr1\t100\tA\tG\tlikely_benign\t0.0431\nchr1\t300\tC\tT\tambiguous\t\n",
        );
        let var = dir.path().join("var.parquet");
        write_synthetic_variation(&var, &[("1", 100, "A/G", 0i8)]);

        let toml = format!(
            r##"
plugin_name = "am"
coordinate_system = "1-based"
ingest_sql = """
SELECT chrom, CAST(pos AS INT) AS start, CAST(pos AS INT) AS end,
       concat(ref, '/', alt) AS allele_string,
       am_class AS am_class,
       CAST(score AS FLOAT) AS am_pathogenicity
FROM plugin_am_src
"""

[[source]]
provider = "tsv"
path = "{}"
  [source.csv]
  delimiter = "\t"
  has_header = false
  compression = "gzip"
  schema = [
    {{ name = "chrom",    type = "Utf8" }},
    {{ name = "pos",      type = "Utf8" }},
    {{ name = "ref",      type = "Utf8" }},
    {{ name = "alt",      type = "Utf8" }},
    {{ name = "am_class", type = "Utf8" }},
    {{ name = "score",    type = "Utf8" }},
  ]

[[value_columns]]
column = "am_class"
csq_field = "am_class"
type = "Utf8"

[[value_columns]]
column = "am_pathogenicity"
csq_field = "am_pathogenicity"
type = "Float32"
"##,
            tsv.display()
        );
        let manifest: SourceManifest = toml::from_str(&toml).unwrap();
        let out = dir.path().join("out");
        let entry = build_plugin_chrom(&manifest, "am.source.toml", &var, &out, "1")
            .await
            .expect("the AlphaMissense types are legitimate and must keep building");
        assert_eq!(entry.rows, 2);
        assert_eq!(entry.warm, 1);
        assert_eq!(entry.cold, 1);
    }

    /// `warm == 0` with `rows > 0` is the shared signature of EVERY silent failure in
    /// this build — a wrong coordinate_system, a `chr1`-vs-`1` contig mismatch, a
    /// mis-shaped allele_string: not one row joined the variation cache. The number was
    /// already computed and printed; nothing acted on it. It must at least warn.
    #[test]
    fn all_cold_chrom_is_detected() {
        assert!(
            no_warm_rows_warning("demo", "chr1", 1_000, 0).is_some(),
            "rows > 0 with warm == 0 must be flagged"
        );
        let msg = no_warm_rows_warning("demo", "chr1", 1_000, 0).unwrap();
        assert!(msg.contains("demo") && msg.contains("chr1"), "{msg}");
        // The usual causes, so the author has somewhere to look.
        assert!(msg.contains("contig"), "must list contig naming: {msg}");
        assert!(
            msg.contains("coordinate"),
            "must list coordinate system: {msg}"
        );
        assert!(
            msg.contains("allele_string"),
            "must list allele-string format: {msg}"
        );
    }

    #[test]
    fn a_healthy_or_empty_chrom_does_not_warn() {
        // At least one row joined -> the key format is demonstrably right.
        assert!(no_warm_rows_warning("demo", "chr1", 1_000, 1).is_none());
        // No rows at all -> a different (already visible) problem, not this one.
        assert!(no_warm_rows_warning("demo", "chr1", 0, 0).is_none());
    }

    /// End-to-end: a `chr1`-vs-`1` contig mismatch in the VARIATION shard (the plugin
    /// side is canonicalized, the variation side here is not) builds cleanly with
    /// rows > 0 and warm == 0 — and the detector fires on the returned entry.
    #[tokio::test(flavor = "multi_thread")]
    async fn contig_mismatch_builds_all_cold_and_is_flagged() {
        let dir = tempfile::tempdir().unwrap();
        let tsv = dir.path().join("src.tsv.gz");
        write_gz(&tsv, "chr1\t100\tA\tG\t0.9\nchr1\t300\tG\tA\t0.7\n");
        let var = dir.path().join("var.parquet");
        // The variation shard spells the contig "chr1"; the plugin side canonicalizes to
        // "1", so NOTHING joins — exactly the failure this warning exists to surface.
        write_synthetic_variation(&var, &[("chr1", 100, "A/G", 0i8)]);

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
            .expect("a contig mismatch is a warning, NOT an error: cold rows stay probeable");
        assert_eq!(entry.rows, 2, "the rows are written");
        assert_eq!(
            entry.warm, 0,
            "…but not one of them joined the variation cache"
        );
        assert!(
            no_warm_rows_warning(&manifest.plugin_name, &entry.chrom, entry.rows, entry.warm)
                .is_some(),
            "the all-cold chrom must be flagged: {entry:?}"
        );
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
