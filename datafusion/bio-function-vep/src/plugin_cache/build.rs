//! End-to-end per-chrom plugin cache build: register sources → ingest view →
//! normalization wrapper → variation-frequency join/tier → two-pass tiered
//! shard write → cache-manifest chrom entry.

use std::path::Path;
use std::sync::Arc;

use datafusion::arrow::array::{Array, BooleanArray, Int8Array};
use datafusion::arrow::compute::{
    cast, concat_batches, filter_record_batch, sort_to_indices, take,
};
use datafusion::arrow::datatypes::SchemaRef;
use datafusion::arrow::record_batch::RecordBatch;
use datafusion::common::{DataFusionError, Result};
use datafusion::prelude::SessionContext;
use futures::StreamExt;

use crate::cache::manifest::canonical_chrom_label;
use crate::plugin_cache::cache_manifest::ChromEntry;
use crate::plugin_cache::join::tiered_stream;
use crate::plugin_cache::normalize::{canonical_contig_udf, wrap_normalization};
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
    af_max_sql: &str,
    output_cache_root: &Path,
    chrom: &str,
) -> Result<ChromEntry> {
    let ctx = SessionContext::new();
    ctx.register_udf(canonical_contig_udf());
    register_sources(&ctx, src).await?;

    // Ingest view (raw column mapping), then normalized view (contig + coords).
    ctx.sql(&format!(
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
    let source_chrom = chrom.strip_prefix("chr").unwrap_or(chrom);
    ctx.sql(&format!(
        "CREATE OR REPLACE VIEW {norm_view} AS SELECT * FROM ({norm_sql}) WHERE chrom = '{source_chrom}'"
    ))
    .await?;

    let out_schema = plugin_output_schema(&src.match_columns, &src.value_columns);
    let plugin_dir = output_cache_root.join("plugin").join(&src.plugin_name);
    std::fs::create_dir_all(&plugin_dir)
        .map_err(|e| DataFusionError::Execution(format!("mkdir {}: {e}", plugin_dir.display())))?;
    let file_name = format!("{}.parquet", canonical_chrom_label(chrom));
    let shard_path = plugin_dir.join(&file_name);

    // Materialize tiered rows.
    let mut stream = tiered_stream(
        &ctx,
        &norm_view,
        variation_shard,
        af_max_sql,
        src.tier.threshold,
    )
    .await?;
    let mut batches = Vec::new();
    while let Some(b) = stream.next().await {
        batches.push(b?);
    }

    // Empty chrom → no shard (matches variation builder cleanup).
    let non_empty: Vec<RecordBatch> = batches.into_iter().filter(|b| b.num_rows() > 0).collect();
    if non_empty.is_empty() {
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
    use datafusion::arrow::array::{Float64Array, StringArray, UInt32Array};
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

    /// Minimal variation-like shard: (chrom, start, allele_string, minor_allele_freq).
    fn write_synthetic_variation(path: &std::path::Path, rows: &[(&str, u32, &str, f64)]) {
        let schema = Arc::new(Schema::new(vec![
            Field::new("chrom", DataType::Utf8, false),
            Field::new("start", DataType::UInt32, false),
            Field::new("allele_string", DataType::Utf8, false),
            Field::new("minor_allele_freq", DataType::Float64, false),
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
                Arc::new(Float64Array::from(
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
        write_synthetic_variation(&var, &[("1", 100, "A/G", 0.5)]);

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

[tier]
threshold = 0.01
unmatched = "cold"
"##,
            tsv.display()
        );
        let manifest: SourceManifest = toml::from_str(&toml).unwrap();
        let out = dir.path().join("out");
        let entry = build_plugin_chrom(
            &manifest,
            "demo.source.toml",
            &var,
            "minor_allele_freq",
            &out,
            "1",
        )
        .await
        .unwrap();
        assert_eq!(entry.rows, 2);
        assert_eq!(entry.warm, 1); // start 100 matched af 0.5
        assert_eq!(entry.cold, 1); // start 300 miss
        assert!(
            out.join("plugin")
                .join("demo")
                .join("chr1.parquet")
                .exists()
        );
    }
}
