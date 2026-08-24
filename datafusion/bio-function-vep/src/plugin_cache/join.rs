//! Variation-tier inheritance: LEFT-joins the normalized plugin rows against the
//! variation shard on `(chrom, start, allele_string)` and inherits the variation
//! record's `tier` (`COALESCE(v.tier, 1)` — no match → cold). Variation columns
//! drop; only the plugin's value columns plus `tier` survive.

use std::path::Path;

use datafusion::common::Result;
use datafusion::physical_plan::SendableRecordBatchStream;
use datafusion::prelude::SessionContext;
use log::info;

/// SQL that LEFT-joins `normalized_view` to a registered `variation_probe`
/// exposing `(chrom, start, allele_string, tier)` and inherits `tier`. The value
/// columns of `normalized_view` pass through (`p.*`); variation columns drop.
pub fn tier_sql(normalized_view: &str, variation_probe: &str) -> String {
    format!(
        "SELECT p.*, CAST(COALESCE(v.tier, 1) AS TINYINT) AS tier \
         FROM {normalized_view} p \
         LEFT JOIN {variation_probe} v \
         ON p.chrom = v.chrom AND p.start = v.start AND p.allele_string = v.allele_string"
    )
}

/// Same join, with an explicit `ORDER BY` so the output is position-ascending
/// regardless of which side of the physical `HashJoinExec` the optimizer put
/// the plugin on (see `build.rs`'s module doc). Used only as the retry path
/// after the cheap, unsorted `tiered_stream` trips `assert_start_monotonic` --
/// the sort is real work, so it's not paid by the common case (a whole-genome
/// plugin bigger than the variation probe, which streams through unreordered
/// already).
pub fn tier_sql_sorted(normalized_view: &str, variation_probe: &str) -> String {
    format!(
        "{} ORDER BY p.start",
        tier_sql(normalized_view, variation_probe)
    )
}

/// Register the variation shard as a `tier`-bearing probe view, then stream the
/// tiered rows (tier inherited from the variation record) in whatever order
/// the join happens to produce.
pub async fn tiered_stream(
    ctx: &SessionContext,
    normalized_view: &str,
    variation_shard: &Path,
) -> Result<SendableRecordBatchStream> {
    tiered_stream_impl(ctx, normalized_view, variation_shard, false).await
}

/// Same as `tiered_stream`, but the output is guaranteed position-ascending on
/// `start` (see `tier_sql_sorted`). Costs a real sort -- use only as the
/// fallback once the unsorted stream has been shown to need it.
pub async fn tiered_stream_sorted(
    ctx: &SessionContext,
    normalized_view: &str,
    variation_shard: &Path,
) -> Result<SendableRecordBatchStream> {
    tiered_stream_impl(ctx, normalized_view, variation_shard, true).await
}

async fn tiered_stream_impl(
    ctx: &SessionContext,
    normalized_view: &str,
    variation_shard: &Path,
    sorted: bool,
) -> Result<SendableRecordBatchStream> {
    let shard = variation_shard.to_string_lossy();
    ctx.register_parquet(
        "plugin_variation_raw",
        shard.as_ref(),
        datafusion::prelude::ParquetReadOptions::default(),
    )
    .await?;
    ctx.sql(&format!(
        // GROUP BY (not DISTINCT): the variation shard can carry multiple
        // source rows for the same (chrom, start, allele_string) (e.g.
        // distinct dbSNP/COSMIC-origin entries for one variant). If those
        // duplicates ever disagree on `tier` (one warm, one cold), a plain
        // `SELECT DISTINCT chrom, start, allele_string, tier` keeps BOTH rows
        // — the join key is then non-unique on the build side, so the tier
        // LEFT JOIN below fans out the plugin row once per surviving variant,
        // silently inflating the cache with duplicate rows for that key
        // (worse: two rows with the same (start, allele_string) but different
        // tier, which the point-lookup path never expects). `MIN(tier)`
        // collapses to exactly one row per key regardless of how many
        // conflicting-tier duplicates the raw shard has — 0 (warm) wins over
        // 1 (cold), i.e. any evidence of a known/frequent variant marks it
        // warm. This also guarantees `HashJoinExec`'s build side can never
        // fan out, which is a (data-dependent) join-skew risk this GROUP BY
        // removes entirely rather than merely reduces.
        //
        // The inner JOIN against a DISTINCT projection of the plugin's own
        // keys is a semi-join emulation that restricts the aggregation to
        // keys the plugin table actually has, computed BEFORE the `GROUP BY`
        // rather than after: for a sparse plugin (far fewer rows than the
        // chromosome's full variation shard), grouping the *entire* shard
        // just to inherit tier for a handful of keys is wasted work the LEFT
        // JOIN below never uses (variation rows with no matching plugin key
        // can never survive a `p.*`-projected LEFT JOIN regardless).
        // Semantically a no-op — this narrows which rows get grouped, not
        // which key wins ties or what a hit resolves to. The DISTINCT on the
        // plugin side keeps this join's build side small and its own key
        // unique, so it can't reintroduce the fan-out risk the GROUP BY
        // above exists to remove.
        "CREATE OR REPLACE VIEW plugin_variation_probe AS \
         SELECT v.chrom, v.start, v.allele_string, MIN(v.tier) AS tier \
         FROM plugin_variation_raw v \
         INNER JOIN (SELECT DISTINCT chrom, start, allele_string FROM {normalized_view}) p \
         ON v.chrom = p.chrom AND v.start = p.start AND v.allele_string = p.allele_string \
         GROUP BY v.chrom, v.start, v.allele_string"
    ))
    .await?;
    let sql = if sorted {
        tier_sql_sorted(normalized_view, "plugin_variation_probe")
    } else {
        tier_sql(normalized_view, "plugin_variation_probe")
    };
    let df = ctx.sql(&sql).await?;
    // Opt-in plan dump. The tier join's memory behaviour depends on which side
    // of each HashJoinExec the optimizer collects, and that is decided from
    // statistics -- so it cannot be reasoned out from the SQL, it has to be
    // read off the physical plan for the actual data. Costs nothing unset.
    if std::env::var("VEPYR_EXPLAIN_TIER_JOIN").is_ok() {
        let plan = df.clone().create_physical_plan().await?;
        info!(
            "plugin_cache: tier join physical plan (sorted={sorted}):\n{}",
            datafusion::physical_plan::displayable(plan.as_ref()).indent(true)
        );
    }
    df.execute_stream().await
}

#[cfg(test)]
mod tests {
    use super::*;
    use datafusion::arrow::array::{Float32Array, Int8Array, StringArray, UInt32Array};
    use datafusion::arrow::datatypes::{DataType, Field, Schema};
    use datafusion::arrow::record_batch::RecordBatch;
    use datafusion::prelude::col;
    use futures::StreamExt;
    use parquet::arrow::ArrowWriter;
    use std::sync::Arc;

    // Codex review, PR #196: for a sparse plugin the pre-fix
    // `plugin_variation_probe` unconditionally grouped the ENTIRE variation
    // shard, including keys no plugin row will ever reference. The semi-join
    // narrowing must be a pure performance change, not a semantic one --
    // verify a large-relative-to-plugin variation shard (rows the plugin
    // never touches, PLUS a conflicting-tier duplicate at a key the plugin
    // DOES touch) still tiers exactly as before.
    #[tokio::test(flavor = "multi_thread")]
    async fn tiered_stream_ignores_variation_rows_no_plugin_key_references() {
        let dir = tempfile::tempdir().unwrap();
        let schema = Arc::new(Schema::new(vec![
            Field::new("chrom", DataType::Utf8, false),
            Field::new("start", DataType::UInt32, false),
            Field::new("allele_string", DataType::Utf8, false),
            Field::new("tier", DataType::Int8, false),
        ]));
        // 100: two conflicting-tier duplicates (must still collapse to warm).
        // 999/1000/1001: irrelevant rows no plugin key ever probes.
        let batch = RecordBatch::try_new(
            schema.clone(),
            vec![
                Arc::new(StringArray::from(vec!["1", "1", "1", "1", "1"])),
                Arc::new(UInt32Array::from(vec![100u32, 100, 999, 1000, 1001])),
                Arc::new(StringArray::from(vec!["A/G", "A/G", "G/T", "C/A", "T/G"])),
                Arc::new(Int8Array::from(vec![1i8, 0, 0, 0, 1])),
            ],
        )
        .unwrap();
        let var_path = dir.path().join("var.parquet");
        let file = std::fs::File::create(&var_path).unwrap();
        let mut w = ArrowWriter::try_new(file, schema, None).unwrap();
        w.write(&batch).unwrap();
        w.close().unwrap();

        let ctx = SessionContext::new();
        let plug = RecordBatch::try_new(
            Arc::new(Schema::new(vec![
                Field::new("chrom", DataType::Utf8, false),
                Field::new("start", DataType::UInt32, false),
                Field::new("end", DataType::UInt32, false),
                Field::new("allele_string", DataType::Utf8, false),
                Field::new("demo_score", DataType::Float32, true),
            ])),
            vec![
                Arc::new(StringArray::from(vec!["1", "1"])),
                Arc::new(UInt32Array::from(vec![100u32, 300])),
                Arc::new(UInt32Array::from(vec![100u32, 300])),
                Arc::new(StringArray::from(vec!["A/G", "G/A"])),
                Arc::new(Float32Array::from(vec![Some(0.9f32), Some(0.7)])),
            ],
        )
        .unwrap();
        ctx.register_batch("plugin_demo_norm", plug).unwrap();

        let mut stream = tiered_stream(&ctx, "plugin_demo_norm", &var_path)
            .await
            .unwrap();
        let mut rows: Vec<(u32, i8)> = Vec::new();
        while let Some(b) = stream.next().await {
            let b = b.unwrap();
            let starts = b
                .column(b.schema().index_of("start").unwrap())
                .as_any()
                .downcast_ref::<UInt32Array>()
                .unwrap();
            let tiers = b
                .column(b.schema().index_of("tier").unwrap())
                .as_any()
                .downcast_ref::<Int8Array>()
                .unwrap();
            for i in 0..b.num_rows() {
                rows.push((starts.value(i), tiers.value(i)));
            }
        }
        rows.sort();
        // 100: conflicting-tier duplicate collapses to warm (0), matching
        // pre-semi-join behavior exactly. 300: no variation match -> cold (1).
        // The three irrelevant variation rows (999/1000/1001) never surface
        // here since only the plugin's own rows are projected -- their
        // absence from `plugin_variation_raw`'s post-semi-join scan is the
        // thing under test, not something directly observable in the output.
        assert_eq!(rows, vec![(100, 0), (300, 1)]);
    }

    #[tokio::test(flavor = "multi_thread")]
    async fn inherits_tier_from_variation_and_miss_is_cold() {
        let ctx = SessionContext::new();

        // variation probe now carries `tier` directly (0=warm, 1=cold).
        let var = RecordBatch::try_new(
            Arc::new(Schema::new(vec![
                Field::new("chrom", DataType::Utf8, false),
                Field::new("start", DataType::UInt32, false),
                Field::new("allele_string", DataType::Utf8, false),
                Field::new("tier", DataType::Int8, false),
            ])),
            vec![
                Arc::new(StringArray::from(vec!["1", "1"])),
                Arc::new(UInt32Array::from(vec![100u32, 200])),
                Arc::new(StringArray::from(vec!["A/G", "C/T"])),
                Arc::new(Int8Array::from(vec![0i8, 1])), // warm, cold
            ],
        )
        .unwrap();
        ctx.register_batch("plugin_variation_probe", var).unwrap();

        // plugin ingest: warm match, cold match, and a miss.
        let plug = RecordBatch::try_new(
            Arc::new(Schema::new(vec![
                Field::new("chrom", DataType::Utf8, false),
                Field::new("start", DataType::UInt32, false),
                Field::new("end", DataType::UInt32, false),
                Field::new("allele_string", DataType::Utf8, false),
                Field::new("demo_score", DataType::Float32, true),
            ])),
            vec![
                Arc::new(StringArray::from(vec!["1", "1", "1"])),
                Arc::new(UInt32Array::from(vec![100u32, 200, 300])),
                Arc::new(UInt32Array::from(vec![100u32, 200, 300])),
                Arc::new(StringArray::from(vec!["A/G", "C/T", "G/A"])),
                Arc::new(Float32Array::from(vec![Some(0.9f32), Some(0.1), Some(0.7)])),
            ],
        )
        .unwrap();
        ctx.register_batch("plugin_demo_norm", plug).unwrap();

        let df = ctx
            .sql(&tier_sql("plugin_demo_norm", "plugin_variation_probe"))
            .await
            .unwrap()
            .sort(vec![col("start").sort(true, true)])
            .unwrap();
        let batches = df.collect().await.unwrap();
        let b = &batches[0];
        let tier = b
            .column(b.schema().index_of("tier").unwrap())
            .as_any()
            .downcast_ref::<Int8Array>()
            .unwrap();
        assert_eq!(tier.value(0), 0); // 100 -> warm (inherited)
        assert_eq!(tier.value(1), 1); // 200 -> cold (inherited)
        assert_eq!(tier.value(2), 1); // 300 -> no match -> cold
        // value column preserved, variation columns dropped
        assert!(b.schema().index_of("demo_score").is_ok());
        assert!(b.schema().index_of("tier").is_ok());
    }
}
