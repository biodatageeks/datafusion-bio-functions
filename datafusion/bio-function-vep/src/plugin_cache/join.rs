//! Variation-tier inheritance: LEFT-joins the normalized plugin rows against the
//! variation shard on `(chrom, start, allele_string)` and inherits the variation
//! record's `tier` (`COALESCE(v.tier, 1)` — no match → cold). Variation columns
//! drop; only the plugin's value columns plus `tier` survive.

use std::path::Path;

use datafusion::common::Result;
use datafusion::physical_plan::SendableRecordBatchStream;
use datafusion::prelude::SessionContext;

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

/// Register the variation shard as a `tier`-bearing probe view, then stream the
/// tiered rows (tier inherited from the variation record).
pub async fn tiered_stream(
    ctx: &SessionContext,
    normalized_view: &str,
    variation_shard: &Path,
) -> Result<SendableRecordBatchStream> {
    let shard = variation_shard.to_string_lossy();
    ctx.register_parquet(
        "plugin_variation_raw",
        shard.as_ref(),
        datafusion::prelude::ParquetReadOptions::default(),
    )
    .await?;
    ctx.sql(
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
        "CREATE OR REPLACE VIEW plugin_variation_probe AS \
         SELECT chrom, start, allele_string, MIN(tier) AS tier \
         FROM plugin_variation_raw GROUP BY chrom, start, allele_string",
    )
    .await?;
    let df = ctx
        .sql(&tier_sql(normalized_view, "plugin_variation_probe"))
        .await?;
    df.execute_stream().await
}

#[cfg(test)]
mod tests {
    use super::*;
    use datafusion::arrow::array::{Float32Array, Int8Array, StringArray, UInt32Array};
    use datafusion::arrow::datatypes::{DataType, Field, Schema};
    use datafusion::arrow::record_batch::RecordBatch;
    use datafusion::prelude::col;
    use std::sync::Arc;

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
