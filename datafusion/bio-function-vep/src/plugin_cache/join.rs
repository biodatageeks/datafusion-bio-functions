//! Variation-frequency join → warm/cold tier derivation (AF discarded).
//!
//! LEFT-joins the normalized plugin rows against the variation shard on the
//! shared key `(chrom, start, allele_string)`, reduces the shard's stored AF
//! columns to a single `af_max`, and appends `tier` (0 = warm if
//! `af_max >= threshold`, else 1 = cold / no match). The variation columns are
//! dropped; only the plugin's value columns plus `tier` survive.

use std::path::Path;

use datafusion::common::Result;
use datafusion::physical_plan::SendableRecordBatchStream;
use datafusion::prelude::SessionContext;

/// SQL that LEFT-joins `normalized_view` to a registered `variation_probe` table
/// exposing `(chrom, start, allele_string, af_max)` and appends `tier`. The
/// value columns of `normalized_view` pass through (`p.*`); variation columns
/// drop.
pub fn tier_sql(normalized_view: &str, variation_probe: &str, threshold: f64) -> String {
    format!(
        "SELECT p.*, \
         CASE WHEN COALESCE(v.af_max, 0.0) >= {threshold} \
         THEN CAST(0 AS TINYINT) ELSE CAST(1 AS TINYINT) END AS tier \
         FROM {normalized_view} p \
         LEFT JOIN {variation_probe} v \
         ON p.chrom = v.chrom AND p.start = v.start AND p.allele_string = v.allele_string"
    )
}

/// Register the variation shard as an `af_max`-bearing probe view, then stream
/// the tiered rows. `af_max_sql` is the expression (over the shard's columns)
/// that reduces the stored AF columns to a single `af_max`.
pub async fn tiered_stream(
    ctx: &SessionContext,
    normalized_view: &str,
    variation_shard: &Path,
    af_max_sql: &str,
    threshold: f64,
) -> Result<SendableRecordBatchStream> {
    let shard = variation_shard.to_string_lossy();
    ctx.register_parquet(
        "plugin_variation_raw",
        shard.as_ref(),
        datafusion::prelude::ParquetReadOptions::default(),
    )
    .await?;
    ctx.sql(&format!(
        "CREATE OR REPLACE VIEW plugin_variation_probe AS \
         SELECT chrom, start, allele_string, {af_max_sql} AS af_max FROM plugin_variation_raw"
    ))
    .await?;
    let df = ctx
        .sql(&tier_sql(
            normalized_view,
            "plugin_variation_probe",
            threshold,
        ))
        .await?;
    df.execute_stream().await
}

#[cfg(test)]
mod tests {
    use super::*;
    use datafusion::arrow::array::{
        Float32Array, Float64Array, Int8Array, StringArray, UInt32Array,
    };
    use datafusion::arrow::datatypes::{DataType, Field, Schema};
    use datafusion::arrow::record_batch::RecordBatch;
    use datafusion::prelude::col;
    use std::sync::Arc;

    #[tokio::test(flavor = "multi_thread")]
    async fn derives_tier_warm_cold_and_miss() {
        let ctx = SessionContext::new();

        // variation probe: (chrom,start,allele_string,af_max)
        let var = RecordBatch::try_new(
            Arc::new(Schema::new(vec![
                Field::new("chrom", DataType::Utf8, false),
                Field::new("start", DataType::UInt32, false),
                Field::new("allele_string", DataType::Utf8, false),
                Field::new("af_max", DataType::Float64, false),
            ])),
            vec![
                Arc::new(StringArray::from(vec!["1", "1"])),
                Arc::new(UInt32Array::from(vec![100u32, 200])),
                Arc::new(StringArray::from(vec!["A/G", "C/T"])),
                Arc::new(Float64Array::from(vec![0.5, 0.001])), // warm, cold
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
            .sql(&tier_sql(
                "plugin_demo_norm",
                "plugin_variation_probe",
                0.01,
            ))
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
        assert_eq!(tier.value(0), 0); // start 100, af 0.5 -> warm
        assert_eq!(tier.value(1), 1); // start 200, af 0.001 -> cold
        assert_eq!(tier.value(2), 1); // start 300, no match -> cold
        // value column preserved, variation columns dropped
        assert!(b.schema().index_of("demo_score").is_ok());
        assert!(b.schema().index_of("af_max").is_err());
    }
}
