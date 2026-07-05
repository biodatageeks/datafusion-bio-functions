//! Shared contig/coordinate normalization applied by the builder over the
//! ingest view. Emits `chrom`/`start`/`end`/`allele_string` in the exact form
//! the variation cache uses (bare Ensembl contig, 1-based coordinates), followed
//! by the plugin's value columns. Written once, applied identically to every
//! plugin, so no per-source SQL re-implements contig/coordinate rules.

use std::sync::Arc;

use datafusion::arrow::array::{ArrayRef, StringArray};
use datafusion::arrow::datatypes::DataType;
use datafusion::common::Result;
use datafusion::logical_expr::{ColumnarValue, ScalarUDF, Volatility};
use datafusion::prelude::create_udf;

use crate::plugin_cache::source_manifest::CoordinateSystem;

/// Bare Ensembl contig form matching the variation `chrom` column
/// (`1`..`22`, `X`, `Y`, `MT`): strip any `chr` prefix, fold the mitochondrial
/// aliases to `MT`, uppercase the rest. Mirrors `cache::key_encoding` prefix
/// handling so plugin and variation contigs agree by construction.
pub fn canonical_contig_str(raw: &str) -> String {
    let bare = raw.strip_prefix("chr").unwrap_or(raw);
    match bare {
        "M" | "MT" => "MT".to_string(),
        other => other.to_ascii_uppercase(),
    }
}

/// A DataFusion scalar UDF wrapping [`canonical_contig_str`].
pub fn canonical_contig_udf() -> ScalarUDF {
    let f = Arc::new(|args: &[ColumnarValue]| -> Result<ColumnarValue> {
        let arrays = ColumnarValue::values_to_arrays(args)?;
        let input = arrays[0]
            .as_any()
            .downcast_ref::<StringArray>()
            .expect("canonical_contig expects a Utf8 argument");
        let out: StringArray = input.iter().map(|v| v.map(canonical_contig_str)).collect();
        Ok(ColumnarValue::Array(Arc::new(out) as ArrayRef))
    });
    create_udf(
        "canonical_contig",
        vec![DataType::Utf8],
        DataType::Utf8,
        Volatility::Immutable,
        f,
    )
}

/// Build the normalization SQL over `inner_view`: canonicalize `chrom`, cast
/// `start`/`end` (with the coordinate shift), keep `allele_string`, then append
/// each value column verbatim. Value columns are enumerated explicitly (no
/// `SELECT *`/`EXCLUDE`) so the projection is stable across DataFusion versions.
pub fn wrap_normalization(
    inner_view: &str,
    coord: CoordinateSystem,
    value_columns: &[String],
) -> String {
    let start_expr = match coord {
        CoordinateSystem::OneBased => "CAST(start AS BIGINT)".to_string(),
        CoordinateSystem::ZeroBasedHalfOpen => "CAST(start AS BIGINT) + 1".to_string(),
    };
    let mut projection = format!(
        "canonical_contig(chrom) AS chrom, {start_expr} AS start, CAST(\"end\" AS BIGINT) AS \"end\", allele_string"
    );
    for col in value_columns {
        projection.push_str(&format!(", {col}"));
    }
    format!("SELECT {projection} FROM {inner_view}")
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn canonicalizes_contig_strings() {
        assert_eq!(canonical_contig_str("chr1"), "1");
        assert_eq!(canonical_contig_str("1"), "1");
        assert_eq!(canonical_contig_str("chrX"), "X");
        assert_eq!(canonical_contig_str("chrM"), "MT");
        assert_eq!(canonical_contig_str("M"), "MT");
        assert_eq!(canonical_contig_str("chrMT"), "MT");
        assert_eq!(canonical_contig_str("MT"), "MT");
    }

    #[test]
    fn one_based_passes_through_zero_based_shifts() {
        let vals = vec!["demo_score".to_string()];
        let one = wrap_normalization("plugin_demo_ingest", CoordinateSystem::OneBased, &vals);
        assert!(one.contains("canonical_contig(chrom) AS chrom"));
        assert!(one.contains("CAST(start AS BIGINT) AS start"));
        assert!(one.contains(", demo_score"));
        assert!(one.ends_with("FROM plugin_demo_ingest"));

        let zero = wrap_normalization(
            "plugin_demo_ingest",
            CoordinateSystem::ZeroBasedHalfOpen,
            &vals,
        );
        assert!(zero.contains("CAST(start AS BIGINT) + 1 AS start"));
    }

    #[tokio::test(flavor = "multi_thread")]
    async fn udf_executes_in_sql() {
        use datafusion::prelude::SessionContext;
        let ctx = SessionContext::new();
        ctx.register_udf(canonical_contig_udf());
        let batches = ctx
            .sql("SELECT canonical_contig('chr22') AS c, canonical_contig('chrM') AS m")
            .await
            .unwrap()
            .collect()
            .await
            .unwrap();
        let b = &batches[0];
        let c = b.column(0).as_any().downcast_ref::<StringArray>().unwrap();
        let m = b.column(1).as_any().downcast_ref::<StringArray>().unwrap();
        assert_eq!(c.value(0), "22");
        assert_eq!(m.value(0), "MT");
    }
}
