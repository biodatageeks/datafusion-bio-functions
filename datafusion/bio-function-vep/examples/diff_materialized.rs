use datafusion::common::{DataFusionError, Result, ScalarValue};
use datafusion::prelude::SessionContext;

fn scalar_to_u64(value: ScalarValue, label: &str) -> Result<u64> {
    match value {
        ScalarValue::UInt64(Some(v)) => Ok(v),
        ScalarValue::Int64(Some(v)) if v >= 0 => Ok(v as u64),
        ScalarValue::UInt32(Some(v)) => Ok(v as u64),
        ScalarValue::Int32(Some(v)) if v >= 0 => Ok(v as u64),
        other => Err(DataFusionError::Execution(format!(
            "unexpected scalar for {label}: {other:?}",
        ))),
    }
}

#[tokio::main]
async fn main() -> Result<()> {
    let args: Vec<String> = std::env::args().collect();
    if args.len() != 3 {
        eprintln!("Usage: {} <pq_dir> <fj_dir>", args[0]);
        std::process::exit(2);
    }

    let ctx = SessionContext::new();
    ctx.register_parquet("pq_out", &args[1], Default::default())
        .await?;
    ctx.register_parquet("fj_out", &args[2], Default::default())
        .await?;

    let summary_sql = "
        WITH pq AS (
          SELECT replace(chrom, 'chr', '') AS chrom_norm, start, \"end\", ref, alt, cache_variation_name, COUNT(*) AS cnt
          FROM pq_out
          GROUP BY replace(chrom, 'chr', ''), start, \"end\", ref, alt, cache_variation_name
        ),
        fj AS (
          SELECT replace(chrom, 'chr', '') AS chrom_norm, start, \"end\", ref, alt, cache_variation_name, COUNT(*) AS cnt
          FROM fj_out
          GROUP BY replace(chrom, 'chr', ''), start, \"end\", ref, alt, cache_variation_name
        )
        SELECT
          COUNT(*) AS total_groups,
          SUM(CASE WHEN pq.cnt IS NULL THEN 1 ELSE 0 END) AS only_fj_groups,
          SUM(CASE WHEN fj.cnt IS NULL THEN 1 ELSE 0 END) AS only_pq_groups,
          SUM(CASE WHEN pq.cnt IS NOT NULL AND fj.cnt IS NOT NULL AND pq.cnt <> fj.cnt THEN 1 ELSE 0 END) AS count_mismatch_groups,
          SUM(ABS(COALESCE(pq.cnt, 0) - COALESCE(fj.cnt, 0))) AS row_delta
        FROM pq
        FULL OUTER JOIN fj
          ON pq.chrom_norm = fj.chrom_norm
         AND pq.start = fj.start
         AND pq.\"end\" = fj.\"end\"
         AND pq.ref = fj.ref
         AND pq.alt = fj.alt
         AND pq.cache_variation_name IS NOT DISTINCT FROM fj.cache_variation_name
    ";

    let summary_batches = ctx.sql(summary_sql).await?.collect().await?;
    let b = summary_batches
        .first()
        .ok_or_else(|| DataFusionError::Execution("empty summary".to_string()))?;

    let total_groups = scalar_to_u64(ScalarValue::try_from_array(b.column(0), 0)?, "total_groups")?;
    let only_fj_groups = scalar_to_u64(
        ScalarValue::try_from_array(b.column(1), 0)?,
        "only_fj_groups",
    )?;
    let only_pq_groups = scalar_to_u64(
        ScalarValue::try_from_array(b.column(2), 0)?,
        "only_pq_groups",
    )?;
    let count_mismatch_groups = scalar_to_u64(
        ScalarValue::try_from_array(b.column(3), 0)?,
        "count_mismatch_groups",
    )?;
    let row_delta = scalar_to_u64(ScalarValue::try_from_array(b.column(4), 0)?, "row_delta")?;

    println!(
        "summary: total_groups={total_groups} only_fj={only_fj_groups} only_pq={only_pq_groups} count_mismatch={count_mismatch_groups} row_delta={row_delta}"
    );

    let sample_sql = "
        WITH pq AS (
          SELECT replace(chrom, 'chr', '') AS chrom_norm, start, \"end\", ref, alt, cache_variation_name, COUNT(*) AS cnt
          FROM pq_out
          GROUP BY replace(chrom, 'chr', ''), start, \"end\", ref, alt, cache_variation_name
        ),
        fj AS (
          SELECT replace(chrom, 'chr', '') AS chrom_norm, start, \"end\", ref, alt, cache_variation_name, COUNT(*) AS cnt
          FROM fj_out
          GROUP BY replace(chrom, 'chr', ''), start, \"end\", ref, alt, cache_variation_name
        )
        SELECT
          COALESCE(pq.chrom_norm, fj.chrom_norm) AS chrom,
          COALESCE(pq.start, fj.start) AS start,
          COALESCE(pq.\"end\", fj.\"end\") AS \"end\",
          COALESCE(pq.ref, fj.ref) AS ref,
          COALESCE(pq.alt, fj.alt) AS alt,
          COALESCE(pq.cache_variation_name, fj.cache_variation_name) AS cache_variation_name,
          COALESCE(pq.cnt, 0) AS pq_cnt,
          COALESCE(fj.cnt, 0) AS fj_cnt,
          ABS(COALESCE(pq.cnt, 0) - COALESCE(fj.cnt, 0)) AS delta
        FROM pq
        FULL OUTER JOIN fj
          ON pq.chrom_norm = fj.chrom_norm
         AND pq.start = fj.start
         AND pq.\"end\" = fj.\"end\"
         AND pq.ref = fj.ref
         AND pq.alt = fj.alt
         AND pq.cache_variation_name IS NOT DISTINCT FROM fj.cache_variation_name
        WHERE COALESCE(pq.cnt, 0) <> COALESCE(fj.cnt, 0)
        ORDER BY delta DESC, chrom, start
        LIMIT 40
    ";

    let sample_batches = ctx.sql(sample_sql).await?.collect().await?;
    println!("sample mismatches:");
    for batch in sample_batches {
        println!("{batch:?}");
    }

    Ok(())
}
