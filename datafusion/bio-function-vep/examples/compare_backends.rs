use std::sync::Arc;
use std::time::Instant;

use datafusion::common::{DataFusionError, Result, ScalarValue};
use datafusion::config::ConfigOptions;
use datafusion::dataframe::DataFrameWriteOptions;
use datafusion::logical_expr::Partitioning;
use datafusion::prelude::{SessionConfig, SessionContext};
use datafusion_bio_format_vcf::table_provider::VcfTableProvider;
use datafusion_bio_function_ranges::{BioConfig, BioSessionExt};
use datafusion_bio_function_vep::register_vep_functions;
use datafusion_bio_function_vep_cache::KvCacheTableProvider;

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

async fn scalar_u64(ctx: &SessionContext, sql: &str, label: &str) -> Result<u64> {
    let batches = ctx.sql(sql).await?.collect().await?;
    let Some(batch) = batches.first() else {
        return Err(DataFusionError::Execution(format!(
            "empty result for {label}",
        )));
    };
    if batch.num_rows() == 0 {
        return Err(DataFusionError::Execution(format!("no rows for {label}")));
    }
    scalar_to_u64(ScalarValue::try_from_array(batch.column(0), 0)?, label)
}

#[tokio::main]
async fn main() -> Result<()> {
    let args: Vec<String> = std::env::args().collect();
    if args.len() < 5 || args.len() > 6 {
        eprintln!(
            "Usage: {} <vcf_path> <parquet_cache_path> <fjall_cache_path> <threads> [columns]",
            args[0]
        );
        std::process::exit(2);
    }

    let vcf_path = &args[1];
    let parquet_cache_path = &args[2];
    let fjall_cache_path = &args[3];
    let threads: usize = args[4]
        .parse()
        .map_err(|e| DataFusionError::Execution(format!("invalid threads: {e}")))?;
    let cols = args
        .get(5)
        .cloned()
        .unwrap_or_else(|| "variation_name".to_string());

    let out_root = "/tmp/lookup_variants_compare";
    let out_pq = format!("{out_root}/from_parquet");
    let out_fj = format!("{out_root}/from_fjall");

    let _ = std::fs::remove_dir_all(out_root);
    std::fs::create_dir_all(out_root)
        .map_err(|e| DataFusionError::Execution(format!("create output root failed: {e}")))?;

    let mut bio = BioConfig::default();
    bio.interval_join_partitioned_left_join = true;

    let config = SessionConfig::from(ConfigOptions::new())
        .with_option_extension(bio)
        .with_information_schema(true)
        .with_target_partitions(threads)
        .with_repartition_file_scans(true)
        .with_repartition_file_min_size(0);
    let ctx = SessionContext::new_with_bio(config);
    register_vep_functions(&ctx);

    let vcf = VcfTableProvider::new(vcf_path.to_string(), Some(vec![]), Some(vec![]), None, true)?;
    ctx.register_table("vcf", Arc::new(vcf))?;
    let vcf_work = ctx
        .table("vcf")
        .await?
        .repartition(Partitioning::RoundRobinBatch(threads))?
        .into_view();
    ctx.register_table("vcf_work", vcf_work)?;

    ctx.register_parquet("cache_parquet", parquet_cache_path, Default::default())
        .await?;
    let kv = KvCacheTableProvider::open_with_cache_size(fjall_cache_path, 512 * 1024 * 1024)?;
    ctx.register_table("cache_fjall", Arc::new(kv))?;

    println!(
        "Comparing with target_partitions={threads}\nvcf={vcf_path}\nparquet={parquet_cache_path}\nfjall={fjall_cache_path}\ncolumns={cols}"
    );

    let t_total = Instant::now();

    let t_pq = Instant::now();
    let parquet_df = ctx
        .sql(&format!(
            "SELECT * FROM lookup_variants('vcf_work', 'cache_parquet', '{cols}')"
        ))
        .await?;
    parquet_df
        .write_parquet(&out_pq, DataFrameWriteOptions::new(), None)
        .await?;
    println!(
        "materialized parquet-backend output -> {out_pq} ({:.2}s)",
        t_pq.elapsed().as_secs_f64()
    );

    let t_fj = Instant::now();
    let fjall_df = ctx
        .sql(&format!(
            "SELECT * FROM lookup_variants('vcf_work', 'cache_fjall', '{cols}')"
        ))
        .await?;
    fjall_df
        .write_parquet(&out_fj, DataFrameWriteOptions::new(), None)
        .await?;
    println!(
        "materialized fjall-backend output -> {out_fj} ({:.2}s)",
        t_fj.elapsed().as_secs_f64()
    );

    let cmp_ctx =
        SessionContext::new_with_config(SessionConfig::new().with_target_partitions(threads));
    cmp_ctx
        .register_parquet("pq_out", &out_pq, Default::default())
        .await?;
    cmp_ctx
        .register_parquet("fj_out", &out_fj, Default::default())
        .await?;

    let pq_rows = scalar_u64(&cmp_ctx, "SELECT COUNT(*) FROM pq_out", "pq_out rows").await?;
    let fj_rows = scalar_u64(&cmp_ctx, "SELECT COUNT(*) FROM fj_out", "fj_out rows").await?;

    let diff_sql = "
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
          SUM(CASE WHEN COALESCE(pq.cnt, 0) <> COALESCE(fj.cnt, 0) THEN 1 ELSE 0 END) AS differing_groups,
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

    let diffs = cmp_ctx.sql(diff_sql).await?.collect().await?;
    let diff_batch = diffs
        .first()
        .ok_or_else(|| DataFusionError::Execution("empty diff output".to_string()))?;
    let differing_groups = scalar_to_u64(
        ScalarValue::try_from_array(diff_batch.column(0), 0)?,
        "differing_groups",
    )?;
    let row_delta = scalar_to_u64(
        ScalarValue::try_from_array(diff_batch.column(1), 0)?,
        "row_delta",
    )?;

    println!("pq rows: {pq_rows}");
    println!("fj rows: {fj_rows}");
    println!("differing grouped rows: {differing_groups}");
    println!("sum absolute row-count delta across groups: {row_delta}");
    println!("total elapsed: {:.2}s", t_total.elapsed().as_secs_f64());

    Ok(())
}
