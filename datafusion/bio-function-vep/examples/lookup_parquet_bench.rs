use std::sync::Arc;
use std::time::Instant;

use datafusion::common::{DataFusionError, Result, ScalarValue};
use datafusion::config::ConfigOptions;
use datafusion::physical_plan::displayable;
use datafusion::prelude::{ParquetReadOptions, SessionConfig, SessionContext};
use datafusion_bio_format_vcf::table_provider::VcfTableProvider;
use datafusion_bio_function::BioSessionExt;
use datafusion_bio_function_ranges::BioConfig;
use datafusion_bio_function_vep::register_vep_functions;

fn scalar_to_u64(value: ScalarValue) -> Result<u64> {
    match value {
        ScalarValue::UInt64(Some(v)) => Ok(v),
        ScalarValue::Int64(Some(v)) if v >= 0 => Ok(v as u64),
        ScalarValue::UInt32(Some(v)) => Ok(v as u64),
        ScalarValue::Int32(Some(v)) if v >= 0 => Ok(v as u64),
        other => Err(DataFusionError::Execution(format!(
            "unexpected COUNT(*) scalar value: {other:?}"
        ))),
    }
}

async fn run_case(
    vcf_path: &str,
    cache_path: &str,
    columns: &str,
    label: &str,
    target_partitions: usize,
) -> Result<()> {
    let config = SessionConfig::from(ConfigOptions::new())
        .with_option_extension(BioConfig::default())
        .with_information_schema(true)
        .with_repartition_joins(false)
        .with_target_partitions(target_partitions);
    let ctx = SessionContext::new_with_bio(config);
    datafusion_bio_function_ranges::table_function::register_ranges_functions(&ctx);
    register_vep_functions(&ctx);

    let vcf = VcfTableProvider::new(vcf_path.to_string(), Some(vec![]), Some(vec![]), None, true)?;
    ctx.register_table("vcf", Arc::new(vcf))?;

    ctx.register_parquet("var_cache", cache_path, ParquetReadOptions::default())
        .await?;

    // Count VCF rows
    let count_batches = ctx
        .sql("SELECT COUNT(*) AS cnt FROM vcf")
        .await?
        .collect()
        .await?;
    let input_rows = if let Some(batch) = count_batches.first() {
        if batch.num_rows() > 0 {
            scalar_to_u64(ScalarValue::try_from_array(batch.column(0), 0)?)?
        } else {
            0
        }
    } else {
        0
    };

    let sql = if columns.is_empty() {
        "SELECT * FROM lookup_variants('vcf', 'var_cache', '')".to_string()
    } else {
        format!("SELECT * FROM lookup_variants('vcf', 'var_cache', '{columns}')")
    };

    // Show the physical plan
    let df = ctx.sql(&sql).await?;
    let plan = df.create_physical_plan().await?;
    let plan_str = displayable(plan.as_ref()).indent(true).to_string();
    println!("\n=== Physical plan for case={label} ===");
    println!("{plan_str}");

    // Execute and time
    let start = Instant::now();
    let out_batches = ctx.sql(&sql).await?.collect().await?;
    let elapsed = start.elapsed().as_secs_f64();
    let output_rows: u64 = out_batches.iter().map(|b| b.num_rows() as u64).sum();

    println!(
        "case={label} target_partitions={target_partitions} input_rows={input_rows} output_rows={output_rows} elapsed_s={elapsed:.3} input_rows_per_s={:.1}",
        if elapsed > 0.0 {
            input_rows as f64 / elapsed
        } else {
            0.0
        }
    );

    Ok(())
}

#[tokio::main]
async fn main() -> Result<()> {
    env_logger::init();

    let args: Vec<String> = std::env::args().collect();
    let vcf_path = args.get(1).map(|s| s.as_str()).unwrap_or(
        "/Users/mwiewior/research/git/polars-bio-vep-benchmark/vep-benchmark/data/HG002_GRCh38_1_22_v4.2.1_benchmark.vcf.gz",
    );
    let cache_path = args
        .get(2)
        .map(|s| s.as_str())
        .unwrap_or("/Users/mwiewior/research/data/vep/115_GRCh38_variants.parquet");
    let target_partitions: usize = args.get(3).and_then(|s| s.parse().ok()).unwrap_or(1);

    println!("vcf={vcf_path}\ncache={cache_path}\ntarget_partitions={target_partitions}");

    run_case(
        vcf_path,
        cache_path,
        "variation_name,allele_string,clin_sig",
        "three_columns",
        target_partitions,
    )
    .await?;

    Ok(())
}
