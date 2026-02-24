use std::sync::Arc;
use std::time::Instant;

use datafusion::common::{DataFusionError, Result, ScalarValue};
use datafusion::prelude::{SessionConfig, SessionContext};
use datafusion_bio_format_vcf::table_provider::VcfTableProvider;
use datafusion_bio_function_vep::register_vep_functions;
#[cfg(feature = "kv-cache")]
use datafusion_bio_function_vep_cache::KvCacheTableProvider;

#[derive(Clone, Debug)]
struct Case {
    label: &'static str,
    columns: Option<&'static str>,
    projection: Option<&'static str>,
}

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

async fn run_vcf_only(vcf_path: &str, iteration: usize) -> Result<(u64, f64)> {
    let config = SessionConfig::new().with_target_partitions(1);
    let ctx = SessionContext::new_with_config(config);

    let vcf = VcfTableProvider::new(vcf_path.to_string(), Some(vec![]), Some(vec![]), None, true)?;
    ctx.register_table("vcf", Arc::new(vcf))?;

    let start = Instant::now();
    let batches = ctx.sql("SELECT * FROM vcf").await?.collect().await?;
    let elapsed = start.elapsed().as_secs_f64();
    let rows: u64 = batches.iter().map(|b| b.num_rows() as u64).sum();
    println!(
        "case=vcf_only iter={} rows={} elapsed_s={:.3} rows_per_s={:.1}",
        iteration,
        rows,
        elapsed,
        if elapsed > 0.0 {
            rows as f64 / elapsed
        } else {
            0.0
        }
    );
    Ok((rows, elapsed))
}

async fn run_case(
    vcf_path: &str,
    cache_path: &str,
    cache_size_mb: u64,
    case: &Case,
    iteration: usize,
) -> Result<(u64, u64, f64)> {
    let config = SessionConfig::new().with_target_partitions(1);
    let ctx = SessionContext::new_with_config(config);
    register_vep_functions(&ctx);

    let vcf = VcfTableProvider::new(vcf_path.to_string(), Some(vec![]), Some(vec![]), None, true)?;
    ctx.register_table("vcf", Arc::new(vcf))?;

    #[cfg(feature = "kv-cache")]
    {
        let cache =
            KvCacheTableProvider::open_with_cache_size(cache_path, cache_size_mb * 1024 * 1024)?;
        ctx.register_table("var_cache", Arc::new(cache))?;
    }

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

    let sql = match case.columns {
        Some(cols) => {
            if let Some(projection) = case.projection {
                format!("SELECT {projection} FROM lookup_variants('vcf', 'var_cache', '{cols}')")
            } else {
                format!("SELECT * FROM lookup_variants('vcf', 'var_cache', '{cols}')")
            }
        }
        None => {
            if let Some(projection) = case.projection {
                format!("SELECT {projection} FROM lookup_variants('vcf', 'var_cache')")
            } else {
                "SELECT * FROM lookup_variants('vcf', 'var_cache')".to_string()
            }
        }
    };

    let start = Instant::now();
    let out_batches = ctx.sql(&sql).await?.collect().await?;
    let elapsed = start.elapsed().as_secs_f64();
    let output_rows: u64 = out_batches.iter().map(|b| b.num_rows() as u64).sum();

    println!(
        "case={} iter={} input_rows={} output_rows={} elapsed_s={:.3} input_rows_per_s={:.1} output_rows_per_s={:.1}",
        case.label,
        iteration,
        input_rows,
        output_rows,
        elapsed,
        if elapsed > 0.0 {
            input_rows as f64 / elapsed
        } else {
            0.0
        },
        if elapsed > 0.0 {
            output_rows as f64 / elapsed
        } else {
            0.0
        },
    );

    Ok((input_rows, output_rows, elapsed))
}

#[tokio::main]
async fn main() -> Result<()> {
    let args: Vec<String> = std::env::args().collect();
    let vcf_path = args.get(1).map(|s| s.as_str()).unwrap_or(
        "/Users/mwiewior/research/git/polars-bio-vep-benchmark/vep-benchmark/data/HG002_chr22.vcf.gz",
    );
    let cache_path = args
        .get(2)
        .map(|s| s.as_str())
        .unwrap_or("/Users/mwiewior/research/data/vep/variation_fjall");
    let iterations: usize = args.get(3).and_then(|s| s.parse().ok()).unwrap_or(3);
    let cache_size_mb: u64 = args.get(4).and_then(|s| s.parse().ok()).unwrap_or(512);

    let cases = vec![
        Case {
            label: "no_cache_columns",
            columns: Some(""),
            projection: None,
        },
        Case {
            label: "all_cache_columns_default",
            columns: None,
            projection: None,
        },
        Case {
            label: "three_columns",
            columns: Some("variation_name,allele_string,clin_sig"),
            projection: None,
        },
        Case {
            label: "three_columns_projected_only",
            columns: Some("variation_name,allele_string,clin_sig"),
            projection: Some("cache_variation_name, cache_allele_string, cache_clin_sig"),
        },
    ];

    println!(
        "vcf={} cache={} target_partitions=1 iterations={} cache_size_mb={}",
        vcf_path, cache_path, iterations, cache_size_mb
    );

    println!("\n=== case=vcf_only ===");
    let mut vcf_total_rows = 0u64;
    let mut vcf_total_elapsed = 0.0f64;
    for iter in 1..=iterations {
        let (rows, elapsed) = run_vcf_only(vcf_path, iter).await?;
        vcf_total_rows += rows;
        vcf_total_elapsed += elapsed;
    }
    println!(
        "avg case=vcf_only avg_rows_per_s={:.1}",
        if vcf_total_elapsed > 0.0 {
            vcf_total_rows as f64 / vcf_total_elapsed
        } else {
            0.0
        }
    );

    for case in &cases {
        println!("\n=== case={} ===", case.label);
        let mut total_input = 0u64;
        let mut total_output = 0u64;
        let mut total_elapsed = 0.0f64;
        for iter in 1..=iterations {
            let (input_rows, output_rows, elapsed) =
                run_case(vcf_path, cache_path, cache_size_mb, case, iter).await?;
            total_input += input_rows;
            total_output += output_rows;
            total_elapsed += elapsed;
        }
        let avg_input_rows_per_s = if total_elapsed > 0.0 {
            total_input as f64 / total_elapsed
        } else {
            0.0
        };
        let avg_output_rows_per_s = if total_elapsed > 0.0 {
            total_output as f64 / total_elapsed
        } else {
            0.0
        };
        println!(
            "avg case={} avg_input_rows_per_s={:.1} avg_output_rows_per_s={:.1}",
            case.label, avg_input_rows_per_s, avg_output_rows_per_s
        );
    }

    Ok(())
}
