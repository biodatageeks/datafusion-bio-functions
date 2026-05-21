use std::sync::Arc;
use std::time::Instant;

use datafusion::common::Result;
use datafusion::physical_plan::{ExecutionPlanProperties, displayable};
use datafusion::prelude::{SessionConfig, SessionContext};
use datafusion_bio_format_vcf::table_provider::VcfTableProvider;

#[tokio::main]
async fn main() -> Result<()> {
    let args: Vec<String> = std::env::args().collect();
    if args.len() < 3 {
        eprintln!(
            "Usage: {} <vcf_path> <cache_path> [backend: fjall|parquet] [options_json]",
            args[0]
        );
        std::process::exit(1);
    }

    let vcf_path = &args[1];
    let cache_path = &args[2];
    let backend = args.get(3).map(String::as_str).unwrap_or("fjall");
    let default_options = if backend == "fjall" {
        r#"{"partitioned":true,"use_fjall":true}"#
    } else {
        r#"{"partitioned":true}"#
    };
    let options = args.get(4).map(String::as_str).unwrap_or(default_options);
    let cache_path_sql = cache_path.replace('\'', "''");
    let options_sql = options.replace('\'', "''");

    println!("vcf={vcf_path}");
    println!("cache={cache_path}");
    println!("backend={backend}");
    println!("options={options}");
    println!();

    for target_partitions in [1, 4, 8] {
        let use_fjall = options.contains("\"use_fjall\":true")
            || options.contains("\"use_fjall\": true")
            || backend == "fjall";
        let annotation_mode = if use_fjall && target_partitions > 1 {
            "parallel-fjall-annotation"
        } else if use_fjall {
            "single-fjall-annotation"
        } else {
            "serial-parquet-annotation"
        };
        let config = SessionConfig::new().with_target_partitions(target_partitions);
        let ctx = SessionContext::new_with_config(config);
        datafusion_bio_function_vep::register_vep_functions(&ctx);

        let vcf_provider = VcfTableProvider::new(vcf_path.to_string(), None, None, None, false)?;
        ctx.register_table("vcf", Arc::new(vcf_provider))?;

        let sql = format!(
            "SELECT chrom, start, `CSQ` FROM annotate_vep('vcf', '{cache_path_sql}', '{backend}', '{options_sql}')"
        );
        let df = ctx.sql(&sql).await?;
        let plan = df.clone().create_physical_plan().await?;
        let output_partitions = plan.output_partitioning().partition_count();
        let plan_text = displayable(plan.as_ref()).indent(true).to_string();

        let started = Instant::now();
        let batches = df.collect().await?;
        let elapsed = started.elapsed();
        let rows: usize = batches.iter().map(|batch| batch.num_rows()).sum();

        println!(
            "mode={annotation_mode} target_partitions={target_partitions} output_partitions={output_partitions} rows={rows} elapsed_s={:.3} rows_per_s={:.1}",
            elapsed.as_secs_f64(),
            rows as f64 / elapsed.as_secs_f64()
        );
        println!("{plan_text}");
    }

    Ok(())
}
