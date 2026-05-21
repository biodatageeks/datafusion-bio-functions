use std::sync::Arc;
use std::time::Instant;

use datafusion::common::Result;
use datafusion::physical_plan::{ExecutionPlanProperties, displayable};
use datafusion::prelude::{SessionConfig, SessionContext};
use datafusion_bio_format_vcf::table_provider::VcfTableProvider;
use futures::StreamExt;

#[tokio::main]
async fn main() -> Result<()> {
    let args: Vec<String> = std::env::args().collect();
    if args.len() < 3 {
        eprintln!(
            "Usage: {} <vcf_path> <cache_path> [backend: fjall|parquet] [options_json] [target_partitions_csv] [limit] [--show-plan]",
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
    let target_partitions_list: Vec<usize> = args
        .get(5)
        .filter(|s| !s.starts_with("--"))
        .map(|s| {
            s.split(',')
                .filter_map(|part| part.trim().parse::<usize>().ok())
                .filter(|part| *part > 0)
                .collect()
        })
        .unwrap_or_else(|| vec![1, 2, 4, 8, 16]);
    let limit = args
        .get(6)
        .filter(|s| !s.starts_with("--"))
        .and_then(|s| s.parse::<usize>().ok());
    let show_plan = args.iter().any(|arg| arg == "--show-plan");
    let cache_path_sql = cache_path.replace('\'', "''");
    let options_sql = options.replace('\'', "''");
    let limit_clause = limit.map(|n| format!(" LIMIT {n}")).unwrap_or_default();

    println!("vcf={vcf_path}");
    println!("cache={cache_path}");
    println!("backend={backend}");
    println!("options={options}");
    println!("target_partitions={target_partitions_list:?}");
    println!(
        "limit={}",
        limit
            .map(|n| n.to_string())
            .unwrap_or_else(|| "none".to_string())
    );
    println!();

    for target_partitions in target_partitions_list {
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
            "SELECT chrom, start, `CSQ` FROM annotate_vep('vcf', '{cache_path_sql}', '{backend}', '{options_sql}'){limit_clause}"
        );
        let df = ctx.sql(&sql).await?;
        let plan = df.clone().create_physical_plan().await?;
        let output_partitions = plan.output_partitioning().partition_count();
        let plan_text = if show_plan {
            Some(displayable(plan.as_ref()).indent(true).to_string())
        } else {
            None
        };
        drop(plan);

        let started = Instant::now();
        let mut stream = df.execute_stream().await?;
        let mut batches = 0usize;
        let mut rows = 0usize;
        while let Some(batch_result) = stream.next().await {
            let batch = batch_result?;
            batches += 1;
            rows += batch.num_rows();
        }
        let elapsed = started.elapsed();

        println!(
            "mode={annotation_mode} target_partitions={target_partitions} output_partitions={output_partitions} batches={batches} rows={rows} elapsed_s={:.3} rows_per_s={:.1}",
            elapsed.as_secs_f64(),
            rows as f64 / elapsed.as_secs_f64()
        );
        if let Some(plan_text) = &plan_text {
            println!("{plan_text}");
        }
    }

    Ok(())
}
