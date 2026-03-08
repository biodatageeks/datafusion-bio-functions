use std::sync::Arc;
use std::time::Instant;

use datafusion::common::Result;
use datafusion::physical_plan::{ExecutionPlan, ExecutionPlanProperties};
use datafusion::prelude::SessionContext;
use datafusion_bio_format_vcf::table_provider::VcfTableProvider;
use futures::StreamExt;

use datafusion_bio_function_vep::kv_cache::cache_exec::{KvLookupExec, KvMatchMode};
use datafusion_bio_function_vep::kv_cache::kv_store::VepKvStore;

fn allele_matches(_vcf_ref: &str, vcf_alt: &str, allele_string: &str) -> bool {
    let mut parts = allele_string.split('/');
    let Some(_cache_ref) = parts.next() else {
        return false;
    };
    parts.any(|a| a == vcf_alt)
}

#[tokio::main]
async fn main() -> Result<()> {
    let args: Vec<String> = std::env::args().collect();
    if args.len() < 3 {
        eprintln!(
            "Usage: {} <fjall_cache_path> <vcf_path> [threads: 1,2,4,8] [cache_size_mb: 256]",
            args[0]
        );
        std::process::exit(1);
    }

    let cache_path = &args[1];
    let vcf_path = &args[2];
    let thread_list: Vec<usize> = if let Some(t) = args.get(3) {
        t.split(',').filter_map(|s| s.trim().parse().ok()).collect()
    } else {
        vec![1, 2, 4, 8]
    };
    let cache_size_mb: u64 = args.get(4).and_then(|s| s.parse().ok()).unwrap_or(256);
    let cache_size_bytes = cache_size_mb * 1024 * 1024;

    let store = Arc::new(VepKvStore::open_with_cache_size(
        cache_path,
        cache_size_bytes,
    )?);
    let cache_schema = store.schema();
    let cache_columns: Vec<String> = cache_schema
        .fields()
        .iter()
        .filter(|f| f.name() != "chrom" && f.name() != "start" && f.name() != "end")
        .map(|f| f.name().clone())
        .collect();

    println!(
        "Cache: format=v{}, columns={}, cache_size={}MB, path={}",
        store.format_version(),
        cache_columns.len(),
        cache_size_mb,
        cache_path,
    );
    println!("VCF: {}", vcf_path);
    println!();
    println!(
        "{:>7} | {:>5} | {:>10} | {:>6} | {:>12}",
        "threads", "parts", "rows", "time", "rows/s"
    );
    println!("{}", "-".repeat(55));

    for &num_threads in &thread_list {
        let config = datafusion::prelude::SessionConfig::new().with_target_partitions(num_threads);
        let ctx = SessionContext::new_with_config(config);

        let table =
            VcfTableProvider::new(vcf_path.to_string(), Some(vec![]), Some(vec![]), None, true)?;
        ctx.register_table("vcf", Arc::new(table))?;

        let df = ctx
            .sql("SELECT `chrom`, `start`, `end`, `ref`, `alt` FROM vcf")
            .await?;
        let physical_plan = df.create_physical_plan().await?;
        let num_partitions = physical_plan.output_partitioning().partition_count();

        let lookup_exec = Arc::new(KvLookupExec::new(
            physical_plan,
            store.clone(),
            cache_columns.clone(),
            KvMatchMode::Exact,
            allele_matches,
            None,
            true,  // vcf_has_chr
            true,  // vcf_zero_based
            false, // cache_zero_based (VEP cache is 1-based)
            false, // extended_probes
        )?);

        let task_ctx = ctx.task_ctx();
        let start_time = Instant::now();

        let mut handles = Vec::new();
        for partition in 0..num_partitions {
            let plan = lookup_exec.clone();
            let ctx = task_ctx.clone();

            let handle = tokio::spawn(async move {
                let mut stream = plan.execute(partition, ctx)?;
                let mut total_rows = 0usize;

                while let Some(batch_result) = stream.next().await {
                    let batch = batch_result?;
                    total_rows += batch.num_rows();
                }

                Ok::<_, datafusion::common::DataFusionError>(total_rows)
            });
            handles.push(handle);
        }

        let mut total_rows = 0usize;
        for handle in handles {
            total_rows += handle.await.unwrap()?;
        }

        let elapsed = start_time.elapsed().as_secs_f64();
        println!(
            "{:>7} | {:>5} | {:>10} | {:>5.1}s | {:>12.0}",
            num_threads,
            num_partitions,
            total_rows,
            elapsed,
            total_rows as f64 / elapsed,
        );
    }

    Ok(())
}
