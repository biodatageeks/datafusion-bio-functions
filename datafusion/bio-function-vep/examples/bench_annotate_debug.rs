use std::sync::Arc;

use datafusion::arrow::array::Array;
use datafusion::common::Result;
use datafusion::physical_plan::ExecutionPlanProperties;
use datafusion::prelude::SessionContext;
use datafusion_bio_format_vcf::table_provider::VcfTableProvider;
use futures::StreamExt;

#[tokio::main]
async fn main() -> Result<()> {
    let vcf_path = std::env::args()
        .nth(1)
        .unwrap_or("/Users/mwiewior/research/git/polars-bio-vep-benchmark/vep-benchmark/data/HG002_chr22.vcf.gz".to_string());

    for partitions in [1, 4, 8] {
        let config = datafusion::prelude::SessionConfig::new().with_target_partitions(partitions);
        let ctx = SessionContext::new_with_config(config);

        let table =
            VcfTableProvider::new(vcf_path.clone(), Some(vec![]), Some(vec![]), None, true)?;
        ctx.register_table("vcf", Arc::new(table))?;

        let df = ctx.sql("SELECT `chrom`, `start` FROM vcf").await?;
        let plan = df.create_physical_plan().await?;
        let num_parts = plan.output_partitioning().partition_count();
        let task_ctx = ctx.task_ctx();

        println!("=== target_partitions={partitions}, actual={num_parts} ===");

        for p in 0..num_parts {
            let mut stream = plan.execute(p, task_ctx.clone())?;
            let mut rows = 0u64;
            let mut min_pos = u32::MAX;
            let mut max_pos = 0u32;

            while let Some(batch_result) = stream.next().await {
                let batch = batch_result?;
                rows += batch.num_rows() as u64;
                let starts = batch.column(1);
                if let Some(arr) = starts
                    .as_any()
                    .downcast_ref::<datafusion::arrow::array::UInt32Array>()
                {
                    for i in 0..arr.len() {
                        let v = arr.value(i);
                        min_pos = min_pos.min(v);
                        max_pos = max_pos.max(v);
                    }
                }
            }
            println!(
                "  partition {p}: {rows} rows, pos range [{min_pos}, {max_pos}], span={}",
                if max_pos >= min_pos {
                    max_pos - min_pos
                } else {
                    0
                }
            );
        }
        println!();
    }

    Ok(())
}
