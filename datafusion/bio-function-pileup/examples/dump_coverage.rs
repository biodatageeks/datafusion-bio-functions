use std::sync::Arc;

use datafusion::arrow::array::{Int16Array, Int32Array, StringArray};
use datafusion::physical_plan::ExecutionPlan;
use datafusion::prelude::*;
use datafusion_bio_format_bam::table_provider::BamTableProvider;
use datafusion_bio_function_pileup::physical_exec::{PileupConfig, PileupExec};
use futures::StreamExt;

/// Dumps coverage output as BED (0-based half-open) to stdout for comparison with mosdepth.
/// Our output is 0-based closed, so we convert pos_end from inclusive to exclusive (pos_end + 1).
#[tokio::main]
async fn main() {
    let args: Vec<String> = std::env::args().collect();
    if args.len() < 2 {
        eprintln!("Usage: dump_coverage <bam_file>");
        std::process::exit(1);
    }
    let bam_path = &args[1];

    let table = BamTableProvider::new(bam_path.clone(), None, true, None, true)
        .await
        .expect("Failed to open BAM file");

    // Use 1 partition for deterministic sorted output
    let config = SessionConfig::new().with_target_partitions(1);
    let ctx = SessionContext::new_with_config(config);
    ctx.register_table("reads", Arc::new(table)).unwrap();

    let df = ctx.table("reads").await.unwrap();
    let plan = df.create_physical_plan().await.unwrap();
    let pileup = PileupExec::new(plan, PileupConfig::default());
    let task_ctx = ctx.task_ctx();

    let stream = pileup.execute(0, task_ctx).unwrap();
    let batches: Vec<_> = stream
        .collect::<Vec<_>>()
        .await
        .into_iter()
        .filter_map(|r| r.ok())
        .collect();

    for batch in &batches {
        let contigs = batch
            .column(0)
            .as_any()
            .downcast_ref::<StringArray>()
            .unwrap();
        let starts = batch
            .column(1)
            .as_any()
            .downcast_ref::<Int32Array>()
            .unwrap();
        let ends = batch
            .column(2)
            .as_any()
            .downcast_ref::<Int32Array>()
            .unwrap();
        let covs = batch
            .column(3)
            .as_any()
            .downcast_ref::<Int16Array>()
            .unwrap();

        for i in 0..batch.num_rows() {
            // Convert from 0-based closed to 0-based half-open (BED format)
            println!(
                "{}\t{}\t{}\t{}",
                contigs.value(i),
                starts.value(i),
                ends.value(i) + 1,
                covs.value(i)
            );
        }
    }
}
