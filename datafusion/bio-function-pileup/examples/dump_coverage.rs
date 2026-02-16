use std::sync::Arc;

use datafusion::arrow::array::{Int16Array, Int32Array, StringArray};
use datafusion::physical_plan::ExecutionPlan;
use datafusion::prelude::*;
use datafusion_bio_format_bam::table_provider::BamTableProvider;
use datafusion_bio_function_pileup::physical_exec::{PileupConfig, PileupExec};
use futures::StreamExt;

/// Dumps coverage output to stdout for comparison with mosdepth / samtools.
///
/// Pass `--zero-based` for 0-based output (default is 1-based).
/// Pass `--per-base` to emit per-position rows (contig, pos, coverage) instead of blocks.
#[tokio::main]
async fn main() {
    let args: Vec<String> = std::env::args().collect();
    if args.len() < 2 {
        eprintln!("Usage: dump_coverage <bam_file> [--zero-based] [--per-base]");
        std::process::exit(1);
    }
    let bam_path = &args[1];
    let zero_based = args.iter().any(|a| a == "--zero-based");
    let per_base = args.iter().any(|a| a == "--per-base");

    let table = BamTableProvider::new(bam_path.clone(), None, zero_based, None, true)
        .await
        .expect("Failed to open BAM file");

    let ctx = SessionContext::new();
    ctx.register_table("reads", Arc::new(table)).unwrap();

    let df = ctx.table("reads").await.unwrap();
    let plan = df.create_physical_plan().await.unwrap();
    let pileup = PileupExec::new(
        plan,
        PileupConfig {
            zero_based,
            per_base,
            ..PileupConfig::default()
        },
    );
    let task_ctx = ctx.task_ctx();

    let stream = pileup.execute(0, task_ctx).unwrap();
    let batches: Vec<_> = stream
        .collect::<Vec<_>>()
        .await
        .into_iter()
        .filter_map(|r| r.ok())
        .collect();

    if per_base {
        // Per-base schema: contig (Utf8), position (Int32), coverage (Int16)
        for batch in &batches {
            let contigs = batch
                .column(0)
                .as_any()
                .downcast_ref::<StringArray>()
                .unwrap();
            let positions = batch
                .column(1)
                .as_any()
                .downcast_ref::<Int32Array>()
                .unwrap();
            let covs = batch
                .column(2)
                .as_any()
                .downcast_ref::<Int16Array>()
                .unwrap();

            for i in 0..batch.num_rows() {
                println!(
                    "{}\t{}\t{}",
                    contigs.value(i),
                    positions.value(i),
                    covs.value(i)
                );
            }
        }
    } else {
        // Block schema: contig (Utf8), pos_start (Int32), pos_end (Int32), coverage (Int16)
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
}
