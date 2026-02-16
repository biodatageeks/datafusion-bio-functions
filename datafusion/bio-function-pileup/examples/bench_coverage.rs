use std::sync::Arc;
use std::time::Instant;

use datafusion::physical_plan::ExecutionPlan;
use datafusion::prelude::*;
use datafusion_bio_format_bam::table_provider::BamTableProvider;
use datafusion_bio_function_pileup::physical_exec::{DenseMode, PileupConfig, PileupExec};
use futures::StreamExt;

#[tokio::main]
async fn main() {
    let args: Vec<String> = std::env::args().collect();
    if args.len() < 2 {
        eprintln!(
            "Usage: bench_coverage <bam_file> [--partitions N] [--dense] [--binary-cigar] [--zero-based]"
        );
        std::process::exit(1);
    }
    let bam_path = &args[1];
    let binary_cigar = args.iter().any(|a| a == "--binary-cigar");
    let zero_based = args.iter().any(|a| a == "--zero-based");

    println!("BAM file: {bam_path}");
    println!("---");

    let start = Instant::now();

    let table = BamTableProvider::new(bam_path.clone(), None, zero_based, None, binary_cigar)
        .await
        .expect("Failed to open BAM file");
    let open_time = start.elapsed();
    println!("BAM open time: {:.3}s", open_time.as_secs_f64());

    // Parse optional --partitions flag (default: use DataFusion's target_partitions)
    let target_partitions: Option<usize> = args
        .iter()
        .position(|a| a == "--partitions")
        .and_then(|i| args.get(i + 1))
        .and_then(|v| v.parse().ok());

    let config = if let Some(tp) = target_partitions {
        SessionConfig::new().with_target_partitions(tp)
    } else {
        SessionConfig::new()
    };
    let ctx = SessionContext::new_with_config(config);
    ctx.register_table("reads", Arc::new(table)).unwrap();

    // Only select columns needed for coverage (avoids decoding sequence, quality, tags)
    let df = ctx
        .sql("SELECT chrom, start, flags, cigar, mapping_quality FROM reads")
        .await
        .unwrap();
    let plan = df.create_physical_plan().await.unwrap();

    let num_partitions = plan.properties().partitioning.partition_count();
    println!("Input partitions: {num_partitions}");

    let dense_mode = if args.iter().any(|a| a == "--dense") {
        DenseMode::Force
    } else {
        DenseMode::default()
    };
    println!("Dense mode: {dense_mode:?}");
    println!("Binary CIGAR: {binary_cigar}");
    println!("Zero-based: {zero_based}");

    let pileup_config = PileupConfig {
        dense_mode,
        binary_cigar,
        zero_based,
        ..PileupConfig::default()
    };
    let pileup = Arc::new(PileupExec::new(plan, pileup_config));
    let task_ctx = ctx.task_ctx();

    let compute_start = Instant::now();

    // Run all partitions concurrently
    let mut handles = Vec::new();
    for partition in 0..num_partitions {
        let pileup = pileup.clone();
        let task_ctx = task_ctx.clone();
        handles.push(tokio::spawn(async move {
            let stream = pileup.execute(partition, task_ctx).unwrap();
            let batches: Vec<_> = stream
                .collect::<Vec<_>>()
                .await
                .into_iter()
                .filter_map(|r| r.ok())
                .collect();

            let mut blocks = 0u64;
            for batch in &batches {
                blocks += batch.num_rows() as u64;
            }
            blocks
        }));
    }

    let mut total_blocks = 0u64;
    for handle in handles {
        total_blocks += handle.await.unwrap();
    }

    let compute_time = compute_start.elapsed();
    let total_time = start.elapsed();

    println!("---");
    println!("Coverage blocks emitted: {total_blocks}");
    println!("Compute time:  {:.3}s", compute_time.as_secs_f64());
    println!("Total time:    {:.3}s", total_time.as_secs_f64());
}
