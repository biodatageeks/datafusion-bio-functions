use std::sync::Arc;
use std::time::Instant;

use datafusion::arrow::array::{Array, AsArray, RecordBatch};
use datafusion::prelude::*;
use datafusion_bio_format_bam::table_provider::BamTableProvider;
use datafusion_bio_function_pileup::cigar;
use futures::StreamExt;

/// Measures time spent in CIGAR string parsing vs total pipeline.
#[tokio::main]
async fn main() {
    let args: Vec<String> = std::env::args().collect();
    if args.len() < 2 {
        eprintln!("Usage: bench_cigar_parse <bam_file>");
        std::process::exit(1);
    }
    let bam_path = &args[1];

    let table = BamTableProvider::new(bam_path.clone(), None, true, None, false)
        .await
        .expect("Failed to open BAM file");

    let config = SessionConfig::new().with_target_partitions(1);
    let ctx = SessionContext::new_with_config(config);
    ctx.register_table("reads", Arc::new(table)).unwrap();

    let df = ctx.sql("SELECT start, cigar FROM reads").await.unwrap();
    let plan = df.create_physical_plan().await.unwrap();
    let task_ctx = ctx.task_ctx();

    let total_start = Instant::now();
    let stream = plan.execute(0, task_ctx).unwrap();
    let batches: Vec<RecordBatch> = stream
        .collect::<Vec<_>>()
        .await
        .into_iter()
        .filter_map(|r| r.ok())
        .collect();
    let io_time = total_start.elapsed();

    // Now measure CIGAR parsing time over the collected batches
    let parse_start = Instant::now();
    let mut total_reads = 0u64;
    let mut total_events = 0u64;
    for batch in &batches {
        let start_arr = batch
            .column(0)
            .as_primitive::<datafusion::arrow::datatypes::UInt32Type>();
        let cigar_arr = batch.column(1).as_string::<i32>();
        for row in 0..batch.num_rows() {
            if cigar_arr.is_null(row) {
                continue;
            }
            let cigar_str = cigar_arr.value(row);
            if cigar_str == "*" {
                continue;
            }
            let start = start_arr.value(row);
            let events = cigar::generate_events(start, cigar_str);
            total_events += events.len() as u64;
            total_reads += 1;
        }
    }
    let parse_time = parse_start.elapsed();

    let total_time = total_start.elapsed();

    println!("Reads processed:    {total_reads}");
    println!("Events generated:   {total_events}");
    println!("BAM I/O time:       {:.3}s", io_time.as_secs_f64());
    println!("CIGAR parse time:   {:.3}s", parse_time.as_secs_f64());
    println!("Total time:         {:.3}s", total_time.as_secs_f64());
    println!(
        "CIGAR parse fraction: {:.1}%",
        parse_time.as_secs_f64() / total_time.as_secs_f64() * 100.0
    );
}
