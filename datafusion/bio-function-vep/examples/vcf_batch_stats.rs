//! Reads a VCF file via the datafusion-bio-formats VCF reader and prints
//! per-batch statistics: number of rows, min/max position, and genomic span.

use std::sync::Arc;

use datafusion::arrow::array::{Array, UInt32Array};
use datafusion::datasource::TableProvider;
use datafusion::prelude::SessionContext;
use datafusion_bio_format_vcf::table_provider::VcfTableProvider;

#[tokio::main]
async fn main() -> datafusion::common::Result<()> {
    let args: Vec<String> = std::env::args().collect();
    let vcf_path = if args.len() > 1 {
        args[1].clone()
    } else {
        "/Users/mwiewior/research/git/polars-bio-vep-benchmark/vep-benchmark/data/HG002_chr22.vcf.gz"
            .to_string()
    };

    println!("VCF file: {vcf_path}");

    // Create the VCF table provider.
    // info_fields=None means include all INFO fields; format_fields=None means all FORMAT fields.
    // Using None/None for simplicity but we only need core columns (chrom, start, end).
    // coordinate_system_zero_based=true => 0-based half-open output.
    let table = VcfTableProvider::new(
        vcf_path.clone(),
        Some(vec![]), // no INFO fields (faster)
        Some(vec![]), // no FORMAT fields (faster)
        None,         // no object storage options
        true,         // 0-based half-open coordinates
    )?;

    println!("Schema: {:?}", table.schema());
    println!("---");

    let ctx = SessionContext::new();
    ctx.register_table("vcf", Arc::new(table))?;

    // Collect all batches via SQL to get start/end per batch.
    let df = ctx.sql("SELECT `chrom`, `start`, `end` FROM vcf").await?;
    let batches = df.collect().await?;

    let num_batches = batches.len();
    let mut total_rows: usize = 0;

    println!(
        "{:>5}  {:>8}  {:>12}  {:>12}  {:>12}  {:>12}",
        "batch", "rows", "min_start", "max_start", "min_end", "span (end-start)"
    );
    println!("{}", "-".repeat(75));

    for (i, batch) in batches.iter().enumerate() {
        let rows = batch.num_rows();
        total_rows += rows;

        let start_col = batch.column_by_name("start").expect("missing start column");
        let end_col = batch.column_by_name("end").expect("missing end column");

        // start and end are UInt32 in the VCF reader
        let starts = start_col
            .as_any()
            .downcast_ref::<UInt32Array>()
            .expect("start column is not UInt32");
        let ends = end_col
            .as_any()
            .downcast_ref::<UInt32Array>()
            .expect("end column is not UInt32");

        let min_start = (0..starts.len())
            .filter(|&j| !starts.is_null(j))
            .map(|j| starts.value(j))
            .min()
            .unwrap_or(0);
        let max_start = (0..starts.len())
            .filter(|&j| !starts.is_null(j))
            .map(|j| starts.value(j))
            .max()
            .unwrap_or(0);
        let min_end = (0..ends.len())
            .filter(|&j| !ends.is_null(j))
            .map(|j| ends.value(j))
            .min()
            .unwrap_or(0);
        let max_end = (0..ends.len())
            .filter(|&j| !ends.is_null(j))
            .map(|j| ends.value(j))
            .max()
            .unwrap_or(0);

        let span = max_end.saturating_sub(min_start);

        println!(
            "{:>5}  {:>8}  {:>12}  {:>12}  {:>12}  {:>12}",
            i, rows, min_start, max_start, min_end, span
        );
    }

    println!("{}", "-".repeat(75));
    println!("Total: {} batches, {} rows", num_batches, total_rows);

    Ok(())
}
