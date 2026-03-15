use datafusion::dataframe::DataFrameWriteOptions;
use datafusion::prelude::SessionContext;

#[tokio::main]
async fn main() -> datafusion::common::Result<()> {
    let args: Vec<String> = std::env::args().collect();
    if args.len() < 3 {
        eprintln!(
            "Usage: {} <input_parquet_path> <output_partitioned_dir>",
            args[0]
        );
        std::process::exit(1);
    }

    let input_path = &args[1];
    let output_dir = &args[2];

    let ctx = SessionContext::new();
    let df = ctx.read_parquet(input_path, Default::default()).await?;

    let output_path = std::path::Path::new(output_dir);
    if output_path.exists() {
        std::fs::remove_dir_all(output_path)?;
    }

    df.write_parquet(
        output_dir,
        DataFrameWriteOptions::new().with_partition_by(vec!["chrom".to_string()]),
        None,
    )
    .await?;

    println!("Rewritten parquet to: {output_dir}");
    Ok(())
}
