use datafusion::common::{DataFusionError, Result};
use datafusion::prelude::SessionContext;
use datafusion_bio_function_vep::register_vortex_cache;

#[tokio::main]
async fn main() -> Result<()> {
    let args: Vec<String> = std::env::args().collect();
    if args.len() != 2 {
        return Err(DataFusionError::Execution(format!(
            "Usage: {} <vortex_path>",
            args[0]
        )));
    }
    let path = &args[1];
    let ctx = SessionContext::new();
    register_vortex_cache(&ctx, "vx", path).await?;
    let batches = ctx.sql("SELECT * FROM vx LIMIT 1").await?.collect().await?;
    let rows: usize = batches.iter().map(|b| b.num_rows()).sum();
    let cols = batches.first().map(|b| b.num_columns()).unwrap_or(0);
    println!("ok: readable vortex path={} rows_sampled={} columns={}", path, rows, cols);
    Ok(())
}
