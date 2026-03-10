use datafusion::common::{DataFusionError, Result};
use datafusion::prelude::SessionContext;
use datafusion_bio_function_vep::register_vortex_cache;

#[tokio::main]
async fn main() -> Result<()> {
    let args: Vec<String> = std::env::args().collect();
    if args.len() < 2 || args.len() > 3 {
        return Err(DataFusionError::Execution(format!(
            "Usage: {} <vortex_path> [count]",
            args[0]
        )));
    }
    let path = &args[1];
    let do_count = matches!(args.get(2).map(String::as_str), Some("count"));
    let ctx = SessionContext::new();
    register_vortex_cache(&ctx, "vx", path).await?;
    if do_count {
        let batches = ctx
            .sql("SELECT COUNT(*) AS cnt FROM vx")
            .await?
            .collect()
            .await?;
        println!(
            "ok: count query succeeded path={} batches={}",
            path,
            batches.len()
        );
    } else {
        let batches = ctx.sql("SELECT * FROM vx LIMIT 1").await?.collect().await?;
        let rows: usize = batches.iter().map(|b| b.num_rows()).sum();
        let cols = batches.first().map(|b| b.num_columns()).unwrap_or(0);
        println!(
            "ok: readable vortex path={} rows_sampled={} columns={}",
            path, rows, cols
        );
    }
    Ok(())
}
