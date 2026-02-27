use std::path::{Path, PathBuf};
use std::sync::Arc;

use datafusion::common::{DataFusionError, Result};
use datafusion::execution::TaskContext;
use datafusion::logical_expr::Partitioning;
use datafusion::physical_plan::{ExecutionPlan, ExecutionPlanProperties};
use datafusion::prelude::{SessionConfig, SessionContext};
use futures::StreamExt;
use vortex::ArrayRef;
use vortex::VortexSessionDefault;
use vortex::arrow::FromArrowArray;
use vortex::compressor::CompactCompressor;
use vortex::dtype::DType;
use vortex::dtype::arrow::FromArrowType;
use vortex::file::{WriteOptionsSessionExt, WriteStrategyBuilder};
use vortex::layout::LayoutStrategy;
use vortex::session::VortexSession;

fn parse_arg<T: std::str::FromStr>(
    args: &[String],
    idx: usize,
    default: T,
    name: &str,
) -> Result<T> {
    args.get(idx)
        .map(|s| {
            s.parse()
                .map_err(|_| DataFusionError::Execution(format!("invalid {name}: '{s}'")))
        })
        .unwrap_or(Ok(default))
}

#[tokio::main]
async fn main() -> Result<()> {
    let args: Vec<String> = std::env::args().collect();
    if args.len() < 3 || args.len() > 9 {
        eprintln!(
            "Usage: {} <parquet_path> <vortex_output_dir> [target_partitions] [zstd_level] [values_per_page] [row_block_size] [max_parallel_writers] [limit_rows]",
            args[0]
        );
        std::process::exit(2);
    }

    let parquet_path = &args[1];
    let vortex_output = PathBuf::from(&args[2]);
    let target_partitions: usize = parse_arg(&args, 3, 8usize, "target_partitions")?;
    let zstd_level: i32 = parse_arg(&args, 4, 9i32, "zstd_level")?;
    let values_per_page: usize = parse_arg(&args, 5, 16_384usize, "values_per_page")?;
    let row_block_size: usize = parse_arg(&args, 6, 16_384usize, "row_block_size")?;
    let max_parallel_writers: usize =
        parse_arg(&args, 7, target_partitions.max(1), "max_parallel_writers")?;
    let limit_rows: Option<usize> = args
        .get(8)
        .map(|s| {
            s.parse().map_err(|_| {
                DataFusionError::Execution(format!("invalid limit_rows argument: '{s}'"))
            })
        })
        .transpose()?;

    if target_partitions == 0 {
        return Err(DataFusionError::Execution(
            "target_partitions must be >= 1".to_string(),
        ));
    }
    if values_per_page == 0 {
        return Err(DataFusionError::Execution(
            "values_per_page must be >= 1".to_string(),
        ));
    }
    if row_block_size == 0 {
        return Err(DataFusionError::Execution(
            "row_block_size must be >= 1".to_string(),
        ));
    }

    if vortex_output.exists() {
        let md = std::fs::metadata(&vortex_output)
            .map_err(|e| DataFusionError::Execution(format!("metadata failed: {e}")))?;
        if md.is_dir() {
            std::fs::remove_dir_all(&vortex_output).map_err(|e| {
                DataFusionError::Execution(format!("failed to remove existing dir: {e}"))
            })?;
        } else {
            std::fs::remove_file(&vortex_output).map_err(|e| {
                DataFusionError::Execution(format!("failed to remove existing file: {e}"))
            })?;
        }
    }
    std::fs::create_dir_all(&vortex_output)
        .map_err(|e| DataFusionError::Execution(format!("failed to create output dir: {e}")))?;

    let ctx = SessionContext::new_with_config(
        SessionConfig::new()
            .with_target_partitions(target_partitions)
            .with_repartition_file_scans(true)
            .with_repartition_file_min_size(0),
    );

    ctx.register_parquet("src", parquet_path, Default::default())
        .await?;

    let mut src_df = ctx.table("src").await?;
    if let Some(limit_rows) = limit_rows {
        src_df = src_df.limit(0, Some(limit_rows))?;
    }
    let src_df = src_df.repartition(Partitioning::RoundRobinBatch(target_partitions))?;

    let plan = src_df.create_physical_plan().await?;
    let partition_count = plan.output_partitioning().partition_count();
    let task_ctx = ctx.task_ctx();

    let session = VortexSession::default();
    let strategy = WriteStrategyBuilder::new()
        .with_compressor(
            CompactCompressor::default()
                .with_zstd_level(zstd_level)
                .with_values_per_page(values_per_page),
        )
        .with_row_block_size(row_block_size)
        .build();
    let dtype = DType::from_arrow(plan.schema());

    let max_parallel_writers = max_parallel_writers.clamp(1, partition_count.max(1));

    println!(
        "start: parquet={} output={} partitions={} zstd_level={} values_per_page={} row_block_size={} max_parallel_writers={} limit_rows={}",
        parquet_path,
        vortex_output.display(),
        partition_count,
        zstd_level,
        values_per_page,
        row_block_size,
        max_parallel_writers,
        limit_rows
            .map(|v| v.to_string())
            .unwrap_or_else(|| "none".to_string())
    );

    let mut total_rows = 0u64;
    let mut total_bytes = 0u64;
    let partition_futures = (0..partition_count).map(|partition_id| {
        write_partition(
            Arc::clone(&plan),
            Arc::clone(&task_ctx),
            session.clone(),
            strategy.clone(),
            dtype.clone(),
            vortex_output.clone(),
            partition_id,
        )
    });

    let mut results =
        futures::stream::iter(partition_futures).buffer_unordered(max_parallel_writers);
    while let Some(partition_result) = results.next().await {
        let (partition_id, rows, bytes) = partition_result?;
        total_rows += rows;
        total_bytes += bytes;
        if rows > 0 {
            println!("partition={partition_id} rows={rows} bytes={bytes}");
        }
    }

    let dir_size = dir_size_bytes(&vortex_output)?;
    println!(
        "ok: wrote vortex output at {} rows={} bytes={} dir_bytes={}",
        vortex_output.display(),
        total_rows,
        total_bytes,
        dir_size
    );
    Ok(())
}

async fn write_partition(
    plan: Arc<dyn ExecutionPlan>,
    task_ctx: Arc<TaskContext>,
    session: VortexSession,
    strategy: Arc<dyn LayoutStrategy>,
    dtype: DType,
    out_dir: PathBuf,
    partition_id: usize,
) -> Result<(usize, u64, u64)> {
    let mut stream = plan.execute(partition_id, task_ctx)?;
    let out_file = out_dir.join(format!("part_{partition_id:03}.vortex"));

    let mut writer = session.write_options().with_strategy(strategy).writer(
        tokio::fs::File::create(&out_file).await.map_err(|e| {
            DataFusionError::Execution(format!("failed to create {}: {e}", out_file.display()))
        })?,
        dtype,
    );

    let mut rows = 0u64;
    while let Some(batch_result) = stream.next().await {
        let batch = batch_result?;
        if batch.num_rows() == 0 {
            continue;
        }
        rows += batch.num_rows() as u64;
        writer
            .push(ArrayRef::from_arrow(batch, false))
            .await
            .map_err(|e| {
                DataFusionError::Execution(format!(
                    "failed to write partition {} chunk: {e}",
                    partition_id
                ))
            })?;
    }

    writer.finish().await.map_err(|e| {
        DataFusionError::Execution(format!(
            "failed to finish partition {} writer: {e}",
            partition_id
        ))
    })?;

    if rows == 0 {
        let _ = tokio::fs::remove_file(&out_file).await;
        return Ok((partition_id, rows, 0u64));
    }

    let bytes = tokio::fs::metadata(&out_file)
        .await
        .map_err(|e| {
            DataFusionError::Execution(format!("failed to stat {}: {e}", out_file.display()))
        })?
        .len();

    Ok((partition_id, rows, bytes))
}

fn dir_size_bytes(path: &Path) -> Result<u64> {
    let mut total = 0u64;
    for entry in std::fs::read_dir(path).map_err(|e| {
        DataFusionError::Execution(format!("failed to read {}: {e}", path.display()))
    })? {
        let entry = entry
            .map_err(|e| DataFusionError::Execution(format!("failed to read dir entry: {e}")))?;
        let md = entry.metadata().map_err(|e| {
            DataFusionError::Execution(format!(
                "failed to read metadata for {:?}: {e}",
                entry.path()
            ))
        })?;
        if md.is_file() {
            total = total.saturating_add(md.len());
        }
    }
    Ok(total)
}
