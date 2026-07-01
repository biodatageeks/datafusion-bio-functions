use std::path::PathBuf;

use anyhow::{Context, Result, bail};
use arrow_array::{Array, RecordBatch};
use clap::Parser;
use lance::Dataset;
use lance::dataset::ProjectionRequest;
use serde_json::Value;

const MAGIC: &[u8; 8] = b"LRSIDX02";

#[derive(Debug, Parser)]
struct Args {
    #[arg(long)]
    dataset_path: PathBuf,
    #[arg(long)]
    positions_file: PathBuf,
    #[arg(long)]
    sidecar_path: PathBuf,
    #[arg(long)]
    report_path: PathBuf,
    #[arg(long, default_value_t = 5000)]
    position_batch_size: usize,
}

#[derive(Debug)]
struct PositionRowIds {
    positions: Vec<u32>,
    row_ids: Vec<u64>,
}

#[derive(Debug, Default, Clone, Copy)]
struct BatchMemory {
    rows: usize,
    array_bytes: usize,
    buffer_bytes: usize,
}

#[tokio::main(flavor = "multi_thread")]
async fn main() -> Result<()> {
    let args = Args::parse();
    let dataset = Dataset::open(
        args.dataset_path
            .to_str()
            .context("dataset path is not valid UTF-8")?,
    )
    .await?;
    let projection = projection_from_report(&args.report_path)?;
    let projection_request =
        ProjectionRequest::from_columns(projection.iter().map(String::as_str), dataset.schema());
    let positions = read_positions(&args.positions_file)?;
    let sidecar = read_sidecar(&args.sidecar_path)?;

    let mut cursor = 0usize;
    let all_row_ids = sidecar.resolve_sorted_positions(&positions, &mut cursor)?;
    let all_batch = dataset
        .take_rows(&all_row_ids, projection_request.clone())
        .await?;
    let all_memory = batch_memory(&all_batch);

    drop(all_batch);

    let mut cursor = 0usize;
    let mut chunk_count = 0usize;
    let mut chunk_rows = 0usize;
    let mut retained_array_bytes = 0usize;
    let mut retained_buffer_bytes = 0usize;
    let mut max_chunk = BatchMemory::default();
    for chunk in positions.chunks(args.position_batch_size) {
        let row_ids = sidecar.resolve_sorted_positions(chunk, &mut cursor)?;
        let batch = dataset
            .take_rows(&row_ids, projection_request.clone())
            .await?;
        let memory = batch_memory(&batch);
        chunk_count += 1;
        chunk_rows += batch.num_rows();
        retained_array_bytes += memory.array_bytes;
        retained_buffer_bytes += memory.buffer_bytes;
        if memory.array_bytes > max_chunk.array_bytes {
            max_chunk = memory;
        }
    }

    println!("positions\t{}", positions.len());
    println!("resolved_row_ids\t{}", all_row_ids.len());
    println!("projection_columns\t{}", projection.len());
    println!("all_rows\t{}", all_memory.rows);
    println!("all_array_bytes\t{}", all_memory.array_bytes);
    println!("all_buffer_bytes\t{}", all_memory.buffer_bytes);
    println!("position_batch_size\t{}", args.position_batch_size);
    println!("chunk_count\t{}", chunk_count);
    println!("chunk_rows\t{}", chunk_rows);
    println!("max_chunk_rows\t{}", max_chunk.rows);
    println!("max_chunk_array_bytes\t{}", max_chunk.array_bytes);
    println!("max_chunk_buffer_bytes\t{}", max_chunk.buffer_bytes);
    println!("retained_chunks_array_bytes\t{}", retained_array_bytes);
    println!("retained_chunks_buffer_bytes\t{}", retained_buffer_bytes);
    Ok(())
}

fn projection_from_report(path: &PathBuf) -> Result<Vec<String>> {
    let value: Value = serde_json::from_slice(
        &std::fs::read(path).with_context(|| format!("failed to read '{}'", path.display()))?,
    )?;
    let fields = value
        .get("resolved_dataset_projection")
        .and_then(Value::as_array)
        .context("benchmark report missing resolved_dataset_projection")?;
    fields
        .iter()
        .map(|field| {
            field
                .as_str()
                .map(str::to_string)
                .context("projection field must be a string")
        })
        .collect()
}

fn read_positions(path: &PathBuf) -> Result<Vec<u64>> {
    let text = std::fs::read_to_string(path)
        .with_context(|| format!("failed to read '{}'", path.display()))?;
    let mut values = text
        .lines()
        .filter(|line| !line.trim().is_empty())
        .map(|line| line.trim().parse::<u64>().context("invalid position"))
        .collect::<Result<Vec<_>>>()?;
    values.sort_unstable();
    values.dedup();
    Ok(values)
}

fn read_sidecar(path: &PathBuf) -> Result<PositionRowIds> {
    let bytes =
        std::fs::read(path).with_context(|| format!("failed to read '{}'", path.display()))?;
    let mut offset = 0usize;
    if bytes.get(..MAGIC.len()) != Some(MAGIC.as_slice()) {
        bail!("invalid sidecar magic");
    }
    offset += MAGIC.len();
    let _unique_positions = read_u64(&bytes, &mut offset)? as usize;
    let row_id_count = read_u64(&bytes, &mut offset)? as usize;
    let positions = (0..row_id_count)
        .map(|_| read_u32(&bytes, &mut offset))
        .collect::<Result<Vec<_>>>()?;
    let row_ids = (0..row_id_count)
        .map(|_| read_u64(&bytes, &mut offset))
        .collect::<Result<Vec<_>>>()?;
    Ok(PositionRowIds { positions, row_ids })
}

impl PositionRowIds {
    fn resolve_sorted_positions(
        &self,
        query_positions: &[u64],
        cursor: &mut usize,
    ) -> Result<Vec<u64>> {
        let mut row_ids = Vec::new();
        for &query_position in query_positions {
            let Ok(position) = u32::try_from(query_position) else {
                continue;
            };
            while *cursor < self.positions.len() && self.positions[*cursor] < position {
                *cursor += 1;
            }
            if *cursor == self.positions.len() {
                break;
            }
            if self.positions[*cursor] == position {
                while *cursor < self.positions.len() && self.positions[*cursor] == position {
                    row_ids.push(self.row_ids[*cursor]);
                    *cursor += 1;
                }
            }
        }
        Ok(row_ids)
    }
}

fn batch_memory(batch: &RecordBatch) -> BatchMemory {
    BatchMemory {
        rows: batch.num_rows(),
        array_bytes: batch
            .columns()
            .iter()
            .map(|array| array.get_array_memory_size())
            .sum(),
        buffer_bytes: batch
            .columns()
            .iter()
            .map(|array| array.get_buffer_memory_size())
            .sum(),
    }
}

fn read_u32(bytes: &[u8], offset: &mut usize) -> Result<u32> {
    let end = *offset + 4;
    let value = u32::from_le_bytes(bytes[*offset..end].try_into()?);
    *offset = end;
    Ok(value)
}

fn read_u64(bytes: &[u8], offset: &mut usize) -> Result<u64> {
    let end = *offset + 8;
    let value = u64::from_le_bytes(bytes[*offset..end].try_into()?);
    *offset = end;
    Ok(value)
}
