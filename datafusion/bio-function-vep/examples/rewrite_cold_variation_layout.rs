//! Rewrite a cold variation parquet file with position-aligned row groups.
//!
//! Usage:
//!   cargo run --release --features kv-cache --example rewrite_cold_variation_layout -- \
//!     --input <chrN_cold.parquet> \
//!     --output-dir <dir-containing-chrN_cold.parquet> \
//!     --index-dir <dir-containing-chrN.posidx> \
//!     --chrom chrN \
//!     --row-group-rows 64000 \
//!     [--data-page-row-count 8192] \
//!     [--batch-size 262144] \
//!     [--bin-bits 20] \
//!     [--max-row-groups-per-file 30000]

use std::fs::File;
use std::path::{Path, PathBuf};
use std::time::Instant;

use datafusion::arrow::array::{Array, RecordBatch, UInt64Array};
use datafusion::arrow::compute::concat_batches;
use datafusion::arrow::datatypes::SchemaRef;
use datafusion::common::{DataFusionError, Result};
use datafusion_bio_function_vep::kv_cache::position_index::{
    PositionIndex, cold_variation_files_for_chrom, position_index_file,
};
use datafusion_bio_function_vep::warm_cache::cold_parquet::{
    ColdParquetLookupSet, cold_parquet_projection_columns,
};
use parquet::arrow::ArrowWriter;
use parquet::arrow::arrow_reader::ParquetRecordBatchReaderBuilder;
use parquet::basic::Compression;
use parquet::file::metadata::KeyValue;
use parquet::file::properties::WriterProperties;

#[derive(Debug)]
struct Args {
    input: PathBuf,
    output_dir: PathBuf,
    index_dir: PathBuf,
    chrom: String,
    row_group_rows: usize,
    data_page_row_count: Option<usize>,
    batch_size: usize,
    bin_bits: Option<u8>,
    max_row_groups_per_file: Option<usize>,
}

#[derive(Debug)]
struct RewriteStats {
    rows: usize,
    row_groups: usize,
    index_positions: usize,
    files: usize,
}

fn parse_args() -> Result<Args> {
    parse_args_from(std::env::args())
}

fn parse_args_from<I, S>(args: I) -> Result<Args>
where
    I: IntoIterator<Item = S>,
    S: Into<String>,
{
    let args: Vec<String> = args.into_iter().map(Into::into).collect();
    let mut input = None;
    let mut output_dir = None;
    let mut index_dir = None;
    let mut chrom = None;
    let mut row_group_rows = None;
    let mut data_page_row_count = None;
    let mut batch_size = 262_144usize;
    let mut bin_bits = None;
    let mut max_row_groups_per_file = None;

    let mut i = 1;
    while i < args.len() {
        match args[i].as_str() {
            "--input" => {
                i += 1;
                input = Some(PathBuf::from(require_arg(&args, i, "--input")?));
            }
            "--output-dir" => {
                i += 1;
                output_dir = Some(PathBuf::from(require_arg(&args, i, "--output-dir")?));
            }
            "--index-dir" => {
                i += 1;
                index_dir = Some(PathBuf::from(require_arg(&args, i, "--index-dir")?));
            }
            "--chrom" => {
                i += 1;
                chrom = Some(require_arg(&args, i, "--chrom")?.to_string());
            }
            "--row-group-rows" => {
                i += 1;
                row_group_rows = Some(
                    require_arg(&args, i, "--row-group-rows")?
                        .parse::<usize>()
                        .map_err(|error| {
                            DataFusionError::Execution(format!(
                                "invalid --row-group-rows value: {error}"
                            ))
                        })?,
                );
            }
            "--data-page-row-count" => {
                i += 1;
                data_page_row_count = Some(
                    require_arg(&args, i, "--data-page-row-count")?
                        .parse::<usize>()
                        .map_err(|error| {
                            DataFusionError::Execution(format!(
                                "invalid --data-page-row-count value: {error}"
                            ))
                        })?
                        .max(1),
                );
            }
            "--batch-size" => {
                i += 1;
                batch_size = require_arg(&args, i, "--batch-size")?
                    .parse::<usize>()
                    .map_err(|error| {
                        DataFusionError::Execution(format!("invalid --batch-size value: {error}"))
                    })?;
            }
            "--bin-bits" => {
                i += 1;
                bin_bits = Some(require_arg(&args, i, "--bin-bits")?.parse::<u8>().map_err(
                    |error| {
                        DataFusionError::Execution(format!("invalid --bin-bits value: {error}"))
                    },
                )?);
            }
            "--max-row-groups-per-file" => {
                i += 1;
                max_row_groups_per_file = Some(
                    require_arg(&args, i, "--max-row-groups-per-file")?
                        .parse::<usize>()
                        .map_err(|error| {
                            DataFusionError::Execution(format!(
                                "invalid --max-row-groups-per-file value: {error}"
                            ))
                        })?
                        .max(1),
                );
            }
            other => {
                return Err(DataFusionError::Execution(format!(
                    "unknown argument: {other}"
                )));
            }
        }
        i += 1;
    }

    Ok(Args {
        input: input.ok_or_else(|| DataFusionError::Execution("--input is required".into()))?,
        output_dir: output_dir
            .ok_or_else(|| DataFusionError::Execution("--output-dir is required".into()))?,
        index_dir: index_dir
            .ok_or_else(|| DataFusionError::Execution("--index-dir is required".into()))?,
        chrom: chrom.ok_or_else(|| DataFusionError::Execution("--chrom is required".into()))?,
        row_group_rows: row_group_rows
            .ok_or_else(|| DataFusionError::Execution("--row-group-rows is required".into()))?
            .max(1),
        data_page_row_count,
        batch_size: batch_size.max(1),
        bin_bits,
        max_row_groups_per_file,
    })
}

fn require_arg<'a>(args: &'a [String], index: usize, flag: &str) -> Result<&'a str> {
    args.get(index)
        .map(String::as_str)
        .ok_or_else(|| DataFusionError::Execution(format!("{flag} requires a value")))
}

fn main() -> Result<()> {
    let args = parse_args()?;
    std::fs::create_dir_all(&args.output_dir).map_err(|error| {
        DataFusionError::Execution(format!(
            "failed to create output dir '{}': {error}",
            args.output_dir.display()
        ))
    })?;
    std::fs::create_dir_all(&args.index_dir).map_err(|error| {
        DataFusionError::Execution(format!(
            "failed to create index dir '{}': {error}",
            args.index_dir.display()
        ))
    })?;

    let output = output_path_for_part(&args, 0);
    eprintln!("input={}", args.input.display());
    eprintln!("output={}", output.display());
    eprintln!("index_dir={}", args.index_dir.display());
    eprintln!("row_group_rows={}", args.row_group_rows);
    eprintln!(
        "data_page_row_count={}",
        args.data_page_row_count
            .map(|value| value.to_string())
            .unwrap_or_else(|| "default".to_string())
    );
    eprintln!("batch_size={}", args.batch_size);
    eprintln!(
        "bin_bits={}",
        args.bin_bits
            .map(|bits| bits.to_string())
            .unwrap_or_else(|| "none".to_string())
    );
    eprintln!(
        "max_row_groups_per_file={}",
        args.max_row_groups_per_file
            .map(|value| value.to_string())
            .unwrap_or_else(|| "none".to_string())
    );

    let started = Instant::now();
    let stats = rewrite(&args)?;
    eprintln!(
        "rewritten rows={} row_groups={} files={} index_positions={} elapsed={:.2}s",
        stats.rows,
        stats.row_groups,
        stats.files,
        stats.index_positions,
        started.elapsed().as_secs_f64()
    );
    Ok(())
}

fn rewrite(args: &Args) -> Result<RewriteStats> {
    let input = File::open(&args.input).map_err(|error| {
        DataFusionError::Execution(format!(
            "failed to open input '{}': {error}",
            args.input.display()
        ))
    })?;
    let builder = ParquetRecordBatchReaderBuilder::try_new(input)?;
    let schema = builder.schema().clone();
    schema.index_of("position_key")?;

    let mut writer = PositionAlignedWriter::new(args, schema, args.row_group_rows, args.bin_bits)?;
    let reader = builder.with_batch_size(args.batch_size).build()?;
    let mut rows = 0usize;

    for batch in reader {
        let batch = batch?;
        rows += batch.num_rows();
        writer.write(&batch)?;
    }

    let row_groups = writer.close()?;
    let output_files = cold_variation_files_for_chrom(&args.output_dir, &args.chrom);
    ColdParquetLookupSet::open(
        output_files.iter(),
        cold_parquet_projection_columns(&[], false),
        args.batch_size,
        1,
    )?;

    let index = PositionIndex::from_parquet_files(output_files.iter(), args.batch_size)?;
    let index_positions = index.len();
    let index_path = position_index_file(&args.index_dir, &args.chrom);
    index.write_to_path(index_path)?;

    Ok(RewriteStats {
        rows,
        row_groups,
        index_positions,
        files: output_files.len(),
    })
}

fn output_path_for_part(args: &Args, part: usize) -> PathBuf {
    if args.max_row_groups_per_file.is_some() {
        args.output_dir
            .join(format!("{}_cold.part{part:04}.parquet", args.chrom))
    } else {
        args.output_dir.join(format!("{}_cold.parquet", args.chrom))
    }
}

fn create_writer(path: &Path, schema: SchemaRef, args: &Args) -> Result<ArrowWriter<File>> {
    let mut metadata = vec![
        KeyValue::new("vepyr.cache_tier".to_string(), "cold".to_string()),
        KeyValue::new(
            "vepyr.row_group_rows".to_string(),
            args.row_group_rows.to_string(),
        ),
    ];
    if let Some(bits) = args.bin_bits {
        metadata.push(KeyValue::new(
            "vepyr.echtvar_bin_bits".to_string(),
            bits.to_string(),
        ));
    }
    if let Some(rows) = args.data_page_row_count {
        metadata.push(KeyValue::new(
            "vepyr.data_page_row_count".to_string(),
            rows.to_string(),
        ));
    }
    let mut props = WriterProperties::builder()
        .set_compression(Compression::ZSTD(Default::default()))
        .set_max_row_group_size(usize::MAX)
        .set_key_value_metadata(Some(metadata));
    if let Some(rows) = args.data_page_row_count {
        props = props
            .set_data_page_row_count_limit(rows)
            .set_write_batch_size(rows);
    }
    let props = props.build();
    let file = File::create(path).map_err(|error| {
        DataFusionError::Execution(format!(
            "failed to create output '{}': {error}",
            path.display()
        ))
    })?;
    Ok(ArrowWriter::try_new(file, schema, Some(props))?)
}

struct PositionAlignedWriter {
    args: Args,
    writer: Option<ArrowWriter<File>>,
    schema: SchemaRef,
    position_idx: usize,
    target_rows: usize,
    bin_bits: Option<u8>,
    pending: Vec<RecordBatch>,
    pending_rows: usize,
    row_groups: usize,
    current_file_row_groups: usize,
    current_part: usize,
}

impl PositionAlignedWriter {
    fn new(
        args: &Args,
        schema: SchemaRef,
        target_rows: usize,
        bin_bits: Option<u8>,
    ) -> Result<Self> {
        let position_idx = schema.index_of("position_key")?;
        let output = output_path_for_part(args, 0);
        let writer = create_writer(&output, schema.clone(), args)?;
        Ok(Self {
            args: Args {
                input: args.input.clone(),
                output_dir: args.output_dir.clone(),
                index_dir: args.index_dir.clone(),
                chrom: args.chrom.clone(),
                row_group_rows: args.row_group_rows,
                data_page_row_count: args.data_page_row_count,
                batch_size: args.batch_size,
                bin_bits: args.bin_bits,
                max_row_groups_per_file: args.max_row_groups_per_file,
            },
            writer: Some(writer),
            schema,
            position_idx,
            target_rows: target_rows.max(1),
            bin_bits,
            pending: Vec::new(),
            pending_rows: 0,
            row_groups: 0,
            current_file_row_groups: 0,
            current_part: 0,
        })
    }

    fn write(&mut self, batch: &RecordBatch) -> Result<()> {
        if batch.num_rows() == 0 {
            return Ok(());
        }
        self.pending_rows += batch.num_rows();
        self.pending.push(batch.clone());
        self.flush_ready(false)
    }

    fn close(mut self) -> Result<usize> {
        self.flush_ready(true)?;
        if let Some(writer) = self.writer.take() {
            writer.close()?;
        }
        Ok(self.row_groups)
    }

    fn flush_ready(&mut self, final_flush: bool) -> Result<()> {
        loop {
            let cutoff = if final_flush {
                self.pending_rows
            } else {
                match self.next_flush_cutoff()? {
                    Some(cutoff) => cutoff,
                    None => return Ok(()),
                }
            };
            if cutoff == 0 {
                return Ok(());
            }
            self.flush_prefix(cutoff)?;
            if final_flush || self.pending_rows < self.target_rows {
                return Ok(());
            }
        }
    }

    fn next_flush_cutoff(&self) -> Result<Option<usize>> {
        if self.pending_rows < self.target_rows && self.bin_bits.is_none() {
            return Ok(None);
        }

        let mut previous_position = None;
        let mut previous_bin = None;
        let mut global_row = 0usize;
        for batch in &self.pending {
            let positions = position_key_array(batch, self.position_idx)?;
            for row in 0..batch.num_rows() {
                if positions.is_null(row) {
                    return Err(DataFusionError::Execution(
                        "position_key must be non-null UInt64".into(),
                    ));
                }
                let current = positions.value(row);
                let current_bin = self.bin_bits.map(|bits| genomic_bin(current, bits));
                if global_row > 0 && previous_bin.zip(current_bin).is_some_and(|(a, b)| a != b) {
                    return Ok(Some(global_row));
                }
                if global_row >= self.target_rows
                    && previous_position.is_some_and(|prev| prev != current)
                {
                    return Ok(Some(global_row));
                }
                previous_position = Some(current);
                previous_bin = current_bin;
                global_row += 1;
            }
        }

        Ok(None)
    }

    fn flush_prefix(&mut self, cutoff: usize) -> Result<()> {
        let mut remaining = cutoff;
        let mut flush_batches = Vec::new();
        let mut keep_batches = Vec::new();

        for batch in self.pending.drain(..) {
            let rows = batch.num_rows();
            if remaining == 0 {
                keep_batches.push(batch);
            } else if rows <= remaining {
                remaining -= rows;
                flush_batches.push(batch);
            } else {
                flush_batches.push(batch.slice(0, remaining));
                keep_batches.push(batch.slice(remaining, rows - remaining));
                remaining = 0;
            }
        }

        if remaining != 0 {
            return Err(DataFusionError::Execution(
                "internal error while splitting position-aligned row group".into(),
            ));
        }

        let batch = concat_batches(&self.schema, flush_batches.iter())
            .map_err(|error| DataFusionError::ArrowError(Box::new(error), None))?;
        self.rotate_if_needed()?;
        let writer = self
            .writer
            .as_mut()
            .ok_or_else(|| DataFusionError::Execution("cold parquet writer closed".into()))?;
        writer.write(&batch)?;
        writer.flush()?;
        self.row_groups += 1;
        self.current_file_row_groups += 1;
        self.pending = keep_batches;
        self.pending_rows -= cutoff;
        Ok(())
    }

    fn rotate_if_needed(&mut self) -> Result<()> {
        let Some(max_row_groups) = self.args.max_row_groups_per_file else {
            return Ok(());
        };
        if self.current_file_row_groups < max_row_groups {
            return Ok(());
        }

        if let Some(writer) = self.writer.take() {
            writer.close()?;
        }
        self.current_part += 1;
        self.current_file_row_groups = 0;
        let output = output_path_for_part(&self.args, self.current_part);
        self.writer = Some(create_writer(&output, self.schema.clone(), &self.args)?);
        Ok(())
    }
}

fn position_key_array(batch: &RecordBatch, idx: usize) -> Result<&UInt64Array> {
    batch
        .column(idx)
        .as_any()
        .downcast_ref::<UInt64Array>()
        .ok_or_else(|| DataFusionError::Execution("position_key must be UInt64".into()))
}

fn genomic_bin(position_key: u64, bits: u8) -> u64 {
    const POSITION_MASK: u64 = (1_u64 << 48) - 1;
    (position_key & POSITION_MASK) >> bits
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn data_page_row_count_arg_is_parsed() {
        let args = parse_args_from([
            "rewrite_cold_variation_layout",
            "--input",
            "chr4_cold.parquet",
            "--output-dir",
            "variation",
            "--index-dir",
            "variation.position_index",
            "--chrom",
            "chr4",
            "--row-group-rows",
            "8192",
            "--data-page-row-count",
            "2048",
        ])
        .unwrap();

        assert_eq!(args.data_page_row_count, Some(2048));
    }
}
