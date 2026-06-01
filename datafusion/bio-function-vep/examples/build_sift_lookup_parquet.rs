//! Build compact transcript-id keyed SIFT/PolyPhen lookup parquet.
//!
//! Usage:
//!   cargo run --release -p datafusion-bio-function-vep --example build_sift_lookup_parquet -- \
//!     <translation_sift.parquet-or-dir> <output-dir> [row_group_size] [zstd|lz4_raw|snappy|none]
//!
//! Defaults are tuned for lookup speed: 16 rows per group and no compression.

use std::path::{Path, PathBuf};
use std::time::Instant;

use datafusion::arrow::array::{
    Array, ArrayRef, BinaryArray, Float32Array, Int32Array, LargeStringArray, ListArray,
    RecordBatch, StringArray, StringViewArray,
};
use datafusion::arrow::datatypes::{DataType, Field, Schema};
use datafusion::common::{DataFusionError, Result};
use datafusion_bio_function_vep::transcript_consequence::{CachedPredictions, CompactPrediction};
use parquet::arrow::ArrowWriter;
use parquet::arrow::ProjectionMask;
use parquet::arrow::arrow_reader::{
    ArrowReaderMetadata, ArrowReaderOptions, ParquetRecordBatchReaderBuilder,
};
use parquet::basic::Compression;
use parquet::file::properties::WriterProperties;

const DEFAULT_ROW_GROUP_SIZE: usize = 16;

fn main() -> Result<()> {
    let mut args = std::env::args().skip(1);
    let input = args.next().map(PathBuf::from).ok_or_else(usage)?;
    let output_dir = args.next().map(PathBuf::from).ok_or_else(usage)?;
    let row_group_size = args
        .next()
        .and_then(|value| value.parse::<usize>().ok())
        .unwrap_or(DEFAULT_ROW_GROUP_SIZE)
        .max(1);
    let compression = args
        .next()
        .as_deref()
        .map(parse_compression)
        .transpose()?
        .unwrap_or(Compression::UNCOMPRESSED);

    std::fs::create_dir_all(&output_dir).map_err(|error| {
        DataFusionError::Execution(format!(
            "failed to create output directory '{}': {error}",
            output_dir.display()
        ))
    })?;

    let input_files = input_files(&input)?;
    if input_files.is_empty() {
        return Err(DataFusionError::Execution(format!(
            "no parquet files found in '{}'",
            input.display()
        )));
    }

    let started = Instant::now();
    let mut total = BuildStats::default();
    for input_file in input_files {
        let output_file = output_dir.join(input_file.file_name().ok_or_else(|| {
            DataFusionError::Execution(format!(
                "input path '{}' has no file name",
                input_file.display()
            ))
        })?);
        let stats = convert_file(&input_file, &output_file, row_group_size, compression)?;
        eprintln!(
            "{} -> {} transcripts={} sift_entries={} polyphen_entries={} bytes={} elapsed_s={:.3}",
            input_file.display(),
            output_file.display(),
            stats.transcripts,
            stats.sift_entries,
            stats.polyphen_entries,
            stats.prediction_bytes,
            stats.elapsed_s,
        );
        total += stats;
    }
    eprintln!(
        "done files={} transcripts={} sift_entries={} polyphen_entries={} bytes={} elapsed_s={:.3}",
        total.files,
        total.transcripts,
        total.sift_entries,
        total.polyphen_entries,
        total.prediction_bytes,
        started.elapsed().as_secs_f64(),
    );
    Ok(())
}

fn usage() -> DataFusionError {
    DataFusionError::Execution(
        "usage: build_sift_lookup_parquet <translation_sift.parquet-or-dir> <output-dir> [row_group_size=16] [zstd|lz4_raw|snappy|none; default=none]"
            .to_string(),
    )
}

fn parse_compression(value: &str) -> Result<Compression> {
    match value.to_ascii_lowercase().as_str() {
        "zstd" => Ok(Compression::ZSTD(Default::default())),
        "lz4" | "lz4_raw" => Ok(Compression::LZ4_RAW),
        "snappy" => Ok(Compression::SNAPPY),
        "none" | "uncompressed" => Ok(Compression::UNCOMPRESSED),
        other => Err(DataFusionError::Execution(format!(
            "unsupported compression '{other}'; expected zstd, lz4_raw, snappy, or none"
        ))),
    }
}

#[derive(Debug, Default, Clone, Copy)]
struct BuildStats {
    files: u64,
    transcripts: u64,
    sift_entries: u64,
    polyphen_entries: u64,
    prediction_bytes: u64,
    elapsed_s: f64,
}

impl std::ops::AddAssign for BuildStats {
    fn add_assign(&mut self, rhs: Self) {
        self.files += rhs.files;
        self.transcripts += rhs.transcripts;
        self.sift_entries += rhs.sift_entries;
        self.polyphen_entries += rhs.polyphen_entries;
        self.prediction_bytes += rhs.prediction_bytes;
        self.elapsed_s += rhs.elapsed_s;
    }
}

fn input_files(path: &Path) -> Result<Vec<PathBuf>> {
    if path.is_file() {
        return Ok(vec![path.to_path_buf()]);
    }
    let mut files = std::fs::read_dir(path)
        .map_err(|error| {
            DataFusionError::Execution(format!(
                "failed to list input directory '{}': {error}",
                path.display()
            ))
        })?
        .filter_map(|entry| entry.ok().map(|entry| entry.path()))
        .filter(|path| {
            path.extension()
                .and_then(|ext| ext.to_str())
                .is_some_and(|ext| ext.eq_ignore_ascii_case("parquet"))
        })
        .collect::<Vec<_>>();
    files.sort();
    Ok(files)
}

fn convert_file(
    input: &Path,
    output: &Path,
    row_group_size: usize,
    compression: Compression,
) -> Result<BuildStats> {
    let started = Instant::now();
    let input_file = std::fs::File::open(input).map_err(|error| {
        DataFusionError::Execution(format!("failed to open '{}': {error}", input.display()))
    })?;
    let metadata =
        ArrowReaderMetadata::load(&input_file, ArrowReaderOptions::default()).map_err(|error| {
            DataFusionError::Execution(format!(
                "failed to load parquet metadata '{}': {error}",
                input.display()
            ))
        })?;
    let projection = projection(
        &metadata,
        &["transcript_id", "sift_predictions", "polyphen_predictions"],
    )?;
    let input_file = std::fs::File::open(input).map_err(|error| {
        DataFusionError::Execution(format!("failed to open '{}': {error}", input.display()))
    })?;
    let reader = ParquetRecordBatchReaderBuilder::new_with_metadata(input_file, metadata)
        .with_projection(projection)
        .build()
        .map_err(|error| {
            DataFusionError::Execution(format!(
                "failed to build parquet reader '{}': {error}",
                input.display()
            ))
        })?;

    let output_schema = std::sync::Arc::new(Schema::new(vec![
        Field::new("transcript_id", DataType::Utf8, false),
        Field::new("predictions", DataType::Binary, false),
    ]));
    let output_file = std::fs::File::create(output).map_err(|error| {
        DataFusionError::Execution(format!("failed to create '{}': {error}", output.display()))
    })?;
    let props = WriterProperties::builder()
        .set_max_row_group_size(row_group_size)
        .set_dictionary_enabled(false)
        .set_compression(compression)
        .build();
    let mut writer = ArrowWriter::try_new(output_file, output_schema.clone(), Some(props))?;

    let mut ids = Vec::with_capacity(row_group_size);
    let mut prediction_blobs = Vec::with_capacity(row_group_size);
    let mut stats = BuildStats {
        files: 1,
        ..BuildStats::default()
    };

    for batch_result in reader {
        let batch = batch_result.map_err(|error| {
            DataFusionError::Execution(format!("failed to read '{}': {error}", input.display()))
        })?;
        let schema = batch.schema();
        let tx_idx = schema.index_of("transcript_id").ok();
        let sift_idx = schema.index_of("sift_predictions").ok();
        let polyphen_idx = schema.index_of("polyphen_predictions").ok();
        let Some(tx_idx) = tx_idx else { continue };

        for row in 0..batch.num_rows() {
            let Some(transcript_id) = string_value(batch.column(tx_idx).as_ref(), row)? else {
                continue;
            };
            let mut predictions = CachedPredictions::default();
            if let Some(idx) = sift_idx {
                predictions.sift = read_compact_predictions(batch.column(idx).as_ref(), row);
            }
            if let Some(idx) = polyphen_idx {
                predictions.polyphen = read_compact_predictions(batch.column(idx).as_ref(), row);
            }
            predictions.sort();
            stats.sift_entries += predictions.sift.len() as u64;
            stats.polyphen_entries += predictions.polyphen.len() as u64;
            let bytes = serialize_predictions(&predictions);
            stats.prediction_bytes += bytes.len() as u64;
            stats.transcripts += 1;
            ids.push(transcript_id.to_string());
            prediction_blobs.push(bytes);

            if ids.len() >= row_group_size {
                flush(
                    &mut writer,
                    output_schema.clone(),
                    &mut ids,
                    &mut prediction_blobs,
                )?;
            }
        }
    }
    flush(&mut writer, output_schema, &mut ids, &mut prediction_blobs)?;
    writer.close()?;
    stats.elapsed_s = started.elapsed().as_secs_f64();
    Ok(stats)
}

fn flush(
    writer: &mut ArrowWriter<std::fs::File>,
    schema: std::sync::Arc<Schema>,
    ids: &mut Vec<String>,
    prediction_blobs: &mut Vec<Vec<u8>>,
) -> Result<()> {
    if ids.is_empty() {
        return Ok(());
    }
    let batch = RecordBatch::try_new(
        schema,
        vec![
            std::sync::Arc::new(StringArray::from(std::mem::take(ids))) as ArrayRef,
            std::sync::Arc::new(BinaryArray::from_iter_values(std::mem::take(
                prediction_blobs,
            ))) as ArrayRef,
        ],
    )
    .map_err(|error| DataFusionError::ArrowError(Box::new(error), None))?;
    writer.write(&batch)?;
    Ok(())
}

fn projection(metadata: &ArrowReaderMetadata, names: &[&str]) -> Result<ProjectionMask> {
    let fields = metadata.schema().fields();
    let roots = names
        .iter()
        .filter_map(|name| fields.iter().position(|field| field.name() == *name))
        .collect::<Vec<_>>();
    Ok(ProjectionMask::roots(metadata.parquet_schema(), roots))
}

fn string_value(array: &dyn Array, row: usize) -> Result<Option<&str>> {
    if row >= array.len() || array.is_null(row) {
        return Ok(None);
    }
    if let Some(array) = array.as_any().downcast_ref::<StringArray>() {
        Ok(Some(array.value(row)))
    } else if let Some(array) = array.as_any().downcast_ref::<StringViewArray>() {
        Ok(Some(array.value(row)))
    } else if let Some(array) = array.as_any().downcast_ref::<LargeStringArray>() {
        Ok(Some(array.value(row)))
    } else {
        Err(DataFusionError::Execution(format!(
            "transcript_id expected string array, got {:?}",
            array.data_type()
        )))
    }
}

fn read_compact_predictions(col: &dyn Array, row: usize) -> Vec<CompactPrediction> {
    if row >= col.len() || col.is_null(row) {
        return Vec::new();
    }
    let Some(list_arr) = col.as_any().downcast_ref::<ListArray>() else {
        return Vec::new();
    };
    let offsets = list_arr.offsets();
    let start_off = offsets[row] as usize;
    let end_off = offsets[row + 1] as usize;
    if start_off == end_off {
        return Vec::new();
    }
    let values = list_arr.values();
    let Some(struct_arr) = values
        .as_any()
        .downcast_ref::<datafusion::arrow::array::StructArray>()
    else {
        return Vec::new();
    };
    let Some(positions) = struct_arr.column_by_name("position") else {
        return Vec::new();
    };
    let Some(amino_acids) = struct_arr.column_by_name("amino_acid") else {
        return Vec::new();
    };
    let Some(predictions) = struct_arr.column_by_name("prediction") else {
        return Vec::new();
    };
    let Some(scores) = struct_arr.column_by_name("score") else {
        return Vec::new();
    };
    let pos_arr = positions.as_any().downcast_ref::<Int32Array>();
    let aa_arr = amino_acids.as_any().downcast_ref::<StringArray>();
    let aa_view = amino_acids.as_any().downcast_ref::<StringViewArray>();
    let aa_large = amino_acids.as_any().downcast_ref::<LargeStringArray>();
    let pred_arr = predictions.as_any().downcast_ref::<StringArray>();
    let pred_view = predictions.as_any().downcast_ref::<StringViewArray>();
    let pred_large = predictions.as_any().downcast_ref::<LargeStringArray>();
    let score_arr = scores.as_any().downcast_ref::<Float32Array>();

    let mut out = Vec::with_capacity(end_off - start_off);
    for i in start_off..end_off {
        let Some(pos) = pos_arr.and_then(|arr| (!arr.is_null(i)).then(|| arr.value(i))) else {
            continue;
        };
        let aa = aa_arr
            .and_then(|arr| (!arr.is_null(i)).then(|| arr.value(i)))
            .or_else(|| aa_view.and_then(|arr| (!arr.is_null(i)).then(|| arr.value(i))))
            .or_else(|| aa_large.and_then(|arr| (!arr.is_null(i)).then(|| arr.value(i))));
        let prediction = pred_arr
            .and_then(|arr| (!arr.is_null(i)).then(|| arr.value(i)))
            .or_else(|| pred_view.and_then(|arr| (!arr.is_null(i)).then(|| arr.value(i))))
            .or_else(|| pred_large.and_then(|arr| (!arr.is_null(i)).then(|| arr.value(i))));
        let score = score_arr.and_then(|arr| (!arr.is_null(i)).then(|| arr.value(i)));
        if let (Some(aa), Some(prediction), Some(score)) = (aa, prediction, score)
            && let Some(amino_acid) = CompactPrediction::encode_amino_acid(aa)
        {
            out.push(CompactPrediction {
                position: pos,
                amino_acid,
                prediction: CompactPrediction::encode_prediction(prediction),
                score,
            });
        }
    }
    out
}

fn serialize_predictions(predictions: &CachedPredictions) -> Vec<u8> {
    let sift_count = predictions.sift.len() as u32;
    let polyphen_count = predictions.polyphen.len() as u32;
    let mut buf = Vec::with_capacity(8 + (sift_count + polyphen_count) as usize * 10);
    buf.extend_from_slice(&sift_count.to_le_bytes());
    buf.extend_from_slice(&polyphen_count.to_le_bytes());
    for prediction in &predictions.sift {
        append_prediction(&mut buf, prediction);
    }
    for prediction in &predictions.polyphen {
        append_prediction(&mut buf, prediction);
    }
    buf
}

fn append_prediction(buf: &mut Vec<u8>, prediction: &CompactPrediction) {
    buf.extend_from_slice(&prediction.position.to_le_bytes());
    buf.push(prediction.amino_acid);
    buf.push(prediction.prediction);
    buf.extend_from_slice(&prediction.score.to_le_bytes());
}
