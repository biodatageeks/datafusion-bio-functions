use std::collections::BTreeSet;
use std::fs::File;
use std::path::{Path, PathBuf};
use std::sync::Arc;
use std::time::Instant;

use datafusion::arrow::array::{
    Array, ArrayRef, BooleanBuilder, Float32Array, Float64Array, Int32Array, Int64Array,
    LargeStringArray, ListBuilder, RecordBatch, StringArray, StringViewArray, UInt32Array,
    UInt64Builder,
};
use datafusion::arrow::compute::filter_record_batch;
use datafusion::arrow::datatypes::{DataType, Field, FieldRef, Schema, SchemaRef};
use datafusion_bio_function_vep::warm_cache::key::{position_key, variant_keys_from_allele_string};
use datafusion_bio_function_vep::warm_cache::reader::projection_for_existing_roots;
use datafusion_bio_function_vep::warm_cache::split::{
    FrequencyFields, WarmPositionCollector, max_af_from_pairs, max_global_af,
};
use parquet::arrow::ArrowWriter;
use parquet::arrow::arrow_reader::ParquetRecordBatchReaderBuilder;
use parquet::basic::Compression;
use parquet::file::metadata::KeyValue;
use parquet::file::properties::WriterProperties;

type Result<T> = std::result::Result<T, Box<dyn std::error::Error + Send + Sync>>;

#[derive(Debug)]
struct Args {
    input: PathBuf,
    output_dir: PathBuf,
    af_threshold: f64,
    position_radius: i64,
    row_group_rows: usize,
    batch_size: usize,
}

fn main() -> Result<()> {
    let args = Args::parse()?;
    std::fs::create_dir_all(&args.output_dir)?;

    let chrom = chrom_from_input(&args.input)?;
    let started = Instant::now();
    let warm_positions = collect_warm_positions(&args)?;
    eprintln!(
        "selected {} warm positions in {:.2}s",
        warm_positions.len(),
        started.elapsed().as_secs_f64()
    );

    let written = write_split_files(&args, &chrom, &warm_positions)?;
    eprintln!(
        "wrote {} warm rows and {} cold rows in {:.2}s",
        written.warm_rows,
        written.cold_rows,
        started.elapsed().as_secs_f64()
    );

    Ok(())
}

#[derive(Debug, Default)]
struct WrittenRows {
    warm_rows: usize,
    cold_rows: usize,
}

impl Args {
    fn parse() -> Result<Self> {
        let mut input = None;
        let mut output_dir = None;
        let mut af_threshold = 0.01;
        let mut position_radius = 1_i64;
        let mut row_group_rows = 250_000_usize;
        let mut batch_size = 65_536_usize;

        let mut args = std::env::args().skip(1);
        while let Some(arg) = args.next() {
            match arg.as_str() {
                "--input" => input = Some(PathBuf::from(require_value(&mut args, "--input")?)),
                "--output-dir" => {
                    output_dir = Some(PathBuf::from(require_value(&mut args, "--output-dir")?))
                }
                "--af-threshold" => {
                    af_threshold = require_value(&mut args, "--af-threshold")?.parse()?
                }
                "--position-radius" => {
                    position_radius = require_value(&mut args, "--position-radius")?.parse()?
                }
                "--row-group-rows" => {
                    row_group_rows = require_value(&mut args, "--row-group-rows")?.parse()?
                }
                "--batch-size" => batch_size = require_value(&mut args, "--batch-size")?.parse()?,
                "--help" | "-h" => {
                    print_usage();
                    std::process::exit(0);
                }
                other => return Err(format!("unknown argument: {other}").into()),
            }
        }

        Ok(Self {
            input: input.ok_or("--input is required")?,
            output_dir: output_dir.ok_or("--output-dir is required")?,
            af_threshold,
            position_radius,
            row_group_rows,
            batch_size,
        })
    }
}

fn require_value(args: &mut impl Iterator<Item = String>, name: &str) -> Result<String> {
    args.next()
        .ok_or_else(|| format!("{name} requires a value").into())
}

fn print_usage() {
    eprintln!(
        "Usage: build_warm_variation_cache --input chrN.parquet --output-dir variation \
         [--af-threshold 0.01] [--position-radius 1] [--row-group-rows 250000]"
    );
}

fn chrom_from_input(input: &Path) -> Result<String> {
    let stem = input.file_stem().and_then(|s| s.to_str()).ok_or_else(|| {
        format!(
            "cannot derive chromosome from input path {}",
            input.display()
        )
    })?;
    Ok(stem
        .strip_suffix("_warm")
        .or_else(|| stem.strip_suffix("_cold"))
        .unwrap_or(stem)
        .to_string())
}

fn collect_warm_positions(args: &Args) -> Result<BTreeSet<i64>> {
    let file = File::open(&args.input)?;
    let builder = ParquetRecordBatchReaderBuilder::try_new(file)?;
    let mask = projection_for_existing_roots(
        builder.schema(),
        builder.parquet_schema(),
        &["start", "minor_allele_freq", "AF", "gnomADg", "gnomADe"],
    );
    let reader = builder
        .with_projection(mask)
        .with_batch_size(args.batch_size)
        .build()?;

    let mut collector = WarmPositionCollector::new(args.af_threshold, args.position_radius);
    for batch in reader {
        let batch = batch?;
        let columns = FrequencyColumnIndices::new(batch.schema())?;
        for row in 0..batch.num_rows() {
            let start = int64_value(batch.column(columns.start).as_ref(), row)
                .ok_or("start must be non-null Int64/Int32/UInt32")?;
            let fields = FrequencyFields {
                minor_allele_freq: columns
                    .minor_allele_freq
                    .and_then(|idx| numeric_value(batch.column(idx).as_ref(), row)),
                af: columns
                    .af
                    .and_then(|idx| string_value(batch.column(idx).as_ref(), row)),
                gnomadg: columns
                    .gnomadg
                    .and_then(|idx| string_value(batch.column(idx).as_ref(), row)),
                gnomade: columns
                    .gnomade
                    .and_then(|idx| string_value(batch.column(idx).as_ref(), row)),
            };
            collector.push(start, max_global_af(&fields))?;
        }
    }

    Ok(collector.finish())
}

fn write_split_files(
    args: &Args,
    chrom: &str,
    warm_positions: &BTreeSet<i64>,
) -> Result<WrittenRows> {
    let file = File::open(&args.input)?;
    let builder = ParquetRecordBatchReaderBuilder::try_new(file)?;
    let source_schema = builder.schema().clone();
    let output_schema = output_schema(&source_schema)?;

    let warm_path = args.output_dir.join(format!("{chrom}_warm.parquet"));
    let cold_path = args.output_dir.join(format!("{chrom}_cold.parquet"));
    let mut warm_writer = create_writer(&warm_path, output_schema.clone(), "warm", args)?;
    let mut cold_writer = create_writer(&cold_path, output_schema.clone(), "cold", args)?;

    let reader = builder.with_batch_size(args.batch_size).build()?;
    let mut written = WrittenRows::default();

    for batch in reader {
        let batch = batch?;
        let enriched = append_key_columns(&batch, &output_schema, chrom)?;
        let warm_mask = warm_mask(&batch, warm_positions)?;
        let cold_mask = invert_mask(&warm_mask);

        let warm_batch = filter_record_batch(&enriched, &warm_mask)?;
        if warm_batch.num_rows() > 0 {
            written.warm_rows += warm_batch.num_rows();
            warm_writer.write(&warm_batch)?;
        }

        let cold_batch = filter_record_batch(&enriched, &cold_mask)?;
        if cold_batch.num_rows() > 0 {
            written.cold_rows += cold_batch.num_rows();
            cold_writer.write(&cold_batch)?;
        }
    }

    warm_writer.close()?;
    cold_writer.close()?;
    Ok(written)
}

#[derive(Debug)]
struct FrequencyColumnIndices {
    start: usize,
    minor_allele_freq: Option<usize>,
    af: Option<usize>,
    gnomadg: Option<usize>,
    gnomade: Option<usize>,
}

impl FrequencyColumnIndices {
    fn new(schema: SchemaRef) -> Result<Self> {
        Ok(Self {
            start: schema
                .index_of("start")
                .map_err(|_| "required column 'start' not found")?,
            minor_allele_freq: schema.index_of("minor_allele_freq").ok(),
            af: schema.index_of("AF").ok(),
            gnomadg: schema.index_of("gnomADg").ok(),
            gnomade: schema.index_of("gnomADe").ok(),
        })
    }
}

fn output_schema(source_schema: &SchemaRef) -> Result<SchemaRef> {
    if source_schema.index_of("position_key").is_ok()
        || source_schema.index_of("variant_keys").is_ok()
    {
        return Err("input already contains position_key or variant_keys".into());
    }

    let mut fields: Vec<FieldRef> = source_schema.fields().iter().cloned().collect();
    fields.push(Arc::new(Field::new(
        "position_key",
        DataType::UInt64,
        false,
    )));
    fields.push(Arc::new(Field::new_list(
        "variant_keys",
        Arc::new(Field::new_list_field(DataType::UInt64, true)),
        false,
    )));
    Ok(Arc::new(Schema::new_with_metadata(
        fields,
        source_schema.metadata().clone(),
    )))
}

fn create_writer(
    path: &Path,
    schema: SchemaRef,
    tier: &str,
    args: &Args,
) -> Result<ArrowWriter<File>> {
    let props = WriterProperties::builder()
        .set_compression(Compression::ZSTD(Default::default()))
        .set_max_row_group_size(args.row_group_rows)
        .set_key_value_metadata(Some(vec![
            KeyValue::new("vepyr.cache_tier".to_string(), tier.to_string()),
            KeyValue::new(
                "vepyr.warm_selector".to_string(),
                format!(
                    "max_global_af>={},+/-{}",
                    args.af_threshold, args.position_radius
                ),
            ),
            KeyValue::new("vepyr.key_version".to_string(), "1".to_string()),
            KeyValue::new(
                "vepyr.row_group_rows".to_string(),
                args.row_group_rows.to_string(),
            ),
        ]))
        .build();
    let file = File::create(path)?;
    Ok(ArrowWriter::try_new(file, schema, Some(props))?)
}

fn append_key_columns(
    batch: &RecordBatch,
    schema: &SchemaRef,
    default_chrom: &str,
) -> Result<RecordBatch> {
    let start_idx = batch.schema().index_of("start")?;
    let allele_idx = batch.schema().index_of("allele_string")?;
    let chrom_idx = batch.schema().index_of("chrom").ok();

    let mut position_builder = UInt64Builder::with_capacity(batch.num_rows());
    let mut variants_builder = ListBuilder::new(UInt64Builder::new());

    for row in 0..batch.num_rows() {
        let start = int64_value(batch.column(start_idx).as_ref(), row)
            .ok_or("start must be non-null Int64/Int32/UInt32")?;
        let chrom = chrom_idx
            .and_then(|idx| string_value(batch.column(idx).as_ref(), row))
            .unwrap_or(default_chrom);
        let allele_string = string_value(batch.column(allele_idx).as_ref(), row).unwrap_or("");

        position_builder.append_value(position_key(chrom, start)?);
        let variant_keys = variant_keys_from_allele_string(chrom, start, allele_string)?;
        for key in variant_keys {
            variants_builder.values().append_value(key);
        }
        variants_builder.append(true);
    }

    let mut columns: Vec<ArrayRef> = batch.columns().to_vec();
    columns.push(Arc::new(position_builder.finish()));
    columns.push(Arc::new(variants_builder.finish()));

    Ok(RecordBatch::try_new(schema.clone(), columns)?)
}

fn warm_mask(
    batch: &RecordBatch,
    warm_positions: &BTreeSet<i64>,
) -> Result<datafusion::arrow::array::BooleanArray> {
    let start_idx = batch.schema().index_of("start")?;
    let mut builder = BooleanBuilder::with_capacity(batch.num_rows());
    for row in 0..batch.num_rows() {
        let start = int64_value(batch.column(start_idx).as_ref(), row)
            .ok_or("start must be non-null Int64/Int32/UInt32")?;
        builder.append_value(warm_positions.contains(&start));
    }
    Ok(builder.finish())
}

fn invert_mask(
    mask: &datafusion::arrow::array::BooleanArray,
) -> datafusion::arrow::array::BooleanArray {
    let mut builder = BooleanBuilder::with_capacity(mask.len());
    for row in 0..mask.len() {
        builder.append_value(!mask.value(row));
    }
    builder.finish()
}

fn int64_value(array: &dyn Array, row: usize) -> Option<i64> {
    if array.is_null(row) {
        return None;
    }
    if let Some(a) = array.as_any().downcast_ref::<Int64Array>() {
        Some(a.value(row))
    } else if let Some(a) = array.as_any().downcast_ref::<Int32Array>() {
        Some(a.value(row) as i64)
    } else {
        array
            .as_any()
            .downcast_ref::<UInt32Array>()
            .map(|a| a.value(row) as i64)
    }
}

fn numeric_value(array: &dyn Array, row: usize) -> Option<f64> {
    if array.is_null(row) {
        return None;
    }
    if let Some(a) = array.as_any().downcast_ref::<Float64Array>() {
        Some(a.value(row))
    } else if let Some(a) = array.as_any().downcast_ref::<Float32Array>() {
        Some(a.value(row) as f64)
    } else {
        string_value(array, row).map(|s| max_af_from_pairs(Some(s)))
    }
}

fn string_value(array: &dyn Array, row: usize) -> Option<&str> {
    if array.is_null(row) {
        return None;
    }
    if let Some(a) = array.as_any().downcast_ref::<StringArray>() {
        Some(a.value(row))
    } else if let Some(a) = array.as_any().downcast_ref::<LargeStringArray>() {
        Some(a.value(row))
    } else if let Some(a) = array.as_any().downcast_ref::<StringViewArray>() {
        Some(a.value(row))
    } else {
        None
    }
}
