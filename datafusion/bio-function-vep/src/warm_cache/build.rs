use std::collections::BTreeSet;
use std::fs::File;
use std::path::{Path, PathBuf};
use std::sync::Arc;

use crate::kv_cache::position_index::{PositionIndex, position_index_file};
use crate::warm_cache::key::{position_key, variant_keys_from_allele_string};
use crate::warm_cache::reader::projection_for_existing_roots;
use crate::warm_cache::split::{
    FrequencyFields, WarmPositionCollector, max_af_from_pairs, max_global_af,
};
use datafusion::arrow::array::{
    Array, ArrayRef, BooleanBuilder, Float32Array, Float64Array, Int32Array, Int64Array,
    LargeStringArray, ListBuilder, RecordBatch, StringArray, StringViewArray, UInt32Array,
    UInt64Array, UInt64Builder,
};
use datafusion::arrow::compute::{concat_batches, filter_record_batch};
use datafusion::arrow::datatypes::{DataType, Field, FieldRef, Schema, SchemaRef};
use parquet::arrow::ArrowWriter;
use parquet::arrow::arrow_reader::ParquetRecordBatchReaderBuilder;
use parquet::basic::Compression;
use parquet::file::metadata::KeyValue;
use parquet::file::properties::WriterProperties;

type Result<T> = std::result::Result<T, Box<dyn std::error::Error + Send + Sync>>;

#[derive(Debug, Clone)]
pub struct WarmVariationTierOptions {
    pub input: PathBuf,
    pub output_dir: PathBuf,
    pub af_threshold: f64,
    pub position_radius: i64,
    pub row_group_rows: usize,
    pub batch_size: usize,
}

impl WarmVariationTierOptions {
    pub fn new(input: impl Into<PathBuf>, output_dir: impl Into<PathBuf>) -> Self {
        Self {
            input: input.into(),
            output_dir: output_dir.into(),
            af_threshold: 0.01,
            position_radius: 1,
            row_group_rows: 500_000,
            batch_size: 65_536,
        }
    }
}

#[derive(Debug, Clone)]
pub struct WarmVariationTierPartsOptions {
    pub inputs: Vec<PathBuf>,
    pub output_dir: PathBuf,
    pub work_dir: PathBuf,
    pub chrom: String,
    pub af_threshold: f64,
    pub position_radius: i64,
    pub row_group_rows: usize,
    pub batch_size: usize,
}

impl WarmVariationTierPartsOptions {
    pub fn new(
        inputs: Vec<PathBuf>,
        output_dir: impl Into<PathBuf>,
        work_dir: impl Into<PathBuf>,
        chrom: impl Into<String>,
    ) -> Self {
        Self {
            inputs,
            output_dir: output_dir.into(),
            work_dir: work_dir.into(),
            chrom: chrom.into(),
            af_threshold: 0.01,
            position_radius: 1,
            row_group_rows: 500_000,
            batch_size: 65_536,
        }
    }
}

#[derive(Debug, Clone, Default)]
pub struct WarmVariationTierStats {
    pub chrom: String,
    pub warm_positions: usize,
    pub warm_rows: usize,
    pub cold_rows: usize,
    pub warm_row_groups: usize,
    pub cold_row_groups: usize,
    pub cold_rows_sharing_warm_positions: usize,
    pub row_group_position_splits: usize,
}

#[derive(Debug)]
struct Args {
    input: PathBuf,
    output_dir: PathBuf,
    af_threshold: f64,
    position_radius: i64,
    row_group_rows: usize,
    batch_size: usize,
}

impl From<WarmVariationTierOptions> for Args {
    fn from(options: WarmVariationTierOptions) -> Self {
        Self {
            input: options.input,
            output_dir: options.output_dir,
            af_threshold: options.af_threshold,
            position_radius: options.position_radius,
            row_group_rows: options.row_group_rows,
            batch_size: options.batch_size,
        }
    }
}

pub fn build_warm_variation_tier(
    options: WarmVariationTierOptions,
) -> Result<WarmVariationTierStats> {
    let args = Args::from(options);
    std::fs::create_dir_all(&args.output_dir)?;

    let chrom = chrom_from_input(&args.input)?;
    let warm_positions = collect_warm_positions(&args)?;
    let written = write_split_files(&args, &chrom, &warm_positions)?;

    Ok(WarmVariationTierStats {
        chrom,
        warm_positions: warm_positions.len(),
        warm_rows: written.warm_rows,
        cold_rows: written.cold_rows,
        warm_row_groups: written.warm_row_groups,
        cold_row_groups: written.cold_row_groups,
        cold_rows_sharing_warm_positions: written.cold_rows_sharing_warm_positions,
        row_group_position_splits: written.row_group_position_splits,
    })
}

pub fn build_warm_variation_tier_from_parts(
    options: WarmVariationTierPartsOptions,
) -> Result<WarmVariationTierStats> {
    if options.inputs.is_empty() {
        return Err("warm/cold variation tier requires at least one input part".into());
    }

    std::fs::create_dir_all(&options.output_dir)?;
    if options.work_dir.exists() {
        std::fs::remove_dir_all(&options.work_dir)?;
    }
    std::fs::create_dir_all(&options.work_dir)?;

    let args_for = |input: PathBuf| Args {
        input,
        output_dir: options.output_dir.clone(),
        af_threshold: options.af_threshold,
        position_radius: options.position_radius,
        row_group_rows: options.row_group_rows,
        batch_size: options.batch_size,
    };

    let mut collect_handles = Vec::with_capacity(options.inputs.len());
    for input in &options.inputs {
        let args = args_for(input.clone());
        collect_handles.push(std::thread::spawn(move || collect_warm_positions(&args)));
    }

    let mut warm_positions = BTreeSet::new();
    for handle in collect_handles {
        let positions = handle
            .join()
            .map_err(|_| "warm position collection thread panicked")??;
        warm_positions.extend(positions);
    }

    let warm_positions = Arc::new(warm_positions);
    let mut split_handles = Vec::with_capacity(options.inputs.len());
    for (idx, input) in options.inputs.iter().enumerate() {
        let args = args_for(input.clone());
        let chrom = options.chrom.clone();
        let warm_positions = Arc::clone(&warm_positions);
        let warm_path = options.work_dir.join(format!("part_{idx}_warm.parquet"));
        let cold_path = options.work_dir.join(format!("part_{idx}_cold.parquet"));
        split_handles.push(std::thread::spawn(move || {
            let written =
                write_split_files_to_paths(&args, &chrom, &warm_positions, &warm_path, &cold_path)?;
            Ok::<_, Box<dyn std::error::Error + Send + Sync>>((idx, warm_path, cold_path, written))
        }));
    }

    let mut split_parts = Vec::with_capacity(options.inputs.len());
    let mut cold_rows_sharing_warm_positions = 0usize;
    for handle in split_handles {
        let (idx, warm_path, cold_path, written) = handle
            .join()
            .map_err(|_| "warm/cold split thread panicked")??;
        cold_rows_sharing_warm_positions += written.cold_rows_sharing_warm_positions;
        split_parts.push((idx, warm_path, cold_path));
    }
    split_parts.sort_by_key(|(idx, _, _)| *idx);

    let warm_part_paths: Vec<PathBuf> = split_parts
        .iter()
        .map(|(_, warm_path, _)| warm_path.clone())
        .collect();
    let cold_part_paths: Vec<PathBuf> = split_parts
        .iter()
        .map(|(_, _, cold_path)| cold_path.clone())
        .collect();

    let merge_args = args_for(options.inputs[0].clone());
    let warm_path = options
        .output_dir
        .join(format!("{}_warm.parquet", options.chrom));
    let cold_path = options
        .output_dir
        .join(format!("{}_cold.parquet", options.chrom));

    let (warm_rows, warm_row_groups) =
        merge_position_aligned_files(&warm_part_paths, &warm_path, "warm", &merge_args)?;
    let (cold_rows, cold_row_groups) =
        merge_position_aligned_files(&cold_part_paths, &cold_path, "cold", &merge_args)?;
    let row_group_position_splits =
        verify_position_aligned_row_groups(&warm_path, options.batch_size)?
            + verify_position_aligned_row_groups(&cold_path, options.batch_size)?;
    write_cold_position_index(
        &options.output_dir,
        &options.chrom,
        &cold_path,
        options.batch_size,
    )?;

    if let Err(error) = std::fs::remove_dir_all(&options.work_dir) {
        eprintln!(
            "warning: failed to remove warm/cold tier work dir {}: {error}",
            options.work_dir.display()
        );
    }

    Ok(WarmVariationTierStats {
        chrom: options.chrom,
        warm_positions: warm_positions.len(),
        warm_rows,
        cold_rows,
        warm_row_groups,
        cold_row_groups,
        cold_rows_sharing_warm_positions,
        row_group_position_splits,
    })
}

#[derive(Debug, Default)]
struct WrittenRows {
    warm_rows: usize,
    cold_rows: usize,
    warm_row_groups: usize,
    cold_row_groups: usize,
    cold_rows_sharing_warm_positions: usize,
    row_group_position_splits: usize,
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
    let warm_path = args.output_dir.join(format!("{chrom}_warm.parquet"));
    let cold_path = args.output_dir.join(format!("{chrom}_cold.parquet"));
    let written = write_split_files_to_paths(args, chrom, warm_positions, &warm_path, &cold_path)?;
    write_cold_position_index(&args.output_dir, chrom, &cold_path, args.batch_size)?;
    Ok(written)
}

fn write_split_files_to_paths(
    args: &Args,
    chrom: &str,
    warm_positions: &BTreeSet<i64>,
    warm_path: &Path,
    cold_path: &Path,
) -> Result<WrittenRows> {
    let file = File::open(&args.input)?;
    let builder = ParquetRecordBatchReaderBuilder::try_new(file)?;
    let source_schema = builder.schema().clone();
    let output_schema = output_schema(&source_schema)?;

    let warm_writer = create_writer(warm_path, output_schema.clone(), "warm", args)?;
    let cold_writer = create_writer(cold_path, output_schema.clone(), "cold", args)?;
    let mut warm_writer =
        PositionAlignedWriter::new(warm_writer, output_schema.clone(), args.row_group_rows)?;
    let mut cold_writer =
        PositionAlignedWriter::new(cold_writer, output_schema.clone(), args.row_group_rows)?;

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
            written.cold_rows_sharing_warm_positions +=
                count_rows_with_warm_positions(&cold_batch, warm_positions)?;
            written.cold_rows += cold_batch.num_rows();
            cold_writer.write(&cold_batch)?;
        }
    }

    written.warm_row_groups = warm_writer.close()?;
    written.cold_row_groups = cold_writer.close()?;
    written.row_group_position_splits =
        verify_position_aligned_row_groups(warm_path, args.batch_size)?
            + verify_position_aligned_row_groups(cold_path, args.batch_size)?;
    Ok(written)
}

fn merge_position_aligned_files(
    inputs: &[PathBuf],
    output_path: &Path,
    tier: &str,
    args: &Args,
) -> Result<(usize, usize)> {
    let first = inputs
        .first()
        .ok_or("position-aligned merge requires at least one input file")?;
    let first_file = File::open(first)?;
    let first_builder = ParquetRecordBatchReaderBuilder::try_new(first_file)?;
    let schema = first_builder.schema().clone();
    drop(first_builder);

    let writer = create_writer(output_path, schema.clone(), tier, args)?;
    let mut writer = PositionAlignedWriter::new(writer, schema, args.row_group_rows)?;
    let mut rows = 0usize;

    for input in inputs {
        let file = File::open(input)?;
        let reader = ParquetRecordBatchReaderBuilder::try_new(file)?
            .with_batch_size(args.batch_size)
            .build()?;
        for batch in reader {
            let batch = batch?;
            if batch.num_rows() == 0 {
                continue;
            }
            rows += batch.num_rows();
            writer.write(&batch)?;
        }
    }

    let row_groups = writer.close()?;
    Ok((rows, row_groups))
}

fn write_cold_position_index(
    output_dir: &Path,
    chrom: &str,
    cold_path: &Path,
    batch_size: usize,
) -> Result<()> {
    let index_dir = position_index_output_dir(output_dir);
    let index_path = position_index_file(&index_dir, chrom);
    let index = PositionIndex::from_parquet(cold_path, batch_size)?;
    index.write_to_path(&index_path)?;
    eprintln!(
        "cold_position_index={} positions={} source={}",
        index_path.display(),
        index.len(),
        cold_path.display()
    );
    Ok(())
}

fn position_index_output_dir(output_dir: &Path) -> PathBuf {
    if output_dir.file_name().and_then(|name| name.to_str()) == Some("variation") {
        output_dir
            .parent()
            .unwrap_or(output_dir)
            .join("variation.position_index")
    } else {
        output_dir.join("variation.position_index")
    }
}

struct PositionAlignedWriter {
    writer: ArrowWriter<File>,
    schema: SchemaRef,
    position_idx: usize,
    target_rows: usize,
    pending: Vec<RecordBatch>,
    pending_rows: usize,
    row_groups: usize,
}

impl PositionAlignedWriter {
    fn new(writer: ArrowWriter<File>, schema: SchemaRef, target_rows: usize) -> Result<Self> {
        let position_idx = schema.index_of("position_key")?;
        Ok(Self {
            writer,
            schema,
            position_idx,
            target_rows: target_rows.max(1),
            pending: Vec::new(),
            pending_rows: 0,
            row_groups: 0,
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
        self.writer.close()?;
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
        if self.pending_rows < self.target_rows {
            return Ok(None);
        }

        let mut previous = None;
        let mut global_row = 0usize;
        for batch in &self.pending {
            let positions = position_key_array(batch, self.position_idx)?;
            for row in 0..batch.num_rows() {
                if positions.is_null(row) {
                    return Err("position_key must be non-null UInt64".into());
                }
                let current = positions.value(row);
                if global_row >= self.target_rows && previous.is_some_and(|prev| prev != current) {
                    return Ok(Some(global_row));
                }
                previous = Some(current);
                global_row += 1;
            }
        }

        Ok(None)
    }

    fn flush_prefix(&mut self, cutoff: usize) -> Result<()> {
        debug_assert!(cutoff <= self.pending_rows);
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
            return Err("internal error while splitting position-aligned row group".into());
        }

        let batch = concat_batches(&self.schema, flush_batches.iter())?;
        self.writer.write(&batch)?;
        self.writer.flush()?;
        self.row_groups += 1;
        self.pending = keep_batches;
        self.pending_rows -= cutoff;
        Ok(())
    }
}

fn position_key_array(batch: &RecordBatch, idx: usize) -> Result<&UInt64Array> {
    batch
        .column(idx)
        .as_any()
        .downcast_ref::<UInt64Array>()
        .ok_or_else(|| "position_key must be UInt64".into())
}

fn count_rows_with_warm_positions(
    batch: &RecordBatch,
    warm_positions: &BTreeSet<i64>,
) -> Result<usize> {
    let start_idx = batch.schema().index_of("start")?;
    let mut count = 0usize;
    for row in 0..batch.num_rows() {
        let start = int64_value(batch.column(start_idx).as_ref(), row)
            .ok_or("start must be non-null Int64/Int32/UInt32")?;
        if warm_positions.contains(&start) {
            count += 1;
        }
    }
    Ok(count)
}

fn verify_position_aligned_row_groups(path: &Path, batch_size: usize) -> Result<usize> {
    let file = File::open(path)?;
    let builder = ParquetRecordBatchReaderBuilder::try_new(file)?;
    let row_groups = builder.metadata().num_row_groups();
    let mut previous_max = None;
    let mut splits = 0usize;

    for row_group in 0..row_groups {
        let (min, max) = row_group_position_range(path, row_group, batch_size)?;
        if previous_max.is_some_and(|previous| previous >= min) {
            splits += 1;
        }
        previous_max = Some(max);
    }

    Ok(splits)
}

fn row_group_position_range(
    path: &Path,
    row_group: usize,
    batch_size: usize,
) -> Result<(u64, u64)> {
    let file = File::open(path)?;
    let builder = ParquetRecordBatchReaderBuilder::try_new(file)?;
    let mask = projection_for_existing_roots(
        builder.schema(),
        builder.parquet_schema(),
        &["position_key"],
    );
    let reader = builder
        .with_projection(mask)
        .with_row_groups(vec![row_group])
        .with_batch_size(batch_size)
        .build()?;

    let mut first = None;
    let mut last = None;
    for batch in reader {
        let batch = batch?;
        let idx = batch.schema().index_of("position_key")?;
        let positions = position_key_array(&batch, idx)?;
        for row in 0..batch.num_rows() {
            if positions.is_null(row) {
                return Err("position_key must be non-null UInt64".into());
            }
            let value = positions.value(row);
            first.get_or_insert(value);
            last = Some(value);
        }
    }

    match (first, last) {
        (Some(first), Some(last)) => Ok((first, last)),
        _ => Err(format!(
            "row group {row_group} in {} has no position_key values",
            path.display()
        )
        .into()),
    }
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
        .set_max_row_group_size(usize::MAX)
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

#[cfg(test)]
mod tests {
    use super::*;

    fn write_input_part(path: &Path, starts: Vec<i64>, mafs: Vec<f64>) {
        let rows = starts.len();
        let schema = Arc::new(Schema::new(vec![
            Field::new("chrom", DataType::Utf8, false),
            Field::new("start", DataType::Int64, false),
            Field::new("allele_string", DataType::Utf8, true),
            Field::new("minor_allele_freq", DataType::Float64, true),
        ]));
        let batch = RecordBatch::try_new(
            schema.clone(),
            vec![
                Arc::new(StringArray::from(vec!["1"; rows])) as ArrayRef,
                Arc::new(Int64Array::from(starts)) as ArrayRef,
                Arc::new(StringArray::from(vec!["A/G"; rows])) as ArrayRef,
                Arc::new(Float64Array::from(mafs)) as ArrayRef,
            ],
        )
        .unwrap();
        let file = File::create(path).unwrap();
        let mut writer = ArrowWriter::try_new(file, schema, None).unwrap();
        writer.write(&batch).unwrap();
        writer.close().unwrap();
    }

    fn count_rows(path: &Path) -> usize {
        let file = File::open(path).unwrap();
        let reader = ParquetRecordBatchReaderBuilder::try_new(file)
            .unwrap()
            .build()
            .unwrap();
        reader.map(|batch| batch.unwrap().num_rows()).sum()
    }

    #[test]
    fn build_warm_variation_tier_from_parts_writes_final_split_files() {
        let dir = tempfile::tempdir().unwrap();
        let output_dir = dir.path().join("variation");
        let work_dir = dir.path().join("work");
        std::fs::create_dir_all(&output_dir).unwrap();

        let part_0 = dir.path().join("part_0.parquet");
        let part_1 = dir.path().join("part_1.parquet");
        write_input_part(&part_0, vec![100, 200], vec![0.02, 0.0]);
        write_input_part(&part_1, vec![101, 300], vec![0.0, 0.0]);

        let mut options = WarmVariationTierPartsOptions::new(
            vec![part_0, part_1],
            &output_dir,
            &work_dir,
            "chr1",
        );
        options.af_threshold = 0.01;
        options.position_radius = 1;
        options.row_group_rows = 2;
        options.batch_size = 2;

        let stats = build_warm_variation_tier_from_parts(options).unwrap();

        let warm_path = output_dir.join("chr1_warm.parquet");
        let cold_path = output_dir.join("chr1_cold.parquet");
        assert_eq!(stats.warm_rows, 2);
        assert_eq!(stats.cold_rows, 2);
        assert_eq!(count_rows(&warm_path), 2);
        assert_eq!(count_rows(&cold_path), 2);
        assert!(
            dir.path()
                .join("variation.position_index/chr1.posidx")
                .exists()
        );
        assert!(!work_dir.exists());
    }
}
