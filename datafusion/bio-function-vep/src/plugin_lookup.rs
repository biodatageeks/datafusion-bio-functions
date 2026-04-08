//! Plugin lookup engine.
//!
//! Runtime prefers plugin fjall stores for point lookups and falls back to the
//! older per-contig parquet-in-memory indexes when no fjall store is present.

use std::collections::HashMap;
use std::sync::Arc;

use datafusion::arrow::array::{
    Array, ArrayRef, Float32Array, Float32Builder, Int32Array, Int32Builder, RecordBatch,
    StringArray, StringBuilder, new_null_array,
};
use datafusion::arrow::datatypes::{DataType, Field, Schema, SchemaRef};
use datafusion::common::Result;
use datafusion::prelude::{ParquetReadOptions, SessionContext};

#[cfg(feature = "kv-cache")]
use crate::kv_cache::VepKvStore;
#[cfg(feature = "kv-cache")]
use crate::kv_cache::key_encoding::chrom_to_code;
#[cfg(feature = "kv-cache")]
use crate::kv_cache::position_entry::PositionEntryReader;
use crate::plugin::{ActivePlugins, PluginKind, PluginSourceKind};

enum PluginBackend {
    Parquet(ParquetPluginIndex),
    #[cfg(feature = "kv-cache")]
    Fjall(FjallPluginIndex),
}

/// Lookup backend for a single plugin on a single contig.
pub struct PluginIndex {
    pub kind: PluginSourceKind,
    backend: PluginBackend,
}

/// In-memory Parquet index for a single plugin on a single contig.
struct ParquetPluginIndex {
    /// Map from (pos, ref, alt) → row index in `data`.
    lookup: HashMap<(i64, Box<str>, Box<str>), u32>,
    /// Plugin value columns only (no chrom/pos/ref/alt).
    data: RecordBatch,
    /// Schema of the value columns.
    value_schema: SchemaRef,
}

impl PluginIndex {
    pub async fn load(
        ctx: &SessionContext,
        kind: PluginSourceKind,
        plugin_dir: &str,
        chrom: &str,
    ) -> Result<Option<Self>> {
        #[cfg(feature = "kv-cache")]
        if let Some(index) = Self::load_fjall(kind, plugin_dir, chrom)? {
            return Ok(Some(index));
        }
        Self::load_parquet(ctx, kind, plugin_dir, chrom).await
    }

    /// Load a plugin index from a per-chromosome parquet file.
    /// Returns `None` if the parquet file does not exist.
    pub async fn load_parquet(
        ctx: &SessionContext,
        kind: PluginSourceKind,
        plugin_dir: &str,
        chrom: &str,
    ) -> Result<Option<Self>> {
        let path = format!("{plugin_dir}/chr{chrom}.parquet");
        if !std::path::Path::new(&path).exists() {
            return Ok(None);
        }

        let df = ctx
            .read_parquet(&path, ParquetReadOptions::default())
            .await?;
        let batches = df.collect().await?;
        if batches.is_empty() {
            return Ok(None);
        }

        let full_schema = batches[0].schema();
        let join_cols = ["chrom", "pos", "ref", "alt"];
        let expected_value_names = kind
            .plugin_kind()
            .output_fields()
            .into_iter()
            .map(|field| field.name().to_string())
            .collect::<Vec<_>>();

        // Identify value columns (everything not in join key).
        let value_fields: Vec<Arc<Field>> = full_schema
            .fields()
            .iter()
            .filter(|f| {
                !join_cols.contains(&f.name().as_str())
                    && expected_value_names.iter().any(|name| name == f.name())
            })
            .cloned()
            .collect();
        let value_schema = Arc::new(Schema::new(value_fields));

        // Build lookup map and collect value columns.
        let mut lookup = HashMap::new();
        let mut value_batches: Vec<RecordBatch> = Vec::new();
        let mut global_offset: u32 = 0;

        for batch in &batches {
            let pos_col = get_i64_column(batch, "pos")?;
            let ref_col = get_string_column(batch, "ref")?;
            let alt_col = get_string_column(batch, "alt")?;

            for row in 0..batch.num_rows() {
                let pos = pos_col[row];
                let ref_allele: Box<str> = ref_col.value(row).into();
                let alt_allele: Box<str> = alt_col.value(row).into();
                lookup.insert((pos, ref_allele, alt_allele), global_offset + row as u32);
            }
            global_offset += batch.num_rows() as u32;

            // Extract value columns.
            let value_cols: Vec<ArrayRef> = value_schema
                .fields()
                .iter()
                .map(|f| {
                    let idx = full_schema.index_of(f.name()).unwrap();
                    batch.column(idx).clone()
                })
                .collect();
            value_batches.push(RecordBatch::try_new(value_schema.clone(), value_cols)?);
        }

        // Concatenate all value batches into one.
        let data = datafusion::arrow::compute::concat_batches(&value_schema, &value_batches)?;

        Ok(Some(Self {
            kind,
            backend: PluginBackend::Parquet(ParquetPluginIndex {
                lookup,
                data,
                value_schema,
            }),
        }))
    }

    #[cfg(feature = "kv-cache")]
    fn load_fjall(kind: PluginSourceKind, plugin_dir: &str, chrom: &str) -> Result<Option<Self>> {
        let path = std::path::Path::new(plugin_dir)
            .parent()
            .unwrap_or(std::path::Path::new(plugin_dir))
            .join(format!(
                "{}.fjall",
                std::path::Path::new(plugin_dir)
                    .file_name()
                    .and_then(|name| name.to_str())
                    .unwrap_or("plugin")
            ));
        if !path.exists() {
            return Ok(None);
        }

        let store = VepKvStore::open(&path)?;
        let value_fields: Vec<Arc<Field>> = store
            .schema()
            .fields()
            .iter()
            .filter(|field| {
                !["chrom", "start", "end", "allele_string"].contains(&field.name().as_str())
            })
            .cloned()
            .collect();
        let value_schema = Arc::new(Schema::new(value_fields));
        let value_col_indices = (0..value_schema.fields().len()).collect();

        Ok(Some(Self {
            kind,
            backend: PluginBackend::Fjall(FjallPluginIndex {
                store,
                chrom_code: chrom_to_code(chrom),
                value_schema,
                value_col_indices,
            }),
        }))
    }

    pub fn append_lookup_columns(
        &self,
        pos_values: &[i64],
        ref_col: &StringArray,
        alt_col: &StringArray,
        num_rows: usize,
        all_columns: &mut Vec<ArrayRef>,
        all_fields: &mut Vec<Arc<Field>>,
    ) -> Result<()> {
        match &self.backend {
            PluginBackend::Parquet(index) => {
                let mut row_indices: Vec<Option<u32>> = Vec::with_capacity(num_rows);
                for row in 0..num_rows {
                    row_indices.push(index.get(
                        pos_values[row],
                        ref_col.value(row),
                        alt_col.value(row),
                    ));
                }
                for (col_idx, field) in index.value_schema.fields().iter().enumerate() {
                    let source = index.data.column(col_idx);
                    all_columns.push(build_lookup_column(source, &row_indices, num_rows)?);
                    all_fields.push(field.clone());
                }
            }
            #[cfg(feature = "kv-cache")]
            PluginBackend::Fjall(index) => {
                let arrays = index.lookup_columns(pos_values, ref_col, alt_col)?;
                for (field, array) in index.value_schema.fields().iter().cloned().zip(arrays) {
                    all_fields.push(field);
                    all_columns.push(array);
                }
            }
        }
        Ok(())
    }

    pub fn value_schema(&self) -> &SchemaRef {
        match &self.backend {
            PluginBackend::Parquet(index) => &index.value_schema,
            #[cfg(feature = "kv-cache")]
            PluginBackend::Fjall(index) => &index.value_schema,
        }
    }

    pub fn csq_values_for_variant(
        &self,
        pos: i64,
        ref_allele: &str,
        alt_allele: &str,
    ) -> Vec<String> {
        match &self.backend {
            PluginBackend::Parquet(index) => {
                index.csq_values_for_variant(pos, ref_allele, alt_allele)
            }
            #[cfg(feature = "kv-cache")]
            PluginBackend::Fjall(index) => {
                index.csq_values_for_variant(pos, ref_allele, alt_allele)
            }
        }
    }
}

impl ParquetPluginIndex {
    fn get(&self, pos: i64, ref_allele: &str, alt_allele: &str) -> Option<u32> {
        let key = (
            pos,
            Box::<str>::from(ref_allele),
            Box::<str>::from(alt_allele),
        );
        self.lookup.get(&key).copied()
    }

    fn csq_values_for_variant(&self, pos: i64, ref_allele: &str, alt_allele: &str) -> Vec<String> {
        let row_idx = self.get(pos, ref_allele, alt_allele);
        (0..self.data.num_columns())
            .map(|col_idx| {
                row_idx
                    .and_then(|row| string_value(self.data.column(col_idx), row as usize))
                    .unwrap_or_default()
            })
            .collect()
    }
}

#[cfg(feature = "kv-cache")]
struct FjallPluginIndex {
    store: VepKvStore,
    chrom_code: u16,
    value_schema: SchemaRef,
    value_col_indices: Vec<usize>,
}

#[cfg(feature = "kv-cache")]
impl FjallPluginIndex {
    fn lookup_columns(
        &self,
        pos_values: &[i64],
        ref_col: &StringArray,
        alt_col: &StringArray,
    ) -> Result<Vec<ArrayRef>> {
        let num_rows = pos_values.len();
        let mut builders = self
            .value_schema
            .fields()
            .iter()
            .map(|field| crate::kv_cache::position_entry::make_builder(field.data_type(), num_rows))
            .collect::<Result<Vec<_>>>()?;

        let mut decompressor = self.store.create_decompressor()?;
        let mut buffer = Vec::new();

        for row in 0..num_rows {
            let allele_row = self.lookup_allele_row(
                pos_values[row],
                ref_col.value(row),
                alt_col.value(row),
                decompressor.as_mut(),
                &mut buffer,
            )?;

            for ((builder, field), col_idx) in builders
                .iter_mut()
                .zip(self.value_schema.fields().iter())
                .zip(self.value_col_indices.iter().copied())
            {
                append_plugin_value(
                    builder.as_mut(),
                    field.data_type(),
                    allele_row,
                    &buffer,
                    col_idx,
                )?;
            }
        }

        Ok(builders
            .into_iter()
            .map(|mut builder| builder.finish())
            .collect())
    }

    fn csq_values_for_variant(&self, pos: i64, ref_allele: &str, alt_allele: &str) -> Vec<String> {
        let mut decompressor = self.store.create_decompressor().ok().flatten();
        let mut buffer = Vec::new();
        let allele_row = self
            .lookup_allele_row(
                pos,
                ref_allele,
                alt_allele,
                decompressor.as_mut(),
                &mut buffer,
            )
            .ok()
            .flatten();

        self.value_col_indices
            .iter()
            .copied()
            .map(|col_idx| match allele_row {
                Some(row) => string_value_from_entry(&buffer, col_idx, row).unwrap_or_default(),
                None => String::new(),
            })
            .collect()
    }

    fn lookup_allele_row(
        &self,
        pos: i64,
        ref_allele: &str,
        alt_allele: &str,
        decompressor: Option<&mut zstd::bulk::Decompressor<'_>>,
        buffer: &mut Vec<u8>,
    ) -> Result<Option<usize>> {
        if !self
            .store
            .get_position_entry_fast(self.chrom_code, pos, decompressor, buffer)?
        {
            return Ok(None);
        }

        let reader = PositionEntryReader::new(buffer)?;
        let target = format!("{ref_allele}/{alt_allele}");
        Ok((0..reader.num_alleles()).find(|&row| reader.allele_string(row) == target))
    }
}

/// Collection of plugin indexes for a single contig.
pub struct ContigPlugins {
    pub active_plugins: ActivePlugins,
    pub indexes: Vec<PluginIndex>,
}

impl ContigPlugins {
    /// Append plugin columns to an annotation batch.
    ///
    /// For each row in `batch`, looks up (pos, ref, alt) in each plugin index.
    /// Matched rows get the plugin's value; unmatched rows get NULL.
    pub fn annotate_batch(&self, batch: RecordBatch) -> Result<RecordBatch> {
        if self.active_plugins.is_empty() {
            return Ok(batch);
        }

        let num_rows = batch.num_rows();
        let schema = batch.schema();

        // Extract lookup keys from the batch.
        let pos_col = get_i64_column(&batch, "start")?;
        let ref_col = get_string_column(&batch, "ref")?;
        let alt_col = get_string_column(&batch, "alt")?;

        // Start with existing columns.
        let mut all_columns: Vec<ArrayRef> = (0..schema.fields().len())
            .map(|i| batch.column(i).clone())
            .collect();
        let mut all_fields: Vec<Arc<Field>> = schema.fields().iter().cloned().collect();

        // For each configured plugin, append either real lookup values or NULLs.
        for plugin_cfg in &self.active_plugins.configs {
            if plugin_cfg.kind == PluginKind::Cadd {
                self.append_cadd_columns(
                    &pos_col,
                    ref_col,
                    alt_col,
                    num_rows,
                    &mut all_columns,
                    &mut all_fields,
                )?;
                continue;
            }

            let maybe_index = self
                .indexes
                .iter()
                .find(|index| index.kind.plugin_kind() == plugin_cfg.kind);

            if let Some(index) = maybe_index {
                index.append_lookup_columns(
                    &pos_col,
                    ref_col,
                    alt_col,
                    num_rows,
                    &mut all_columns,
                    &mut all_fields,
                )?;
            } else {
                for field in plugin_cfg.kind.output_fields() {
                    all_columns.push(new_null_array(field.data_type(), num_rows));
                    all_fields.push(Arc::new(field));
                }
            }
        }

        let new_schema = Arc::new(Schema::new(all_fields));
        Ok(RecordBatch::try_new(new_schema, all_columns)?)
    }

    /// Return a pipe-delimited CSQ suffix for the given variant.
    pub fn csq_suffix_for_variant(&self, pos: i64, ref_allele: &str, alt_allele: &str) -> String {
        let mut values = Vec::new();
        for plugin_cfg in &self.active_plugins.configs {
            if plugin_cfg.kind == PluginKind::Cadd {
                values.extend(self.cadd_values_for_variant(pos, ref_allele, alt_allele));
                continue;
            }

            if let Some(index) = self
                .indexes
                .iter()
                .find(|index| index.kind.plugin_kind() == plugin_cfg.kind)
            {
                values.extend(index.csq_values_for_variant(pos, ref_allele, alt_allele));
            } else {
                for _ in 0..plugin_cfg.kind.output_fields().len() {
                    values.push(String::new());
                }
            }
        }
        values.join("|")
    }

    fn append_cadd_columns(
        &self,
        pos_values: &[i64],
        ref_col: &StringArray,
        alt_col: &StringArray,
        num_rows: usize,
        all_columns: &mut Vec<ArrayRef>,
        all_fields: &mut Vec<Arc<Field>>,
    ) -> Result<()> {
        let mut raw_builder = Float32Builder::with_capacity(num_rows);
        let mut phred_builder = Float32Builder::with_capacity(num_rows);

        for row in 0..num_rows {
            let values = self.cadd_values_for_variant(
                pos_values[row],
                ref_col.value(row),
                alt_col.value(row),
            );
            append_optional_f32(&mut raw_builder, values.first().map(String::as_str));
            append_optional_f32(&mut phred_builder, values.get(1).map(String::as_str));
        }

        all_columns.push(Arc::new(raw_builder.finish()));
        all_fields.push(Arc::new(Field::new("raw_score", DataType::Float32, true)));
        all_columns.push(Arc::new(phred_builder.finish()));
        all_fields.push(Arc::new(Field::new("phred_score", DataType::Float32, true)));
        Ok(())
    }

    fn cadd_values_for_variant(&self, pos: i64, ref_allele: &str, alt_allele: &str) -> Vec<String> {
        self.indexes
            .iter()
            .find(|index| index.kind == PluginSourceKind::Cadd)
            .map(|index| index.csq_values_for_variant(pos, ref_allele, alt_allele))
            .unwrap_or_else(|| vec![String::new(), String::new()])
    }
}

fn append_optional_f32(builder: &mut Float32Builder, value: Option<&str>) {
    match value.and_then(|raw| (!raw.is_empty()).then_some(raw)) {
        Some(raw) => match raw.parse::<f32>() {
            Ok(parsed) => builder.append_value(parsed),
            Err(_) => builder.append_null(),
        },
        None => builder.append_null(),
    }
}

#[cfg(feature = "kv-cache")]
fn append_plugin_value(
    builder: &mut dyn datafusion::arrow::array::ArrayBuilder,
    data_type: &DataType,
    allele_row: Option<usize>,
    buffer: &[u8],
    col_idx: usize,
) -> Result<()> {
    match allele_row {
        Some(row) => {
            let reader = PositionEntryReader::new(buffer)?;
            match data_type {
                DataType::Float32 => {
                    let builder = builder
                        .as_any_mut()
                        .downcast_mut::<Float32Builder>()
                        .ok_or_else(|| {
                            datafusion::common::DataFusionError::Execution(
                                "expected Float32Builder".into(),
                            )
                        })?;
                    if let Some(value) = reader.read_f32_value(col_idx, row) {
                        builder.append_value(value);
                    } else {
                        builder.append_null();
                    }
                }
                DataType::Int32 => {
                    let builder = builder
                        .as_any_mut()
                        .downcast_mut::<Int32Builder>()
                        .ok_or_else(|| {
                            datafusion::common::DataFusionError::Execution(
                                "expected Int32Builder".into(),
                            )
                        })?;
                    if let Some(value) = reader.read_i64_value(col_idx, row) {
                        builder.append_value(value as i32);
                    } else {
                        builder.append_null();
                    }
                }
                DataType::Utf8 => {
                    let builder = builder
                        .as_any_mut()
                        .downcast_mut::<StringBuilder>()
                        .ok_or_else(|| {
                            datafusion::common::DataFusionError::Execution(
                                "expected StringBuilder".into(),
                            )
                        })?;
                    if let Some(value) = reader.read_string_value(col_idx, row) {
                        builder.append_value(value);
                    } else {
                        builder.append_null();
                    }
                }
                other => {
                    return Err(datafusion::common::DataFusionError::Execution(format!(
                        "Unsupported plugin fjall output type: {other}"
                    )));
                }
            }
        }
        None => match data_type {
            DataType::Float32 => builder
                .as_any_mut()
                .downcast_mut::<Float32Builder>()
                .ok_or_else(|| {
                    datafusion::common::DataFusionError::Execution("expected Float32Builder".into())
                })?
                .append_null(),
            DataType::Int32 => builder
                .as_any_mut()
                .downcast_mut::<Int32Builder>()
                .ok_or_else(|| {
                    datafusion::common::DataFusionError::Execution("expected Int32Builder".into())
                })?
                .append_null(),
            DataType::Utf8 => builder
                .as_any_mut()
                .downcast_mut::<StringBuilder>()
                .ok_or_else(|| {
                    datafusion::common::DataFusionError::Execution("expected StringBuilder".into())
                })?
                .append_null(),
            other => {
                return Err(datafusion::common::DataFusionError::Execution(format!(
                    "Unsupported plugin fjall output type: {other}"
                )));
            }
        },
    }
    Ok(())
}

#[cfg(feature = "kv-cache")]
fn string_value_from_entry(buffer: &[u8], col_idx: usize, row: usize) -> Option<String> {
    let reader = PositionEntryReader::new(buffer).ok()?;
    reader
        .read_string_value(col_idx, row)
        .or_else(|| {
            reader
                .read_i64_value(col_idx, row)
                .map(|value| value.to_string())
        })
        .or_else(|| {
            reader
                .read_f32_value(col_idx, row)
                .map(|value| value.to_string())
        })
}

/// Build an output column by looking up values from `source` at the given indices.
fn build_lookup_column(
    source: &ArrayRef,
    row_indices: &[Option<u32>],
    num_rows: usize,
) -> Result<ArrayRef> {
    match source.data_type() {
        DataType::Float32 => {
            let src = source.as_any().downcast_ref::<Float32Array>().unwrap();
            let mut builder = Float32Builder::with_capacity(num_rows);
            for idx in row_indices {
                match idx {
                    Some(i) if !src.is_null(*i as usize) => {
                        builder.append_value(src.value(*i as usize))
                    }
                    _ => builder.append_null(),
                }
            }
            Ok(Arc::new(builder.finish()))
        }
        DataType::Int32 => {
            let src = source.as_any().downcast_ref::<Int32Array>().unwrap();
            let mut builder = Int32Builder::with_capacity(num_rows);
            for idx in row_indices {
                match idx {
                    Some(i) if !src.is_null(*i as usize) => {
                        builder.append_value(src.value(*i as usize))
                    }
                    _ => builder.append_null(),
                }
            }
            Ok(Arc::new(builder.finish()))
        }
        DataType::Utf8 => {
            let src = source.as_any().downcast_ref::<StringArray>().unwrap();
            let mut builder = StringBuilder::with_capacity(num_rows, num_rows * 16);
            for idx in row_indices {
                match idx {
                    Some(i) if !src.is_null(*i as usize) => {
                        builder.append_value(src.value(*i as usize))
                    }
                    _ => builder.append_null(),
                }
            }
            Ok(Arc::new(builder.finish()))
        }
        dt => {
            // Fallback: create null array for unsupported types.
            Ok(new_null_array(dt, num_rows))
        }
    }
}

fn string_value(source: &ArrayRef, row_idx: usize) -> Option<String> {
    match source.data_type() {
        DataType::Float32 => source
            .as_any()
            .downcast_ref::<Float32Array>()
            .and_then(|array| (!array.is_null(row_idx)).then(|| array.value(row_idx).to_string())),
        DataType::Int32 => source
            .as_any()
            .downcast_ref::<Int32Array>()
            .and_then(|array| (!array.is_null(row_idx)).then(|| array.value(row_idx).to_string())),
        DataType::Utf8 => source
            .as_any()
            .downcast_ref::<StringArray>()
            .and_then(|array| (!array.is_null(row_idx)).then(|| array.value(row_idx).to_string())),
        _ => None,
    }
}

fn get_i64_column(batch: &RecordBatch, name: &str) -> Result<Vec<i64>> {
    let schema = batch.schema();
    let idx = schema.index_of(name).map_err(|e| {
        datafusion::common::DataFusionError::Execution(format!(
            "Column '{name}' not found in batch: {e}"
        ))
    })?;
    let col = batch.column(idx);

    // Try Int64 first, then UInt32 (VCF 'start' is often UInt32).
    if let Some(arr) = col
        .as_any()
        .downcast_ref::<datafusion::arrow::array::Int64Array>()
    {
        Ok(arr.values().to_vec())
    } else if let Some(arr) = col
        .as_any()
        .downcast_ref::<datafusion::arrow::array::UInt32Array>()
    {
        Ok(arr.values().iter().map(|v| *v as i64).collect())
    } else if let Some(arr) = col
        .as_any()
        .downcast_ref::<datafusion::arrow::array::Int32Array>()
    {
        Ok(arr.values().iter().map(|v| *v as i64).collect())
    } else {
        Err(datafusion::common::DataFusionError::Execution(format!(
            "Column '{name}' has unsupported type: {:?}",
            col.data_type()
        )))
    }
}

fn get_string_column<'a>(batch: &'a RecordBatch, name: &str) -> Result<&'a StringArray> {
    let schema = batch.schema();
    let idx = schema.index_of(name).map_err(|e| {
        datafusion::common::DataFusionError::Execution(format!(
            "Column '{name}' not found in batch: {e}"
        ))
    })?;
    batch
        .column(idx)
        .as_any()
        .downcast_ref::<StringArray>()
        .ok_or_else(|| {
            datafusion::common::DataFusionError::Execution(format!(
                "Column '{name}' is not StringArray"
            ))
        })
}

#[cfg(test)]
mod tests {
    use super::*;
    use datafusion::arrow::array::{Float32Array, Int64Array};
    use datafusion::arrow::datatypes::{Field, Schema};
    use std::sync::Arc;

    #[cfg(feature = "kv-cache")]
    #[test]
    fn load_fjall_plugin_and_lookup_exact_allele() {
        use crate::kv_cache::VepKvStore;
        use crate::kv_cache::position_entry::serialize_position_entry;
        use datafusion::arrow::array::StringArray;

        let temp = tempfile::tempdir().expect("tempdir");
        let plugin_dir = temp.path().join("cadd");
        std::fs::create_dir_all(&plugin_dir).expect("plugin dir");
        let fjall_dir = temp.path().join("cadd.fjall");

        let schema = Arc::new(Schema::new(vec![
            Field::new("chrom", DataType::Utf8, false),
            Field::new("start", DataType::Int64, false),
            Field::new("end", DataType::Int64, false),
            Field::new("allele_string", DataType::Utf8, false),
            Field::new("phred_score", DataType::Float32, true),
        ]));
        let batch = RecordBatch::try_new(
            schema.clone(),
            vec![
                Arc::new(StringArray::from(vec!["1"])),
                Arc::new(Int64Array::from(vec![123_i64])),
                Arc::new(Int64Array::from(vec![123_i64])),
                Arc::new(StringArray::from(vec!["A/G"])),
                Arc::new(Float32Array::from(vec![17.5_f32])),
            ],
        )
        .expect("batch");

        let store = VepKvStore::create(&fjall_dir, schema).expect("store");
        let value = serialize_position_entry(&[0], &batch, &[4], 3).expect("serialize");
        store.put_position_entry("1", 123, &value).expect("put");
        store.persist().expect("persist");
        drop(store);

        let index =
            PluginIndex::load_fjall(PluginSourceKind::Cadd, plugin_dir.to_str().unwrap(), "1")
                .expect("load")
                .expect("present");

        assert_eq!(index.csq_values_for_variant(123, "A", "G"), vec!["17.5"]);
        assert_eq!(index.csq_values_for_variant(123, "A", "T"), vec![""]);
    }

    #[test]
    fn cadd_values_are_loaded_from_single_backend() {
        let schema = Arc::new(Schema::new(vec![Field::new(
            "phred_score",
            DataType::Float32,
            true,
        )]));
        let cadd_index = PluginIndex {
            kind: PluginSourceKind::Cadd,
            backend: PluginBackend::Parquet(ParquetPluginIndex {
                lookup: HashMap::from([
                    ((101_i64, Box::<str>::from("A"), Box::<str>::from("G")), 0),
                    ((202_i64, Box::<str>::from("A"), Box::<str>::from("AT")), 1),
                ]),
                data: RecordBatch::try_new(
                    schema.clone(),
                    vec![Arc::new(Float32Array::from(vec![11.0_f32, 22.0_f32]))],
                )
                .expect("cadd batch"),
                value_schema: schema.clone(),
            }),
        };

        let plugins = ContigPlugins {
            active_plugins: ActivePlugins::default(),
            indexes: vec![cadd_index],
        };

        assert_eq!(plugins.cadd_values_for_variant(101, "A", "G"), vec!["11"]);
        assert_eq!(plugins.cadd_values_for_variant(202, "A", "AT"), vec!["22"]);
    }
}
