//! Plugin lookup engine.
//!
//! Runtime prefers plugin fjall stores for point lookups and falls back to the
//! older per-contig parquet-in-memory indexes when no fjall store is present.

use std::collections::{HashMap, HashSet};
use std::sync::Arc;

use datafusion::arrow::array::{
    Array, ArrayRef, Float32Array, Float32Builder, Int32Array, Int32Builder, LargeStringArray,
    RecordBatch, StringArray, StringBuilder, StringViewArray, new_null_array,
};
use datafusion::arrow::datatypes::{DataType, Field, Schema, SchemaRef};
use datafusion::common::Result;
use datafusion::prelude::{Expr, ParquetReadOptions, SessionContext, col, lit};

#[cfg(feature = "kv-cache")]
use crate::kv_cache::VepKvStore;
#[cfg(feature = "kv-cache")]
use crate::kv_cache::key_encoding::chrom_to_code;
#[cfg(feature = "kv-cache")]
use crate::kv_cache::position_entry::PositionEntryReader;
use crate::plugin::{ActivePlugins, PluginKind, PluginSourceKind};

pub type PluginTargetKey = (i64, String, String);

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
    /// Map from (pos, ref, alt) → row indices in `data`.
    lookup: HashMap<(i64, Box<str>, Box<str>), Vec<u32>>,
    /// Plugin value columns only (no chrom/pos/ref/alt).
    data: RecordBatch,
    /// Schema of the value columns.
    value_schema: SchemaRef,
    dbnsfp_match_data: Option<DbnsfpMatchData>,
    spliceai_match_data: Option<SpliceAiMatchData>,
}

struct DbnsfpMatchData {
    aapos: StringArray,
    aaref: StringArray,
    aaalt: StringArray,
}

struct SpliceAiMatchData {
    symbol: StringArray,
}

impl PluginIndex {
    pub async fn load(
        ctx: &SessionContext,
        kind: PluginSourceKind,
        plugin_dir: &str,
        chrom: &str,
        target_keys: Option<&HashSet<PluginTargetKey>>,
    ) -> Result<Option<Self>> {
        #[cfg(feature = "kv-cache")]
        if let Some(index) = Self::load_fjall(kind, plugin_dir, chrom)? {
            return Ok(Some(index));
        }
        Self::load_parquet(ctx, kind, plugin_dir, chrom, target_keys).await
    }

    /// Load a plugin index from a per-chromosome parquet file.
    /// Returns `None` if the parquet file does not exist.
    pub async fn load_parquet(
        ctx: &SessionContext,
        kind: PluginSourceKind,
        plugin_dir: &str,
        chrom: &str,
        target_keys: Option<&HashSet<PluginTargetKey>>,
    ) -> Result<Option<Self>> {
        let path = format!("{plugin_dir}/chr{chrom}.parquet");
        if !std::path::Path::new(&path).exists() {
            return Ok(None);
        }
        if target_keys.is_some_and(|keys| keys.is_empty()) {
            return Ok(None);
        }

        let mut df = ctx
            .read_parquet(&path, ParquetReadOptions::default())
            .await?;
        if let Some(keys) = target_keys {
            let pos_exprs = keys
                .iter()
                .filter_map(|(pos, _, _)| u32::try_from(*pos).ok())
                .collect::<HashSet<_>>()
                .into_iter()
                .map(|pos| lit(pos))
                .collect::<Vec<Expr>>();
            if pos_exprs.is_empty() {
                return Ok(None);
            }
            df = df.filter(col("pos").in_list(pos_exprs, false))?;
        }
        let batches = df.collect().await?;
        if batches.is_empty() {
            return Ok(None);
        }

        let full_schema = batches[0].schema();
        let join_cols = ["chrom", "pos", "ref", "alt"];
        let expected_fields = kind
            .plugin_kind()
            .output_fields()
            .into_iter()
            .map(|field| (field.name().to_string(), Arc::new(field)))
            .collect::<HashMap<_, _>>();

        // Keep runtime schema stable even if parquet was inferred as Utf8View/LargeUtf8.
        let value_fields: Vec<Arc<Field>> = full_schema
            .fields()
            .iter()
            .filter_map(|field| {
                (!join_cols.contains(&field.name().as_str()))
                    .then(|| expected_fields.get(field.name()).cloned())
                    .flatten()
            })
            .collect();
        let value_schema = Arc::new(Schema::new(value_fields));

        // Build lookup map and collect value columns.
        let mut lookup = HashMap::new();
        let mut value_batches: Vec<RecordBatch> = Vec::new();
        let dbnsfp_match_schema = (kind == PluginSourceKind::DbNSFP).then(|| {
            Arc::new(Schema::new(vec![
                Field::new("__vepyr_dbnsfp_aapos", DataType::Utf8, true),
                Field::new("__vepyr_dbnsfp_aaref", DataType::Utf8, true),
                Field::new("__vepyr_dbnsfp_aaalt", DataType::Utf8, true),
            ]))
        });
        let mut dbnsfp_match_batches: Vec<RecordBatch> = Vec::new();
        let spliceai_match_schema = (kind == PluginSourceKind::SpliceAI).then(|| {
            Arc::new(Schema::new(vec![Field::new(
                "symbol",
                DataType::Utf8,
                true,
            )]))
        });
        let mut spliceai_match_batches: Vec<RecordBatch> = Vec::new();
        let mut global_offset: u32 = 0;

        for batch in &batches {
            let pos_col = get_i64_column(batch, "pos")?;
            let ref_col = get_string_column(batch, "ref")?;
            let alt_col = get_string_column(batch, "alt")?;

            for row in 0..batch.num_rows() {
                let pos = pos_col[row];
                if target_keys.is_some_and(|keys| {
                    !keys.contains(&(
                        pos,
                        ref_col.value(row).to_string(),
                        alt_col.value(row).to_string(),
                    ))
                }) {
                    continue;
                }
                let ref_allele: Box<str> = ref_col.value(row).into();
                let alt_allele: Box<str> = alt_col.value(row).into();
                lookup
                    .entry((pos, ref_allele, alt_allele))
                    .or_insert_with(Vec::new)
                    .push(global_offset + row as u32);
            }
            global_offset += batch.num_rows() as u32;

            // Extract value columns.
            let value_cols: Vec<ArrayRef> = value_schema
                .fields()
                .iter()
                .map(|f| {
                    let idx = full_schema.index_of(f.name()).unwrap();
                    normalize_array_to_field(batch.column(idx), f)
                })
                .collect::<Result<_>>()?;
            value_batches.push(RecordBatch::try_new(value_schema.clone(), value_cols)?);

            if let Some(match_schema) = &dbnsfp_match_schema {
                let match_cols: Vec<ArrayRef> = match_schema
                    .fields()
                    .iter()
                    .map(|f| {
                        let idx = full_schema.index_of(f.name()).unwrap();
                        normalize_array_to_field(batch.column(idx), f)
                    })
                    .collect::<Result<_>>()?;
                dbnsfp_match_batches.push(RecordBatch::try_new(match_schema.clone(), match_cols)?);
            }
            if let Some(match_schema) = &spliceai_match_schema {
                let match_cols: Vec<ArrayRef> = match_schema
                    .fields()
                    .iter()
                    .map(|f| {
                        let idx = full_schema.index_of(f.name()).unwrap();
                        normalize_array_to_field(batch.column(idx), f)
                    })
                    .collect::<Result<_>>()?;
                spliceai_match_batches
                    .push(RecordBatch::try_new(match_schema.clone(), match_cols)?);
            }
        }

        // Concatenate all value batches into one.
        let data = datafusion::arrow::compute::concat_batches(&value_schema, &value_batches)?;
        let dbnsfp_match_data = if let Some(match_schema) = &dbnsfp_match_schema {
            let batch =
                datafusion::arrow::compute::concat_batches(match_schema, &dbnsfp_match_batches)?;
            Some(DbnsfpMatchData {
                aapos: get_string_column(&batch, "__vepyr_dbnsfp_aapos")?,
                aaref: get_string_column(&batch, "__vepyr_dbnsfp_aaref")?,
                aaalt: get_string_column(&batch, "__vepyr_dbnsfp_aaalt")?,
            })
        } else {
            None
        };
        let spliceai_match_data = if let Some(match_schema) = &spliceai_match_schema {
            let batch =
                datafusion::arrow::compute::concat_batches(match_schema, &spliceai_match_batches)?;
            Some(SpliceAiMatchData {
                symbol: get_string_column(&batch, "symbol")?,
            })
        } else {
            None
        };

        Ok(Some(Self {
            kind,
            backend: PluginBackend::Parquet(ParquetPluginIndex {
                lookup,
                data,
                value_schema,
                dbnsfp_match_data,
                spliceai_match_data,
            }),
        }))
    }

    #[cfg(feature = "kv-cache")]
    fn load_fjall(kind: PluginSourceKind, plugin_dir: &str, chrom: &str) -> Result<Option<Self>> {
        let plugin_path = std::path::Path::new(plugin_dir);
        let plugin_name = plugin_path
            .file_name()
            .and_then(|name| name.to_str())
            .unwrap_or("plugin");
        let parent = plugin_path.parent().unwrap_or(plugin_path);
        let path = if parent
            .parent()
            .and_then(|path| path.file_name())
            .and_then(|name| name.to_str())
            .is_some_and(|name| name == "parquet")
        {
            parent
                .parent()
                .and_then(|path| path.parent())
                .unwrap_or(parent)
                .join(format!("{plugin_name}.fjall"))
        } else {
            parent.join(format!("{plugin_name}.fjall"))
        };
        if !path.exists() {
            return Ok(None);
        }

        let store = VepKvStore::open(&path)?;
        let stored_entry_indices = store
            .schema()
            .fields()
            .iter()
            .enumerate()
            .filter(|(_, field)| !["chrom", "start"].contains(&field.name().as_str()))
            .map(|(idx, field)| (field.name().to_string(), idx - 2))
            .collect::<HashMap<_, _>>();
        let value_fields: Vec<Arc<Field>> = kind
            .plugin_kind()
            .output_fields()
            .into_iter()
            .map(Arc::new)
            .filter(|field| stored_entry_indices.contains_key(field.name()))
            .collect();
        let value_schema = Arc::new(Schema::new(value_fields));
        let value_col_indices = value_schema
            .fields()
            .iter()
            .map(|field| {
                stored_entry_indices
                    .get(field.name())
                    .copied()
                    .ok_or_else(|| {
                        datafusion::common::DataFusionError::Execution(format!(
                            "plugin fjall store missing field '{}'",
                            field.name()
                        ))
                    })
            })
            .collect::<Result<Vec<_>>>()?;

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

    pub fn csq_values_for_consequence(
        &self,
        pos: i64,
        ref_allele: &str,
        alt_allele: &str,
        consequence_symbol: Option<&str>,
        protein_position: Option<&str>,
        amino_acids: Option<&str>,
    ) -> Vec<String> {
        match &self.backend {
            PluginBackend::Parquet(index) => index.csq_values_for_consequence(
                pos,
                ref_allele,
                alt_allele,
                consequence_symbol,
                protein_position,
                amino_acids,
            ),
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
        self.lookup.get(&key).and_then(|rows| rows.last()).copied()
    }

    fn get_candidates(&self, pos: i64, ref_allele: &str, alt_allele: &str) -> Option<&[u32]> {
        let key = (
            pos,
            Box::<str>::from(ref_allele),
            Box::<str>::from(alt_allele),
        );
        self.lookup.get(&key).map(Vec::as_slice)
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

    fn csq_values_for_consequence(
        &self,
        pos: i64,
        ref_allele: &str,
        alt_allele: &str,
        consequence_symbol: Option<&str>,
        protein_position: Option<&str>,
        amino_acids: Option<&str>,
    ) -> Vec<String> {
        let row_idx = if self.spliceai_match_data.is_some() {
            self.matching_spliceai_row(pos, ref_allele, alt_allele, consequence_symbol)
        } else {
            self.matching_dbnsfp_row(pos, ref_allele, alt_allele, protein_position, amino_acids)
                .or_else(|| self.get(pos, ref_allele, alt_allele))
        };
        (0..self.data.num_columns())
            .map(|col_idx| {
                row_idx
                    .and_then(|row| string_value(self.data.column(col_idx), row as usize))
                    .unwrap_or_default()
            })
            .collect()
    }

    fn matching_dbnsfp_row(
        &self,
        pos: i64,
        ref_allele: &str,
        alt_allele: &str,
        protein_position: Option<&str>,
        amino_acids: Option<&str>,
    ) -> Option<u32> {
        let match_data = self.dbnsfp_match_data.as_ref()?;
        let candidates = self.get_candidates(pos, ref_allele, alt_allele)?;
        if candidates.len() <= 1 {
            return candidates.first().copied();
        }
        let (aa_ref, aa_alt) = amino_acids.and_then(|value| value.split_once('/'))?;
        let protein_position = protein_position.filter(|value| !value.is_empty())?;
        candidates.iter().copied().find(|row| {
            dbnsfp_row_matches_consequence(
                match_data.aapos.value(*row as usize),
                match_data.aaref.value(*row as usize),
                match_data.aaalt.value(*row as usize),
                protein_position,
                aa_ref,
                aa_alt,
            )
        })
    }

    fn matching_spliceai_row(
        &self,
        pos: i64,
        ref_allele: &str,
        alt_allele: &str,
        consequence_symbol: Option<&str>,
    ) -> Option<u32> {
        let match_data = self.spliceai_match_data.as_ref()?;
        let consequence_symbol = consequence_symbol.filter(|value| !value.is_empty())?;
        self.get_candidates(pos, ref_allele, alt_allele)?
            .iter()
            .copied()
            .find(|row| match_data.symbol.value(*row as usize) == consequence_symbol)
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
                    &ref_col,
                    &alt_col,
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
                    &ref_col,
                    &alt_col,
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
        self.csq_suffix_for_consequence(pos, ref_allele, alt_allele, true, true, None, None, None)
    }

    /// Return a pipe-delimited CSQ suffix for a specific consequence row.
    ///
    /// Some plugins are consequence-type-specific in VEP output even when the
    /// source lookup itself is variant-level. AlphaMissense is emitted only on
    /// matched missense rows. dbNSFP's VEP plugin defaults to a consequence
    /// filter for missense, stop-lost, stop-gained, and start-lost rows.
    pub fn csq_suffix_for_consequence(
        &self,
        pos: i64,
        ref_allele: &str,
        alt_allele: &str,
        include_alphamissense: bool,
        include_dbnsfp: bool,
        consequence_symbol: Option<&str>,
        protein_position: Option<&str>,
        amino_acids: Option<&str>,
    ) -> String {
        let mut values = Vec::new();
        for plugin_cfg in &self.active_plugins.configs {
            if plugin_cfg.kind == PluginKind::Cadd {
                values.extend(
                    self.cadd_values_for_variant(pos, ref_allele, alt_allele)
                        .into_iter()
                        .map(|value| sanitize_csq_plugin_value(&value)),
                );
                continue;
            }

            if plugin_cfg.kind == PluginKind::AlphaMissense && !include_alphamissense {
                for _ in 0..plugin_cfg.kind.output_fields().len() {
                    values.push(String::new());
                }
                continue;
            }

            if plugin_cfg.kind == PluginKind::DbNSFP && !include_dbnsfp {
                for _ in 0..plugin_cfg.kind.output_fields().len() {
                    values.push(String::new());
                }
                continue;
            }

            if let Some(index) = self
                .indexes
                .iter()
                .find(|index| index.kind.plugin_kind() == plugin_cfg.kind)
            {
                let plugin_values = if plugin_cfg.kind == PluginKind::DbNSFP
                    || plugin_cfg.kind == PluginKind::SpliceAI
                {
                    index.csq_values_for_consequence(
                        pos,
                        ref_allele,
                        alt_allele,
                        consequence_symbol,
                        protein_position,
                        amino_acids,
                    )
                } else {
                    index.csq_values_for_variant(pos, ref_allele, alt_allele)
                };
                values.extend(
                    plugin_values
                        .into_iter()
                        .map(|value| sanitize_plugin_value_for_kind(plugin_cfg.kind, &value)),
                );
            } else {
                for _ in 0..plugin_cfg.kind.output_fields().len() {
                    values.push(String::new());
                }
            }
        }
        values.join("|")
    }

    pub fn alphamissense_matches_transcript_consequence(
        &self,
        pos: i64,
        ref_allele: &str,
        alt_allele: &str,
        protein_position: Option<&str>,
        amino_acids: Option<&str>,
    ) -> bool {
        let Some(index) = self
            .indexes
            .iter()
            .find(|index| index.kind == PluginSourceKind::AlphaMissense)
        else {
            return false;
        };
        let values = index.csq_values_for_variant(pos, ref_allele, alt_allele);
        let Some(plugin_protein_variant) = values.get(3).map(String::as_str) else {
            return false;
        };
        protein_variant_matches(plugin_protein_variant, protein_position, amino_acids)
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

fn sanitize_csq_plugin_value(value: &str) -> String {
    value.replace('|', "&")
}

fn sanitize_dbnsfp_csq_plugin_value(value: &str) -> String {
    value
        .chars()
        .map(|ch| if ch == '|' || ch == ';' { '&' } else { ch })
        .collect()
}

fn sanitize_plugin_value_for_kind(kind: PluginKind, value: &str) -> String {
    if kind == PluginKind::DbNSFP {
        sanitize_dbnsfp_csq_plugin_value(value)
    } else {
        sanitize_csq_plugin_value(value)
    }
}

fn dbnsfp_row_matches_consequence(
    aapos_values: &str,
    aaref_values: &str,
    aaalt_values: &str,
    protein_position: &str,
    aa_ref: &str,
    aa_alt: &str,
) -> bool {
    let mut positions = aapos_values.split(';');
    let mut refs = aaref_values.split(';');
    let mut alts = aaalt_values.split(';');
    loop {
        let (Some(pos), Some(row_ref), Some(row_alt)) =
            (positions.next(), refs.next(), alts.next())
        else {
            return false;
        };
        if pos == protein_position && row_ref == aa_ref && row_alt == aa_alt {
            return true;
        }
    }
}

fn protein_variant_matches(
    plugin_protein_variant: &str,
    protein_position: Option<&str>,
    amino_acids: Option<&str>,
) -> bool {
    let Some((ref_aa, alt_aa)) = amino_acids.and_then(|value| value.split_once('/')) else {
        return false;
    };
    let Some(position) = protein_position.filter(|value| !value.is_empty()) else {
        return false;
    };
    if plugin_protein_variant.len() < 3 {
        return false;
    }
    let mut chars = plugin_protein_variant.chars();
    let Some(plugin_ref) = chars.next() else {
        return false;
    };
    let Some(plugin_alt) = plugin_protein_variant.chars().last() else {
        return false;
    };
    let plugin_pos = &plugin_protein_variant[1..plugin_protein_variant.len() - 1];
    plugin_pos == position
        && ref_aa.len() == 1
        && alt_aa.len() == 1
        && ref_aa.starts_with(plugin_ref)
        && alt_aa.starts_with(plugin_alt)
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
                DataType::Utf8 | DataType::Utf8View | DataType::LargeUtf8 => {
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
            DataType::Utf8 | DataType::Utf8View | DataType::LargeUtf8 => builder
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
        DataType::Utf8View => {
            let src = source.as_any().downcast_ref::<StringViewArray>().unwrap();
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
        DataType::LargeUtf8 => {
            let src = source.as_any().downcast_ref::<LargeStringArray>().unwrap();
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

fn normalize_array_to_field(source: &ArrayRef, field: &Field) -> Result<ArrayRef> {
    match (source.data_type(), field.data_type()) {
        (DataType::Utf8, DataType::Utf8) => Ok(source.clone()),
        (DataType::Utf8View, DataType::Utf8) => {
            let array = source
                .as_any()
                .downcast_ref::<StringViewArray>()
                .ok_or_else(|| {
                    datafusion::common::DataFusionError::Execution(
                        "expected StringViewArray".into(),
                    )
                })?;
            Ok(Arc::new(StringArray::from_iter(
                (0..array.len()).map(|i| (!array.is_null(i)).then(|| array.value(i))),
            )))
        }
        (DataType::LargeUtf8, DataType::Utf8) => {
            let array = source
                .as_any()
                .downcast_ref::<LargeStringArray>()
                .ok_or_else(|| {
                    datafusion::common::DataFusionError::Execution(
                        "expected LargeStringArray".into(),
                    )
                })?;
            Ok(Arc::new(StringArray::from_iter(
                (0..array.len()).map(|i| (!array.is_null(i)).then(|| array.value(i))),
            )))
        }
        (source_ty, field_ty) if source_ty == field_ty => Ok(source.clone()),
        (source_ty, field_ty) => Err(datafusion::common::DataFusionError::Execution(format!(
            "plugin field '{}' expected {field_ty} but source column had {source_ty}",
            field.name()
        ))),
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
        DataType::Utf8View => source
            .as_any()
            .downcast_ref::<StringViewArray>()
            .and_then(|array| (!array.is_null(row_idx)).then(|| array.value(row_idx).to_string())),
        DataType::LargeUtf8 => source
            .as_any()
            .downcast_ref::<LargeStringArray>()
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

fn get_string_column(batch: &RecordBatch, name: &str) -> Result<StringArray> {
    let schema = batch.schema();
    let idx = schema.index_of(name).map_err(|e| {
        datafusion::common::DataFusionError::Execution(format!(
            "Column '{name}' not found in batch: {e}"
        ))
    })?;
    let column = batch.column(idx);
    match column.data_type() {
        DataType::Utf8 => column
            .as_any()
            .downcast_ref::<StringArray>()
            .cloned()
            .ok_or_else(|| {
                datafusion::common::DataFusionError::Execution(format!(
                    "Column '{name}' is not StringArray"
                ))
            }),
        DataType::Utf8View => column
            .as_any()
            .downcast_ref::<StringViewArray>()
            .map(|array| {
                StringArray::from_iter(
                    (0..array.len()).map(|i| (!array.is_null(i)).then(|| array.value(i))),
                )
            })
            .ok_or_else(|| {
                datafusion::common::DataFusionError::Execution(format!(
                    "Column '{name}' is not StringViewArray"
                ))
            }),
        DataType::LargeUtf8 => column
            .as_any()
            .downcast_ref::<LargeStringArray>()
            .map(|array| {
                StringArray::from_iter(
                    (0..array.len()).map(|i| (!array.is_null(i)).then(|| array.value(i))),
                )
            })
            .ok_or_else(|| {
                datafusion::common::DataFusionError::Execution(format!(
                    "Column '{name}' is not LargeStringArray"
                ))
            }),
        other => Err(datafusion::common::DataFusionError::Execution(format!(
            "Column '{name}' has unsupported string type: {other}"
        ))),
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::plugin::PluginConfig;
    use datafusion::arrow::array::{
        Float32Array, Int32Array, Int64Array, StringArray, UInt32Array,
    };
    use datafusion::arrow::datatypes::{Field, Schema};
    use datafusion::parquet::arrow::ArrowWriter;
    use std::sync::Arc;

    #[tokio::test]
    async fn load_parquet_plugin_filters_to_requested_target_keys() {
        let temp = tempfile::tempdir().expect("tempdir");
        let plugin_dir = temp.path().join("spliceai");
        std::fs::create_dir_all(&plugin_dir).expect("plugin dir");
        let path = plugin_dir.join("chr1.parquet");

        let schema = Arc::new(Schema::new(vec![
            Field::new("chrom", DataType::Utf8, false),
            Field::new("pos", DataType::UInt32, false),
            Field::new("ref", DataType::Utf8, false),
            Field::new("alt", DataType::Utf8, false),
            Field::new("symbol", DataType::Utf8, true),
            Field::new("ds_ag", DataType::Float32, true),
            Field::new("dp_ag", DataType::Int32, true),
        ]));
        let batch = RecordBatch::try_new(
            schema.clone(),
            vec![
                Arc::new(StringArray::from(vec!["1", "1"])),
                Arc::new(UInt32Array::from(vec![65_420_u32, 65_421_u32])),
                Arc::new(StringArray::from(vec!["C", "G"])),
                Arc::new(StringArray::from(vec!["A", "T"])),
                Arc::new(StringArray::from(vec!["OR4F5", "SKIP"])),
                Arc::new(Float32Array::from(vec![0.12_f32, 0.99_f32])),
                Arc::new(Int32Array::from(vec![7_i32, 99_i32])),
            ],
        )
        .expect("batch");
        let file = std::fs::File::create(path).expect("create parquet");
        let mut writer = ArrowWriter::try_new(file, schema, None).expect("writer");
        writer.write(&batch).expect("write");
        writer.close().expect("close");

        let mut target_keys = HashSet::new();
        target_keys.insert((65_420_i64, "C".to_string(), "A".to_string()));

        let ctx = SessionContext::new();
        let index = PluginIndex::load_parquet(
            &ctx,
            PluginSourceKind::SpliceAI,
            plugin_dir.to_str().unwrap(),
            "1",
            Some(&target_keys),
        )
        .await
        .expect("load")
        .expect("present");

        assert_eq!(
            index.csq_values_for_variant(65_420, "C", "A"),
            vec!["OR4F5", "0.12", "7"]
        );
        assert_eq!(
            index.csq_values_for_variant(65_421, "G", "T"),
            vec!["", "", ""]
        );
    }

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

    #[cfg(feature = "kv-cache")]
    #[test]
    fn load_fjall_plugin_from_cache_root_for_partitioned_layout() {
        use crate::kv_cache::VepKvStore;
        use crate::kv_cache::position_entry::serialize_position_entry;
        use datafusion::arrow::array::StringArray;

        let temp = tempfile::tempdir().expect("tempdir");
        let plugin_dir = temp
            .path()
            .join("parquet")
            .join("115_GRCh38_vep")
            .join("cadd");
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
                    (
                        (101_i64, Box::<str>::from("A"), Box::<str>::from("G")),
                        vec![0],
                    ),
                    (
                        (202_i64, Box::<str>::from("A"), Box::<str>::from("AT")),
                        vec![1],
                    ),
                ]),
                data: RecordBatch::try_new(
                    schema.clone(),
                    vec![Arc::new(Float32Array::from(vec![11.0_f32, 22.0_f32]))],
                )
                .expect("cadd batch"),
                value_schema: schema.clone(),
                dbnsfp_match_data: None,
                spliceai_match_data: None,
            }),
        };

        let plugins = ContigPlugins {
            active_plugins: ActivePlugins::default(),
            indexes: vec![cadd_index],
        };

        assert_eq!(plugins.cadd_values_for_variant(101, "A", "G"), vec!["11"]);
        assert_eq!(plugins.cadd_values_for_variant(202, "A", "AT"), vec!["22"]);
    }

    #[test]
    fn csq_suffix_escapes_inner_pipe_characters_in_plugin_values() {
        let schema = Arc::new(Schema::new(vec![
            Field::new("ClinVar", DataType::Utf8, true),
            Field::new("ClinVar_CLNSIG", DataType::Utf8, true),
            Field::new("ClinVar_CLNREVSTAT", DataType::Utf8, true),
            Field::new("ClinVar_CLNDN", DataType::Utf8, true),
            Field::new("ClinVar_CLNVC", DataType::Utf8, true),
            Field::new("ClinVar_CLNVI", DataType::Utf8, true),
        ]));
        let clinvar_index = PluginIndex {
            kind: PluginSourceKind::ClinVar,
            backend: PluginBackend::Parquet(ParquetPluginIndex {
                lookup: HashMap::from([(
                    (101_i64, Box::<str>::from("A"), Box::<str>::from("G")),
                    vec![0],
                )]),
                data: RecordBatch::try_new(
                    schema.clone(),
                    vec![
                        Arc::new(StringArray::from(vec!["12345"])),
                        Arc::new(StringArray::from(vec!["Benign"])),
                        Arc::new(StringArray::from(vec![
                            "criteria_provided|_multiple_submitters|_no_conflicts",
                        ])),
                        Arc::new(StringArray::from(vec!["Disease_A|Disease_B"])),
                        Arc::new(StringArray::from(vec!["single_nucleotide_variant"])),
                        Arc::new(StringArray::from(vec!["ClinGen:CA1|Other:2"])),
                    ],
                )
                .expect("clinvar batch"),
                value_schema: schema,
                dbnsfp_match_data: None,
                spliceai_match_data: None,
            }),
        };

        let plugins = ContigPlugins {
            active_plugins: ActivePlugins {
                configs: vec![PluginConfig {
                    kind: PluginKind::ClinVar,
                    source_dirs: vec![],
                }],
            },
            indexes: vec![clinvar_index],
        };

        assert_eq!(
            plugins.csq_suffix_for_variant(101, "A", "G"),
            "12345|Benign|criteria_provided&_multiple_submitters&_no_conflicts|Disease_A&Disease_B|single_nucleotide_variant|ClinGen:CA1&Other:2"
        );
    }

    #[test]
    fn dbnsfp_csq_suffix_escapes_vep_style_inner_separators() {
        let schema = Arc::new(Schema::new(vec![
            Field::new("sift4g_score", DataType::Utf8, true),
            Field::new("sift4g_pred", DataType::Utf8, true),
        ]));
        let dbnsfp_index = PluginIndex {
            kind: PluginSourceKind::DbNSFP,
            backend: PluginBackend::Parquet(ParquetPluginIndex {
                lookup: HashMap::from([(
                    (101_i64, Box::<str>::from("A"), Box::<str>::from("G")),
                    vec![0],
                )]),
                data: RecordBatch::try_new(
                    schema.clone(),
                    vec![
                        Arc::new(StringArray::from(vec![".;0.354;0.471"])),
                        Arc::new(StringArray::from(vec!["T|D"])),
                    ],
                )
                .expect("dbnsfp batch"),
                value_schema: schema,
                dbnsfp_match_data: None,
                spliceai_match_data: None,
            }),
        };

        let plugins = ContigPlugins {
            active_plugins: ActivePlugins {
                configs: vec![PluginConfig {
                    kind: PluginKind::DbNSFP,
                    source_dirs: vec![],
                }],
            },
            indexes: vec![dbnsfp_index],
        };

        assert_eq!(
            plugins.csq_suffix_for_variant(101, "A", "G"),
            ".&0.354&0.471|T&D"
        );
    }

    #[test]
    fn dbnsfp_csq_suffix_selects_duplicate_row_by_protein_change() {
        let schema = Arc::new(Schema::new(vec![
            Field::new("sift4g_score", DataType::Utf8, true),
            Field::new("revel_score", DataType::Utf8, true),
        ]));
        let dbnsfp_index = PluginIndex {
            kind: PluginSourceKind::DbNSFP,
            backend: PluginBackend::Parquet(ParquetPluginIndex {
                lookup: HashMap::from([(
                    (254484_i64, Box::<str>::from("A"), Box::<str>::from("C")),
                    vec![0, 1],
                )]),
                data: RecordBatch::try_new(
                    schema.clone(),
                    vec![
                        Arc::new(StringArray::from(vec!["0.021", ".;.;0.136"])),
                        Arc::new(StringArray::from(vec!["0.238", ".;.;0.533"])),
                    ],
                )
                .expect("dbnsfp batch"),
                value_schema: schema,
                dbnsfp_match_data: Some(DbnsfpMatchData {
                    aapos: StringArray::from(vec!["112", "450;581;647"]),
                    aaref: StringArray::from(vec!["M", "Y;Y;Y"]),
                    aaalt: StringArray::from(vec!["L", "S;S;S"]),
                }),
                spliceai_match_data: None,
            }),
        };

        let plugins = ContigPlugins {
            active_plugins: ActivePlugins {
                configs: vec![PluginConfig {
                    kind: PluginKind::DbNSFP,
                    source_dirs: vec![],
                }],
            },
            indexes: vec![dbnsfp_index],
        };

        assert_eq!(
            plugins.csq_suffix_for_consequence(
                254484,
                "A",
                "C",
                true,
                true,
                None,
                Some("112"),
                Some("M/L"),
            ),
            "0.021|0.238"
        );
        assert_eq!(
            plugins.csq_suffix_for_consequence(
                254484,
                "A",
                "C",
                true,
                true,
                None,
                Some("647"),
                Some("Y/S"),
            ),
            ".&.&0.136|.&.&0.533"
        );
    }

    #[test]
    fn spliceai_csq_suffix_selects_row_by_consequence_symbol() {
        let schema = Arc::new(Schema::new(vec![
            Field::new("symbol", DataType::Utf8, true),
            Field::new("ds_ag", DataType::Utf8, true),
        ]));
        let spliceai_index = PluginIndex {
            kind: PluginSourceKind::SpliceAI,
            backend: PluginBackend::Parquet(ParquetPluginIndex {
                lookup: HashMap::from([(
                    (1029699_i64, Box::<str>::from("A"), Box::<str>::from("T")),
                    vec![0],
                )]),
                data: RecordBatch::try_new(
                    schema.clone(),
                    vec![
                        Arc::new(StringArray::from(vec!["C7orf50"])),
                        Arc::new(StringArray::from(vec!["0.12"])),
                    ],
                )
                .expect("spliceai batch"),
                value_schema: schema,
                dbnsfp_match_data: None,
                spliceai_match_data: Some(SpliceAiMatchData {
                    symbol: StringArray::from(vec!["C7orf50"]),
                }),
            }),
        };

        let plugins = ContigPlugins {
            active_plugins: ActivePlugins {
                configs: vec![PluginConfig {
                    kind: PluginKind::SpliceAI,
                    source_dirs: vec![],
                }],
            },
            indexes: vec![spliceai_index],
        };

        assert_eq!(
            plugins.csq_suffix_for_consequence(
                1029699,
                "A",
                "T",
                true,
                true,
                Some("C7orf50"),
                None,
                None,
            ),
            "C7orf50|0.12"
        );
        assert_eq!(
            plugins.csq_suffix_for_consequence(
                1029699,
                "A",
                "T",
                true,
                true,
                Some("CHL1"),
                None,
                None,
            ),
            "|"
        );
    }
}
