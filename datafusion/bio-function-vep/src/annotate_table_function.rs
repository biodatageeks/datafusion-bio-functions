//! Table function registration for consequence annotation.
//!
//! `annotate_vep()` is the high-level consequence annotation entrypoint.

use std::sync::Arc;

use datafusion::arrow::datatypes::Schema;
use datafusion::catalog::{CatalogProviderList, SchemaProvider, TableFunctionImpl};
use datafusion::common::{DataFusionError, Result, ScalarValue};
use datafusion::datasource::TableProvider;
use datafusion::logical_expr::Expr;
use datafusion::prelude::SessionContext;

use crate::annotate_provider::AnnotateProvider;
use crate::annotation_store::AnnotationBackend;

/// Table function implementing
/// `annotate_vep(vcf_table, cache_source, backend [, options_json])`.
pub struct AnnotateFunction {
    session: Arc<SessionContext>,
    /// Catalog list captured at registration time to avoid acquiring
    /// SessionState locks during `call()` (planning time).
    catalog_list: Arc<dyn CatalogProviderList>,
    /// Default catalog name captured at registration time.
    default_catalog: String,
    /// Default schema name captured at registration time.
    default_schema: String,
}

impl AnnotateFunction {
    pub fn new(session: Arc<SessionContext>) -> Self {
        // Capture catalog metadata at registration time (before any planner
        // locks exist) so resolve_schema can look up tables without touching
        // the SessionState RwLock.
        let state = session.state();
        let catalog_list = Arc::clone(state.catalog_list());
        let default_catalog = state.config_options().catalog.default_catalog.clone();
        let default_schema = state.config_options().catalog.default_schema.clone();
        Self {
            session,
            catalog_list,
            default_catalog,
            default_schema,
        }
    }
}

impl std::fmt::Debug for AnnotateFunction {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(f, "AnnotateFunction")
    }
}

impl TableFunctionImpl for AnnotateFunction {
    fn call(&self, args: &[Expr]) -> Result<Arc<dyn TableProvider>> {
        if args.len() < 3 {
            return Err(DataFusionError::Plan(
                "annotate_vep() requires at least 3 arguments: vcf_table, cache_source, backend"
                    .to_string(),
            ));
        }

        let vcf_table = extract_string_arg(&args[0], "vcf_table", "annotate_vep")?;
        let cache_source = extract_string_arg(&args[1], "cache_source", "annotate_vep")?;
        let backend_raw = extract_string_arg(&args[2], "backend", "annotate_vep")?;
        let backend = AnnotationBackend::parse(&backend_raw)?;

        let options_json = if args.len() > 3 {
            Some(extract_string_arg(
                &args[3],
                "options_json",
                "annotate_vep",
            )?)
        } else {
            None
        };

        let vcf_schema = resolve_schema_from_catalog(
            &*self.catalog_list,
            &self.default_catalog,
            &self.default_schema,
            &vcf_table,
        )?;

        Ok(Arc::new(AnnotateProvider::new(
            Arc::clone(&self.session),
            vcf_table,
            cache_source,
            backend,
            options_json,
            vcf_schema,
        )))
    }
}

/// Extract a string literal from an expression.
fn extract_string_arg(arg: &Expr, name: &str, fn_name: &str) -> Result<String> {
    match arg {
        Expr::Literal(ScalarValue::Utf8(Some(val)), _) => {
            if val.contains('`') {
                return Err(DataFusionError::Plan(format!(
                    "{fn_name}() {name} must not contain backtick characters, got: {val}"
                )));
            }
            Ok(val.clone())
        }
        other => Err(DataFusionError::Plan(format!(
            "{fn_name}() {name} must be a string literal, got: {other}"
        ))),
    }
}

/// Resolve the Arrow schema of a registered table using a pre-captured
/// `CatalogProviderList`, completely bypassing `SessionContext` and its
/// `SessionState` RwLock.
///
/// This is critical for vepyr / polars-bio integration: DataFusion's SQL
/// planner holds a **write lock** on `SessionState` while resolving table
/// functions. Any call to `session.state()`, `session.catalog()`, or
/// `session.table()` from within `TableFunctionImpl::call()` will deadlock
/// because they all acquire a read lock on the same RwLock.
///
/// By capturing the `CatalogProviderList` at registration time (before any
/// planner locks exist), we can look up tables without touching the session.
fn resolve_schema_from_catalog(
    catalog_list: &dyn CatalogProviderList,
    default_catalog: &str,
    default_schema: &str,
    table_name: &str,
) -> Result<Schema> {
    // Support bare names ("vcf"), schema-qualified ("public.vcf"), and
    // fully-qualified ("datafusion.public.vcf") table references.
    let parts: Vec<&str> = table_name.split('.').collect();
    let (cat_name, schema_name, bare_name) = match parts.len() {
        3 => (parts[0], parts[1], parts[2]),
        2 => (default_catalog, parts[0], parts[1]),
        _ => (default_catalog, default_schema, table_name),
    };

    let catalog = catalog_list
        .catalog(cat_name)
        .ok_or_else(|| DataFusionError::Plan(format!("Catalog '{cat_name}' not found")))?;
    let schema_provider = catalog.schema(schema_name).ok_or_else(|| {
        DataFusionError::Plan(format!(
            "Schema '{schema_name}' not found in catalog '{cat_name}'"
        ))
    })?;

    let table_provider = resolve_table_sync(&*schema_provider, bare_name)?;
    Ok(table_provider.schema().as_ref().clone())
}

/// Run `SchemaProvider::table()` synchronously, handling both tokio-context
/// and no-tokio-context cases.
fn resolve_table_sync(
    schema_provider: &dyn SchemaProvider,
    table_name: &str,
) -> Result<Arc<dyn TableProvider>> {
    let result = match tokio::runtime::Handle::try_current() {
        Ok(handle) => {
            tokio::task::block_in_place(|| handle.block_on(schema_provider.table(table_name)))
        }
        Err(_) => {
            let rt = tokio::runtime::Runtime::new()
                .map_err(|e| DataFusionError::External(Box::new(e)))?;
            rt.block_on(schema_provider.table(table_name))
        }
    };
    result
        .map_err(|e| DataFusionError::External(Box::new(e)))?
        .ok_or_else(|| DataFusionError::Plan(format!("Table '{table_name}' not found")))
}

#[cfg(test)]
mod tests {
    use crate::create_vep_session;
    #[cfg(feature = "kv-cache")]
    use crate::kv_cache::{VepKvStore, position_entry::serialize_position_entry};
    use crate::so_terms::SoTerm;
    use datafusion::arrow::array::{Array, Float64Array, Int64Array, RecordBatch, StringArray};
    use datafusion::arrow::datatypes::{DataType, Field, Schema};
    use datafusion::datasource::MemTable;
    use parquet::arrow::ArrowWriter;
    use std::collections::BTreeSet;
    use std::fs::File;
    use std::sync::Arc;
    use tempfile::TempDir;

    /// Write record batches into the partitioned cache layout under
    /// `{tmpdir}/{table_type}/{chrom}.parquet`, grouping rows by the "chrom"
    /// column in each batch.
    fn write_partitioned_cache(tmpdir: &TempDir, table_type: &str, batches: &[RecordBatch]) {
        let type_dir = tmpdir.path().join(table_type);
        std::fs::create_dir_all(&type_dir).expect("create type dir");

        // Collect all rows grouped by chrom value.
        let mut chrom_batches: std::collections::HashMap<String, Vec<RecordBatch>> =
            std::collections::HashMap::new();

        for batch in batches {
            let schema = batch.schema();
            let chrom_idx = schema
                .index_of("chrom")
                .expect("batch must have a 'chrom' column");
            let chrom_arr = batch
                .column(chrom_idx)
                .as_any()
                .downcast_ref::<StringArray>()
                .expect("chrom column must be StringArray");

            // Get unique chroms in this batch.
            let unique_chroms: BTreeSet<String> = (0..chrom_arr.len())
                .map(|i| chrom_arr.value(i).to_string())
                .collect();

            for chrom in unique_chroms {
                // Build boolean mask for this chrom.
                let indices: Vec<u32> = (0..chrom_arr.len())
                    .filter(|&i| chrom_arr.value(i) == chrom)
                    .map(|i| i as u32)
                    .collect();
                let indices_arr = datafusion::arrow::array::UInt32Array::from(indices);
                let filtered_columns: Vec<Arc<dyn Array>> = (0..batch.num_columns())
                    .map(|c| {
                        datafusion::arrow::compute::take(
                            batch.column(c).as_ref(),
                            &indices_arr,
                            None,
                        )
                        .expect("take should succeed")
                    })
                    .collect();
                let filtered_batch =
                    RecordBatch::try_new(schema.clone(), filtered_columns).expect("filtered batch");
                chrom_batches.entry(chrom).or_default().push(filtered_batch);
            }
        }

        // Write per-chrom parquet files.
        for (chrom, batches) in &chrom_batches {
            let path = type_dir.join(format!("{chrom}.parquet"));
            let file = File::create(&path).expect("create parquet file");
            let mut writer = ArrowWriter::try_new(file, batches[0].schema(), None)
                .expect("create parquet writer");
            for b in batches {
                writer.write(b).expect("write parquet batch");
            }
            writer.close().expect("close parquet writer");
        }
    }

    /// Convenience: write a single RecordBatch to the partitioned cache.
    fn write_batch_to_cache(tmpdir: &TempDir, table_type: &str, batch: &RecordBatch) {
        write_partitioned_cache(tmpdir, table_type, &[batch.clone()]);
    }

    /// Write a RecordBatch directly into a specific chrom parquet file,
    /// for tables that do not have a "chrom" column (e.g. exons, translations).
    fn write_batch_to_chrom(tmpdir: &TempDir, table_type: &str, chrom: &str, batch: &RecordBatch) {
        let type_dir = tmpdir.path().join(table_type);
        std::fs::create_dir_all(&type_dir).expect("create type dir");
        let path = type_dir.join(format!("{chrom}.parquet"));
        let file = File::create(&path).expect("create parquet file");
        let mut writer =
            ArrowWriter::try_new(file, batch.schema(), None).expect("create parquet writer");
        writer.write(batch).expect("write parquet batch");
        writer.close().expect("close parquet writer");
    }

    #[cfg(feature = "kv-cache")]
    fn write_batch_to_fjall(tmpdir: &TempDir, batch: &RecordBatch) {
        let fjall_dir = tmpdir.path().join("variation.fjall");
        let schema = batch.schema();
        let store = VepKvStore::create(&fjall_dir, schema.clone()).expect("create fjall cache");

        let chrom_idx = schema.index_of("chrom").expect("chrom column");
        let start_idx = schema.index_of("start").expect("start column");
        let allele_idx = schema
            .index_of("allele_string")
            .expect("allele_string column");
        let stored_col_indices: Vec<usize> = (0..schema.fields().len())
            .filter(|&idx| idx != chrom_idx && idx != start_idx)
            .collect();

        let chroms = batch
            .column(chrom_idx)
            .as_any()
            .downcast_ref::<StringArray>()
            .expect("chrom must be StringArray");
        let starts = batch
            .column(start_idx)
            .as_any()
            .downcast_ref::<Int64Array>()
            .expect("start must be Int64Array");

        for row in 0..batch.num_rows() {
            let entry = serialize_position_entry(&[row], batch, &stored_col_indices, allele_idx)
                .expect("serialize position entry");
            store
                .put_position_entry(chroms.value(row), starts.value(row), &entry)
                .expect("write fjall position entry");
        }

        store.persist().expect("persist fjall cache");
    }

    fn vcf_table() -> MemTable {
        let schema = Arc::new(Schema::new(vec![
            Field::new("chrom", DataType::Utf8, false),
            Field::new("start", DataType::Int64, false),
            Field::new("end", DataType::Int64, false),
            Field::new("ref", DataType::Utf8, false),
            Field::new("alt", DataType::Utf8, false),
        ]));
        let batch = RecordBatch::try_new(
            schema.clone(),
            vec![
                Arc::new(StringArray::from(vec!["1", "2"])),
                Arc::new(Int64Array::from(vec![100, 200])),
                Arc::new(Int64Array::from(vec![101, 201])),
                Arc::new(StringArray::from(vec!["A", "C"])),
                Arc::new(StringArray::from(vec!["G", "T"])),
            ],
        )
        .expect("valid vcf batch");
        MemTable::try_new(schema, vec![vec![batch]]).expect("valid vcf memtable")
    }

    fn cache_batch() -> RecordBatch {
        use crate::annotate_provider::cache_lookup_column_names;

        // Core columns required for the join.
        // Two rows: chrom "1" at pos 100 (matches vcf_table row 1) and
        // chrom "2" at pos 200 (matches vcf_table row 2).
        let mut fields = vec![
            Field::new("chrom", DataType::Utf8, false),
            Field::new("start", DataType::Int64, false),
            Field::new("end", DataType::Int64, false),
            Field::new("allele_string", DataType::Utf8, false),
            Field::new("failed", DataType::Int64, false),
        ];
        let mut columns: Vec<Arc<dyn datafusion::arrow::array::Array>> = vec![
            Arc::new(StringArray::from(vec!["1", "2"])),
            Arc::new(Int64Array::from(vec![100, 200])),
            Arc::new(Int64Array::from(vec![101, 201])),
            Arc::new(StringArray::from(vec!["A/G", "C/T"])),
            Arc::new(Int64Array::from(vec![0, 0])),
        ];
        // Add all cache lookup columns as nullable Utf8.
        for col in cache_lookup_column_names() {
            fields.push(Field::new(col, DataType::Utf8, true));
            let val: Option<&str> = match col {
                "variation_name" => Some("rs100"),
                "clin_sig" => Some("benign"),
                _ => None,
            };
            columns.push(Arc::new(StringArray::from(vec![val, None])));
        }
        let schema = Arc::new(Schema::new(fields));
        RecordBatch::try_new(schema, columns).expect("valid cache batch")
    }

    fn transcripts_batch() -> RecordBatch {
        let schema = Arc::new(Schema::new(vec![
            Field::new("transcript_id", DataType::Utf8, false),
            Field::new("chrom", DataType::Utf8, false),
            Field::new("start", DataType::Int64, false),
            Field::new("end", DataType::Int64, false),
            Field::new("strand", DataType::Int64, false),
            Field::new("biotype", DataType::Utf8, false),
            Field::new("cds_start", DataType::Int64, true),
            Field::new("cds_end", DataType::Int64, true),
        ]));
        RecordBatch::try_new(
            schema,
            vec![
                Arc::new(StringArray::from(vec![
                    "ENST00000000010",
                    "ENST00000000011",
                ])),
                Arc::new(StringArray::from(vec!["1", "2"])),
                Arc::new(Int64Array::from(vec![50, 150])),
                Arc::new(Int64Array::from(vec![200, 250])),
                Arc::new(Int64Array::from(vec![1, 1])),
                Arc::new(StringArray::from(vec!["protein_coding", "lincRNA"])),
                Arc::new(Int64Array::from(vec![80, 0])),
                Arc::new(Int64Array::from(vec![180, 0])),
            ],
        )
        .expect("valid transcript batch")
    }

    fn exons_batch() -> RecordBatch {
        let schema = Arc::new(Schema::new(vec![
            Field::new("transcript_id", DataType::Utf8, false),
            Field::new("exon_number", DataType::Int64, false),
            Field::new("start", DataType::Int64, false),
            Field::new("end", DataType::Int64, false),
        ]));
        RecordBatch::try_new(
            schema,
            vec![
                Arc::new(StringArray::from(vec![
                    "ENST00000000010",
                    "ENST00000000011",
                ])),
                Arc::new(Int64Array::from(vec![1, 1])),
                Arc::new(Int64Array::from(vec![50, 150])),
                Arc::new(Int64Array::from(vec![200, 250])),
            ],
        )
        .expect("valid exon batch")
    }

    fn refseq_vcf_table() -> MemTable {
        let schema = Arc::new(Schema::new(vec![
            Field::new("chrom", DataType::Utf8, false),
            Field::new("start", DataType::Int64, false),
            Field::new("end", DataType::Int64, false),
            Field::new("ref", DataType::Utf8, false),
            Field::new("alt", DataType::Utf8, false),
        ]));
        let batch = RecordBatch::try_new(
            schema.clone(),
            vec![
                Arc::new(StringArray::from(vec!["1"])),
                Arc::new(Int64Array::from(vec![100])),
                Arc::new(Int64Array::from(vec![101])),
                Arc::new(StringArray::from(vec!["A"])),
                Arc::new(StringArray::from(vec!["G"])),
            ],
        )
        .expect("valid refseq vcf batch");
        MemTable::try_new(schema, vec![vec![batch]]).expect("valid refseq vcf memtable")
    }

    fn refseq_cache_batch() -> RecordBatch {
        use crate::annotate_provider::cache_lookup_column_names;

        let mut fields = vec![
            Field::new("chrom", DataType::Utf8, false),
            Field::new("start", DataType::Int64, false),
            Field::new("end", DataType::Int64, false),
            Field::new("allele_string", DataType::Utf8, false),
            Field::new("failed", DataType::Int64, false),
        ];
        let mut columns: Vec<Arc<dyn datafusion::arrow::array::Array>> = vec![
            Arc::new(StringArray::from(vec!["1"])),
            Arc::new(Int64Array::from(vec![100])),
            Arc::new(Int64Array::from(vec![101])),
            Arc::new(StringArray::from(vec!["A/G"])),
            Arc::new(Int64Array::from(vec![0_i64])),
        ];
        for col in cache_lookup_column_names() {
            fields.push(Field::new(col, DataType::Utf8, true));
            columns.push(Arc::new(StringArray::from(vec![None::<&str>])));
        }
        let schema = Arc::new(Schema::new(fields));
        RecordBatch::try_new(schema, columns).expect("valid refseq cache batch")
    }

    fn refseq_transcripts_batch() -> RecordBatch {
        let schema = Arc::new(Schema::new(vec![
            Field::new("transcript_id", DataType::Utf8, false),
            Field::new("chrom", DataType::Utf8, false),
            Field::new("start", DataType::Int64, false),
            Field::new("end", DataType::Int64, false),
            Field::new("strand", DataType::Int64, false),
            Field::new("biotype", DataType::Utf8, false),
            Field::new("cds_start", DataType::Int64, true),
            Field::new("cds_end", DataType::Int64, true),
            Field::new("source", DataType::Utf8, true),
            Field::new("bam_edit_status", DataType::Utf8, true),
            Field::new("raw_object_json", DataType::Utf8, true),
        ]));
        let raw = r#"{"__class":"Bio::EnsEMBL::Transcript","__value":{"_source_cache":"RefSeq","display_xref":{"display_id":"NM_000001"},"attributes":[{"__class":"Bio::EnsEMBL::Attribute","__value":{"code":"enst_refseq_compare","value":"ENST00000332831:cds_only"}},{"__class":"Bio::EnsEMBL::Attribute","__value":{"code":"rseq_ens_match_cds","value":"1"}}]}}"#;
        RecordBatch::try_new(
            schema,
            vec![
                Arc::new(StringArray::from(vec!["NM_000001"])),
                Arc::new(StringArray::from(vec!["1"])),
                Arc::new(Int64Array::from(vec![50])),
                Arc::new(Int64Array::from(vec![200])),
                Arc::new(Int64Array::from(vec![1])),
                Arc::new(StringArray::from(vec!["lncRNA"])),
                Arc::new(Int64Array::from(vec![None::<i64>])),
                Arc::new(Int64Array::from(vec![None::<i64>])),
                Arc::new(StringArray::from(vec![Some("BestRefSeq")])),
                Arc::new(StringArray::from(vec![Some("ok")])),
                Arc::new(StringArray::from(vec![Some(raw)])),
            ],
        )
        .expect("valid refseq transcript batch")
    }

    fn refseq_exons_batch() -> RecordBatch {
        let schema = Arc::new(Schema::new(vec![
            Field::new("transcript_id", DataType::Utf8, false),
            Field::new("exon_number", DataType::Int64, false),
            Field::new("start", DataType::Int64, false),
            Field::new("end", DataType::Int64, false),
        ]));
        RecordBatch::try_new(
            schema,
            vec![
                Arc::new(StringArray::from(vec!["NM_000001"])),
                Arc::new(Int64Array::from(vec![1])),
                Arc::new(Int64Array::from(vec![50])),
                Arc::new(Int64Array::from(vec![200])),
            ],
        )
        .expect("valid refseq exon batch")
    }

    fn synonymous_vcf_table() -> MemTable {
        let schema = Arc::new(Schema::new(vec![
            Field::new("chrom", DataType::Utf8, false),
            Field::new("start", DataType::Int64, false),
            Field::new("end", DataType::Int64, false),
            Field::new("ref", DataType::Utf8, false),
            Field::new("alt", DataType::Utf8, false),
        ]));
        let batch = RecordBatch::try_new(
            schema.clone(),
            vec![
                Arc::new(StringArray::from(vec!["1"])),
                Arc::new(Int64Array::from(vec![100])),
                Arc::new(Int64Array::from(vec![102])),
                Arc::new(StringArray::from(vec!["GCT"])),
                Arc::new(StringArray::from(vec!["GCC"])),
            ],
        )
        .expect("valid synonymous vcf batch");
        MemTable::try_new(schema, vec![vec![batch]]).expect("valid synonymous vcf memtable")
    }

    fn synonymous_cache_batch() -> RecordBatch {
        let schema = Arc::new(Schema::new(vec![
            Field::new("chrom", DataType::Utf8, false),
            Field::new("start", DataType::Int64, false),
            Field::new("end", DataType::Int64, false),
            Field::new("variation_name", DataType::Utf8, true),
            Field::new("allele_string", DataType::Utf8, false),
            Field::new("clin_sig", DataType::Utf8, true),
            Field::new("AF", DataType::Float64, true),
            Field::new("failed", DataType::Int64, false),
        ]));
        RecordBatch::try_new(
            schema,
            vec![
                Arc::new(StringArray::from(vec!["1"])),
                Arc::new(Int64Array::from(vec![100])),
                Arc::new(Int64Array::from(vec![102])),
                Arc::new(StringArray::from(vec!["rs_syn"])),
                Arc::new(StringArray::from(vec!["GCT/GCC"])),
                Arc::new(StringArray::from(vec!["benign"])),
                Arc::new(Float64Array::from(vec![0.01_f64])),
                Arc::new(Int64Array::from(vec![0])),
            ],
        )
        .expect("valid synonymous cache batch")
    }

    fn synonymous_transcripts_batch() -> RecordBatch {
        let schema = Arc::new(Schema::new(vec![
            Field::new("transcript_id", DataType::Utf8, false),
            Field::new("chrom", DataType::Utf8, false),
            Field::new("start", DataType::Int64, false),
            Field::new("end", DataType::Int64, false),
            Field::new("strand", DataType::Int64, false),
            Field::new("biotype", DataType::Utf8, false),
            Field::new("cds_start", DataType::Int64, true),
            Field::new("cds_end", DataType::Int64, true),
        ]));
        RecordBatch::try_new(
            schema,
            vec![
                Arc::new(StringArray::from(vec!["ENST00000000001"])),
                Arc::new(StringArray::from(vec!["1"])),
                Arc::new(Int64Array::from(vec![90])),
                Arc::new(Int64Array::from(vec![140])),
                Arc::new(Int64Array::from(vec![1])),
                Arc::new(StringArray::from(vec!["protein_coding"])),
                Arc::new(Int64Array::from(vec![97])),
                Arc::new(Int64Array::from(vec![105])),
            ],
        )
        .expect("valid synonymous transcript batch")
    }

    fn synonymous_exons_batch() -> RecordBatch {
        let schema = Arc::new(Schema::new(vec![
            Field::new("transcript_id", DataType::Utf8, false),
            Field::new("exon_number", DataType::Int64, false),
            Field::new("start", DataType::Int64, false),
            Field::new("end", DataType::Int64, false),
        ]));
        RecordBatch::try_new(
            schema,
            vec![
                Arc::new(StringArray::from(vec!["ENST00000000001"])),
                Arc::new(Int64Array::from(vec![1])),
                Arc::new(Int64Array::from(vec![90])),
                Arc::new(Int64Array::from(vec![140])),
            ],
        )
        .expect("valid synonymous exon batch")
    }

    fn synonymous_translations_batch() -> RecordBatch {
        let schema = Arc::new(Schema::new(vec![
            Field::new("transcript_id", DataType::Utf8, false),
            Field::new("cds_len", DataType::Int64, true),
            Field::new("protein_len", DataType::Int64, true),
            Field::new("translation_seq", DataType::Utf8, true),
            Field::new("cds_sequence", DataType::Utf8, true),
        ]));
        // CDS from 97..105 => ATG GCT TAA
        RecordBatch::try_new(
            schema,
            vec![
                Arc::new(StringArray::from(vec!["ENST00000000001"])),
                Arc::new(Int64Array::from(vec![9])),
                Arc::new(Int64Array::from(vec![3])),
                Arc::new(StringArray::from(vec!["MA*"])),
                Arc::new(StringArray::from(vec!["ATGGCTTAA"])),
            ],
        )
        .expect("valid synonymous translations batch")
    }

    fn context_vcf_table() -> MemTable {
        let schema = Arc::new(Schema::new(vec![
            Field::new("chrom", DataType::Utf8, false),
            Field::new("start", DataType::Int64, false),
            Field::new("end", DataType::Int64, false),
            Field::new("ref", DataType::Utf8, false),
            Field::new("alt", DataType::Utf8, false),
        ]));
        let batch = RecordBatch::try_new(
            schema.clone(),
            vec![
                Arc::new(StringArray::from(vec!["1"])),
                Arc::new(Int64Array::from(vec![155])),
                Arc::new(Int64Array::from(vec![155])),
                Arc::new(StringArray::from(vec!["A"])),
                Arc::new(StringArray::from(vec!["G"])),
            ],
        )
        .expect("valid context vcf batch");
        MemTable::try_new(schema, vec![vec![batch]]).expect("valid context vcf memtable")
    }

    fn context_cache_batch() -> RecordBatch {
        let schema = Arc::new(Schema::new(vec![
            Field::new("chrom", DataType::Utf8, false),
            Field::new("start", DataType::Int64, false),
            Field::new("end", DataType::Int64, false),
            Field::new("variation_name", DataType::Utf8, true),
            Field::new("allele_string", DataType::Utf8, false),
            Field::new("clin_sig", DataType::Utf8, true),
            Field::new("AF", DataType::Float64, true),
            Field::new("failed", DataType::Int64, false),
        ]));
        RecordBatch::try_new(
            schema,
            vec![
                Arc::new(StringArray::from(vec!["1"])),
                Arc::new(Int64Array::from(vec![155])),
                Arc::new(Int64Array::from(vec![155])),
                Arc::new(StringArray::from(vec!["rs_ctx"])),
                Arc::new(StringArray::from(vec!["A/G"])),
                Arc::new(StringArray::from(vec!["pathogenic"])),
                Arc::new(Float64Array::from(vec![0.42_f64])),
                Arc::new(Int64Array::from(vec![0])),
            ],
        )
        .expect("valid context cache batch")
    }

    fn context_transcripts_batch() -> RecordBatch {
        let schema = Arc::new(Schema::new(vec![
            Field::new("transcript_id", DataType::Utf8, false),
            Field::new("chrom", DataType::Utf8, false),
            Field::new("start", DataType::Int64, false),
            Field::new("end", DataType::Int64, false),
            Field::new("strand", DataType::Int64, false),
            Field::new("biotype", DataType::Utf8, false),
            Field::new("cds_start", DataType::Int64, true),
            Field::new("cds_end", DataType::Int64, true),
        ]));
        RecordBatch::try_new(
            schema,
            vec![
                Arc::new(StringArray::from(vec!["ENST00000000002"])),
                Arc::new(StringArray::from(vec!["1"])),
                Arc::new(Int64Array::from(vec![100])),
                Arc::new(Int64Array::from(vec![250])),
                Arc::new(Int64Array::from(vec![1])),
                Arc::new(StringArray::from(vec!["protein_coding"])),
                Arc::new(Int64Array::from(vec![120])),
                Arc::new(Int64Array::from(vec![240])),
            ],
        )
        .expect("valid context transcript batch")
    }

    fn context_exons_batch() -> RecordBatch {
        let schema = Arc::new(Schema::new(vec![
            Field::new("transcript_id", DataType::Utf8, false),
            Field::new("exon_number", DataType::Int64, false),
            Field::new("start", DataType::Int64, false),
            Field::new("end", DataType::Int64, false),
        ]));
        RecordBatch::try_new(
            schema,
            vec![
                Arc::new(StringArray::from(vec!["ENST00000000002"])),
                Arc::new(Int64Array::from(vec![1])),
                Arc::new(Int64Array::from(vec![100])),
                Arc::new(Int64Array::from(vec![250])),
            ],
        )
        .expect("valid context exon batch")
    }

    fn context_regulatory_batch() -> RecordBatch {
        let schema = Arc::new(Schema::new(vec![
            Field::new("stable_id", DataType::Utf8, false),
            Field::new("chrom", DataType::Utf8, false),
            Field::new("start", DataType::Int64, false),
            Field::new("end", DataType::Int64, false),
        ]));
        RecordBatch::try_new(
            schema,
            vec![
                Arc::new(StringArray::from(vec!["reg_ctx"])),
                Arc::new(StringArray::from(vec!["1"])),
                Arc::new(Int64Array::from(vec![150])),
                Arc::new(Int64Array::from(vec![160])),
            ],
        )
        .expect("valid context regulatory batch")
    }

    fn context_motif_batch() -> RecordBatch {
        let schema = Arc::new(Schema::new(vec![
            Field::new("motif_id", DataType::Utf8, false),
            Field::new("chrom", DataType::Utf8, false),
            Field::new("start", DataType::Int64, false),
            Field::new("end", DataType::Int64, false),
        ]));
        RecordBatch::try_new(
            schema,
            vec![
                Arc::new(StringArray::from(vec!["motif_ctx"])),
                Arc::new(StringArray::from(vec!["1"])),
                Arc::new(Int64Array::from(vec![150])),
                Arc::new(Int64Array::from(vec![160])),
            ],
        )
        .expect("valid context motif batch")
    }

    fn golden_context_vcf_table() -> MemTable {
        let schema = Arc::new(Schema::new(vec![
            Field::new("chrom", DataType::Utf8, false),
            Field::new("start", DataType::Int64, false),
            Field::new("end", DataType::Int64, false),
            Field::new("ref", DataType::Utf8, false),
            Field::new("alt", DataType::Utf8, false),
        ]));
        let batch = RecordBatch::try_new(
            schema.clone(),
            vec![
                Arc::new(StringArray::from(vec!["1"])),
                Arc::new(Int64Array::from(vec![500])),
                Arc::new(Int64Array::from(vec![500])),
                Arc::new(StringArray::from(vec!["A"])),
                Arc::new(StringArray::from(vec!["G"])),
            ],
        )
        .expect("valid golden context vcf batch");
        MemTable::try_new(schema, vec![vec![batch]]).expect("valid golden context vcf memtable")
    }

    fn golden_context_cache_batch() -> RecordBatch {
        let schema = Arc::new(Schema::new(vec![
            Field::new("chrom", DataType::Utf8, false),
            Field::new("start", DataType::Int64, false),
            Field::new("end", DataType::Int64, false),
            Field::new("variation_name", DataType::Utf8, true),
            Field::new("allele_string", DataType::Utf8, false),
            Field::new("clin_sig", DataType::Utf8, true),
            Field::new("AF", DataType::Float64, true),
            Field::new("failed", DataType::Int64, false),
        ]));
        RecordBatch::try_new(
            schema,
            vec![
                Arc::new(StringArray::from(vec!["1"])),
                Arc::new(Int64Array::from(vec![500])),
                Arc::new(Int64Array::from(vec![500])),
                Arc::new(StringArray::from(vec!["rs_golden_ctx"])),
                Arc::new(StringArray::from(vec!["A/G"])),
                Arc::new(StringArray::from(vec!["pathogenic"])),
                Arc::new(Float64Array::from(vec![0.11_f64])),
                Arc::new(Int64Array::from(vec![0])),
            ],
        )
        .expect("valid golden context cache batch")
    }

    fn golden_context_transcripts_batch() -> RecordBatch {
        let schema = Arc::new(Schema::new(vec![
            Field::new("transcript_id", DataType::Utf8, false),
            Field::new("chrom", DataType::Utf8, false),
            Field::new("start", DataType::Int64, false),
            Field::new("end", DataType::Int64, false),
            Field::new("strand", DataType::Int64, false),
            Field::new("biotype", DataType::Utf8, false),
            Field::new("cds_start", DataType::Int64, true),
            Field::new("cds_end", DataType::Int64, true),
        ]));
        RecordBatch::try_new(
            schema,
            vec![
                Arc::new(StringArray::from(vec!["ENST00000000003"])),
                Arc::new(StringArray::from(vec!["1"])),
                Arc::new(Int64Array::from(vec![20_000])),
                Arc::new(Int64Array::from(vec![20_100])),
                Arc::new(Int64Array::from(vec![1])),
                Arc::new(StringArray::from(vec!["protein_coding"])),
                Arc::new(Int64Array::from(vec![20_020])),
                Arc::new(Int64Array::from(vec![20_080])),
            ],
        )
        .expect("valid golden context transcript batch")
    }

    fn golden_context_exons_batch() -> RecordBatch {
        let schema = Arc::new(Schema::new(vec![
            Field::new("transcript_id", DataType::Utf8, false),
            Field::new("exon_number", DataType::Int64, false),
            Field::new("start", DataType::Int64, false),
            Field::new("end", DataType::Int64, false),
        ]));
        RecordBatch::try_new(
            schema,
            vec![
                Arc::new(StringArray::from(vec!["ENST00000000003"])),
                Arc::new(Int64Array::from(vec![1])),
                Arc::new(Int64Array::from(vec![20_000])),
                Arc::new(Int64Array::from(vec![20_100])),
            ],
        )
        .expect("valid golden context exon batch")
    }

    fn golden_context_regulatory_batch() -> RecordBatch {
        let schema = Arc::new(Schema::new(vec![
            Field::new("stable_id", DataType::Utf8, false),
            Field::new("chrom", DataType::Utf8, false),
            Field::new("start", DataType::Int64, false),
            Field::new("end", DataType::Int64, false),
        ]));
        RecordBatch::try_new(
            schema,
            vec![
                Arc::new(StringArray::from(vec!["reg_golden"])),
                Arc::new(StringArray::from(vec!["1"])),
                Arc::new(Int64Array::from(vec![495])),
                Arc::new(Int64Array::from(vec![505])),
            ],
        )
        .expect("valid golden context regulatory batch")
    }

    fn golden_context_motif_batch() -> RecordBatch {
        let schema = Arc::new(Schema::new(vec![
            Field::new("motif_id", DataType::Utf8, false),
            Field::new("chrom", DataType::Utf8, false),
            Field::new("start", DataType::Int64, false),
            Field::new("end", DataType::Int64, false),
        ]));
        RecordBatch::try_new(
            schema,
            vec![
                Arc::new(StringArray::from(vec!["motif_golden"])),
                Arc::new(StringArray::from(vec!["1"])),
                Arc::new(Int64Array::from(vec![495])),
                Arc::new(Int64Array::from(vec![505])),
            ],
        )
        .expect("valid golden context motif batch")
    }

    fn splice_vcf_table() -> MemTable {
        let schema = Arc::new(Schema::new(vec![
            Field::new("chrom", DataType::Utf8, false),
            Field::new("start", DataType::Int64, false),
            Field::new("end", DataType::Int64, false),
            Field::new("ref", DataType::Utf8, false),
            Field::new("alt", DataType::Utf8, false),
        ]));
        let batch = RecordBatch::try_new(
            schema.clone(),
            vec![
                Arc::new(StringArray::from(vec!["1"])),
                Arc::new(Int64Array::from(vec![151])),
                Arc::new(Int64Array::from(vec![151])),
                Arc::new(StringArray::from(vec!["A"])),
                Arc::new(StringArray::from(vec!["G"])),
            ],
        )
        .expect("valid splice vcf batch");
        MemTable::try_new(schema, vec![vec![batch]]).expect("valid splice vcf memtable")
    }

    fn splice_cache_batch() -> RecordBatch {
        let schema = Arc::new(Schema::new(vec![
            Field::new("chrom", DataType::Utf8, false),
            Field::new("start", DataType::Int64, false),
            Field::new("end", DataType::Int64, false),
            Field::new("variation_name", DataType::Utf8, true),
            Field::new("allele_string", DataType::Utf8, false),
            Field::new("clin_sig", DataType::Utf8, true),
            Field::new("AF", DataType::Float64, true),
            Field::new("failed", DataType::Int64, false),
        ]));
        RecordBatch::try_new(
            schema,
            vec![
                Arc::new(StringArray::from(vec!["1"])),
                Arc::new(Int64Array::from(vec![151])),
                Arc::new(Int64Array::from(vec![151])),
                Arc::new(StringArray::from(vec!["rs_splice"])),
                Arc::new(StringArray::from(vec!["A/G"])),
                Arc::new(StringArray::from(vec!["uncertain_significance"])),
                Arc::new(Float64Array::from(vec![0.2_f64])),
                Arc::new(Int64Array::from(vec![0])),
            ],
        )
        .expect("valid splice cache batch")
    }

    fn splice_transcripts_batch() -> RecordBatch {
        let schema = Arc::new(Schema::new(vec![
            Field::new("transcript_id", DataType::Utf8, false),
            Field::new("chrom", DataType::Utf8, false),
            Field::new("start", DataType::Int64, false),
            Field::new("end", DataType::Int64, false),
            Field::new("strand", DataType::Int64, false),
            Field::new("biotype", DataType::Utf8, false),
            Field::new("cds_start", DataType::Int64, true),
            Field::new("cds_end", DataType::Int64, true),
        ]));
        RecordBatch::try_new(
            schema,
            vec![
                Arc::new(StringArray::from(vec!["ENST00000000004"])),
                Arc::new(StringArray::from(vec!["1"])),
                Arc::new(Int64Array::from(vec![100])),
                Arc::new(Int64Array::from(vec![300])),
                Arc::new(Int64Array::from(vec![1])),
                Arc::new(StringArray::from(vec!["protein_coding"])),
                Arc::new(Int64Array::from(vec![120])),
                Arc::new(Int64Array::from(vec![280])),
            ],
        )
        .expect("valid splice transcript batch")
    }

    fn splice_exons_batch() -> RecordBatch {
        let schema = Arc::new(Schema::new(vec![
            Field::new("transcript_id", DataType::Utf8, false),
            Field::new("exon_number", DataType::Int64, false),
            Field::new("start", DataType::Int64, false),
            Field::new("end", DataType::Int64, false),
        ]));
        RecordBatch::try_new(
            schema,
            vec![
                Arc::new(StringArray::from(vec![
                    "ENST00000000004",
                    "ENST00000000004",
                ])),
                Arc::new(Int64Array::from(vec![1, 2])),
                Arc::new(Int64Array::from(vec![100, 250])),
                Arc::new(Int64Array::from(vec![150, 300])),
            ],
        )
        .expect("valid splice exon batch")
    }

    fn repeat_shift_vcf_table() -> MemTable {
        let schema = Arc::new(Schema::new(vec![
            Field::new("chrom", DataType::Utf8, false),
            Field::new("start", DataType::Int64, false),
            Field::new("end", DataType::Int64, false),
            Field::new("ref", DataType::Utf8, false),
            Field::new("alt", DataType::Utf8, false),
        ]));
        let batch = RecordBatch::try_new(
            schema.clone(),
            vec![
                Arc::new(StringArray::from(vec!["1"])),
                Arc::new(Int64Array::from(vec![55])),
                Arc::new(Int64Array::from(vec![70])),
                Arc::new(StringArray::from(vec!["GAAGAAGAAGAAGAA"])),
                Arc::new(StringArray::from(vec!["G"])),
            ],
        )
        .expect("valid repeat-shift vcf batch");
        MemTable::try_new(schema, vec![vec![batch]]).expect("valid repeat-shift vcf memtable")
    }

    fn repeat_shift_cache_batch() -> RecordBatch {
        let schema = Arc::new(Schema::new(vec![
            Field::new("chrom", DataType::Utf8, false),
            Field::new("start", DataType::Int64, false),
            Field::new("end", DataType::Int64, false),
            Field::new("variation_name", DataType::Utf8, true),
            Field::new("allele_string", DataType::Utf8, false),
            Field::new("clin_sig", DataType::Utf8, true),
            Field::new("AF", DataType::Float64, true),
            Field::new("failed", DataType::Int64, false),
        ]));
        RecordBatch::try_new(
            schema,
            vec![
                Arc::new(StringArray::from(vec!["1"])),
                Arc::new(Int64Array::from(vec![66])),
                Arc::new(Int64Array::from(vec![79])),
                Arc::new(StringArray::from(vec!["rs_repeat_shift"])),
                Arc::new(StringArray::from(vec!["AAGAAGAAGAAGAA/-"])),
                Arc::new(StringArray::from(vec!["likely_pathogenic"])),
                Arc::new(Float64Array::from(vec![0.005_f64])),
                Arc::new(Int64Array::from(vec![0])),
            ],
        )
        .expect("valid repeat-shift cache batch")
    }

    fn repeat_shift_transcripts_batch() -> RecordBatch {
        let schema = Arc::new(Schema::new(vec![
            Field::new("transcript_id", DataType::Utf8, false),
            Field::new("chrom", DataType::Utf8, false),
            Field::new("start", DataType::Int64, false),
            Field::new("end", DataType::Int64, false),
            Field::new("strand", DataType::Int64, false),
            Field::new("biotype", DataType::Utf8, false),
            Field::new("cds_start", DataType::Int64, true),
            Field::new("cds_end", DataType::Int64, true),
        ]));
        RecordBatch::try_new(
            schema,
            vec![
                Arc::new(StringArray::from(vec!["ENST00000000005"])),
                Arc::new(StringArray::from(vec!["1"])),
                Arc::new(Int64Array::from(vec![40])),
                Arc::new(Int64Array::from(vec![120])),
                Arc::new(Int64Array::from(vec![1])),
                Arc::new(StringArray::from(vec!["protein_coding"])),
                Arc::new(Int64Array::from(vec![45])),
                Arc::new(Int64Array::from(vec![110])),
            ],
        )
        .expect("valid repeat-shift transcript batch")
    }

    fn repeat_shift_exons_batch() -> RecordBatch {
        let schema = Arc::new(Schema::new(vec![
            Field::new("transcript_id", DataType::Utf8, false),
            Field::new("exon_number", DataType::Int64, false),
            Field::new("start", DataType::Int64, false),
            Field::new("end", DataType::Int64, false),
        ]));
        RecordBatch::try_new(
            schema,
            vec![
                Arc::new(StringArray::from(vec!["ENST00000000005"])),
                Arc::new(Int64Array::from(vec![1])),
                Arc::new(Int64Array::from(vec![40])),
                Arc::new(Int64Array::from(vec![120])),
            ],
        )
        .expect("valid repeat-shift exon batch")
    }

    fn neg_strand_vcf_table() -> MemTable {
        let schema = Arc::new(Schema::new(vec![
            Field::new("chrom", DataType::Utf8, false),
            Field::new("start", DataType::Int64, false),
            Field::new("end", DataType::Int64, false),
            Field::new("ref", DataType::Utf8, false),
            Field::new("alt", DataType::Utf8, false),
        ]));
        let batch = RecordBatch::try_new(
            schema.clone(),
            vec![
                Arc::new(StringArray::from(vec!["1"])),
                Arc::new(Int64Array::from(vec![103])),
                Arc::new(Int64Array::from(vec![103])),
                Arc::new(StringArray::from(vec!["A"])),
                Arc::new(StringArray::from(vec!["G"])),
            ],
        )
        .expect("valid negative-strand vcf batch");
        MemTable::try_new(schema, vec![vec![batch]]).expect("valid negative-strand vcf memtable")
    }

    fn neg_strand_cache_batch() -> RecordBatch {
        let schema = Arc::new(Schema::new(vec![
            Field::new("chrom", DataType::Utf8, false),
            Field::new("start", DataType::Int64, false),
            Field::new("end", DataType::Int64, false),
            Field::new("variation_name", DataType::Utf8, true),
            Field::new("allele_string", DataType::Utf8, false),
            Field::new("clin_sig", DataType::Utf8, true),
            Field::new("AF", DataType::Float64, true),
            Field::new("failed", DataType::Int64, false),
        ]));
        RecordBatch::try_new(
            schema,
            vec![
                Arc::new(StringArray::from(vec!["1"])),
                Arc::new(Int64Array::from(vec![103])),
                Arc::new(Int64Array::from(vec![103])),
                Arc::new(StringArray::from(vec!["rs_neg_syn"])),
                Arc::new(StringArray::from(vec!["A/G"])),
                Arc::new(StringArray::from(vec!["benign"])),
                Arc::new(Float64Array::from(vec![0.03_f64])),
                Arc::new(Int64Array::from(vec![0])),
            ],
        )
        .expect("valid negative-strand cache batch")
    }

    fn neg_strand_transcripts_batch() -> RecordBatch {
        let schema = Arc::new(Schema::new(vec![
            Field::new("transcript_id", DataType::Utf8, false),
            Field::new("chrom", DataType::Utf8, false),
            Field::new("start", DataType::Int64, false),
            Field::new("end", DataType::Int64, false),
            Field::new("strand", DataType::Int64, false),
            Field::new("biotype", DataType::Utf8, false),
            Field::new("cds_start", DataType::Int64, true),
            Field::new("cds_end", DataType::Int64, true),
        ]));
        RecordBatch::try_new(
            schema,
            vec![
                Arc::new(StringArray::from(vec!["ENST00000000006"])),
                Arc::new(StringArray::from(vec!["1"])),
                Arc::new(Int64Array::from(vec![90])),
                Arc::new(Int64Array::from(vec![140])),
                Arc::new(Int64Array::from(vec![-1])),
                Arc::new(StringArray::from(vec!["protein_coding"])),
                Arc::new(Int64Array::from(vec![100])),
                Arc::new(Int64Array::from(vec![108])),
            ],
        )
        .expect("valid negative-strand transcript batch")
    }

    fn neg_strand_exons_batch() -> RecordBatch {
        let schema = Arc::new(Schema::new(vec![
            Field::new("transcript_id", DataType::Utf8, false),
            Field::new("exon_number", DataType::Int64, false),
            Field::new("start", DataType::Int64, false),
            Field::new("end", DataType::Int64, false),
        ]));
        RecordBatch::try_new(
            schema,
            vec![
                Arc::new(StringArray::from(vec!["ENST00000000006"])),
                Arc::new(Int64Array::from(vec![1])),
                Arc::new(Int64Array::from(vec![90])),
                Arc::new(Int64Array::from(vec![140])),
            ],
        )
        .expect("valid negative-strand exon batch")
    }

    fn neg_strand_translations_batch() -> RecordBatch {
        let schema = Arc::new(Schema::new(vec![
            Field::new("transcript_id", DataType::Utf8, false),
            Field::new("cds_len", DataType::Int64, true),
            Field::new("protein_len", DataType::Int64, true),
            Field::new("translation_seq", DataType::Utf8, true),
            Field::new("cds_sequence", DataType::Utf8, true),
        ]));
        RecordBatch::try_new(
            schema,
            vec![
                Arc::new(StringArray::from(vec!["ENST00000000006"])),
                Arc::new(Int64Array::from(vec![9])),
                Arc::new(Int64Array::from(vec![3])),
                Arc::new(StringArray::from(vec!["MA*"])),
                Arc::new(StringArray::from(vec!["ATGGCTTAA"])),
            ],
        )
        .expect("valid negative-strand translation batch")
    }

    fn stop_loss_vcf_table() -> MemTable {
        let schema = Arc::new(Schema::new(vec![
            Field::new("chrom", DataType::Utf8, false),
            Field::new("start", DataType::Int64, false),
            Field::new("end", DataType::Int64, false),
            Field::new("ref", DataType::Utf8, false),
            Field::new("alt", DataType::Utf8, false),
        ]));
        let batch = RecordBatch::try_new(
            schema.clone(),
            vec![
                Arc::new(StringArray::from(vec!["1"])),
                Arc::new(Int64Array::from(vec![106])),
                Arc::new(Int64Array::from(vec![108])),
                Arc::new(StringArray::from(vec!["TAA"])),
                Arc::new(StringArray::from(vec!["-"])),
            ],
        )
        .expect("valid stop-loss vcf batch");
        MemTable::try_new(schema, vec![vec![batch]]).expect("valid stop-loss vcf memtable")
    }

    fn stop_loss_cache_batch() -> RecordBatch {
        let schema = Arc::new(Schema::new(vec![
            Field::new("chrom", DataType::Utf8, false),
            Field::new("start", DataType::Int64, false),
            Field::new("end", DataType::Int64, false),
            Field::new("variation_name", DataType::Utf8, true),
            Field::new("allele_string", DataType::Utf8, false),
            Field::new("clin_sig", DataType::Utf8, true),
            Field::new("AF", DataType::Float64, true),
            Field::new("failed", DataType::Int64, false),
        ]));
        RecordBatch::try_new(
            schema,
            vec![
                Arc::new(StringArray::from(vec!["1"])),
                Arc::new(Int64Array::from(vec![106])),
                Arc::new(Int64Array::from(vec![108])),
                Arc::new(StringArray::from(vec!["rs_stop_loss"])),
                Arc::new(StringArray::from(vec!["TAA/-"])),
                Arc::new(StringArray::from(vec!["pathogenic"])),
                Arc::new(Float64Array::from(vec![0.001_f64])),
                Arc::new(Int64Array::from(vec![0])),
            ],
        )
        .expect("valid stop-loss cache batch")
    }

    fn stop_loss_transcripts_batch() -> RecordBatch {
        let schema = Arc::new(Schema::new(vec![
            Field::new("transcript_id", DataType::Utf8, false),
            Field::new("chrom", DataType::Utf8, false),
            Field::new("start", DataType::Int64, false),
            Field::new("end", DataType::Int64, false),
            Field::new("strand", DataType::Int64, false),
            Field::new("biotype", DataType::Utf8, false),
            Field::new("cds_start", DataType::Int64, true),
            Field::new("cds_end", DataType::Int64, true),
        ]));
        RecordBatch::try_new(
            schema,
            vec![
                Arc::new(StringArray::from(vec!["ENST00000000007"])),
                Arc::new(StringArray::from(vec!["1"])),
                Arc::new(Int64Array::from(vec![90])),
                Arc::new(Int64Array::from(vec![140])),
                Arc::new(Int64Array::from(vec![1])),
                Arc::new(StringArray::from(vec!["protein_coding"])),
                Arc::new(Int64Array::from(vec![100])),
                Arc::new(Int64Array::from(vec![108])),
            ],
        )
        .expect("valid stop-loss transcript batch")
    }

    fn stop_loss_exons_batch() -> RecordBatch {
        let schema = Arc::new(Schema::new(vec![
            Field::new("transcript_id", DataType::Utf8, false),
            Field::new("exon_number", DataType::Int64, false),
            Field::new("start", DataType::Int64, false),
            Field::new("end", DataType::Int64, false),
        ]));
        RecordBatch::try_new(
            schema,
            vec![
                Arc::new(StringArray::from(vec!["ENST00000000007"])),
                Arc::new(Int64Array::from(vec![1])),
                Arc::new(Int64Array::from(vec![90])),
                Arc::new(Int64Array::from(vec![140])),
            ],
        )
        .expect("valid stop-loss exon batch")
    }

    fn stop_loss_translations_batch() -> RecordBatch {
        let schema = Arc::new(Schema::new(vec![
            Field::new("transcript_id", DataType::Utf8, false),
            Field::new("cds_len", DataType::Int64, true),
            Field::new("protein_len", DataType::Int64, true),
            Field::new("translation_seq", DataType::Utf8, true),
            Field::new("cds_sequence", DataType::Utf8, true),
        ]));
        RecordBatch::try_new(
            schema,
            vec![
                Arc::new(StringArray::from(vec!["ENST00000000007"])),
                Arc::new(Int64Array::from(vec![9])),
                Arc::new(Int64Array::from(vec![3])),
                Arc::new(StringArray::from(vec!["MA*"])),
                Arc::new(StringArray::from(vec!["ATGGCTTAA"])),
            ],
        )
        .expect("valid stop-loss translation batch")
    }

    fn hgvs_vcf_table() -> MemTable {
        let schema = Arc::new(Schema::new(vec![
            Field::new("chrom", DataType::Utf8, false),
            Field::new("start", DataType::Int64, false),
            Field::new("end", DataType::Int64, false),
            Field::new("ref", DataType::Utf8, false),
            Field::new("alt", DataType::Utf8, false),
        ]));
        let batch = RecordBatch::try_new(
            schema.clone(),
            vec![
                Arc::new(StringArray::from(vec!["1"])),
                Arc::new(Int64Array::from(vec![103])),
                Arc::new(Int64Array::from(vec![103])),
                Arc::new(StringArray::from(vec!["G"])),
                Arc::new(StringArray::from(vec!["A"])),
            ],
        )
        .expect("valid hgvs vcf batch");
        MemTable::try_new(schema, vec![vec![batch]]).expect("valid hgvs vcf memtable")
    }

    fn hgvs_cache_batch() -> RecordBatch {
        let schema = Arc::new(Schema::new(vec![
            Field::new("chrom", DataType::Utf8, false),
            Field::new("start", DataType::Int64, false),
            Field::new("end", DataType::Int64, false),
            Field::new("variation_name", DataType::Utf8, true),
            Field::new("allele_string", DataType::Utf8, false),
            Field::new("clin_sig", DataType::Utf8, true),
            Field::new("AF", DataType::Float64, true),
            Field::new("failed", DataType::Int64, false),
        ]));
        RecordBatch::try_new(
            schema,
            vec![
                Arc::new(StringArray::from(vec!["1"])),
                Arc::new(Int64Array::from(vec![103])),
                Arc::new(Int64Array::from(vec![103])),
                Arc::new(StringArray::from(vec!["rs_hgvs"])),
                Arc::new(StringArray::from(vec!["G/A"])),
                Arc::new(StringArray::from(vec!["pathogenic"])),
                Arc::new(Float64Array::from(vec![0.02_f64])),
                Arc::new(Int64Array::from(vec![0])),
            ],
        )
        .expect("valid hgvs cache batch")
    }

    fn hgvs_transcripts_batch() -> RecordBatch {
        let schema = Arc::new(Schema::new(vec![
            Field::new("transcript_id", DataType::Utf8, false),
            Field::new("chrom", DataType::Utf8, false),
            Field::new("start", DataType::Int64, false),
            Field::new("end", DataType::Int64, false),
            Field::new("strand", DataType::Int64, false),
            Field::new("biotype", DataType::Utf8, false),
            Field::new("cds_start", DataType::Int64, true),
            Field::new("cds_end", DataType::Int64, true),
            Field::new("version", DataType::Int64, true),
            Field::new("translation_stable_id", DataType::Utf8, true),
        ]));
        RecordBatch::try_new(
            schema,
            vec![
                Arc::new(StringArray::from(vec!["ENSTHGVS000001"])),
                Arc::new(StringArray::from(vec!["1"])),
                Arc::new(Int64Array::from(vec![90])),
                Arc::new(Int64Array::from(vec![140])),
                Arc::new(Int64Array::from(vec![1])),
                Arc::new(StringArray::from(vec!["protein_coding"])),
                Arc::new(Int64Array::from(vec![100])),
                Arc::new(Int64Array::from(vec![108])),
                Arc::new(Int64Array::from(vec![1])),
                Arc::new(StringArray::from(vec!["ENSPHGVS000001"])),
            ],
        )
        .expect("valid hgvs transcript batch")
    }

    fn hgvs_exons_batch() -> RecordBatch {
        let schema = Arc::new(Schema::new(vec![
            Field::new("transcript_id", DataType::Utf8, false),
            Field::new("exon_number", DataType::Int64, false),
            Field::new("start", DataType::Int64, false),
            Field::new("end", DataType::Int64, false),
        ]));
        RecordBatch::try_new(
            schema,
            vec![
                Arc::new(StringArray::from(vec!["ENSTHGVS000001"])),
                Arc::new(Int64Array::from(vec![1])),
                Arc::new(Int64Array::from(vec![90])),
                Arc::new(Int64Array::from(vec![140])),
            ],
        )
        .expect("valid hgvs exon batch")
    }

    fn hgvs_translations_batch() -> RecordBatch {
        let schema = Arc::new(Schema::new(vec![
            Field::new("transcript_id", DataType::Utf8, false),
            Field::new("cds_len", DataType::Int64, true),
            Field::new("protein_len", DataType::Int64, true),
            Field::new("translation_seq", DataType::Utf8, true),
            Field::new("cds_sequence", DataType::Utf8, true),
            Field::new("stable_id", DataType::Utf8, true),
            Field::new("version", DataType::Int64, true),
        ]));
        // CDS from 100..108 => ATG GCT TAA, variant 103 G>A gives Ala2Thr.
        RecordBatch::try_new(
            schema,
            vec![
                Arc::new(StringArray::from(vec!["ENSTHGVS000001"])),
                Arc::new(Int64Array::from(vec![9])),
                Arc::new(Int64Array::from(vec![3])),
                Arc::new(StringArray::from(vec!["MA*"])),
                Arc::new(StringArray::from(vec!["ATGGCTTAA"])),
                Arc::new(StringArray::from(vec!["ENSPHGVS000001"])),
                Arc::new(Int64Array::from(vec![1])),
            ],
        )
        .expect("valid hgvs translation batch")
    }

    fn write_test_indexed_fasta(
        sequence_name: &str,
        sequence: &str,
    ) -> (tempfile::TempDir, String) {
        let dir = tempfile::tempdir().expect("create temp fasta dir");
        let fasta_path = dir.path().join("test.fa");
        let fasta_contents = format!(">{sequence_name}\n{sequence}\n");
        std::fs::write(&fasta_path, fasta_contents).expect("write temp fasta");
        let fasta_index = format!(
            "{sequence_name}\t{}\t{}\t{}\t{}\n",
            sequence.len(),
            sequence_name.len() + 2,
            sequence.len(),
            sequence.len() + 1
        );
        std::fs::write(fasta_path.with_extension("fa.fai"), fasta_index)
            .expect("write temp fasta index");
        (
            dir,
            fasta_path
                .to_str()
                .expect("temp fasta path should be UTF-8")
                .to_string(),
        )
    }

    fn distance_vcf_table() -> MemTable {
        let schema = Arc::new(Schema::new(vec![
            Field::new("chrom", DataType::Utf8, false),
            Field::new("start", DataType::Int64, false),
            Field::new("end", DataType::Int64, false),
            Field::new("ref", DataType::Utf8, false),
            Field::new("alt", DataType::Utf8, false),
        ]));
        let batch = RecordBatch::try_new(
            schema.clone(),
            vec![
                Arc::new(StringArray::from(vec!["1"])),
                Arc::new(Int64Array::from(vec![30_000])),
                Arc::new(Int64Array::from(vec![30_000])),
                Arc::new(StringArray::from(vec!["A"])),
                Arc::new(StringArray::from(vec!["G"])),
            ],
        )
        .expect("valid distance vcf batch");
        MemTable::try_new(schema, vec![vec![batch]]).expect("valid distance vcf memtable")
    }

    fn distance_cache_batch() -> RecordBatch {
        let schema = Arc::new(Schema::new(vec![
            Field::new("chrom", DataType::Utf8, false),
            Field::new("start", DataType::Int64, false),
            Field::new("end", DataType::Int64, false),
            Field::new("variation_name", DataType::Utf8, true),
            Field::new("allele_string", DataType::Utf8, false),
            Field::new("clin_sig", DataType::Utf8, true),
            Field::new("AF", DataType::Float64, true),
            Field::new("failed", DataType::Int64, false),
        ]));
        RecordBatch::try_new(
            schema,
            vec![
                Arc::new(StringArray::from(vec!["1"])),
                Arc::new(Int64Array::from(vec![30_000])),
                Arc::new(Int64Array::from(vec![30_000])),
                Arc::new(StringArray::from(vec!["rs_distance"])),
                Arc::new(StringArray::from(vec!["A/G"])),
                Arc::new(StringArray::from(vec![None::<&str>])),
                Arc::new(Float64Array::from(vec![None::<f64>])),
                Arc::new(Int64Array::from(vec![0])),
            ],
        )
        .expect("valid distance cache batch")
    }

    fn distance_transcripts_batch() -> RecordBatch {
        let schema = Arc::new(Schema::new(vec![
            Field::new("transcript_id", DataType::Utf8, false),
            Field::new("chrom", DataType::Utf8, false),
            Field::new("start", DataType::Int64, false),
            Field::new("end", DataType::Int64, false),
            Field::new("strand", DataType::Int64, false),
            Field::new("biotype", DataType::Utf8, false),
            Field::new("cds_start", DataType::Int64, true),
            Field::new("cds_end", DataType::Int64, true),
        ]));
        RecordBatch::try_new(
            schema,
            vec![
                Arc::new(StringArray::from(vec!["ENSTDISTUP", "ENSTDISTDOWN"])),
                Arc::new(StringArray::from(vec!["1", "1"])),
                Arc::new(Int64Array::from(vec![37_079, 14_000])),
                Arc::new(Int64Array::from(vec![37_200, 15_000])),
                Arc::new(Int64Array::from(vec![1, 1])),
                Arc::new(StringArray::from(vec!["protein_coding", "protein_coding"])),
                Arc::new(Int64Array::from(vec![None::<i64>, None])),
                Arc::new(Int64Array::from(vec![None::<i64>, None])),
            ],
        )
        .expect("valid distance transcript batch")
    }

    fn distance_exons_batch() -> RecordBatch {
        let schema = Arc::new(Schema::new(vec![
            Field::new("transcript_id", DataType::Utf8, false),
            Field::new("exon_number", DataType::Int64, false),
            Field::new("start", DataType::Int64, false),
            Field::new("end", DataType::Int64, false),
        ]));
        RecordBatch::try_new(
            schema,
            vec![
                Arc::new(StringArray::from(vec!["ENSTDISTUP", "ENSTDISTDOWN"])),
                Arc::new(Int64Array::from(vec![1, 1])),
                Arc::new(Int64Array::from(vec![37_079, 14_000])),
                Arc::new(Int64Array::from(vec![37_200, 15_000])),
            ],
        )
        .expect("valid distance exon batch")
    }

    fn string_values(col: &Arc<dyn datafusion::arrow::array::Array>) -> Vec<Option<String>> {
        if let Some(arr) = col.as_any().downcast_ref::<StringArray>() {
            (0..arr.len())
                .map(|i| {
                    if arr.is_null(i) {
                        None
                    } else {
                        Some(arr.value(i).to_string())
                    }
                })
                .collect()
        } else if let Some(arr) = col
            .as_any()
            .downcast_ref::<datafusion::arrow::array::StringViewArray>()
        {
            (0..arr.len())
                .map(|i| {
                    if arr.is_null(i) {
                        None
                    } else {
                        Some(arr.value(i).to_string())
                    }
                })
                .collect()
        } else {
            panic!("expected string or string view");
        }
    }

    fn csq_entries(csq: &str) -> Vec<Vec<&str>> {
        csq.split(',')
            .map(|entry| entry.split('|').collect())
            .collect()
    }

    fn find_csq_entry<'a>(csq: &'a str, feature_type: &str, feature: &str) -> Vec<&'a str> {
        csq_entries(csq)
            .into_iter()
            .find(|fields| fields.len() == 74 && fields[5] == feature_type && fields[6] == feature)
            .expect("expected CSQ entry to exist")
    }

    fn regulatory_feature_batch(rows: &[(&str, &str, i64, i64, Option<&str>)]) -> RecordBatch {
        let schema = Arc::new(Schema::new(vec![
            Field::new("stable_id", DataType::Utf8, false),
            Field::new("chrom", DataType::Utf8, false),
            Field::new("start", DataType::Int64, false),
            Field::new("end", DataType::Int64, false),
            Field::new("feature_type", DataType::Utf8, true),
        ]));
        let stable_ids: Vec<&str> = rows.iter().map(|row| row.0).collect();
        let chroms: Vec<&str> = rows.iter().map(|row| row.1).collect();
        let starts: Vec<i64> = rows.iter().map(|row| row.2).collect();
        let ends: Vec<i64> = rows.iter().map(|row| row.3).collect();
        let feature_types: Vec<Option<&str>> = rows.iter().map(|row| row.4).collect();
        RecordBatch::try_new(
            schema,
            vec![
                Arc::new(StringArray::from(stable_ids)),
                Arc::new(StringArray::from(chroms)),
                Arc::new(Int64Array::from(starts)),
                Arc::new(Int64Array::from(ends)),
                Arc::new(StringArray::from(feature_types)),
            ],
        )
        .expect("valid regulatory feature batch")
    }

    fn extract_terms(csq: &str) -> Vec<String> {
        let mut parts = csq.split('|');
        let _allele = parts.next();
        let terms = parts.next().unwrap_or("");
        terms.split('&').map(|s| s.to_string()).collect()
    }

    fn extract_term_set(csq: &str) -> BTreeSet<String> {
        let mut out = BTreeSet::new();
        for ann in csq.split(',') {
            let mut parts = ann.split('|');
            let _allele = parts.next();
            let terms = parts.next().unwrap_or("");
            for term in terms.split('&').filter(|t| !t.is_empty()) {
                out.insert(term.to_string());
            }
        }
        out
    }

    fn assert_term_set_exact(csq: &str, expected_terms: &[&str]) {
        let expected: BTreeSet<String> = expected_terms.iter().map(|t| (*t).to_string()).collect();
        let observed = extract_term_set(csq);
        assert_eq!(
            observed, expected,
            "CSQ term set mismatch.\nobserved={observed:?}\nexpected={expected:?}"
        );
    }

    fn assert_terms_sorted_by_rank(csq: &str) {
        let terms = extract_terms(csq);
        let ranks = terms
            .iter()
            .map(|t| SoTerm::from_str(t).expect("term should parse").rank())
            .collect::<Vec<_>>();
        let mut sorted = ranks.clone();
        sorted.sort_unstable();
        assert_eq!(ranks, sorted, "terms should be rank-sorted in CSQ");
    }

    #[tokio::test(flavor = "multi_thread")]
    async fn test_annotate_vep_appends_annotation_columns() {
        let backend = "parquet";
        let tmpdir = TempDir::new().expect("create temp dir");
        write_batch_to_cache(&tmpdir, "variation", &cache_batch());

        let ctx = create_vep_session();
        ctx.register_table("vcf_data", Arc::new(vcf_table()))
            .expect("register vcf table");

        let cache_path = tmpdir.path().to_str().expect("utf8 path");
        let sql = format!(
            "SELECT * FROM annotate_vep('vcf_data', '{cache_path}', '{backend}', '{{\"partitioned\":true}}')"
        );
        let df = ctx
            .sql(&sql)
            .await
            .expect("annotate_vep query should parse");

        let batches = df.collect().await.expect("collect annotate_vep");
        let total_rows: usize = batches.iter().map(|b| b.num_rows()).sum();
        assert_eq!(total_rows, 2);

        let mut csq_values = Vec::new();
        let mut most_values = Vec::new();
        for batch in &batches {
            assert!(batch.column_by_name("CSQ").is_some());
            assert!(batch.column_by_name("most_severe_consequence").is_some());
            csq_values.extend(string_values(
                batch.column_by_name("CSQ").expect("csq column exists"),
            ));
            most_values.extend(string_values(
                batch
                    .column_by_name("most_severe_consequence")
                    .expect("most_severe_consequence column exists"),
            ));
        }
        // Partitioned path uses transcript engine — without context tables,
        // the variant is correctly classified as intergenic_variant.
        assert!(
            csq_values
                .iter()
                .any(|v| v.as_ref().is_some_and(|s| s.contains("intergenic_variant"))),
        );
        assert!(
            most_values
                .iter()
                .any(|v| v.as_ref() == Some(&"intergenic_variant".to_string())),
        );
    }

    #[tokio::test(flavor = "multi_thread")]
    async fn test_annotate_vep_projection_includes_null_placeholder_fields() {
        let backend = "parquet";
        let tmpdir = TempDir::new().expect("create temp dir");
        write_batch_to_cache(&tmpdir, "variation", &cache_batch());

        let ctx = create_vep_session();
        ctx.register_table("vcf_data", Arc::new(vcf_table()))
            .expect("register vcf table");

        let cache_path = tmpdir.path().to_str().expect("utf8 path");
        let sql = format!(
            "SELECT chrom, \"CSQ\" FROM annotate_vep('vcf_data', '{cache_path}', '{backend}', '{{\"partitioned\":true}}')"
        );
        let df = ctx.sql(&sql).await.expect("projection query should parse");

        let batches = df.collect().await.expect("collect projected annotate_vep");
        let total_rows: usize = batches.iter().map(|b| b.num_rows()).sum();
        assert_eq!(total_rows, 2);
        for batch in &batches {
            assert_eq!(batch.num_columns(), 2);
            assert_eq!(batch.schema().field(0).name(), "chrom");
            assert_eq!(batch.schema().field(1).name(), "CSQ");
        }
        let mut csq_values = Vec::new();
        for batch in &batches {
            csq_values.extend(string_values(
                batch.column_by_name("CSQ").expect("csq column exists"),
            ));
        }
        // Both rows have CSQ (intergenic_variant when no context tables).
        assert!(csq_values.iter().all(|v| v.is_some()));
    }

    #[tokio::test(flavor = "multi_thread")]
    async fn test_annotate_vep_uses_transcript_context_tables_when_available() {
        let backend = "parquet";
        let tmpdir = TempDir::new().expect("create temp dir");
        write_batch_to_cache(&tmpdir, "variation", &cache_batch());
        write_batch_to_cache(&tmpdir, "transcript", &transcripts_batch());
        // exons_batch() has no chrom column; write per-chrom explicitly
        write_batch_to_chrom(&tmpdir, "exon", "1", &exons_batch());
        write_batch_to_chrom(&tmpdir, "exon", "2", &exons_batch());

        let ctx = create_vep_session();
        ctx.register_table("vcf_data", Arc::new(vcf_table()))
            .expect("register vcf table");

        let cache_path = tmpdir.path().to_str().expect("utf8 path");
        let sql = format!(
            "SELECT chrom, \"CSQ\", most_severe_consequence \
             FROM annotate_vep('vcf_data', '{cache_path}', '{backend}', '{{\"partitioned\":true}}') \
             ORDER BY chrom"
        );
        let df = ctx.sql(&sql).await.expect("query should parse");

        let batches = df
            .collect()
            .await
            .expect("collect transcript-aware annotate_vep");
        let total_rows: usize = batches.iter().map(|b| b.num_rows()).sum();
        assert_eq!(total_rows, 2);

        let mut chrom = Vec::new();
        let mut csq = Vec::new();
        let mut most = Vec::new();
        for batch in &batches {
            chrom.extend(string_values(
                batch.column_by_name("chrom").expect("chrom column exists"),
            ));
            csq.extend(string_values(
                batch.column_by_name("CSQ").expect("csq column exists"),
            ));
            most.extend(string_values(
                batch
                    .column_by_name("most_severe_consequence")
                    .expect("most_severe_consequence column exists"),
            ));
        }

        assert_eq!(chrom, vec![Some("1".to_string()), Some("2".to_string())],);
        assert!(csq.iter().all(|v| v.is_some()));
        assert!(most.iter().all(|v| v.is_some()));
        assert!(
            csq[0]
                .as_ref()
                .is_some_and(|s| s.contains("coding_sequence_variant")),
        );
        assert!(
            csq[1]
                .as_ref()
                .is_some_and(|s| s.contains("non_coding_transcript_exon_variant")),
        );
        assert_eq!(
            most,
            vec![
                Some("coding_sequence_variant".to_string()),
                Some("non_coding_transcript_exon_variant".to_string())
            ],
        );
    }

    #[tokio::test(flavor = "multi_thread")]
    async fn test_annotate_vep_uses_partitioned_context_for_transcript_annotation() {
        let backend = "parquet";
        let tmpdir = TempDir::new().expect("create temp dir");
        write_batch_to_cache(&tmpdir, "variation", &cache_batch());
        write_batch_to_cache(&tmpdir, "transcript", &transcripts_batch());
        write_batch_to_chrom(&tmpdir, "exon", "1", &exons_batch());
        write_batch_to_chrom(&tmpdir, "exon", "2", &exons_batch());

        let ctx = create_vep_session();
        ctx.register_table("vcf_data", Arc::new(vcf_table()))
            .expect("register vcf table");

        let cache_path = tmpdir.path().to_str().expect("utf8 path");
        let sql = format!(
            "SELECT \"CSQ\", most_severe_consequence \
             FROM annotate_vep( \
               'vcf_data', \
               '{cache_path}', \
               '{backend}', \
               '{{\"partitioned\":true}}' \
             )"
        );
        let df = ctx.sql(&sql).await.expect("query should parse");

        let batches = df
            .collect()
            .await
            .expect("collect transcript-aware annotate_vep");
        let mut csq = Vec::new();
        let mut most = Vec::new();
        for batch in &batches {
            csq.extend(string_values(
                batch.column_by_name("CSQ").expect("csq column exists"),
            ));
            most.extend(string_values(
                batch
                    .column_by_name("most_severe_consequence")
                    .expect("most_severe_consequence column exists"),
            ));
        }
        assert!(csq.iter().all(|v| v.is_some()));
        assert!(most.iter().all(|v| v.is_some()));
    }

    // Mirrors Ensembl VEP Runner distance coverage:
    // - https://github.com/Ensembl/ensembl-vep/blob/release/115/t/Runner.t#L535-L571
    #[tokio::test(flavor = "multi_thread")]
    async fn test_annotate_vep_respects_options_json_distance_for_upstream_and_downstream() {
        let backend = "parquet";
        let tmpdir = TempDir::new().expect("create temp dir");
        write_batch_to_cache(&tmpdir, "variation", &distance_cache_batch());
        write_batch_to_cache(&tmpdir, "transcript", &distance_transcripts_batch());
        write_batch_to_chrom(&tmpdir, "exon", "1", &distance_exons_batch());

        let ctx = create_vep_session();
        ctx.register_table("vcf_distance", Arc::new(distance_vcf_table()))
            .expect("register distance vcf");

        let cache_path = tmpdir.path().to_str().expect("utf8 path");

        let default_batches = ctx
            .sql(&format!(
                "SELECT \"CSQ\" FROM annotate_vep('vcf_distance', '{cache_path}', '{backend}', '{{\"partitioned\":true}}')"
            ))
            .await
            .expect("default distance query should parse")
            .collect()
            .await
            .expect("collect default distance annotate_vep");
        let default_csq = string_values(
            default_batches[0]
                .column_by_name("CSQ")
                .expect("default csq column exists"),
        );
        let default_csq0 = default_csq[0]
            .as_ref()
            .expect("default csq should be present");
        assert!(default_csq0.contains("intergenic_variant"));
        assert!(!default_csq0.contains("ENSTDISTUP"));
        assert!(!default_csq0.contains("ENSTDISTDOWN"));

        let numeric_batches = ctx
            .sql(&format!(
                "SELECT \"CSQ\" FROM annotate_vep( \
                   'vcf_distance', \
                   '{cache_path}', \
                   '{backend}', \
                   '{{\"partitioned\":true,\"distance\":10000}}' \
                 )"
            ))
            .await
            .expect("numeric distance query should parse")
            .collect()
            .await
            .expect("collect numeric distance annotate_vep");
        let numeric_csq = string_values(
            numeric_batches[0]
                .column_by_name("CSQ")
                .expect("numeric csq column exists"),
        );
        let numeric_csq0 = numeric_csq[0]
            .as_ref()
            .expect("numeric csq should be present");
        let numeric_entries = csq_entries(numeric_csq0);
        assert_eq!(numeric_entries.len(), 1);
        let upstream_entry = find_csq_entry(numeric_csq0, "Transcript", "ENSTDISTUP");
        assert_eq!(upstream_entry[1], "upstream_gene_variant");
        assert_eq!(upstream_entry[18], "7079");
        assert!(!numeric_csq0.contains("ENSTDISTDOWN"));

        let pair_batches = ctx
            .sql(&format!(
                "SELECT \"CSQ\" FROM annotate_vep( \
                   'vcf_distance', \
                   '{cache_path}', \
                   '{backend}', \
                   '{{\"partitioned\":true,\"distance\":\"10000,20000\"}}' \
                 )"
            ))
            .await
            .expect("pair distance query should parse")
            .collect()
            .await
            .expect("collect pair distance annotate_vep");
        let pair_csq = string_values(
            pair_batches[0]
                .column_by_name("CSQ")
                .expect("pair csq column exists"),
        );
        let pair_csq0 = pair_csq[0].as_ref().expect("pair csq should be present");
        let pair_entries = csq_entries(pair_csq0);
        assert_eq!(pair_entries.len(), 2);
        let upstream_entry = find_csq_entry(pair_csq0, "Transcript", "ENSTDISTUP");
        let downstream_entry = find_csq_entry(pair_csq0, "Transcript", "ENSTDISTDOWN");
        assert_eq!(upstream_entry[1], "upstream_gene_variant");
        assert_eq!(upstream_entry[18], "7079");
        assert_eq!(downstream_entry[1], "downstream_gene_variant");
        assert_eq!(downstream_entry[18], "15000");
    }

    // Mirrors Ensembl VEP offline HGVS gating and output toggles:
    // - https://github.com/Ensembl/ensembl-vep/blob/release/115/modules/Bio/EnsEMBL/VEP/Runner.pm#L726-L738
    // - https://github.com/Ensembl/ensembl-vep/blob/release/115/modules/Bio/EnsEMBL/VEP/OutputFactory.pm#L1698-L1715
    #[tokio::test(flavor = "multi_thread")]
    async fn test_annotate_vep_hgvs_fields_require_flag_and_reference_fasta() {
        let backend = "parquet";
        let tmpdir = TempDir::new().expect("create temp dir");
        write_batch_to_cache(&tmpdir, "variation", &hgvs_cache_batch());
        write_batch_to_cache(&tmpdir, "transcript", &hgvs_transcripts_batch());
        write_batch_to_chrom(&tmpdir, "exon", "1", &hgvs_exons_batch());
        write_batch_to_chrom(&tmpdir, "translation_core", "1", &hgvs_translations_batch());

        let ctx = create_vep_session();
        ctx.register_table("vcf_hgvs", Arc::new(hgvs_vcf_table()))
            .expect("register hgvs vcf");

        let cache_path = tmpdir.path().to_str().expect("utf8 path");

        let default_batches = ctx
            .sql(&format!(
                "SELECT \"CSQ\" FROM annotate_vep('vcf_hgvs', '{cache_path}', '{backend}', '{{\"partitioned\":true}}')"
            ))
            .await
            .expect("default hgvs query should parse")
            .collect()
            .await
            .expect("collect default hgvs annotate_vep");
        let default_csq = string_values(
            default_batches[0]
                .column_by_name("CSQ")
                .expect("default csq column exists"),
        );
        let default_entry = csq_entries(default_csq[0].as_ref().expect("default csq present"))
            .into_iter()
            .find(|fields| fields.len() == 74 && fields[5] == "Transcript")
            .expect("transcript entry should exist");
        assert_eq!(default_entry[10], "");
        assert_eq!(default_entry[11], "");

        let missing_fasta_err = ctx
            .sql(&format!(
                "SELECT \"CSQ\" FROM annotate_vep( \
                   'vcf_hgvs', \
                   '{cache_path}', \
                   '{backend}', \
                   '{{\"partitioned\":true,\"hgvs\":true}}' \
                 )"
            ))
            .await
            .expect("hgvs query should parse")
            .collect()
            .await
            .expect_err("hgvs without reference fasta should fail")
            .to_string();
        assert!(missing_fasta_err.contains("Cannot generate HGVS coordinates"));

        let (_fasta_dir, fasta_path) = write_test_indexed_fasta("1", "NNNNNNNNNN");

        let hgvs_batches = ctx
            .sql(&format!(
                "SELECT \"CSQ\" FROM annotate_vep( \
                   'vcf_hgvs', \
                   '{cache_path}', \
                   '{backend}', \
                   '{{\"partitioned\":true,\"hgvs\":true,\"reference_fasta_path\":\"{}\"}}' \
                 )",
                fasta_path.replace('\'', "''")
            ))
            .await
            .expect("hgvs enabled query should parse")
            .collect()
            .await
            .expect("collect hgvs enabled annotate_vep");
        let hgvs_csq = string_values(
            hgvs_batches[0]
                .column_by_name("CSQ")
                .expect("hgvs csq column exists"),
        );
        let hgvs_entry = csq_entries(hgvs_csq[0].as_ref().expect("hgvs csq present"))
            .into_iter()
            .find(|fields| fields.len() == 74 && fields[5] == "Transcript")
            .expect("transcript hgvs entry should exist");
        assert_eq!(hgvs_entry[10], "ENSTHGVS000001.1:c.4G>A");
        assert_eq!(hgvs_entry[11], "ENSPHGVS000001.1:p.Ala2Thr");

        let hgvs_prediction_batches = ctx
            .sql(&format!(
                "SELECT \"CSQ\" FROM annotate_vep( \
                   'vcf_hgvs', \
                   '{cache_path}', \
                   '{backend}', \
                   '{{\"partitioned\":true,\"hgvsp\":true,\"hgvsp_use_prediction\":true,\"reference_fasta_path\":\"{}\"}}' \
                 )",
                fasta_path.replace('\'', "''")
            ))
            .await
            .expect("hgvs prediction query should parse")
            .collect()
            .await
            .expect("collect hgvs prediction annotate_vep");
        let hgvs_prediction_csq = string_values(
            hgvs_prediction_batches[0]
                .column_by_name("CSQ")
                .expect("hgvs prediction csq column exists"),
        );
        let hgvs_prediction_entry = csq_entries(
            hgvs_prediction_csq[0]
                .as_ref()
                .expect("hgvs prediction csq present"),
        )
        .into_iter()
        .find(|fields| fields.len() == 74 && fields[5] == "Transcript")
        .expect("transcript hgvs prediction entry should exist");
        assert_eq!(hgvs_prediction_entry[10], "");
        assert_eq!(hgvs_prediction_entry[11], "ENSPHGVS000001.1:p.(Ala2Thr)",);
    }

    #[tokio::test(flavor = "multi_thread")]
    async fn test_annotate_vep_uses_translation_context_for_synonymous_classification() {
        let backend = "parquet";
        let tmpdir = TempDir::new().expect("create temp dir");
        write_batch_to_cache(&tmpdir, "variation", &synonymous_cache_batch());
        write_batch_to_cache(&tmpdir, "transcript", &synonymous_transcripts_batch());
        write_batch_to_chrom(&tmpdir, "exon", "1", &synonymous_exons_batch());
        write_batch_to_chrom(
            &tmpdir,
            "translation_core",
            "1",
            &synonymous_translations_batch(),
        );

        let ctx = create_vep_session();
        ctx.register_table("vcf_syn", Arc::new(synonymous_vcf_table()))
            .expect("register synonymous vcf table");

        let cache_path = tmpdir.path().to_str().expect("utf8 path");
        let sql = format!(
            "SELECT \"CSQ\", most_severe_consequence \
             FROM annotate_vep('vcf_syn', '{cache_path}', '{backend}', '{{\"partitioned\":true}}')"
        );
        let df = ctx.sql(&sql).await.expect("query should parse");
        let batches = df.collect().await.expect("collect annotate_vep");
        let mut csq = Vec::new();
        let mut most = Vec::new();
        for batch in &batches {
            csq.extend(string_values(
                batch.column_by_name("CSQ").expect("csq column exists"),
            ));
            most.extend(string_values(
                batch
                    .column_by_name("most_severe_consequence")
                    .expect("most_severe_consequence column exists"),
            ));
        }
        assert_eq!(csq.len(), 1);
        assert_eq!(most.len(), 1);
        let csq0 = csq[0].as_ref().expect("csq should be present");
        assert!(csq0.contains("synonymous_variant"));
        assert!(!csq0.contains("missense_variant"));
        assert_eq!(
            most[0],
            Some("synonymous_variant".to_string()),
            "translation-aware codon classification should produce synonymous_variant",
        );
    }

    #[tokio::test(flavor = "multi_thread")]
    async fn test_annotate_vep_csq_terms_deterministic_with_context_tables() {
        let backend = "parquet";
        let tmpdir = TempDir::new().expect("create temp dir");
        write_batch_to_cache(&tmpdir, "variation", &context_cache_batch());
        write_batch_to_cache(&tmpdir, "transcript", &context_transcripts_batch());
        write_batch_to_chrom(&tmpdir, "exon", "1", &context_exons_batch());
        write_batch_to_cache(&tmpdir, "regulatory", &context_regulatory_batch());
        write_batch_to_cache(&tmpdir, "motif", &context_motif_batch());

        let ctx = create_vep_session();
        ctx.register_table("vcf_ctx", Arc::new(context_vcf_table()))
            .expect("register context vcf table");

        let cache_path = tmpdir.path().to_str().expect("utf8 path");
        let sql = format!(
            "SELECT \"CSQ\", most_severe_consequence \
             FROM annotate_vep('vcf_ctx', '{cache_path}', '{backend}', '{{\"partitioned\":true}}')"
        );
        let batches = ctx
            .sql(&sql)
            .await
            .expect("query should parse")
            .collect()
            .await
            .expect("collect annotate_vep");
        let csq = string_values(batches[0].column_by_name("CSQ").expect("csq column exists"));
        let most = string_values(
            batches[0]
                .column_by_name("most_severe_consequence")
                .expect("most_severe_consequence column exists"),
        );
        let csq0 = csq[0].as_ref().expect("csq should be present").to_string();
        let most0 = most[0]
            .as_ref()
            .expect("most severe should be present")
            .to_string();

        assert_terms_sorted_by_rank(&csq0);
        assert!(csq0.contains("regulatory_region_variant"));
        assert!(csq0.contains("TF_binding_site_variant"));
        // miRNA and SV features are not yet partitioned, so those terms
        // are not present in the partitioned path output.
        assert_eq!(most0, "coding_sequence_variant");

        // Verify determinism: run the same query again.
        let batches2 = ctx
            .sql(&sql)
            .await
            .expect("second query should parse")
            .collect()
            .await
            .expect("collect second annotate_vep");
        let csq2 = string_values(
            batches2[0]
                .column_by_name("CSQ")
                .expect("csq column exists"),
        );
        let csq2_0 = csq2[0].as_ref().expect("csq should be present").to_string();
        assert_eq!(csq0, csq2_0, "CSQ should be deterministic across runs");
    }

    // Verifies that CSQ entries are sorted lexicographically by Feature ID
    // within each Feature_type group, matching Ensembl VEP behaviour.
    // See: https://github.com/biodatageeks/datafusion-bio-functions/issues/83
    #[tokio::test(flavor = "multi_thread")]
    async fn test_csq_entries_sorted_by_feature_id_within_feature_type() {
        let backend = "parquet";
        let tmpdir = TempDir::new().expect("create temp dir");
        write_batch_to_cache(&tmpdir, "variation", &context_cache_batch());

        // Three transcripts overlapping position 155 on chr 1, deliberately
        // supplied in non-lexicographic order (C > A > B).
        let tx_schema = Arc::new(Schema::new(vec![
            Field::new("transcript_id", DataType::Utf8, false),
            Field::new("chrom", DataType::Utf8, false),
            Field::new("start", DataType::Int64, false),
            Field::new("end", DataType::Int64, false),
            Field::new("strand", DataType::Int64, false),
            Field::new("biotype", DataType::Utf8, false),
            Field::new("cds_start", DataType::Int64, true),
            Field::new("cds_end", DataType::Int64, true),
        ]));
        let tx_batch = RecordBatch::try_new(
            tx_schema,
            vec![
                Arc::new(StringArray::from(vec![
                    "ENST00000900000",
                    "ENST00000100000",
                    "ENST00000500000",
                ])),
                Arc::new(StringArray::from(vec!["1", "1", "1"])),
                Arc::new(Int64Array::from(vec![100, 100, 100])),
                Arc::new(Int64Array::from(vec![250, 250, 250])),
                Arc::new(Int64Array::from(vec![1, 1, 1])),
                Arc::new(StringArray::from(vec![
                    "protein_coding",
                    "protein_coding",
                    "protein_coding",
                ])),
                Arc::new(Int64Array::from(vec![120, 120, 120])),
                Arc::new(Int64Array::from(vec![240, 240, 240])),
            ],
        )
        .expect("valid multi-transcript batch");
        write_batch_to_cache(&tmpdir, "transcript", &tx_batch);

        // Exons for all three transcripts.
        let exon_schema = Arc::new(Schema::new(vec![
            Field::new("transcript_id", DataType::Utf8, false),
            Field::new("exon_number", DataType::Int64, false),
            Field::new("start", DataType::Int64, false),
            Field::new("end", DataType::Int64, false),
        ]));
        let exon_batch = RecordBatch::try_new(
            exon_schema,
            vec![
                Arc::new(StringArray::from(vec![
                    "ENST00000900000",
                    "ENST00000100000",
                    "ENST00000500000",
                ])),
                Arc::new(Int64Array::from(vec![1, 1, 1])),
                Arc::new(Int64Array::from(vec![100, 100, 100])),
                Arc::new(Int64Array::from(vec![250, 250, 250])),
            ],
        )
        .expect("valid multi-exon batch");
        write_batch_to_chrom(&tmpdir, "exon", "1", &exon_batch);

        let ctx = create_vep_session();
        ctx.register_table("vcf_sort", Arc::new(context_vcf_table()))
            .expect("register vcf table");

        let cache_path = tmpdir.path().to_str().expect("utf8 path");
        let sql = format!(
            "SELECT \"CSQ\" \
             FROM annotate_vep('vcf_sort', '{cache_path}', '{backend}', '{{\"partitioned\":true}}')"
        );
        let batches = ctx
            .sql(&sql)
            .await
            .expect("query should parse")
            .collect()
            .await
            .expect("collect annotate_vep");
        let csq = string_values(batches[0].column_by_name("CSQ").expect("csq column"));
        let csq0 = csq[0].as_ref().expect("csq should be present");
        let entries = csq_entries(csq0);
        // 74 = number of pipe-delimited CSQ FORMAT fields per entry.
        let feature_ids: Vec<&str> = entries
            .iter()
            .filter(|f| f.len() == 74 && f[5] == "Transcript")
            .map(|f| f[6])
            .collect();
        assert!(
            !feature_ids.is_empty(),
            "expected at least one Transcript CSQ entry"
        );
        assert_eq!(
            feature_ids,
            vec!["ENST00000100000", "ENST00000500000", "ENST00000900000"],
            "CSQ transcript entries must be sorted lexicographically by Feature ID"
        );
    }

    // Verifies that CSQ entries are grouped by Feature type in VEP order
    // (Transcript → RegulatoryFeature → MotifFeature) and sorted within
    // each group by feature ID.
    // See: https://github.com/biodatageeks/datafusion-bio-functions/issues/83
    #[tokio::test(flavor = "multi_thread")]
    async fn test_csq_entries_grouped_by_feature_type_then_sorted_by_id() {
        let backend = "parquet";
        let tmpdir = TempDir::new().expect("create temp dir");
        write_batch_to_cache(&tmpdir, "variation", &context_cache_batch());
        write_batch_to_cache(&tmpdir, "transcript", &context_transcripts_batch());
        write_batch_to_chrom(&tmpdir, "exon", "1", &context_exons_batch());

        // Two regulatory features in reverse order.
        write_batch_to_cache(
            &tmpdir,
            "regulatory",
            &regulatory_feature_batch(&[
                ("ENSR_BBB", "1", 150, 160, Some("promoter")),
                ("ENSR_AAA", "1", 150, 160, Some("enhancer")),
            ]),
        );
        write_batch_to_cache(&tmpdir, "motif", &context_motif_batch());

        let ctx = create_vep_session();
        ctx.register_table("vcf_group", Arc::new(context_vcf_table()))
            .expect("register vcf table");

        let cache_path = tmpdir.path().to_str().expect("utf8 path");
        let sql = format!(
            "SELECT \"CSQ\" \
             FROM annotate_vep('vcf_group', '{cache_path}', '{backend}', '{{\"partitioned\":true}}')"
        );
        let batches = ctx
            .sql(&sql)
            .await
            .expect("query should parse")
            .collect()
            .await
            .expect("collect annotate_vep");
        let csq = string_values(batches[0].column_by_name("CSQ").expect("csq column"));
        let csq0 = csq[0].as_ref().expect("csq should be present");
        let entries = csq_entries(csq0);
        // 74 = number of pipe-delimited CSQ FORMAT fields per entry.
        let feature_types: Vec<&str> = entries
            .iter()
            .filter(|f| f.len() == 74)
            .map(|f| f[5])
            .collect();
        assert!(
            !feature_types.is_empty(),
            "expected at least one CSQ entry with 74 fields"
        );

        // Transcript entries must come before RegulatoryFeature entries,
        // which must come before MotifFeature entries.
        let mut seen_reg = false;
        let mut seen_motif = false;
        for ft in &feature_types {
            match *ft {
                "Transcript" => {
                    assert!(
                        !seen_reg && !seen_motif,
                        "Transcript must appear before Regulatory and Motif"
                    );
                }
                "RegulatoryFeature" => {
                    assert!(!seen_motif, "RegulatoryFeature must appear before Motif");
                    seen_reg = true;
                }
                "MotifFeature" => {
                    seen_motif = true;
                }
                _ => {}
            }
        }
        assert!(seen_reg, "expected at least one RegulatoryFeature entry");
        assert!(seen_motif, "expected at least one MotifFeature entry");

        // Regulatory entries should be sorted by stable_id.
        // 74 = number of pipe-delimited CSQ FORMAT fields per entry.
        let reg_ids: Vec<&str> = entries
            .iter()
            .filter(|f| f.len() == 74 && f[5] == "RegulatoryFeature")
            .map(|f| f[6])
            .collect();
        let mut sorted_reg = reg_ids.clone();
        sorted_reg.sort();
        assert_eq!(
            reg_ids, sorted_reg,
            "Regulatory CSQ entries must be sorted by stable_id"
        );
    }

    // Mirrors the cache-backed regulatory source coverage from:
    // - https://github.com/Ensembl/ensembl-vep/blob/release/115/t/AnnotationSource_Cache_RegFeat.t#L166-L205
    #[tokio::test(flavor = "multi_thread")]
    async fn test_annotate_vep_deduplicates_duplicate_regulatory_rows_from_context_table() {
        let backend = "parquet";
        let tmpdir = TempDir::new().expect("create temp dir");
        write_batch_to_cache(&tmpdir, "variation", &context_cache_batch());

        let duplicate_batch = regulatory_feature_batch(&[
            ("ENSR_DUP", "1", 150, 160, Some("promoter")),
            ("ENSR_DUP", "1", 150, 160, Some("promoter")),
        ]);
        write_batch_to_cache(&tmpdir, "regulatory", &duplicate_batch);

        let ctx = create_vep_session();
        ctx.register_table("vcf_reg_dup", Arc::new(context_vcf_table()))
            .expect("register duplicate regulatory vcf");

        let cache_path = tmpdir.path().to_str().expect("utf8 path");
        let sql = format!(
            "SELECT \"CSQ\", most_severe_consequence \
             FROM annotate_vep( \
               'vcf_reg_dup', \
               '{cache_path}', \
               '{backend}', \
               '{{\"partitioned\":true}}' \
             )"
        );
        let batches = ctx
            .sql(&sql)
            .await
            .expect("duplicate regulatory query should parse")
            .collect()
            .await
            .expect("collect duplicate regulatory annotate_vep");

        let csq = string_values(batches[0].column_by_name("CSQ").expect("csq column exists"));
        let most = string_values(
            batches[0]
                .column_by_name("most_severe_consequence")
                .expect("most severe column exists"),
        );
        let csq0 = csq[0].as_ref().expect("csq should be present");
        let entries = csq_entries(csq0);
        let regulatory_entries: Vec<&Vec<&str>> = entries
            .iter()
            .filter(|fields| fields.len() == 74 && fields[5] == "RegulatoryFeature")
            .collect();
        assert_eq!(regulatory_entries.len(), 1);
        assert_eq!(regulatory_entries[0][6], "ENSR_DUP");
        assert_eq!(regulatory_entries[0][7], "promoter");
        assert_eq!(most[0], Some("regulatory_region_variant".to_string()),);
    }

    // Mirrors Ensembl VEP regulatory output serialization coverage:
    // - https://github.com/Ensembl/ensembl-vep/blob/release/115/t/OutputFactory.t#L1207-L1236
    // - https://github.com/Ensembl/ensembl-vep/blob/release/115/modules/Bio/EnsEMBL/VEP/OutputFactory.pm#L1802-L1829
    #[tokio::test(flavor = "multi_thread")]
    async fn test_annotate_vep_regulatory_csq_serializer_matches_vep_output_factory_shape() {
        let backend = "parquet";
        let tmpdir = TempDir::new().expect("create temp dir");
        write_batch_to_cache(&tmpdir, "variation", &context_cache_batch());
        let reg_batch =
            regulatory_feature_batch(&[("ENSR00000140763", "1", 150, 160, Some("promoter"))]);
        write_batch_to_cache(&tmpdir, "regulatory", &reg_batch);

        let ctx = create_vep_session();
        ctx.register_table("vcf_reg_ser", Arc::new(context_vcf_table()))
            .expect("register regulatory serializer vcf");

        let cache_path = tmpdir.path().to_str().expect("utf8 path");
        let sql = format!(
            "SELECT \"CSQ\", most_severe_consequence \
             FROM annotate_vep('vcf_reg_ser', '{cache_path}', '{backend}', '{{\"partitioned\":true}}')"
        );
        let batches = ctx
            .sql(&sql)
            .await
            .expect("regulatory serializer query should parse")
            .collect()
            .await
            .expect("collect regulatory serializer annotate_vep");

        let csq = string_values(batches[0].column_by_name("CSQ").expect("csq column exists"));
        let most = string_values(
            batches[0]
                .column_by_name("most_severe_consequence")
                .expect("most severe column exists"),
        );
        let csq0 = csq[0].as_ref().expect("csq should be present");
        let entries = csq_entries(csq0);
        assert_eq!(entries.len(), 2);

        let entry = find_csq_entry(csq0, "RegulatoryFeature", "ENSR00000140763");
        assert_eq!(entry.len(), 74);
        assert_eq!(entry[0], "G");
        assert_eq!(entry[1], "regulatory_region_variant");
        assert_eq!(entry[2], "MODIFIER");
        assert_eq!(entry[3], "");
        assert_eq!(entry[4], "");
        assert_eq!(entry[5], "RegulatoryFeature");
        assert_eq!(entry[6], "ENSR00000140763");
        assert_eq!(entry[7], "promoter");
        assert_eq!(entry[8], "");
        assert_eq!(entry[9], "");
        assert_eq!(entry[18], "");
        assert!(
            entries
                .iter()
                .any(|fields| fields.len() == 74 && fields[1] == "intergenic_variant"),
            "regulatory-only output should retain the orthogonal intergenic entry"
        );
        assert_eq!(most[0], Some("regulatory_region_variant".to_string()),);
    }

    #[tokio::test(flavor = "multi_thread")]
    async fn test_annotate_vep_golden_fixture_term_parity_for_context_reg_motif() {
        let backend = "parquet";
        let tmpdir = TempDir::new().expect("create temp dir");
        write_batch_to_cache(&tmpdir, "variation", &golden_context_cache_batch());
        write_batch_to_cache(&tmpdir, "transcript", &golden_context_transcripts_batch());
        write_batch_to_chrom(&tmpdir, "exon", "1", &golden_context_exons_batch());
        write_batch_to_cache(&tmpdir, "regulatory", &golden_context_regulatory_batch());
        write_batch_to_cache(&tmpdir, "motif", &golden_context_motif_batch());

        let ctx = create_vep_session();
        ctx.register_table("vcf_golden_ctx", Arc::new(golden_context_vcf_table()))
            .expect("register golden context vcf table");

        // miRNA and structural features are not yet partitioned, so those
        // terms are excluded from the expected set.
        let expected_terms = [
            "TF_binding_site_variant",
            "regulatory_region_variant",
            "intergenic_variant",
        ];

        let cache_path = tmpdir.path().to_str().expect("utf8 path");
        let sql = format!(
            "SELECT \"CSQ\", most_severe_consequence \
             FROM annotate_vep('vcf_golden_ctx', '{cache_path}', '{backend}', '{{\"partitioned\":true}}')"
        );
        let batches = ctx
            .sql(&sql)
            .await
            .expect("query should parse")
            .collect()
            .await
            .expect("collect annotate_vep");
        let csq = string_values(batches[0].column_by_name("CSQ").expect("csq column exists"));
        let most = string_values(
            batches[0]
                .column_by_name("most_severe_consequence")
                .expect("most_severe_consequence column exists"),
        );
        let csq0 = csq[0].as_ref().expect("csq should be present").to_string();
        let most0 = most[0]
            .as_ref()
            .expect("most severe should be present")
            .to_string();
        assert_terms_sorted_by_rank(&csq0);
        assert_term_set_exact(&csq0, &expected_terms);
        assert_eq!(most0, "TF_binding_site_variant");
    }

    #[tokio::test(flavor = "multi_thread")]
    async fn test_annotate_vep_golden_fixture_term_parity_for_splice_boundary_case() {
        let ctx = create_vep_session();
        ctx.register_table("vcf_splice", Arc::new(splice_vcf_table()))
            .expect("register splice vcf");
        let tmpdir = TempDir::new().expect("create tmpdir");
        write_batch_to_cache(&tmpdir, "variation", &splice_cache_batch());
        write_batch_to_cache(&tmpdir, "transcript", &splice_transcripts_batch());
        write_batch_to_chrom(&tmpdir, "exon", "1", &splice_exons_batch());
        let cache_path = tmpdir.path().display().to_string();

        let expected_terms = ["splice_donor_variant"];

        {
            let sql = format!(
                "SELECT \"CSQ\", most_severe_consequence \
                 FROM annotate_vep('vcf_splice', '{cache_path}', 'parquet', '{{\"partitioned\":true}}')"
            );
            let batches = ctx
                .sql(&sql)
                .await
                .expect("query should parse")
                .collect()
                .await
                .expect("collect annotate_vep");
            let csq = string_values(batches[0].column_by_name("CSQ").expect("csq column exists"));
            let most = string_values(
                batches[0]
                    .column_by_name("most_severe_consequence")
                    .expect("most_severe_consequence column exists"),
            );
            let csq0 = csq[0].as_ref().expect("csq should be present");
            assert_terms_sorted_by_rank(csq0);
            assert_term_set_exact(csq0, &expected_terms);
            assert_eq!(most[0], Some("splice_donor_variant".to_string()));
        }
    }

    #[tokio::test(flavor = "multi_thread")]
    async fn test_annotate_vep_golden_fixture_term_parity_for_repeat_shifted_deletion() {
        let ctx = create_vep_session();
        ctx.register_table("vcf_repeat", Arc::new(repeat_shift_vcf_table()))
            .expect("register repeat vcf");
        let tmpdir = TempDir::new().expect("create tmpdir");
        write_batch_to_cache(&tmpdir, "variation", &repeat_shift_cache_batch());
        write_batch_to_cache(&tmpdir, "transcript", &repeat_shift_transcripts_batch());
        write_batch_to_chrom(&tmpdir, "exon", "1", &repeat_shift_exons_batch());
        let cache_path = tmpdir.path().display().to_string();

        let expected_terms = ["frameshift_variant"];

        {
            let sql = format!(
                "SELECT \"CSQ\", most_severe_consequence \
                 FROM annotate_vep('vcf_repeat', '{cache_path}', 'parquet', '{{\"partitioned\":true}}')"
            );
            let batches = ctx
                .sql(&sql)
                .await
                .expect("query should parse")
                .collect()
                .await
                .expect("collect annotate_vep");
            let csq = string_values(batches[0].column_by_name("CSQ").expect("csq column exists"));
            let most = string_values(
                batches[0]
                    .column_by_name("most_severe_consequence")
                    .expect("most_severe_consequence column exists"),
            );
            let csq0 = csq[0].as_ref().expect("csq should be present");
            // variation_name is no longer emitted in CSQ (VEP Existing_variation is
            // empty for non-merged caches). The cache lookup still works; we just
            // don't put variation_name in the Allele|...|SOURCE format string.
            assert_terms_sorted_by_rank(csq0);
            assert_term_set_exact(csq0, &expected_terms);
            assert_eq!(most[0], Some("frameshift_variant".to_string()));
        }
    }

    #[tokio::test(flavor = "multi_thread")]
    async fn test_annotate_vep_translation_parity_handles_negative_strand_synonymous_case() {
        {
            let ctx = create_vep_session();
            ctx.register_table("vcf_neg", Arc::new(neg_strand_vcf_table()))
                .expect("register negative-strand vcf");
            let tmpdir = TempDir::new().expect("create tmpdir");
            write_batch_to_cache(&tmpdir, "variation", &neg_strand_cache_batch());
            write_batch_to_cache(&tmpdir, "transcript", &neg_strand_transcripts_batch());
            write_batch_to_chrom(&tmpdir, "exon", "1", &neg_strand_exons_batch());
            write_batch_to_chrom(
                &tmpdir,
                "translation_core",
                "1",
                &neg_strand_translations_batch(),
            );
            let cache_path = tmpdir.path().display().to_string();

            let sql = format!(
                "SELECT \"CSQ\", most_severe_consequence \
                 FROM annotate_vep('vcf_neg', '{cache_path}', 'parquet', '{{\"partitioned\":true}}')"
            );
            let df = ctx.sql(&sql).await.expect("query should parse");
            let batches = df.collect().await.expect("collect annotate_vep");
            let csq = string_values(batches[0].column_by_name("CSQ").expect("csq column exists"));
            let most = string_values(
                batches[0]
                    .column_by_name("most_severe_consequence")
                    .expect("most_severe_consequence column exists"),
            );
            let csq0 = csq[0].as_ref().expect("csq should be present");
            assert_terms_sorted_by_rank(csq0);
            assert_term_set_exact(csq0, &["synonymous_variant"]);
            assert_eq!(
                most[0],
                Some("synonymous_variant".to_string()),
                "backend=parquet"
            );
        }
    }

    #[tokio::test(flavor = "multi_thread")]
    async fn test_annotate_vep_translation_parity_handles_inframe_stop_loss_case() {
        {
            let ctx = create_vep_session();
            ctx.register_table("vcf_stop_loss", Arc::new(stop_loss_vcf_table()))
                .expect("register stop-loss vcf");
            let tmpdir = TempDir::new().expect("create tmpdir");
            write_batch_to_cache(&tmpdir, "variation", &stop_loss_cache_batch());
            write_batch_to_cache(&tmpdir, "transcript", &stop_loss_transcripts_batch());
            write_batch_to_chrom(&tmpdir, "exon", "1", &stop_loss_exons_batch());
            write_batch_to_chrom(
                &tmpdir,
                "translation_core",
                "1",
                &stop_loss_translations_batch(),
            );
            let cache_path = tmpdir.path().display().to_string();

            let sql = format!(
                "SELECT \"CSQ\", most_severe_consequence \
                 FROM annotate_vep('vcf_stop_loss', '{cache_path}', 'parquet', '{{\"partitioned\":true}}')"
            );
            let df = ctx.sql(&sql).await.expect("query should parse");
            let batches = df.collect().await.expect("collect annotate_vep");
            let csq = string_values(batches[0].column_by_name("CSQ").expect("csq column exists"));
            let most = string_values(
                batches[0]
                    .column_by_name("most_severe_consequence")
                    .expect("most_severe_consequence column exists"),
            );
            let csq0 = csq[0].as_ref().expect("csq should be present");
            assert_terms_sorted_by_rank(csq0);
            assert_term_set_exact(csq0, &["inframe_deletion"]);
            assert_eq!(
                most[0],
                Some("inframe_deletion".to_string()),
                "backend=parquet"
            );
        }
    }

    #[tokio::test(flavor = "multi_thread")]
    async fn test_annotate_vep_csq_has_74_pipe_delimited_fields_per_transcript() {
        {
            let ctx = create_vep_session();
            ctx.register_table("vcf_data", Arc::new(vcf_table()))
                .expect("register vcf table");
            let tmpdir = TempDir::new().expect("create tmpdir");
            write_batch_to_cache(&tmpdir, "variation", &cache_batch());
            write_batch_to_cache(&tmpdir, "transcript", &transcripts_batch());
            write_batch_to_chrom(&tmpdir, "exon", "1", &exons_batch());
            write_batch_to_chrom(&tmpdir, "exon", "2", &exons_batch());
            let cache_path = tmpdir.path().display().to_string();
            let backend = "parquet";

            let sql = format!(
                "SELECT \"CSQ\" FROM annotate_vep('vcf_data', '{cache_path}', 'parquet', '{{\"partitioned\":true}}') \
                 ORDER BY chrom"
            );
            let df = ctx.sql(&sql).await.expect("query should parse");
            let batches = df.collect().await.expect("collect annotate_vep");
            let csq_values = string_values(batches[0].column_by_name("CSQ").expect("csq column"));

            for csq_opt in &csq_values {
                let Some(csq) = csq_opt else { continue };
                for entry in csq.split(',') {
                    let fields: Vec<&str> = entry.split('|').collect();
                    assert_eq!(
                        fields.len(),
                        74,
                        "backend={backend}: CSQ entry should have 74 pipe-delimited fields, got {}: {:?}",
                        fields.len(),
                        entry
                    );
                    assert!(
                        !fields[0].is_empty(),
                        "backend={backend}: Allele field should not be empty"
                    );
                    assert!(
                        !fields[1].is_empty(),
                        "backend={backend}: Consequence field should not be empty"
                    );
                    assert!(
                        ["HIGH", "MODERATE", "LOW", "MODIFIER"].contains(&fields[2]),
                        "backend={backend}: IMPACT field should be a standard value, got: {}",
                        fields[2]
                    );
                    if fields[5] == "Transcript" {
                        assert!(
                            !fields[6].is_empty(),
                            "backend={backend}: Feature should have transcript_id when Feature_type is Transcript"
                        );
                    }
                }
            }
        }
    }

    #[tokio::test(flavor = "multi_thread")]
    async fn test_annotate_vep_refseq_and_merged_modes_emit_refseq_specific_csq_fields() {
        let ctx = create_vep_session();
        ctx.register_table("vcf_refseq", Arc::new(refseq_vcf_table()))
            .expect("register refseq vcf table");
        let tmpdir = TempDir::new().expect("create tmpdir");
        write_batch_to_cache(&tmpdir, "variation", &refseq_cache_batch());
        write_batch_to_cache(&tmpdir, "transcript", &refseq_transcripts_batch());
        write_batch_to_chrom(&tmpdir, "exon", "1", &refseq_exons_batch());
        let cache_path = tmpdir.path().display().to_string();

        let default_sql = format!(
            "SELECT \"CSQ\" FROM annotate_vep('vcf_refseq', '{cache_path}', 'parquet', '{{\"partitioned\":true}}')"
        );
        let default_batches = ctx
            .sql(&default_sql)
            .await
            .expect("default query should parse")
            .collect()
            .await
            .expect("collect default annotate_vep");
        let default_csq = string_values(
            default_batches[0]
                .column_by_name("CSQ")
                .expect("default csq column"),
        );
        let default_entry = csq_entries(default_csq[0].as_deref().expect("default csq present"))
            .into_iter()
            .next()
            .expect("default CSQ entry");
        assert_eq!(default_entry.len(), 74);
        assert_eq!(default_entry[1], "intergenic_variant");

        let refseq_sql = format!(
            "SELECT \"CSQ\" FROM annotate_vep('vcf_refseq', '{cache_path}', 'parquet', '{{\"partitioned\":true,\"refseq\":true}}')"
        );
        let refseq_batches = ctx
            .sql(&refseq_sql)
            .await
            .expect("refseq query should parse")
            .collect()
            .await
            .expect("collect refseq annotate_vep");
        let refseq_csq = string_values(
            refseq_batches[0]
                .column_by_name("CSQ")
                .expect("refseq csq column"),
        );
        let refseq_entry = csq_entries(refseq_csq[0].as_deref().expect("refseq csq present"))
            .into_iter()
            .find(|fields| fields.len() == 78 && fields[5] == "Transcript")
            .expect("expected transcript CSQ entry in refseq mode");
        assert_eq!(refseq_entry[6], "NM_000001");
        assert_eq!(refseq_entry[28], "rseq_ens_match_cds");
        assert_eq!(refseq_entry[29], "");
        assert_eq!(refseq_entry[30], "A");
        assert_eq!(refseq_entry[31], "A");
        assert_eq!(refseq_entry[32], "OK");

        let merged_sql = format!(
            "SELECT \"CSQ\" FROM annotate_vep('vcf_refseq', '{cache_path}', 'parquet', '{{\"partitioned\":true,\"merged\":true}}')"
        );
        let merged_batches = ctx
            .sql(&merged_sql)
            .await
            .expect("merged query should parse")
            .collect()
            .await
            .expect("collect merged annotate_vep");
        let merged_csq = string_values(
            merged_batches[0]
                .column_by_name("CSQ")
                .expect("merged csq column"),
        );
        let merged_entry = csq_entries(merged_csq[0].as_deref().expect("merged csq present"))
            .into_iter()
            .find(|fields| fields.len() == 79 && fields[5] == "Transcript")
            .expect("expected transcript CSQ entry in merged mode");
        assert_eq!(merged_entry[6], "NM_000001");
        assert_eq!(merged_entry[28], "rseq_ens_match_cds");
        assert_eq!(merged_entry[29], "RefSeq");
        assert_eq!(merged_entry[30], "");
        assert_eq!(merged_entry[31], "A");
        assert_eq!(merged_entry[32], "A");
        assert_eq!(merged_entry[33], "OK");
    }

    #[tokio::test(flavor = "multi_thread")]
    async fn test_annotate_vep_rejects_unknown_backend() {
        let ctx = create_vep_session();
        ctx.register_table("vcf_data", Arc::new(vcf_table()))
            .expect("register vcf table");

        let err = ctx
            .sql("SELECT * FROM annotate_vep('vcf_data', '/tmp/vep_cache', 'bad_backend')")
            .await
            .expect_err("unknown backend should fail")
            .to_string();

        assert!(err.contains("annotate_vep() backend must be one of"));
    }

    /// VCF table with 10 rows on the same contig for LIMIT tests.
    fn multi_row_vcf_table() -> MemTable {
        let schema = Arc::new(Schema::new(vec![
            Field::new("chrom", DataType::Utf8, false),
            Field::new("start", DataType::Int64, false),
            Field::new("end", DataType::Int64, false),
            Field::new("ref", DataType::Utf8, false),
            Field::new("alt", DataType::Utf8, false),
        ]));
        let n = 10;
        let batch = RecordBatch::try_new(
            schema.clone(),
            vec![
                Arc::new(StringArray::from(vec!["1"; n])),
                Arc::new(Int64Array::from(
                    (0..n as i64).map(|i| 100 + i * 10).collect::<Vec<_>>(),
                )),
                Arc::new(Int64Array::from(
                    (0..n as i64).map(|i| 101 + i * 10).collect::<Vec<_>>(),
                )),
                Arc::new(StringArray::from(vec!["A"; n])),
                Arc::new(StringArray::from(vec!["G"; n])),
            ],
        )
        .expect("valid multi-row vcf batch");
        MemTable::try_new(schema, vec![vec![batch]]).expect("valid multi-row vcf memtable")
    }

    /// Variation cache with 10 rows matching multi_row_vcf_table.
    fn multi_row_cache_batch() -> RecordBatch {
        use crate::annotate_provider::cache_lookup_column_names;

        let n = 10;
        let mut fields = vec![
            Field::new("chrom", DataType::Utf8, false),
            Field::new("start", DataType::Int64, false),
            Field::new("end", DataType::Int64, false),
            Field::new("allele_string", DataType::Utf8, false),
            Field::new("failed", DataType::Int64, false),
        ];
        let mut columns: Vec<Arc<dyn datafusion::arrow::array::Array>> = vec![
            Arc::new(StringArray::from(vec!["1"; n])),
            Arc::new(Int64Array::from(
                (0..n as i64).map(|i| 100 + i * 10).collect::<Vec<_>>(),
            )),
            Arc::new(Int64Array::from(
                (0..n as i64).map(|i| 101 + i * 10).collect::<Vec<_>>(),
            )),
            Arc::new(StringArray::from(vec!["A/G"; n])),
            Arc::new(Int64Array::from(vec![0_i64; n])),
        ];
        for col in cache_lookup_column_names() {
            fields.push(Field::new(col, DataType::Utf8, true));
            columns.push(Arc::new(StringArray::from(vec![None::<&str>; n])));
        }
        let schema = Arc::new(Schema::new(fields));
        RecordBatch::try_new(schema, columns).expect("valid multi-row cache batch")
    }

    #[tokio::test(flavor = "multi_thread")]
    async fn test_annotate_vep_limit_returns_exact_count() {
        let tmpdir = TempDir::new().expect("create temp dir");
        write_batch_to_cache(&tmpdir, "variation", &multi_row_cache_batch());

        let ctx = create_vep_session();
        ctx.register_table("vcf_data", Arc::new(multi_row_vcf_table()))
            .expect("register vcf table");

        let cache_path = tmpdir.path().to_str().expect("utf8 path");
        let sql = format!(
            "SELECT * FROM annotate_vep('vcf_data', '{cache_path}', 'parquet', \
             '{{\"partitioned\":true}}') LIMIT 3"
        );
        let batches = ctx
            .sql(&sql)
            .await
            .expect("annotate_vep with LIMIT")
            .collect()
            .await
            .expect("collect");

        let total_rows: usize = batches.iter().map(|b| b.num_rows()).sum();
        assert_eq!(total_rows, 3, "LIMIT 3 should return exactly 3 rows");
    }

    #[tokio::test(flavor = "multi_thread")]
    async fn test_annotate_vep_limit_one() {
        let tmpdir = TempDir::new().expect("create temp dir");
        write_batch_to_cache(&tmpdir, "variation", &multi_row_cache_batch());

        let ctx = create_vep_session();
        ctx.register_table("vcf_data", Arc::new(multi_row_vcf_table()))
            .expect("register vcf table");

        let cache_path = tmpdir.path().to_str().expect("utf8 path");
        let sql = format!(
            "SELECT * FROM annotate_vep('vcf_data', '{cache_path}', 'parquet', \
             '{{\"partitioned\":true}}') LIMIT 1"
        );
        let batches = ctx
            .sql(&sql)
            .await
            .expect("annotate_vep with LIMIT 1")
            .collect()
            .await
            .expect("collect");

        let total_rows: usize = batches.iter().map(|b| b.num_rows()).sum();
        assert_eq!(total_rows, 1, "LIMIT 1 should return exactly 1 row");
    }

    #[tokio::test(flavor = "multi_thread")]
    async fn test_annotate_vep_limit_exceeds_total_returns_all() {
        let tmpdir = TempDir::new().expect("create temp dir");
        write_batch_to_cache(&tmpdir, "variation", &multi_row_cache_batch());

        let ctx = create_vep_session();
        ctx.register_table("vcf_data", Arc::new(multi_row_vcf_table()))
            .expect("register vcf table");

        let cache_path = tmpdir.path().to_str().expect("utf8 path");
        let sql = format!(
            "SELECT * FROM annotate_vep('vcf_data', '{cache_path}', 'parquet', \
             '{{\"partitioned\":true}}') LIMIT 100"
        );
        let batches = ctx
            .sql(&sql)
            .await
            .expect("annotate_vep with LIMIT 100")
            .collect()
            .await
            .expect("collect");

        let total_rows: usize = batches.iter().map(|b| b.num_rows()).sum();
        assert_eq!(
            total_rows, 10,
            "LIMIT 100 with 10 input rows should return all 10"
        );
    }

    #[tokio::test(flavor = "multi_thread")]
    async fn test_annotate_vep_no_limit_returns_all() {
        let tmpdir = TempDir::new().expect("create temp dir");
        write_batch_to_cache(&tmpdir, "variation", &multi_row_cache_batch());

        let ctx = create_vep_session();
        ctx.register_table("vcf_data", Arc::new(multi_row_vcf_table()))
            .expect("register vcf table");

        let cache_path = tmpdir.path().to_str().expect("utf8 path");
        let sql = format!(
            "SELECT * FROM annotate_vep('vcf_data', '{cache_path}', 'parquet', \
             '{{\"partitioned\":true}}')"
        );
        let batches = ctx
            .sql(&sql)
            .await
            .expect("annotate_vep without LIMIT")
            .collect()
            .await
            .expect("collect");

        let total_rows: usize = batches.iter().map(|b| b.num_rows()).sum();
        assert_eq!(total_rows, 10, "No LIMIT should return all 10 rows");
    }

    #[cfg(feature = "kv-cache")]
    #[tokio::test(flavor = "multi_thread")]
    async fn test_annotate_vep_limit_with_fjall_releases_lock() {
        let tmpdir = TempDir::new().expect("create temp dir");
        let cache_batch = multi_row_cache_batch();
        write_batch_to_cache(&tmpdir, "variation", &cache_batch);
        write_batch_to_fjall(&tmpdir, &cache_batch);

        let ctx = create_vep_session();
        ctx.register_table("vcf_data", Arc::new(multi_row_vcf_table()))
            .expect("register vcf table");

        let cache_path = tmpdir.path().to_str().expect("utf8 path");
        let sql = format!(
            "SELECT * FROM annotate_vep('vcf_data', '{cache_path}', 'parquet', \
             '{{\"partitioned\":true,\"use_fjall\":true}}') LIMIT 3"
        );
        let batches = ctx
            .sql(&sql)
            .await
            .expect("annotate_vep with LIMIT 3 and fjall")
            .collect()
            .await
            .expect("collect");

        let total_rows: usize = batches.iter().map(|b| b.num_rows()).sum();
        assert_eq!(total_rows, 3, "LIMIT 3 should return exactly 3 rows");

        let reopen_err = VepKvStore::open(tmpdir.path().join("variation.fjall"))
            .err()
            .map(|e| e.to_string());
        assert!(
            reopen_err.is_none(),
            "fjall cache should reopen after LIMIT without lingering lock: {}",
            reopen_err.unwrap_or_default()
        );
    }
}
