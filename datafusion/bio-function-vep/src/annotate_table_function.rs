//! Table function registration for consequence annotation.
//!
//! `annotate_vep()` is the high-level consequence annotation entrypoint.

use std::sync::Arc;

use datafusion::arrow::datatypes::Schema;
use datafusion::catalog::TableFunctionImpl;
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
}

impl AnnotateFunction {
    pub fn new(session: Arc<SessionContext>) -> Self {
        Self { session }
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

        let (vcf_schema, _) = resolve_schema(&self.session, &vcf_table)?;

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

fn resolve_schema(session: &SessionContext, vcf_table: &str) -> Result<(Schema, String)> {
    match tokio::runtime::Handle::try_current() {
        Ok(handle) => tokio::task::block_in_place(|| {
            let vcf = handle.block_on(session.table(vcf_table))?;
            Ok::<_, DataFusionError>((vcf.schema().as_arrow().clone(), vcf_table.to_string()))
        }),
        Err(_) => {
            let rt = tokio::runtime::Runtime::new()
                .map_err(|e| DataFusionError::External(Box::new(e)))?;
            let vcf = rt.block_on(session.table(vcf_table))?;
            Ok((vcf.schema().as_arrow().clone(), vcf_table.to_string()))
        }
    }
}

#[cfg(test)]
mod tests {
    use crate::create_vep_session;
    use crate::so_terms::SoTerm;
    use datafusion::arrow::array::{Array, Float64Array, Int64Array, RecordBatch, StringArray};
    use datafusion::arrow::datatypes::{DataType, Field, Schema};
    use datafusion::datasource::MemTable;
    use datafusion::prelude::ParquetReadOptions;
    use parquet::arrow::ArrowWriter;
    use std::collections::BTreeSet;
    use std::fs::File;
    use std::sync::Arc;
    use tempfile::TempDir;

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

    fn cache_table() -> MemTable {
        use crate::annotate_provider::CACHE_OUTPUT_COLUMNS;

        // Core columns required for the join.
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
            Arc::new(Int64Array::from(vec![0])),
        ];
        // Add all CACHE_OUTPUT_COLUMNS as nullable Utf8.
        for &col in CACHE_OUTPUT_COLUMNS {
            fields.push(Field::new(col, DataType::Utf8, true));
            let val: Option<&str> = match col {
                "variation_name" => Some("rs100"),
                "clin_sig" => Some("benign"),
                _ => None,
            };
            columns.push(Arc::new(StringArray::from(vec![val])));
        }
        let schema = Arc::new(Schema::new(fields));
        let batch = RecordBatch::try_new(schema.clone(), columns).expect("valid cache batch");
        MemTable::try_new(schema, vec![vec![batch]]).expect("valid cache memtable")
    }

    fn transcripts_table() -> MemTable {
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
        let batch = RecordBatch::try_new(
            schema.clone(),
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
        .expect("valid transcript batch");
        MemTable::try_new(schema, vec![vec![batch]]).expect("valid transcript memtable")
    }

    fn exons_table() -> MemTable {
        let schema = Arc::new(Schema::new(vec![
            Field::new("transcript_id", DataType::Utf8, false),
            Field::new("exon_number", DataType::Int64, false),
            Field::new("start", DataType::Int64, false),
            Field::new("end", DataType::Int64, false),
        ]));
        let batch = RecordBatch::try_new(
            schema.clone(),
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
        .expect("valid exon batch");
        MemTable::try_new(schema, vec![vec![batch]]).expect("valid exon memtable")
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

    fn synonymous_cache_table() -> MemTable {
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
        let batch = RecordBatch::try_new(
            schema.clone(),
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
        .expect("valid synonymous cache batch");
        MemTable::try_new(schema, vec![vec![batch]]).expect("valid synonymous cache memtable")
    }

    fn synonymous_transcripts_table() -> MemTable {
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
        let batch = RecordBatch::try_new(
            schema.clone(),
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
        .expect("valid synonymous transcript batch");
        MemTable::try_new(schema, vec![vec![batch]]).expect("valid synonymous transcript memtable")
    }

    fn synonymous_exons_table() -> MemTable {
        let schema = Arc::new(Schema::new(vec![
            Field::new("transcript_id", DataType::Utf8, false),
            Field::new("exon_number", DataType::Int64, false),
            Field::new("start", DataType::Int64, false),
            Field::new("end", DataType::Int64, false),
        ]));
        let batch = RecordBatch::try_new(
            schema.clone(),
            vec![
                Arc::new(StringArray::from(vec!["ENST00000000001"])),
                Arc::new(Int64Array::from(vec![1])),
                Arc::new(Int64Array::from(vec![90])),
                Arc::new(Int64Array::from(vec![140])),
            ],
        )
        .expect("valid synonymous exon batch");
        MemTable::try_new(schema, vec![vec![batch]]).expect("valid synonymous exon memtable")
    }

    fn synonymous_translations_table() -> MemTable {
        let schema = Arc::new(Schema::new(vec![
            Field::new("transcript_id", DataType::Utf8, false),
            Field::new("cds_len", DataType::Int64, true),
            Field::new("protein_len", DataType::Int64, true),
            Field::new("translation_seq", DataType::Utf8, true),
            Field::new("cds_sequence", DataType::Utf8, true),
        ]));
        // CDS from 97..105 => ATG GCT TAA
        let batch = RecordBatch::try_new(
            schema.clone(),
            vec![
                Arc::new(StringArray::from(vec!["ENST00000000001"])),
                Arc::new(Int64Array::from(vec![9])),
                Arc::new(Int64Array::from(vec![3])),
                Arc::new(StringArray::from(vec!["MA*"])),
                Arc::new(StringArray::from(vec!["ATGGCTTAA"])),
            ],
        )
        .expect("valid synonymous translations batch");
        MemTable::try_new(schema, vec![vec![batch]])
            .expect("valid synonymous translations memtable")
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

    fn context_cache_table() -> MemTable {
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
        let batch = RecordBatch::try_new(
            schema.clone(),
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
        .expect("valid context cache batch");
        MemTable::try_new(schema, vec![vec![batch]]).expect("valid context cache memtable")
    }

    fn context_transcripts_table() -> MemTable {
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
        let batch = RecordBatch::try_new(
            schema.clone(),
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
        .expect("valid context transcript batch");
        MemTable::try_new(schema, vec![vec![batch]]).expect("valid context transcript memtable")
    }

    fn context_exons_table() -> MemTable {
        let schema = Arc::new(Schema::new(vec![
            Field::new("transcript_id", DataType::Utf8, false),
            Field::new("exon_number", DataType::Int64, false),
            Field::new("start", DataType::Int64, false),
            Field::new("end", DataType::Int64, false),
        ]));
        let batch = RecordBatch::try_new(
            schema.clone(),
            vec![
                Arc::new(StringArray::from(vec!["ENST00000000002"])),
                Arc::new(Int64Array::from(vec![1])),
                Arc::new(Int64Array::from(vec![100])),
                Arc::new(Int64Array::from(vec![250])),
            ],
        )
        .expect("valid context exon batch");
        MemTable::try_new(schema, vec![vec![batch]]).expect("valid context exon memtable")
    }

    fn context_regulatory_table() -> MemTable {
        let schema = Arc::new(Schema::new(vec![
            Field::new("stable_id", DataType::Utf8, false),
            Field::new("chrom", DataType::Utf8, false),
            Field::new("start", DataType::Int64, false),
            Field::new("end", DataType::Int64, false),
        ]));
        let batch = RecordBatch::try_new(
            schema.clone(),
            vec![
                Arc::new(StringArray::from(vec!["reg_ctx"])),
                Arc::new(StringArray::from(vec!["1"])),
                Arc::new(Int64Array::from(vec![150])),
                Arc::new(Int64Array::from(vec![160])),
            ],
        )
        .expect("valid context regulatory batch");
        MemTable::try_new(schema, vec![vec![batch]]).expect("valid context regulatory memtable")
    }

    fn context_motif_table() -> MemTable {
        let schema = Arc::new(Schema::new(vec![
            Field::new("motif_id", DataType::Utf8, false),
            Field::new("chrom", DataType::Utf8, false),
            Field::new("start", DataType::Int64, false),
            Field::new("end", DataType::Int64, false),
        ]));
        let batch = RecordBatch::try_new(
            schema.clone(),
            vec![
                Arc::new(StringArray::from(vec!["motif_ctx"])),
                Arc::new(StringArray::from(vec!["1"])),
                Arc::new(Int64Array::from(vec![150])),
                Arc::new(Int64Array::from(vec![160])),
            ],
        )
        .expect("valid context motif batch");
        MemTable::try_new(schema, vec![vec![batch]]).expect("valid context motif memtable")
    }

    fn context_mirna_table() -> MemTable {
        let schema = Arc::new(Schema::new(vec![
            Field::new("mirna_id", DataType::Utf8, false),
            Field::new("chrom", DataType::Utf8, false),
            Field::new("start", DataType::Int64, false),
            Field::new("end", DataType::Int64, false),
        ]));
        let batch = RecordBatch::try_new(
            schema.clone(),
            vec![
                Arc::new(StringArray::from(vec!["mir_ctx"])),
                Arc::new(StringArray::from(vec!["1"])),
                Arc::new(Int64Array::from(vec![150])),
                Arc::new(Int64Array::from(vec![160])),
            ],
        )
        .expect("valid context miRNA batch");
        MemTable::try_new(schema, vec![vec![batch]]).expect("valid context miRNA memtable")
    }

    fn context_sv_table() -> MemTable {
        let schema = Arc::new(Schema::new(vec![
            Field::new("feature_id", DataType::Utf8, false),
            Field::new("feature_kind", DataType::Utf8, false),
            Field::new("event_type", DataType::Utf8, false),
            Field::new("chrom", DataType::Utf8, false),
            Field::new("start", DataType::Int64, false),
            Field::new("end", DataType::Int64, false),
        ]));
        let batch = RecordBatch::try_new(
            schema.clone(),
            vec![
                Arc::new(StringArray::from(vec![
                    "sv_tx_ablation",
                    "sv_reg_amp",
                    "sv_tfbs_ablation",
                ])),
                Arc::new(StringArray::from(vec!["transcript", "regulatory", "tfbs"])),
                Arc::new(StringArray::from(vec![
                    "ablation",
                    "amplification",
                    "ablation",
                ])),
                Arc::new(StringArray::from(vec!["1", "1", "1"])),
                Arc::new(Int64Array::from(vec![150, 150, 150])),
                Arc::new(Int64Array::from(vec![160, 160, 160])),
            ],
        )
        .expect("valid context sv batch");
        MemTable::try_new(schema, vec![vec![batch]]).expect("valid context sv memtable")
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

    fn golden_context_cache_table() -> MemTable {
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
        let batch = RecordBatch::try_new(
            schema.clone(),
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
        .expect("valid golden context cache batch");
        MemTable::try_new(schema, vec![vec![batch]]).expect("valid golden context cache memtable")
    }

    fn golden_context_transcripts_table() -> MemTable {
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
        let batch = RecordBatch::try_new(
            schema.clone(),
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
        .expect("valid golden context transcript batch");
        MemTable::try_new(schema, vec![vec![batch]])
            .expect("valid golden context transcript memtable")
    }

    fn golden_context_exons_table() -> MemTable {
        let schema = Arc::new(Schema::new(vec![
            Field::new("transcript_id", DataType::Utf8, false),
            Field::new("exon_number", DataType::Int64, false),
            Field::new("start", DataType::Int64, false),
            Field::new("end", DataType::Int64, false),
        ]));
        let batch = RecordBatch::try_new(
            schema.clone(),
            vec![
                Arc::new(StringArray::from(vec!["ENST00000000003"])),
                Arc::new(Int64Array::from(vec![1])),
                Arc::new(Int64Array::from(vec![20_000])),
                Arc::new(Int64Array::from(vec![20_100])),
            ],
        )
        .expect("valid golden context exon batch");
        MemTable::try_new(schema, vec![vec![batch]]).expect("valid golden context exon memtable")
    }

    fn golden_context_regulatory_table() -> MemTable {
        let schema = Arc::new(Schema::new(vec![
            Field::new("stable_id", DataType::Utf8, false),
            Field::new("chrom", DataType::Utf8, false),
            Field::new("start", DataType::Int64, false),
            Field::new("end", DataType::Int64, false),
        ]));
        let batch = RecordBatch::try_new(
            schema.clone(),
            vec![
                Arc::new(StringArray::from(vec!["reg_golden"])),
                Arc::new(StringArray::from(vec!["1"])),
                Arc::new(Int64Array::from(vec![495])),
                Arc::new(Int64Array::from(vec![505])),
            ],
        )
        .expect("valid golden context regulatory batch");
        MemTable::try_new(schema, vec![vec![batch]])
            .expect("valid golden context regulatory memtable")
    }

    fn golden_context_motif_table() -> MemTable {
        let schema = Arc::new(Schema::new(vec![
            Field::new("motif_id", DataType::Utf8, false),
            Field::new("chrom", DataType::Utf8, false),
            Field::new("start", DataType::Int64, false),
            Field::new("end", DataType::Int64, false),
        ]));
        let batch = RecordBatch::try_new(
            schema.clone(),
            vec![
                Arc::new(StringArray::from(vec!["motif_golden"])),
                Arc::new(StringArray::from(vec!["1"])),
                Arc::new(Int64Array::from(vec![495])),
                Arc::new(Int64Array::from(vec![505])),
            ],
        )
        .expect("valid golden context motif batch");
        MemTable::try_new(schema, vec![vec![batch]]).expect("valid golden context motif memtable")
    }

    fn golden_context_mirna_table() -> MemTable {
        let schema = Arc::new(Schema::new(vec![
            Field::new("mirna_id", DataType::Utf8, false),
            Field::new("chrom", DataType::Utf8, false),
            Field::new("start", DataType::Int64, false),
            Field::new("end", DataType::Int64, false),
        ]));
        let batch = RecordBatch::try_new(
            schema.clone(),
            vec![
                Arc::new(StringArray::from(vec!["mirna_golden"])),
                Arc::new(StringArray::from(vec!["1"])),
                Arc::new(Int64Array::from(vec![495])),
                Arc::new(Int64Array::from(vec![505])),
            ],
        )
        .expect("valid golden context miRNA batch");
        MemTable::try_new(schema, vec![vec![batch]]).expect("valid golden context miRNA memtable")
    }

    fn golden_context_sv_table() -> MemTable {
        let schema = Arc::new(Schema::new(vec![
            Field::new("feature_id", DataType::Utf8, false),
            Field::new("feature_kind", DataType::Utf8, false),
            Field::new("event_type", DataType::Utf8, false),
            Field::new("chrom", DataType::Utf8, false),
            Field::new("start", DataType::Int64, false),
            Field::new("end", DataType::Int64, false),
        ]));
        let batch = RecordBatch::try_new(
            schema.clone(),
            vec![
                Arc::new(StringArray::from(vec![
                    "sv_tx_ablation",
                    "sv_tx_amplification",
                    "sv_tx_elongation",
                    "sv_tx_truncation",
                    "sv_reg_ablation",
                    "sv_reg_amplification",
                    "sv_tfbs_ablation",
                    "sv_tfbs_amplification",
                ])),
                Arc::new(StringArray::from(vec![
                    "transcript",
                    "transcript",
                    "transcript",
                    "transcript",
                    "regulatory",
                    "regulatory",
                    "tfbs",
                    "tfbs",
                ])),
                Arc::new(StringArray::from(vec![
                    "ablation",
                    "amplification",
                    "elongation",
                    "truncation",
                    "ablation",
                    "amplification",
                    "ablation",
                    "amplification",
                ])),
                Arc::new(StringArray::from(vec![
                    "1", "1", "1", "1", "1", "1", "1", "1",
                ])),
                Arc::new(Int64Array::from(vec![
                    495, 495, 495, 495, 495, 495, 495, 495,
                ])),
                Arc::new(Int64Array::from(vec![
                    505, 505, 505, 505, 505, 505, 505, 505,
                ])),
            ],
        )
        .expect("valid golden context sv batch");
        MemTable::try_new(schema, vec![vec![batch]]).expect("valid golden context sv memtable")
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

    fn splice_cache_table() -> MemTable {
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
        let batch = RecordBatch::try_new(
            schema.clone(),
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
        .expect("valid splice cache batch");
        MemTable::try_new(schema, vec![vec![batch]]).expect("valid splice cache memtable")
    }

    fn splice_transcripts_table() -> MemTable {
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
        let batch = RecordBatch::try_new(
            schema.clone(),
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
        .expect("valid splice transcript batch");
        MemTable::try_new(schema, vec![vec![batch]]).expect("valid splice transcript memtable")
    }

    fn splice_exons_table() -> MemTable {
        let schema = Arc::new(Schema::new(vec![
            Field::new("transcript_id", DataType::Utf8, false),
            Field::new("exon_number", DataType::Int64, false),
            Field::new("start", DataType::Int64, false),
            Field::new("end", DataType::Int64, false),
        ]));
        let batch = RecordBatch::try_new(
            schema.clone(),
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
        .expect("valid splice exon batch");
        MemTable::try_new(schema, vec![vec![batch]]).expect("valid splice exon memtable")
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

    fn repeat_shift_cache_table() -> MemTable {
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
        let batch = RecordBatch::try_new(
            schema.clone(),
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
        .expect("valid repeat-shift cache batch");
        MemTable::try_new(schema, vec![vec![batch]]).expect("valid repeat-shift cache memtable")
    }

    fn repeat_shift_transcripts_table() -> MemTable {
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
        let batch = RecordBatch::try_new(
            schema.clone(),
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
        .expect("valid repeat-shift transcript batch");
        MemTable::try_new(schema, vec![vec![batch]])
            .expect("valid repeat-shift transcript memtable")
    }

    fn repeat_shift_exons_table() -> MemTable {
        let schema = Arc::new(Schema::new(vec![
            Field::new("transcript_id", DataType::Utf8, false),
            Field::new("exon_number", DataType::Int64, false),
            Field::new("start", DataType::Int64, false),
            Field::new("end", DataType::Int64, false),
        ]));
        let batch = RecordBatch::try_new(
            schema.clone(),
            vec![
                Arc::new(StringArray::from(vec!["ENST00000000005"])),
                Arc::new(Int64Array::from(vec![1])),
                Arc::new(Int64Array::from(vec![40])),
                Arc::new(Int64Array::from(vec![120])),
            ],
        )
        .expect("valid repeat-shift exon batch");
        MemTable::try_new(schema, vec![vec![batch]]).expect("valid repeat-shift exon memtable")
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

    fn neg_strand_cache_table() -> MemTable {
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
        let batch = RecordBatch::try_new(
            schema.clone(),
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
        .expect("valid negative-strand cache batch");
        MemTable::try_new(schema, vec![vec![batch]]).expect("valid negative-strand cache memtable")
    }

    fn neg_strand_transcripts_table() -> MemTable {
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
        let batch = RecordBatch::try_new(
            schema.clone(),
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
        .expect("valid negative-strand transcript batch");
        MemTable::try_new(schema, vec![vec![batch]])
            .expect("valid negative-strand transcript memtable")
    }

    fn neg_strand_exons_table() -> MemTable {
        let schema = Arc::new(Schema::new(vec![
            Field::new("transcript_id", DataType::Utf8, false),
            Field::new("exon_number", DataType::Int64, false),
            Field::new("start", DataType::Int64, false),
            Field::new("end", DataType::Int64, false),
        ]));
        let batch = RecordBatch::try_new(
            schema.clone(),
            vec![
                Arc::new(StringArray::from(vec!["ENST00000000006"])),
                Arc::new(Int64Array::from(vec![1])),
                Arc::new(Int64Array::from(vec![90])),
                Arc::new(Int64Array::from(vec![140])),
            ],
        )
        .expect("valid negative-strand exon batch");
        MemTable::try_new(schema, vec![vec![batch]]).expect("valid negative-strand exon memtable")
    }

    fn neg_strand_translations_table() -> MemTable {
        let schema = Arc::new(Schema::new(vec![
            Field::new("transcript_id", DataType::Utf8, false),
            Field::new("cds_len", DataType::Int64, true),
            Field::new("protein_len", DataType::Int64, true),
            Field::new("translation_seq", DataType::Utf8, true),
            Field::new("cds_sequence", DataType::Utf8, true),
        ]));
        let batch = RecordBatch::try_new(
            schema.clone(),
            vec![
                Arc::new(StringArray::from(vec!["ENST00000000006"])),
                Arc::new(Int64Array::from(vec![9])),
                Arc::new(Int64Array::from(vec![3])),
                Arc::new(StringArray::from(vec!["MA*"])),
                Arc::new(StringArray::from(vec!["ATGGCTTAA"])),
            ],
        )
        .expect("valid negative-strand translation batch");
        MemTable::try_new(schema, vec![vec![batch]])
            .expect("valid negative-strand translation memtable")
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

    fn stop_loss_cache_table() -> MemTable {
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
        let batch = RecordBatch::try_new(
            schema.clone(),
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
        .expect("valid stop-loss cache batch");
        MemTable::try_new(schema, vec![vec![batch]]).expect("valid stop-loss cache memtable")
    }

    fn stop_loss_transcripts_table() -> MemTable {
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
        let batch = RecordBatch::try_new(
            schema.clone(),
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
        .expect("valid stop-loss transcript batch");
        MemTable::try_new(schema, vec![vec![batch]]).expect("valid stop-loss transcript memtable")
    }

    fn stop_loss_exons_table() -> MemTable {
        let schema = Arc::new(Schema::new(vec![
            Field::new("transcript_id", DataType::Utf8, false),
            Field::new("exon_number", DataType::Int64, false),
            Field::new("start", DataType::Int64, false),
            Field::new("end", DataType::Int64, false),
        ]));
        let batch = RecordBatch::try_new(
            schema.clone(),
            vec![
                Arc::new(StringArray::from(vec!["ENST00000000007"])),
                Arc::new(Int64Array::from(vec![1])),
                Arc::new(Int64Array::from(vec![90])),
                Arc::new(Int64Array::from(vec![140])),
            ],
        )
        .expect("valid stop-loss exon batch");
        MemTable::try_new(schema, vec![vec![batch]]).expect("valid stop-loss exon memtable")
    }

    fn stop_loss_translations_table() -> MemTable {
        let schema = Arc::new(Schema::new(vec![
            Field::new("transcript_id", DataType::Utf8, false),
            Field::new("cds_len", DataType::Int64, true),
            Field::new("protein_len", DataType::Int64, true),
            Field::new("translation_seq", DataType::Utf8, true),
            Field::new("cds_sequence", DataType::Utf8, true),
        ]));
        let batch = RecordBatch::try_new(
            schema.clone(),
            vec![
                Arc::new(StringArray::from(vec!["ENST00000000007"])),
                Arc::new(Int64Array::from(vec![9])),
                Arc::new(Int64Array::from(vec![3])),
                Arc::new(StringArray::from(vec!["MA*"])),
                Arc::new(StringArray::from(vec!["ATGGCTTAA"])),
            ],
        )
        .expect("valid stop-loss translation batch");
        MemTable::try_new(schema, vec![vec![batch]]).expect("valid stop-loss translation memtable")
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

    fn hgvs_cache_table() -> MemTable {
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
        let batch = RecordBatch::try_new(
            schema.clone(),
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
        .expect("valid hgvs cache batch");
        MemTable::try_new(schema, vec![vec![batch]]).expect("valid hgvs cache memtable")
    }

    fn hgvs_transcripts_table() -> MemTable {
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
        let batch = RecordBatch::try_new(
            schema.clone(),
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
        .expect("valid hgvs transcript batch");
        MemTable::try_new(schema, vec![vec![batch]]).expect("valid hgvs transcript memtable")
    }

    fn hgvs_exons_table() -> MemTable {
        let schema = Arc::new(Schema::new(vec![
            Field::new("transcript_id", DataType::Utf8, false),
            Field::new("exon_number", DataType::Int64, false),
            Field::new("start", DataType::Int64, false),
            Field::new("end", DataType::Int64, false),
        ]));
        let batch = RecordBatch::try_new(
            schema.clone(),
            vec![
                Arc::new(StringArray::from(vec!["ENSTHGVS000001"])),
                Arc::new(Int64Array::from(vec![1])),
                Arc::new(Int64Array::from(vec![90])),
                Arc::new(Int64Array::from(vec![140])),
            ],
        )
        .expect("valid hgvs exon batch");
        MemTable::try_new(schema, vec![vec![batch]]).expect("valid hgvs exon memtable")
    }

    fn hgvs_translations_table() -> MemTable {
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
        let batch = RecordBatch::try_new(
            schema.clone(),
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
        .expect("valid hgvs translation batch");
        MemTable::try_new(schema, vec![vec![batch]]).expect("valid hgvs translation memtable")
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

    fn distance_cache_table() -> MemTable {
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
        let batch = RecordBatch::try_new(
            schema.clone(),
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
        .expect("valid distance cache batch");
        MemTable::try_new(schema, vec![vec![batch]]).expect("valid distance cache memtable")
    }

    fn distance_transcripts_table() -> MemTable {
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
        let batch = RecordBatch::try_new(
            schema.clone(),
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
        .expect("valid distance transcript batch");
        MemTable::try_new(schema, vec![vec![batch]]).expect("valid distance transcript memtable")
    }

    fn distance_exons_table() -> MemTable {
        let schema = Arc::new(Schema::new(vec![
            Field::new("transcript_id", DataType::Utf8, false),
            Field::new("exon_number", DataType::Int64, false),
            Field::new("start", DataType::Int64, false),
            Field::new("end", DataType::Int64, false),
        ]));
        let batch = RecordBatch::try_new(
            schema.clone(),
            vec![
                Arc::new(StringArray::from(vec!["ENSTDISTUP", "ENSTDISTDOWN"])),
                Arc::new(Int64Array::from(vec![1, 1])),
                Arc::new(Int64Array::from(vec![37_079, 14_000])),
                Arc::new(Int64Array::from(vec![37_200, 15_000])),
            ],
        )
        .expect("valid distance exon batch");
        MemTable::try_new(schema, vec![vec![batch]]).expect("valid distance exon memtable")
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

    fn regulatory_table_with_rows(rows: &[(&str, &str, i64, i64, Option<&str>)]) -> MemTable {
        let batch = regulatory_feature_batch(rows);
        MemTable::try_new(batch.schema(), vec![vec![batch]]).expect("valid regulatory memtable")
    }

    async fn register_single_batch_parquet(
        ctx: &datafusion::prelude::SessionContext,
        table_name: &str,
        batch: &RecordBatch,
    ) -> TempDir {
        let dir = tempfile::tempdir().expect("create temp parquet dir");
        let path = dir.path().join(format!("{table_name}.parquet"));
        let file = File::create(&path).expect("create parquet file");
        let mut writer =
            ArrowWriter::try_new(file, batch.schema(), None).expect("create parquet writer");
        writer.write(batch).expect("write parquet batch");
        writer.close().expect("close parquet writer");
        ctx.register_parquet(
            table_name,
            path.to_str().expect("utf8 parquet path"),
            ParquetReadOptions::default(),
        )
        .await
        .expect("register parquet table");
        dir
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
        for backend in ["parquet", "fjall"] {
            let ctx = create_vep_session();
            ctx.register_table("vcf_data", Arc::new(vcf_table()))
                .expect("register vcf table");
            ctx.register_table("var_cache", Arc::new(cache_table()))
                .expect("register cache table");

            let sql = format!("SELECT * FROM annotate_vep('vcf_data', 'var_cache', '{backend}')");
            let df = ctx
                .sql(&sql)
                .await
                .expect("annotate_vep query should parse");

            let batches = df.collect().await.expect("collect annotate_vep");
            let total_rows: usize = batches.iter().map(|b| b.num_rows()).sum();
            assert_eq!(total_rows, 2, "backend={backend}");

            let mut csq_values = Vec::new();
            let mut most_values = Vec::new();
            for batch in &batches {
                assert!(batch.column_by_name("csq").is_some());
                assert!(batch.column_by_name("most_severe_consequence").is_some());
                csq_values.extend(string_values(
                    batch.column_by_name("csq").expect("csq column exists"),
                ));
                most_values.extend(string_values(
                    batch
                        .column_by_name("most_severe_consequence")
                        .expect("most_severe_consequence column exists"),
                ));
            }
            assert!(
                csq_values
                    .iter()
                    .any(|v| v.as_ref().is_some_and(|s| s.contains("sequence_variant"))),
                "backend={backend}"
            );
            assert!(csq_values.iter().any(|v| v.is_none()), "backend={backend}");
            assert!(
                most_values
                    .iter()
                    .any(|v| v.as_ref() == Some(&"sequence_variant".to_string())),
                "backend={backend}"
            );
            assert!(most_values.iter().any(|v| v.is_none()), "backend={backend}");
        }
    }

    #[tokio::test(flavor = "multi_thread")]
    async fn test_annotate_vep_projection_includes_null_placeholder_fields() {
        for backend in ["parquet", "fjall"] {
            let ctx = create_vep_session();
            ctx.register_table("vcf_data", Arc::new(vcf_table()))
                .expect("register vcf table");
            ctx.register_table("var_cache", Arc::new(cache_table()))
                .expect("register cache table");

            let sql = format!(
                "SELECT chrom, csq FROM annotate_vep('vcf_data', 'var_cache', '{backend}')"
            );
            let df = ctx.sql(&sql).await.expect("projection query should parse");

            let batches = df.collect().await.expect("collect projected annotate_vep");
            let total_rows: usize = batches.iter().map(|b| b.num_rows()).sum();
            assert_eq!(total_rows, 2, "backend={backend}");
            for batch in &batches {
                assert_eq!(batch.num_columns(), 2, "backend={backend}");
                assert_eq!(batch.schema().field(0).name(), "chrom");
                assert_eq!(batch.schema().field(1).name(), "csq");
            }
            let mut csq_values = Vec::new();
            for batch in &batches {
                csq_values.extend(string_values(
                    batch.column_by_name("csq").expect("csq column exists"),
                ));
            }
            assert!(csq_values.iter().any(|v| v.is_some()), "backend={backend}");
            assert!(csq_values.iter().any(|v| v.is_none()), "backend={backend}");
        }
    }

    #[tokio::test(flavor = "multi_thread")]
    async fn test_annotate_vep_uses_transcript_context_tables_when_available() {
        for backend in ["parquet", "fjall"] {
            let ctx = create_vep_session();
            ctx.register_table("vcf_data", Arc::new(vcf_table()))
                .expect("register vcf table");
            ctx.register_table("var_cache", Arc::new(cache_table()))
                .expect("register cache table");
            ctx.register_table("var_cache_transcripts", Arc::new(transcripts_table()))
                .expect("register transcripts table");
            ctx.register_table("var_cache_exons", Arc::new(exons_table()))
                .expect("register exons table");

            let sql = format!(
                "SELECT chrom, csq, most_severe_consequence \
                 FROM annotate_vep('vcf_data', 'var_cache', '{backend}') \
                 ORDER BY chrom"
            );
            let df = ctx.sql(&sql).await.expect("query should parse");

            let batches = df
                .collect()
                .await
                .expect("collect transcript-aware annotate_vep");
            let total_rows: usize = batches.iter().map(|b| b.num_rows()).sum();
            assert_eq!(total_rows, 2, "backend={backend}");

            let mut chrom = Vec::new();
            let mut csq = Vec::new();
            let mut most = Vec::new();
            for batch in &batches {
                chrom.extend(string_values(
                    batch.column_by_name("chrom").expect("chrom column exists"),
                ));
                csq.extend(string_values(
                    batch.column_by_name("csq").expect("csq column exists"),
                ));
                most.extend(string_values(
                    batch
                        .column_by_name("most_severe_consequence")
                        .expect("most_severe_consequence column exists"),
                ));
            }

            assert_eq!(
                chrom,
                vec![Some("1".to_string()), Some("2".to_string())],
                "backend={backend}"
            );
            assert!(csq.iter().all(|v| v.is_some()), "backend={backend}");
            assert!(most.iter().all(|v| v.is_some()), "backend={backend}");
            assert!(
                csq[0]
                    .as_ref()
                    .is_some_and(|s| s.contains("coding_sequence_variant")),
                "backend={backend}"
            );
            assert!(
                csq[1]
                    .as_ref()
                    .is_some_and(|s| s.contains("non_coding_transcript_exon_variant")),
                "backend={backend}"
            );
            assert_eq!(
                most,
                vec![
                    Some("coding_sequence_variant".to_string()),
                    Some("non_coding_transcript_exon_variant".to_string())
                ],
                "backend={backend}"
            );
        }
    }

    #[tokio::test(flavor = "multi_thread")]
    async fn test_annotate_vep_uses_options_json_table_overrides_for_transcript_context() {
        for backend in ["parquet", "fjall"] {
            let ctx = create_vep_session();
            ctx.register_table("vcf_data", Arc::new(vcf_table()))
                .expect("register vcf table");
            ctx.register_table("var_cache", Arc::new(cache_table()))
                .expect("register cache table");
            ctx.register_table("tx_ctx", Arc::new(transcripts_table()))
                .expect("register transcript context table");
            ctx.register_table("ex_ctx", Arc::new(exons_table()))
                .expect("register exon context table");

            let sql = format!(
                "SELECT csq, most_severe_consequence \
                 FROM annotate_vep( \
                   'vcf_data', \
                   'var_cache', \
                   '{backend}', \
                   '{{\"transcripts_table\":\"tx_ctx\",\"exons_table\":\"ex_ctx\"}}' \
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
                    batch.column_by_name("csq").expect("csq column exists"),
                ));
                most.extend(string_values(
                    batch
                        .column_by_name("most_severe_consequence")
                        .expect("most_severe_consequence column exists"),
                ));
            }
            assert!(csq.iter().all(|v| v.is_some()), "backend={backend}");
            assert!(most.iter().all(|v| v.is_some()), "backend={backend}");
        }
    }

    // Mirrors Ensembl VEP Runner distance coverage:
    // - https://github.com/Ensembl/ensembl-vep/blob/release/115/t/Runner.t#L535-L571
    #[tokio::test(flavor = "multi_thread")]
    async fn test_annotate_vep_respects_options_json_distance_for_upstream_and_downstream() {
        for backend in ["parquet", "fjall"] {
            let ctx = create_vep_session();
            ctx.register_table("vcf_distance", Arc::new(distance_vcf_table()))
                .expect("register distance vcf");
            ctx.register_table("var_distance_cache", Arc::new(distance_cache_table()))
                .expect("register distance cache");
            ctx.register_table(
                "var_distance_cache_transcripts",
                Arc::new(distance_transcripts_table()),
            )
            .expect("register distance transcripts");
            ctx.register_table("var_distance_cache_exons", Arc::new(distance_exons_table()))
                .expect("register distance exons");

            let default_batches = ctx
                .sql(&format!(
                    "SELECT csq FROM annotate_vep('vcf_distance', 'var_distance_cache', '{backend}')"
                ))
                .await
                .expect("default distance query should parse")
                .collect()
                .await
                .expect("collect default distance annotate_vep");
            let default_csq = string_values(
                default_batches[0]
                    .column_by_name("csq")
                    .expect("default csq column exists"),
            );
            let default_csq0 = default_csq[0]
                .as_ref()
                .expect("default csq should be present");
            assert!(
                default_csq0.contains("intergenic_variant"),
                "backend={backend}"
            );
            assert!(!default_csq0.contains("ENSTDISTUP"), "backend={backend}");
            assert!(!default_csq0.contains("ENSTDISTDOWN"), "backend={backend}");

            let numeric_batches = ctx
                .sql(&format!(
                    "SELECT csq FROM annotate_vep( \
                       'vcf_distance', \
                       'var_distance_cache', \
                       '{backend}', \
                       '{{\"distance\":10000}}' \
                     )"
                ))
                .await
                .expect("numeric distance query should parse")
                .collect()
                .await
                .expect("collect numeric distance annotate_vep");
            let numeric_csq = string_values(
                numeric_batches[0]
                    .column_by_name("csq")
                    .expect("numeric csq column exists"),
            );
            let numeric_csq0 = numeric_csq[0]
                .as_ref()
                .expect("numeric csq should be present");
            let numeric_entries = csq_entries(numeric_csq0);
            assert_eq!(numeric_entries.len(), 1, "backend={backend}");
            let upstream_entry = find_csq_entry(numeric_csq0, "Transcript", "ENSTDISTUP");
            assert_eq!(
                upstream_entry[1], "upstream_gene_variant",
                "backend={backend}"
            );
            assert_eq!(upstream_entry[18], "7079", "backend={backend}");
            assert!(!numeric_csq0.contains("ENSTDISTDOWN"), "backend={backend}");

            let pair_batches = ctx
                .sql(&format!(
                    "SELECT csq FROM annotate_vep( \
                       'vcf_distance', \
                       'var_distance_cache', \
                       '{backend}', \
                       '{{\"distance\":\"10000,20000\"}}' \
                     )"
                ))
                .await
                .expect("pair distance query should parse")
                .collect()
                .await
                .expect("collect pair distance annotate_vep");
            let pair_csq = string_values(
                pair_batches[0]
                    .column_by_name("csq")
                    .expect("pair csq column exists"),
            );
            let pair_csq0 = pair_csq[0].as_ref().expect("pair csq should be present");
            let pair_entries = csq_entries(pair_csq0);
            assert_eq!(pair_entries.len(), 2, "backend={backend}");
            let upstream_entry = find_csq_entry(pair_csq0, "Transcript", "ENSTDISTUP");
            let downstream_entry = find_csq_entry(pair_csq0, "Transcript", "ENSTDISTDOWN");
            assert_eq!(
                upstream_entry[1], "upstream_gene_variant",
                "backend={backend}"
            );
            assert_eq!(upstream_entry[18], "7079", "backend={backend}");
            assert_eq!(
                downstream_entry[1], "downstream_gene_variant",
                "backend={backend}"
            );
            assert_eq!(downstream_entry[18], "15000", "backend={backend}");
        }
    }

    // Mirrors Ensembl VEP offline HGVS gating and output toggles:
    // - https://github.com/Ensembl/ensembl-vep/blob/release/115/modules/Bio/EnsEMBL/VEP/Runner.pm#L726-L738
    // - https://github.com/Ensembl/ensembl-vep/blob/release/115/modules/Bio/EnsEMBL/VEP/OutputFactory.pm#L1698-L1715
    #[tokio::test(flavor = "multi_thread")]
    async fn test_annotate_vep_hgvs_fields_require_flag_and_reference_fasta() {
        for backend in ["parquet", "fjall"] {
            let ctx = create_vep_session();
            ctx.register_table("vcf_hgvs", Arc::new(hgvs_vcf_table()))
                .expect("register hgvs vcf");
            ctx.register_table("var_hgvs_cache", Arc::new(hgvs_cache_table()))
                .expect("register hgvs cache");
            ctx.register_table(
                "var_hgvs_cache_transcripts",
                Arc::new(hgvs_transcripts_table()),
            )
            .expect("register hgvs transcripts");
            ctx.register_table("var_hgvs_cache_exons", Arc::new(hgvs_exons_table()))
                .expect("register hgvs exons");
            ctx.register_table(
                "var_hgvs_cache_translations",
                Arc::new(hgvs_translations_table()),
            )
            .expect("register hgvs translations");

            let default_batches = ctx
                .sql(&format!(
                    "SELECT csq FROM annotate_vep('vcf_hgvs', 'var_hgvs_cache', '{backend}')"
                ))
                .await
                .expect("default hgvs query should parse")
                .collect()
                .await
                .expect("collect default hgvs annotate_vep");
            let default_csq = string_values(
                default_batches[0]
                    .column_by_name("csq")
                    .expect("default csq column exists"),
            );
            let default_entry = csq_entries(default_csq[0].as_ref().expect("default csq present"))
                .into_iter()
                .find(|fields| fields.len() == 74 && fields[5] == "Transcript")
                .expect("transcript entry should exist");
            assert_eq!(default_entry[10], "", "backend={backend}");
            assert_eq!(default_entry[11], "", "backend={backend}");

            let missing_fasta_err = ctx
                .sql(&format!(
                    "SELECT csq FROM annotate_vep( \
                       'vcf_hgvs', \
                       'var_hgvs_cache', \
                       '{backend}', \
                       '{{\"hgvs\":true}}' \
                     )"
                ))
                .await
                .expect("hgvs query should parse")
                .collect()
                .await
                .expect_err("hgvs without reference fasta should fail")
                .to_string();
            assert!(
                missing_fasta_err.contains("Cannot generate HGVS coordinates"),
                "backend={backend}"
            );

            let (_fasta_dir, fasta_path) = write_test_indexed_fasta("1", "NNNNNNNNNN");

            let hgvs_batches = ctx
                .sql(&format!(
                    "SELECT csq FROM annotate_vep( \
                       'vcf_hgvs', \
                       'var_hgvs_cache', \
                       '{backend}', \
                       '{{\"hgvs\":true,\"reference_fasta_path\":\"{}\"}}' \
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
                    .column_by_name("csq")
                    .expect("hgvs csq column exists"),
            );
            let hgvs_entry = csq_entries(hgvs_csq[0].as_ref().expect("hgvs csq present"))
                .into_iter()
                .find(|fields| fields.len() == 74 && fields[5] == "Transcript")
                .expect("transcript hgvs entry should exist");
            assert_eq!(
                hgvs_entry[10], "ENSTHGVS000001.1:c.4G>A",
                "backend={backend}"
            );
            assert_eq!(
                hgvs_entry[11], "ENSPHGVS000001.1:p.Ala2Thr",
                "backend={backend}"
            );

            let hgvs_prediction_batches = ctx
                .sql(&format!(
                    "SELECT csq FROM annotate_vep( \
                       'vcf_hgvs', \
                       'var_hgvs_cache', \
                       '{backend}', \
                       '{{\"hgvsp\":true,\"hgvsp_use_prediction\":true,\"reference_fasta_path\":\"{}\"}}' \
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
                    .column_by_name("csq")
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
            assert_eq!(hgvs_prediction_entry[10], "", "backend={backend}");
            assert_eq!(
                hgvs_prediction_entry[11], "ENSPHGVS000001.1:p.(Ala2Thr)",
                "backend={backend}"
            );
        }
    }

    #[tokio::test(flavor = "multi_thread")]
    async fn test_annotate_vep_uses_translation_context_for_synonymous_classification() {
        for backend in ["parquet", "fjall"] {
            let ctx = create_vep_session();
            ctx.register_table("vcf_syn", Arc::new(synonymous_vcf_table()))
                .expect("register synonymous vcf table");
            ctx.register_table("var_syn_cache", Arc::new(synonymous_cache_table()))
                .expect("register synonymous cache table");
            ctx.register_table(
                "var_syn_cache_transcripts",
                Arc::new(synonymous_transcripts_table()),
            )
            .expect("register transcripts table");
            ctx.register_table("var_syn_cache_exons", Arc::new(synonymous_exons_table()))
                .expect("register exons table");
            ctx.register_table(
                "var_syn_cache_translations",
                Arc::new(synonymous_translations_table()),
            )
            .expect("register translations table");

            let sql = format!(
                "SELECT csq, most_severe_consequence \
                 FROM annotate_vep('vcf_syn', 'var_syn_cache', '{backend}')"
            );
            let df = ctx.sql(&sql).await.expect("query should parse");
            let batches = df.collect().await.expect("collect annotate_vep");
            let mut csq = Vec::new();
            let mut most = Vec::new();
            for batch in &batches {
                csq.extend(string_values(
                    batch.column_by_name("csq").expect("csq column exists"),
                ));
                most.extend(string_values(
                    batch
                        .column_by_name("most_severe_consequence")
                        .expect("most_severe_consequence column exists"),
                ));
            }
            assert_eq!(csq.len(), 1, "backend={backend}");
            assert_eq!(most.len(), 1, "backend={backend}");
            let csq0 = csq[0].as_ref().expect("csq should be present");
            assert!(csq0.contains("synonymous_variant"), "backend={backend}");
            assert!(!csq0.contains("missense_variant"), "backend={backend}");
            assert_eq!(
                most[0],
                Some("synonymous_variant".to_string()),
                "backend={backend}: translation-aware codon classification should produce synonymous_variant",
            );
        }
    }

    #[tokio::test(flavor = "multi_thread")]
    async fn test_annotate_vep_csq_terms_deterministic_and_backend_consistent_with_context_tables()
    {
        let ctx = create_vep_session();
        ctx.register_table("vcf_ctx", Arc::new(context_vcf_table()))
            .expect("register context vcf table");
        ctx.register_table("var_ctx_cache", Arc::new(context_cache_table()))
            .expect("register context cache table");
        ctx.register_table(
            "var_ctx_cache_transcripts",
            Arc::new(context_transcripts_table()),
        )
        .expect("register context transcript table");
        ctx.register_table("var_ctx_cache_exons", Arc::new(context_exons_table()))
            .expect("register context exon table");
        ctx.register_table(
            "var_ctx_cache_regulatory_features",
            Arc::new(context_regulatory_table()),
        )
        .expect("register context regulatory table");
        ctx.register_table(
            "var_ctx_cache_motif_features",
            Arc::new(context_motif_table()),
        )
        .expect("register context motif table");
        ctx.register_table(
            "var_ctx_cache_mirna_features",
            Arc::new(context_mirna_table()),
        )
        .expect("register context mirna table");
        ctx.register_table("var_ctx_cache_sv_features", Arc::new(context_sv_table()))
            .expect("register context sv table");

        let parquet_df = ctx
            .sql(
                "SELECT csq, most_severe_consequence \
                 FROM annotate_vep('vcf_ctx', 'var_ctx_cache', 'parquet')",
            )
            .await
            .expect("parquet query should parse");
        let parquet_batches = parquet_df
            .collect()
            .await
            .expect("collect parquet annotate_vep");
        let parquet_csq = string_values(
            parquet_batches[0]
                .column_by_name("csq")
                .expect("csq column exists"),
        );
        let parquet_most = string_values(
            parquet_batches[0]
                .column_by_name("most_severe_consequence")
                .expect("most_severe_consequence column exists"),
        );
        let parquet_csq0 = parquet_csq[0]
            .as_ref()
            .expect("parquet csq should be present")
            .to_string();
        let parquet_most0 = parquet_most[0]
            .as_ref()
            .expect("parquet most severe should be present")
            .to_string();

        let fjall_df = ctx
            .sql(
                "SELECT csq, most_severe_consequence \
                 FROM annotate_vep('vcf_ctx', 'var_ctx_cache', 'fjall')",
            )
            .await
            .expect("fjall query should parse");
        let fjall_batches = fjall_df
            .collect()
            .await
            .expect("collect fjall annotate_vep");
        let fjall_csq = string_values(
            fjall_batches[0]
                .column_by_name("csq")
                .expect("csq column exists"),
        );
        let fjall_most = string_values(
            fjall_batches[0]
                .column_by_name("most_severe_consequence")
                .expect("most_severe_consequence column exists"),
        );
        let fjall_csq0 = fjall_csq[0]
            .as_ref()
            .expect("fjall csq should be present")
            .to_string();
        let fjall_most0 = fjall_most[0]
            .as_ref()
            .expect("fjall most severe should be present")
            .to_string();

        assert_eq!(
            parquet_csq0, fjall_csq0,
            "CSQ should be deterministic across backend modes for identical logical context"
        );
        assert_eq!(
            parquet_most0, fjall_most0,
            "most_severe_consequence should be deterministic across backend modes"
        );
        assert_terms_sorted_by_rank(&parquet_csq0);
        assert!(parquet_csq0.contains("transcript_ablation"));
        assert!(parquet_csq0.contains("regulatory_region_variant"));
        assert!(parquet_csq0.contains("TF_binding_site_variant"));
        assert!(parquet_csq0.contains("mature_miRNA_variant"));
        assert_eq!(parquet_most0, "transcript_ablation");
    }

    // Mirrors the cache-backed regulatory source coverage from:
    // - https://github.com/Ensembl/ensembl-vep/blob/release/115/t/AnnotationSource_Cache_RegFeat.t#L166-L205
    #[tokio::test(flavor = "multi_thread")]
    async fn test_annotate_vep_deduplicates_duplicate_regulatory_rows_from_context_table() {
        for backend in ["parquet", "fjall"] {
            let ctx = create_vep_session();
            ctx.register_table("vcf_reg_dup", Arc::new(context_vcf_table()))
                .expect("register duplicate regulatory vcf");
            ctx.register_table("var_reg_dup_cache", Arc::new(context_cache_table()))
                .expect("register duplicate regulatory cache");

            let duplicate_batch = regulatory_feature_batch(&[
                ("ENSR_DUP", "1", 150, 160, Some("promoter")),
                ("ENSR_DUP", "1", 150, 160, Some("promoter")),
            ]);
            // For parquet backend, test via Parquet file; for fjall, use MemTable.
            // Both exercise the same dedup logic.
            let _parquet_dir = if backend == "parquet" {
                Some(register_single_batch_parquet(&ctx, "reg_dup_table", &duplicate_batch).await)
            } else {
                let dup_schema = duplicate_batch.schema();
                let dup_mem =
                    MemTable::try_new(dup_schema, vec![vec![duplicate_batch.clone()]]).unwrap();
                ctx.register_table("reg_dup_table", Arc::new(dup_mem))
                    .unwrap();
                None
            };

            let sql = format!(
                "SELECT csq, most_severe_consequence \
                 FROM annotate_vep( \
                   'vcf_reg_dup', \
                   'var_reg_dup_cache', \
                   '{backend}', \
                   '{{\"regulatory_table\":\"reg_dup_table\"}}' \
                 )"
            );
            let batches = ctx
                .sql(&sql)
                .await
                .expect("duplicate regulatory query should parse")
                .collect()
                .await
                .expect("collect duplicate regulatory annotate_vep");

            let csq = string_values(batches[0].column_by_name("csq").expect("csq column exists"));
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
            assert_eq!(regulatory_entries.len(), 1, "backend={backend}");
            assert_eq!(regulatory_entries[0][6], "ENSR_DUP", "backend={backend}");
            assert_eq!(regulatory_entries[0][7], "promoter", "backend={backend}");
            assert_eq!(
                most[0],
                Some("regulatory_region_variant".to_string()),
                "backend={backend}"
            );
        }
    }

    // Mirrors Ensembl VEP regulatory output serialization coverage:
    // - https://github.com/Ensembl/ensembl-vep/blob/release/115/t/OutputFactory.t#L1207-L1236
    // - https://github.com/Ensembl/ensembl-vep/blob/release/115/modules/Bio/EnsEMBL/VEP/OutputFactory.pm#L1802-L1829
    #[tokio::test(flavor = "multi_thread")]
    async fn test_annotate_vep_regulatory_csq_serializer_matches_vep_output_factory_shape() {
        for backend in ["parquet", "fjall"] {
            let ctx = create_vep_session();
            ctx.register_table("vcf_reg_ser", Arc::new(context_vcf_table()))
                .expect("register regulatory serializer vcf");
            ctx.register_table("var_reg_ser_cache", Arc::new(context_cache_table()))
                .expect("register regulatory serializer cache");
            ctx.register_table(
                "var_reg_ser_cache_regulatory_features",
                Arc::new(regulatory_table_with_rows(&[(
                    "ENSR00000140763",
                    "1",
                    150,
                    160,
                    Some("promoter"),
                )])),
            )
            .expect("register regulatory serializer table");

            let sql = format!(
                "SELECT csq, most_severe_consequence \
                 FROM annotate_vep('vcf_reg_ser', 'var_reg_ser_cache', '{backend}')"
            );
            let batches = ctx
                .sql(&sql)
                .await
                .expect("regulatory serializer query should parse")
                .collect()
                .await
                .expect("collect regulatory serializer annotate_vep");

            let csq = string_values(batches[0].column_by_name("csq").expect("csq column exists"));
            let most = string_values(
                batches[0]
                    .column_by_name("most_severe_consequence")
                    .expect("most severe column exists"),
            );
            let csq0 = csq[0].as_ref().expect("csq should be present");
            let entries = csq_entries(csq0);
            assert_eq!(entries.len(), 2, "backend={backend}");

            let entry = find_csq_entry(csq0, "RegulatoryFeature", "ENSR00000140763");
            assert_eq!(entry.len(), 74, "backend={backend}");
            assert_eq!(entry[0], "G", "backend={backend}");
            assert_eq!(entry[1], "regulatory_region_variant", "backend={backend}");
            assert_eq!(entry[2], "MODIFIER", "backend={backend}");
            assert_eq!(entry[3], "", "backend={backend}");
            assert_eq!(entry[4], "", "backend={backend}");
            assert_eq!(entry[5], "RegulatoryFeature", "backend={backend}");
            assert_eq!(entry[6], "ENSR00000140763", "backend={backend}");
            assert_eq!(entry[7], "promoter", "backend={backend}");
            assert_eq!(entry[8], "", "backend={backend}");
            assert_eq!(entry[9], "", "backend={backend}");
            assert_eq!(entry[18], "", "backend={backend}");
            assert!(
                entries
                    .iter()
                    .any(|fields| fields.len() == 74 && fields[1] == "intergenic_variant"),
                "backend={backend}: regulatory-only output should retain the orthogonal intergenic entry"
            );
            assert_eq!(
                most[0],
                Some("regulatory_region_variant".to_string()),
                "backend={backend}"
            );
        }
    }

    #[tokio::test(flavor = "multi_thread")]
    async fn test_annotate_vep_golden_fixture_term_parity_for_context_reg_motif_mirna_sv() {
        let ctx = create_vep_session();
        ctx.register_table("vcf_golden_ctx", Arc::new(golden_context_vcf_table()))
            .expect("register golden context vcf table");
        ctx.register_table(
            "var_golden_ctx_cache",
            Arc::new(golden_context_cache_table()),
        )
        .expect("register golden context cache table");
        ctx.register_table(
            "var_golden_ctx_cache_transcripts",
            Arc::new(golden_context_transcripts_table()),
        )
        .expect("register golden context transcripts");
        ctx.register_table(
            "var_golden_ctx_cache_exons",
            Arc::new(golden_context_exons_table()),
        )
        .expect("register golden context exons");
        ctx.register_table(
            "var_golden_ctx_cache_regulatory_features",
            Arc::new(golden_context_regulatory_table()),
        )
        .expect("register golden context regulatory");
        ctx.register_table(
            "var_golden_ctx_cache_motif_features",
            Arc::new(golden_context_motif_table()),
        )
        .expect("register golden context motif");
        ctx.register_table(
            "var_golden_ctx_cache_mirna_features",
            Arc::new(golden_context_mirna_table()),
        )
        .expect("register golden context mirna");
        ctx.register_table(
            "var_golden_ctx_cache_sv_features",
            Arc::new(golden_context_sv_table()),
        )
        .expect("register golden context sv");

        let expected_terms = [
            "transcript_ablation",
            "feature_elongation",
            "feature_truncation",
            "transcript_amplification",
            "TFBS_ablation",
            "TFBS_amplification",
            "regulatory_region_ablation",
            "regulatory_region_amplification",
            "mature_miRNA_variant",
            "TF_binding_site_variant",
            "regulatory_region_variant",
            "intergenic_variant",
        ];

        let mut observed = Vec::new();
        for backend in ["parquet", "fjall"] {
            let sql = format!(
                "SELECT csq, most_severe_consequence \
                 FROM annotate_vep('vcf_golden_ctx', 'var_golden_ctx_cache', '{backend}')"
            );
            let batches = ctx
                .sql(&sql)
                .await
                .expect("query should parse")
                .collect()
                .await
                .expect("collect annotate_vep");
            let csq = string_values(batches[0].column_by_name("csq").expect("csq column exists"));
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
            assert_eq!(most0, "transcript_ablation");
            observed.push((csq0, most0));
        }

        assert_eq!(
            observed[0], observed[1],
            "golden context term parity should be backend-consistent"
        );
    }

    #[tokio::test(flavor = "multi_thread")]
    async fn test_annotate_vep_golden_fixture_term_parity_for_splice_boundary_case() {
        let ctx = create_vep_session();
        ctx.register_table("vcf_splice", Arc::new(splice_vcf_table()))
            .expect("register splice vcf");
        ctx.register_table("var_splice_cache", Arc::new(splice_cache_table()))
            .expect("register splice cache");
        ctx.register_table(
            "var_splice_cache_transcripts",
            Arc::new(splice_transcripts_table()),
        )
        .expect("register splice transcripts");
        ctx.register_table("var_splice_cache_exons", Arc::new(splice_exons_table()))
            .expect("register splice exons");

        let expected_terms = ["splice_donor_variant"];

        for backend in ["parquet", "fjall"] {
            let sql = format!(
                "SELECT csq, most_severe_consequence \
                 FROM annotate_vep('vcf_splice', 'var_splice_cache', '{backend}')"
            );
            let batches = ctx
                .sql(&sql)
                .await
                .expect("query should parse")
                .collect()
                .await
                .expect("collect annotate_vep");
            let csq = string_values(batches[0].column_by_name("csq").expect("csq column exists"));
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
        ctx.register_table("var_repeat_cache", Arc::new(repeat_shift_cache_table()))
            .expect("register repeat cache");
        ctx.register_table(
            "var_repeat_cache_transcripts",
            Arc::new(repeat_shift_transcripts_table()),
        )
        .expect("register repeat transcripts");
        ctx.register_table(
            "var_repeat_cache_exons",
            Arc::new(repeat_shift_exons_table()),
        )
        .expect("register repeat exons");

        let expected_terms = ["frameshift_variant"];

        for backend in ["parquet", "fjall"] {
            let sql = format!(
                "SELECT csq, most_severe_consequence \
                 FROM annotate_vep('vcf_repeat', 'var_repeat_cache', '{backend}')"
            );
            let batches = ctx
                .sql(&sql)
                .await
                .expect("query should parse")
                .collect()
                .await
                .expect("collect annotate_vep");
            let csq = string_values(batches[0].column_by_name("csq").expect("csq column exists"));
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
        for backend in ["parquet", "fjall"] {
            let ctx = create_vep_session();
            ctx.register_table("vcf_neg", Arc::new(neg_strand_vcf_table()))
                .expect("register negative-strand vcf");
            ctx.register_table("var_neg_cache", Arc::new(neg_strand_cache_table()))
                .expect("register negative-strand cache");
            ctx.register_table(
                "var_neg_cache_transcripts",
                Arc::new(neg_strand_transcripts_table()),
            )
            .expect("register negative-strand transcripts");
            ctx.register_table("var_neg_cache_exons", Arc::new(neg_strand_exons_table()))
                .expect("register negative-strand exons");
            ctx.register_table(
                "var_neg_cache_translations",
                Arc::new(neg_strand_translations_table()),
            )
            .expect("register negative-strand translations");

            let sql = format!(
                "SELECT csq, most_severe_consequence \
                 FROM annotate_vep('vcf_neg', 'var_neg_cache', '{backend}')"
            );
            let df = ctx.sql(&sql).await.expect("query should parse");
            let batches = df.collect().await.expect("collect annotate_vep");
            let csq = string_values(batches[0].column_by_name("csq").expect("csq column exists"));
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
                "backend={backend}"
            );
        }
    }

    #[tokio::test(flavor = "multi_thread")]
    async fn test_annotate_vep_translation_parity_handles_inframe_stop_loss_case() {
        for backend in ["parquet", "fjall"] {
            let ctx = create_vep_session();
            ctx.register_table("vcf_stop_loss", Arc::new(stop_loss_vcf_table()))
                .expect("register stop-loss vcf");
            ctx.register_table("var_stop_loss_cache", Arc::new(stop_loss_cache_table()))
                .expect("register stop-loss cache");
            ctx.register_table(
                "var_stop_loss_cache_transcripts",
                Arc::new(stop_loss_transcripts_table()),
            )
            .expect("register stop-loss transcripts");
            ctx.register_table(
                "var_stop_loss_cache_exons",
                Arc::new(stop_loss_exons_table()),
            )
            .expect("register stop-loss exons");
            ctx.register_table(
                "var_stop_loss_cache_translations",
                Arc::new(stop_loss_translations_table()),
            )
            .expect("register stop-loss translations");

            let sql = format!(
                "SELECT csq, most_severe_consequence \
                 FROM annotate_vep('vcf_stop_loss', 'var_stop_loss_cache', '{backend}')"
            );
            let df = ctx.sql(&sql).await.expect("query should parse");
            let batches = df.collect().await.expect("collect annotate_vep");
            let csq = string_values(batches[0].column_by_name("csq").expect("csq column exists"));
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
                "backend={backend}"
            );
        }
    }

    #[tokio::test(flavor = "multi_thread")]
    async fn test_annotate_vep_csq_has_74_pipe_delimited_fields_per_transcript() {
        for backend in ["parquet", "fjall"] {
            let ctx = create_vep_session();
            ctx.register_table("vcf_data", Arc::new(vcf_table()))
                .expect("register vcf table");
            ctx.register_table("var_cache", Arc::new(cache_table()))
                .expect("register cache table");
            ctx.register_table("var_cache_transcripts", Arc::new(transcripts_table()))
                .expect("register transcripts table");
            ctx.register_table("var_cache_exons", Arc::new(exons_table()))
                .expect("register exons table");

            let sql = format!(
                "SELECT csq FROM annotate_vep('vcf_data', 'var_cache', '{backend}') \
                 ORDER BY chrom"
            );
            let df = ctx.sql(&sql).await.expect("query should parse");
            let batches = df.collect().await.expect("collect annotate_vep");
            let csq_values = string_values(batches[0].column_by_name("csq").expect("csq column"));

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
}
