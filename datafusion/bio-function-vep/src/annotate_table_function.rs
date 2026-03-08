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
    use std::collections::BTreeSet;
    use std::sync::Arc;

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
        let schema = Arc::new(Schema::new(vec![
            Field::new("chrom", DataType::Utf8, false),
            Field::new("start", DataType::Int64, false),
            Field::new("end", DataType::Int64, false),
            Field::new("variation_name", DataType::Utf8, true),
            Field::new("allele_string", DataType::Utf8, false),
            Field::new("clin_sig", DataType::Utf8, true),
            Field::new("AF", DataType::Float64, true),
        ]));
        let batch = RecordBatch::try_new(
            schema.clone(),
            vec![
                Arc::new(StringArray::from(vec!["1"])),
                Arc::new(Int64Array::from(vec![100])),
                Arc::new(Int64Array::from(vec![101])),
                Arc::new(StringArray::from(vec!["rs100"])),
                Arc::new(StringArray::from(vec!["A/G"])),
                Arc::new(StringArray::from(vec!["benign"])),
                Arc::new(Float64Array::from(vec![0.12_f64])),
            ],
        )
        .expect("valid cache batch");
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
                Arc::new(StringArray::from(vec!["ENST00000000010", "ENST00000000011"])),
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
                Arc::new(StringArray::from(vec!["ENST00000000010", "ENST00000000011"])),
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
                Arc::new(StringArray::from(vec!["ENST00000000004", "ENST00000000004"])),
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
        let ctx = create_vep_session();
        ctx.register_table("vcf_data", Arc::new(vcf_table()))
            .expect("register vcf table");
        ctx.register_table("var_cache", Arc::new(cache_table()))
            .expect("register cache table");

        let df = ctx
            .sql("SELECT * FROM annotate_vep('vcf_data', 'var_cache', 'parquet')")
            .await
            .expect("annotate_vep query should parse");

        let batches = df.collect().await.expect("collect annotate_vep");
        let total_rows: usize = batches.iter().map(|b| b.num_rows()).sum();
        assert_eq!(total_rows, 2);

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
                .any(|v| v.as_ref().is_some_and(|s| s.contains("sequence_variant")))
        );
        assert!(csq_values.iter().any(|v| v.is_none()));
        assert!(
            most_values
                .iter()
                .any(|v| v.as_ref() == Some(&"sequence_variant".to_string()))
        );
        assert!(most_values.iter().any(|v| v.is_none()));
    }

    #[tokio::test(flavor = "multi_thread")]
    async fn test_annotate_vep_projection_includes_null_placeholder_fields() {
        let ctx = create_vep_session();
        ctx.register_table("vcf_data", Arc::new(vcf_table()))
            .expect("register vcf table");
        ctx.register_table("var_cache", Arc::new(cache_table()))
            .expect("register cache table");

        let df = ctx
            .sql("SELECT chrom, csq FROM annotate_vep('vcf_data', 'var_cache', 'parquet')")
            .await
            .expect("projection query should parse");

        let batches = df.collect().await.expect("collect projected annotate_vep");
        let total_rows: usize = batches.iter().map(|b| b.num_rows()).sum();
        assert_eq!(total_rows, 2);
        for batch in &batches {
            assert_eq!(batch.num_columns(), 2);
            assert_eq!(batch.schema().field(0).name(), "chrom");
            assert_eq!(batch.schema().field(1).name(), "csq");
        }
        let mut csq_values = Vec::new();
        for batch in &batches {
            csq_values.extend(string_values(
                batch.column_by_name("csq").expect("csq column exists"),
            ));
        }
        assert!(csq_values.iter().any(|v| v.is_some()));
        assert!(csq_values.iter().any(|v| v.is_none()));
    }

    #[tokio::test(flavor = "multi_thread")]
    async fn test_annotate_vep_uses_transcript_context_tables_when_available() {
        let ctx = create_vep_session();
        ctx.register_table("vcf_data", Arc::new(vcf_table()))
            .expect("register vcf table");
        ctx.register_table("var_cache", Arc::new(cache_table()))
            .expect("register cache table");
        ctx.register_table("var_cache_transcripts", Arc::new(transcripts_table()))
            .expect("register transcripts table");
        ctx.register_table("var_cache_exons", Arc::new(exons_table()))
            .expect("register exons table");

        let df = ctx
            .sql(
                "SELECT chrom, csq, most_severe_consequence \
                 FROM annotate_vep('vcf_data', 'var_cache', 'parquet') \
                 ORDER BY chrom",
            )
            .await
            .expect("query should parse");

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
                batch.column_by_name("csq").expect("csq column exists"),
            ));
            most.extend(string_values(
                batch
                    .column_by_name("most_severe_consequence")
                    .expect("most_severe_consequence column exists"),
            ));
        }

        assert_eq!(chrom, vec![Some("1".to_string()), Some("2".to_string())]);
        assert!(csq.iter().all(|v| v.is_some()));
        assert!(most.iter().all(|v| v.is_some()));
        // Without translation tables, SNV in CDS produces
        // coding_sequence_variant (no codon evidence for missense).
        assert!(
            csq[0]
                .as_ref()
                .is_some_and(|s| s.contains("coding_sequence_variant"))
        );
        assert!(csq[0].as_ref().is_some_and(|s| s.contains("rs100")));
        assert!(
            csq[1]
                .as_ref()
                .is_some_and(|s| s.contains("non_coding_transcript_exon_variant"))
        );
        assert_eq!(
            most,
            vec![
                Some("coding_sequence_variant".to_string()),
                Some("non_coding_transcript_exon_variant".to_string())
            ]
        );
    }

    #[tokio::test(flavor = "multi_thread")]
    async fn test_annotate_vep_uses_options_json_table_overrides_for_transcript_context() {
        let ctx = create_vep_session();
        ctx.register_table("vcf_data", Arc::new(vcf_table()))
            .expect("register vcf table");
        ctx.register_table("var_cache", Arc::new(cache_table()))
            .expect("register cache table");
        ctx.register_table("tx_ctx", Arc::new(transcripts_table()))
            .expect("register transcript context table");
        ctx.register_table("ex_ctx", Arc::new(exons_table()))
            .expect("register exon context table");

        let df = ctx
            .sql(
                "SELECT csq, most_severe_consequence \
                 FROM annotate_vep( \
                   'vcf_data', \
                   'var_cache', \
                   'parquet', \
                   '{\"transcripts_table\":\"tx_ctx\",\"exons_table\":\"ex_ctx\"}' \
                 )",
            )
            .await
            .expect("query should parse");

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
        assert!(csq.iter().all(|v| v.is_some()));
        assert!(most.iter().all(|v| v.is_some()));
    }

    #[tokio::test(flavor = "multi_thread")]
    async fn test_annotate_vep_uses_translation_context_for_synonymous_classification() {
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

        let df = ctx
            .sql(
                "SELECT csq, most_severe_consequence \
                 FROM annotate_vep('vcf_syn', 'var_syn_cache', 'parquet')",
            )
            .await
            .expect("query should parse");
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

        let expected_terms = [
            "splice_donor_variant",
            "intron_variant",
        ];

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

        let expected_terms = [
            "frameshift_variant",
        ];

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
            assert!(
                csq0.contains("|rs_repeat_shift|"),
                "repeat-shifted allele should resolve to cache variation_name via extended probes, got: {csq0}"
            );
            assert_terms_sorted_by_rank(csq0);
            assert_term_set_exact(csq0, &expected_terms);
            assert_eq!(most[0], Some("frameshift_variant".to_string()));
        }
    }

    #[tokio::test(flavor = "multi_thread")]
    async fn test_annotate_vep_translation_parity_handles_negative_strand_synonymous_case() {
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

        let df = ctx
            .sql(
                "SELECT csq, most_severe_consequence \
                 FROM annotate_vep('vcf_neg', 'var_neg_cache', 'parquet')",
            )
            .await
            .expect("query should parse");
        let batches = df.collect().await.expect("collect annotate_vep");
        let csq = string_values(batches[0].column_by_name("csq").expect("csq column exists"));
        let most = string_values(
            batches[0]
                .column_by_name("most_severe_consequence")
                .expect("most_severe_consequence column exists"),
        );
        let csq0 = csq[0].as_ref().expect("csq should be present");
        assert_terms_sorted_by_rank(csq0);
        assert_term_set_exact(
            csq0,
            &[
                "synonymous_variant",
            ],
        );
        assert_eq!(most[0], Some("synonymous_variant".to_string()));
    }

    #[tokio::test(flavor = "multi_thread")]
    async fn test_annotate_vep_translation_parity_handles_inframe_stop_loss_case() {
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

        let df = ctx
            .sql(
                "SELECT csq, most_severe_consequence \
                 FROM annotate_vep('vcf_stop_loss', 'var_stop_loss_cache', 'parquet')",
            )
            .await
            .expect("query should parse");
        let batches = df.collect().await.expect("collect annotate_vep");
        let csq = string_values(batches[0].column_by_name("csq").expect("csq column exists"));
        let most = string_values(
            batches[0]
                .column_by_name("most_severe_consequence")
                .expect("most_severe_consequence column exists"),
        );
        let csq0 = csq[0].as_ref().expect("csq should be present");
        assert_terms_sorted_by_rank(csq0);
        assert_term_set_exact(
            csq0,
            &[
                "stop_lost",
                "inframe_deletion",
            ],
        );
        assert_eq!(most[0], Some("stop_lost".to_string()));
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
