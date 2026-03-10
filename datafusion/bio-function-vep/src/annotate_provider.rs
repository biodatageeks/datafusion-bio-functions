//! Provider for `annotate_vep()` table function.
//!
//! Runtime behavior:
//! - always starts from `lookup_variants()` for known-variant metadata,
//! - when transcript/exon tables are available, computes transcript-driven
//!   consequence terms and most-severe ranking,
//! - otherwise falls back to phase-1.5 known-variant CSQ placeholders.

use std::any::Any;
use std::collections::HashSet;
use std::collections::hash_map::DefaultHasher;
use std::fmt::{Debug, Formatter};
use std::hash::{Hash, Hasher};
use std::sync::Arc;

use async_trait::async_trait;
use datafusion::arrow::array::{
    Array, AsArray, BooleanArray, Float32Array, Float64Array, Int8Array, Int16Array, Int32Array,
    Int64Array, LargeStringArray, ListArray, RecordBatch, StringArray, StringBuilder,
    StringViewArray, UInt8Array, UInt16Array, UInt32Array, UInt64Array, new_null_array,
};
use datafusion::arrow::datatypes::{DataType, Field, Schema, SchemaRef};
use datafusion::catalog::Session;
use datafusion::common::{DataFusionError, Result};
use datafusion::datasource::{MemTable, TableProvider, TableType};
use datafusion::physical_plan::ExecutionPlan;
use datafusion::prelude::{Expr, ParquetReadOptions, SessionContext};

use crate::allele::vcf_to_vep_allele;
use crate::annotation_store::{AnnotationBackend, build_store};
#[cfg(feature = "kv-cache")]
use crate::kv_cache::KvCacheTableProvider;
use crate::so_terms::{SoImpact, SoTerm, most_severe_term};
use crate::transcript_consequence::{
    ExonFeature, MirnaFeature, MotifFeature, PreparedContext, RegulatoryFeature, StructuralFeature,
    SvEventKind, SvFeatureKind, TranscriptConsequenceEngine, TranscriptFeature, TranslationFeature,
    VariantInput, is_vep_transcript,
};

/// Known variation cache annotation columns exposed as top-level output fields.
/// All are nullable Utf8. Columns not present in the actual cache emit NULLs.
pub const CACHE_OUTPUT_COLUMNS: &[&str] = &[
    // Variant identity
    "variation_name",
    // Clinical
    "clin_sig",
    "clin_sig_allele",
    "clinical_impact",
    "phenotype_or_disease",
    "pubmed",
    // Flags
    "somatic",
    "minor_allele",
    "minor_allele_freq",
    // 1000 Genomes
    "AF",
    "AFR",
    "AMR",
    "EAS",
    "EUR",
    "SAS",
    // gnomAD exome
    "gnomADe",
    "gnomADe_AFR",
    "gnomADe_AMR",
    "gnomADe_ASJ",
    "gnomADe_EAS",
    "gnomADe_FIN",
    "gnomADe_NFE",
    "gnomADe_SAS",
    "gnomADe_MID",
    "gnomADe_REMAINING",
    // gnomAD genome
    "gnomADg",
    "gnomADg_AFR",
    "gnomADg_AMI",
    "gnomADg_AMR",
    "gnomADg_ASJ",
    "gnomADg_EAS",
    "gnomADg_FIN",
    "gnomADg_MID",
    "gnomADg_NFE",
    "gnomADg_SAS",
    "gnomADg_REMAINING",
    // Cross-reference IDs
    "clinvar_ids",
    "cosmic_ids",
    "dbsnp_ids",
];

/// AF column definition: how to read, emit, and name each frequency population.
struct AfColumn {
    /// Column name in the variation cache parquet (e.g. `"gnomADg_FIN"`).
    cache_col: &'static str,
    /// Flag group: 0 = `--af`, 1 = `--af_1kg`, 2 = `--af_gnomade`, 3 = `--af_gnomadg`.
    flag_group: u8,
    /// Whether VEP emits this field's frequency in the individual CSQ slot.
    /// VEP's offline cache mode only emits global AF + 1000G sub-pops + gnomAD global;
    /// gnomAD sub-population frequencies are NOT emitted in individual CSQ fields.
    emit_in_csq: bool,
    /// Population name for MAX_AF_POPS (VEP-internal naming convention).
    /// `None` means this entry is excluded from MAX_AF computation (globals).
    max_af_pop: Option<&'static str>,
}

const AF_COLUMNS: &[AfColumn] = &[
    // --af (global 1000 Genomes) — emitted in CSQ, excluded from MAX_AF_POPS
    AfColumn { cache_col: "AF", flag_group: 0, emit_in_csq: true, max_af_pop: None },
    // --af_1kg (continental) — emitted, MAX_AF uses short names (AFR not AFR_AF)
    AfColumn { cache_col: "AFR", flag_group: 1, emit_in_csq: true, max_af_pop: Some("AFR") },
    AfColumn { cache_col: "AMR", flag_group: 1, emit_in_csq: true, max_af_pop: Some("AMR") },
    AfColumn { cache_col: "EAS", flag_group: 1, emit_in_csq: true, max_af_pop: Some("EAS") },
    AfColumn { cache_col: "EUR", flag_group: 1, emit_in_csq: true, max_af_pop: Some("EUR") },
    AfColumn { cache_col: "SAS", flag_group: 1, emit_in_csq: true, max_af_pop: Some("SAS") },
    // --af_gnomade — only global emitted in CSQ; sub-pops used for MAX_AF only
    AfColumn { cache_col: "gnomADe", flag_group: 2, emit_in_csq: true, max_af_pop: None },
    AfColumn { cache_col: "gnomADe_AFR", flag_group: 2, emit_in_csq: false, max_af_pop: Some("gnomADe_AFR") },
    AfColumn { cache_col: "gnomADe_AMR", flag_group: 2, emit_in_csq: false, max_af_pop: Some("gnomADe_AMR") },
    AfColumn { cache_col: "gnomADe_ASJ", flag_group: 2, emit_in_csq: false, max_af_pop: Some("gnomADe_ASJ") },
    AfColumn { cache_col: "gnomADe_EAS", flag_group: 2, emit_in_csq: false, max_af_pop: Some("gnomADe_EAS") },
    AfColumn { cache_col: "gnomADe_FIN", flag_group: 2, emit_in_csq: false, max_af_pop: Some("gnomADe_FIN") },
    AfColumn { cache_col: "gnomADe_MID", flag_group: 2, emit_in_csq: false, max_af_pop: Some("gnomADe_MID") },
    AfColumn { cache_col: "gnomADe_NFE", flag_group: 2, emit_in_csq: false, max_af_pop: Some("gnomADe_NFE") },
    AfColumn { cache_col: "gnomADe_REMAINING", flag_group: 2, emit_in_csq: false, max_af_pop: Some("gnomADe_REMAINING") },
    AfColumn { cache_col: "gnomADe_SAS", flag_group: 2, emit_in_csq: false, max_af_pop: Some("gnomADe_SAS") },
    // --af_gnomadg — only global emitted in CSQ; sub-pops used for MAX_AF only
    AfColumn { cache_col: "gnomADg", flag_group: 3, emit_in_csq: true, max_af_pop: None },
    AfColumn { cache_col: "gnomADg_AFR", flag_group: 3, emit_in_csq: false, max_af_pop: Some("gnomADg_AFR") },
    AfColumn { cache_col: "gnomADg_AMI", flag_group: 3, emit_in_csq: false, max_af_pop: Some("gnomADg_AMI") },
    AfColumn { cache_col: "gnomADg_AMR", flag_group: 3, emit_in_csq: false, max_af_pop: Some("gnomADg_AMR") },
    AfColumn { cache_col: "gnomADg_ASJ", flag_group: 3, emit_in_csq: false, max_af_pop: Some("gnomADg_ASJ") },
    AfColumn { cache_col: "gnomADg_EAS", flag_group: 3, emit_in_csq: false, max_af_pop: Some("gnomADg_EAS") },
    AfColumn { cache_col: "gnomADg_FIN", flag_group: 3, emit_in_csq: false, max_af_pop: Some("gnomADg_FIN") },
    AfColumn { cache_col: "gnomADg_MID", flag_group: 3, emit_in_csq: false, max_af_pop: Some("gnomADg_MID") },
    AfColumn { cache_col: "gnomADg_NFE", flag_group: 3, emit_in_csq: false, max_af_pop: Some("gnomADg_NFE") },
    AfColumn { cache_col: "gnomADg_REMAINING", flag_group: 3, emit_in_csq: false, max_af_pop: Some("gnomADg_REMAINING") },
    AfColumn { cache_col: "gnomADg_SAS", flag_group: 3, emit_in_csq: false, max_af_pop: Some("gnomADg_SAS") },
];

/// Parsed VEP option flags controlling which Batch 3 CSQ fields are emitted.
/// Flag names match Ensembl VEP CLI: `--check_existing`, `--af`, `--af_1kg`,
/// `--af_gnomade`, `--af_gnomadg`, `--max_af`, `--pubmed`.
#[derive(Debug, Clone)]
struct VepFlags {
    check_existing: bool,
    af: bool,
    af_1kg: bool,
    af_gnomade: bool,
    af_gnomadg: bool,
    max_af: bool,
    pubmed: bool,
}

impl VepFlags {
    fn from_options_json(options_json: Option<&str>) -> Self {
        let parse = |key| {
            options_json
                .and_then(|opts| AnnotateProvider::parse_json_bool_option(opts, key))
                .unwrap_or(false)
        };
        let af = parse("af");
        let af_1kg = parse("af_1kg");
        let af_gnomade = parse("af_gnomade");
        let af_gnomadg = parse("af_gnomadg");
        let max_af = parse("max_af");
        let pubmed = parse("pubmed");
        // VEP behavior: AF flags imply --check_existing.
        let check_existing =
            parse("check_existing") || af || af_1kg || af_gnomade || af_gnomadg || max_af || pubmed;
        Self {
            check_existing,
            af,
            af_1kg,
            af_gnomade,
            af_gnomadg,
            max_af,
            pubmed,
        }
    }

    /// Whether this AF column's flag group is enabled.
    fn af_group_enabled(&self, group: u8) -> bool {
        match group {
            0 => self.af,
            1 => self.af_1kg,
            2 => self.af_gnomade,
            3 => self.af_gnomadg,
            _ => false,
        }
    }
}

/// Parse a VEP cache `"allele:freq"` string and extract the frequency for the
/// specified VEP-minimized allele.
///
/// Cache format examples:
///   `"T:0.9301"`           — single allele
///   `"A:0.006,G:0.994"`    — multi-allele (comma-separated)
///   `"-:0.001"`            — deletion allele
fn extract_af_for_allele<'a>(cache_af_str: &'a str, vep_allele: &str) -> &'a str {
    if cache_af_str.is_empty() {
        return "";
    }
    for entry in cache_af_str.split(',') {
        if let Some((allele, freq)) = entry.split_once(':') {
            if allele == vep_allele {
                return freq;
            }
        }
    }
    ""
}

/// Compute MAX_AF and MAX_AF_POPS from collected (population_name, frequency_str) pairs.
///
/// Returns `(max_af_str, max_af_pops_str)` where `max_af_pops_str` uses `&` as
/// separator when multiple populations tie for the maximum (matching VEP format).
fn compute_max_af(af_entries: &[(&str, &str)]) -> (String, String) {
    let mut max_val: f64 = f64::NEG_INFINITY;
    let mut max_pops: Vec<&str> = Vec::new();
    let mut found_any = false;

    for &(pop_name, freq_str) in af_entries {
        if freq_str.is_empty() {
            continue;
        }
        let Ok(freq) = freq_str.parse::<f64>() else {
            continue;
        };
        found_any = true;
        if freq > max_val {
            max_val = freq;
            max_pops.clear();
            max_pops.push(pop_name);
        } else if (freq - max_val).abs() < f64::EPSILON {
            max_pops.push(pop_name);
        }
    }

    if !found_any {
        return (String::new(), String::new());
    }
    // Format MAX_AF the same way VEP does: fixed decimal, no trailing zeros.
    let max_af_str = format!("{max_val}");
    let max_af_pops_str = max_pops.join("&");
    (max_af_str, max_af_pops_str)
}

/// Table provider implementing `annotate_vep(...)`.
pub struct AnnotateProvider {
    session: Arc<SessionContext>,
    vcf_table: String,
    cache_source: String,
    backend: AnnotationBackend,
    options_json: Option<String>,
    schema: SchemaRef,
}

impl AnnotateProvider {
    pub fn new(
        session: Arc<SessionContext>,
        vcf_table: String,
        cache_source: String,
        backend: AnnotationBackend,
        options_json: Option<String>,
        vcf_schema: Schema,
    ) -> Self {
        // Output schema starts with all VCF columns and appends annotation fields.
        let mut fields: Vec<Arc<Field>> = vcf_schema
            .fields()
            .iter()
            .map(|field| {
                Arc::new(Field::new(
                    field.name(),
                    field.data_type().clone(),
                    field.is_nullable(),
                ))
            })
            .collect();

        fields.push(Arc::new(Field::new("csq", DataType::Utf8, true)));
        fields.push(Arc::new(Field::new(
            "most_severe_consequence",
            DataType::Utf8,
            true,
        )));
        for &col_name in CACHE_OUTPUT_COLUMNS {
            fields.push(Arc::new(Field::new(col_name, DataType::Utf8, true)));
        }

        Self {
            session,
            vcf_table,
            cache_source,
            backend,
            options_json,
            schema: Arc::new(Schema::new(fields)),
        }
    }

    async fn resolve_cache_table_name(&self) -> Result<String> {
        if self.session.table(&self.cache_source).await.is_ok() {
            return Ok(self.cache_source.clone());
        }

        let table_name = self.generated_cache_table_name();
        if self.session.table(&table_name).await.is_ok() {
            return Ok(table_name);
        }

        match self.backend {
            AnnotationBackend::Parquet => {
                self.session
                    .register_parquet(
                        &table_name,
                        &self.cache_source,
                        ParquetReadOptions::default(),
                    )
                    .await?;
            }
            AnnotationBackend::Fjall => {
                #[cfg(feature = "kv-cache")]
                {
                    let provider = KvCacheTableProvider::open(&self.cache_source).map_err(|e| {
                        DataFusionError::Execution(format!(
                            "annotate_vep(): failed to open fjall cache '{}': {e}",
                            self.cache_source
                        ))
                    })?;
                    self.session
                        .register_table(&table_name, Arc::new(provider))?;
                }
                #[cfg(not(feature = "kv-cache"))]
                {
                    return Err(DataFusionError::Execution(
                        "annotate_vep(): fjall backend requires kv-cache feature".to_string(),
                    ));
                }
            }
        }

        Ok(table_name)
    }

    fn generated_cache_table_name(&self) -> String {
        let mut hasher = DefaultHasher::new();
        self.backend.as_str().hash(&mut hasher);
        self.cache_source.hash(&mut hasher);
        format!(
            "__annotate_cache_{}_{:x}",
            self.backend.as_str(),
            hasher.finish()
        )
    }

    fn escaped_sql_literal(value: &str) -> String {
        value.replace('\'', "''")
    }

    fn vcf_field_count(&self) -> usize {
        self.schema
            .fields()
            .len()
            .saturating_sub(2 + CACHE_OUTPUT_COLUMNS.len())
    }

    fn vcf_field_names(&self) -> Vec<String> {
        (0..self.vcf_field_count())
            .map(|idx| self.schema.field(idx).name().clone())
            .collect()
    }

    fn parse_json_string_option(json: &str, key: &str) -> Option<String> {
        let needle = format!("\"{key}\"");
        let start = json.find(&needle)?;
        let rest = &json[start + needle.len()..];
        let colon = rest.find(':')?;
        let after_colon = rest[colon + 1..].trim_start();
        let after_quote = after_colon.strip_prefix('"')?;
        let end_quote = after_quote.find('"')?;
        let value = &after_quote[..end_quote];
        if value.is_empty() || value.contains('`') {
            return None;
        }
        Some(value.to_string())
    }

    fn parse_json_bool_option(json: &str, key: &str) -> Option<bool> {
        let needle = format!("\"{key}\"");
        let start = json.find(&needle)?;
        let rest = &json[start + needle.len()..];
        let colon = rest.find(':')?;
        let after_colon = rest[colon + 1..].trim_start();
        if after_colon.starts_with("true") {
            return Some(true);
        }
        if after_colon.starts_with("false") {
            return Some(false);
        }
        None
    }

    async fn resolve_transcript_context_tables(
        &self,
        cache_table: &str,
    ) -> Result<Option<(String, String)>> {
        let mut candidates: Vec<(String, String)> = Vec::new();

        if let Some(options) = self.options_json.as_deref() {
            let tx = Self::parse_json_string_option(options, "transcripts_table");
            let ex = Self::parse_json_string_option(options, "exons_table");
            if let (Some(tx_table), Some(exon_table)) = (tx, ex) {
                candidates.push((tx_table, exon_table));
            }
        }

        candidates.push((
            format!("{cache_table}_transcripts"),
            format!("{cache_table}_exons"),
        ));

        if cache_table != self.cache_source {
            candidates.push((
                format!("{}_transcripts", self.cache_source),
                format!("{}_exons", self.cache_source),
            ));
        }

        for (tx_table, exon_table) in candidates {
            if self.session.table(&tx_table).await.is_ok()
                && self.session.table(&exon_table).await.is_ok()
            {
                return Ok(Some((tx_table, exon_table)));
            }
        }
        Ok(None)
    }

    async fn resolve_optional_context_table(
        &self,
        options_key: &str,
        cache_table: &str,
        suffix: &str,
    ) -> Result<Option<String>> {
        if let Some(options) = self.options_json.as_deref() {
            if let Some(name) = Self::parse_json_string_option(options, options_key) {
                if self.session.table(&name).await.is_ok() {
                    return Ok(Some(name));
                }
            }
        }

        let primary = format!("{cache_table}_{suffix}");
        if self.session.table(&primary).await.is_ok() {
            return Ok(Some(primary));
        }
        if cache_table != self.cache_source {
            let secondary = format!("{}_{}", self.cache_source, suffix);
            if self.session.table(&secondary).await.is_ok() {
                return Ok(Some(secondary));
            }
        }
        Ok(None)
    }

    fn chrom_filter_clause(chroms: &HashSet<String>) -> String {
        if chroms.is_empty() {
            return String::new();
        }
        // Include both "chr"-prefixed and bare variants so mismatched naming
        // between VCF (chr22) and context tables (22) still matches.
        let mut expanded: HashSet<String> = HashSet::new();
        for c in chroms {
            let escaped = c.replace('\'', "''");
            expanded.insert(escaped.clone());
            if let Some(bare) = escaped.strip_prefix("chr") {
                expanded.insert(bare.to_string());
            } else {
                expanded.insert(format!("chr{escaped}"));
            }
        }
        let literals: Vec<String> = expanded.iter().map(|c| format!("'{c}'")).collect();
        format!(" WHERE chrom IN ({})", literals.join(", "))
    }

    async fn load_transcripts(
        &self,
        table: &str,
        chroms: &HashSet<String>,
    ) -> Result<Vec<TranscriptFeature>> {
        let filter = Self::chrom_filter_clause(chroms);
        let query = format!("SELECT * FROM `{table}`{filter}");
        let batches = self.session.sql(&query).await?.collect().await?;
        let mut out = Vec::new();

        for batch in &batches {
            let schema = batch.schema();
            let tx_idx = schema
                .index_of("transcript_id")
                .or_else(|_| schema.index_of("stable_id"))
                .map_err(|_| {
                    DataFusionError::Execution(format!(
                        "annotate_vep(): transcript table '{table}' is missing required column transcript_id (or stable_id)"
                    ))
                })?;
            let chrom_idx = schema.index_of("chrom").map_err(|_| {
                DataFusionError::Execution(format!(
                    "annotate_vep(): transcript table '{table}' is missing required column chrom"
                ))
            })?;
            let start_idx = schema.index_of("start").map_err(|_| {
                DataFusionError::Execution(format!(
                    "annotate_vep(): transcript table '{table}' is missing required column start"
                ))
            })?;
            let end_idx = schema.index_of("end").map_err(|_| {
                DataFusionError::Execution(format!(
                    "annotate_vep(): transcript table '{table}' is missing required column end"
                ))
            })?;
            let strand_idx = schema.index_of("strand").map_err(|_| {
                DataFusionError::Execution(format!(
                    "annotate_vep(): transcript table '{table}' is missing required column strand"
                ))
            })?;
            let biotype_idx = schema.index_of("biotype").map_err(|_| {
                DataFusionError::Execution(format!(
                    "annotate_vep(): transcript table '{table}' is missing required column biotype"
                ))
            })?;
            let cds_start_idx = schema.index_of("cds_start").ok();
            let cds_end_idx = schema.index_of("cds_end").ok();
            let gene_stable_id_idx = schema.index_of("gene_stable_id").ok();
            let gene_symbol_idx = schema.index_of("gene_symbol").ok();
            let gene_symbol_source_idx = schema.index_of("gene_symbol_source").ok();
            let gene_hgnc_id_idx = schema.index_of("gene_hgnc_id").ok();
            let source_idx = schema.index_of("source").ok();
            let version_idx = schema.index_of("version").ok();
            let raw_json_idx = schema.index_of("raw_object_json").ok();
            let cds_start_nf_idx = schema.index_of("cds_start_nf").ok();
            let cds_end_nf_idx = schema.index_of("cds_end_nf").ok();
            let mirna_regions_idx = schema.index_of("mature_mirna_regions").ok();
            // Batch 1 columns.
            let is_canonical_idx = schema.index_of("is_canonical").ok();
            let tsl_idx = schema.index_of("tsl").ok();
            let mane_select_idx = schema.index_of("mane_select").ok();
            let mane_plus_clinical_idx = schema.index_of("mane_plus_clinical").ok();
            let translation_stable_id_idx = schema.index_of("translation_stable_id").ok();
            let gene_phenotype_idx = schema.index_of("gene_phenotype").ok();
            let ccds_idx = schema.index_of("ccds").ok();
            let swissprot_idx = schema.index_of("swissprot").ok();
            let trembl_idx = schema.index_of("trembl").ok();
            let uniparc_idx = schema.index_of("uniparc").ok();
            let uniprot_isoform_idx = schema.index_of("uniprot_isoform").ok();

            for row in 0..batch.num_rows() {
                let Some(transcript_id) = string_at(batch.column(tx_idx).as_ref(), row) else {
                    continue;
                };
                let Some(chrom) = string_at(batch.column(chrom_idx).as_ref(), row) else {
                    continue;
                };
                let Some(start) = int64_at(batch.column(start_idx).as_ref(), row) else {
                    continue;
                };
                let Some(end) = int64_at(batch.column(end_idx).as_ref(), row) else {
                    continue;
                };
                let Some(strand_raw) = int64_at(batch.column(strand_idx).as_ref(), row) else {
                    continue;
                };
                let Some(biotype) = string_at(batch.column(biotype_idx).as_ref(), row) else {
                    continue;
                };

                let strand = if strand_raw >= 0 { 1 } else { -1 };
                let cds_start = cds_start_idx
                    .and_then(|idx| int64_at(batch.column(idx).as_ref(), row))
                    .filter(|v| *v > 0);
                let cds_end = cds_end_idx
                    .and_then(|idx| int64_at(batch.column(idx).as_ref(), row))
                    .filter(|v| *v > 0);

                // Mature miRNA regions from promoted List<Struct<start,end>>
                // column (already genomic coordinates).
                let mature_mirna_regions = if biotype == "miRNA" {
                    mirna_regions_idx
                        .and_then(|idx| read_mirna_regions(batch, idx, row))
                        .unwrap_or_default()
                } else {
                    Vec::new()
                };

                // CDS incompleteness flags from promoted boolean columns.
                let cds_start_nf = cds_start_nf_idx
                    .and_then(|idx| bool_at(batch.column(idx).as_ref(), row))
                    .unwrap_or(false);
                let cds_end_nf = cds_end_nf_idx
                    .and_then(|idx| bool_at(batch.column(idx).as_ref(), row))
                    .unwrap_or(false);
                let flags_str = flags_str_from_bools(cds_start_nf, cds_end_nf);

                let gene_stable_id =
                    gene_stable_id_idx.and_then(|idx| string_at(batch.column(idx).as_ref(), row));
                let gene_symbol =
                    gene_symbol_idx.and_then(|idx| string_at(batch.column(idx).as_ref(), row));
                let gene_symbol_source = gene_symbol_source_idx
                    .and_then(|idx| string_at(batch.column(idx).as_ref(), row));
                let gene_hgnc_id =
                    gene_hgnc_id_idx.and_then(|idx| string_at(batch.column(idx).as_ref(), row));
                let source = source_idx.and_then(|idx| string_at(batch.column(idx).as_ref(), row));
                let version = version_idx
                    .and_then(|idx| int64_at(batch.column(idx).as_ref(), row))
                    .and_then(|v| i32::try_from(v).ok());

                // Batch 1 fields.
                let is_canonical = is_canonical_idx
                    .and_then(|idx| bool_at(batch.column(idx).as_ref(), row))
                    .unwrap_or(false);
                let tsl = tsl_idx
                    .and_then(|idx| int64_at(batch.column(idx).as_ref(), row))
                    .and_then(|v| i32::try_from(v).ok());
                let mane_select = mane_select_idx
                    .and_then(|idx| string_at(batch.column(idx).as_ref(), row));
                let mane_plus_clinical = mane_plus_clinical_idx
                    .and_then(|idx| string_at(batch.column(idx).as_ref(), row));
                let translation_stable_id = translation_stable_id_idx
                    .and_then(|idx| string_at(batch.column(idx).as_ref(), row));
                let gene_phenotype = gene_phenotype_idx
                    .and_then(|idx| bool_at(batch.column(idx).as_ref(), row))
                    .unwrap_or(false);
                let ccds =
                    ccds_idx.and_then(|idx| string_at(batch.column(idx).as_ref(), row));
                let swissprot =
                    swissprot_idx.and_then(|idx| string_at(batch.column(idx).as_ref(), row));
                let trembl =
                    trembl_idx.and_then(|idx| string_at(batch.column(idx).as_ref(), row));
                let uniparc =
                    uniparc_idx.and_then(|idx| string_at(batch.column(idx).as_ref(), row));
                let uniprot_isoform = uniprot_isoform_idx
                    .and_then(|idx| string_at(batch.column(idx).as_ref(), row));

                out.push(TranscriptFeature {
                    transcript_id,
                    chrom,
                    start,
                    end,
                    strand,
                    biotype,
                    cds_start,
                    cds_end,
                    mature_mirna_regions,
                    gene_stable_id,
                    gene_symbol,
                    gene_symbol_source,
                    gene_hgnc_id,
                    source,
                    version,
                    cds_start_nf,
                    cds_end_nf,
                    flags_str,
                    is_canonical,
                    tsl,
                    mane_select,
                    mane_plus_clinical,
                    translation_stable_id,
                    gene_phenotype,
                    ccds,
                    swissprot,
                    trembl,
                    uniparc,
                    uniprot_isoform,
                });
            }
        }

        Ok(out)
    }

    async fn load_exons(&self, table: &str, chroms: &HashSet<String>) -> Result<Vec<ExonFeature>> {
        let has_chrom = self
            .session
            .table(table)
            .await
            .ok()
            .map(|t| {
                t.schema()
                    .as_arrow()
                    .fields()
                    .iter()
                    .any(|f| f.name() == "chrom")
            })
            .unwrap_or(false);
        let filter = if has_chrom {
            Self::chrom_filter_clause(chroms)
        } else {
            String::new()
        };
        let query = format!("SELECT * FROM `{table}`{filter}");
        let batches = self.session.sql(&query).await?.collect().await?;
        let mut out = Vec::new();

        for batch in &batches {
            let schema = batch.schema();
            let tx_idx = schema
                .index_of("transcript_id")
                .or_else(|_| schema.index_of("stable_id"))
                .map_err(|_| {
                    DataFusionError::Execution(format!(
                        "annotate_vep(): exon table '{table}' is missing required column transcript_id (or stable_id)"
                    ))
                })?;
            let exon_idx = schema.index_of("exon_number").map_err(|_| {
                DataFusionError::Execution(format!(
                    "annotate_vep(): exon table '{table}' is missing required column exon_number"
                ))
            })?;
            let start_idx = schema.index_of("start").map_err(|_| {
                DataFusionError::Execution(format!(
                    "annotate_vep(): exon table '{table}' is missing required column start"
                ))
            })?;
            let end_idx = schema.index_of("end").map_err(|_| {
                DataFusionError::Execution(format!(
                    "annotate_vep(): exon table '{table}' is missing required column end"
                ))
            })?;

            for row in 0..batch.num_rows() {
                let Some(transcript_id) = string_at(batch.column(tx_idx).as_ref(), row) else {
                    continue;
                };
                let Some(exon_number_raw) = int64_at(batch.column(exon_idx).as_ref(), row) else {
                    continue;
                };
                let Some(start) = int64_at(batch.column(start_idx).as_ref(), row) else {
                    continue;
                };
                let Some(end) = int64_at(batch.column(end_idx).as_ref(), row) else {
                    continue;
                };

                out.push(ExonFeature {
                    transcript_id,
                    exon_number: i32::try_from(exon_number_raw).unwrap_or(i32::MAX),
                    start,
                    end,
                });
            }
        }

        Ok(out)
    }

    async fn load_translations(
        &self,
        table: &str,
        chroms: &HashSet<String>,
    ) -> Result<Vec<TranslationFeature>> {
        let has_chrom = self
            .session
            .table(table)
            .await
            .ok()
            .map(|t| {
                t.schema()
                    .as_arrow()
                    .fields()
                    .iter()
                    .any(|f| f.name() == "chrom")
            })
            .unwrap_or(false);
        let filter = if has_chrom {
            Self::chrom_filter_clause(chroms)
        } else {
            String::new()
        };
        let query = format!("SELECT * FROM `{table}`{filter}");
        let batches = self.session.sql(&query).await?.collect().await?;
        let mut out = Vec::new();

        for batch in &batches {
            let schema = batch.schema();
            let tx_idx = schema
                .index_of("transcript_id")
                .or_else(|_| schema.index_of("stable_id"))
                .map_err(|_| {
                    DataFusionError::Execution(format!(
                        "annotate_vep(): translation table '{table}' is missing required column transcript_id (or stable_id)"
                    ))
                })?;
            let cds_len_idx = schema
                .index_of("cds_len")
                .or_else(|_| schema.index_of("cds_length"))
                .ok();
            let protein_len_idx = schema.index_of("protein_len").ok();
            let translation_seq_idx = schema.index_of("translation_seq").ok();
            let cds_seq_idx = schema
                .index_of("cds_sequence")
                .or_else(|_| schema.index_of("cds_seq"))
                .or_else(|_| schema.index_of("coding_sequence"))
                .ok();
            let tl_stable_id_idx = schema.index_of("stable_id").ok();
            let tl_version_idx = schema.index_of("version").ok();
            for row in 0..batch.num_rows() {
                let Some(transcript_id) = string_at(batch.column(tx_idx).as_ref(), row) else {
                    continue;
                };
                let cds_len = cds_len_idx
                    .and_then(|idx| int64_at(batch.column(idx).as_ref(), row))
                    .and_then(|v| usize::try_from(v).ok());
                let protein_len = protein_len_idx
                    .and_then(|idx| int64_at(batch.column(idx).as_ref(), row))
                    .and_then(|v| usize::try_from(v).ok());
                let translation_seq =
                    translation_seq_idx.and_then(|idx| string_at(batch.column(idx).as_ref(), row));
                let cds_sequence =
                    cds_seq_idx.and_then(|idx| string_at(batch.column(idx).as_ref(), row));
                let tl_stable_id =
                    tl_stable_id_idx.and_then(|idx| string_at(batch.column(idx).as_ref(), row));
                let tl_version = tl_version_idx
                    .and_then(|idx| int64_at(batch.column(idx).as_ref(), row))
                    .and_then(|v| i32::try_from(v).ok());

                out.push(TranslationFeature {
                    transcript_id,
                    cds_len,
                    protein_len,
                    translation_seq,
                    cds_sequence,
                    stable_id: tl_stable_id,
                    version: tl_version,
                });
            }
        }

        Ok(out)
    }

    async fn load_regulatory_features(
        &self,
        table: &str,
        chroms: &HashSet<String>,
    ) -> Result<Vec<RegulatoryFeature>> {
        let filter = Self::chrom_filter_clause(chroms);
        let query = format!("SELECT * FROM `{table}`{filter}");
        let batches = self.session.sql(&query).await?.collect().await?;
        let mut out = Vec::new();

        for batch in &batches {
            let schema = batch.schema();
            let id_idx = schema
                .index_of("stable_id")
                .or_else(|_| schema.index_of("feature_id"))
                .ok();
            let ft_idx = schema.index_of("feature_type").ok();
            let chrom_idx = schema.index_of("chrom").map_err(|_| {
                DataFusionError::Execution(format!(
                    "annotate_vep(): regulatory table '{table}' is missing required column chrom"
                ))
            })?;
            let start_idx = schema.index_of("start").map_err(|_| {
                DataFusionError::Execution(format!(
                    "annotate_vep(): regulatory table '{table}' is missing required column start"
                ))
            })?;
            let end_idx = schema.index_of("end").map_err(|_| {
                DataFusionError::Execution(format!(
                    "annotate_vep(): regulatory table '{table}' is missing required column end"
                ))
            })?;

            for row in 0..batch.num_rows() {
                let Some(chrom) = string_at(batch.column(chrom_idx).as_ref(), row) else {
                    continue;
                };
                let Some(start) = int64_at(batch.column(start_idx).as_ref(), row) else {
                    continue;
                };
                let Some(end) = int64_at(batch.column(end_idx).as_ref(), row) else {
                    continue;
                };
                let feature_id = id_idx
                    .and_then(|idx| string_at(batch.column(idx).as_ref(), row))
                    .unwrap_or_else(|| "reg".to_string());
                let feature_type =
                    ft_idx.and_then(|idx| string_at(batch.column(idx).as_ref(), row));
                out.push(RegulatoryFeature {
                    feature_id,
                    chrom,
                    start,
                    end,
                    feature_type,
                });
            }
        }

        Ok(out)
    }

    async fn load_motif_features(
        &self,
        table: &str,
        chroms: &HashSet<String>,
    ) -> Result<Vec<MotifFeature>> {
        let filter = Self::chrom_filter_clause(chroms);
        let query = format!("SELECT * FROM `{table}`{filter}");
        let batches = self.session.sql(&query).await?.collect().await?;
        let mut out = Vec::new();

        for batch in &batches {
            let schema = batch.schema();
            let id_idx = schema
                .index_of("motif_id")
                .or_else(|_| schema.index_of("feature_id"))
                .ok();
            let chrom_idx = schema.index_of("chrom").map_err(|_| {
                DataFusionError::Execution(format!(
                    "annotate_vep(): motif table '{table}' is missing required column chrom"
                ))
            })?;
            let start_idx = schema.index_of("start").map_err(|_| {
                DataFusionError::Execution(format!(
                    "annotate_vep(): motif table '{table}' is missing required column start"
                ))
            })?;
            let end_idx = schema.index_of("end").map_err(|_| {
                DataFusionError::Execution(format!(
                    "annotate_vep(): motif table '{table}' is missing required column end"
                ))
            })?;

            for row in 0..batch.num_rows() {
                let Some(chrom) = string_at(batch.column(chrom_idx).as_ref(), row) else {
                    continue;
                };
                let Some(start) = int64_at(batch.column(start_idx).as_ref(), row) else {
                    continue;
                };
                let Some(end) = int64_at(batch.column(end_idx).as_ref(), row) else {
                    continue;
                };
                let motif_id = id_idx
                    .and_then(|idx| string_at(batch.column(idx).as_ref(), row))
                    .unwrap_or_else(|| "motif".to_string());
                out.push(MotifFeature {
                    motif_id,
                    chrom,
                    start,
                    end,
                });
            }
        }

        Ok(out)
    }

    async fn load_mirna_features(
        &self,
        table: &str,
        chroms: &HashSet<String>,
    ) -> Result<Vec<MirnaFeature>> {
        let filter = Self::chrom_filter_clause(chroms);
        let query = format!("SELECT * FROM `{table}`{filter}");
        let batches = self.session.sql(&query).await?.collect().await?;
        let mut out = Vec::new();

        for batch in &batches {
            let schema = batch.schema();
            let id_idx = schema
                .index_of("mirna_id")
                .or_else(|_| schema.index_of("feature_id"))
                .ok();
            let chrom_idx = schema.index_of("chrom").map_err(|_| {
                DataFusionError::Execution(format!(
                    "annotate_vep(): miRNA table '{table}' is missing required column chrom"
                ))
            })?;
            let start_idx = schema.index_of("start").map_err(|_| {
                DataFusionError::Execution(format!(
                    "annotate_vep(): miRNA table '{table}' is missing required column start"
                ))
            })?;
            let end_idx = schema.index_of("end").map_err(|_| {
                DataFusionError::Execution(format!(
                    "annotate_vep(): miRNA table '{table}' is missing required column end"
                ))
            })?;

            for row in 0..batch.num_rows() {
                let Some(chrom) = string_at(batch.column(chrom_idx).as_ref(), row) else {
                    continue;
                };
                let Some(start) = int64_at(batch.column(start_idx).as_ref(), row) else {
                    continue;
                };
                let Some(end) = int64_at(batch.column(end_idx).as_ref(), row) else {
                    continue;
                };
                let mirna_id = id_idx
                    .and_then(|idx| string_at(batch.column(idx).as_ref(), row))
                    .unwrap_or_else(|| "mirna".to_string());
                out.push(MirnaFeature {
                    mirna_id,
                    chrom,
                    start,
                    end,
                });
            }
        }

        Ok(out)
    }

    async fn load_structural_features(
        &self,
        table: &str,
        chroms: &HashSet<String>,
    ) -> Result<Vec<StructuralFeature>> {
        let filter = Self::chrom_filter_clause(chroms);
        let query = format!("SELECT * FROM `{table}`{filter}");
        let batches = self.session.sql(&query).await?.collect().await?;
        let mut out = Vec::new();

        for batch in &batches {
            let schema = batch.schema();
            let id_idx = schema
                .index_of("feature_id")
                .or_else(|_| schema.index_of("stable_id"))
                .ok();
            let kind_idx = schema.index_of("feature_kind").map_err(|_| {
                DataFusionError::Execution(format!(
                    "annotate_vep(): structural table '{table}' is missing required column feature_kind"
                ))
            })?;
            let event_idx = schema.index_of("event_type").map_err(|_| {
                DataFusionError::Execution(format!(
                    "annotate_vep(): structural table '{table}' is missing required column event_type"
                ))
            })?;
            let chrom_idx = schema.index_of("chrom").map_err(|_| {
                DataFusionError::Execution(format!(
                    "annotate_vep(): structural table '{table}' is missing required column chrom"
                ))
            })?;
            let start_idx = schema.index_of("start").map_err(|_| {
                DataFusionError::Execution(format!(
                    "annotate_vep(): structural table '{table}' is missing required column start"
                ))
            })?;
            let end_idx = schema.index_of("end").map_err(|_| {
                DataFusionError::Execution(format!(
                    "annotate_vep(): structural table '{table}' is missing required column end"
                ))
            })?;

            for row in 0..batch.num_rows() {
                let Some(chrom) = string_at(batch.column(chrom_idx).as_ref(), row) else {
                    continue;
                };
                let Some(start) = int64_at(batch.column(start_idx).as_ref(), row) else {
                    continue;
                };
                let Some(end) = int64_at(batch.column(end_idx).as_ref(), row) else {
                    continue;
                };
                let Some(kind_raw) = string_at(batch.column(kind_idx).as_ref(), row) else {
                    continue;
                };
                let Some(event_raw) = string_at(batch.column(event_idx).as_ref(), row) else {
                    continue;
                };
                let Some(feature_kind) = parse_sv_feature_kind(&kind_raw) else {
                    continue;
                };
                let Some(event_kind) = parse_sv_event_kind(&event_raw) else {
                    continue;
                };
                let feature_id = id_idx
                    .and_then(|idx| string_at(batch.column(idx).as_ref(), row))
                    .unwrap_or_else(|| "sv".to_string());
                out.push(StructuralFeature {
                    feature_id,
                    chrom,
                    start,
                    end,
                    feature_kind,
                    event_kind,
                });
            }
        }

        Ok(out)
    }

    async fn scan_with_transcript_engine(
        &self,
        state: &dyn Session,
        projection: Option<&Vec<usize>>,
        lookup_sql: &str,
        requested_columns: &[&str],
        transcripts_table: Option<&str>,
        exons_table: Option<&str>,
        translations_table: Option<&str>,
        regulatory_table: Option<&str>,
        motif_table: Option<&str>,
        mirna_table: Option<&str>,
        sv_table: Option<&str>,
    ) -> Result<Arc<dyn ExecutionPlan>> {
        let vcf_names = self.vcf_field_names();
        let mut select_exprs = Vec::new();
        for name in &vcf_names {
            select_exprs.push(format!("l.`{name}` AS `{name}`"));
        }
        for name in requested_columns {
            select_exprs.push(format!("l.`cache_{name}` AS `cache_{name}`"));
        }
        let base_query = format!(
            "SELECT {} FROM {} AS l",
            select_exprs.join(", "),
            lookup_sql
        );
        let base_batches = self.session.sql(&base_query).await?.collect().await?;

        // Check if there are any cache misses (rows without cached most_severe).
        // Only load context tables if we actually need the transcript engine.
        let has_cache_misses = base_batches.iter().any(|batch| {
            let Ok(idx) = batch.schema().index_of("cache_most_severe_consequence") else {
                return true; // Column missing — all rows are misses.
            };
            let col = batch.column(idx);
            (0..batch.num_rows()).any(|row| string_at(col.as_ref(), row).is_none())
        });

        let (transcripts, exons, translations, regulatory, motifs, mirnas, structural) =
            if has_cache_misses {
                // Extract distinct chrom values from VCF batches for predicate pushdown
                let vcf_chroms: HashSet<String> = {
                    let mut set = HashSet::new();
                    for batch in &base_batches {
                        if let Ok(idx) = batch.schema().index_of("chrom") {
                            let col = batch.column(idx);
                            for row in 0..batch.num_rows() {
                                if let Some(v) = string_at(col.as_ref(), row) {
                                    set.insert(v);
                                }
                            }
                        }
                    }
                    set
                };

                let merged = self
                    .options_json
                    .as_deref()
                    .and_then(|opts| Self::parse_json_bool_option(opts, "merged"))
                    .unwrap_or(false);

                let tx = if let Some(table) = transcripts_table {
                    self.load_transcripts(table, &vcf_chroms)
                        .await?
                        .into_iter()
                        .filter(|t| is_vep_transcript(&t.transcript_id, merged))
                        .collect::<Vec<_>>()
                } else {
                    Vec::new()
                };
                let ex = if let Some(table) = exons_table {
                    self.load_exons(table, &vcf_chroms).await?
                } else {
                    Vec::new()
                };
                let tl = if let Some(table) = translations_table {
                    self.load_translations(table, &vcf_chroms).await?
                } else {
                    Vec::new()
                };
                let rg = if let Some(table) = regulatory_table {
                    self.load_regulatory_features(table, &vcf_chroms).await?
                } else {
                    Vec::new()
                };
                let mt = if let Some(table) = motif_table {
                    self.load_motif_features(table, &vcf_chroms).await?
                } else {
                    Vec::new()
                };
                let mi = if let Some(table) = mirna_table {
                    self.load_mirna_features(table, &vcf_chroms).await?
                } else {
                    Vec::new()
                };
                let st = if let Some(table) = sv_table {
                    self.load_structural_features(table, &vcf_chroms).await?
                } else {
                    Vec::new()
                };
                (tx, ex, tl, rg, mt, mi, st)
            } else {
                (
                    Vec::new(),
                    Vec::new(),
                    Vec::new(),
                    Vec::new(),
                    Vec::new(),
                    Vec::new(),
                    Vec::new(),
                )
            };
        let engine = TranscriptConsequenceEngine::default();
        let ctx = PreparedContext::new(
            &transcripts,
            &exons,
            &translations,
            &regulatory,
            &motifs,
            &mirnas,
            &structural,
        );

        let mut annotated_batches = Vec::with_capacity(base_batches.len());
        for batch in &base_batches {
            annotated_batches
                .push(self.annotate_batch_with_transcript_engine(batch, &engine, &ctx)?);
        }

        let projected_batches = if let Some(indices) = projection {
            let mut out = Vec::with_capacity(annotated_batches.len());
            for batch in annotated_batches {
                out.push(batch.project(indices)?);
            }
            out
        } else {
            annotated_batches
        };
        let projected_schema = if let Some(indices) = projection {
            Arc::new(self.schema.project(indices)?)
        } else {
            self.schema.clone()
        };

        let mem = MemTable::try_new(projected_schema, vec![projected_batches])?;
        mem.scan(state, None, &[], None).await
    }

    fn annotate_batch_with_transcript_engine(
        &self,
        batch: &RecordBatch,
        engine: &TranscriptConsequenceEngine,
        ctx: &PreparedContext<'_>,
    ) -> Result<RecordBatch> {
        let schema = batch.schema();
        let chrom_idx = schema.index_of("chrom").map_err(|_| {
            DataFusionError::Execution(
                "annotate_vep(): input VCF row is missing required chrom column".to_string(),
            )
        })?;
        let start_idx = schema.index_of("start").map_err(|_| {
            DataFusionError::Execution(
                "annotate_vep(): input VCF row is missing required start column".to_string(),
            )
        })?;
        let end_idx = schema.index_of("end").map_err(|_| {
            DataFusionError::Execution(
                "annotate_vep(): input VCF row is missing required end column".to_string(),
            )
        })?;
        let ref_idx = schema.index_of("ref").map_err(|_| {
            DataFusionError::Execution(
                "annotate_vep(): input VCF row is missing required ref column".to_string(),
            )
        })?;
        let alt_idx = schema.index_of("alt").map_err(|_| {
            DataFusionError::Execution(
                "annotate_vep(): input VCF row is missing required alt column".to_string(),
            )
        })?;
        let variation_name_idx = schema.index_of("cache_variation_name").ok();
        let cached_csq_idx = schema.index_of("cache_consequence_types").ok();
        let cached_most_idx = schema.index_of("cache_most_severe_consequence").ok();
        let merged = self
            .options_json
            .as_deref()
            .and_then(|opts| Self::parse_json_bool_option(opts, "merged"))
            .unwrap_or(false);
        let flags = VepFlags::from_options_json(self.options_json.as_deref());

        // Resolve cache column indices for Batch 3 AF columns.
        let af_col_indices: Vec<Option<usize>> = AF_COLUMNS
            .iter()
            .map(|col| schema.index_of(&format!("cache_{}", col.cache_col)).ok())
            .collect();
        // Clinical/frequency metadata indices.
        let clin_sig_idx = schema.index_of("cache_clin_sig").ok();
        let somatic_idx = schema.index_of("cache_somatic").ok();
        let pheno_idx = schema.index_of("cache_phenotype_or_disease").ok();
        let pubmed_idx = schema.index_of("cache_pubmed").ok();

        let mut csq_builder = StringBuilder::with_capacity(batch.num_rows(), batch.num_rows() * 40);
        let mut most_builder =
            StringBuilder::with_capacity(batch.num_rows(), batch.num_rows() * 16);

        for row in 0..batch.num_rows() {
            let Some(chrom) = string_at(batch.column(chrom_idx).as_ref(), row) else {
                csq_builder.append_null();
                most_builder.append_null();
                continue;
            };
            let Some(alt_allele) = string_at(batch.column(alt_idx).as_ref(), row) else {
                csq_builder.append_null();
                most_builder.append_null();
                continue;
            };

            // VEP-style allele minimization: strip shared prefix and suffix between REF and ALT.
            let ref_al = string_at(batch.column(ref_idx).as_ref(), row).unwrap_or_default();
            let (vep_ref, vep_allele) = vcf_to_vep_allele(&ref_al, &alt_allele);
            let variant_class = classify_variant(&vep_ref, &vep_allele);

            // Cache-hit fast path: use pre-computed consequence from variation cache.
            let cached_most =
                cached_most_idx.and_then(|idx| string_at(batch.column(idx).as_ref(), row));
            let cached_csq =
                cached_csq_idx.and_then(|idx| string_at(batch.column(idx).as_ref(), row));

            let variation_name = variation_name_idx
                .and_then(|idx| string_at(batch.column(idx).as_ref(), row))
                .unwrap_or_default();

            // --- Batch 3: per-variant fields (same for every transcript entry) ---
            let existing_var = if flags.check_existing {
                &variation_name
            } else {
                ""
            };

            // Extract allele frequencies from cache columns.
            let af_raw: Vec<String> = af_col_indices
                .iter()
                .map(|opt_idx| {
                    opt_idx
                        .and_then(|idx| string_at(batch.column(idx).as_ref(), row))
                        .unwrap_or_default()
                })
                .collect();
            // Parse allele:freq for each AF column. Build two vectors:
            // 1. af_csq_values: what goes into CSQ fields (empty for gnomAD sub-pops)
            // 2. max_af_entries: sub-population (name, freq) pairs for MAX_AF computation
            let mut af_csq_values: Vec<String> = Vec::with_capacity(AF_COLUMNS.len());
            let mut max_af_entries: Vec<(&str, &str)> = Vec::new();
            for (i, col) in AF_COLUMNS.iter().enumerate() {
                let freq = if flags.af_group_enabled(col.flag_group) {
                    extract_af_for_allele(&af_raw[i], &vep_allele)
                } else {
                    ""
                };
                // CSQ output: only emit for columns VEP emits in offline/cache mode.
                let csq_val = if col.emit_in_csq { freq } else { "" };
                af_csq_values.push(csq_val.to_string());
                // MAX_AF: only include sub-populations (not globals).
                if let Some(pop_name) = col.max_af_pop {
                    if !freq.is_empty() {
                        max_af_entries.push((pop_name, freq));
                    }
                }
            }
            let (max_af, max_af_pops) = if flags.max_af {
                compute_max_af(&max_af_entries)
            } else {
                (String::new(), String::new())
            };

            // Clinical / phenotype fields.
            let clin_sig = if flags.check_existing {
                clin_sig_idx
                    .and_then(|idx| string_at(batch.column(idx).as_ref(), row))
                    .unwrap_or_default()
            } else {
                String::new()
            };
            let somatic_val = if flags.check_existing {
                somatic_idx
                    .and_then(|idx| int64_at(batch.column(idx).as_ref(), row))
                    .map(|v| if v != 0 { "1" } else { "" })
                    .unwrap_or("")
            } else {
                ""
            };
            let pheno_val = if flags.check_existing {
                pheno_idx
                    .and_then(|idx| int64_at(batch.column(idx).as_ref(), row))
                    .map(|v| if v != 0 { "1" } else { "" })
                    .unwrap_or("")
            } else {
                ""
            };
            let pubmed_val = if flags.pubmed {
                pubmed_idx
                    .and_then(|idx| string_at(batch.column(idx).as_ref(), row))
                    .unwrap_or_default()
            } else {
                String::new()
            };

            // Build the 33-field Batch 3 suffix (positions 41-73) shared across all transcripts.
            let batch3_suffix = format!(
                "{}|{}|{}|{clin_sig}|{somatic_val}|{pheno_val}|{pubmed_val}",
                af_csq_values.join("|"),
                max_af,
                max_af_pops,
            );

            let (csq_string, most_str) = if let Some(most_val) = &cached_most {
                // Cache hit — produce single CSQ entry with empty transcript fields.
                let csq_val = cached_csq.unwrap_or_default();
                let impact = SoTerm::from_str(most_val)
                    .map(|t| impact_label(t.impact()))
                    .unwrap_or_else(|| impact_label(SoImpact::Modifier));
                // 74-field CSQ: 29 base + 12 Batch 1 + 33 Batch 3.
                let csq_entry = format!(
                    "{vep_allele}|{csq_val}|{impact}|||||||||||||||||{existing_var}||||||||||\
                     {variant_class}||||||||||||{batch3_suffix}"
                );
                (csq_entry, most_val.clone())
            } else {
                // Cache miss — compute via transcript engine and produce per-transcript CSQ.
                let Some(start) = int64_at(batch.column(start_idx).as_ref(), row) else {
                    csq_builder.append_null();
                    most_builder.append_null();
                    continue;
                };
                let Some(end) = int64_at(batch.column(end_idx).as_ref(), row) else {
                    csq_builder.append_null();
                    most_builder.append_null();
                    continue;
                };
                let Some(ref_allele) = string_at(batch.column(ref_idx).as_ref(), row) else {
                    csq_builder.append_null();
                    most_builder.append_null();
                    continue;
                };

                let variant = VariantInput::from_vcf(
                    chrom.clone(),
                    start,
                    end,
                    ref_allele,
                    alt_allele.clone(),
                );
                let assignments = engine.evaluate_variant_prepared(&variant, ctx);

                // Derive most_severe from all assignments.
                let mut all_terms =
                    TranscriptConsequenceEngine::collapse_variant_terms(&assignments);
                if all_terms.is_empty() {
                    all_terms.push(SoTerm::SequenceVariant);
                }
                let most = most_severe_term(all_terms.iter()).unwrap_or(SoTerm::SequenceVariant);
                let most_str = most.as_str().to_string();

                // Build per-transcript CSQ entries, comma-separated.
                let mut csq_entries: Vec<String> = Vec::with_capacity(assignments.len());
                for tc in &assignments {
                    let terms_str = tc
                        .terms
                        .iter()
                        .map(|t| t.as_str())
                        .collect::<Vec<_>>()
                        .join("&");
                    let tc_impact = most_severe_term(tc.terms.iter())
                        .map(|t| impact_label(t.impact()))
                        .unwrap_or_else(|| impact_label(SoImpact::Modifier));
                    let feature_type = tc.feature_type.as_str();
                    let feature = tc.transcript_id.as_deref().unwrap_or("");
                    // Look up transcript metadata via index (zero-copy).
                    let tx_opt = tc.transcript_idx.map(|idx| &ctx.transcripts[idx]);
                    let (symbol, gene, biotype_tx, strand_str, symbol_source, hgnc_id, source) =
                        if let Some(tx) = tx_opt {
                            (
                                tx.gene_symbol.as_deref().unwrap_or(""),
                                tx.gene_stable_id.as_deref().unwrap_or(""),
                                tx.biotype.as_str(),
                                if tx.strand >= 0 { "1" } else { "-1" },
                                tx.gene_symbol_source.as_deref().unwrap_or(""),
                                tx.gene_hgnc_id.as_deref().unwrap_or(""),
                                tx.source.as_deref().unwrap_or(""),
                            )
                        } else {
                            ("", "", "", "", "", "", "")
                        };
                    let biotype = tc.biotype_override.as_deref().unwrap_or(biotype_tx);
                    let exon = tc.exon_str.as_deref().unwrap_or("");
                    let intron = tc.intron_str.as_deref().unwrap_or("");
                    let cdna_pos = tc.cdna_position.as_deref().unwrap_or("");
                    let cds_pos = tc.cds_position.as_deref().unwrap_or("");
                    let protein_pos = tc.protein_position.as_deref().unwrap_or("");
                    let amino_acids = tc.amino_acids.as_deref().unwrap_or("");
                    let codons_str = tc.codons.as_deref().unwrap_or("");
                    let distance = tc.distance.map(|d| d.to_string()).unwrap_or_default();
                    let tc_flags = tc.flags.as_deref().unwrap_or("");
                    let hgvsc = "";
                    let hgvsp = "";
                    let source_val = if merged { source } else { "" };

                    // Batch 1 fields from transcript metadata.
                    let canonical = tx_opt
                        .map(|tx| if tx.is_canonical { "YES" } else { "" })
                        .unwrap_or("");
                    let tsl_str = tx_opt
                        .and_then(|tx| tx.tsl)
                        .map(|v| v.to_string())
                        .unwrap_or_default();
                    let mane_select = tx_opt
                        .and_then(|tx| tx.mane_select.as_deref())
                        .unwrap_or("");
                    let mane_plus = tx_opt
                        .and_then(|tx| tx.mane_plus_clinical.as_deref())
                        .unwrap_or("");
                    let ensp = tx_opt
                        .and_then(|tx| tx.translation_stable_id.as_deref())
                        .unwrap_or("");
                    let gene_pheno = tx_opt
                        .map(|tx| if tx.gene_phenotype { "1" } else { "" })
                        .unwrap_or("");
                    let ccds = tx_opt
                        .and_then(|tx| tx.ccds.as_deref())
                        .unwrap_or("");
                    let swissprot = tx_opt
                        .and_then(|tx| tx.swissprot.as_deref())
                        .unwrap_or("");
                    let trembl = tx_opt
                        .and_then(|tx| tx.trembl.as_deref())
                        .unwrap_or("");
                    let uniparc = tx_opt
                        .and_then(|tx| tx.uniparc.as_deref())
                        .unwrap_or("");
                    let uniprot_isoform = tx_opt
                        .and_then(|tx| tx.uniprot_isoform.as_deref())
                        .unwrap_or("");

                    // 74-field CSQ: 29 base + 12 Batch 1 + 33 Batch 3.
                    csq_entries.push(format!(
                        "{vep_allele}|{terms_str}|{tc_impact}|{symbol}|{gene}|{feature_type}|{feature}|{biotype}|\
                         {exon}|{intron}|{hgvsc}|{hgvsp}|\
                         {cdna_pos}|{cds_pos}|{protein_pos}|{amino_acids}|{codons_str}|\
                         {existing_var}|{distance}|{strand_str}|{tc_flags}|{symbol_source}|{hgnc_id}|\
                         |||||{source_val}|\
                         {variant_class}|{canonical}|{tsl_str}|{mane_select}|{mane_plus}|\
                         {ensp}|{gene_pheno}|{ccds}|{swissprot}|{trembl}|{uniparc}|{uniprot_isoform}|\
                         {batch3_suffix}"
                    ));
                }
                if csq_entries.is_empty() {
                    let impact = impact_label(SoImpact::Modifier);
                    csq_entries.push(format!(
                        "{vep_allele}|sequence_variant|{impact}|||||||||||||||||{existing_var}||||||||||\
                         {variant_class}||||||||||||{batch3_suffix}"
                    ));
                }
                (csq_entries.join(","), most_str)
            };

            csq_builder.append_value(&csq_string);
            most_builder.append_value(&most_str);
        }

        let mut out_cols =
            Vec::with_capacity(self.vcf_field_count() + 2 + CACHE_OUTPUT_COLUMNS.len());
        for name in self.vcf_field_names() {
            let idx = schema.index_of(&name).map_err(|_| {
                DataFusionError::Execution(format!(
                    "annotate_vep(): expected VCF output column '{name}' missing from intermediate lookup output"
                ))
            })?;
            out_cols.push(batch.column(idx).clone());
        }
        out_cols.push(Arc::new(csq_builder.finish()));
        out_cols.push(Arc::new(most_builder.finish()));

        // Pass through extra cache annotation columns.
        for &col_name in CACHE_OUTPUT_COLUMNS {
            let col_name_in_batch = format!("cache_{col_name}");
            if let Ok(idx) = schema.index_of(&col_name_in_batch) {
                let col = batch.column(idx);
                let mut builder =
                    StringBuilder::with_capacity(batch.num_rows(), batch.num_rows() * 8);
                for row in 0..batch.num_rows() {
                    match string_at(col.as_ref(), row) {
                        Some(v) => builder.append_value(&v),
                        None => builder.append_null(),
                    }
                }
                out_cols.push(Arc::new(builder.finish()));
            } else {
                out_cols.push(new_null_array(&DataType::Utf8, batch.num_rows()));
            }
        }

        Ok(RecordBatch::try_new(self.schema.clone(), out_cols)?)
    }
}

fn impact_label(impact: SoImpact) -> &'static str {
    match impact {
        SoImpact::High => "HIGH",
        SoImpact::Moderate => "MODERATE",
        SoImpact::Low => "LOW",
        SoImpact::Modifier => "MODIFIER",
    }
}

/// Classify a variant from VEP-minimized REF/ALT alleles.
fn classify_variant(vep_ref: &str, vep_alt: &str) -> &'static str {
    match (vep_ref, vep_alt) {
        ("-", a) if !a.is_empty() => "insertion",
        (r, "-") if !r.is_empty() => "deletion",
        (r, a) if r.len() == 1 && a.len() == 1 => "SNV",
        (r, a) if r.len() == a.len() => "substitution",
        _ => "indel",
    }
}

fn parse_sv_feature_kind(value: &str) -> Option<SvFeatureKind> {
    match value.to_ascii_lowercase().as_str() {
        "transcript" | "tx" => Some(SvFeatureKind::Transcript),
        "regulatory" | "reg" => Some(SvFeatureKind::Regulatory),
        "tfbs" | "motif" => Some(SvFeatureKind::Tfbs),
        "feature" | "generic" => Some(SvFeatureKind::Generic),
        _ => None,
    }
}

fn parse_sv_event_kind(value: &str) -> Option<SvEventKind> {
    match value.to_ascii_lowercase().as_str() {
        "ablation" | "deletion" | "del" => Some(SvEventKind::Ablation),
        "amplification" | "duplication" | "dup" | "amp" => Some(SvEventKind::Amplification),
        "elongation" | "elongate" => Some(SvEventKind::Elongation),
        "truncation" | "truncate" => Some(SvEventKind::Truncation),
        _ => None,
    }
}

/// Parse mature miRNA genomic regions from the `raw_object_json` transcript
/// attribute.  VEP stores miRNA cDNA coordinates in the transcript's attribute
/// array as `{code: "miRNA", value: "42-59"}`.  We map those cDNA coords to
/// genomic coordinates using the strand and transcript boundaries.
///
/// miRNA transcripts are almost always single-exon, so the mapping is trivial:
/// - Plus strand:  `genomic = tx.start + cdna - 1`
/// - Minus strand: `genomic_start = tx.end - cdna_end + 1`, `genomic_end = tx.end - cdna_start + 1`
/// Reconstruct `FLAGS` string from promoted boolean columns using a canonical
/// order (cds_start_NF before cds_end_NF).  VEP's encounter order varies per
/// transcript but our golden benchmark already normalizes term order, so
/// canonical ordering is correct.
fn flags_str_from_bools(cds_start_nf: bool, cds_end_nf: bool) -> Option<String> {
    match (cds_start_nf, cds_end_nf) {
        (true, true) => Some("cds_start_NF&cds_end_NF".to_string()),
        (true, false) => Some("cds_start_NF".to_string()),
        (false, true) => Some("cds_end_NF".to_string()),
        (false, false) => None,
    }
}

/// Read mature miRNA genomic regions from a promoted `List<Struct<start,end>>`
/// column.  Returns `None` if the cell is NULL (letting the caller fall back
/// to JSON parsing if needed).
fn read_mirna_regions(batch: &RecordBatch, col_idx: usize, row: usize) -> Option<Vec<(i64, i64)>> {
    let col = batch.column(col_idx);
    if col.is_null(row) {
        return None;
    }
    let list_arr = col.as_any().downcast_ref::<ListArray>()?;
    let offsets = list_arr.offsets();
    let start_off = offsets[row] as usize;
    let end_off = offsets[row + 1] as usize;
    if start_off == end_off {
        return Some(Vec::new());
    }
    let values = list_arr.values();
    let struct_arr = values.as_struct();
    let starts = struct_arr.column_by_name("start")?;
    let ends = struct_arr.column_by_name("end")?;

    let mut regions = Vec::with_capacity(end_off - start_off);
    for i in start_off..end_off {
        let s = int64_at(starts.as_ref(), i)?;
        let e = int64_at(ends.as_ref(), i)?;
        regions.push((s, e));
    }
    Some(regions)
}

fn bool_at(array: &dyn Array, row: usize) -> Option<bool> {
    if array.is_null(row) {
        return None;
    }
    if let Some(arr) = array.as_any().downcast_ref::<BooleanArray>() {
        return Some(arr.value(row));
    }
    None
}

fn int64_at(array: &dyn Array, row: usize) -> Option<i64> {
    if array.is_null(row) {
        return None;
    }
    if let Some(arr) = array.as_any().downcast_ref::<Int64Array>() {
        return Some(arr.value(row));
    }
    if let Some(arr) = array.as_any().downcast_ref::<Int32Array>() {
        return Some(arr.value(row) as i64);
    }
    if let Some(arr) = array.as_any().downcast_ref::<Int16Array>() {
        return Some(arr.value(row) as i64);
    }
    if let Some(arr) = array.as_any().downcast_ref::<Int8Array>() {
        return Some(arr.value(row) as i64);
    }
    if let Some(arr) = array.as_any().downcast_ref::<UInt64Array>() {
        return i64::try_from(arr.value(row)).ok();
    }
    if let Some(arr) = array.as_any().downcast_ref::<UInt32Array>() {
        return Some(arr.value(row) as i64);
    }
    if let Some(arr) = array.as_any().downcast_ref::<UInt16Array>() {
        return Some(arr.value(row) as i64);
    }
    if let Some(arr) = array.as_any().downcast_ref::<UInt8Array>() {
        return Some(arr.value(row) as i64);
    }
    None
}

fn string_at(array: &dyn Array, row: usize) -> Option<String> {
    if array.is_null(row) {
        return None;
    }
    if let Some(arr) = array.as_any().downcast_ref::<StringArray>() {
        return Some(arr.value(row).to_string());
    }
    if let Some(arr) = array.as_any().downcast_ref::<StringViewArray>() {
        return Some(arr.value(row).to_string());
    }
    if let Some(arr) = array.as_any().downcast_ref::<LargeStringArray>() {
        return Some(arr.value(row).to_string());
    }
    if let Some(arr) = array.as_any().downcast_ref::<Float64Array>() {
        return Some(arr.value(row).to_string());
    }
    if let Some(arr) = array.as_any().downcast_ref::<Float32Array>() {
        return Some(arr.value(row).to_string());
    }
    if let Some(arr) = array.as_any().downcast_ref::<Int64Array>() {
        return Some(arr.value(row).to_string());
    }
    if let Some(arr) = array.as_any().downcast_ref::<Int32Array>() {
        return Some(arr.value(row).to_string());
    }
    if let Some(arr) = array.as_any().downcast_ref::<UInt64Array>() {
        return Some(arr.value(row).to_string());
    }
    if let Some(arr) = array.as_any().downcast_ref::<UInt32Array>() {
        return Some(arr.value(row).to_string());
    }
    None
}

impl Debug for AnnotateProvider {
    fn fmt(&self, f: &mut Formatter<'_>) -> std::fmt::Result {
        write!(
            f,
            "AnnotateProvider {{ vcf: {}, cache_source: {}, backend: {} }}",
            self.vcf_table,
            self.cache_source,
            self.backend.as_str()
        )
    }
}

#[async_trait]
impl TableProvider for AnnotateProvider {
    fn as_any(&self) -> &dyn Any {
        self
    }

    fn schema(&self) -> SchemaRef {
        self.schema.clone()
    }

    fn table_type(&self) -> TableType {
        TableType::Temporary
    }

    async fn scan(
        &self,
        state: &dyn Session,
        projection: Option<&Vec<usize>>,
        _filters: &[Expr],
        _limit: Option<usize>,
    ) -> Result<Arc<dyn ExecutionPlan>> {
        let _store = build_store(self.backend, self.cache_source.clone());
        let cache_table = self.resolve_cache_table_name().await?;

        let cache_schema = self
            .session
            .table(&cache_table)
            .await?
            .schema()
            .as_arrow()
            .clone();
        let available_cache_columns: HashSet<String> = cache_schema
            .fields()
            .iter()
            .map(|f| f.name().clone())
            .collect();

        let mut preferred_columns = vec!["consequence_types", "most_severe_consequence"];
        for &col in CACHE_OUTPUT_COLUMNS {
            if !preferred_columns.contains(&col) {
                preferred_columns.push(col);
            }
        }
        let requested_columns: Vec<&str> = preferred_columns
            .iter()
            .copied()
            .filter(|name| available_cache_columns.contains(*name))
            .collect();
        let requested_columns_sql = requested_columns.join(",");

        let vcf_table_lit = Self::escaped_sql_literal(&self.vcf_table);
        let cache_table_lit = Self::escaped_sql_literal(&cache_table);
        let columns_lit = Self::escaped_sql_literal(&requested_columns_sql);
        // Extended probes use interval-overlap join to handle VEP-style
        // indel coordinate encodings. Enabled by default for backward
        // compatibility. When input is pre-normalized (e.g. bcftools norm),
        // set "extended_probes":false in options_json for faster equi-join.
        let extended_probes = self
            .options_json
            .as_deref()
            .and_then(|opts| Self::parse_json_bool_option(opts, "extended_probes"))
            .unwrap_or(true);
        let lookup_sql = format!(
            "lookup_variants('{vcf_table_lit}', '{cache_table_lit}', '{columns_lit}', 'exact', {extended_probes})"
        );

        let transcript_pair = self.resolve_transcript_context_tables(&cache_table).await?;
        let translations_table = self
            .resolve_optional_context_table("translations_table", &cache_table, "translations")
            .await?;
        let regulatory_table = self
            .resolve_optional_context_table("regulatory_table", &cache_table, "regulatory_features")
            .await?;
        let motif_table = self
            .resolve_optional_context_table("motif_table", &cache_table, "motif_features")
            .await?;
        let mirna_table = self
            .resolve_optional_context_table("mirna_table", &cache_table, "mirna_features")
            .await?;
        let sv_table = self
            .resolve_optional_context_table("sv_table", &cache_table, "sv_features")
            .await?;

        if transcript_pair.is_some()
            || translations_table.is_some()
            || regulatory_table.is_some()
            || motif_table.is_some()
            || mirna_table.is_some()
            || sv_table.is_some()
        {
            let (tx_table, ex_table) = transcript_pair
                .as_ref()
                .map(|(tx, ex)| (Some(tx.as_str()), Some(ex.as_str())))
                .unwrap_or((None, None));
            return self
                .scan_with_transcript_engine(
                    state,
                    projection,
                    &lookup_sql,
                    &requested_columns,
                    tx_table,
                    ex_table,
                    translations_table.as_deref(),
                    regulatory_table.as_deref(),
                    motif_table.as_deref(),
                    mirna_table.as_deref(),
                    sv_table.as_deref(),
                )
                .await;
        }

        let has_col = |name: &str| requested_columns.contains(&name);
        let present_checks = [
            ("variation_name", has_col("variation_name")),
            ("clin_sig", has_col("clin_sig")),
            ("AF", has_col("AF")),
            ("somatic", has_col("somatic")),
            ("phenotype_or_disease", has_col("phenotype_or_disease")),
            ("pubmed", has_col("pubmed")),
        ]
        .iter()
        .filter_map(|(name, present)| {
            if *present {
                Some(format!("l.`cache_{name}` IS NOT NULL"))
            } else {
                None
            }
        })
        .collect::<Vec<_>>();
        let present_condition_sql = if present_checks.is_empty() {
            "FALSE".to_string()
        } else {
            present_checks.join(" OR ")
        };

        let var_expr = if has_col("variation_name") {
            "COALESCE(CAST(l.`cache_variation_name` AS VARCHAR), '')"
        } else {
            "''"
        };
        let clin_expr = if has_col("clin_sig") {
            "COALESCE(CAST(l.`cache_clin_sig` AS VARCHAR), '')"
        } else {
            "''"
        };
        let af_expr = if has_col("AF") {
            "COALESCE(CAST(l.`cache_AF` AS VARCHAR), '')"
        } else {
            "''"
        };

        let csq_expr = format!(
            "CASE WHEN {present_condition_sql} THEN \
             CONCAT(COALESCE(CAST(l.`alt` AS VARCHAR), ''), \
                    '|sequence_variant|MODIFIER|', {var_expr}, '|', {clin_expr}, '|', {af_expr}) \
             ELSE CAST(NULL AS VARCHAR) END"
        );
        let most_severe_expr = format!(
            "CASE WHEN {present_condition_sql} THEN 'sequence_variant' ELSE CAST(NULL AS VARCHAR) END"
        );

        let total_fields = self.schema.fields().len();
        let vcf_fields = self.vcf_field_count();
        let csq_idx = vcf_fields;
        let most_severe_idx = vcf_fields + 1;
        let cache_col_start = vcf_fields + 2;
        let projected_indices: Vec<usize> = projection
            .cloned()
            .unwrap_or_else(|| (0..total_fields).collect());

        let mut projected_exprs = Vec::with_capacity(projected_indices.len());

        for idx in projected_indices {
            if idx >= total_fields {
                return Err(DataFusionError::Execution(format!(
                    "annotate_vep(): projection index {idx} out of bounds for schema with {total_fields} fields"
                )));
            }

            if idx < vcf_fields {
                let name = self.schema.field(idx).name().clone();
                projected_exprs.push(format!("l.`{name}` AS `{name}`"));
            } else if idx == csq_idx {
                projected_exprs.push(format!("{csq_expr} AS `csq`"));
            } else if idx == most_severe_idx {
                projected_exprs.push(format!("{most_severe_expr} AS `most_severe_consequence`"));
            } else if idx >= cache_col_start {
                let col_name = CACHE_OUTPUT_COLUMNS[idx - cache_col_start];
                projected_exprs.push(format!(
                    "CAST(l.`cache_{col_name}` AS VARCHAR) AS `{col_name}`"
                ));
            }
        }

        // Fallback path: when transcript/exon context tables are unavailable.
        let _ = &self.options_json;
        let query = format!(
            "SELECT {} FROM {} AS l",
            projected_exprs.join(", "),
            lookup_sql
        );

        let df = self.session.sql(&query).await?;
        df.create_physical_plan().await
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_flags_str_both_true() {
        assert_eq!(
            flags_str_from_bools(true, true),
            Some("cds_start_NF&cds_end_NF".to_string())
        );
    }

    #[test]
    fn test_flags_str_start_only() {
        assert_eq!(
            flags_str_from_bools(true, false),
            Some("cds_start_NF".to_string())
        );
    }

    #[test]
    fn test_flags_str_end_only() {
        assert_eq!(
            flags_str_from_bools(false, true),
            Some("cds_end_NF".to_string())
        );
    }

    #[test]
    fn test_flags_str_neither() {
        assert_eq!(flags_str_from_bools(false, false), None);
    }

    // --- extract_af_for_allele tests ---

    #[test]
    fn test_extract_af_single_allele() {
        assert_eq!(extract_af_for_allele("T:0.9301", "T"), "0.9301");
    }

    #[test]
    fn test_extract_af_multi_allele() {
        assert_eq!(extract_af_for_allele("A:0.006,G:0.994", "G"), "0.994");
        assert_eq!(extract_af_for_allele("A:0.006,G:0.994", "A"), "0.006");
    }

    #[test]
    fn test_extract_af_deletion_allele() {
        assert_eq!(extract_af_for_allele("-:0.001,T:0.999", "-"), "0.001");
    }

    #[test]
    fn test_extract_af_no_match() {
        assert_eq!(extract_af_for_allele("A:0.5", "G"), "");
    }

    #[test]
    fn test_extract_af_empty_input() {
        assert_eq!(extract_af_for_allele("", "T"), "");
    }

    #[test]
    fn test_extract_af_zero_freq() {
        assert_eq!(extract_af_for_allele("T:0", "T"), "0");
    }

    // --- compute_max_af tests ---

    #[test]
    fn test_compute_max_af_single_pop() {
        let entries = vec![("gnomADg_NFE", "0.05")];
        let (max_af, max_pops) = compute_max_af(&entries);
        assert_eq!(max_af, "0.05");
        assert_eq!(max_pops, "gnomADg_NFE");
    }

    #[test]
    fn test_compute_max_af_multiple_pops() {
        let entries = vec![
            ("AFR_AF", "0.01"),
            ("EUR_AF", "0.05"),
            ("gnomADg_NFE", "0.03"),
        ];
        let (max_af, max_pops) = compute_max_af(&entries);
        assert_eq!(max_af, "0.05");
        assert_eq!(max_pops, "EUR_AF");
    }

    #[test]
    fn test_compute_max_af_tie() {
        let entries = vec![("AFR_AF", "0.05"), ("EUR_AF", "0.05")];
        let (max_af, max_pops) = compute_max_af(&entries);
        assert_eq!(max_af, "0.05");
        assert_eq!(max_pops, "AFR_AF&EUR_AF");
    }

    #[test]
    fn test_compute_max_af_empty() {
        let entries: Vec<(&str, &str)> = vec![];
        let (max_af, max_pops) = compute_max_af(&entries);
        assert_eq!(max_af, "");
        assert_eq!(max_pops, "");
    }

    #[test]
    fn test_compute_max_af_skips_empty_freq() {
        let entries = vec![("AF", ""), ("gnomADg_NFE", "0.02")];
        let (max_af, max_pops) = compute_max_af(&entries);
        assert_eq!(max_af, "0.02");
        assert_eq!(max_pops, "gnomADg_NFE");
    }

    // --- VepFlags tests ---

    #[test]
    fn test_vep_flags_none() {
        let flags = VepFlags::from_options_json(None);
        assert!(!flags.check_existing);
        assert!(!flags.af);
    }

    #[test]
    fn test_vep_flags_af_implies_check_existing() {
        let flags = VepFlags::from_options_json(Some(r#"{"af":true}"#));
        assert!(flags.af);
        assert!(flags.check_existing);
    }

    #[test]
    fn test_vep_flags_pubmed_implies_check_existing() {
        let flags = VepFlags::from_options_json(Some(r#"{"pubmed":true}"#));
        assert!(flags.pubmed);
        assert!(flags.check_existing);
    }

    #[test]
    fn test_vep_flags_check_existing_standalone() {
        let flags = VepFlags::from_options_json(Some(r#"{"check_existing":true}"#));
        assert!(flags.check_existing);
        assert!(!flags.af);
        assert!(!flags.pubmed);
    }

    #[test]
    fn test_vep_flags_af_group_enabled() {
        let flags = VepFlags::from_options_json(Some(r#"{"af_gnomadg":true}"#));
        assert!(flags.af_group_enabled(3)); // gnomADg
        assert!(!flags.af_group_enabled(0)); // global AF
        assert!(!flags.af_group_enabled(2)); // gnomADe
    }
}
