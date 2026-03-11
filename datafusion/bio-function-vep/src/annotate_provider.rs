//! Provider for `annotate_vep()` table function.
//!
//! Runtime behavior:
//! - always starts from `lookup_variants()` for known-variant metadata,
//! - when transcript/exon tables are available, computes transcript-driven
//!   consequence terms and most-severe ranking,
//! - otherwise falls back to phase-1.5 known-variant CSQ placeholders.

use std::any::Any;
use std::collections::HashMap;
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
    /// Apply `sprintf("%.4f")` formatting (VEP does this for the global AF field only).
    format_4f: bool,
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
    // --af (global 1000 Genomes) — emitted in CSQ, excluded from MAX_AF_POPS, formatted %.4f
    AfColumn { cache_col: "AF", format_4f: true, flag_group: 0, emit_in_csq: true, max_af_pop: None },
    // --af_1kg (continental) — emitted, MAX_AF uses short names (AFR not AFR_AF)
    AfColumn { cache_col: "AFR", format_4f: false, flag_group: 1, emit_in_csq: true, max_af_pop: Some("AFR") },
    AfColumn { cache_col: "AMR", format_4f: false, flag_group: 1, emit_in_csq: true, max_af_pop: Some("AMR") },
    AfColumn { cache_col: "EAS", format_4f: false, flag_group: 1, emit_in_csq: true, max_af_pop: Some("EAS") },
    AfColumn { cache_col: "EUR", format_4f: false, flag_group: 1, emit_in_csq: true, max_af_pop: Some("EUR") },
    AfColumn { cache_col: "SAS", format_4f: false, flag_group: 1, emit_in_csq: true, max_af_pop: Some("SAS") },
    // --af_gnomade — only global emitted in CSQ; sub-pops used for MAX_AF only
    AfColumn { cache_col: "gnomADe", format_4f: false, flag_group: 2, emit_in_csq: true, max_af_pop: None },
    AfColumn { cache_col: "gnomADe_AFR", format_4f: false, flag_group: 2, emit_in_csq: false, max_af_pop: Some("gnomADe_AFR") },
    AfColumn { cache_col: "gnomADe_AMR", format_4f: false, flag_group: 2, emit_in_csq: false, max_af_pop: Some("gnomADe_AMR") },
    AfColumn { cache_col: "gnomADe_ASJ", format_4f: false, flag_group: 2, emit_in_csq: false, max_af_pop: Some("gnomADe_ASJ") },
    AfColumn { cache_col: "gnomADe_EAS", format_4f: false, flag_group: 2, emit_in_csq: false, max_af_pop: Some("gnomADe_EAS") },
    AfColumn { cache_col: "gnomADe_FIN", format_4f: false, flag_group: 2, emit_in_csq: false, max_af_pop: Some("gnomADe_FIN") },
    AfColumn { cache_col: "gnomADe_MID", format_4f: false, flag_group: 2, emit_in_csq: false, max_af_pop: Some("gnomADe_MID") },
    AfColumn { cache_col: "gnomADe_NFE", format_4f: false, flag_group: 2, emit_in_csq: false, max_af_pop: Some("gnomADe_NFE") },
    AfColumn { cache_col: "gnomADe_REMAINING", format_4f: false, flag_group: 2, emit_in_csq: false, max_af_pop: Some("gnomADe_REMAINING") },
    AfColumn { cache_col: "gnomADe_SAS", format_4f: false, flag_group: 2, emit_in_csq: false, max_af_pop: Some("gnomADe_SAS") },
    // --af_gnomadg — only global emitted in CSQ; sub-pops used for MAX_AF only
    AfColumn { cache_col: "gnomADg", format_4f: false, flag_group: 3, emit_in_csq: true, max_af_pop: None },
    AfColumn { cache_col: "gnomADg_AFR", format_4f: false, flag_group: 3, emit_in_csq: false, max_af_pop: Some("gnomADg_AFR") },
    AfColumn { cache_col: "gnomADg_AMI", format_4f: false, flag_group: 3, emit_in_csq: false, max_af_pop: Some("gnomADg_AMI") },
    AfColumn { cache_col: "gnomADg_AMR", format_4f: false, flag_group: 3, emit_in_csq: false, max_af_pop: Some("gnomADg_AMR") },
    AfColumn { cache_col: "gnomADg_ASJ", format_4f: false, flag_group: 3, emit_in_csq: false, max_af_pop: Some("gnomADg_ASJ") },
    AfColumn { cache_col: "gnomADg_EAS", format_4f: false, flag_group: 3, emit_in_csq: false, max_af_pop: Some("gnomADg_EAS") },
    AfColumn { cache_col: "gnomADg_FIN", format_4f: false, flag_group: 3, emit_in_csq: false, max_af_pop: Some("gnomADg_FIN") },
    AfColumn { cache_col: "gnomADg_MID", format_4f: false, flag_group: 3, emit_in_csq: false, max_af_pop: Some("gnomADg_MID") },
    AfColumn { cache_col: "gnomADg_NFE", format_4f: false, flag_group: 3, emit_in_csq: false, max_af_pop: Some("gnomADg_NFE") },
    AfColumn { cache_col: "gnomADg_REMAINING", format_4f: false, flag_group: 3, emit_in_csq: false, max_af_pop: Some("gnomADg_REMAINING") },
    AfColumn { cache_col: "gnomADg_SAS", format_4f: false, flag_group: 3, emit_in_csq: false, max_af_pop: Some("gnomADg_SAS") },
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

/// A single co-located variant entry with allele and clinical metadata.
#[derive(Debug, Clone)]
struct ColocatedEntry {
    variation_name: String,
    allele_string: String,
    somatic: i64,
    pheno: i64,
    clin_sig_allele: Option<String>,
    pubmed: Option<String>,
}

/// Aggregated data from ALL co-located variants at the same position.
///
/// VEP collects every variant at the same chromosomal position (regardless of
/// allele matching) and aggregates their metadata into `&`-separated strings.
/// This struct holds the pre-aggregated values for a given (chrom, start, end).
#[derive(Debug, Default, Clone)]
struct ColocatedData {
    /// All variation names joined with `&`, ordered non-somatic first.
    existing_variation: String,
    /// Raw entries for allele-specific filtering of SOMATIC, PHENO, CLIN_SIG, PUBMED.
    entries: Vec<ColocatedEntry>,
}

impl ColocatedData {
    /// Build allele-filtered SOMATIC flags.
    ///
    /// VEP includes only allele-matching variants (plus COSMIC_MUTATION/HGMD_MUTATION)
    /// for SOMATIC/PHENO, while Existing_variation includes ALL variants at the position.
    fn somatic_flags_for_allele(&self, vcf_ref: &str, vcf_alt: &str) -> String {
        let filtered = self.allele_filtered_entries(vcf_ref, vcf_alt);
        let strs: Vec<&str> = filtered
            .iter()
            .map(|e| if e.somatic != 0 { "1" } else { "0" })
            .collect();
        if strs.iter().all(|s| *s == "0") {
            String::new()
        } else {
            strs.join("&")
        }
    }

    /// Build allele-filtered PHENO flags.
    fn pheno_flags_for_allele(&self, vcf_ref: &str, vcf_alt: &str) -> String {
        let filtered = self.allele_filtered_entries(vcf_ref, vcf_alt);
        let strs: Vec<&str> = filtered
            .iter()
            .map(|e| if e.pheno != 0 { "1" } else { "0" })
            .collect();
        if strs.iter().all(|s| *s == "0") {
            String::new()
        } else {
            strs.join("&")
        }
    }

    /// Build allele-filtered CLIN_SIG from clin_sig_allele column.
    ///
    /// VEP extracts per-allele clinical significance from clin_sig_allele entries
    /// like `"C:benign;C:benign/likely_benign"` — filtering for the VEP allele and
    /// preserving compound terms like `benign/likely_benign`.
    fn clin_sig_for_allele(&self, vcf_ref: &str, vcf_alt: &str) -> String {
        use crate::allele::vcf_to_vep_allele;
        let (_, vep_allele) = vcf_to_vep_allele(vcf_ref, vcf_alt);
        let filtered = self.allele_filtered_entries(vcf_ref, vcf_alt);
        let mut terms: Vec<String> = Vec::new();
        for entry in &filtered {
            if let Some(csa) = &entry.clin_sig_allele {
                for pair in csa.split(';') {
                    let pair = pair.trim();
                    if let Some(colon) = pair.find(':') {
                        let allele_part = &pair[..colon];
                        let sig_part = &pair[colon + 1..];
                        if allele_part == vep_allele && !sig_part.is_empty() {
                            // Preserve compound terms (e.g. "benign/likely_benign") as-is.
                            if !terms.contains(&sig_part.to_string()) {
                                terms.push(sig_part.to_string());
                            }
                        }
                    }
                }
            }
        }
        terms.join("&")
    }

    /// Build allele-filtered PUBMED by aggregating across all allele-matched variants.
    fn pubmed_for_allele(&self, vcf_ref: &str, vcf_alt: &str) -> String {
        let filtered = self.allele_filtered_entries(vcf_ref, vcf_alt);
        let mut ids: Vec<String> = Vec::new();
        for entry in &filtered {
            if let Some(pm) = &entry.pubmed {
                for id in pm.split(',') {
                    let id = id.trim();
                    if !id.is_empty() && !ids.iter().any(|x| x == id) {
                        ids.push(id.to_string());
                    }
                }
            }
        }
        ids.join("&")
    }

    /// Filter entries by allele compatibility.
    ///
    /// Always includes COSMIC_MUTATION and HGMD_MUTATION entries (no allele check).
    /// For others, uses `allele_matches` to check VCF REF/ALT against cache allele_string.
    ///
    /// Returns entries sorted to match VEP's co-located ordering for SOMATIC/PHENO:
    /// non-somatic first, then somatic, with rs* before others within each group.
    fn allele_filtered_entries(&self, vcf_ref: &str, vcf_alt: &str) -> Vec<&ColocatedEntry> {
        use crate::allele::allele_matches;
        let mut filtered: Vec<&ColocatedEntry> = self
            .entries
            .iter()
            .filter(|e| {
                // COSMIC and HGMD mutations are always included (no allele check).
                if e.allele_string.starts_with("COSMIC_MUTATION")
                    || e.allele_string.starts_with("HGMD_MUTATION")
                {
                    return true;
                }
                allele_matches(vcf_ref, vcf_alt, &e.allele_string)
            })
            .collect();
        // Sort: non-somatic (0) before somatic (1), then rs* before others, then alphabetical.
        filtered.sort_by(|a, b| {
            let a_som = if a.somatic != 0 { 1 } else { 0 };
            let b_som = if b.somatic != 0 { 1 } else { 0 };
            a_som.cmp(&b_som).then_with(|| {
                let a_rs = a.variation_name.starts_with("rs");
                let b_rs = b.variation_name.starts_with("rs");
                b_rs.cmp(&a_rs)
                    .then_with(|| a.variation_name.cmp(&b.variation_name))
            })
        });
        filtered
    }
}

/// Build co-located variant aggregation from cache entries.
///
/// Groups cache rows by (chrom, start, end) and produces one `ColocatedData`
/// per position. Variants are ordered: rs* IDs first, then others (HGMD, COSMIC),
/// alphabetical within each group.
fn build_colocated_map(
    coloc_batches: &[RecordBatch],
) -> HashMap<(String, i64, i64), ColocatedData> {
    // Collect raw entries keyed by position.
    let mut raw: HashMap<(String, i64, i64), Vec<ColocatedEntry>> = HashMap::new();

    for batch in coloc_batches {
        let schema = batch.schema();
        let chrom_idx = schema.index_of("coloc_chrom").unwrap_or(0);
        let start_idx = schema.index_of("coloc_start").unwrap_or(1);
        let end_idx = schema.index_of("coloc_end").unwrap_or(2);
        let var_name_idx = schema.index_of("coloc_variation_name").unwrap_or(3);
        let somatic_idx = schema.index_of("coloc_somatic").unwrap_or(4);
        let pheno_idx = schema.index_of("coloc_pheno").unwrap_or(5);
        let allele_str_idx = schema.index_of("coloc_allele_string").ok();
        let clin_sig_allele_idx = schema.index_of("coloc_clin_sig_allele").ok();
        let pubmed_col_idx = schema.index_of("coloc_pubmed").ok();

        for row in 0..batch.num_rows() {
            let Some(chrom) = string_at(batch.column(chrom_idx).as_ref(), row) else {
                continue;
            };
            let Some(start) = int64_at(batch.column(start_idx).as_ref(), row) else {
                continue;
            };
            let end = int64_at(batch.column(end_idx).as_ref(), row).unwrap_or(start);
            let var_name = string_at(batch.column(var_name_idx).as_ref(), row)
                .unwrap_or_default();
            let somatic = int64_at(batch.column(somatic_idx).as_ref(), row).unwrap_or(0);
            let pheno = int64_at(batch.column(pheno_idx).as_ref(), row).unwrap_or(0);
            let allele_string = allele_str_idx
                .and_then(|idx| string_at(batch.column(idx).as_ref(), row))
                .unwrap_or_default();
            let clin_sig_allele = clin_sig_allele_idx
                .and_then(|idx| string_at(batch.column(idx).as_ref(), row));
            let pubmed = pubmed_col_idx
                .and_then(|idx| string_at(batch.column(idx).as_ref(), row));

            if var_name.is_empty() {
                continue;
            }
            raw.entry((chrom, start, end)).or_default().push(ColocatedEntry {
                variation_name: var_name,
                allele_string,
                somatic,
                pheno,
                clin_sig_allele,
                pubmed,
            });
        }
    }

    // Aggregate into ColocatedData.
    let mut map = HashMap::with_capacity(raw.len());
    for (key, mut entries) in raw {
        // Sort to match VEP's co-located variant ordering:
        // rs* IDs (dbSNP) first, then others (HGMD, COSMIC), alphabetical within.
        entries.sort_by(|a, b| {
            let a_rs = a.variation_name.starts_with("rs");
            let b_rs = b.variation_name.starts_with("rs");
            b_rs.cmp(&a_rs).then_with(|| a.variation_name.cmp(&b.variation_name))
        });

        let names: Vec<&str> = entries.iter().map(|e| e.variation_name.as_str()).collect();
        let existing_variation = names.join("&");

        map.insert(
            key,
            ColocatedData {
                existing_variation,
                entries,
            },
        );
    }
    map
}

/// Replace commas with `&` in a CSQ field value.
///
/// CSQ entries are comma-separated, so any literal comma in a field value
/// (e.g. TREMBL `"B0QYZ8.91,X5DR28.81"` or PUBMED `"123,456"`) would split
/// the entry and corrupt field alignment.  VEP uses `&` as the multi-value
/// separator within CSQ fields.
#[inline]
fn csq_escape(val: &str) -> std::borrow::Cow<'_, str> {
    if val.contains(',') {
        std::borrow::Cow::Owned(val.replace(',', "&"))
    } else {
        std::borrow::Cow::Borrowed(val)
    }
}

/// Format an AF value with Perl's `sprintf("%.4f", $freq)` — 4 decimal places.
/// VEP applies this only to the global `AF` field for backward compatibility.
fn format_af_4f(raw: &str) -> String {
    if raw.is_empty() {
        return String::new();
    }
    let Ok(val) = raw.parse::<f64>() else {
        return raw.to_string();
    };
    format!("{val:.4}")
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

    /// Query the variation cache for ALL variants at positions found in base_batches.
    ///
    /// Aggregates co-located variant data (Existing_variation, SOMATIC, PHENO)
    /// by (chrom, start, end), including non-allele-matching entries like COSMIC.
    async fn build_colocated_variant_map(
        &self,
        base_batches: &[RecordBatch],
        cache_table: &str,
    ) -> Result<HashMap<(String, i64, i64), ColocatedData>> {
        use crate::allele::{vep_norm_end, vep_norm_start};

        // Extract distinct (chrom, vep_start, vep_end) from base_batches.
        // Use VEP-normalized coordinates so they match the cache's coordinate convention.
        let mut positions: HashSet<(String, i64, i64)> = HashSet::new();
        for batch in base_batches {
            let chrom_idx = batch.schema().index_of("chrom").ok();
            let start_idx = batch.schema().index_of("start").ok();
            let ref_idx = batch.schema().index_of("ref").ok();
            let alt_idx = batch.schema().index_of("alt").ok();
            let (Some(ci), Some(si)) = (chrom_idx, start_idx) else {
                continue;
            };
            for row in 0..batch.num_rows() {
                if let (Some(chrom), Some(start)) = (
                    string_at(batch.column(ci).as_ref(), row),
                    int64_at(batch.column(si).as_ref(), row),
                ) {
                    let ref_allele = ref_idx
                        .and_then(|i| string_at(batch.column(i).as_ref(), row));
                    let alt_allele = alt_idx
                        .and_then(|i| string_at(batch.column(i).as_ref(), row));
                    let (vep_start, vep_end) =
                        if let (Some(ref_a), Some(alt_a)) = (&ref_allele, &alt_allele) {
                            (
                                vep_norm_start(start, ref_a, alt_a),
                                vep_norm_end(start, ref_a, alt_a),
                            )
                        } else {
                            (start, start)
                        };
                    let chrom_norm = chrom
                        .strip_prefix("chr")
                        .unwrap_or(&chrom)
                        .to_string();
                    positions.insert((chrom_norm, vep_start, vep_end));
                }
            }
        }
        if positions.is_empty() {
            return Ok(HashMap::new());
        }

        // Register a temporary table of positions for an efficient semi-join.
        let chroms: Vec<&str> = positions.iter().map(|(c, _, _)| c.as_str()).collect();
        let starts: Vec<i64> = positions.iter().map(|(_, s, _)| *s).collect();
        let ends: Vec<i64> = positions.iter().map(|(_, _, e)| *e).collect();
        let pos_batch = RecordBatch::try_new(
            Arc::new(Schema::new(vec![
                Field::new("p_chrom", DataType::Utf8, false),
                Field::new("p_start", DataType::Int64, false),
                Field::new("p_end", DataType::Int64, false),
            ])),
            vec![
                Arc::new(StringArray::from(chroms)),
                Arc::new(Int64Array::from(starts)),
                Arc::new(Int64Array::from(ends)),
            ],
        )?;
        let pos_table = MemTable::try_new(pos_batch.schema(), vec![vec![pos_batch]])?;
        self.session
            .register_table("__vep_coloc_positions", Arc::new(pos_table))?;

        let coloc_sql = format!(
            "SELECT c.`chrom` AS coloc_chrom, \
                    CAST(c.`start` AS BIGINT) AS coloc_start, \
                    CAST(c.`end` AS BIGINT) AS coloc_end, \
                    c.`variation_name` AS coloc_variation_name, \
                    CAST(COALESCE(c.`somatic`, 0) AS BIGINT) AS coloc_somatic, \
                    CAST(COALESCE(c.`phenotype_or_disease`, 0) AS BIGINT) AS coloc_pheno, \
                    c.`allele_string` AS coloc_allele_string, \
                    CAST(c.`clin_sig_allele` AS VARCHAR) AS coloc_clin_sig_allele, \
                    CAST(c.`pubmed` AS VARCHAR) AS coloc_pubmed \
             FROM `{cache_table}` AS c \
             INNER JOIN __vep_coloc_positions AS p \
             ON c.`chrom` = p.`p_chrom` \
             AND CAST(c.`start` AS BIGINT) = p.`p_start` \
             AND CAST(c.`end` AS BIGINT) = p.`p_end` \
             WHERE c.`variation_name` IS NOT NULL"
        );
        let coloc_batches = self.session.sql(&coloc_sql).await?.collect().await?;

        // Clean up temp table (best-effort).
        let _ = self.session.deregister_table("__vep_coloc_positions");

        Ok(build_colocated_map(&coloc_batches))
    }

    async fn scan_with_transcript_engine(
        &self,
        state: &dyn Session,
        projection: Option<&Vec<usize>>,
        lookup_sql: &str,
        requested_columns: &[&str],
        cache_table: &str,
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

        // Build co-located variant aggregation.
        // VEP includes ALL variants at the same position (regardless of allele)
        // for Existing_variation, SOMATIC, and PHENO fields.
        let flags = VepFlags::from_options_json(self.options_json.as_deref());
        let colocated_map = if flags.check_existing {
            self.build_colocated_variant_map(&base_batches, cache_table)
                .await
                .unwrap_or_default()
        } else {
            HashMap::new()
        };

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
            annotated_batches.push(self.annotate_batch_with_transcript_engine(
                batch,
                &engine,
                &ctx,
                &colocated_map,
            )?);
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
        colocated_map: &HashMap<(String, i64, i64), ColocatedData>,
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
            // Look up co-located variant aggregation (all variants at same position).
            // Use VEP-normalized coordinates to match the cache's coordinate convention.
            let start_val = int64_at(batch.column(start_idx).as_ref(), row).unwrap_or(0);
            let chrom_norm = chrom.strip_prefix("chr").unwrap_or(&chrom);
            let (vep_start, vep_end) = {
                use crate::allele::{vep_norm_end, vep_norm_start};
                (
                    vep_norm_start(start_val, &ref_al, &alt_allele),
                    vep_norm_end(start_val, &ref_al, &alt_allele),
                )
            };
            let coloc = colocated_map.get(&(chrom_norm.to_string(), vep_start, vep_end));

            let coloc_existing = coloc.map(|c| c.existing_variation.as_str());
            let existing_var = if flags.check_existing {
                coloc_existing.unwrap_or(&variation_name)
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
            let mut af_raw_formatted: Vec<String> = Vec::with_capacity(AF_COLUMNS.len());
            for (i, col) in AF_COLUMNS.iter().enumerate() {
                let raw_freq = if flags.af_group_enabled(col.flag_group) {
                    extract_af_for_allele(&af_raw[i], &vep_allele)
                } else {
                    ""
                };
                // VEP applies sprintf("%.4f") to the global AF field only.
                let freq = if col.format_4f && !raw_freq.is_empty() {
                    format_af_4f(raw_freq)
                } else {
                    raw_freq.to_string()
                };
                // CSQ output: only emit for columns VEP emits in offline/cache mode.
                let csq_val = if col.emit_in_csq { freq.clone() } else { String::new() };
                af_csq_values.push(csq_val);
                af_raw_formatted.push(freq);
            }
            // MAX_AF: collect sub-population (name, freq) pairs.
            let mut max_af_entries: Vec<(&str, &str)> = Vec::new();
            for (i, col) in AF_COLUMNS.iter().enumerate() {
                if let Some(pop_name) = col.max_af_pop {
                    if !af_raw_formatted[i].is_empty() {
                        max_af_entries.push((pop_name, &af_raw_formatted[i]));
                    }
                }
            }
            let (max_af, max_af_pops) = if flags.max_af {
                compute_max_af(&max_af_entries)
            } else {
                (String::new(), String::new())
            };

            // Clinical / phenotype fields — allele-filtered from co-located data.
            // SOMATIC/PHENO: filter by allele compatibility (always include COSMIC/HGMD).
            // CLIN_SIG: extract from clin_sig_allele column, preserving compound terms.
            // PUBMED: aggregate across all allele-matched co-located variants.
            let allele_filtered_somatic: String;
            let allele_filtered_pheno: String;
            let allele_filtered_clin_sig: String;
            let allele_filtered_pubmed: String;

            let somatic_val: &str = if flags.check_existing {
                if let Some(c) = coloc {
                    allele_filtered_somatic = c.somatic_flags_for_allele(&ref_al, &alt_allele);
                    &allele_filtered_somatic
                } else {
                    somatic_idx
                        .and_then(|idx| int64_at(batch.column(idx).as_ref(), row))
                        .map(|v| if v != 0 { "1" } else { "" })
                        .unwrap_or("")
                }
            } else {
                ""
            };
            let pheno_val: &str = if flags.check_existing {
                if let Some(c) = coloc {
                    allele_filtered_pheno = c.pheno_flags_for_allele(&ref_al, &alt_allele);
                    &allele_filtered_pheno
                } else {
                    pheno_idx
                        .and_then(|idx| int64_at(batch.column(idx).as_ref(), row))
                        .map(|v| if v != 0 { "1" } else { "" })
                        .unwrap_or("")
                }
            } else {
                ""
            };
            let clin_sig: String = if flags.check_existing {
                if let Some(c) = coloc {
                    allele_filtered_clin_sig = c.clin_sig_for_allele(&ref_al, &alt_allele);
                    allele_filtered_clin_sig
                } else {
                    clin_sig_idx
                        .and_then(|idx| string_at(batch.column(idx).as_ref(), row))
                        .map(|v| csq_escape(&v).into_owned())
                        .unwrap_or_default()
                }
            } else {
                String::new()
            };
            let pubmed_val = if flags.pubmed {
                if let Some(c) = coloc {
                    allele_filtered_pubmed = c.pubmed_for_allele(&ref_al, &alt_allele);
                    allele_filtered_pubmed
                } else {
                    pubmed_idx
                        .and_then(|idx| string_at(batch.column(idx).as_ref(), row))
                        .map(|v| csq_escape(&v).into_owned())
                        .unwrap_or_default()
                }
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
                    "{vep_allele}|{csq_val}|{impact}|||||||||||||||{existing_var}||||||||||||\
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
                    let trembl_raw = tx_opt
                        .and_then(|tx| tx.trembl.as_deref())
                        .unwrap_or("");
                    let trembl = csq_escape(trembl_raw);
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
                        "{vep_allele}|sequence_variant|{impact}|||||||||||||||{existing_var}||||||||||||\
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
                    &cache_table,
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

    // =======================================================================
    // flags_str_from_bools
    // =======================================================================

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

    // =======================================================================
    // csq_escape
    // =======================================================================

    #[test]
    fn test_csq_escape_no_comma() {
        let result = csq_escape("benign");
        assert_eq!(result.as_ref(), "benign");
        assert!(matches!(result, std::borrow::Cow::Borrowed(_)));
    }

    #[test]
    fn test_csq_escape_with_comma() {
        let result = csq_escape("benign,likely_benign");
        assert_eq!(result.as_ref(), "benign&likely_benign");
        assert!(matches!(result, std::borrow::Cow::Owned(_)));
    }

    #[test]
    fn test_csq_escape_multiple_commas() {
        assert_eq!(csq_escape("a,b,c").as_ref(), "a&b&c");
    }

    #[test]
    fn test_csq_escape_empty() {
        let result = csq_escape("");
        assert_eq!(result.as_ref(), "");
        assert!(matches!(result, std::borrow::Cow::Borrowed(_)));
    }

    // =======================================================================
    // extract_af_for_allele
    // =======================================================================

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

    #[test]
    fn test_extract_af_no_colon_in_entry() {
        // Malformed entry without colon — should be skipped, return empty.
        assert_eq!(extract_af_for_allele("T0.5", "T"), "");
    }

    #[test]
    fn test_extract_af_insertion_allele() {
        assert_eq!(extract_af_for_allele("AGT:0.01,-:0.99", "AGT"), "0.01");
    }

    // =======================================================================
    // format_af_4f
    // =======================================================================

    #[test]
    fn test_format_af_4f_rounds_to_4_decimals() {
        assert_eq!(format_af_4f("0.03015"), "0.0301");
        assert_eq!(format_af_4f("0.007987"), "0.0080");
        assert_eq!(format_af_4f("0.000599"), "0.0006");
    }

    #[test]
    fn test_format_af_4f_preserves_short_values() {
        assert_eq!(format_af_4f("0.5"), "0.5000");
        assert_eq!(format_af_4f("0.1769"), "0.1769");
        assert_eq!(format_af_4f("0"), "0.0000");
    }

    #[test]
    fn test_format_af_4f_empty() {
        assert_eq!(format_af_4f(""), "");
    }

    #[test]
    fn test_format_af_4f_unparseable() {
        // Non-numeric input — returned as-is.
        assert_eq!(format_af_4f("not_a_number"), "not_a_number");
    }

    #[test]
    fn test_format_af_4f_one() {
        assert_eq!(format_af_4f("1"), "1.0000");
    }

    // =======================================================================
    // compute_max_af
    // =======================================================================

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

    #[test]
    fn test_compute_max_af_skips_unparseable_freq() {
        let entries = vec![("AF", "bad"), ("gnomADg_NFE", "0.03")];
        let (max_af, max_pops) = compute_max_af(&entries);
        assert_eq!(max_af, "0.03");
        assert_eq!(max_pops, "gnomADg_NFE");
    }

    #[test]
    fn test_compute_max_af_all_empty_freqs() {
        let entries = vec![("AF", ""), ("gnomADg_NFE", "")];
        let (max_af, max_pops) = compute_max_af(&entries);
        assert_eq!(max_af, "");
        assert_eq!(max_pops, "");
    }

    // =======================================================================
    // classify_variant
    // =======================================================================

    #[test]
    fn test_classify_variant_snv() {
        assert_eq!(classify_variant("A", "G"), "SNV");
        assert_eq!(classify_variant("C", "T"), "SNV");
    }

    #[test]
    fn test_classify_variant_insertion() {
        assert_eq!(classify_variant("-", "T"), "insertion");
        assert_eq!(classify_variant("-", "AGT"), "insertion");
    }

    #[test]
    fn test_classify_variant_deletion() {
        assert_eq!(classify_variant("A", "-"), "deletion");
        assert_eq!(classify_variant("AGT", "-"), "deletion");
    }

    #[test]
    fn test_classify_variant_substitution() {
        assert_eq!(classify_variant("AT", "GC"), "substitution");
        assert_eq!(classify_variant("ATG", "GCA"), "substitution");
    }

    #[test]
    fn test_classify_variant_indel() {
        assert_eq!(classify_variant("AT", "G"), "indel");
        assert_eq!(classify_variant("A", "GC"), "indel");
    }

    // =======================================================================
    // impact_label
    // =======================================================================

    #[test]
    fn test_impact_label_all_variants() {
        assert_eq!(impact_label(SoImpact::High), "HIGH");
        assert_eq!(impact_label(SoImpact::Moderate), "MODERATE");
        assert_eq!(impact_label(SoImpact::Low), "LOW");
        assert_eq!(impact_label(SoImpact::Modifier), "MODIFIER");
    }

    // =======================================================================
    // parse_sv_feature_kind
    // =======================================================================

    #[test]
    fn test_parse_sv_feature_kind_known() {
        assert_eq!(parse_sv_feature_kind("transcript"), Some(SvFeatureKind::Transcript));
        assert_eq!(parse_sv_feature_kind("tx"), Some(SvFeatureKind::Transcript));
        assert_eq!(parse_sv_feature_kind("regulatory"), Some(SvFeatureKind::Regulatory));
        assert_eq!(parse_sv_feature_kind("reg"), Some(SvFeatureKind::Regulatory));
        assert_eq!(parse_sv_feature_kind("tfbs"), Some(SvFeatureKind::Tfbs));
        assert_eq!(parse_sv_feature_kind("motif"), Some(SvFeatureKind::Tfbs));
        assert_eq!(parse_sv_feature_kind("feature"), Some(SvFeatureKind::Generic));
        assert_eq!(parse_sv_feature_kind("generic"), Some(SvFeatureKind::Generic));
    }

    #[test]
    fn test_parse_sv_feature_kind_case_insensitive() {
        assert_eq!(parse_sv_feature_kind("TRANSCRIPT"), Some(SvFeatureKind::Transcript));
        assert_eq!(parse_sv_feature_kind("Regulatory"), Some(SvFeatureKind::Regulatory));
    }

    #[test]
    fn test_parse_sv_feature_kind_unknown() {
        assert_eq!(parse_sv_feature_kind("unknown"), None);
        assert_eq!(parse_sv_feature_kind(""), None);
    }

    // =======================================================================
    // parse_sv_event_kind
    // =======================================================================

    #[test]
    fn test_parse_sv_event_kind_known() {
        assert_eq!(parse_sv_event_kind("ablation"), Some(SvEventKind::Ablation));
        assert_eq!(parse_sv_event_kind("deletion"), Some(SvEventKind::Ablation));
        assert_eq!(parse_sv_event_kind("del"), Some(SvEventKind::Ablation));
        assert_eq!(parse_sv_event_kind("amplification"), Some(SvEventKind::Amplification));
        assert_eq!(parse_sv_event_kind("duplication"), Some(SvEventKind::Amplification));
        assert_eq!(parse_sv_event_kind("dup"), Some(SvEventKind::Amplification));
        assert_eq!(parse_sv_event_kind("amp"), Some(SvEventKind::Amplification));
        assert_eq!(parse_sv_event_kind("elongation"), Some(SvEventKind::Elongation));
        assert_eq!(parse_sv_event_kind("elongate"), Some(SvEventKind::Elongation));
        assert_eq!(parse_sv_event_kind("truncation"), Some(SvEventKind::Truncation));
        assert_eq!(parse_sv_event_kind("truncate"), Some(SvEventKind::Truncation));
    }

    #[test]
    fn test_parse_sv_event_kind_case_insensitive() {
        assert_eq!(parse_sv_event_kind("ABLATION"), Some(SvEventKind::Ablation));
        assert_eq!(parse_sv_event_kind("Truncation"), Some(SvEventKind::Truncation));
    }

    #[test]
    fn test_parse_sv_event_kind_unknown() {
        assert_eq!(parse_sv_event_kind("unknown"), None);
        assert_eq!(parse_sv_event_kind(""), None);
    }

    // =======================================================================
    // VepFlags
    // =======================================================================

    #[test]
    fn test_vep_flags_none() {
        let flags = VepFlags::from_options_json(None);
        assert!(!flags.check_existing);
        assert!(!flags.af);
        assert!(!flags.af_1kg);
        assert!(!flags.af_gnomade);
        assert!(!flags.af_gnomadg);
        assert!(!flags.max_af);
        assert!(!flags.pubmed);
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

    #[test]
    fn test_vep_flags_af_group_enabled_out_of_range() {
        let flags = VepFlags::from_options_json(Some(r#"{"af":true}"#));
        assert!(!flags.af_group_enabled(99));
    }

    #[test]
    fn test_vep_flags_all_af_flags() {
        let flags = VepFlags::from_options_json(Some(
            r#"{"af":true,"af_1kg":true,"af_gnomade":true,"af_gnomadg":true,"max_af":true}"#,
        ));
        assert!(flags.af);
        assert!(flags.af_1kg);
        assert!(flags.af_gnomade);
        assert!(flags.af_gnomadg);
        assert!(flags.max_af);
        assert!(flags.check_existing);
    }

    #[test]
    fn test_vep_flags_explicit_false() {
        let flags = VepFlags::from_options_json(Some(r#"{"af":false}"#));
        assert!(!flags.af);
        assert!(!flags.check_existing);
    }

    #[test]
    fn test_vep_flags_empty_json() {
        let flags = VepFlags::from_options_json(Some("{}"));
        assert!(!flags.check_existing);
        assert!(!flags.af);
    }

    #[test]
    fn test_vep_flags_max_af_implies_check_existing() {
        let flags = VepFlags::from_options_json(Some(r#"{"max_af":true}"#));
        assert!(flags.max_af);
        assert!(flags.check_existing);
    }

    #[test]
    fn test_vep_flags_af_1kg_implies_check_existing() {
        let flags = VepFlags::from_options_json(Some(r#"{"af_1kg":true}"#));
        assert!(flags.af_1kg);
        assert!(flags.check_existing);
    }

    #[test]
    fn test_vep_flags_af_gnomade_implies_check_existing() {
        let flags = VepFlags::from_options_json(Some(r#"{"af_gnomade":true}"#));
        assert!(flags.af_gnomade);
        assert!(flags.check_existing);
    }

    // =======================================================================
    // parse_json_string_option / parse_json_bool_option
    // =======================================================================

    #[test]
    fn test_parse_json_string_option_basic() {
        let json = r#"{"key":"value"}"#;
        assert_eq!(
            AnnotateProvider::parse_json_string_option(json, "key"),
            Some("value".to_string())
        );
    }

    #[test]
    fn test_parse_json_string_option_missing_key() {
        let json = r#"{"other":"value"}"#;
        assert_eq!(AnnotateProvider::parse_json_string_option(json, "key"), None);
    }

    #[test]
    fn test_parse_json_string_option_empty_value() {
        let json = r#"{"key":""}"#;
        assert_eq!(AnnotateProvider::parse_json_string_option(json, "key"), None);
    }

    #[test]
    fn test_parse_json_string_option_backtick_rejected() {
        let json = r#"{"key":"val`ue"}"#;
        assert_eq!(AnnotateProvider::parse_json_string_option(json, "key"), None);
    }

    #[test]
    fn test_parse_json_string_option_with_spaces() {
        let json = r#"{"key" : "value"}"#;
        assert_eq!(
            AnnotateProvider::parse_json_string_option(json, "key"),
            Some("value".to_string())
        );
    }

    #[test]
    fn test_parse_json_string_option_multiple_keys() {
        let json = r#"{"a":"first","b":"second"}"#;
        assert_eq!(
            AnnotateProvider::parse_json_string_option(json, "b"),
            Some("second".to_string())
        );
    }

    #[test]
    fn test_parse_json_bool_option_true() {
        let json = r#"{"flag":true}"#;
        assert_eq!(AnnotateProvider::parse_json_bool_option(json, "flag"), Some(true));
    }

    #[test]
    fn test_parse_json_bool_option_false() {
        let json = r#"{"flag":false}"#;
        assert_eq!(
            AnnotateProvider::parse_json_bool_option(json, "flag"),
            Some(false)
        );
    }

    #[test]
    fn test_parse_json_bool_option_missing_key() {
        let json = r#"{"other":true}"#;
        assert_eq!(AnnotateProvider::parse_json_bool_option(json, "flag"), None);
    }

    #[test]
    fn test_parse_json_bool_option_non_bool_value() {
        let json = r#"{"flag":"yes"}"#;
        assert_eq!(AnnotateProvider::parse_json_bool_option(json, "flag"), None);
    }

    #[test]
    fn test_parse_json_bool_option_with_spaces() {
        let json = r#"{"flag" :  true}"#;
        assert_eq!(AnnotateProvider::parse_json_bool_option(json, "flag"), Some(true));
    }

    // =======================================================================
    // escaped_sql_literal
    // =======================================================================

    #[test]
    fn test_escaped_sql_literal_no_quotes() {
        assert_eq!(AnnotateProvider::escaped_sql_literal("table_name"), "table_name");
    }

    #[test]
    fn test_escaped_sql_literal_single_quote() {
        assert_eq!(AnnotateProvider::escaped_sql_literal("it's"), "it''s");
    }

    #[test]
    fn test_escaped_sql_literal_multiple_quotes() {
        assert_eq!(AnnotateProvider::escaped_sql_literal("a'b'c"), "a''b''c");
    }

    // =======================================================================
    // ColocatedData — allele filtering methods
    // =======================================================================

    fn make_entry(
        name: &str,
        allele_string: &str,
        somatic: i64,
        pheno: i64,
        clin_sig_allele: Option<&str>,
        pubmed: Option<&str>,
    ) -> ColocatedEntry {
        ColocatedEntry {
            variation_name: name.to_string(),
            allele_string: allele_string.to_string(),
            somatic,
            pheno,
            clin_sig_allele: clin_sig_allele.map(|s| s.to_string()),
            pubmed: pubmed.map(|s| s.to_string()),
        }
    }

    fn make_coloc_data(entries: Vec<ColocatedEntry>) -> ColocatedData {
        let names: Vec<&str> = entries.iter().map(|e| e.variation_name.as_str()).collect();
        ColocatedData {
            existing_variation: names.join("&"),
            entries,
        }
    }

    #[test]
    fn test_allele_filtered_entries_exact_match() {
        let data = make_coloc_data(vec![
            make_entry("rs1", "T/C", 0, 0, None, None),
            make_entry("rs2", "T/A", 0, 0, None, None),
        ]);
        // VCF T>C should match T/C but not T/A.
        let filtered = data.allele_filtered_entries("T", "C");
        assert_eq!(filtered.len(), 1);
        assert_eq!(filtered[0].variation_name, "rs1");
    }

    #[test]
    fn test_allele_filtered_entries_cosmic_always_included() {
        let data = make_coloc_data(vec![
            make_entry("rs1", "T/C", 0, 0, None, None),
            make_entry("COSV123", "COSMIC_MUTATION", 1, 0, None, None),
        ]);
        let filtered = data.allele_filtered_entries("T", "A");
        // rs1 (T/C) doesn't match T/A, but COSMIC always included.
        assert_eq!(filtered.len(), 1);
        assert_eq!(filtered[0].variation_name, "COSV123");
    }

    #[test]
    fn test_allele_filtered_entries_hgmd_always_included() {
        let data = make_coloc_data(vec![
            make_entry("HGMD1", "HGMD_MUTATION", 0, 1, None, None),
        ]);
        let filtered = data.allele_filtered_entries("A", "G");
        assert_eq!(filtered.len(), 1);
        assert_eq!(filtered[0].variation_name, "HGMD1");
    }

    #[test]
    fn test_allele_filtered_entries_sort_non_somatic_first() {
        let data = make_coloc_data(vec![
            make_entry("COSV1", "COSMIC_MUTATION", 1, 0, None, None),
            make_entry("rs1", "T/C", 0, 0, None, None),
        ]);
        let filtered = data.allele_filtered_entries("T", "C");
        assert_eq!(filtered.len(), 2);
        // Non-somatic (rs1) should come before somatic (COSV1).
        assert_eq!(filtered[0].variation_name, "rs1");
        assert_eq!(filtered[1].variation_name, "COSV1");
    }

    #[test]
    fn test_allele_filtered_entries_empty_result() {
        let data = make_coloc_data(vec![
            make_entry("rs1", "T/A", 0, 0, None, None),
        ]);
        let filtered = data.allele_filtered_entries("G", "C");
        assert!(filtered.is_empty());
    }

    #[test]
    fn test_allele_filtered_entries_rs_before_others_within_group() {
        let data = make_coloc_data(vec![
            make_entry("COSV1", "COSMIC_MUTATION", 0, 0, None, None),
            make_entry("rs1", "T/C", 0, 0, None, None),
            make_entry("rs2", "T/C", 0, 0, None, None),
        ]);
        let filtered = data.allele_filtered_entries("T", "C");
        assert_eq!(filtered.len(), 3);
        // All non-somatic: rs before COSV, alphabetical within.
        assert_eq!(filtered[0].variation_name, "rs1");
        assert_eq!(filtered[1].variation_name, "rs2");
        assert_eq!(filtered[2].variation_name, "COSV1");
    }

    // --- somatic_flags_for_allele ---

    #[test]
    fn test_somatic_flags_all_non_somatic() {
        let data = make_coloc_data(vec![
            make_entry("rs1", "T/C", 0, 0, None, None),
            make_entry("rs2", "T/C", 0, 0, None, None),
        ]);
        // All zeros → empty string.
        assert_eq!(data.somatic_flags_for_allele("T", "C"), "");
    }

    #[test]
    fn test_somatic_flags_mixed() {
        let data = make_coloc_data(vec![
            make_entry("rs1", "T/C", 0, 0, None, None),
            make_entry("COSV1", "COSMIC_MUTATION", 1, 0, None, None),
        ]);
        // Non-somatic first: rs1(0), COSV1(1).
        assert_eq!(data.somatic_flags_for_allele("T", "C"), "0&1");
    }

    #[test]
    fn test_somatic_flags_single_somatic() {
        let data = make_coloc_data(vec![
            make_entry("COSV1", "COSMIC_MUTATION", 1, 0, None, None),
        ]);
        assert_eq!(data.somatic_flags_for_allele("A", "G"), "1");
    }

    #[test]
    fn test_somatic_flags_no_matching_entries() {
        let data = make_coloc_data(vec![
            make_entry("rs1", "T/A", 1, 0, None, None),
        ]);
        // No allele match, no COSMIC → empty (no entries).
        assert_eq!(data.somatic_flags_for_allele("G", "C"), "");
    }

    // --- pheno_flags_for_allele ---

    #[test]
    fn test_pheno_flags_all_non_pheno() {
        let data = make_coloc_data(vec![
            make_entry("rs1", "T/C", 0, 0, None, None),
        ]);
        assert_eq!(data.pheno_flags_for_allele("T", "C"), "");
    }

    #[test]
    fn test_pheno_flags_mixed() {
        let data = make_coloc_data(vec![
            make_entry("rs1", "T/C", 0, 0, None, None),
            make_entry("rs2", "T/C", 0, 1, None, None),
        ]);
        assert_eq!(data.pheno_flags_for_allele("T", "C"), "0&1");
    }

    // --- clin_sig_for_allele ---

    #[test]
    fn test_clin_sig_for_allele_simple() {
        let data = make_coloc_data(vec![
            make_entry("rs1", "T/C", 0, 0, Some("C:benign"), None),
        ]);
        assert_eq!(data.clin_sig_for_allele("T", "C"), "benign");
    }

    #[test]
    fn test_clin_sig_for_allele_compound_term() {
        let data = make_coloc_data(vec![
            make_entry("rs1", "T/C", 0, 0, Some("C:benign;C:benign/likely_benign"), None),
        ]);
        assert_eq!(
            data.clin_sig_for_allele("T", "C"),
            "benign&benign/likely_benign"
        );
    }

    #[test]
    fn test_clin_sig_for_allele_wrong_allele_filtered() {
        let data = make_coloc_data(vec![
            make_entry("rs1", "T/C", 0, 0, Some("A:pathogenic"), None),
        ]);
        // VEP allele for T>C is "C", but clin_sig_allele has "A:pathogenic" → no match.
        assert_eq!(data.clin_sig_for_allele("T", "C"), "");
    }

    #[test]
    fn test_clin_sig_for_allele_no_clin_sig() {
        let data = make_coloc_data(vec![
            make_entry("rs1", "T/C", 0, 0, None, None),
        ]);
        assert_eq!(data.clin_sig_for_allele("T", "C"), "");
    }

    #[test]
    fn test_clin_sig_for_allele_deduplicates() {
        let data = make_coloc_data(vec![
            make_entry("rs1", "T/C", 0, 0, Some("C:benign"), None),
            make_entry("rs2", "T/C", 0, 0, Some("C:benign"), None),
        ]);
        assert_eq!(data.clin_sig_for_allele("T", "C"), "benign");
    }

    #[test]
    fn test_clin_sig_for_allele_multiple_entries() {
        let data = make_coloc_data(vec![
            make_entry("rs1", "T/C", 0, 0, Some("C:benign"), None),
            make_entry("rs2", "T/C", 0, 0, Some("C:pathogenic"), None),
        ]);
        assert_eq!(data.clin_sig_for_allele("T", "C"), "benign&pathogenic");
    }

    #[test]
    fn test_clin_sig_for_allele_missing_colon() {
        // Malformed clin_sig_allele without colon → skipped.
        let data = make_coloc_data(vec![
            make_entry("rs1", "T/C", 0, 0, Some("benign_only"), None),
        ]);
        assert_eq!(data.clin_sig_for_allele("T", "C"), "");
    }

    #[test]
    fn test_clin_sig_for_allele_empty_sig_part() {
        // Empty significance after colon → skipped.
        let data = make_coloc_data(vec![
            make_entry("rs1", "T/C", 0, 0, Some("C:"), None),
        ]);
        assert_eq!(data.clin_sig_for_allele("T", "C"), "");
    }

    // --- pubmed_for_allele ---

    #[test]
    fn test_pubmed_for_allele_single() {
        let data = make_coloc_data(vec![
            make_entry("rs1", "T/C", 0, 0, None, Some("12345")),
        ]);
        assert_eq!(data.pubmed_for_allele("T", "C"), "12345");
    }

    #[test]
    fn test_pubmed_for_allele_aggregation() {
        let data = make_coloc_data(vec![
            make_entry("rs1", "T/C", 0, 0, None, Some("111")),
            make_entry("rs2", "T/C", 0, 0, None, Some("222")),
        ]);
        assert_eq!(data.pubmed_for_allele("T", "C"), "111&222");
    }

    #[test]
    fn test_pubmed_for_allele_deduplication() {
        let data = make_coloc_data(vec![
            make_entry("rs1", "T/C", 0, 0, None, Some("111,222")),
            make_entry("rs2", "T/C", 0, 0, None, Some("222,333")),
        ]);
        assert_eq!(data.pubmed_for_allele("T", "C"), "111&222&333");
    }

    #[test]
    fn test_pubmed_for_allele_no_pubmed() {
        let data = make_coloc_data(vec![
            make_entry("rs1", "T/C", 0, 0, None, None),
        ]);
        assert_eq!(data.pubmed_for_allele("T", "C"), "");
    }

    #[test]
    fn test_pubmed_for_allele_non_matching_allele_excluded() {
        let data = make_coloc_data(vec![
            make_entry("rs1", "T/A", 0, 0, None, Some("12345")),
        ]);
        // T/A doesn't match VCF T>C.
        assert_eq!(data.pubmed_for_allele("T", "C"), "");
    }

    #[test]
    fn test_pubmed_for_allele_cosmic_included() {
        let data = make_coloc_data(vec![
            make_entry("COSV1", "COSMIC_MUTATION", 1, 0, None, Some("99999")),
        ]);
        assert_eq!(data.pubmed_for_allele("A", "G"), "99999");
    }

    // =======================================================================
    // build_colocated_map
    // =======================================================================

    fn make_coloc_batch(
        rows: &[(&str, i64, i64, &str, i64, i64, &str, Option<&str>, Option<&str>)],
    ) -> RecordBatch {
        let chroms: Vec<&str> = rows.iter().map(|r| r.0).collect();
        let starts: Vec<i64> = rows.iter().map(|r| r.1).collect();
        let ends: Vec<i64> = rows.iter().map(|r| r.2).collect();
        let names: Vec<&str> = rows.iter().map(|r| r.3).collect();
        let somatics: Vec<i64> = rows.iter().map(|r| r.4).collect();
        let phenos: Vec<i64> = rows.iter().map(|r| r.5).collect();
        let allele_strings: Vec<&str> = rows.iter().map(|r| r.6).collect();
        let clin_sig_alleles: Vec<Option<&str>> = rows.iter().map(|r| r.7).collect();
        let pubmeds: Vec<Option<&str>> = rows.iter().map(|r| r.8).collect();

        RecordBatch::try_new(
            Arc::new(Schema::new(vec![
                Field::new("coloc_chrom", DataType::Utf8, false),
                Field::new("coloc_start", DataType::Int64, false),
                Field::new("coloc_end", DataType::Int64, false),
                Field::new("coloc_variation_name", DataType::Utf8, false),
                Field::new("coloc_somatic", DataType::Int64, false),
                Field::new("coloc_pheno", DataType::Int64, false),
                Field::new("coloc_allele_string", DataType::Utf8, false),
                Field::new("coloc_clin_sig_allele", DataType::Utf8, true),
                Field::new("coloc_pubmed", DataType::Utf8, true),
            ])),
            vec![
                Arc::new(StringArray::from(chroms)),
                Arc::new(Int64Array::from(starts)),
                Arc::new(Int64Array::from(ends)),
                Arc::new(StringArray::from(names)),
                Arc::new(Int64Array::from(somatics)),
                Arc::new(Int64Array::from(phenos)),
                Arc::new(StringArray::from(allele_strings)),
                Arc::new(StringArray::from(clin_sig_alleles)),
                Arc::new(StringArray::from(pubmeds)),
            ],
        )
        .unwrap()
    }

    #[test]
    fn test_build_colocated_map_single_entry() {
        let batch = make_coloc_batch(&[
            ("22", 100, 100, "rs123", 0, 0, "T/C", None, None),
        ]);
        let map = build_colocated_map(&[batch]);
        assert_eq!(map.len(), 1);
        let data = &map[&("22".to_string(), 100, 100)];
        assert_eq!(data.existing_variation, "rs123");
        assert_eq!(data.entries.len(), 1);
    }

    #[test]
    fn test_build_colocated_map_rs_sorted_first() {
        let batch = make_coloc_batch(&[
            ("22", 100, 100, "COSV123", 1, 0, "COSMIC_MUTATION", None, None),
            ("22", 100, 100, "rs456", 0, 0, "T/C", None, None),
            ("22", 100, 100, "rs123", 0, 0, "T/A", None, None),
        ]);
        let map = build_colocated_map(&[batch]);
        let data = &map[&("22".to_string(), 100, 100)];
        // rs* first (alphabetical), then COSV.
        assert_eq!(data.existing_variation, "rs123&rs456&COSV123");
    }

    #[test]
    fn test_build_colocated_map_skips_empty_variation_name() {
        let batch = make_coloc_batch(&[
            ("22", 100, 100, "", 0, 0, "T/C", None, None),
            ("22", 100, 100, "rs1", 0, 0, "T/C", None, None),
        ]);
        let map = build_colocated_map(&[batch]);
        let data = &map[&("22".to_string(), 100, 100)];
        assert_eq!(data.existing_variation, "rs1");
        assert_eq!(data.entries.len(), 1);
    }

    #[test]
    fn test_build_colocated_map_multiple_positions() {
        let batch = make_coloc_batch(&[
            ("22", 100, 100, "rs1", 0, 0, "T/C", None, None),
            ("22", 200, 200, "rs2", 0, 0, "A/G", None, None),
        ]);
        let map = build_colocated_map(&[batch]);
        assert_eq!(map.len(), 2);
        assert!(map.contains_key(&("22".to_string(), 100, 100)));
        assert!(map.contains_key(&("22".to_string(), 200, 200)));
    }

    #[test]
    fn test_build_colocated_map_preserves_clinical_metadata() {
        let batch = make_coloc_batch(&[
            ("22", 100, 100, "rs1", 0, 1, "T/C", Some("C:benign"), Some("12345")),
        ]);
        let map = build_colocated_map(&[batch]);
        let data = &map[&("22".to_string(), 100, 100)];
        let entry = &data.entries[0];
        assert_eq!(entry.clin_sig_allele.as_deref(), Some("C:benign"));
        assert_eq!(entry.pubmed.as_deref(), Some("12345"));
        assert_eq!(entry.pheno, 1);
    }

    #[test]
    fn test_build_colocated_map_empty_batches() {
        let map = build_colocated_map(&[]);
        assert!(map.is_empty());
    }

    #[test]
    fn test_build_colocated_map_multiple_batches() {
        let batch1 = make_coloc_batch(&[
            ("22", 100, 100, "rs1", 0, 0, "T/C", None, None),
        ]);
        let batch2 = make_coloc_batch(&[
            ("22", 100, 100, "rs2", 0, 0, "T/A", None, None),
        ]);
        let map = build_colocated_map(&[batch1, batch2]);
        let data = &map[&("22".to_string(), 100, 100)];
        assert_eq!(data.entries.len(), 2);
        assert_eq!(data.existing_variation, "rs1&rs2");
    }

    // =======================================================================
    // int64_at — type coercion
    // =======================================================================

    #[test]
    fn test_int64_at_int64() {
        let arr = Int64Array::from(vec![42]);
        assert_eq!(int64_at(&arr, 0), Some(42));
    }

    #[test]
    fn test_int64_at_int32() {
        let arr = Int32Array::from(vec![42]);
        assert_eq!(int64_at(&arr, 0), Some(42));
    }

    #[test]
    fn test_int64_at_int16() {
        let arr = Int16Array::from(vec![42]);
        assert_eq!(int64_at(&arr, 0), Some(42));
    }

    #[test]
    fn test_int64_at_int8() {
        let arr = Int8Array::from(vec![42]);
        assert_eq!(int64_at(&arr, 0), Some(42));
    }

    #[test]
    fn test_int64_at_uint64() {
        let arr = UInt64Array::from(vec![42u64]);
        assert_eq!(int64_at(&arr, 0), Some(42));
    }

    #[test]
    fn test_int64_at_uint64_overflow() {
        let arr = UInt64Array::from(vec![u64::MAX]);
        // u64::MAX can't fit in i64 → None.
        assert_eq!(int64_at(&arr, 0), None);
    }

    #[test]
    fn test_int64_at_uint32() {
        let arr = UInt32Array::from(vec![42u32]);
        assert_eq!(int64_at(&arr, 0), Some(42));
    }

    #[test]
    fn test_int64_at_uint16() {
        let arr = UInt16Array::from(vec![42u16]);
        assert_eq!(int64_at(&arr, 0), Some(42));
    }

    #[test]
    fn test_int64_at_uint8() {
        let arr = UInt8Array::from(vec![42u8]);
        assert_eq!(int64_at(&arr, 0), Some(42));
    }

    #[test]
    fn test_int64_at_null() {
        let arr = Int64Array::from(vec![None as Option<i64>]);
        assert_eq!(int64_at(&arr, 0), None);
    }

    #[test]
    fn test_int64_at_unsupported_type() {
        let arr = Float64Array::from(vec![42.0]);
        assert_eq!(int64_at(&arr, 0), None);
    }

    // =======================================================================
    // string_at — type coercion
    // =======================================================================

    #[test]
    fn test_string_at_string_array() {
        let arr = StringArray::from(vec!["hello"]);
        assert_eq!(string_at(&arr, 0), Some("hello".to_string()));
    }

    #[test]
    fn test_string_at_large_string_array() {
        let arr = LargeStringArray::from(vec!["large"]);
        assert_eq!(string_at(&arr, 0), Some("large".to_string()));
    }

    #[test]
    fn test_string_at_string_view_array() {
        let arr = StringViewArray::from(vec!["view"]);
        assert_eq!(string_at(&arr, 0), Some("view".to_string()));
    }

    #[test]
    fn test_string_at_float64() {
        let arr = Float64Array::from(vec![3.14]);
        assert_eq!(string_at(&arr, 0), Some("3.14".to_string()));
    }

    #[test]
    fn test_string_at_float32() {
        let arr = Float32Array::from(vec![2.5f32]);
        assert_eq!(string_at(&arr, 0), Some("2.5".to_string()));
    }

    #[test]
    fn test_string_at_int64() {
        let arr = Int64Array::from(vec![42]);
        assert_eq!(string_at(&arr, 0), Some("42".to_string()));
    }

    #[test]
    fn test_string_at_int32() {
        let arr = Int32Array::from(vec![42]);
        assert_eq!(string_at(&arr, 0), Some("42".to_string()));
    }

    #[test]
    fn test_string_at_uint64() {
        let arr = UInt64Array::from(vec![42u64]);
        assert_eq!(string_at(&arr, 0), Some("42".to_string()));
    }

    #[test]
    fn test_string_at_uint32() {
        let arr = UInt32Array::from(vec![42u32]);
        assert_eq!(string_at(&arr, 0), Some("42".to_string()));
    }

    #[test]
    fn test_string_at_null() {
        let arr = StringArray::from(vec![None as Option<&str>]);
        assert_eq!(string_at(&arr, 0), None);
    }

    #[test]
    fn test_string_at_unsupported_type() {
        let arr = BooleanArray::from(vec![true]);
        assert_eq!(string_at(&arr, 0), None);
    }

    // =======================================================================
    // bool_at
    // =======================================================================

    #[test]
    fn test_bool_at_true() {
        let arr = BooleanArray::from(vec![true]);
        assert_eq!(bool_at(&arr, 0), Some(true));
    }

    #[test]
    fn test_bool_at_false() {
        let arr = BooleanArray::from(vec![false]);
        assert_eq!(bool_at(&arr, 0), Some(false));
    }

    #[test]
    fn test_bool_at_null() {
        let arr = BooleanArray::from(vec![None as Option<bool>]);
        assert_eq!(bool_at(&arr, 0), None);
    }

    #[test]
    fn test_bool_at_unsupported_type() {
        let arr = Int64Array::from(vec![1]);
        assert_eq!(bool_at(&arr, 0), None);
    }

    // =======================================================================
    // read_mirna_regions
    // =======================================================================

    #[test]
    fn test_read_mirna_regions_null_cell() {
        use datafusion::arrow::datatypes::Fields;

        let struct_fields = Fields::from(vec![
            Field::new("start", DataType::Int64, false),
            Field::new("end", DataType::Int64, false),
        ]);
        let list_field = Field::new("item", DataType::Struct(struct_fields.clone()), true);
        let list_dt = DataType::List(Arc::new(list_field));

        let null_arr = new_null_array(&list_dt, 1);
        let batch = RecordBatch::try_new(
            Arc::new(Schema::new(vec![Field::new("mirna", list_dt, true)])),
            vec![null_arr],
        )
        .unwrap();

        assert_eq!(read_mirna_regions(&batch, 0, 0), None);
    }

    // =======================================================================
    // AF_COLUMNS constant validation
    // =======================================================================

    #[test]
    fn test_af_columns_count() {
        // 1 global + 5 continental + 1 gnomADe global + 9 gnomADe sub + 1 gnomADg global + 10 gnomADg sub = 27
        assert_eq!(AF_COLUMNS.len(), 27);
    }

    #[test]
    fn test_af_columns_only_global_af_has_format_4f() {
        let formatted: Vec<&str> = AF_COLUMNS
            .iter()
            .filter(|c| c.format_4f)
            .map(|c| c.cache_col)
            .collect();
        assert_eq!(formatted, vec!["AF"]);
    }

    #[test]
    fn test_af_columns_emit_in_csq_count() {
        // Emitted in CSQ: AF, AFR, AMR, EAS, EUR, SAS, gnomADe, gnomADg = 8
        let emitted: Vec<&str> = AF_COLUMNS
            .iter()
            .filter(|c| c.emit_in_csq)
            .map(|c| c.cache_col)
            .collect();
        assert_eq!(emitted.len(), 8);
        assert!(emitted.contains(&"AF"));
        assert!(emitted.contains(&"gnomADe"));
        assert!(emitted.contains(&"gnomADg"));
    }

    #[test]
    fn test_af_columns_max_af_pop_excludes_globals() {
        // Globals (AF, gnomADe, gnomADg) should have max_af_pop = None.
        for col in AF_COLUMNS {
            if col.cache_col == "AF" || col.cache_col == "gnomADe" || col.cache_col == "gnomADg" {
                assert!(col.max_af_pop.is_none(), "{} should not have max_af_pop", col.cache_col);
            }
        }
    }

    #[test]
    fn test_af_columns_flag_groups() {
        // Group 0 = af (1 entry), 1 = af_1kg (5), 2 = af_gnomade (10), 3 = af_gnomadg (11)
        assert_eq!(AF_COLUMNS.iter().filter(|c| c.flag_group == 0).count(), 1);
        assert_eq!(AF_COLUMNS.iter().filter(|c| c.flag_group == 1).count(), 5);
        assert_eq!(AF_COLUMNS.iter().filter(|c| c.flag_group == 2).count(), 10);
        assert_eq!(AF_COLUMNS.iter().filter(|c| c.flag_group == 3).count(), 11);
    }

    // =======================================================================
    // CSQ format string pipe-count correctness
    //
    // The 74-field CSQ format is critical — if any format string has the wrong
    // number of pipes, all downstream fields shift and the annotation is corrupt.
    // These tests replicate the exact format!() patterns from
    // annotate_batch_with_transcript_engine() and verify field count + positions.
    // =======================================================================

    /// Build the batch3 suffix exactly as the production code does.
    fn build_test_batch3_suffix(
        af_vals: &[&str],
        max_af: &str,
        max_af_pops: &str,
        clin_sig: &str,
        somatic: &str,
        pheno: &str,
        pubmed: &str,
    ) -> String {
        let af_joined: Vec<String> = af_vals.iter().map(|s| s.to_string()).collect();
        format!(
            "{}|{}|{}|{clin_sig}|{somatic}|{pheno}|{pubmed}",
            af_joined.join("|"),
            max_af,
            max_af_pops,
        )
    }

    #[test]
    fn test_batch3_suffix_field_count() {
        let af_vals: Vec<&str> = vec![""; AF_COLUMNS.len()];
        let suffix = build_test_batch3_suffix(&af_vals, "", "", "", "", "", "");
        // batch3 has 33 fields (27 AF + max_af + max_af_pops + clin_sig + somatic + pheno + pubmed).
        let field_count = suffix.split('|').count();
        assert_eq!(field_count, 33, "batch3_suffix must have 33 fields, got {field_count}");
    }

    #[test]
    fn test_cache_hit_csq_has_74_fields() {
        let vep_allele = "C";
        let csq_val = "missense_variant";
        let impact = "MODERATE";
        let existing_var = "rs123";
        let variant_class = "SNV";
        let af_vals: Vec<&str> = vec![""; AF_COLUMNS.len()];
        let batch3_suffix = build_test_batch3_suffix(&af_vals, "", "", "", "", "", "");

        // This is the exact format from line 1841-1843.
        let csq = format!(
            "{vep_allele}|{csq_val}|{impact}|||||||||||||||{existing_var}||||||||||||\
             {variant_class}||||||||||||{batch3_suffix}"
        );
        let field_count = csq.split('|').count();
        assert_eq!(field_count, 74, "cache-hit CSQ must have 74 fields, got {field_count}");

        // Verify key field positions (0-indexed).
        let fields: Vec<&str> = csq.split('|').collect();
        assert_eq!(fields[0], "C", "field 0 = Allele");
        assert_eq!(fields[1], "missense_variant", "field 1 = Consequence");
        assert_eq!(fields[2], "MODERATE", "field 2 = IMPACT");
        assert_eq!(fields[17], "rs123", "field 17 = Existing_variation");
        assert_eq!(fields[29], "SNV", "field 29 = VARIANT_CLASS");
    }

    #[test]
    fn test_cache_miss_csq_has_74_fields() {
        let vep_allele = "G";
        let terms_str = "missense_variant&splice_region_variant";
        let tc_impact = "MODERATE";
        let symbol = "BRCA1";
        let gene = "ENSG00000012048";
        let feature_type = "Transcript";
        let feature = "ENST00000357654";
        let biotype = "protein_coding";
        let exon = "10/23";
        let intron = "";
        let hgvsc = "";
        let hgvsp = "";
        let cdna_pos = "1234";
        let cds_pos = "1001";
        let protein_pos = "334";
        let amino_acids = "V/G";
        let codons_str = "gTc/gGc";
        let existing_var = "rs12345";
        let distance = "";
        let strand_str = "1";
        let tc_flags = "";
        let symbol_source = "HGNC";
        let hgnc_id = "HGNC:1100";
        let source_val = "";
        let variant_class = "SNV";
        let canonical = "YES";
        let tsl_str = "1";
        let mane_select = "NM_007294.4";
        let mane_plus = "";
        let ensp = "ENSP00000350283";
        let gene_pheno = "1";
        let ccds = "CCDS11453.1";
        let swissprot = "P38398";
        let trembl = "";
        let uniparc = "UPI000002ED67";
        let uniprot_isoform = "";
        let af_vals: Vec<&str> = vec![""; AF_COLUMNS.len()];
        let batch3_suffix = build_test_batch3_suffix(&af_vals, "", "", "benign", "", "1", "");

        // This is the exact format from lines 1964-1972.
        let csq = format!(
            "{vep_allele}|{terms_str}|{tc_impact}|{symbol}|{gene}|{feature_type}|{feature}|{biotype}|\
             {exon}|{intron}|{hgvsc}|{hgvsp}|\
             {cdna_pos}|{cds_pos}|{protein_pos}|{amino_acids}|{codons_str}|\
             {existing_var}|{distance}|{strand_str}|{tc_flags}|{symbol_source}|{hgnc_id}|\
             |||||{source_val}|\
             {variant_class}|{canonical}|{tsl_str}|{mane_select}|{mane_plus}|\
             {ensp}|{gene_pheno}|{ccds}|{swissprot}|{trembl}|{uniparc}|{uniprot_isoform}|\
             {batch3_suffix}"
        );
        let fields: Vec<&str> = csq.split('|').collect();
        assert_eq!(
            fields.len(),
            74,
            "cache-miss CSQ must have 74 fields, got {}",
            fields.len()
        );

        // Verify field positions match VEP CSQ spec.
        assert_eq!(fields[0], "G", "field 0 = Allele");
        assert_eq!(fields[1], "missense_variant&splice_region_variant", "field 1 = Consequence");
        assert_eq!(fields[2], "MODERATE", "field 2 = IMPACT");
        assert_eq!(fields[3], "BRCA1", "field 3 = SYMBOL");
        assert_eq!(fields[4], "ENSG00000012048", "field 4 = Gene");
        assert_eq!(fields[5], "Transcript", "field 5 = Feature_type");
        assert_eq!(fields[6], "ENST00000357654", "field 6 = Feature");
        assert_eq!(fields[7], "protein_coding", "field 7 = BIOTYPE");
        assert_eq!(fields[8], "10/23", "field 8 = EXON");
        assert_eq!(fields[9], "", "field 9 = INTRON");
        assert_eq!(fields[10], "", "field 10 = HGVSc");
        assert_eq!(fields[11], "", "field 11 = HGVSp");
        assert_eq!(fields[12], "1234", "field 12 = cDNA_position");
        assert_eq!(fields[13], "1001", "field 13 = CDS_position");
        assert_eq!(fields[14], "334", "field 14 = Protein_position");
        assert_eq!(fields[15], "V/G", "field 15 = Amino_acids");
        assert_eq!(fields[16], "gTc/gGc", "field 16 = Codons");
        assert_eq!(fields[17], "rs12345", "field 17 = Existing_variation");
        assert_eq!(fields[18], "", "field 18 = DISTANCE");
        assert_eq!(fields[19], "1", "field 19 = STRAND");
        assert_eq!(fields[20], "", "field 20 = FLAGS");
        assert_eq!(fields[21], "HGNC", "field 21 = SYMBOL_SOURCE");
        assert_eq!(fields[22], "HGNC:1100", "field 22 = HGNC_ID");
        // Fields 23-27: MOTIF_NAME, MOTIF_POS, HIGH_INF_POS, MOTIF_SCORE_CHANGE, TRANSCRIPTION_FACTORS
        for i in 23..28 {
            assert_eq!(fields[i], "", "field {i} should be empty (motif/TF placeholder)");
        }
        assert_eq!(fields[28], "", "field 28 = SOURCE (empty when not merged)");
        assert_eq!(fields[29], "SNV", "field 29 = VARIANT_CLASS");
        assert_eq!(fields[30], "YES", "field 30 = CANONICAL");
        assert_eq!(fields[31], "1", "field 31 = TSL");
        assert_eq!(fields[32], "NM_007294.4", "field 32 = MANE_SELECT");
        assert_eq!(fields[33], "", "field 33 = MANE_PLUS_CLINICAL");
        assert_eq!(fields[34], "ENSP00000350283", "field 34 = ENSP");
        assert_eq!(fields[35], "1", "field 35 = GENE_PHENO");
        assert_eq!(fields[36], "CCDS11453.1", "field 36 = CCDS");
        assert_eq!(fields[37], "P38398", "field 37 = SWISSPROT");
        assert_eq!(fields[38], "", "field 38 = TREMBL");
        assert_eq!(fields[39], "UPI000002ED67", "field 39 = UNIPARC");
        assert_eq!(fields[40], "", "field 40 = UNIPROT_ISOFORM");
        // Fields 41-67: 27 AF values (all empty in this test).
        for i in 41..68 {
            assert_eq!(fields[i], "", "field {i} should be empty (AF placeholder)");
        }
        assert_eq!(fields[68], "", "field 68 = MAX_AF");
        assert_eq!(fields[69], "", "field 69 = MAX_AF_POPS");
        assert_eq!(fields[70], "benign", "field 70 = CLIN_SIG");
        assert_eq!(fields[71], "", "field 71 = SOMATIC");
        assert_eq!(fields[72], "1", "field 72 = PHENO");
        assert_eq!(fields[73], "", "field 73 = PUBMED");
    }

    #[test]
    fn test_fallback_csq_has_74_fields() {
        let vep_allele = "T";
        let impact = impact_label(SoImpact::Modifier);
        let existing_var = "";
        let variant_class = "SNV";
        let af_vals: Vec<&str> = vec![""; AF_COLUMNS.len()];
        let batch3_suffix = build_test_batch3_suffix(&af_vals, "", "", "", "", "", "");

        // This is the exact format from lines 1977-1979.
        let csq = format!(
            "{vep_allele}|sequence_variant|{impact}|||||||||||||||{existing_var}||||||||||||\
             {variant_class}||||||||||||{batch3_suffix}"
        );
        let field_count = csq.split('|').count();
        assert_eq!(field_count, 74, "fallback CSQ must have 74 fields, got {field_count}");

        let fields: Vec<&str> = csq.split('|').collect();
        assert_eq!(fields[0], "T", "field 0 = Allele");
        assert_eq!(fields[1], "sequence_variant", "field 1 = Consequence");
        assert_eq!(fields[2], "MODIFIER", "field 2 = IMPACT");
        assert_eq!(fields[29], "SNV", "field 29 = VARIANT_CLASS");
    }

    #[test]
    fn test_cache_hit_csq_and_cache_miss_csq_field_alignment() {
        // Verify that cache-hit and cache-miss CSQ format strings agree on field positions
        // for the fields they both populate: Allele(0), Consequence(1), IMPACT(2),
        // Existing_variation(17), VARIANT_CLASS(29), and all Batch 3 fields (41-73).
        let af_vals: Vec<&str> = vec![""; AF_COLUMNS.len()];
        let batch3 = build_test_batch3_suffix(&af_vals, "0.05", "AFR", "benign", "0&1", "1", "999");

        let hit = format!(
            "C|missense_variant|MODERATE|||||||||||||||rs1||||||||||||SNV||||||||||||{batch3}"
        );
        let miss = format!(
            "C|missense_variant|MODERATE|SYM|ENSG|Transcript|ENST|pc|\
             1/2|||||||||rs1||1||HGNC|H1|||||||\
             SNV|YES|1|NM|NP|ENSP|1|CCDS|SP||UP|UI|{batch3}"
        );
        let hit_fields: Vec<&str> = hit.split('|').collect();
        let miss_fields: Vec<&str> = miss.split('|').collect();
        assert_eq!(hit_fields.len(), 74);
        assert_eq!(miss_fields.len(), 74);

        // Shared fields must be at the same position.
        assert_eq!(hit_fields[0], miss_fields[0], "Allele");
        assert_eq!(hit_fields[1], miss_fields[1], "Consequence");
        assert_eq!(hit_fields[2], miss_fields[2], "IMPACT");
        assert_eq!(hit_fields[17], miss_fields[17], "Existing_variation");
        assert_eq!(hit_fields[29], miss_fields[29], "VARIANT_CLASS");
        // All batch3 fields (41-73) must be identical.
        for i in 41..74 {
            assert_eq!(hit_fields[i], miss_fields[i], "batch3 field {i}");
        }
    }

    #[test]
    fn test_batch3_suffix_with_populated_af_fields() {
        // Verify AF fields land in correct positions within batch3.
        let mut af_vals = vec![""; AF_COLUMNS.len()];
        af_vals[0] = "0.0301"; // AF (global, format_4f)
        af_vals[1] = "0.05";   // AFR
        af_vals[6] = "0.03";   // gnomADe (global)
        af_vals[16] = "0.02";  // gnomADg (global)

        let suffix = build_test_batch3_suffix(&af_vals, "0.05", "AFR", "benign", "0&1", "1", "123");
        let fields: Vec<&str> = suffix.split('|').collect();
        assert_eq!(fields.len(), 33);
        assert_eq!(fields[0], "0.0301", "AF");
        assert_eq!(fields[1], "0.05", "AFR");
        assert_eq!(fields[6], "0.03", "gnomADe");
        assert_eq!(fields[16], "0.02", "gnomADg");
        assert_eq!(fields[27], "0.05", "MAX_AF");
        assert_eq!(fields[28], "AFR", "MAX_AF_POPS");
        assert_eq!(fields[29], "benign", "CLIN_SIG");
        assert_eq!(fields[30], "0&1", "SOMATIC");
        assert_eq!(fields[31], "1", "PHENO");
        assert_eq!(fields[32], "123", "PUBMED");
    }

    #[test]
    fn test_source_field_only_populated_when_merged() {
        let af_vals: Vec<&str> = vec![""; AF_COLUMNS.len()];
        let batch3 = build_test_batch3_suffix(&af_vals, "", "", "", "", "", "");

        // Not merged: source_val = ""
        let csq_not_merged = format!(
            "C|x|HIGH|SYM|G|Transcript|ENST|pc||||||||||\
             rs1||1||HGNC|H1|||||||\
             SNV||||||||||||{batch3}"
        );
        let fields: Vec<&str> = csq_not_merged.split('|').collect();
        assert_eq!(fields.len(), 74);
        assert_eq!(fields[28], "", "SOURCE should be empty when not merged");

        // Merged: source_val = "Ensembl"
        let csq_merged = format!(
            "C|x|HIGH|SYM|G|Transcript|ENST|pc||||||||||\
             rs1||1||HGNC|H1||||||Ensembl|\
             SNV||||||||||||{batch3}"
        );
        let fields_merged: Vec<&str> = csq_merged.split('|').collect();
        assert_eq!(fields_merged.len(), 74);
        assert_eq!(fields_merged[28], "Ensembl", "SOURCE should be 'Ensembl' when merged");
    }

    #[test]
    fn test_csq_multiple_transcript_entries_comma_joined() {
        // When multiple transcripts match, entries are joined with comma.
        let entry1 = "C|missense_variant|MODERATE|A|B|Transcript|ENST1|pc||||||||||||||||||||||SNV||||||||||||||||||||||||||||||||||||||||||||";
        let entry2 = "C|synonymous_variant|LOW|A|B|Transcript|ENST2|pc||||||||||||||||||||||SNV||||||||||||||||||||||||||||||||||||||||||||";
        let joined = [entry1, entry2].join(",");
        let entries: Vec<&str> = joined.split(',').collect();
        assert_eq!(entries.len(), 2);
        // Each entry should have 74 fields.
        // Note: this breaks if any field value contains a comma — csq_escape handles that.
        assert_eq!(entries[0].split('|').count(), 74);
        assert_eq!(entries[1].split('|').count(), 74);
    }

    #[test]
    fn test_csq_escape_prevents_comma_splitting() {
        // TREMBL values like "B0QYZ8.91,X5DR28.81" would corrupt CSV-like CSQ splitting.
        // csq_escape replaces commas with '&'.
        let trembl_raw = "B0QYZ8.91,X5DR28.81";
        let trembl = csq_escape(trembl_raw);
        assert_eq!(trembl.as_ref(), "B0QYZ8.91&X5DR28.81");
        assert!(!trembl.contains(','), "escaped value must not contain commas");
    }

    // =======================================================================
    // Existing_variation fallback logic
    // =======================================================================

    #[test]
    fn test_existing_var_prefers_coloc_over_variation_name() {
        // When co-located data exists, use it; otherwise fall back to variation_name.
        let coloc = ColocatedData {
            existing_variation: "rs1&rs2".to_string(),
            entries: vec![],
        };
        let variation_name = "rs1";

        // With co-located data → use co-located.
        let existing = Some(&coloc).map(|c| c.existing_variation.as_str()).unwrap_or(variation_name);
        assert_eq!(existing, "rs1&rs2");

        // Without co-located data → fall back to variation_name.
        let existing_fallback: Option<&ColocatedData> = None;
        let existing2 = existing_fallback
            .map(|c| c.existing_variation.as_str())
            .unwrap_or(variation_name);
        assert_eq!(existing2, "rs1");
    }

    #[test]
    fn test_existing_var_empty_when_check_existing_false() {
        // When check_existing=false, Existing_variation should be empty.
        let flags = VepFlags::from_options_json(None);
        let existing = if flags.check_existing { "rs1" } else { "" };
        assert_eq!(existing, "");
    }

    // =======================================================================
    // chrom_filter_clause
    // =======================================================================

    #[test]
    fn test_chrom_filter_clause_empty_set() {
        let chroms: HashSet<String> = HashSet::new();
        let clause = AnnotateProvider::chrom_filter_clause(&chroms);
        assert_eq!(clause, "");
    }

    #[test]
    fn test_chrom_filter_clause_with_chr_prefix() {
        let mut chroms = HashSet::new();
        chroms.insert("chr22".to_string());
        let clause = AnnotateProvider::chrom_filter_clause(&chroms);
        // Should produce WHERE clause with both chr22 and 22.
        assert!(clause.contains("'chr22'"), "should include chr22");
        assert!(clause.contains("'22'"), "should include 22 without prefix");
    }

    #[test]
    fn test_chrom_filter_clause_without_chr_prefix() {
        let mut chroms = HashSet::new();
        chroms.insert("22".to_string());
        let clause = AnnotateProvider::chrom_filter_clause(&chroms);
        assert!(clause.contains("'22'"), "should include 22");
        assert!(clause.contains("'chr22'"), "should include chr22 with prefix");
    }
}
