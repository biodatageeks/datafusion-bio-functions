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
    Array, Float32Array, Float64Array, Int8Array, Int16Array, Int32Array, Int64Array,
    LargeStringArray, RecordBatch, StringArray, StringBuilder, StringViewArray, UInt8Array,
    UInt16Array, UInt32Array, UInt64Array, new_null_array,
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
            let raw_json_idx = schema.index_of("raw_object_json").ok();
            let gene_stable_id_idx = schema.index_of("gene_stable_id").ok();
            let gene_symbol_idx = schema.index_of("gene_symbol").ok();
            let gene_symbol_source_idx = schema.index_of("gene_symbol_source").ok();
            let gene_hgnc_id_idx = schema.index_of("gene_hgnc_id").ok();
            let source_idx = schema.index_of("source").ok();
            let version_idx = schema.index_of("version").ok();

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

                // Extract raw_object_json once — used for miRNA regions and FLAGS.
                let raw_json =
                    raw_json_idx.and_then(|idx| string_at(batch.column(idx).as_ref(), row));

                // For miRNA transcripts, extract mature miRNA regions from
                // raw_object_json attributes (cDNA coords mapped to genomic).
                let mature_mirna_regions = if biotype == "miRNA" {
                    parse_mirna_regions_from_json(raw_json.as_deref(), start, end, strand)
                } else {
                    Vec::new()
                };

                // Parse CDS incompleteness flags from raw_object_json attributes.
                let (cds_start_nf, cds_end_nf, flags_str) =
                    parse_transcript_flags(raw_json.as_deref());

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
            let vep_allele = {
                let ref_al = string_at(batch.column(ref_idx).as_ref(), row).unwrap_or_default();
                let (_vep_ref, vep_alt) = vcf_to_vep_allele(&ref_al, &alt_allele);
                vep_alt
            };

            // Cache-hit fast path: use pre-computed consequence from variation cache.
            let cached_most =
                cached_most_idx.and_then(|idx| string_at(batch.column(idx).as_ref(), row));
            let cached_csq =
                cached_csq_idx.and_then(|idx| string_at(batch.column(idx).as_ref(), row));

            let variation_name = variation_name_idx
                .and_then(|idx| string_at(batch.column(idx).as_ref(), row))
                .unwrap_or_default();

            let (csq_string, most_str) = if let Some(most_val) = &cached_most {
                // Cache hit — produce single CSQ entry with empty transcript fields.
                let csq_val = cached_csq.unwrap_or_default();
                let impact = SoTerm::from_str(most_val)
                    .map(|t| impact_label(t.impact()))
                    .unwrap_or_else(|| impact_label(SoImpact::Modifier));
                // 29-field CSQ: Allele|Consequence|IMPACT|SYMBOL|Gene|Feature_type|Feature|BIOTYPE|
                // EXON|INTRON|HGVSc|HGVSp|cDNA_pos|CDS_pos|Protein_pos|AA|Codons|
                // Existing_variation|DISTANCE|STRAND|FLAGS|SYMBOL_SOURCE|HGNC_ID|
                // MOTIF_NAME|MOTIF_POS|HIGH_INF_POS|MOTIF_SCORE_CHANGE|TRANSCRIPTION_FACTORS|SOURCE
                let csq_entry =
                    format!("{vep_allele}|{csq_val}|{impact}|||||||||||||||||{variation_name}||||||||||");
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
                    let (symbol, gene, biotype_tx, strand_str, symbol_source, hgnc_id, source) =
                        if let Some(idx) = tc.transcript_idx {
                            let tx = ctx.transcripts[idx];
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
                    let flags = tc.flags.as_deref().unwrap_or("");
                    // HGVSc/HGVSp require VEP's --hgvs flag + FASTA; leave
                    // empty by default to preserve schema compatibility.
                    let hgvsc = "";
                    let hgvsp = "";
                    // VEP only emits SOURCE when using --merged cache.
                    let source_val = if merged { source } else { "" };
                    csq_entries.push(format!(
                        "{vep_allele}|{terms_str}|{tc_impact}|{symbol}|{gene}|{feature_type}|{feature}|{biotype}|\
                         {exon}|{intron}|{hgvsc}|{hgvsp}|\
                         {cdna_pos}|{cds_pos}|{protein_pos}|{amino_acids}|{codons_str}|\
                         |{distance}|{strand_str}|{flags}|{symbol_source}|{hgnc_id}|\
                         |||||{source_val}"
                    ));
                }
                if csq_entries.is_empty() {
                    let impact = impact_label(SoImpact::Modifier);
                    csq_entries.push(format!(
                        "{vep_allele}|sequence_variant|{impact}|||||||||||||||||||||||||||"
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
fn parse_mirna_regions_from_json(
    raw_json: Option<&str>,
    tx_start: i64,
    tx_end: i64,
    strand: i8,
) -> Vec<(i64, i64)> {
    let Some(json_str) = raw_json else {
        return Vec::new();
    };

    // Look for patterns like {"code":"miRNA","value":"42-59"} in the JSON.
    // We do a lightweight parse to avoid pulling in a full JSON library just
    // for this single use case.
    let mut regions = Vec::new();
    // Find all occurrences of "miRNA" code entries and extract their values.
    let needle = "\"miRNA\"";
    let mut search_from = 0;
    while let Some(pos) = json_str[search_from..].find(needle) {
        let abs_pos = search_from + pos + needle.len();
        search_from = abs_pos;

        // Look for the "value" key nearby (within the same JSON object).
        let window_end = (abs_pos + 200).min(json_str.len());
        let window = &json_str[abs_pos..window_end];

        // Find "value" key
        if let Some(val_pos) = window.find("\"value\"") {
            let after_key = &window[val_pos + 7..]; // skip "value"
            // Find the value string: skip colon and whitespace, then quoted string
            if let Some(quote_start) = after_key.find('"') {
                let val_start = quote_start + 1;
                if let Some(quote_end) = after_key[val_start..].find('"') {
                    let val_str = &after_key[val_start..val_start + quote_end];
                    // Parse "N-M" pattern
                    if let Some((cdna_start, cdna_end)) = parse_cdna_range(val_str) {
                        let (gstart, gend) = if strand >= 0 {
                            (tx_start + cdna_start - 1, tx_start + cdna_end - 1)
                        } else {
                            (tx_end - cdna_end + 1, tx_end - cdna_start + 1)
                        };
                        regions.push((gstart, gend));
                    }
                }
            }
        }
    }
    regions
}

/// Parse `cds_start_NF` and `cds_end_NF` flags from `raw_object_json` attributes.
/// These indicate incomplete CDS (5' or 3' truncation).
/// Only returns true when the attribute has `"value":"1"` (not just presence of the key).
/// Returns (start_nf, end_nf, flags_str) where flags_str preserves VEP's encounter order.
fn parse_transcript_flags(raw_json: Option<&str>) -> (bool, bool, Option<String>) {
    let Some(json_str) = raw_json else {
        return (false, false, None);
    };
    let Ok(val) = serde_json::from_str::<serde_json::Value>(json_str) else {
        return (false, false, None);
    };
    let attrs = val
        .pointer("/__value/attributes")
        .and_then(|a| a.as_array());
    let Some(attrs) = attrs else {
        return (false, false, None);
    };
    let mut start_nf = false;
    let mut end_nf = false;
    let mut flags: Vec<&str> = Vec::new();
    for attr in attrs {
        // Each attribute may be wrapped: {"__class":"...","__value":{"code":"...","value":"..."}}
        // Unwrap to the inner object if present, otherwise use the attribute directly.
        let inner = attr.get("__value").unwrap_or(attr);
        let code = inner.get("code").and_then(|c| c.as_str()).unwrap_or("");
        let value = inner.get("value").and_then(|v| v.as_str()).unwrap_or("");
        if code == "cds_start_NF" && value == "1" {
            start_nf = true;
            flags.push("cds_start_NF");
        }
        if code == "cds_end_NF" && value == "1" {
            end_nf = true;
            flags.push("cds_end_NF");
        }
    }
    let flags_str = if flags.is_empty() {
        None
    } else {
        Some(flags.join("&"))
    };
    (start_nf, end_nf, flags_str)
}

fn parse_cdna_range(s: &str) -> Option<(i64, i64)> {
    let parts: Vec<&str> = s.split('-').collect();
    if parts.len() != 2 {
        return None;
    }
    let start: i64 = parts[0].trim().parse().ok()?;
    let end: i64 = parts[1].trim().parse().ok()?;
    Some((start, end))
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
    fn test_parse_flags_cds_start_nf_value_1() {
        let json = r#"{"__value":{"attributes":[{"code":"cds_start_NF","value":"1"}]}}"#;
        assert_eq!(
            parse_transcript_flags(Some(json)),
            (true, false, Some("cds_start_NF".to_string()))
        );
    }

    #[test]
    fn test_parse_flags_cds_start_nf_value_0() {
        let json = r#"{"__value":{"attributes":[{"code":"cds_start_NF","value":"0"}]}}"#;
        assert_eq!(parse_transcript_flags(Some(json)), (false, false, None));
    }

    #[test]
    fn test_parse_flags_both_nf_start_first() {
        let json = r#"{"__value":{"attributes":[{"code":"cds_start_NF","value":"1"},{"code":"cds_end_NF","value":"1"}]}}"#;
        assert_eq!(
            parse_transcript_flags(Some(json)),
            (true, true, Some("cds_start_NF&cds_end_NF".to_string()))
        );
    }

    #[test]
    fn test_parse_flags_both_nf_end_first() {
        // VEP preserves encounter order: end before start.
        let json = r#"{"__value":{"attributes":[{"code":"cds_end_NF","value":"1"},{"code":"cds_start_NF","value":"1"}]}}"#;
        assert_eq!(
            parse_transcript_flags(Some(json)),
            (true, true, Some("cds_end_NF&cds_start_NF".to_string()))
        );
    }

    #[test]
    fn test_parse_flags_no_attributes() {
        let json = r#"{"__value":{}}"#;
        assert_eq!(parse_transcript_flags(Some(json)), (false, false, None));
    }

    #[test]
    fn test_parse_flags_none_input() {
        assert_eq!(parse_transcript_flags(None), (false, false, None));
    }

    #[test]
    fn test_parse_flags_nested_attribute_format() {
        let json = r#"{"__value":{"attributes":[{"__class":"Bio::EnsEMBL::Attribute","__value":{"code":"cds_start_NF","value":"1"}},{"__class":"Bio::EnsEMBL::Attribute","__value":{"code":"cds_end_NF","value":"1"}}]}}"#;
        assert_eq!(
            parse_transcript_flags(Some(json)),
            (true, true, Some("cds_start_NF&cds_end_NF".to_string()))
        );
    }

    #[test]
    fn test_parse_flags_nested_single_flag() {
        let json = r#"{"__value":{"attributes":[{"__class":"Bio::EnsEMBL::Attribute","__value":{"code":"cds_end_NF","value":"1"}}]}}"#;
        assert_eq!(
            parse_transcript_flags(Some(json)),
            (false, true, Some("cds_end_NF".to_string()))
        );
    }

    #[test]
    fn test_parse_flags_mixed_wrapped_unwrapped() {
        // Mix of wrapped (__class/__value) and unwrapped attributes
        let json = r#"{"__value":{"attributes":[{"code":"cds_start_NF","value":"1"},{"__class":"Bio::EnsEMBL::Attribute","__value":{"code":"cds_end_NF","value":"1"}}]}}"#;
        assert_eq!(
            parse_transcript_flags(Some(json)),
            (true, true, Some("cds_start_NF&cds_end_NF".to_string()))
        );
    }

    #[test]
    fn test_parse_flags_non_1_value_rejected() {
        // value="true" should NOT be accepted — only "1" counts
        let json = r#"{"__value":{"attributes":[{"code":"cds_start_NF","value":"true"}]}}"#;
        assert_eq!(parse_transcript_flags(Some(json)), (false, false, None));
    }

    #[test]
    fn test_parse_flags_ignores_non_cds_attributes() {
        // Non-cds_start_NF/cds_end_NF attributes should be ignored
        let json = r#"{"__value":{"attributes":[{"code":"some_other","value":"1"},{"code":"cds_start_NF","value":"1"}]}}"#;
        assert_eq!(
            parse_transcript_flags(Some(json)),
            (true, false, Some("cds_start_NF".to_string()))
        );
    }

    #[test]
    fn test_parse_flags_empty_attributes_array() {
        let json = r#"{"__value":{"attributes":[]}}"#;
        assert_eq!(parse_transcript_flags(Some(json)), (false, false, None));
    }

    #[test]
    fn test_parse_flags_malformed_json() {
        assert_eq!(
            parse_transcript_flags(Some("not valid json")),
            (false, false, None)
        );
    }
}
