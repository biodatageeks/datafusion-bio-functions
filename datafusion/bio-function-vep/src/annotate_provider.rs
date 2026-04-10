//! Provider for `annotate_vep()` table function.
//!
//! Runtime behavior:
//! - always starts from `lookup_variants()` for known-variant metadata,
//! - when transcript/exon tables are available, computes transcript-driven
//!   consequence terms and most-severe ranking,
//! - otherwise falls back to phase-1.5 known-variant CSQ placeholders.

use std::any::Any;
use std::cmp::Ordering;
use std::collections::HashMap;
use std::collections::HashSet;
use std::collections::VecDeque;
use std::collections::hash_map::DefaultHasher;
use std::fmt::{Debug, Formatter};
use std::hash::{Hash, Hasher};
use std::io::{BufRead, Seek};
use std::pin::Pin;
use std::sync::{Arc, Mutex};
use std::task::{Context as TaskCtx, Poll};
use std::time::Instant;

/// Returns true when VEP_PROFILE env var is set (any value).
fn profiling_enabled() -> bool {
    std::env::var("VEP_PROFILE").is_ok()
}

macro_rules! profile_start {
    () => {
        Instant::now()
    };
}

macro_rules! profile_end {
    ($label:expr, $start:expr) => {
        if profiling_enabled() {
            let elapsed = $start.elapsed();
            eprintln!(
                "[VEP_PROFILE] {:.<50} {:>8.1}ms",
                $label,
                elapsed.as_secs_f64() * 1000.0
            );
        }
    };
    ($label:expr, $start:expr, $extra:expr) => {
        if profiling_enabled() {
            let elapsed = $start.elapsed();
            eprintln!(
                "[VEP_PROFILE] {:.<50} {:>8.1}ms  {}",
                $label,
                elapsed.as_secs_f64() * 1000.0,
                $extra
            );
        }
    };
}

use async_trait::async_trait;
use datafusion::arrow::array::{
    Array, AsArray, BooleanArray, Float32Array, Float32Builder, Float64Array, Int8Array,
    Int8Builder, Int16Array, Int32Array, Int64Array, Int64Builder, LargeStringArray, ListArray,
    ListBuilder, RecordBatch, StringArray, StringBuilder, StringViewArray, UInt8Array, UInt16Array,
    UInt32Array, UInt64Array, new_null_array,
};
use datafusion::arrow::datatypes::{DataType, Field, Schema, SchemaRef};
use datafusion::catalog::Session;
use datafusion::common::{DataFusionError, Result};
use datafusion::datasource::{MemTable, TableProvider, TableType};
use datafusion::execution::{RecordBatchStream, SendableRecordBatchStream, TaskContext};
use datafusion::physical_expr::EquivalenceProperties;
use datafusion::physical_plan::ExecutionPlan;
use datafusion::physical_plan::execution_plan::{Boundedness, EmissionType};
use datafusion::physical_plan::{
    DisplayAs, DisplayFormatType, ExecutionPlanProperties, PlanProperties,
};
use datafusion::prelude::{Expr, ParquetReadOptions, SessionContext, col, lit};
use futures::{Future, Stream};
use noodles_core::{Position, Region};
use noodles_fasta as fasta;
use serde_json::Value;

use crate::allele::{
    MatchedVariantAllele, vcf_to_vep_allele, vcf_to_vep_input_allele, vep_norm_end, vep_norm_start,
};
use crate::annotation_store::{AnnotationBackend, build_store};
use crate::config;
#[cfg(feature = "kv-cache")]
use crate::kv_cache::KvCacheTableProvider;
use crate::lookup_provider::LookupProvider;
use crate::miss_worklist::{MissWorklist, collect_miss_worklist};
use crate::partitioned_cache::PartitionedParquetCache;
use crate::so_terms::{SoImpact, SoTerm, most_severe_term};
use crate::transcript_consequence::{
    CachedPredictions, CompactPrediction, ExonFeature, MirnaFeature, MotifFeature, PreparedContext,
    ProteinDomainFeature, RefSeqEdit, RegulatoryFeature, SiftPolyphenCache, StructuralFeature,
    SvEventKind, SvFeatureKind, TranscriptCdnaMapperSegment, TranscriptConsequence,
    TranscriptConsequenceEngine, TranscriptFeature, TranslationFeature, VariantInput,
    infer_refseq_deletion_edits_from_sequences, refseq_edit_offset_delta,
};
use crate::variant_lookup_exec::{
    ColocatedCacheEntry, ColocatedKey, ColocatedSink, ColocatedSinkValue,
};

/// Column categories for typed non-meta annotation columns.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
enum AnnotationCategory {
    /// Transcript-level fields from the most-severe transcript (42 columns).
    Transcript,
    /// Population allele frequency fields, resolved for the matching allele (29 columns).
    Frequency,
    /// Variant-level annotation columns from co-located data (9 columns).
    Variant,
    /// Cache-only columns not present in the CSQ string (7 columns).
    CacheOnly,
}

/// Metadata for a single annotation output column.
struct AnnotationColumnDef {
    /// Output column name as it appears in the Arrow schema.
    name: &'static str,
    /// Arrow data type for this column.
    data_type: DataType,
    /// Which category this column belongs to.
    category: AnnotationCategory,
    /// For cache-only/frequency columns: the column name in the variation cache
    /// parquet (prefixed with `cache_` in the intermediate batch).
    /// `None` for transcript-level columns populated from engine results.
    cache_col: Option<&'static str>,
}

fn list_utf8_data_type() -> DataType {
    DataType::List(Arc::new(Field::new("item", DataType::Utf8, true)))
}

fn list_i8_data_type() -> DataType {
    DataType::List(Arc::new(Field::new("item", DataType::Int8, true)))
}

fn list_i64_data_type() -> DataType {
    DataType::List(Arc::new(Field::new("item", DataType::Int64, true)))
}

/// Full default annotation output column definitions (87 columns, excluding csq and most_severe_consequence).
///
/// Order matches the README schema: transcript-level (42), frequency (29), variant-level (9), cache-only (7).
fn annotation_column_defs() -> Vec<AnnotationColumnDef> {
    use AnnotationCategory::*;
    vec![
        // ── Transcript-level (42) ──
        AnnotationColumnDef {
            name: "Allele",
            data_type: DataType::Utf8,
            category: Transcript,
            cache_col: None,
        },
        AnnotationColumnDef {
            name: "Consequence",
            data_type: list_utf8_data_type(),
            category: Transcript,
            cache_col: None,
        },
        AnnotationColumnDef {
            name: "IMPACT",
            data_type: list_utf8_data_type(),
            category: Transcript,
            cache_col: None,
        },
        AnnotationColumnDef {
            name: "SYMBOL",
            data_type: list_utf8_data_type(),
            category: Transcript,
            cache_col: None,
        },
        AnnotationColumnDef {
            name: "Gene",
            data_type: list_utf8_data_type(),
            category: Transcript,
            cache_col: None,
        },
        AnnotationColumnDef {
            name: "Feature_type",
            data_type: list_utf8_data_type(),
            category: Transcript,
            cache_col: None,
        },
        AnnotationColumnDef {
            name: "Feature",
            data_type: list_utf8_data_type(),
            category: Transcript,
            cache_col: None,
        },
        AnnotationColumnDef {
            name: "BIOTYPE",
            data_type: list_utf8_data_type(),
            category: Transcript,
            cache_col: None,
        },
        AnnotationColumnDef {
            name: "EXON",
            data_type: list_utf8_data_type(),
            category: Transcript,
            cache_col: None,
        },
        AnnotationColumnDef {
            name: "INTRON",
            data_type: list_utf8_data_type(),
            category: Transcript,
            cache_col: None,
        },
        AnnotationColumnDef {
            name: "HGVSc",
            data_type: list_utf8_data_type(),
            category: Transcript,
            cache_col: None,
        },
        AnnotationColumnDef {
            name: "HGVSp",
            data_type: list_utf8_data_type(),
            category: Transcript,
            cache_col: None,
        },
        AnnotationColumnDef {
            name: "cDNA_position",
            data_type: list_utf8_data_type(),
            category: Transcript,
            cache_col: None,
        },
        AnnotationColumnDef {
            name: "CDS_position",
            data_type: list_utf8_data_type(),
            category: Transcript,
            cache_col: None,
        },
        AnnotationColumnDef {
            name: "Protein_position",
            data_type: list_utf8_data_type(),
            category: Transcript,
            cache_col: None,
        },
        AnnotationColumnDef {
            name: "Amino_acids",
            data_type: list_utf8_data_type(),
            category: Transcript,
            cache_col: None,
        },
        AnnotationColumnDef {
            name: "Codons",
            data_type: list_utf8_data_type(),
            category: Transcript,
            cache_col: None,
        },
        AnnotationColumnDef {
            name: "Existing_variation",
            data_type: list_utf8_data_type(),
            category: Transcript,
            cache_col: None,
        },
        AnnotationColumnDef {
            name: "DISTANCE",
            data_type: list_i64_data_type(),
            category: Transcript,
            cache_col: None,
        },
        AnnotationColumnDef {
            name: "STRAND",
            data_type: list_i8_data_type(),
            category: Transcript,
            cache_col: None,
        },
        AnnotationColumnDef {
            name: "FLAGS",
            data_type: list_utf8_data_type(),
            category: Transcript,
            cache_col: None,
        },
        AnnotationColumnDef {
            name: "VARIANT_CLASS",
            data_type: DataType::Utf8,
            category: Transcript,
            cache_col: None,
        },
        AnnotationColumnDef {
            name: "SYMBOL_SOURCE",
            data_type: list_utf8_data_type(),
            category: Transcript,
            cache_col: None,
        },
        AnnotationColumnDef {
            name: "HGNC_ID",
            data_type: list_utf8_data_type(),
            category: Transcript,
            cache_col: None,
        },
        AnnotationColumnDef {
            name: "CANONICAL",
            data_type: list_utf8_data_type(),
            category: Transcript,
            cache_col: None,
        },
        AnnotationColumnDef {
            name: "MANE",
            data_type: list_utf8_data_type(),
            category: Transcript,
            cache_col: None,
        },
        AnnotationColumnDef {
            name: "MANE_SELECT",
            data_type: list_utf8_data_type(),
            category: Transcript,
            cache_col: None,
        },
        AnnotationColumnDef {
            name: "MANE_PLUS_CLINICAL",
            data_type: list_utf8_data_type(),
            category: Transcript,
            cache_col: None,
        },
        AnnotationColumnDef {
            name: "TSL",
            data_type: list_i8_data_type(),
            category: Transcript,
            cache_col: None,
        },
        AnnotationColumnDef {
            name: "APPRIS",
            data_type: list_utf8_data_type(),
            category: Transcript,
            cache_col: None,
        },
        AnnotationColumnDef {
            name: "CCDS",
            data_type: list_utf8_data_type(),
            category: Transcript,
            cache_col: None,
        },
        AnnotationColumnDef {
            name: "ENSP",
            data_type: list_utf8_data_type(),
            category: Transcript,
            cache_col: None,
        },
        AnnotationColumnDef {
            name: "SWISSPROT",
            data_type: list_utf8_data_type(),
            category: Transcript,
            cache_col: None,
        },
        AnnotationColumnDef {
            name: "TREMBL",
            data_type: list_utf8_data_type(),
            category: Transcript,
            cache_col: None,
        },
        AnnotationColumnDef {
            name: "UNIPARC",
            data_type: list_utf8_data_type(),
            category: Transcript,
            cache_col: None,
        },
        AnnotationColumnDef {
            name: "UNIPROT_ISOFORM",
            data_type: list_utf8_data_type(),
            category: Transcript,
            cache_col: None,
        },
        AnnotationColumnDef {
            name: "GENE_PHENO",
            data_type: list_utf8_data_type(),
            category: Transcript,
            cache_col: None,
        },
        AnnotationColumnDef {
            name: "SIFT",
            data_type: list_utf8_data_type(),
            category: Transcript,
            cache_col: None,
        },
        AnnotationColumnDef {
            name: "PolyPhen",
            data_type: list_utf8_data_type(),
            category: Transcript,
            cache_col: None,
        },
        AnnotationColumnDef {
            name: "DOMAINS",
            data_type: list_utf8_data_type(),
            category: Transcript,
            cache_col: None,
        },
        AnnotationColumnDef {
            name: "miRNA",
            data_type: list_utf8_data_type(),
            category: Transcript,
            cache_col: None,
        },
        AnnotationColumnDef {
            name: "HGVS_OFFSET",
            data_type: list_i64_data_type(),
            category: Transcript,
            cache_col: None,
        },
        // ── Frequency (29) ──
        AnnotationColumnDef {
            name: "AF",
            data_type: DataType::Float32,
            category: Frequency,
            cache_col: Some("AF"),
        },
        AnnotationColumnDef {
            name: "AFR_AF",
            data_type: DataType::Float32,
            category: Frequency,
            cache_col: Some("AFR"),
        },
        AnnotationColumnDef {
            name: "AMR_AF",
            data_type: DataType::Float32,
            category: Frequency,
            cache_col: Some("AMR"),
        },
        AnnotationColumnDef {
            name: "EAS_AF",
            data_type: DataType::Float32,
            category: Frequency,
            cache_col: Some("EAS"),
        },
        AnnotationColumnDef {
            name: "EUR_AF",
            data_type: DataType::Float32,
            category: Frequency,
            cache_col: Some("EUR"),
        },
        AnnotationColumnDef {
            name: "SAS_AF",
            data_type: DataType::Float32,
            category: Frequency,
            cache_col: Some("SAS"),
        },
        AnnotationColumnDef {
            name: "gnomADe_AF",
            data_type: DataType::Float32,
            category: Frequency,
            cache_col: Some("gnomADe"),
        },
        AnnotationColumnDef {
            name: "gnomADe_AFR_AF",
            data_type: DataType::Float32,
            category: Frequency,
            cache_col: Some("gnomADe_AFR"),
        },
        AnnotationColumnDef {
            name: "gnomADe_AMR_AF",
            data_type: DataType::Float32,
            category: Frequency,
            cache_col: Some("gnomADe_AMR"),
        },
        AnnotationColumnDef {
            name: "gnomADe_ASJ_AF",
            data_type: DataType::Float32,
            category: Frequency,
            cache_col: Some("gnomADe_ASJ"),
        },
        AnnotationColumnDef {
            name: "gnomADe_EAS_AF",
            data_type: DataType::Float32,
            category: Frequency,
            cache_col: Some("gnomADe_EAS"),
        },
        AnnotationColumnDef {
            name: "gnomADe_FIN_AF",
            data_type: DataType::Float32,
            category: Frequency,
            cache_col: Some("gnomADe_FIN"),
        },
        AnnotationColumnDef {
            name: "gnomADe_MID_AF",
            data_type: DataType::Float32,
            category: Frequency,
            cache_col: Some("gnomADe_MID"),
        },
        AnnotationColumnDef {
            name: "gnomADe_NFE_AF",
            data_type: DataType::Float32,
            category: Frequency,
            cache_col: Some("gnomADe_NFE"),
        },
        AnnotationColumnDef {
            name: "gnomADe_REMAINING_AF",
            data_type: DataType::Float32,
            category: Frequency,
            cache_col: Some("gnomADe_REMAINING"),
        },
        AnnotationColumnDef {
            name: "gnomADe_SAS_AF",
            data_type: DataType::Float32,
            category: Frequency,
            cache_col: Some("gnomADe_SAS"),
        },
        AnnotationColumnDef {
            name: "gnomADg_AF",
            data_type: DataType::Float32,
            category: Frequency,
            cache_col: Some("gnomADg"),
        },
        AnnotationColumnDef {
            name: "gnomADg_AFR_AF",
            data_type: DataType::Float32,
            category: Frequency,
            cache_col: Some("gnomADg_AFR"),
        },
        AnnotationColumnDef {
            name: "gnomADg_AMI_AF",
            data_type: DataType::Float32,
            category: Frequency,
            cache_col: Some("gnomADg_AMI"),
        },
        AnnotationColumnDef {
            name: "gnomADg_AMR_AF",
            data_type: DataType::Float32,
            category: Frequency,
            cache_col: Some("gnomADg_AMR"),
        },
        AnnotationColumnDef {
            name: "gnomADg_ASJ_AF",
            data_type: DataType::Float32,
            category: Frequency,
            cache_col: Some("gnomADg_ASJ"),
        },
        AnnotationColumnDef {
            name: "gnomADg_EAS_AF",
            data_type: DataType::Float32,
            category: Frequency,
            cache_col: Some("gnomADg_EAS"),
        },
        AnnotationColumnDef {
            name: "gnomADg_FIN_AF",
            data_type: DataType::Float32,
            category: Frequency,
            cache_col: Some("gnomADg_FIN"),
        },
        AnnotationColumnDef {
            name: "gnomADg_MID_AF",
            data_type: DataType::Float32,
            category: Frequency,
            cache_col: Some("gnomADg_MID"),
        },
        AnnotationColumnDef {
            name: "gnomADg_NFE_AF",
            data_type: DataType::Float32,
            category: Frequency,
            cache_col: Some("gnomADg_NFE"),
        },
        AnnotationColumnDef {
            name: "gnomADg_REMAINING_AF",
            data_type: DataType::Float32,
            category: Frequency,
            cache_col: Some("gnomADg_REMAINING"),
        },
        AnnotationColumnDef {
            name: "gnomADg_SAS_AF",
            data_type: DataType::Float32,
            category: Frequency,
            cache_col: Some("gnomADg_SAS"),
        },
        AnnotationColumnDef {
            name: "MAX_AF",
            data_type: DataType::Float32,
            category: Frequency,
            cache_col: None,
        },
        AnnotationColumnDef {
            name: "MAX_AF_POPS",
            data_type: DataType::Utf8,
            category: Frequency,
            cache_col: None,
        },
        // ── Variant-level (9) ──
        AnnotationColumnDef {
            name: "CLIN_SIG",
            data_type: list_utf8_data_type(),
            category: Variant,
            cache_col: None,
        },
        AnnotationColumnDef {
            name: "SOMATIC",
            data_type: DataType::Utf8,
            category: Variant,
            cache_col: None,
        },
        AnnotationColumnDef {
            name: "PHENO",
            data_type: DataType::Utf8,
            category: Variant,
            cache_col: None,
        },
        AnnotationColumnDef {
            name: "PUBMED",
            data_type: list_utf8_data_type(),
            category: Variant,
            cache_col: None,
        },
        AnnotationColumnDef {
            name: "MOTIF_NAME",
            data_type: DataType::Utf8,
            category: Variant,
            cache_col: None,
        },
        AnnotationColumnDef {
            name: "MOTIF_POS",
            data_type: DataType::Utf8,
            category: Variant,
            cache_col: None,
        },
        AnnotationColumnDef {
            name: "HIGH_INF_POS",
            data_type: DataType::Utf8,
            category: Variant,
            cache_col: None,
        },
        AnnotationColumnDef {
            name: "MOTIF_SCORE_CHANGE",
            data_type: DataType::Float32,
            category: Variant,
            cache_col: None,
        },
        AnnotationColumnDef {
            name: "TRANSCRIPTION_FACTORS",
            data_type: list_utf8_data_type(),
            category: Variant,
            cache_col: None,
        },
        // ── Cache-only (7) ──
        AnnotationColumnDef {
            name: "clin_sig_allele",
            data_type: list_utf8_data_type(),
            category: CacheOnly,
            cache_col: Some("clin_sig_allele"),
        },
        AnnotationColumnDef {
            name: "clinical_impact",
            data_type: DataType::Utf8,
            category: CacheOnly,
            cache_col: Some("clinical_impact"),
        },
        AnnotationColumnDef {
            name: "minor_allele",
            data_type: DataType::Utf8,
            category: CacheOnly,
            cache_col: Some("minor_allele"),
        },
        AnnotationColumnDef {
            name: "minor_allele_freq",
            data_type: DataType::Float32,
            category: CacheOnly,
            cache_col: Some("minor_allele_freq"),
        },
        AnnotationColumnDef {
            name: "clinvar_ids",
            data_type: list_utf8_data_type(),
            category: CacheOnly,
            cache_col: Some("clinvar_ids"),
        },
        AnnotationColumnDef {
            name: "cosmic_ids",
            data_type: list_utf8_data_type(),
            category: CacheOnly,
            cache_col: Some("cosmic_ids"),
        },
        AnnotationColumnDef {
            name: "dbsnp_ids",
            data_type: list_utf8_data_type(),
            category: CacheOnly,
            cache_col: Some("dbsnp_ids"),
        },
    ]
}

fn refseq_annotation_column_defs(include_source: bool) -> Vec<AnnotationColumnDef> {
    use AnnotationCategory::Transcript;

    let mut defs = annotation_column_defs();
    let insert_at = defs
        .iter()
        .position(|def| def.name == "GENE_PHENO")
        .unwrap_or(defs.len());
    let mut refseq_defs = vec![AnnotationColumnDef {
        name: "REFSEQ_MATCH",
        data_type: list_utf8_data_type(),
        category: Transcript,
        cache_col: None,
    }];
    if include_source {
        refseq_defs.push(AnnotationColumnDef {
            name: "SOURCE",
            data_type: list_utf8_data_type(),
            category: Transcript,
            cache_col: None,
        });
    }
    refseq_defs.extend([
        AnnotationColumnDef {
            name: "REFSEQ_OFFSET",
            data_type: list_i64_data_type(),
            category: Transcript,
            cache_col: None,
        },
        AnnotationColumnDef {
            name: "GIVEN_REF",
            data_type: list_utf8_data_type(),
            category: Transcript,
            cache_col: None,
        },
        AnnotationColumnDef {
            name: "USED_REF",
            data_type: list_utf8_data_type(),
            category: Transcript,
            cache_col: None,
        },
        AnnotationColumnDef {
            name: "BAM_EDIT",
            data_type: list_utf8_data_type(),
            category: Transcript,
            cache_col: None,
        },
    ]);
    defs.splice(insert_at..insert_at, refseq_defs);
    defs
}

fn annotation_column_defs_for_selection(
    transcript_selection: TranscriptSelectionFlags,
) -> Vec<AnnotationColumnDef> {
    if transcript_selection.refseq_fields() {
        refseq_annotation_column_defs(transcript_selection.source_field())
    } else {
        annotation_column_defs()
    }
}

/// Returns the list of cache column names needed for the variation lookup query.
///
/// This is the backward-compatible list of column names that the `requested_columns`
/// logic uses to select columns from the variation cache parquet/fjall store.
pub fn cache_lookup_column_names() -> Vec<&'static str> {
    vec![
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
    ]
}

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
    AfColumn {
        cache_col: "AF",
        format_4f: true,
        flag_group: 0,
        emit_in_csq: true,
        max_af_pop: None,
    },
    // --af_1kg (continental) — emitted, MAX_AF uses short names (AFR not AFR_AF)
    AfColumn {
        cache_col: "AFR",
        format_4f: false,
        flag_group: 1,
        emit_in_csq: true,
        max_af_pop: Some("AFR"),
    },
    AfColumn {
        cache_col: "AMR",
        format_4f: false,
        flag_group: 1,
        emit_in_csq: true,
        max_af_pop: Some("AMR"),
    },
    AfColumn {
        cache_col: "EAS",
        format_4f: false,
        flag_group: 1,
        emit_in_csq: true,
        max_af_pop: Some("EAS"),
    },
    AfColumn {
        cache_col: "EUR",
        format_4f: false,
        flag_group: 1,
        emit_in_csq: true,
        max_af_pop: Some("EUR"),
    },
    AfColumn {
        cache_col: "SAS",
        format_4f: false,
        flag_group: 1,
        emit_in_csq: true,
        max_af_pop: Some("SAS"),
    },
    // --af_gnomade — only global emitted in CSQ; sub-pops used for MAX_AF only
    AfColumn {
        cache_col: "gnomADe",
        format_4f: false,
        flag_group: 2,
        emit_in_csq: true,
        max_af_pop: None,
    },
    AfColumn {
        cache_col: "gnomADe_AFR",
        format_4f: false,
        flag_group: 2,
        emit_in_csq: false,
        max_af_pop: Some("gnomADe_AFR"),
    },
    AfColumn {
        cache_col: "gnomADe_AMR",
        format_4f: false,
        flag_group: 2,
        emit_in_csq: false,
        max_af_pop: Some("gnomADe_AMR"),
    },
    AfColumn {
        cache_col: "gnomADe_ASJ",
        format_4f: false,
        flag_group: 2,
        emit_in_csq: false,
        max_af_pop: Some("gnomADe_ASJ"),
    },
    AfColumn {
        cache_col: "gnomADe_EAS",
        format_4f: false,
        flag_group: 2,
        emit_in_csq: false,
        max_af_pop: Some("gnomADe_EAS"),
    },
    AfColumn {
        cache_col: "gnomADe_FIN",
        format_4f: false,
        flag_group: 2,
        emit_in_csq: false,
        max_af_pop: Some("gnomADe_FIN"),
    },
    AfColumn {
        cache_col: "gnomADe_MID",
        format_4f: false,
        flag_group: 2,
        emit_in_csq: false,
        max_af_pop: Some("gnomADe_MID"),
    },
    AfColumn {
        cache_col: "gnomADe_NFE",
        format_4f: false,
        flag_group: 2,
        emit_in_csq: false,
        max_af_pop: Some("gnomADe_NFE"),
    },
    AfColumn {
        cache_col: "gnomADe_REMAINING",
        format_4f: false,
        flag_group: 2,
        emit_in_csq: false,
        max_af_pop: Some("gnomADe_REMAINING"),
    },
    AfColumn {
        cache_col: "gnomADe_SAS",
        format_4f: false,
        flag_group: 2,
        emit_in_csq: false,
        max_af_pop: Some("gnomADe_SAS"),
    },
    // --af_gnomadg — only global emitted in CSQ; sub-pops used for MAX_AF only
    AfColumn {
        cache_col: "gnomADg",
        format_4f: false,
        flag_group: 3,
        emit_in_csq: true,
        max_af_pop: None,
    },
    AfColumn {
        cache_col: "gnomADg_AFR",
        format_4f: false,
        flag_group: 3,
        emit_in_csq: false,
        max_af_pop: Some("gnomADg_AFR"),
    },
    AfColumn {
        cache_col: "gnomADg_AMI",
        format_4f: false,
        flag_group: 3,
        emit_in_csq: false,
        max_af_pop: Some("gnomADg_AMI"),
    },
    AfColumn {
        cache_col: "gnomADg_AMR",
        format_4f: false,
        flag_group: 3,
        emit_in_csq: false,
        max_af_pop: Some("gnomADg_AMR"),
    },
    AfColumn {
        cache_col: "gnomADg_ASJ",
        format_4f: false,
        flag_group: 3,
        emit_in_csq: false,
        max_af_pop: Some("gnomADg_ASJ"),
    },
    AfColumn {
        cache_col: "gnomADg_EAS",
        format_4f: false,
        flag_group: 3,
        emit_in_csq: false,
        max_af_pop: Some("gnomADg_EAS"),
    },
    AfColumn {
        cache_col: "gnomADg_FIN",
        format_4f: false,
        flag_group: 3,
        emit_in_csq: false,
        max_af_pop: Some("gnomADg_FIN"),
    },
    AfColumn {
        cache_col: "gnomADg_MID",
        format_4f: false,
        flag_group: 3,
        emit_in_csq: false,
        max_af_pop: Some("gnomADg_MID"),
    },
    AfColumn {
        cache_col: "gnomADg_NFE",
        format_4f: false,
        flag_group: 3,
        emit_in_csq: false,
        max_af_pop: Some("gnomADg_NFE"),
    },
    AfColumn {
        cache_col: "gnomADg_REMAINING",
        format_4f: false,
        flag_group: 3,
        emit_in_csq: false,
        max_af_pop: Some("gnomADg_REMAINING"),
    },
    AfColumn {
        cache_col: "gnomADg_SAS",
        format_4f: false,
        flag_group: 3,
        emit_in_csq: false,
        max_af_pop: Some("gnomADg_SAS"),
    },
];

/// Parsed VEP option flags controlling which Batch 3 CSQ fields are emitted.
///
/// Flag names match Ensembl VEP CLI: `--check_existing`, `--af`, `--af_1kg`,
/// `--af_gnomade`, `--af_gnomadg`, `--max_af`, `--pubmed`, `--everything`.
///
/// Traceability:
/// - Ensembl VEP `Config.pm` `--everything` expansion
/// Cached parquet metadata for direct sift/polyphen window reads, bypassing DataFusion SQL.
struct SiftDirectReader {
    path: String,
    arrow_meta: parquet::arrow::arrow_reader::ArrowReaderMetadata,
    projection: parquet::arrow::ProjectionMask,
    rg_ranges: Vec<(i64, i64)>,
}

impl SiftDirectReader {
    /// Read a single genomic window directly from parquet, skipping DataFusion SQL planning.
    fn load_window(
        &self,
        chrom: &str,
        win_start: i64,
        win_end: i64,
        cache: &mut SiftPolyphenCache,
    ) -> Result<()> {
        use parquet::arrow::arrow_reader::ParquetRecordBatchReaderBuilder;

        let matching_rgs: Vec<usize> = self
            .rg_ranges
            .iter()
            .enumerate()
            .filter(|(_, (s, e))| *s <= win_end && *e >= win_start)
            .map(|(i, _)| i)
            .collect();

        if matching_rgs.is_empty() {
            return Ok(());
        }

        let file = std::fs::File::open(&self.path).map_err(|e| {
            DataFusionError::Execution(format!("failed to open sift file '{}': {e}", self.path))
        })?;
        let reader =
            ParquetRecordBatchReaderBuilder::new_with_metadata(file, self.arrow_meta.clone())
                .with_projection(self.projection.clone())
                .with_row_groups(matching_rgs)
                .build()
                .map_err(|e| {
                    DataFusionError::Execution(format!("failed to build sift reader: {e}"))
                })?;

        let chrom_norm = chrom.strip_prefix("chr").unwrap_or(chrom);

        for batch_result in reader {
            let batch = batch_result
                .map_err(|e| DataFusionError::Execution(format!("sift batch read error: {e}")))?;
            let schema = batch.schema();
            let tx_idx = schema
                .index_of("transcript_id")
                .or_else(|_| schema.index_of("stable_id"))
                .ok();
            let end_col_idx = schema.index_of("end").ok();
            let chrom_col_idx = schema.index_of("chrom").ok();
            let sift_col_idx = schema.index_of("sift_predictions").ok();
            let pp_col_idx = schema.index_of("polyphen_predictions").ok();

            let Some(tx_idx) = tx_idx else { continue };

            for row in 0..batch.num_rows() {
                // Filter by chromosome (the file may contain multiple chroms).
                if let Some(ci) = chrom_col_idx {
                    if let Some(row_chrom) = string_at(batch.column(ci).as_ref(), row) {
                        let row_norm = row_chrom.strip_prefix("chr").unwrap_or(&row_chrom);
                        if row_norm != chrom_norm {
                            continue;
                        }
                    }
                }
                let Some(tx_id) = string_at(batch.column(tx_idx).as_ref(), row) else {
                    continue;
                };
                if cache.get(&tx_id).is_some() {
                    continue;
                }
                let genomic_end = end_col_idx
                    .and_then(|idx| int64_at(batch.column(idx).as_ref(), row))
                    .unwrap_or(i64::MAX);
                let mut preds = CachedPredictions::default();
                if let Some(idx) = sift_col_idx {
                    preds.sift = read_compact_predictions(batch.column(idx).as_ref(), row);
                }
                if let Some(idx) = pp_col_idx {
                    preds.polyphen = read_compact_predictions(batch.column(idx).as_ref(), row);
                }
                preds.sort();
                cache.insert(tx_id, preds, genomic_end);
            }
        }
        Ok(())
    }
}

///   <https://github.com/Ensembl/ensembl-vep/blob/release/115/modules/Bio/EnsEMBL/VEP/Config.pm#L346-L374>
#[derive(Debug, Clone)]
struct VepFlags {
    check_existing: bool,
    af: bool,
    af_1kg: bool,
    af_gnomade: bool,
    af_gnomadg: bool,
    max_af: bool,
    pubmed: bool,
    /// When true, all VEP features are enabled and 80-field CSQ schema is used.
    everything: bool,
}

impl VepFlags {
    fn from_options_json(options_json: Option<&str>) -> Self {
        let parse = |key| {
            options_json
                .and_then(|opts| AnnotateProvider::parse_json_bool_option(opts, key))
                .unwrap_or(false)
        };
        let everything = parse("everything");
        // --everything implies all sub-flags per Config.pm#L346-L374.
        let af = everything || parse("af");
        let af_1kg = everything || parse("af_1kg");
        let af_gnomade = everything || parse("af_gnomade");
        let af_gnomadg = everything || parse("af_gnomadg");
        let max_af = everything || parse("max_af");
        let pubmed = everything || parse("pubmed");
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
            everything,
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

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
enum PickCriterion {
    ManeSelect,
    ManePlusClinical,
    Canonical,
    Appris,
    Tsl,
    Biotype,
    Ccds,
    Rank,
    Length,
    Ensembl,
    Refseq,
}

impl PickCriterion {
    fn from_str(raw: &str) -> Result<Self> {
        match raw.trim() {
            "mane_select" => Ok(Self::ManeSelect),
            "mane_plus_clinical" => Ok(Self::ManePlusClinical),
            "canonical" => Ok(Self::Canonical),
            "appris" => Ok(Self::Appris),
            "tsl" => Ok(Self::Tsl),
            "biotype" => Ok(Self::Biotype),
            "ccds" => Ok(Self::Ccds),
            "rank" => Ok(Self::Rank),
            "length" => Ok(Self::Length),
            "ensembl" => Ok(Self::Ensembl),
            "refseq" => Ok(Self::Refseq),
            other => Err(DataFusionError::Execution(format!(
                "annotate_vep(): unsupported pick_order criterion '{other}'. Supported criteria: mane_select, mane_plus_clinical, canonical, appris, tsl, biotype, ccds, rank, length, ensembl, refseq"
            ))),
        }
    }
}

#[derive(Debug, Clone)]
struct PickFlags {
    flag_pick_allele_gene: bool,
    pick_order: Vec<PickCriterion>,
}

impl Default for PickFlags {
    fn default() -> Self {
        Self {
            flag_pick_allele_gene: false,
            // Traceability:
            // - Ensembl VEP release 115 default `pick_order`
            //   <https://github.com/Ensembl/ensembl-vep/blob/release/115/modules/Bio/EnsEMBL/VEP/Config.pm#L300-L307>
            pick_order: vec![
                PickCriterion::ManeSelect,
                PickCriterion::ManePlusClinical,
                PickCriterion::Canonical,
                PickCriterion::Appris,
                PickCriterion::Tsl,
                PickCriterion::Biotype,
                PickCriterion::Ccds,
                PickCriterion::Rank,
                PickCriterion::Length,
                PickCriterion::Ensembl,
                PickCriterion::Refseq,
            ],
        }
    }
}

impl PickFlags {
    fn from_options_json(options_json: Option<&str>) -> Result<Self> {
        let parse = |key| {
            options_json
                .and_then(|opts| AnnotateProvider::parse_json_bool_option(opts, key))
                .unwrap_or(false)
        };

        let mut out = Self {
            flag_pick_allele_gene: parse("flag_pick_allele_gene"),
            ..Self::default()
        };

        if let Some(raw_order) = options_json
            .and_then(|opts| AnnotateProvider::parse_json_string_option(opts, "pick_order"))
        {
            let parsed: Vec<PickCriterion> = raw_order
                .split(',')
                .map(PickCriterion::from_str)
                .collect::<Result<Vec<_>>>()?;
            if parsed.is_empty() {
                return Err(DataFusionError::Execution(
                    "annotate_vep(): pick_order must contain at least one criterion".to_string(),
                ));
            }
            out.pick_order = parsed;
        }

        Ok(out)
    }

    fn requires_transcript_annotations(&self, skip_csq: bool, skip_typed_cols: bool) -> bool {
        self.flag_pick_allele_gene && (!skip_csq || !skip_typed_cols)
    }
}

/// Parsed HGVS-related flags controlling HGVSc/HGVSp emission.
///
/// Traceability:
/// - Ensembl VEP `Config.pm` HGVS-related flags
///   <https://github.com/Ensembl/ensembl-vep/blob/release/115/modules/Bio/EnsEMBL/VEP/Config.pm#L195-L200>
/// - Ensembl VEP `Config.pm` `shift_hgvs`
///   <https://github.com/Ensembl/ensembl-vep/blob/release/115/modules/Bio/EnsEMBL/VEP/Config.pm#L353-L381>
/// - Ensembl VEP `Runner::post_setup_checks()`
///   <https://github.com/Ensembl/ensembl-vep/blob/release/115/modules/Bio/EnsEMBL/VEP/Runner.pm#L771-L773>
/// - Ensembl VEP `OutputFactory::TranscriptVariationAllele_to_output_hash()`
///   <https://github.com/Ensembl/ensembl-vep/blob/release/115/modules/Bio/EnsEMBL/VEP/OutputFactory.pm#L1698-L1715>
#[derive(Debug, Clone, Copy, Default)]
struct HgvsFlags {
    hgvsc: bool,
    hgvsp: bool,
    shift_hgvs: bool,
    no_escape: bool,
    remove_hgvsp_version: bool,
    hgvsp_use_prediction: bool,
}

impl HgvsFlags {
    fn from_options_json(options_json: Option<&str>) -> Self {
        let parse = |key| {
            options_json
                .and_then(|opts| AnnotateProvider::parse_json_bool_option(opts, key))
                .unwrap_or(false)
        };
        let everything = parse("everything");
        // --everything implies --hgvs per Config.pm#L346-L374.
        let hgvs = everything || parse("hgvs");
        let hgvsc = hgvs || parse("hgvsc");
        let hgvsp = hgvs || parse("hgvsp");
        let shift_hgvs = options_json
            .and_then(|opts| AnnotateProvider::parse_json_bool_option(opts, "shift_hgvs"))
            .unwrap_or(hgvsc || hgvsp);
        Self {
            hgvsc,
            hgvsp,
            shift_hgvs,
            no_escape: parse("no_escape"),
            remove_hgvsp_version: parse("remove_hgvsp_version"),
            hgvsp_use_prediction: parse("hgvsp_use_prediction"),
        }
    }

    fn any(self) -> bool {
        self.hgvsc || self.hgvsp
    }
}

#[derive(Debug, Clone, Copy, PartialEq, Eq, Default)]
enum TranscriptSourceMode {
    #[default]
    Ensembl,
    RefSeq,
    Merged,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq, Default)]
struct TranscriptSelectionFlags {
    source_mode: TranscriptSourceMode,
    gencode_basic: bool,
    gencode_primary: bool,
    all_refseq: bool,
    exclude_predicted: bool,
}

impl TranscriptSelectionFlags {
    fn from_options_json(options_json: Option<&str>) -> Result<Self> {
        let parse = |key| {
            options_json
                .and_then(|opts| AnnotateProvider::parse_json_bool_option(opts, key))
                .unwrap_or(false)
        };

        let refseq = parse("refseq");
        let merged = parse("merged");
        let gencode_basic = parse("gencode_basic");
        let gencode_primary = parse("gencode_primary");
        let all_refseq = parse("all_refseq");
        let exclude_predicted = parse("exclude_predicted");

        if refseq && merged {
            return Err(DataFusionError::Execution(
                "annotate_vep(): --refseq and --merged are mutually exclusive".to_string(),
            ));
        }
        if refseq && gencode_basic {
            return Err(DataFusionError::Execution(
                "annotate_vep(): --refseq and --gencode_basic are mutually exclusive".to_string(),
            ));
        }
        if refseq && gencode_primary {
            return Err(DataFusionError::Execution(
                "annotate_vep(): --refseq and --gencode_primary are mutually exclusive".to_string(),
            ));
        }
        if gencode_basic && gencode_primary {
            return Err(DataFusionError::Execution(
                "annotate_vep(): --gencode_basic and --gencode_primary are mutually exclusive"
                    .to_string(),
            ));
        }

        let source_mode = if merged {
            TranscriptSourceMode::Merged
        } else if refseq {
            TranscriptSourceMode::RefSeq
        } else {
            TranscriptSourceMode::Ensembl
        };

        if source_mode == TranscriptSourceMode::Ensembl && all_refseq {
            return Err(DataFusionError::Execution(
                "annotate_vep(): --all_refseq requires --refseq or --merged".to_string(),
            ));
        }
        if source_mode == TranscriptSourceMode::Ensembl && exclude_predicted {
            return Err(DataFusionError::Execution(
                "annotate_vep(): --exclude_predicted requires --refseq or --merged".to_string(),
            ));
        }

        Ok(Self {
            source_mode,
            gencode_basic,
            gencode_primary,
            all_refseq,
            exclude_predicted,
        })
    }

    fn merged(self) -> bool {
        self.source_mode == TranscriptSourceMode::Merged
    }

    fn refseq_fields(self) -> bool {
        matches!(
            self.source_mode,
            TranscriptSourceMode::RefSeq | TranscriptSourceMode::Merged
        )
    }

    fn source_field(self) -> bool {
        self.source_mode == TranscriptSourceMode::Merged
    }
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
enum CsqPlaceholderField {
    Empty,
    Allele,
    Consequence,
    Impact,
    ExistingVariation,
    VariantClass,
    AfValue(usize),
    MaxAf,
    MaxAfPops,
    ClinSig,
    Somatic,
    Pheno,
    Pubmed,
}

impl CsqPlaceholderField {
    fn from_name(name: &str) -> Self {
        match name {
            "Allele" => Self::Allele,
            "Consequence" => Self::Consequence,
            "IMPACT" => Self::Impact,
            "Existing_variation" => Self::ExistingVariation,
            "VARIANT_CLASS" => Self::VariantClass,
            "AF" => Self::AfValue(0),
            "AFR_AF" => Self::AfValue(1),
            "AMR_AF" => Self::AfValue(2),
            "EAS_AF" => Self::AfValue(3),
            "EUR_AF" => Self::AfValue(4),
            "SAS_AF" => Self::AfValue(5),
            "gnomADe_AF" => Self::AfValue(6),
            "gnomADe_AFR" | "gnomADe_AFR_AF" => Self::AfValue(7),
            "gnomADe_AMR" | "gnomADe_AMR_AF" => Self::AfValue(8),
            "gnomADe_ASJ" | "gnomADe_ASJ_AF" => Self::AfValue(9),
            "gnomADe_EAS" | "gnomADe_EAS_AF" => Self::AfValue(10),
            "gnomADe_FIN" | "gnomADe_FIN_AF" => Self::AfValue(11),
            "gnomADe_MID" | "gnomADe_MID_AF" => Self::AfValue(12),
            "gnomADe_NFE" | "gnomADe_NFE_AF" => Self::AfValue(13),
            "gnomADe_REMAINING" | "gnomADe_REMAINING_AF" => Self::AfValue(14),
            "gnomADe_SAS" | "gnomADe_SAS_AF" => Self::AfValue(15),
            "gnomADg_AF" => Self::AfValue(16),
            "gnomADg_AFR" | "gnomADg_AFR_AF" => Self::AfValue(17),
            "gnomADg_AMI" | "gnomADg_AMI_AF" => Self::AfValue(18),
            "gnomADg_AMR" | "gnomADg_AMR_AF" => Self::AfValue(19),
            "gnomADg_ASJ" | "gnomADg_ASJ_AF" => Self::AfValue(20),
            "gnomADg_EAS" | "gnomADg_EAS_AF" => Self::AfValue(21),
            "gnomADg_FIN" | "gnomADg_FIN_AF" => Self::AfValue(22),
            "gnomADg_MID" | "gnomADg_MID_AF" => Self::AfValue(23),
            "gnomADg_NFE" | "gnomADg_NFE_AF" => Self::AfValue(24),
            "gnomADg_REMAINING" | "gnomADg_REMAINING_AF" => Self::AfValue(25),
            "gnomADg_SAS" | "gnomADg_SAS_AF" => Self::AfValue(26),
            "MAX_AF" => Self::MaxAf,
            "MAX_AF_POPS" => Self::MaxAfPops,
            "CLIN_SIG" => Self::ClinSig,
            "SOMATIC" => Self::Somatic,
            "PHENO" => Self::Pheno,
            "PUBMED" => Self::Pubmed,
            _ => Self::Empty,
        }
    }

    fn value<'a>(self, entry: &'a CsqPlaceholderEntry<'a>) -> &'a str {
        match self {
            Self::Empty => "",
            Self::Allele => entry.allele,
            Self::Consequence => entry.consequence,
            Self::Impact => entry.impact,
            Self::ExistingVariation => entry.existing_variation,
            Self::VariantClass => entry.variant_class,
            Self::AfValue(idx) => entry
                .frequency_fields
                .af_values
                .get(idx)
                .map(String::as_str)
                .unwrap_or(""),
            Self::MaxAf => entry.frequency_fields.max_af.as_str(),
            Self::MaxAfPops => entry.frequency_fields.max_af_pops.as_str(),
            Self::ClinSig => entry.variant_fields.clin_sig.as_str(),
            Self::Somatic => entry.variant_fields.somatic.as_str(),
            Self::Pheno => entry.variant_fields.pheno.as_str(),
            Self::Pubmed => entry.variant_fields.pubmed.as_str(),
        }
    }
}

#[derive(Debug, Clone)]
struct CsqPlaceholderLayout {
    fields: Vec<CsqPlaceholderField>,
}

impl CsqPlaceholderLayout {
    fn for_mode(everything: bool, transcript_selection: TranscriptSelectionFlags) -> Self {
        let fields = crate::golden_benchmark::csq_field_names_for_mode(
            everything,
            transcript_selection.source_mode == TranscriptSourceMode::RefSeq,
            transcript_selection.source_mode == TranscriptSourceMode::Merged,
        )
        .into_iter()
        .map(CsqPlaceholderField::from_name)
        .collect();
        Self { fields }
    }

    fn append_entry(&self, buf: &mut String, entry: &CsqPlaceholderEntry<'_>) {
        for (idx, field) in self.fields.iter().enumerate() {
            if idx > 0 {
                buf.push('|');
            }
            buf.push_str(field.value(entry));
        }
    }
}

#[derive(Debug)]
struct CsqPlaceholderEntry<'a> {
    allele: &'a str,
    consequence: &'a str,
    impact: &'a str,
    existing_variation: &'a str,
    variant_class: &'a str,
    frequency_fields: &'a ColocatedFrequencyFields,
    variant_fields: &'a ColocatedVariantFields,
}

#[derive(Debug, Default)]
struct TranscriptRawMetadata {
    display_xref_id: Option<String>,
    source: Option<String>,
    gene_hgnc_id_native: Option<String>,
    refseq_match: Option<String>,
    refseq_edits: Vec<RefSeqEdit>,
    cdna_mapper_segments: Vec<TranscriptCdnaMapperSegment>,
    spliced_seq: Option<String>,
    five_prime_utr_seq: Option<String>,
    three_prime_utr_seq: Option<String>,
    translateable_seq: Option<String>,
    flags_str: Option<String>,
    is_gencode_basic: bool,
    is_gencode_primary: bool,
}

/// A single co-located variant entry with allele and clinical metadata.
#[derive(Debug, Clone)]
struct ColocatedEntry {
    variation_name: String,
    allele_string: String,
    matched_alleles: Vec<MatchedVariantAllele>,
    somatic: i64,
    pheno: i64,
    clin_sig: Option<String>,
    clin_sig_allele: Option<String>,
    pubmed: Option<String>,
    /// Raw AF column values (same order as `AF_COL_NAMES` / `AF_COLUMNS`).
    af_values: Vec<String>,
}

#[derive(Debug, Default, Clone)]
struct ColocatedData {
    entries: Vec<ColocatedEntry>,
    compare_output_allele: Option<String>,
    unshifted_output_allele: Option<String>,
}

#[derive(Debug, Default)]
struct ColocatedVariantFields {
    existing_variation: String,
    clin_sig: String,
    somatic: String,
    pheno: String,
    pubmed: String,
}

#[derive(Debug, Default)]
struct ColocatedFrequencyFields {
    af_values: Vec<String>,
    max_af: String,
    max_af_pops: String,
}

fn variant_prefix_rank(variation_name: &str) -> u8 {
    match variation_name
        .get(..2)
        .unwrap_or(variation_name)
        .to_ascii_lowercase()
        .as_str()
    {
        "rs" => 1,
        "cm" | "ci" | "cd" => 2,
        "co" => 3,
        _ => 100,
    }
}

fn push_unique_value(values: &mut Vec<String>, value: impl Into<String>) {
    let value = value.into();
    if !values.iter().any(|existing| existing == &value) {
        values.push(value);
    }
}

impl ColocatedEntry {
    /// Traceability:
    /// - Ensembl VEP `add_colocated_variant_info()`
    ///   <https://github.com/Ensembl/ensembl-vep/blob/2beada0d57ca6234f467b14a6c60280f4a082717/modules/Bio/EnsEMBL/VEP/OutputFactory.pm#L1012-L1035>
    /// - Ensembl VEP `add_colocated_frequency_data()`
    ///   <https://github.com/Ensembl/ensembl-vep/blob/2beada0d57ca6234f467b14a6c60280f4a082717/modules/Bio/EnsEMBL/VEP/OutputFactory.pm#L1150-L1157>
    ///
    /// The colocated output path must resolve the matched existing-variant
    /// allele against the live CSQ allele, with optional fallback to the
    /// retained unshifted allele when VEP preserved original shift metadata.
    fn matching_allele<'a>(
        &'a self,
        output_allele: &str,
        output_allele_unshifted: Option<&str>,
    ) -> Option<&'a MatchedVariantAllele> {
        self.matched_alleles.iter().find(|matched| {
            matched.a_allele == output_allele
                || output_allele_unshifted.is_some_and(|allele| matched.a_allele == allele)
        })
    }

    /// Traceability:
    /// - Ensembl VEP `add_colocated_variant_info()`
    ///   <https://github.com/Ensembl/ensembl-vep/blob/2beada0d57ca6234f467b14a6c60280f4a082717/modules/Bio/EnsEMBL/VEP/OutputFactory.pm#L1012-L1035>
    ///
    /// Existing variants without `matched_alleles` remain visible, but once
    /// a matched-allele map exists the output must be filtered by the active
    /// CSQ allele exactly as OutputFactory does.
    fn matches_output_allele(
        &self,
        output_allele: &str,
        output_allele_unshifted: Option<&str>,
    ) -> bool {
        self.matched_alleles.is_empty()
            || self
                .matching_allele(output_allele, output_allele_unshifted)
                .is_some()
    }
}

impl ColocatedData {
    /// Traceability:
    /// - Ensembl VEP `add_colocated_variant_info()`
    ///   <https://github.com/Ensembl/ensembl-vep/blob/2beada0d57ca6234f467b14a6c60280f4a082717/modules/Bio/EnsEMBL/VEP/OutputFactory.pm#L1012-L1035>
    ///
    /// Rust stores the active compare-space allele and any retained original
    /// compare-space allele separately on the colocated sink. For the live CSQ
    /// path, `Existing_variation` must prefer the active compare-space allele
    /// and only fall back to the retained original allele when the output
    /// allele already equals the active representation.
    fn variant_match_output_allele<'a>(&'a self, output_allele: &str) -> Option<&'a str> {
        self.compare_output_allele
            .as_deref()
            .filter(|allele| *allele != output_allele)
            .or_else(|| {
                self.unshifted_output_allele
                    .as_deref()
                    .filter(|allele| *allele != output_allele)
            })
    }

    /// Traceability:
    /// - Ensembl VEP `add_colocated_frequency_data()`
    ///   <https://github.com/Ensembl/ensembl-vep/blob/2beada0d57ca6234f467b14a6c60280f4a082717/modules/Bio/EnsEMBL/VEP/OutputFactory.pm#L1150-L1157>
    ///
    /// Frequency output matches the current CSQ allele first, then
    /// `alt_orig_allele_string` when VEP retained original shift metadata.
    /// Because Rust stores both active and retained compare-space alleles on
    /// the sink, the live path must prefer the retained original allele here.
    fn frequency_match_output_allele<'a>(&'a self, output_allele: &str) -> Option<&'a str> {
        self.unshifted_output_allele
            .as_deref()
            .filter(|allele| *allele != output_allele)
            .or_else(|| {
                self.compare_output_allele
                    .as_deref()
                    .filter(|allele| *allele != output_allele)
            })
    }

    /// Traceability:
    /// - Ensembl VEP `add_colocated_variant_info()`
    ///   <https://github.com/Ensembl/ensembl-vep/blob/2beada0d57ca6234f467b14a6c60280f4a082717/modules/Bio/EnsEMBL/VEP/OutputFactory.pm#L1005-L1011>
    ///
    /// OutputFactory sorts co-located variants by somatic status and then by
    /// variation-name class before emitting IDs and metadata.
    fn sorted_entries(&self) -> Vec<&ColocatedEntry> {
        let mut entries: Vec<&ColocatedEntry> = self.entries.iter().collect();
        entries.sort_by(|a, b| {
            (a.somatic != 0).cmp(&(b.somatic != 0)).then_with(|| {
                variant_prefix_rank(&a.variation_name).cmp(&variant_prefix_rank(&b.variation_name))
            })
        });
        entries
    }

    /// Traceability:
    /// - Ensembl VEP `add_colocated_variant_info()`
    ///   <https://github.com/Ensembl/ensembl-vep/blob/2beada0d57ca6234f467b14a6c60280f4a082717/modules/Bio/EnsEMBL/VEP/OutputFactory.pm#L1005-L1120>
    ///
    /// This mirrors OutputFactory's co-located ID and clinical-field assembly:
    /// sort existing variants the same way, filter by `matched_alleles` using
    /// the current output allele plus `alt_orig_allele_string` only when
    /// upstream shift state exists, and preserve the `clin_sig` vs
    /// `clin_sig_allele` split.
    fn variant_fields(
        &self,
        output_allele: &str,
        output_allele_unshifted: Option<&str>,
        include_pubmed: bool,
    ) -> ColocatedVariantFields {
        let mut fields = ColocatedVariantFields::default();
        let mut clin_sig_values: Vec<String> = Vec::new();
        let mut clin_sig_allele_values: Vec<String> = Vec::new();
        let mut clin_sig_allele_exists = false;
        let mut somatic_values: Vec<&str> = Vec::new();
        let mut pheno_values: Vec<&str> = Vec::new();
        let mut pubmed_values: Vec<String> = Vec::new();

        for entry in self.sorted_entries() {
            if !entry.matches_output_allele(output_allele, output_allele_unshifted) {
                continue;
            }

            if !entry.variation_name.is_empty() {
                if !fields.existing_variation.is_empty() {
                    fields.existing_variation.push('&');
                }
                fields.existing_variation.push_str(&entry.variation_name);
            }

            if let Some(clin_sig_allele) = &entry.clin_sig_allele {
                let mut allele_terms: HashMap<String, String> = HashMap::new();
                for chunk in clin_sig_allele.split(';') {
                    let Some((allele, value)) = chunk.split_once(':') else {
                        continue;
                    };
                    let slot = allele_terms.entry(allele.to_string()).or_default();
                    if !slot.is_empty() {
                        slot.push(',');
                    }
                    slot.push_str(value);
                }
                if let Some(value) = allele_terms.get(output_allele) {
                    push_unique_value(&mut clin_sig_allele_values, value.clone());
                }
                clin_sig_allele_exists = true;
            }

            if let Some(clin_sig) = &entry.clin_sig {
                if !clin_sig_allele_exists {
                    for value in clin_sig.split(',') {
                        if !value.is_empty() {
                            clin_sig_values.push(value.to_string());
                        }
                    }
                }
            }

            somatic_values.push(if entry.somatic != 0 { "1" } else { "0" });
            pheno_values.push(if entry.pheno != 0 { "1" } else { "0" });

            if include_pubmed {
                if let Some(pubmed) = &entry.pubmed {
                    for value in pubmed.split(',') {
                        if !value.is_empty() {
                            pubmed_values.push(value.to_string());
                        }
                    }
                }
            }
        }

        if somatic_values.iter().any(|value| *value == "1") {
            fields.somatic = somatic_values.join("&");
        }
        if pheno_values.iter().any(|value| *value == "1") {
            fields.pheno = pheno_values.join("&");
        }
        if !pubmed_values.is_empty() {
            fields.pubmed = csq_escape(&pubmed_values.join("&")).into_owned();
        }
        if !clin_sig_allele_values.is_empty() {
            fields.clin_sig = csq_escape(&clin_sig_allele_values.join(";")).into_owned();
        } else if !clin_sig_values.is_empty() {
            fields.clin_sig = csq_escape(&clin_sig_values.join("&")).into_owned();
        }

        fields
    }

    /// Traceability:
    /// - Ensembl VEP `add_colocated_frequency_data()`
    ///   <https://github.com/Ensembl/ensembl-vep/blob/2beada0d57ca6234f467b14a6c60280f4a082717/modules/Bio/EnsEMBL/VEP/OutputFactory.pm#L1139-L1232>
    /// - Ensembl VEP `get_frequency_data()`
    ///   <https://github.com/Ensembl/ensembl-vep/blob/2beada0d57ca6234f467b14a6c60280f4a082717/modules/Bio/EnsEMBL/VEP/AnnotationSource/Cache/BaseCacheVariation.pm#L179-L255>
    ///
    /// This mirrors VEP's per-existing-variant frequency projection: build
    /// allele-to-frequency maps from the matched existing variant, select the
    /// `b_allele` named in `matched_alleles`, use `alt_orig_allele_string`
    /// only when present on the output allele, and only interpolate the global
    /// `AF` field in the same biallelic case VEP allows.
    fn frequency_fields(
        &self,
        output_allele: &str,
        output_allele_unshifted: Option<&str>,
        flags: &VepFlags,
    ) -> ColocatedFrequencyFields {
        let mut per_column: Vec<Vec<String>> = vec![Vec::new(); AF_COLUMNS.len()];
        let mut max_af: Option<(f64, String)> = None;
        let mut max_af_pops: Vec<String> = Vec::new();

        for entry in self.sorted_entries() {
            let Some(matched_allele) =
                entry.matching_allele(output_allele, output_allele_unshifted)
            else {
                continue;
            };

            let existing_alleles: Vec<&str> = entry.allele_string.split('/').collect();
            let mut entry_max_af: Option<(f64, String)> = None;
            let mut entry_max_af_pops: Vec<String> = Vec::new();

            for (idx, column) in AF_COLUMNS.iter().enumerate() {
                let should_process = flags.max_af || flags.af_group_enabled(column.flag_group);
                if !should_process || idx >= entry.af_values.len() {
                    continue;
                }

                let raw = &entry.af_values[idx];
                if raw.is_empty() {
                    continue;
                }

                let mut freq_data: HashMap<String, String> = HashMap::new();
                let mut remaining: HashSet<String> = existing_alleles
                    .iter()
                    .map(|allele| (*allele).to_string())
                    .collect();
                let mut total = 0.0_f64;

                for pair in raw.split(',') {
                    let Some((allele, freq)) = pair.split_once(':') else {
                        continue;
                    };
                    let formatted = if column.format_4f {
                        format_af_4f(freq)
                    } else {
                        freq.to_string()
                    };
                    freq_data.insert(allele.to_string(), formatted);
                    total += freq.parse::<f64>().unwrap_or(0.0);
                    remaining.remove(allele);
                }

                let mut interpolated = false;
                if existing_alleles.len() == 2 && remaining.len() == 1 && column.cache_col == "AF" {
                    let remaining_allele = remaining.into_iter().next().unwrap();
                    freq_data.insert(remaining_allele, format!("{}", 1.0 - total));
                    interpolated = true;
                }

                let chosen = if let Some(value) = freq_data.get(&matched_allele.b_allele) {
                    Some(value.clone())
                } else if interpolated {
                    freq_data.get(output_allele).cloned()
                } else {
                    None
                };
                let Some(chosen) = chosen else {
                    continue;
                };

                if flags.af_group_enabled(column.flag_group) {
                    push_unique_value(&mut per_column[idx], chosen.clone());
                }

                if flags.max_af {
                    if let Some(pop_name) = column.max_af_pop {
                        if let Ok(freq) = chosen.parse::<f64>() {
                            match entry_max_af {
                                None => {
                                    entry_max_af = Some((freq, chosen.clone()));
                                    entry_max_af_pops.clear();
                                    entry_max_af_pops.push(pop_name.to_string());
                                }
                                Some((current, _)) if freq > current => {
                                    entry_max_af = Some((freq, chosen.clone()));
                                    entry_max_af_pops.clear();
                                    entry_max_af_pops.push(pop_name.to_string());
                                }
                                Some((current, _)) if (freq - current).abs() < f64::EPSILON => {
                                    push_unique_value(&mut entry_max_af_pops, pop_name.to_string());
                                }
                                _ => {}
                            }
                        }
                    }
                }
            }

            if flags.max_af && !entry_max_af_pops.is_empty() {
                let (entry_max, entry_max_str) = entry_max_af.unwrap_or((0.0, String::new()));
                let current_max = max_af.as_ref().map(|(value, _)| *value).unwrap_or(0.0);

                if entry_max > current_max {
                    max_af = Some((entry_max, entry_max_str.clone()));
                    max_af_pops.clear();
                }

                if entry_max >= current_max {
                    if max_af.is_none() {
                        max_af = Some((entry_max, entry_max_str));
                    }
                    max_af_pops.extend(entry_max_af_pops);
                }
            }
        }

        let af_values = AF_COLUMNS
            .iter()
            .enumerate()
            .map(|(idx, column)| {
                if column.emit_in_csq || flags.everything {
                    per_column[idx].join("&")
                } else {
                    String::new()
                }
            })
            .collect();

        ColocatedFrequencyFields {
            af_values,
            max_af: max_af.map(|(_, raw)| raw).unwrap_or_default(),
            max_af_pops: max_af_pops.join("&"),
        }
    }
}

/// Build co-located variant aggregation from the piggybacked collection sink.
///
/// Converts `ColocatedCacheEntry` entries (collected during `VariantLookupExec`
/// probe phase) into the same `ColocatedData` format used by the CSQ assembler,
/// preserving the per-input-allele separation VEP keeps between different
/// parser/decomposed alleles at the same locus.
///
/// Traceability:
/// - Ensembl VEP `compare_existing()`
///   <https://github.com/Ensembl/ensembl-vep/blob/2beada0d57ca6234f467b14a6c60280f4a082717/modules/Bio/EnsEMBL/VEP/AnnotationType/Variation.pm#L199-L206>
/// - Ensembl VEP `add_colocated_variant_info()`
///   <https://github.com/Ensembl/ensembl-vep/blob/2beada0d57ca6234f467b14a6c60280f4a082717/modules/Bio/EnsEMBL/VEP/OutputFactory.pm#L1032-L1049>
///
/// Each existing variant must retain the `matched_alleles` attached during
/// `compare_existing()`, merged across duplicate probe hits but without
/// re-synthesizing allele matches from local heuristics. The active compare
/// allele plus any retained original compare allele must also survive this
/// merge so the live CSQ path can reproduce OutputFactory's shifted-vs-original
/// allele filtering.
fn build_colocated_map_from_sink(
    sink: &HashMap<ColocatedKey, ColocatedSinkValue>,
) -> HashMap<ColocatedKey, ColocatedData> {
    let mut map = HashMap::with_capacity(sink.len());
    for (key, sink_value) in sink {
        let cache_entries = &sink_value.entries;
        let mut entries: Vec<ColocatedEntry> = Vec::new();
        let mut seen: HashMap<String, usize> = HashMap::new();
        for ce in cache_entries {
            if ce.variation_name.is_empty() {
                continue;
            }
            if let Some(existing_idx) = seen.get(&ce.variation_name) {
                let existing = &mut entries[*existing_idx];
                for matched in &ce.matched_alleles {
                    if !existing
                        .matched_alleles
                        .iter()
                        .any(|entry| entry == matched)
                    {
                        existing.matched_alleles.push(matched.clone());
                    }
                }
                continue;
            }
            seen.insert(ce.variation_name.clone(), entries.len());
            entries.push(ColocatedEntry {
                variation_name: ce.variation_name.clone(),
                allele_string: ce.allele_string.clone(),
                matched_alleles: ce.matched_alleles.clone(),
                somatic: ce.somatic,
                pheno: ce.pheno,
                clin_sig: ce.clin_sig.clone(),
                clin_sig_allele: ce.clin_sig_allele.clone(),
                pubmed: ce.pubmed.clone(),
                af_values: ce.af_values.clone(),
            });
        }

        map.insert(
            key.clone(),
            ColocatedData {
                entries,
                compare_output_allele: sink_value.compare_output_allele.clone(),
                unshifted_output_allele: sink_value.unshifted_output_allele.clone(),
            },
        );
    }
    map
}

/// Traceability:
/// - Ensembl VEP `output_hash_to_vcf_info_chunk()`
///   <https://github.com/Ensembl/ensembl-vep/blob/2beada0d57ca6234f467b14a6c60280f4a082717/modules/Bio/EnsEMBL/VEP/OutputFactory/VCF.pm#L379-L405>
///
/// This applies the same CSQ field escaping VEP uses for VCF output:
/// array separators become `&`, semicolons are percent-encoded, spaces become
/// underscores, pipes become `&`, and `-` is serialized as an empty field.
#[inline]
fn csq_escape(val: &str) -> std::borrow::Cow<'_, str> {
    if val == "-" {
        return std::borrow::Cow::Borrowed("");
    }

    let mut changed = false;
    let mut escaped = String::with_capacity(val.len());
    for ch in val.chars() {
        match ch {
            ',' | '|' => {
                escaped.push('&');
                changed = true;
            }
            ';' => {
                escaped.push_str("%3B");
                changed = true;
            }
            ch if ch.is_whitespace() => {
                escaped.push('_');
                changed = true;
            }
            _ => escaped.push(ch),
        }
    }

    if changed {
        std::borrow::Cow::Owned(escaped)
    } else {
        std::borrow::Cow::Borrowed(val)
    }
}

/// Format APPRIS annotation code for CSQ output.
///
/// Traceability:
/// - Ensembl VEP `OutputFactory::TranscriptVariationAllele_to_output_hash()`
///   <https://github.com/Ensembl/ensembl-vep/blob/release/115/modules/Bio/EnsEMBL/VEP/OutputFactory.pm#L1563-L1570>
///
/// VEP abbreviates: `principal1` → `P1`, `alternative2` → `A2`.
fn format_appris(raw: &str) -> String {
    raw.replace("principal", "P").replace("alternative", "A")
}

fn parse_appris_pick_rank(raw: Option<&str>) -> u8 {
    let Some(raw) = raw else {
        return 100;
    };

    // Traceability:
    // - Ensembl VEP release 115 `pick_worst_VariationFeatureOverlapAllele()`
    //   parses APPRIS as `([A-Za-z]).+(\\d+)`, then adds 10 for
    //   `alternative*` values
    //   <https://github.com/Ensembl/ensembl-vep/blob/release/115/modules/Bio/EnsEMBL/VEP/OutputFactory.pm#L756-L767>
    let Some(first) = raw
        .chars()
        .next()
        .filter(|ch| ch.is_ascii_alphabetic())
        .map(|ch| ch.to_ascii_lowercase())
    else {
        return 100;
    };
    let trailing_digits_len = raw
        .chars()
        .rev()
        .take_while(|ch| ch.is_ascii_digit())
        .count();
    if trailing_digits_len == 0 {
        return 100;
    }
    let digits = &raw[raw.len() - trailing_digits_len..];
    let Some(mut grade) = digits.parse::<u8>().ok() else {
        return 100;
    };
    if grade == 0 {
        return 100;
    }
    if first == 'a' {
        grade = grade.saturating_add(10);
    }
    grade
}

#[derive(Debug, Clone)]
struct PickCandidateInfo {
    assignment_idx: usize,
    mane_select: u8,
    mane_plus_clinical: u8,
    canonical: u8,
    ccds: u8,
    length: i64,
    biotype: u8,
    tsl: i32,
    appris: u8,
    ensembl: u8,
    refseq: u8,
    rank: u32,
}

fn build_pick_candidate_info(
    assignment_idx: usize,
    assignments: &[TranscriptConsequence],
    ctx: &PreparedContext<'_>,
) -> PickCandidateInfo {
    let tc = &assignments[assignment_idx];
    let tx = tc.transcript_idx.map(|idx| ctx.transcripts[idx]);

    let (mane_select, mane_plus_clinical, canonical, ccds, length, biotype, tsl, appris) =
        if let Some(tx) = tx {
            // Traceability:
            // - Ensembl VEP release 115 `pick_worst_VariationFeatureOverlapAllele()`
            //   initializes transcript criteria from transcript attributes,
            //   translation/CDS length, and `_source_cache`
            //   <https://github.com/Ensembl/ensembl-vep/blob/release/115/modules/Bio/EnsEMBL/VEP/OutputFactory.pm#L702-L767>
            let length =
                if let Some(translation) = ctx.translation_by_tx.get(tx.transcript_id.as_str()) {
                    let cds_len = translation
                        .cds_sequence
                        .as_ref()
                        .map(|sequence| sequence.len())
                        .or(translation.cds_len)
                        .unwrap_or(0);
                    -i64::try_from(cds_len).unwrap_or(i64::MAX)
                } else {
                    -tx.end.saturating_sub(tx.start).saturating_add(1)
                };
            (
                if tx
                    .mane_select
                    .as_deref()
                    .is_some_and(|value| !value.is_empty())
                {
                    0
                } else {
                    1
                },
                if tx
                    .mane_plus_clinical
                    .as_deref()
                    .is_some_and(|value| !value.is_empty())
                {
                    0
                } else {
                    1
                },
                if tx.is_canonical { 0 } else { 1 },
                if tx
                    .ccds
                    .as_deref()
                    .is_some_and(|value| !value.is_empty() && value != "-")
                {
                    0
                } else {
                    1
                },
                length,
                if tx.biotype == "protein_coding" { 0 } else { 1 },
                tx.tsl.unwrap_or(100),
                parse_appris_pick_rank(tx.appris.as_deref()),
            )
        } else {
            (1, 1, 1, 1, 0, 1, 100, 100)
        };

    PickCandidateInfo {
        assignment_idx,
        mane_select,
        mane_plus_clinical,
        canonical,
        ccds,
        length,
        biotype,
        tsl,
        appris,
        ensembl: pick_assignment_source_rank(tx, "ensembl"),
        refseq: pick_assignment_source_rank(tx, "refseq"),
        rank: pick_assignment_rank(tc),
    }
}

fn pick_assignment_gene_key(
    tc: &TranscriptConsequence,
    ctx: &PreparedContext<'_>,
) -> Option<String> {
    // Traceability:
    // - Ensembl VEP release 115 groups transcript consequences by
    //   `transcript->{_gene_stable_id}`, only falling back to a gene-adaptor
    //   lookup when that cache slot is absent
    //   <https://github.com/Ensembl/ensembl-vep/blob/release/115/modules/Bio/EnsEMBL/VEP/OutputFactory.pm#L845-L855>
    //
    // We intentionally require the cache-provided `gene_stable_id` here and do
    // not substitute `gene_symbol` or `transcript_id`, because those are local
    // heuristics rather than VEP pick-group semantics.
    let tx = ctx.transcripts.get(tc.transcript_idx?)?;
    tx.gene_stable_id.clone()
}

fn pick_assignment_source_rank(tx: Option<&TranscriptFeature>, wanted: &str) -> u8 {
    let Some(tx) = tx else {
        return 1;
    };

    // Traceability:
    // - Ensembl VEP release 115 sets `$info->{lc($tr->{_source_cache})} = 0`
    //   and only exact `ensembl` / `refseq` keys influence those categories
    //   <https://github.com/Ensembl/ensembl-vep/blob/release/115/modules/Bio/EnsEMBL/VEP/OutputFactory.pm#L734-L741>
    if tx
        .source
        .as_deref()
        .is_some_and(|source| source.eq_ignore_ascii_case(wanted))
    {
        0
    } else {
        1
    }
}

fn pick_assignment_rank(tc: &TranscriptConsequence) -> u32 {
    most_severe_term(tc.terms.iter())
        .map(|term| term.rank() as u32)
        .unwrap_or(1000)
}

fn retain_best_pick_candidates<T, F>(
    candidates: &mut Vec<PickCandidateInfo>,
    mut key: F,
) -> Option<usize>
where
    T: Ord + Copy,
    F: FnMut(&PickCandidateInfo) -> T,
{
    candidates.sort_by_key(|candidate| key(candidate));
    let best = key(&candidates[0]);
    let keep = candidates
        .iter()
        .take_while(|candidate| key(candidate) == best)
        .count();
    if keep == 1 {
        return Some(candidates[0].assignment_idx);
    }
    candidates.truncate(keep);
    None
}

fn pick_worst_assignment(
    candidate_indices: &[usize],
    assignments: &[TranscriptConsequence],
    ctx: &PreparedContext<'_>,
    pick_order: &[PickCriterion],
) -> Option<usize> {
    // Traceability:
    // - Ensembl VEP release 115
    //   `OutputFactory::pick_worst_VariationFeatureOverlapAllele()`
    //   <https://github.com/Ensembl/ensembl-vep/blob/release/115/modules/Bio/EnsEMBL/VEP/OutputFactory.pm#L679-L809>
    //
    // VEP computes category values for all candidates, stably sorts on one
    // `pick_order` criterion at a time, keeps only the tied best candidates,
    // and then repeats with the next criterion.
    let mut candidates: Vec<PickCandidateInfo> = candidate_indices
        .iter()
        .copied()
        .map(|assignment_idx| build_pick_candidate_info(assignment_idx, assignments, ctx))
        .collect();

    if candidates.is_empty() {
        return None;
    }
    if candidates.len() == 1 {
        return Some(candidates[0].assignment_idx);
    }

    for criterion in pick_order {
        let winner = match criterion {
            PickCriterion::ManeSelect => {
                retain_best_pick_candidates(&mut candidates, |candidate| candidate.mane_select)
            }
            PickCriterion::ManePlusClinical => {
                retain_best_pick_candidates(&mut candidates, |candidate| {
                    candidate.mane_plus_clinical
                })
            }
            PickCriterion::Canonical => {
                retain_best_pick_candidates(&mut candidates, |candidate| candidate.canonical)
            }
            PickCriterion::Appris => {
                retain_best_pick_candidates(&mut candidates, |candidate| candidate.appris)
            }
            PickCriterion::Tsl => {
                retain_best_pick_candidates(&mut candidates, |candidate| candidate.tsl)
            }
            PickCriterion::Biotype => {
                retain_best_pick_candidates(&mut candidates, |candidate| candidate.biotype)
            }
            PickCriterion::Ccds => {
                retain_best_pick_candidates(&mut candidates, |candidate| candidate.ccds)
            }
            PickCriterion::Rank => {
                retain_best_pick_candidates(&mut candidates, |candidate| candidate.rank)
            }
            PickCriterion::Length => {
                retain_best_pick_candidates(&mut candidates, |candidate| candidate.length)
            }
            PickCriterion::Ensembl => {
                retain_best_pick_candidates(&mut candidates, |candidate| candidate.ensembl)
            }
            PickCriterion::Refseq => {
                retain_best_pick_candidates(&mut candidates, |candidate| candidate.refseq)
            }
        };
        if winner.is_some() {
            return winner;
        }
    }

    candidates.first().map(|candidate| candidate.assignment_idx)
}

fn append_pick_flag(flags: &mut Option<String>) {
    if flags
        .as_deref()
        .is_some_and(|value| value.split('&').any(|part| part == "PICK"))
    {
        return;
    }
    match flags {
        Some(existing) if !existing.is_empty() => existing.push_str("&PICK"),
        _ => *flags = Some("PICK".to_string()),
    }
}

fn mark_flag_pick_allele_gene(
    assignments: &mut [TranscriptConsequence],
    ctx: &PreparedContext<'_>,
    pick_flags: &PickFlags,
) {
    if !pick_flags.flag_pick_allele_gene {
        return;
    }

    let mut grouped: HashMap<String, Vec<usize>> = HashMap::new();
    for (idx, tc) in assignments.iter().enumerate() {
        let Some(gene_key) = pick_assignment_gene_key(tc, ctx) else {
            continue;
        };
        grouped.entry(gene_key).or_default().push(idx);
    }

    for candidate_indices in grouped.values() {
        let Some(winner) =
            pick_worst_assignment(candidate_indices, assignments, ctx, &pick_flags.pick_order)
        else {
            continue;
        };
        append_pick_flag(&mut assignments[winner].flags);
    }
}

/// Compute miRNA CSQ field from ncRNA secondary structure and variant cDNA position.
///
/// Traceability:
/// - Ensembl VEP `OutputFactory::TranscriptVariationAllele_to_output_hash()`
///   <https://github.com/Ensembl/ensembl-vep/blob/release/115/modules/Bio/EnsEMBL/VEP/OutputFactory.pm#L1572-L1612>
///
/// The ncRNA attribute value has format `"start:end structure_string"` where:
/// - `start`/`end` are 1-based cDNA positions of the structure
/// - structure string uses dot-bracket notation with optional RLE counts:
///   `(19` = 19 open-parens (stem), `.6` = 6 dots (loop), bare char = count 1
///
/// VEP maps variant cDNA positions into the expanded structure array:
///   `struct_index = cdna_pos - ncrna_start`
/// Characters `(` and `)` → `miRNA_stem`, `.` → `miRNA_loop`.
/// Output is sorted, `&`-joined SO terms after deduplicating raw structure
/// characters. Because `(` and `)` are distinct before mapping, VEP can emit
/// `miRNA_stem` twice for a single overlapped region.
fn mirna_structure_field(
    ncrna_structure: Option<&str>,
    biotype: &str,
    cdna_start: Option<usize>,
    cdna_end: Option<usize>,
) -> String {
    if biotype != "miRNA" {
        return String::new();
    }
    let Some(raw) = ncrna_structure else {
        return String::new();
    };
    let Some(cdna_s) = cdna_start else {
        return String::new();
    };
    let Some(cdna_e) = cdna_end else {
        return String::new();
    };

    // Parse ncRNA structure. Two formats supported:
    // 1. Full attribute: "start:end structure_string" (e.g. "1:81 (19.(6...")
    // 2. Structure only: "(19.(6..." (from parquet ncrna_structure column)
    //
    // When start:end prefix is missing, assume structure starts at cDNA position 1.
    let parts: Vec<&str> = raw
        .splitn(3, |c: char| c.is_whitespace() || c == ':')
        .collect();
    let (struct_start, struct_str) = if parts.len() >= 3 {
        if let (Ok(s), Ok(_e)) = (parts[0].parse::<usize>(), parts[1].parse::<usize>()) {
            (s, parts[2])
        } else {
            (1usize, raw)
        }
    } else {
        (1usize, raw)
    };

    let (cs, ce) = if cdna_s <= cdna_e {
        (cdna_s, cdna_e)
    } else {
        (cdna_e, cdna_s)
    };

    // Expand RLE structure: "(19" → 19 '(' chars, ".6" → 6 '.' chars, bare char → 1.
    let mut expanded: Vec<u8> = Vec::new();
    let bytes = struct_str.as_bytes();
    let mut i = 0;
    while i < bytes.len() {
        let ch = bytes[i];
        if ch == b'(' || ch == b')' || ch == b'.' {
            // Read optional count after the character.
            let mut count = 0usize;
            let mut j = i + 1;
            while j < bytes.len() && bytes[j].is_ascii_digit() {
                count = count * 10 + (bytes[j] - b'0') as usize;
                j += 1;
            }
            if count == 0 {
                count = 1;
            }
            for _ in 0..count {
                expanded.push(ch);
            }
            i = j;
        } else {
            i += 1;
        }
    }

    // Map variant cDNA positions to structure indices and collect the distinct
    // raw structure characters overlapped by the variant, matching VEP.
    let mut has_open_stem = false;
    let mut has_close_stem = false;
    let mut has_loop = false;
    for pos in cs..=ce {
        if pos < struct_start {
            continue;
        }
        let idx = pos - struct_start;
        if idx >= expanded.len() {
            continue;
        }
        match expanded[idx] {
            b'(' => has_open_stem = true,
            b')' => has_close_stem = true,
            b'.' => has_loop = true,
            _ => {}
        }
    }

    let mut terms = Vec::new();
    if has_open_stem {
        terms.push("miRNA_stem");
    }
    if has_close_stem {
        terms.push("miRNA_stem");
    }
    if has_loop {
        terms.push("miRNA_loop");
    }

    terms.sort_unstable();
    terms.join("&")
}

/// Look up SIFT and PolyPhen predictions from the per-transcript LRU cache.
///
/// Traceability:
/// - Ensembl VEP `OutputFactory::add_sift_polyphen()`
///   <https://github.com/Ensembl/ensembl-vep/blob/release/115/modules/Bio/EnsEMBL/VEP/OutputFactory.pm#L1746-L1799>
/// - Ensembl Variation `TranscriptVariationAllele::_protein_function_prediction()`
///   <https://github.com/Ensembl/ensembl-variation/blob/release/115/modules/Bio/EnsEMBL/Variation/TranscriptVariationAllele.pm>
///
/// VEP requires a single amino acid substitution (`A/B` pattern in `pep_allele_string`).
/// Lookup is by (protein_position, alt_amino_acid). Output uses `--sift b` / `--polyphen b`
/// format: `prediction(score)` with spaces replaced by underscores.
///
/// The cache is populated by `load_sift_polyphen_cache()` with one SQL point
/// query per transcript. Each transcript's predictions are stored in a
/// HashMap for O(1) lookup. See biodatageeks/datafusion-bio-functions#38.
fn lookup_sift_polyphen(
    transcript_id: Option<&str>,
    protein_position: Option<&str>,
    amino_acids: Option<&str>,
    cache: &mut SiftPolyphenCache,
    #[cfg(feature = "kv-cache")] sift_kv: &Option<crate::kv_cache::SiftKvStore>,
    #[cfg(not(feature = "kv-cache"))] _sift_kv: &Option<()>,
) -> (String, String) {
    let empty = || (String::new(), String::new());

    // VEP only produces SIFT/PolyPhen for single amino acid substitutions.
    let Some(aa) = amino_acids else {
        return empty();
    };
    let parts: Vec<&str> = aa.split('/').collect();
    if parts.len() != 2 || parts[0].len() != 1 || parts[1].len() != 1 {
        return empty();
    }
    let alt_aa = parts[1];

    // Parse protein position (may be a range like "42-43" for indels — skip those).
    let Some(pos_str) = protein_position else {
        return empty();
    };
    let Ok(pos) = pos_str.parse::<i32>() else {
        return empty();
    };

    let Some(tx_id) = transcript_id else {
        return empty();
    };

    // Lazy load from fjall sift keyspace on cache miss.
    if cache.get(tx_id).is_none() {
        #[cfg(feature = "kv-cache")]
        if let Some(kv) = sift_kv {
            if let Ok(Some(preds)) = kv.get(tx_id) {
                cache.insert(tx_id.to_string(), preds, i64::MAX);
            }
        }
    }

    let Some(preds) = cache.get(tx_id) else {
        return empty();
    };

    let sift = preds
        .lookup_sift(pos, alt_aa)
        .map(|(pred, score)| format_prediction(pred, score))
        .unwrap_or_default();
    let polyphen = preds
        .lookup_polyphen(pos, alt_aa)
        .map(|(pred, score)| format_prediction(pred, score))
        .unwrap_or_default();

    (sift, polyphen)
}

/// Format a SIFT/PolyPhen prediction as `prediction(score)` with spaces→underscores.
///
/// VEP uses `--sift b` / `--polyphen b` format (both prediction and score).
fn format_prediction(prediction: &str, score: f32) -> String {
    let pred = prediction.replace(' ', "_").replace("_-_", "_");
    format!("{pred}({score})")
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
    transcript_selection: TranscriptSelectionFlags,
    annotation_column_defs: Vec<AnnotationColumnDef>,
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
    ) -> Result<Self> {
        let transcript_selection =
            TranscriptSelectionFlags::from_options_json(options_json.as_deref())?;
        let annotation_column_defs = annotation_column_defs_for_selection(transcript_selection);

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

        fields.push(Arc::new(Field::new("CSQ", DataType::Utf8, true)));
        fields.push(Arc::new(Field::new(
            "most_severe_consequence",
            DataType::Utf8,
            true,
        )));
        for col_def in &annotation_column_defs {
            fields.push(Arc::new(Field::new(
                col_def.name,
                col_def.data_type.clone(),
                true,
            )));
        }

        Ok(Self {
            session,
            vcf_table,
            cache_source,
            backend,
            options_json,
            transcript_selection,
            annotation_column_defs,
            schema: Arc::new(Schema::new(fields)),
        })
    }

    fn escaped_sql_literal(value: &str) -> String {
        value.replace('\'', "''")
    }

    fn vcf_field_count(&self) -> usize {
        self.schema
            .fields()
            .len()
            .saturating_sub(self.annotation_column_count())
    }

    fn annotation_column_count(&self) -> usize {
        self.annotation_column_defs.len() + 2
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

    fn parse_json_i64_option(json: &str, key: &str) -> Option<i64> {
        let needle = format!("\"{key}\"");
        let start = json.find(&needle)?;
        let rest = &json[start + needle.len()..];
        let colon = rest.find(':')?;
        let after_colon = rest[colon + 1..].trim_start();

        let digits_len = after_colon
            .chars()
            .take_while(|ch| *ch == '-' || ch.is_ascii_digit())
            .count();
        if digits_len == 0 {
            return None;
        }

        after_colon[..digits_len].parse().ok()
    }

    /// Traceability:
    /// - Ensembl VEP `Config.pm` `distance` option
    ///   <https://github.com/Ensembl/ensembl-vep/blob/release/115/modules/Bio/EnsEMBL/VEP/Config.pm#L145-L155>
    /// - Ensembl VEP `BaseRunner::_set_package_variables()`
    ///   <https://github.com/Ensembl/ensembl-vep/blob/release/115/modules/Bio/EnsEMBL/VEP/BaseRunner.pm#L499-L511>
    ///
    /// Ensembl accepts `--distance N` and `--distance U,D`, applying the
    /// first value to upstream and the second to downstream when provided.
    fn parse_json_distance_option(json: &str) -> Option<(i64, i64)> {
        if let Some(distance) = Self::parse_json_i64_option(json, "distance") {
            if distance >= 0 {
                return Some((distance, distance));
            }
        }

        let raw = Self::parse_json_string_option(json, "distance")?;
        let parts: Vec<&str> = raw.split(',').map(str::trim).collect();
        let parse_part = |value: &str| value.parse::<i64>().ok().filter(|parsed| *parsed >= 0);

        match parts.as_slice() {
            [single] => parse_part(single).map(|distance| (distance, distance)),
            [upstream, downstream] => Some((parse_part(upstream)?, parse_part(downstream)?)),
            _ => None,
        }
    }

    fn transcript_distance_config(&self) -> (i64, i64) {
        self.options_json
            .as_deref()
            .and_then(Self::parse_json_distance_option)
            .unwrap_or((5000, 5000))
    }

    /// Traceability:
    /// - Ensembl VEP `OutputFactory::TranscriptVariationAllele_to_output_hash()`
    ///   <https://github.com/Ensembl/ensembl-vep/blob/release/115/modules/Bio/EnsEMBL/VEP/OutputFactory.pm#L1706-L1715>
    ///
    /// Ensembl strips the translation version only at output time when
    /// `remove_hgvsp_version` is enabled, and URI-escapes `=` only when
    /// `no_escape` is false.
    fn format_hgvsp_output(
        raw_hgvsp: &str,
        remove_hgvsp_version: bool,
        no_escape: bool,
        prediction_format: bool,
    ) -> String {
        let mut out = raw_hgvsp.to_string();

        if remove_hgvsp_version {
            if let Some((ref_name, suffix)) = out.split_once(":p.") {
                let stripped_ref = ref_name
                    .rsplit_once('.')
                    .filter(|(_, version)| version.chars().all(|ch| ch.is_ascii_digit()))
                    .map(|(base, _)| base)
                    .unwrap_or(ref_name);
                out = format!("{stripped_ref}:p.{suffix}");
            }
        }

        if prediction_format {
            if let Some((ref_name, suffix)) = out.split_once(":p.") {
                out = format!("{ref_name}:p.({suffix})");
            }
        }

        if !no_escape {
            out = out.replace('=', "%3D");
        }

        out
    }

    /// Traceability:
    /// - Ensembl VEP `Runner::post_setup_checks()`
    ///   <https://github.com/Ensembl/ensembl-vep/blob/release/115/modules/Bio/EnsEMBL/VEP/Runner.pm#L726-L738>
    ///
    /// VEP refuses offline HGVS output without an available FASTA. Our table
    /// function is also offline/cache-backed, so require `reference_fasta_path`
    /// whenever HGVS output is requested. Our runtime uses indexed FASTA
    /// access, so validate that the indexed reader can actually be opened here
    /// rather than failing later during execution.
    fn validate_hgvs_reference_fasta(
        hgvs_flags: HgvsFlags,
        reference_fasta_path: Option<&str>,
    ) -> Result<()> {
        if !hgvs_flags.any() {
            return Ok(());
        }
        let Some(path) = reference_fasta_path else {
            return Err(DataFusionError::Execution(
                "annotate_vep(): Cannot generate HGVS coordinates (--hgvs/--hgvsc/--hgvsp) without reference_fasta_path (VEP --fasta)".to_string(),
            ));
        };
        if !std::path::Path::new(path).exists() {
            return Err(DataFusionError::Execution(format!(
                "annotate_vep(): reference_fasta_path does not exist: {path}"
            )));
        }
        fasta::io::indexed_reader::Builder::default()
            .build_from_path(path)
            .map_err(|e| {
                DataFusionError::Execution(format!(
                    "annotate_vep(): failed to open indexed reference FASTA '{path}': {e}"
                ))
            })?;
        Ok(())
    }

    async fn projected_columns_for_table(
        session: &SessionContext,
        table: &str,
        wanted: &[&str],
    ) -> String {
        let Ok(provider) = session.table(table).await else {
            return "*".to_string();
        };
        let schema = provider.schema();
        let arrow = schema.as_arrow();
        let field_names: HashSet<&str> = arrow.fields().iter().map(|f| f.name().as_str()).collect();
        let mut cols: Vec<&str> = Vec::new();
        for &col in wanted {
            let bare = col.trim_matches('"');
            if field_names.contains(bare) {
                cols.push(col);
            }
        }
        if cols.is_empty() {
            "*".to_string()
        } else {
            cols.join(", ")
        }
    }

    async fn load_transcripts(
        &self,
        table: &str,
        worklist: &MissWorklist,
    ) -> Result<(Vec<TranscriptFeature>, HashMap<String, String>)> {
        let filter = worklist.chrom_filter_clause();
        let cols = Self::projected_columns_for_table(
            &self.session,
            table,
            &[
                "transcript_id",
                "stable_id",
                "chrom",
                "start",
                "\"end\"",
                "strand",
                "biotype",
                "cds_start",
                "cds_end",
                "cdna_coding_start",
                "cdna_coding_end",
                "gene_stable_id",
                "gene_symbol",
                "gene_symbol_source",
                "gene_hgnc_id_native",
                "gene_hgnc_id",
                "source",
                "version",
                "raw_object_json",
                "cds_start_nf",
                "cds_end_nf",
                "mature_mirna_regions",
                "cdna_seq",
                "bam_edit_status",
                "has_non_polya_rna_edit",
                "spliced_seq",
                "translateable_seq",
                "flags_str",
                "cdna_mapper_segments",
                "is_canonical",
                "tsl",
                "mane_select",
                "mane_plus_clinical",
                "translation_stable_id",
                "gene_phenotype",
                "ccds",
                "swissprot",
                "trembl",
                "uniparc",
                "uniprot_isoform",
                "appris",
                "ncrna_structure",
            ],
        )
        .await;
        let query = format!("SELECT {cols} FROM `{table}`{filter}");
        let batches = self.session.sql(&query).await?.collect().await?;
        let mut out = Vec::new();
        let mut translateable_seq_by_tx = HashMap::new();
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
            let cdna_coding_start_idx = schema.index_of("cdna_coding_start").ok();
            let cdna_coding_end_idx = schema.index_of("cdna_coding_end").ok();
            let gene_stable_id_idx = schema.index_of("gene_stable_id").ok();
            let gene_symbol_idx = schema.index_of("gene_symbol").ok();
            let gene_symbol_source_idx = schema.index_of("gene_symbol_source").ok();
            let gene_hgnc_id_native_idx = schema.index_of("gene_hgnc_id_native").ok();
            let gene_hgnc_id_idx = schema.index_of("gene_hgnc_id").ok();
            let source_idx = schema.index_of("source").ok();
            let version_idx = schema.index_of("version").ok();
            let raw_object_json_idx = schema.index_of("raw_object_json").ok();
            let cds_start_nf_idx = schema.index_of("cds_start_nf").ok();
            let cds_end_nf_idx = schema.index_of("cds_end_nf").ok();
            let mirna_regions_idx = schema.index_of("mature_mirna_regions").ok();
            let cdna_seq_idx = schema.index_of("cdna_seq").ok();
            // Promoted columns (previously extracted from raw_object_json).
            let bam_edit_status_idx = schema.index_of("bam_edit_status").ok();
            let has_non_polya_rna_edit_idx = schema.index_of("has_non_polya_rna_edit").ok();
            let spliced_seq_idx = schema.index_of("spliced_seq").ok();
            let translateable_seq_idx = schema.index_of("translateable_seq").ok();
            let five_prime_utr_seq_idx = schema.index_of("five_prime_utr_seq").ok();
            let three_prime_utr_seq_idx = schema.index_of("three_prime_utr_seq").ok();
            let flags_str_idx = schema.index_of("flags_str").ok();
            let cdna_mapper_segments_idx = schema.index_of("cdna_mapper_segments").ok();
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
            let appris_idx = schema.index_of("appris").ok();
            let ncrna_structure_idx = schema.index_of("ncrna_structure").ok();

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
                let cdna_coding_start = cdna_coding_start_idx
                    .and_then(|idx| int64_at(batch.column(idx).as_ref(), row))
                    .filter(|v| *v > 0)
                    .and_then(|v| usize::try_from(v).ok());
                let cdna_coding_end = cdna_coding_end_idx
                    .and_then(|idx| int64_at(batch.column(idx).as_ref(), row))
                    .filter(|v| *v > 0)
                    .and_then(|v| usize::try_from(v).ok());

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
                let cdna_mapper_segments = cdna_mapper_segments_idx
                    .map(|idx| {
                        cdna_mapper_segments_from_list_column(batch.column(idx).as_ref(), row)
                    })
                    .unwrap_or_default();
                let raw_metadata = raw_object_json_idx
                    .and_then(|idx| string_at(batch.column(idx).as_ref(), row))
                    .map(|raw| parse_transcript_raw_metadata(&raw))
                    .unwrap_or_default();
                let TranscriptRawMetadata {
                    display_xref_id,
                    source: raw_source,
                    gene_hgnc_id_native: raw_gene_hgnc_id_native,
                    refseq_match,
                    refseq_edits,
                    cdna_mapper_segments: raw_cdna_mapper_segments,
                    spliced_seq: raw_spliced_seq,
                    five_prime_utr_seq,
                    three_prime_utr_seq,
                    translateable_seq: raw_translateable_seq,
                    flags_str: raw_flags_str,
                    is_gencode_basic,
                    is_gencode_primary,
                } = raw_metadata;
                let cdna_mapper_segments = if cdna_mapper_segments.is_empty() {
                    raw_cdna_mapper_segments
                } else {
                    cdna_mapper_segments
                };
                // Ensembl release/115 computes alternate CDS from the live
                // transcript object's `_translateable_seq()` / 3'UTR rather
                // than reconstructing from genomic exons when that state is
                // already cached on the transcript object.
                // https://github.com/Ensembl/ensembl-variation/blob/release/115/modules/Bio/EnsEMBL/Variation/TranscriptVariationAllele.pm#L2470-L2481
                let translateable_seq = translateable_seq_idx
                    .and_then(|idx| string_at(batch.column(idx).as_ref(), row))
                    .or(raw_translateable_seq);
                let five_prime_utr_seq = five_prime_utr_seq_idx
                    .and_then(|idx| string_at(batch.column(idx).as_ref(), row))
                    .or(five_prime_utr_seq);
                let three_prime_utr_seq = three_prime_utr_seq_idx
                    .and_then(|idx| string_at(batch.column(idx).as_ref(), row))
                    .or(three_prime_utr_seq);
                if let Some(seq) = translateable_seq.as_ref() {
                    translateable_seq_by_tx.insert(transcript_id.clone(), seq.clone());
                }
                let flags_str = flags_str_idx
                    .and_then(|idx| string_at(batch.column(idx).as_ref(), row))
                    .or(raw_flags_str)
                    .or_else(|| flags_str_from_bools(cds_start_nf, cds_end_nf));

                let gene_stable_id =
                    gene_stable_id_idx.and_then(|idx| string_at(batch.column(idx).as_ref(), row));
                let gene_symbol =
                    gene_symbol_idx.and_then(|idx| string_at(batch.column(idx).as_ref(), row));
                let gene_symbol_source = gene_symbol_source_idx
                    .and_then(|idx| string_at(batch.column(idx).as_ref(), row));
                let gene_hgnc_id_native = gene_hgnc_id_native_idx
                    .and_then(|idx| string_at(batch.column(idx).as_ref(), row))
                    .or(raw_gene_hgnc_id_native);
                let promoted_gene_hgnc_id =
                    gene_hgnc_id_idx.and_then(|idx| string_at(batch.column(idx).as_ref(), row));
                let gene_hgnc_id = gene_hgnc_id_native.clone().or(promoted_gene_hgnc_id);
                let source = raw_source.or_else(|| {
                    source_idx
                        .and_then(|idx| string_at(batch.column(idx).as_ref(), row))
                        .and_then(|value| normalize_source_label(&value))
                });
                let bam_edit_status =
                    bam_edit_status_idx.and_then(|idx| string_at(batch.column(idx).as_ref(), row));
                let has_non_polya_rna_edit = has_non_polya_rna_edit_idx
                    .and_then(|idx| bool_at(batch.column(idx).as_ref(), row))
                    .unwrap_or(false);
                let cdna_seq =
                    cdna_seq_idx.and_then(|idx| string_at(batch.column(idx).as_ref(), row));
                let spliced_seq = spliced_seq_idx
                    .and_then(|idx| string_at(batch.column(idx).as_ref(), row))
                    .or(raw_spliced_seq)
                    .or_else(|| {
                        synthesize_spliced_seq(
                            five_prime_utr_seq.as_deref(),
                            translateable_seq.as_deref(),
                            three_prime_utr_seq.as_deref(),
                            cdna_coding_start,
                            cdna_coding_end,
                            cdna_seq.as_deref(),
                        )
                    });
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
                let mane_select =
                    mane_select_idx.and_then(|idx| string_at(batch.column(idx).as_ref(), row));
                let mane_plus_clinical = mane_plus_clinical_idx
                    .and_then(|idx| string_at(batch.column(idx).as_ref(), row));
                let translation_stable_id = translation_stable_id_idx
                    .and_then(|idx| string_at(batch.column(idx).as_ref(), row));
                let gene_phenotype = gene_phenotype_idx
                    .and_then(|idx| bool_at(batch.column(idx).as_ref(), row))
                    .unwrap_or(false);
                let ccds = ccds_idx.and_then(|idx| string_at(batch.column(idx).as_ref(), row));
                let swissprot =
                    swissprot_idx.and_then(|idx| string_at(batch.column(idx).as_ref(), row));
                let trembl = trembl_idx.and_then(|idx| string_at(batch.column(idx).as_ref(), row));
                let uniparc =
                    uniparc_idx.and_then(|idx| string_at(batch.column(idx).as_ref(), row));
                let uniprot_isoform =
                    uniprot_isoform_idx.and_then(|idx| string_at(batch.column(idx).as_ref(), row));
                let appris = appris_idx.and_then(|idx| string_at(batch.column(idx).as_ref(), row));
                let ncrna_structure =
                    ncrna_structure_idx.and_then(|idx| string_at(batch.column(idx).as_ref(), row));

                out.push(TranscriptFeature {
                    transcript_id,
                    chrom,
                    start,
                    end,
                    strand,
                    biotype,
                    cds_start,
                    cds_end,
                    cdna_coding_start,
                    cdna_coding_end,
                    cdna_mapper_segments,
                    mature_mirna_regions,
                    gene_stable_id,
                    gene_symbol,
                    gene_symbol_source,
                    gene_hgnc_id_native,
                    gene_hgnc_id,
                    display_xref_id,
                    source,
                    refseq_match,
                    refseq_edits,
                    is_gencode_basic,
                    is_gencode_primary,
                    bam_edit_status,
                    has_non_polya_rna_edit,
                    spliced_seq,
                    five_prime_utr_seq,
                    three_prime_utr_seq,
                    translateable_seq,
                    cdna_seq,
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
                    appris,
                    ncrna_structure,
                });
            }
        }

        Ok((out, translateable_seq_by_tx))
    }

    async fn load_exons(&self, table: &str, worklist: &MissWorklist) -> Result<Vec<ExonFeature>> {
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
            worklist.chrom_filter_clause()
        } else {
            String::new()
        };
        let cols = Self::projected_columns_for_table(
            &self.session,
            table,
            &[
                "transcript_id",
                "stable_id",
                "exon_number",
                "start",
                "\"end\"",
                "chrom",
            ],
        )
        .await;
        let query = format!("SELECT {cols} FROM `{table}`{filter}");
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
        worklist: &MissWorklist,
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
            worklist.chrom_filter_clause()
        } else {
            String::new()
        };
        let cols = Self::projected_columns_for_table(
            &self.session,
            table,
            &[
                "transcript_id",
                "stable_id",
                "chrom",
                "start",
                "\"end\"",
                "cds_len",
                "cds_length",
                "protein_len",
                "translation_seq",
                "cds_sequence",
                "cds_seq",
                "coding_sequence",
                // Upstream d26e370+fixup added canonical (pre-BAM-edit)
                // counterparts that drive HGVSp flanking / dup / 3'-shift
                // checks against VEP's canonical reference rather than the
                // BAM-edited peptide. Must be explicitly projected here or
                // the row-level `schema.index_of(...)` below misses them
                // and `translation_seq_canonical` ends up None for every
                // translation — the whole canonical pipeline silently
                // degrades back to BAM-edited behavior.
                "translation_seq_canonical",
                "cds_sequence_canonical",
                "version",
                "protein_features",
            ],
        )
        .await;
        let query = format!("SELECT {cols} FROM `{table}`{filter}");
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
            // Upstream rev d26e370 added `translation_seq_canonical` and
            // `cds_sequence_canonical`: canonical (pre-BAM-edit) values that
            // VEP uses for HGVSp. Fall back to the BAM-edited columns when
            // the canonical columns are absent (older parquet caches).
            let translation_seq_canonical_idx = schema.index_of("translation_seq_canonical").ok();
            let cds_seq_idx = schema
                .index_of("cds_sequence")
                .or_else(|_| schema.index_of("cds_seq"))
                .or_else(|_| schema.index_of("coding_sequence"))
                .ok();
            let cds_seq_canonical_idx = schema.index_of("cds_sequence_canonical").ok();
            let tl_stable_id_idx = schema.index_of("stable_id").ok();
            let tl_version_idx = schema.index_of("version").ok();
            // SIFT/PolyPhen predictions are NOT loaded here — they are loaded
            // lazily per-transcript via SiftPolyphenCache to avoid ~20GB memory.
            // See biodatageeks/datafusion-bio-functions#38.
            let pf_idx = schema.index_of("protein_features").ok();
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
                // Canonical columns (upstream d26e370). Strict: if the parquet
                // doesn't carry them, `TranslationFeature.translation_seq_canonical`
                // stays `None` and the HGVSp helper falls through to its
                // caller-supplied fallback — no cache-level clone of the
                // BAM-edited value. Legacy parquet caches must be regenerated.
                let translation_seq_canonical = translation_seq_canonical_idx
                    .and_then(|idx| string_at(batch.column(idx).as_ref(), row));
                let cds_sequence_canonical = cds_seq_canonical_idx
                    .and_then(|idx| string_at(batch.column(idx).as_ref(), row));
                let tl_stable_id =
                    tl_stable_id_idx.and_then(|idx| string_at(batch.column(idx).as_ref(), row));
                let tl_version = tl_version_idx
                    .and_then(|idx| int64_at(batch.column(idx).as_ref(), row))
                    .and_then(|v| i32::try_from(v).ok());
                let protein_features = pf_idx
                    .map(|idx| read_protein_features(batch.column(idx).as_ref(), row))
                    .unwrap_or_default();

                out.push(TranslationFeature {
                    transcript_id,
                    cds_len,
                    protein_len,
                    translation_seq,
                    cds_sequence,
                    translation_seq_canonical,
                    cds_sequence_canonical,
                    stable_id: tl_stable_id,
                    version: tl_version,
                    protein_features,
                });
            }
        }

        Ok(out)
    }

    /// Window size for sliding-window SIFT/PolyPhen loading (5MB).
    const SIFT_WINDOW_SIZE: i64 = 5_000_000;

    /// Try to build a direct parquet reader for sift windows, bypassing DataFusion SQL.
    /// Returns cached metadata + projection + RG ranges if the file path can be resolved.
    fn build_sift_direct_reader(path: &str) -> Option<SiftDirectReader> {
        use parquet::arrow::ProjectionMask;
        use parquet::arrow::arrow_reader::{ArrowReaderMetadata, ArrowReaderOptions};
        use parquet::file::statistics::Statistics;

        let file = std::fs::File::open(path).ok()?;
        let arrow_meta = ArrowReaderMetadata::load(&file, ArrowReaderOptions::default()).ok()?;
        let parquet_schema = arrow_meta.metadata().file_metadata().schema_descr_ptr();

        // Find root column indices for projection
        let arrow_schema = arrow_meta.schema();
        let fields = arrow_schema.fields();
        let find_idx = |name: &str| fields.iter().position(|f| f.name() == name);

        let tid_root = find_idx("transcript_id")?;
        let end_root = find_idx("end")?;
        let sift_root = find_idx("sift_predictions")?;
        let poly_root = find_idx("polyphen_predictions")?;
        let chrom_root = find_idx("chrom");

        let mut proj_indices = vec![tid_root, end_root, sift_root, poly_root];
        if let Some(ci) = chrom_root {
            proj_indices.push(ci);
        }
        let projection = ProjectionMask::roots(&parquet_schema, proj_indices);

        // Pre-compute RG position ranges from column statistics
        let num_rgs = arrow_meta.metadata().num_row_groups();
        // Find physical column index for "start" and "end"
        let leaf_cols = parquet_schema.columns();
        let start_leaf = leaf_cols.iter().position(|c| c.name() == "start");
        let end_leaf = leaf_cols.iter().position(|c| c.name() == "end");

        let rg_ranges: Vec<(i64, i64)> = (0..num_rgs)
            .map(|i| {
                let rg = arrow_meta.metadata().row_group(i);
                let min_start = start_leaf
                    .and_then(|idx| rg.column(idx).statistics())
                    .and_then(|s| match s {
                        Statistics::Int64(v) => v.min_opt().copied(),
                        _ => None,
                    })
                    .unwrap_or(i64::MIN);
                let max_end = end_leaf
                    .and_then(|idx| rg.column(idx).statistics())
                    .and_then(|s| match s {
                        Statistics::Int64(v) => v.max_opt().copied(),
                        _ => None,
                    })
                    .unwrap_or(i64::MAX);
                (min_start, max_end)
            })
            .collect();

        Some(SiftDirectReader {
            path: path.to_string(),
            arrow_meta,
            projection,
            rg_ranges,
        })
    }

    /// Load SIFT/PolyPhen predictions for a single genomic window into the cache.
    ///
    /// Queries translations whose CDS overlaps the window `[win_start, win_end)`:
    /// ```sql
    /// SELECT transcript_id, "end", sift_predictions, polyphen_predictions
    /// FROM translations
    /// WHERE chrom = '1' AND start <= win_end AND "end" >= win_start
    /// ```
    ///
    /// Each window typically returns ~20-50 translations (~500K prediction
    /// entries ~20MB).
    ///
    /// With sorted parquet + small row groups (bio-formats#129), DataFusion
    /// uses row-group min/max statistics to skip non-matching row groups,
    /// reading only 1-2 row groups per window query instead of all.
    ///
    /// See biodatageeks/datafusion-bio-functions#38.
    async fn load_sift_window(
        &self,
        table: &str,
        chrom: &str,
        win_start: i64,
        win_end: i64,
        cache: &mut SiftPolyphenCache,
    ) -> Result<()> {
        let escaped_chrom = Self::escaped_sql_literal(chrom);
        let query = format!(
            "SELECT transcript_id, \"end\", sift_predictions, polyphen_predictions \
             FROM `{table}` \
             WHERE chrom = '{escaped_chrom}' \
               AND start <= {win_end} AND \"end\" >= {win_start}"
        );
        let batches = self.session.sql(&query).await?.collect().await?;

        for batch in &batches {
            let schema = batch.schema();
            let tx_idx = schema
                .index_of("transcript_id")
                .or_else(|_| schema.index_of("stable_id"))
                .ok();
            let end_col_idx = schema.index_of("end").ok();
            let sift_col_idx = schema.index_of("sift_predictions").ok();
            let pp_col_idx = schema.index_of("polyphen_predictions").ok();

            let Some(tx_idx) = tx_idx else { continue };

            for row in 0..batch.num_rows() {
                let Some(tx_id) = string_at(batch.column(tx_idx).as_ref(), row) else {
                    continue;
                };
                // Skip if already cached (window overlap with previous window).
                if cache.get(&tx_id).is_some() {
                    continue;
                }
                let genomic_end = end_col_idx
                    .and_then(|idx| int64_at(batch.column(idx).as_ref(), row))
                    .unwrap_or(i64::MAX);
                let mut preds = CachedPredictions::default();
                if let Some(idx) = sift_col_idx {
                    preds.sift = read_compact_predictions(batch.column(idx).as_ref(), row);
                }
                if let Some(idx) = pp_col_idx {
                    preds.polyphen = read_compact_predictions(batch.column(idx).as_ref(), row);
                }
                preds.sort();
                cache.insert(tx_id, preds, genomic_end);
            }
        }

        Ok(())
    }

    async fn load_regulatory_features(
        &self,
        table: &str,
        worklist: &MissWorklist,
    ) -> Result<Vec<RegulatoryFeature>> {
        let filter = worklist.interval_filter_sql();
        let cols = Self::projected_columns_for_table(
            &self.session,
            table,
            &[
                "stable_id",
                "feature_id",
                "feature_type",
                "chrom",
                "start",
                "\"end\"",
            ],
        )
        .await;
        let query = format!("SELECT {cols} FROM `{table}`{filter}");
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
        worklist: &MissWorklist,
    ) -> Result<Vec<MotifFeature>> {
        let filter = worklist.interval_filter_sql();
        let cols = Self::projected_columns_for_table(
            &self.session,
            table,
            &["motif_id", "feature_id", "chrom", "start", "\"end\""],
        )
        .await;
        let query = format!("SELECT {cols} FROM `{table}`{filter}");
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
        worklist: &MissWorklist,
    ) -> Result<Vec<MirnaFeature>> {
        let filter = worklist.interval_filter_sql();
        let cols = Self::projected_columns_for_table(
            &self.session,
            table,
            &["mirna_id", "feature_id", "chrom", "start", "\"end\""],
        )
        .await;
        let query = format!("SELECT {cols} FROM `{table}`{filter}");
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
        worklist: &MissWorklist,
    ) -> Result<Vec<StructuralFeature>> {
        let filter = worklist.interval_filter_sql();
        let cols = Self::projected_columns_for_table(
            &self.session,
            table,
            &[
                "feature_id",
                "stable_id",
                "feature_kind",
                "event_type",
                "chrom",
                "start",
                "\"end\"",
            ],
        )
        .await;
        let query = format!("SELECT {cols} FROM `{table}`{filter}");
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

    /// Discover contigs present in the VCF input.
    ///
    /// Prefers `"bio.vcf.contigs.indexed"` (TBI-derived, data-bearing only,
    /// zero cost), then falls back to `SELECT DISTINCT chrom` which scans
    /// the VCF but returns only contigs with actual data.
    ///
    /// Does NOT use `"bio.vcf.contigs"` (all VCF header contigs) because
    /// GRCh38 headers list ~195 contigs even when only 1–22 have data.
    async fn discover_vcf_contigs(&self) -> Result<Vec<String>> {
        let table = self.session.table(&self.vcf_table).await?;
        let schema = table.schema();
        let arrow_schema = schema.as_arrow();
        let metadata = arrow_schema.metadata();

        // Priority 1: TBI-indexed contigs (only data-bearing contigs, zero cost).
        if let Some(indexed_json) = metadata.get("bio.vcf.contigs.indexed") {
            if let Ok(contigs) = serde_json::from_str::<Vec<String>>(indexed_json) {
                if !contigs.is_empty() {
                    return Ok(contigs);
                }
            }
        }

        // Priority 2: scan chrom column for actual data-bearing contigs.
        let query = format!(
            "SELECT DISTINCT chrom FROM `{}`",
            Self::escaped_sql_literal(&self.vcf_table)
        );
        let batches = self.session.sql(&query).await?.collect().await?;
        let mut contigs = Vec::new();
        for batch in &batches {
            let col = batch.column(0);
            if let Some(arr) = col.as_any().downcast_ref::<StringArray>() {
                for i in 0..arr.len() {
                    if !arr.is_null(i) {
                        contigs.push(arr.value(i).to_string());
                    }
                }
            } else if let Some(arr) = col.as_any().downcast_ref::<LargeStringArray>() {
                for i in 0..arr.len() {
                    if !arr.is_null(i) {
                        contigs.push(arr.value(i).to_string());
                    }
                }
            } else if let Some(arr) = col.as_any().downcast_ref::<StringViewArray>() {
                for i in 0..arr.len() {
                    if !arr.is_null(i) {
                        contigs.push(arr.value(i).to_string());
                    }
                }
            }
        }
        Ok(contigs)
    }

    /// Partitioned per-contig annotation pipeline.
    ///
    /// For each contig discovered in the VCF:
    /// 1. Register per-chrom variation parquet, run lookup
    /// 2. Register per-chrom context parquet (transcript, exon, etc.)
    /// 3. Load context, build PreparedContext, annotate
    /// 4. Deregister ephemeral tables, free memory
    /// 5. Collect annotated batches
    #[allow(clippy::too_many_arguments)]
    async fn scan_with_transcript_engine_partitioned(
        &self,
        _state: &dyn Session,
        projection: Option<&Vec<usize>>,
        requested_columns: &[&str],
        extended_probes: bool,
        cache: &PartitionedParquetCache,
        translations_sift_table: Option<&str>,
        #[cfg(feature = "kv-cache")] kv_store: Option<Arc<crate::kv_cache::VepKvStore>>,
        fetch_limit: Option<usize>,
    ) -> Result<Arc<dyn ExecutionPlan>> {
        if profiling_enabled() {
            eprintln!("[VEP_PROFILE] ====== scan_with_transcript_engine_partitioned START ======");
        }

        let flags = VepFlags::from_options_json(self.options_json.as_deref());
        let hgvs_flags = HgvsFlags::from_options_json(self.options_json.as_deref());
        let transcript_selection = self.transcript_selection;
        let pick_flags = PickFlags::from_options_json(self.options_json.as_deref())?;
        let allowed_failed = self
            .options_json
            .as_deref()
            .and_then(|opts| Self::parse_json_i64_option(opts, "failed"))
            .unwrap_or(0);
        let reference_fasta_path = self
            .options_json
            .as_deref()
            .and_then(|opts| Self::parse_json_string_option(opts, "reference_fasta_path"));
        let input_buffer_size = self
            .options_json
            .as_deref()
            .and_then(|opts| Self::parse_json_i64_option(opts, "buffer_size"))
            .and_then(|value| usize::try_from(value).ok())
            .filter(|value| *value > 0)
            .unwrap_or(VEP_INPUT_BUFFER_SIZE);
        Self::validate_hgvs_reference_fasta(hgvs_flags, reference_fasta_path.as_deref())?;
        let (upstream_distance, downstream_distance) = self.transcript_distance_config();

        // Discover contigs from VCF.
        let t_contigs = profile_start!();
        let vcf_contigs = self.discover_vcf_contigs().await?;
        // Build expanded cache chrom set with both bare and chr-prefixed forms
        // so that VCF "chr1" matches cache "1" and vice versa.
        let mut cache_chroms: HashSet<String> = HashSet::new();
        for c in cache.available_chroms() {
            cache_chroms.insert(c.clone());
            if let Some(bare) = c.strip_prefix("chr") {
                cache_chroms.insert(bare.to_string());
            } else {
                cache_chroms.insert(format!("chr{c}"));
            }
        }
        let contigs: Vec<String> = vcf_contigs
            .iter()
            .filter(|c| cache_chroms.contains(c.as_str()))
            .cloned()
            .collect();
        profile_end!(
            "0. discover_contigs",
            t_contigs,
            format!(
                "{} VCF contigs, {} in cache, {} to process",
                vcf_contigs.len(),
                cache.available_chroms().len(),
                contigs.len()
            )
        );

        let cache_columns: Vec<String> = requested_columns.iter().map(|s| s.to_string()).collect();

        let projected_schema = if let Some(indices) = projection {
            Arc::new(self.schema.project(indices)?)
        } else {
            self.schema.clone()
        };

        let config = ContigAnnotationConfig {
            vcf_table: self.vcf_table.clone(),
            options_json: self.options_json.clone(),
            cache_columns,
            extended_probes,
            translations_sift_table: translations_sift_table.map(|s| s.to_string()),
            flags,
            hgvs_flags,
            transcript_selection,
            pick_flags,
            allowed_failed,
            reference_fasta_path,
            upstream_distance,
            downstream_distance,
            input_buffer_size,
            projection: projection.cloned(),
            annotation_column_count: self.annotation_column_count(),
            fetch_limit,
            #[cfg(feature = "kv-cache")]
            use_fjall: kv_store.is_some(),
            #[cfg(feature = "kv-cache")]
            sift_kv_store: kv_store.as_ref().and_then(|store| {
                let parent = store.root_path().parent()?;
                let sift_path = parent.join("translation_sift.fjall");
                crate::kv_cache::SiftKvStore::open_path(&sift_path)
                    .ok()
                    .flatten()
            }),
            #[cfg(feature = "kv-cache")]
            kv_store,
        };

        let exec = ContigAnnotationExec::new(
            projected_schema,
            self.schema.clone(),
            contigs,
            Arc::clone(&self.session),
            Arc::new(cache.clone()),
            config,
        );

        if profiling_enabled() {
            eprintln!(
                "[VEP_PROFILE] ====== scan_with_transcript_engine_partitioned: returning ContigAnnotationExec ======"
            );
        }

        Ok(Arc::new(exec))
    }

    #[allow(clippy::too_many_arguments)]
    fn annotate_batch_with_transcript_engine(
        &self,
        batch: &RecordBatch,
        engine: &TranscriptConsequenceEngine,
        ctx: &PreparedContext<'_>,
        colocated_map: &HashMap<ColocatedKey, ColocatedData>,
        sift_cache: &mut SiftPolyphenCache,
        #[cfg(feature = "kv-cache")] sift_kv: &Option<crate::kv_cache::SiftKvStore>,
        #[cfg(not(feature = "kv-cache"))] _sift_kv: &Option<()>,
        skip_csq: bool,
        skip_typed_cols: bool,
        flags: &VepFlags,
        hgvs_flags: &HgvsFlags,
        transcript_selection: TranscriptSelectionFlags,
        pick_flags: &PickFlags,
        hgvs_reference_reader: &mut Option<FastaReader>,
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

        let mut csq_builder = if skip_csq {
            StringBuilder::new()
        } else {
            StringBuilder::with_capacity(batch.num_rows(), batch.num_rows() * 40)
        };
        let mut most_builder =
            StringBuilder::with_capacity(batch.num_rows(), batch.num_rows() * 16);

        // --- Annotation column builders (87 columns) ---
        // Transcript-level (42)
        let mut b_allele = StringBuilder::with_capacity(batch.num_rows(), batch.num_rows() * 4);
        let mut b_consequence = ListBuilder::new(StringBuilder::new());
        let mut b_impact = ListBuilder::new(StringBuilder::new());
        let mut b_symbol = ListBuilder::new(StringBuilder::new());
        let mut b_gene = ListBuilder::new(StringBuilder::new());
        let mut b_feature_type = ListBuilder::new(StringBuilder::new());
        let mut b_feature = ListBuilder::new(StringBuilder::new());
        let mut b_biotype = ListBuilder::new(StringBuilder::new());
        let mut b_exon = ListBuilder::new(StringBuilder::new());
        let mut b_intron = ListBuilder::new(StringBuilder::new());
        let mut b_hgvsc = ListBuilder::new(StringBuilder::new());
        let mut b_hgvsp = ListBuilder::new(StringBuilder::new());
        let mut b_cdna_position = ListBuilder::new(StringBuilder::new());
        let mut b_cds_position = ListBuilder::new(StringBuilder::new());
        let mut b_protein_position = ListBuilder::new(StringBuilder::new());
        let mut b_amino_acids = ListBuilder::new(StringBuilder::new());
        let mut b_codons = ListBuilder::new(StringBuilder::new());
        let mut b_existing_variation = ListBuilder::new(StringBuilder::new());
        let mut b_distance = ListBuilder::new(Int64Builder::new());
        let mut b_strand = ListBuilder::new(Int8Builder::new());
        let mut b_flags = ListBuilder::new(StringBuilder::new());
        let mut b_variant_class =
            StringBuilder::with_capacity(batch.num_rows(), batch.num_rows() * 8);
        let mut b_symbol_source = ListBuilder::new(StringBuilder::new());
        let mut b_hgnc_id = ListBuilder::new(StringBuilder::new());
        let mut b_canonical = ListBuilder::new(StringBuilder::new());
        let mut b_mane = ListBuilder::new(StringBuilder::new());
        let mut b_mane_select = ListBuilder::new(StringBuilder::new());
        let mut b_mane_plus_clinical = ListBuilder::new(StringBuilder::new());
        let mut b_tsl = ListBuilder::new(Int8Builder::new());
        let mut b_appris = ListBuilder::new(StringBuilder::new());
        let mut b_ccds = ListBuilder::new(StringBuilder::new());
        let mut b_ensp = ListBuilder::new(StringBuilder::new());
        let mut b_swissprot = ListBuilder::new(StringBuilder::new());
        let mut b_trembl = ListBuilder::new(StringBuilder::new());
        let mut b_uniparc = ListBuilder::new(StringBuilder::new());
        let mut b_uniprot_isoform = ListBuilder::new(StringBuilder::new());
        let mut b_refseq_match = ListBuilder::new(StringBuilder::new());
        let mut b_source = ListBuilder::new(StringBuilder::new());
        let mut b_refseq_offset = ListBuilder::new(Int64Builder::new());
        let mut b_given_ref = ListBuilder::new(StringBuilder::new());
        let mut b_used_ref = ListBuilder::new(StringBuilder::new());
        let mut b_bam_edit = ListBuilder::new(StringBuilder::new());
        let mut b_gene_pheno = ListBuilder::new(StringBuilder::new());
        let mut b_sift = ListBuilder::new(StringBuilder::new());
        let mut b_polyphen = ListBuilder::new(StringBuilder::new());
        let mut b_domains = ListBuilder::new(StringBuilder::new());
        let mut b_mirna = ListBuilder::new(StringBuilder::new());
        let mut b_hgvs_offset = ListBuilder::new(Int64Builder::new());
        // Frequency (29)
        let mut b_af: Vec<Float32Builder> = (0..27)
            .map(|_| Float32Builder::with_capacity(batch.num_rows()))
            .collect();
        let mut b_max_af = Float32Builder::with_capacity(batch.num_rows());
        let mut b_max_af_pops =
            StringBuilder::with_capacity(batch.num_rows(), batch.num_rows() * 8);
        // Variant-level (9)
        let mut b_clin_sig = ListBuilder::new(StringBuilder::new());
        let mut b_somatic = StringBuilder::with_capacity(batch.num_rows(), batch.num_rows() * 4);
        let mut b_pheno = StringBuilder::with_capacity(batch.num_rows(), batch.num_rows() * 4);
        let mut b_pubmed = ListBuilder::new(StringBuilder::new());
        let mut b_motif_name = StringBuilder::with_capacity(batch.num_rows(), batch.num_rows() * 8);
        let mut b_motif_pos = StringBuilder::with_capacity(batch.num_rows(), batch.num_rows() * 4);
        let mut b_high_inf_pos =
            StringBuilder::with_capacity(batch.num_rows(), batch.num_rows() * 4);
        let mut b_motif_score_change = Float32Builder::with_capacity(batch.num_rows());
        let mut b_transcription_factors = ListBuilder::new(StringBuilder::new());
        // Cache-only (7)
        let mut b_clin_sig_allele = ListBuilder::new(StringBuilder::new());
        let mut b_clinical_impact =
            StringBuilder::with_capacity(batch.num_rows(), batch.num_rows() * 8);
        let mut b_minor_allele =
            StringBuilder::with_capacity(batch.num_rows(), batch.num_rows() * 4);
        let mut b_minor_allele_freq = Float32Builder::with_capacity(batch.num_rows());
        let mut b_clinvar_ids = ListBuilder::new(StringBuilder::new());
        let mut b_cosmic_ids = ListBuilder::new(StringBuilder::new());
        let mut b_dbsnp_ids = ListBuilder::new(StringBuilder::new());
        let include_refseq_fields = transcript_selection.refseq_fields();
        let include_source_field = transcript_selection.source_field();

        /// Append NULL to all annotation column builders for the selected mode.
        macro_rules! append_null_annotation_row {
            () => {
                b_allele.append_null();
                b_consequence.append(false);
                b_impact.append(false);
                b_symbol.append(false);
                b_gene.append(false);
                b_feature_type.append(false);
                b_feature.append(false);
                b_biotype.append(false);
                b_exon.append(false);
                b_intron.append(false);
                b_hgvsc.append(false);
                b_hgvsp.append(false);
                b_cdna_position.append(false);
                b_cds_position.append(false);
                b_protein_position.append(false);
                b_amino_acids.append(false);
                b_codons.append(false);
                b_existing_variation.append(false);
                b_distance.append(false);
                b_strand.append(false);
                b_flags.append(false);
                b_variant_class.append_null();
                b_symbol_source.append(false);
                b_hgnc_id.append(false);
                b_canonical.append(false);
                b_mane.append(false);
                b_mane_select.append(false);
                b_mane_plus_clinical.append(false);
                b_tsl.append(false);
                b_appris.append(false);
                b_ccds.append(false);
                b_ensp.append(false);
                b_swissprot.append(false);
                b_trembl.append(false);
                b_uniparc.append(false);
                b_uniprot_isoform.append(false);
                if include_refseq_fields {
                    b_refseq_match.append(false);
                    if include_source_field {
                        b_source.append(false);
                    }
                    b_refseq_offset.append(false);
                    b_given_ref.append(false);
                    b_used_ref.append(false);
                    b_bam_edit.append(false);
                }
                b_gene_pheno.append(false);
                b_sift.append(false);
                b_polyphen.append(false);
                b_domains.append(false);
                b_mirna.append(false);
                b_hgvs_offset.append(false);
                for af_b in b_af.iter_mut() {
                    af_b.append_null();
                }
                b_max_af.append_null();
                b_max_af_pops.append_null();
                b_clin_sig.append(false);
                b_somatic.append_null();
                b_pheno.append_null();
                b_pubmed.append(false);
                b_motif_name.append_null();
                b_motif_pos.append_null();
                b_high_inf_pos.append_null();
                b_motif_score_change.append_null();
                b_transcription_factors.append(false);
                b_clin_sig_allele.append(false);
                b_clinical_impact.append_null();
                b_minor_allele.append_null();
                b_minor_allele_freq.append_null();
                b_clinvar_ids.append(false);
                b_cosmic_ids.append(false);
                b_dbsnp_ids.append(false);
            };
        }

        // Reusable buffers to avoid per-row/per-CSQ-entry String allocations.
        let mut csq_buf = if skip_csq {
            String::new()
        } else {
            String::with_capacity(4096)
        };
        let mut terms_buf = String::with_capacity(128);
        // Reusable permutation index for VEP-compatible CSQ ordering.
        // Allocated once, reused across all rows in the batch.
        let mut sorted_indices: Vec<usize> = Vec::new();
        let placeholder_layout =
            CsqPlaceholderLayout::for_mode(flags.everything, transcript_selection);

        for row in 0..batch.num_rows() {
            let Some(chrom) = string_at(batch.column(chrom_idx).as_ref(), row) else {
                csq_builder.append_null();
                most_builder.append_null();
                append_null_annotation_row!();
                continue;
            };
            let Some(alt_allele) = string_at(batch.column(alt_idx).as_ref(), row) else {
                csq_builder.append_null();
                most_builder.append_null();
                append_null_annotation_row!();
                continue;
            };

            // VEP skips star alleles entirely — no CSQ produced.
            if alt_allele == "*" {
                csq_builder.append_null();
                most_builder.append_null();
                append_null_annotation_row!();
                continue;
            }

            // VEP-style allele minimization: strip shared prefix and suffix between REF and ALT.
            let ref_al = string_at(batch.column(ref_idx).as_ref(), row).unwrap_or_default();
            let (vep_ref, vep_allele) = vcf_to_vep_allele(&ref_al, &alt_allele);
            let variant_class = classify_variant(&vep_ref, &vep_allele);

            // Cache-hit fast path: use pre-computed consequence from variation cache.
            let cached_most =
                cached_most_idx.and_then(|idx| string_at(batch.column(idx).as_ref(), row));
            let cached_csq =
                cached_csq_idx.and_then(|idx| string_at(batch.column(idx).as_ref(), row));

            let _variation_name = variation_name_idx
                .and_then(|idx| string_at(batch.column(idx).as_ref(), row))
                .unwrap_or_default();

            // --- Batch 3: per-variant fields (same for every transcript entry) ---
            // Look up co-located variant aggregation (all variants at same position).
            // Traceability:
            // - Ensembl VEP `Parser::VCF::create_VariationFeatures()`
            //   https://github.com/Ensembl/ensembl-vep/blob/2beada0d57ca6234f467b14a6c60280f4a082717/modules/Bio/EnsEMBL/VEP/Parser/VCF.pm#L321-L345
            // - Ensembl VEP `compare_existing()`
            //   https://github.com/Ensembl/ensembl-vep/blob/2beada0d57ca6234f467b14a6c60280f4a082717/modules/Bio/EnsEMBL/VEP/AnnotationType/Variation.pm#L146-L206
            //
            // VEP keys the existing-variant overlap/matching flow in parser/input
            // coordinate space, not the fully minimized VEP-normalized allele space.
            let start_val = int64_at(batch.column(start_idx).as_ref(), row).unwrap_or(0);
            let end_val = int64_at(batch.column(end_idx).as_ref(), row).unwrap_or(0);
            let chrom_norm = chrom.strip_prefix("chr").unwrap_or(&chrom);
            let (input_ref, input_alt, input_start) =
                vcf_to_vep_input_allele(start_val, &ref_al, &alt_allele);
            let input_allele_string = format!("{input_ref}/{input_alt}");
            let coloc = colocated_map.get(&(
                chrom_norm.to_string(),
                input_start,
                end_val,
                input_allele_string,
            ));
            let (variant_fields, frequency_fields) = if flags.check_existing {
                if let Some(data) = coloc {
                    (
                        data.variant_fields(
                            &vep_allele,
                            data.variant_match_output_allele(&vep_allele),
                            flags.pubmed,
                        ),
                        data.frequency_fields(
                            &vep_allele,
                            data.frequency_match_output_allele(&vep_allele),
                            flags,
                        ),
                    )
                } else {
                    (
                        ColocatedVariantFields::default(),
                        ColocatedFrequencyFields {
                            af_values: vec![String::new(); AF_COLUMNS.len()],
                            max_af: String::new(),
                            max_af_pops: String::new(),
                        },
                    )
                }
            } else {
                (
                    ColocatedVariantFields::default(),
                    ColocatedFrequencyFields {
                        af_values: vec![String::new(); AF_COLUMNS.len()],
                        max_af: String::new(),
                        max_af_pops: String::new(),
                    },
                )
            };
            let existing_var = variant_fields.existing_variation.as_str();

            // Build the 33-field Batch 3 suffix (positions 41-73) shared across all transcripts.
            let batch3_suffix = format!(
                "{}|{}|{}|{}|{}|{}|{}",
                frequency_fields.af_values.join("|"),
                frequency_fields.max_af,
                frequency_fields.max_af_pops,
                variant_fields.clin_sig,
                variant_fields.somatic,
                variant_fields.pheno,
                variant_fields.pubmed,
            );

            let most_str;
            // Store assignment results from cache-miss path for annotation column population.
            let mut row_assignments: Vec<TranscriptConsequence> = Vec::new();
            // Store the VariantInput for HGVS_OFFSET extraction in annotation columns.
            let mut row_variant: Option<VariantInput> = None;
            let require_transcript_annotations =
                pick_flags.requires_transcript_annotations(skip_csq, skip_typed_cols);
            if !skip_csq {
                csq_buf.clear();
            }
            let use_cached_fast_path = cached_most.is_some() && !require_transcript_annotations;
            if use_cached_fast_path {
                use std::fmt::Write;
                let most_val = cached_most.as_deref().unwrap_or_default();
                if !skip_csq {
                    let csq_val = cached_csq.unwrap_or_default();
                    let impact = SoTerm::from_str(most_val)
                        .map(|t| impact_label(t.impact()))
                        .unwrap_or_else(|| impact_label(SoImpact::Modifier));
                    let entry = CsqPlaceholderEntry {
                        allele: &vep_allele,
                        consequence: csq_val.as_str(),
                        impact,
                        existing_variation: existing_var,
                        variant_class,
                        frequency_fields: &frequency_fields,
                        variant_fields: &variant_fields,
                    };
                    placeholder_layout.append_entry(&mut csq_buf, &entry);
                }
                most_str = most_val.to_string();
            } else {
                use std::fmt::Write;
                // Cache miss — compute via transcript engine and produce per-transcript CSQ.
                let Some(start) = int64_at(batch.column(start_idx).as_ref(), row) else {
                    csq_builder.append_null();
                    most_builder.append_null();
                    append_null_annotation_row!();
                    continue;
                };
                let Some(end) = int64_at(batch.column(end_idx).as_ref(), row) else {
                    csq_builder.append_null();
                    most_builder.append_null();
                    append_null_annotation_row!();
                    continue;
                };
                let Some(ref_allele) = string_at(batch.column(ref_idx).as_ref(), row) else {
                    csq_builder.append_null();
                    most_builder.append_null();
                    append_null_annotation_row!();
                    continue;
                };

                // VEP skips star alleles entirely — no CSQ produced.
                if alt_allele == "*" {
                    csq_builder.append_null();
                    most_builder.append_null();
                    append_null_annotation_row!();
                    continue;
                }

                let mut variant = VariantInput::from_vcf(
                    chrom.clone(),
                    start,
                    end,
                    ref_allele,
                    alt_allele.clone(),
                );
                // Only compute genomic shift for indels (ref != alt length).
                // SNVs/MNVs don't shift and skipping avoids allele normalization overhead.
                if let Some(reader) = hgvs_reference_reader.as_mut() {
                    if ref_al.len() != alt_allele.len() {
                        let chrom_norm = chrom.strip_prefix("chr").unwrap_or(&chrom);
                        let (vep_ref_norm, vep_alt_norm) = vcf_to_vep_allele(&ref_al, &alt_allele);
                        let vep_start = vep_norm_start(start, &ref_al, &alt_allele);
                        let vep_end = vep_norm_end(start, &ref_al, &alt_allele);
                        variant.hgvs_shift_forward = crate::hgvs::build_hgvs_genomic_shift(
                            reader,
                            chrom_norm,
                            &vep_ref_norm,
                            &vep_alt_norm,
                            vep_start,
                            vep_end,
                            1,
                        )?;
                        variant.hgvs_shift_reverse = crate::hgvs::build_hgvs_genomic_shift(
                            reader,
                            chrom_norm,
                            &vep_ref_norm,
                            &vep_alt_norm,
                            vep_start,
                            vep_end,
                            -1,
                        )?;
                    }
                }
                let assignments = engine.evaluate_variant_prepared(&variant, ctx);

                // Derive most_severe from all assignments.
                let mut all_terms =
                    TranscriptConsequenceEngine::collapse_variant_terms(&assignments);
                if all_terms.is_empty() {
                    all_terms.push(SoTerm::SequenceVariant);
                }
                let most = most_severe_term(all_terms.iter()).unwrap_or(SoTerm::SequenceVariant);
                most_str = most.as_str().to_string();
                row_assignments = assignments;
                row_variant = Some(variant);
                mark_flag_pick_allele_gene(&mut row_assignments, ctx, pick_flags);

                // Build VEP-compatible sorted permutation index.
                // Used by both CSQ serialization and typed annotation columns
                // so the Nth CSQ entry matches the Nth typed column element.
                //
                // Sort order: Transcript → Regulatory → Motif → Intergenic,
                // then lexicographically by feature stable_id within each
                // group. See ensembl-variation VariationFeature.pm lines 855-864.
                //
                // Source arrays are pre-sorted by feature ID in
                // PreparedContext, so transcript_idx order equals
                // lexicographic transcript_id order. We compare
                // integer indices instead of heap-allocated strings.
                // Non-transcript features (regulatory, motif) are
                // already emitted in ID order by collect_overlapping_indices.
                sorted_indices.clear();
                sorted_indices.extend(0..row_assignments.len());
                if row_assignments.len() > 1 {
                    sorted_indices.sort_unstable_by(|&i, &j| {
                        let a = &row_assignments[i];
                        let b = &row_assignments[j];
                        a.feature_type
                            .rank()
                            .cmp(&b.feature_type.rank())
                            .then_with(|| match (a.transcript_idx, b.transcript_idx) {
                                (Some(ai), Some(bj)) => ai.cmp(&bj),
                                _ => i.cmp(&j),
                            })
                    });
                }

                // Build per-transcript CSQ entries into reusable buffer (already cleared above).
                // Skip the entire CSQ formatting when the csq column is not projected.
                if !skip_csq {
                    for &si in &sorted_indices {
                        let tc = &row_assignments[si];
                        terms_buf.clear();
                        for (i, t) in tc.terms.iter().enumerate() {
                            if i > 0 {
                                terms_buf.push('&');
                            }
                            terms_buf.push_str(t.as_str());
                        }
                        let terms_str = terms_buf.as_str();
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
                        // Write comma separator between CSQ entries.
                        if !csq_buf.is_empty() {
                            csq_buf.push(',');
                        }
                        let distance = tc.distance.map(|d| d.to_string()).unwrap_or_default();
                        let tc_flags = tc.flags.as_deref().unwrap_or("");
                        let hgvsc = if hgvs_flags.hgvsc {
                            tc.hgvsc.as_deref().unwrap_or("")
                        } else {
                            ""
                        };
                        let hgvsp = if hgvs_flags.hgvsp {
                            tc.hgvsp
                                .as_deref()
                                .map(|value| {
                                    Self::format_hgvsp_output(
                                        value,
                                        hgvs_flags.remove_hgvsp_version,
                                        hgvs_flags.no_escape,
                                        hgvs_flags.hgvsp_use_prediction,
                                    )
                                })
                                .unwrap_or_default()
                        } else {
                            String::new()
                        };
                        let refseq_match = tx_opt
                            .and_then(|tx| tx.refseq_match.as_deref())
                            .unwrap_or("");
                        let bam_edit = tx_opt
                            .and_then(|tx| tx.bam_edit_status.as_deref())
                            .map(str::to_ascii_uppercase)
                            .unwrap_or_default();
                        let source_val = if include_source_field { source } else { "" };
                        let refseq_offset_value = tx_opt
                            .filter(|_| hgvs_flags.hgvsc && tc.hgvsc.is_some())
                            .and_then(|tx| {
                                tc.cdna_position
                                    .as_deref()
                                    .and_then(parse_cdna_position_start)
                                    .and_then(|cdna_start| {
                                        refseq_misalignment_offset(tx, cdna_start)
                                    })
                            });
                        let refseq_offset = refseq_offset_value
                            .map(|offset| offset.to_string())
                            .unwrap_or_default();
                        let given_ref = tc.given_ref.as_deref().unwrap_or("");
                        let used_ref = tc.used_ref.as_deref().unwrap_or("");

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
                        let ccds = tx_opt.and_then(|tx| tx.ccds.as_deref()).unwrap_or("");
                        let swissprot_raw =
                            tx_opt.and_then(|tx| tx.swissprot.as_deref()).unwrap_or("");
                        let swissprot = csq_escape(swissprot_raw);
                        let trembl_raw = tx_opt.and_then(|tx| tx.trembl.as_deref()).unwrap_or("");
                        let trembl = csq_escape(trembl_raw);
                        let uniparc = tx_opt.and_then(|tx| tx.uniparc.as_deref()).unwrap_or("");
                        let uniprot_isoform = tx_opt
                            .and_then(|tx| tx.uniprot_isoform.as_deref())
                            .unwrap_or("");

                        if flags.everything {
                            // HGVS_OFFSET mirrors the transcript-level HGVSc shift
                            // decision. RefSeq rows with failed BAM edit replay can
                            // still emit transcript-space HGVSc, but VEP suppresses
                            // the exposed genomic shift for those transcripts.
                            let hgvs_offset = if hgvs_flags.hgvsc {
                                tx_opt
                                    .zip(row_variant.as_ref())
                                    .and_then(|(tx, variant)| {
                                        let ref_allele = tc
                                            .used_ref
                                            .as_deref()
                                            .unwrap_or(variant.ref_allele.as_str());
                                        crate::hgvs::hgvsc_offset_for_output(
                                            tx,
                                            variant,
                                            ref_allele,
                                            tc.hgvsc.as_deref(),
                                        )
                                    })
                                    .map(|offset| offset.to_string())
                                    .unwrap_or_default()
                            } else {
                                String::new()
                            };
                            // MANE generic: VEP emits "MANE_Select" or "MANE_Plus_Clinical"
                            // depending on the transcript's MANE annotation.
                            // Traceability:
                            // - VEP OutputFactory.pm MANE output
                            //   https://github.com/Ensembl/ensembl-vep/blob/release/115/modules/Bio/EnsEMBL/VEP/OutputFactory.pm#L1548-L1560
                            let mane = if tx_opt.and_then(|tx| tx.mane_select.as_deref()).is_some()
                            {
                                "MANE_Select"
                            } else if tx_opt
                                .and_then(|tx| tx.mane_plus_clinical.as_deref())
                                .is_some()
                            {
                                "MANE_Plus_Clinical"
                            } else {
                                ""
                            };
                            // APPRIS: abbreviate principal1→P1, alternative2→A2.
                            // Traceability:
                            // - VEP OutputFactory.pm APPRIS output
                            //   https://github.com/Ensembl/ensembl-vep/blob/release/115/modules/Bio/EnsEMBL/VEP/OutputFactory.pm#L1563-L1570
                            let appris_str = tx_opt
                                .and_then(|tx| tx.appris.as_deref())
                                .map(format_appris)
                                .unwrap_or_default();
                            // SIFT/PolyPhen: lookup by (protein_position, alt_amino_acid).
                            // Traceability:
                            // - VEP OutputFactory.pm SIFT/PolyPhen output
                            //   https://github.com/Ensembl/ensembl-vep/blob/release/115/modules/Bio/EnsEMBL/VEP/OutputFactory.pm#L1746-L1799
                            let (sift_str, polyphen_str) = lookup_sift_polyphen(
                                tc.transcript_id.as_deref(),
                                tc.protein_position.as_deref(),
                                tc.amino_acids.as_deref(),
                                sift_cache,
                                #[cfg(feature = "kv-cache")]
                                sift_kv,
                                #[cfg(not(feature = "kv-cache"))]
                                _sift_kv,
                            );
                            // DOMAINS: overlapping protein domain features.
                            // VEP gates DOMAINS on $pre->{coding} which requires
                            // a valid CDS coordinate mapping. Our cds_position is
                            // only set when the variant falls within the CDS region.
                            // Traceability:
                            // - VEP OutputFactory.pm line 1434: if($self->{domains} && $pre->{coding})
                            // - VEP BaseVariationFeatureOverlapAllele.pm _bvfo_preds lines 449-471
                            //   https://github.com/Ensembl/ensembl-variation/blob/release/115/modules/Bio/EnsEMBL/Variation/BaseVariationFeatureOverlapAllele.pm
                            let is_coding =
                                tc.cds_position.as_deref().is_some_and(|s| !s.is_empty());
                            let domains = if is_coding {
                                lookup_domains(
                                    tc.transcript_id.as_deref(),
                                    tc.protein_position.as_deref(),
                                    tc.amino_acids.as_deref(),
                                    ctx,
                                )
                            } else {
                                String::new()
                            };
                            // miRNA: ncRNA secondary structure overlap.
                            let mirna_str = {
                                let ncrna = tx_opt.and_then(|tx| tx.ncrna_structure.as_deref());
                                // Parse cDNA position range from the "N" or "N-M" string.
                                let (cs, ce) = tc
                                    .cdna_position
                                    .as_deref()
                                    .and_then(|p| {
                                        if let Some((a, b)) = p.split_once('-') {
                                            Some((
                                                a.parse::<usize>().ok()?,
                                                b.parse::<usize>().ok()?,
                                            ))
                                        } else {
                                            let v = p.parse::<usize>().ok()?;
                                            Some((v, v))
                                        }
                                    })
                                    .unwrap_or((0, 0));
                                if cs > 0 {
                                    mirna_structure_field(ncrna, biotype, Some(cs), Some(ce))
                                } else {
                                    String::new()
                                }
                            };
                            // 80-field CSQ: 22 base + 20 Batch 1 + 33 Batch 3 + 5 motif.
                            // Traceability:
                            // - VEP Constants.pm CSQ field order for --everything
                            //   https://github.com/Ensembl/ensembl-vep/blob/release/115/modules/Bio/EnsEMBL/VEP/Constants.pm#L66-L138
                            if include_source_field {
                                let _ = write!(
                                    csq_buf,
                                    "{vep_allele}|{terms_str}|{tc_impact}|{symbol}|{gene}|{feature_type}|{feature}|{biotype}|\
                                 {exon}|{intron}|{hgvsc}|{hgvsp}|\
                                 {cdna_pos}|{cds_pos}|{protein_pos}|{amino_acids}|{codons_str}|\
                                 {existing_var}|{distance}|{strand_str}|{tc_flags}|\
                                 {variant_class}|{symbol_source}|{hgnc_id}|\
                                 {canonical}|{mane}|{mane_select}|{mane_plus}|{tsl_str}|{appris_str}|{ccds}|{ensp}|\
                                 {swissprot}|{trembl}|{uniparc}|{uniprot_isoform}|{refseq_match}|{source_val}|{refseq_offset}|{given_ref}|{used_ref}|{bam_edit}|{gene_pheno}|\
                                 {sift_str}|{polyphen_str}|{domains}|{mirna_str}|\
                                 {hgvs_offset}|\
                                 {batch3_suffix}|||||"
                                );
                            } else if include_refseq_fields {
                                let _ = write!(
                                    csq_buf,
                                    "{vep_allele}|{terms_str}|{tc_impact}|{symbol}|{gene}|{feature_type}|{feature}|{biotype}|\
                                 {exon}|{intron}|{hgvsc}|{hgvsp}|\
                                 {cdna_pos}|{cds_pos}|{protein_pos}|{amino_acids}|{codons_str}|\
                                 {existing_var}|{distance}|{strand_str}|{tc_flags}|\
                                 {variant_class}|{symbol_source}|{hgnc_id}|\
                                 {canonical}|{mane}|{mane_select}|{mane_plus}|{tsl_str}|{appris_str}|{ccds}|{ensp}|\
                                 {swissprot}|{trembl}|{uniparc}|{uniprot_isoform}|{refseq_match}|{refseq_offset}|{given_ref}|{used_ref}|{bam_edit}|{gene_pheno}|\
                                 {sift_str}|{polyphen_str}|{domains}|{mirna_str}|\
                                 {hgvs_offset}|\
                                 {batch3_suffix}|||||"
                                );
                            } else {
                                let _ = write!(
                                    csq_buf,
                                    "{vep_allele}|{terms_str}|{tc_impact}|{symbol}|{gene}|{feature_type}|{feature}|{biotype}|\
                                 {exon}|{intron}|{hgvsc}|{hgvsp}|\
                                 {cdna_pos}|{cds_pos}|{protein_pos}|{amino_acids}|{codons_str}|\
                                 {existing_var}|{distance}|{strand_str}|{tc_flags}|\
                                 {variant_class}|{symbol_source}|{hgnc_id}|\
                                 {canonical}|{mane}|{mane_select}|{mane_plus}|{tsl_str}|{appris_str}|{ccds}|{ensp}|\
                                 {swissprot}|{trembl}|{uniparc}|{uniprot_isoform}|{gene_pheno}|\
                                 {sift_str}|{polyphen_str}|{domains}|{mirna_str}|\
                                 {hgvs_offset}|\
                                 {batch3_suffix}|||||"
                                );
                            }
                        } else {
                            // 74-field CSQ: 29 base + 12 Batch 1 + 33 Batch 3.
                            if include_source_field {
                                let _ = write!(
                                    csq_buf,
                                    "{vep_allele}|{terms_str}|{tc_impact}|{symbol}|{gene}|{feature_type}|{feature}|{biotype}|\
                                 {exon}|{intron}|{hgvsc}|{hgvsp}|\
                                 {cdna_pos}|{cds_pos}|{protein_pos}|{amino_acids}|{codons_str}|\
                                 {existing_var}|{distance}|{strand_str}|{tc_flags}|{symbol_source}|{hgnc_id}|\
                                 |||||{refseq_match}|{source_val}|{refseq_offset}|{given_ref}|{used_ref}|{bam_edit}|\
                                 {variant_class}|{canonical}|{tsl_str}|{mane_select}|{mane_plus}|\
                                 {ensp}|{gene_pheno}|{ccds}|{swissprot}|{trembl}|{uniparc}|{uniprot_isoform}|\
                                 {batch3_suffix}"
                                );
                            } else if include_refseq_fields {
                                let _ = write!(
                                    csq_buf,
                                    "{vep_allele}|{terms_str}|{tc_impact}|{symbol}|{gene}|{feature_type}|{feature}|{biotype}|\
                                 {exon}|{intron}|{hgvsc}|{hgvsp}|\
                                 {cdna_pos}|{cds_pos}|{protein_pos}|{amino_acids}|{codons_str}|\
                                 {existing_var}|{distance}|{strand_str}|{tc_flags}|{symbol_source}|{hgnc_id}|\
                                 |||||{refseq_match}|{refseq_offset}|{given_ref}|{used_ref}|{bam_edit}|\
                                 {variant_class}|{canonical}|{tsl_str}|{mane_select}|{mane_plus}|\
                                 {ensp}|{gene_pheno}|{ccds}|{swissprot}|{trembl}|{uniparc}|{uniprot_isoform}|\
                                 {batch3_suffix}"
                                );
                            } else {
                                let _ = write!(
                                    csq_buf,
                                    "{vep_allele}|{terms_str}|{tc_impact}|{symbol}|{gene}|{feature_type}|{feature}|{biotype}|\
                                 {exon}|{intron}|{hgvsc}|{hgvsp}|\
                                 {cdna_pos}|{cds_pos}|{protein_pos}|{amino_acids}|{codons_str}|\
                                 {existing_var}|{distance}|{strand_str}|{tc_flags}|{symbol_source}|{hgnc_id}|\
                                 |||||{source_val}|\
                                 {variant_class}|{canonical}|{tsl_str}|{mane_select}|{mane_plus}|\
                                 {ensp}|{gene_pheno}|{ccds}|{swissprot}|{trembl}|{uniparc}|{uniprot_isoform}|\
                                 {batch3_suffix}"
                                );
                            }
                        }
                    }
                    if csq_buf.is_empty() {
                        let impact = impact_label(SoImpact::Modifier);
                        let entry = CsqPlaceholderEntry {
                            allele: &vep_allele,
                            consequence: "sequence_variant",
                            impact,
                            existing_variation: existing_var,
                            variant_class,
                            frequency_fields: &frequency_fields,
                            variant_fields: &variant_fields,
                        };
                        placeholder_layout.append_entry(&mut csq_buf, &entry);
                    }
                } // end if !skip_csq (cache-miss CSQ formatting)
            };

            if skip_csq {
                csq_builder.append_null();
            } else {
                csq_builder.append_value(&csq_buf);
            }
            most_builder.append_value(&most_str);

            // --- Populate 87 annotation column builders for this row ---
            // Skip all typed column work when they're not in the projection.
            if skip_typed_cols {
                append_null_annotation_row!();
            } else {
                // -- Transcript-level columns (42) --
                // Allele (scalar, same for all transcripts)
                b_allele.append_value(&vep_allele);
                // VARIANT_CLASS (scalar)
                b_variant_class.append_value(variant_class);
                // Existing_variation (variant-level, not per-transcript)
                if !existing_var.is_empty() {
                    let vals = b_existing_variation.values();
                    for v in existing_var.split(',') {
                        vals.append_value(v.trim());
                    }
                    b_existing_variation.append(true);
                } else {
                    b_existing_variation.append(false);
                }

                if !row_assignments.is_empty() {
                    // Cache-miss: iterate consequence entries in the same
                    // sorted order used for CSQ serialization so that the
                    // Nth typed column element matches the Nth CSQ entry.
                    for &si in &sorted_indices {
                        let tc = &row_assignments[si];
                        let tx_opt = tc.transcript_idx.map(|idx| &ctx.transcripts[idx]);
                        let refseq_match = tx_opt
                            .and_then(|tx| tx.refseq_match.as_deref())
                            .unwrap_or("");
                        let source_val = tx_opt.and_then(|tx| tx.source.as_deref()).unwrap_or("");
                        let refseq_offset_value = tx_opt
                            .filter(|_| hgvs_flags.hgvsc && tc.hgvsc.is_some())
                            .and_then(|tx| {
                                tc.cdna_position
                                    .as_deref()
                                    .and_then(parse_cdna_position_start)
                                    .and_then(|cdna_start| {
                                        refseq_misalignment_offset(tx, cdna_start)
                                    })
                            });
                        let given_ref = tc.given_ref.as_deref().unwrap_or("");
                        let used_ref = tc.used_ref.as_deref().unwrap_or("");
                        let bam_edit = tx_opt
                            .and_then(|tx| tx.bam_edit_status.as_deref())
                            .map(str::to_ascii_uppercase)
                            .unwrap_or_default();

                        // Consequence: "&"-joined terms for this transcript
                        {
                            let terms: Vec<&str> = tc.terms.iter().map(|t| t.as_str()).collect();
                            b_consequence.values().append_value(terms.join("&"));
                        }

                        // IMPACT
                        {
                            let tc_impact = most_severe_term(tc.terms.iter())
                                .map(|t| impact_label(t.impact()))
                                .unwrap_or_else(|| impact_label(SoImpact::Modifier));
                            b_impact.values().append_value(tc_impact);
                        }

                        // SYMBOL
                        append_opt_str(
                            b_symbol.values(),
                            tx_opt.and_then(|tx| tx.gene_symbol.as_deref()),
                        );

                        // Gene
                        append_opt_str(
                            b_gene.values(),
                            tx_opt.and_then(|tx| tx.gene_stable_id.as_deref()),
                        );

                        // Feature_type
                        {
                            let ft = tc.feature_type.as_str();
                            append_opt_str(
                                b_feature_type.values(),
                                if ft.is_empty() { None } else { Some(ft) },
                            );
                        }

                        // Feature
                        append_opt_str(b_feature.values(), tc.transcript_id.as_deref());

                        // BIOTYPE
                        {
                            let biotype = tc
                                .biotype_override
                                .as_deref()
                                .unwrap_or(tx_opt.map(|tx| tx.biotype.as_str()).unwrap_or(""));
                            append_opt_str(
                                b_biotype.values(),
                                if biotype.is_empty() {
                                    None
                                } else {
                                    Some(biotype)
                                },
                            );
                        }

                        // EXON, INTRON
                        append_opt_str(b_exon.values(), tc.exon_str.as_deref());
                        append_opt_str(b_intron.values(), tc.intron_str.as_deref());

                        // HGVSc
                        if hgvs_flags.hgvsc {
                            append_opt_str(b_hgvsc.values(), tc.hgvsc.as_deref());
                        } else {
                            b_hgvsc.values().append_null();
                        }

                        // HGVSp
                        if hgvs_flags.hgvsp {
                            let hgvsp_val = tc.hgvsp.as_deref().map(|value| {
                                Self::format_hgvsp_output(
                                    value,
                                    hgvs_flags.remove_hgvsp_version,
                                    hgvs_flags.no_escape,
                                    hgvs_flags.hgvsp_use_prediction,
                                )
                            });
                            match hgvsp_val {
                                Some(v) if !v.is_empty() => {
                                    b_hgvsp.values().append_value(&v);
                                }
                                _ => b_hgvsp.values().append_null(),
                            }
                        } else {
                            b_hgvsp.values().append_null();
                        }

                        // cDNA_position, CDS_position, Protein_position, Amino_acids, Codons
                        append_opt_str(b_cdna_position.values(), tc.cdna_position.as_deref());
                        append_opt_str(b_cds_position.values(), tc.cds_position.as_deref());
                        append_opt_str(b_protein_position.values(), tc.protein_position.as_deref());
                        append_opt_str(b_amino_acids.values(), tc.amino_acids.as_deref());
                        append_opt_str(b_codons.values(), tc.codons.as_deref());

                        // DISTANCE (Int64)
                        match tc.distance {
                            Some(d) => b_distance.values().append_value(d),
                            None => b_distance.values().append_null(),
                        }

                        // STRAND (Int8)
                        match tx_opt {
                            Some(tx) => {
                                b_strand
                                    .values()
                                    .append_value(if tx.strand >= 0 { 1 } else { -1 })
                            }
                            None => b_strand.values().append_null(),
                        }

                        // FLAGS
                        append_opt_str(b_flags.values(), tc.flags.as_deref());

                        // SYMBOL_SOURCE, HGNC_ID
                        append_opt_str(
                            b_symbol_source.values(),
                            tx_opt.and_then(|tx| tx.gene_symbol_source.as_deref()),
                        );
                        append_opt_str(
                            b_hgnc_id.values(),
                            tx_opt.and_then(|tx| tx.gene_hgnc_id.as_deref()),
                        );

                        // CANONICAL
                        match tx_opt {
                            Some(tx) if tx.is_canonical => {
                                b_canonical.values().append_value("YES");
                            }
                            _ => b_canonical.values().append_null(),
                        }

                        // MANE
                        if flags.everything {
                            let mane_val =
                                if tx_opt.and_then(|tx| tx.mane_select.as_deref()).is_some() {
                                    Some("MANE_Select")
                                } else if tx_opt
                                    .and_then(|tx| tx.mane_plus_clinical.as_deref())
                                    .is_some()
                                {
                                    Some("MANE_Plus_Clinical")
                                } else {
                                    None
                                };
                            append_opt_str(b_mane.values(), mane_val);
                        } else {
                            b_mane.values().append_null();
                        }

                        // MANE_SELECT, MANE_PLUS_CLINICAL
                        append_opt_str(
                            b_mane_select.values(),
                            tx_opt.and_then(|tx| tx.mane_select.as_deref()),
                        );
                        append_opt_str(
                            b_mane_plus_clinical.values(),
                            tx_opt.and_then(|tx| tx.mane_plus_clinical.as_deref()),
                        );

                        // TSL (Int8)
                        match tx_opt.and_then(|tx| tx.tsl) {
                            Some(v) => b_tsl.values().append_value(v as i8),
                            None => b_tsl.values().append_null(),
                        }

                        // APPRIS
                        if flags.everything {
                            match tx_opt.and_then(|tx| tx.appris.as_deref()) {
                                Some(raw) => {
                                    b_appris.values().append_value(format_appris(raw));
                                }
                                None => b_appris.values().append_null(),
                            }
                        } else {
                            b_appris.values().append_null();
                        }

                        // CCDS, ENSP, SWISSPROT, TREMBL, UNIPARC, UNIPROT_ISOFORM
                        append_opt_str(b_ccds.values(), tx_opt.and_then(|tx| tx.ccds.as_deref()));
                        append_opt_str(
                            b_ensp.values(),
                            tx_opt.and_then(|tx| tx.translation_stable_id.as_deref()),
                        );
                        append_opt_str(
                            b_swissprot.values(),
                            tx_opt.and_then(|tx| tx.swissprot.as_deref()),
                        );
                        append_opt_str(
                            b_trembl.values(),
                            tx_opt.and_then(|tx| tx.trembl.as_deref()),
                        );
                        append_opt_str(
                            b_uniparc.values(),
                            tx_opt.and_then(|tx| tx.uniparc.as_deref()),
                        );
                        append_opt_str(
                            b_uniprot_isoform.values(),
                            tx_opt.and_then(|tx| tx.uniprot_isoform.as_deref()),
                        );

                        if include_refseq_fields {
                            append_opt_str(b_refseq_match.values(), Some(refseq_match));
                            if include_source_field {
                                append_opt_str(b_source.values(), Some(source_val));
                            }
                            match refseq_offset_value {
                                Some(offset) => b_refseq_offset.values().append_value(offset),
                                None => b_refseq_offset.values().append_null(),
                            }
                            append_opt_str(b_given_ref.values(), Some(given_ref));
                            append_opt_str(b_used_ref.values(), Some(used_ref));
                            append_opt_str(b_bam_edit.values(), Some(bam_edit.as_str()));
                        }

                        // GENE_PHENO
                        match tx_opt {
                            Some(tx) if tx.gene_phenotype => {
                                b_gene_pheno.values().append_value("1");
                            }
                            _ => b_gene_pheno.values().append_null(),
                        }

                        // SIFT, PolyPhen
                        if flags.everything {
                            let (sift_str, polyphen_str) = lookup_sift_polyphen(
                                tc.transcript_id.as_deref(),
                                tc.protein_position.as_deref(),
                                tc.amino_acids.as_deref(),
                                sift_cache,
                                #[cfg(feature = "kv-cache")]
                                sift_kv,
                                #[cfg(not(feature = "kv-cache"))]
                                _sift_kv,
                            );
                            append_opt_str(
                                b_sift.values(),
                                if sift_str.is_empty() {
                                    None
                                } else {
                                    Some(&sift_str)
                                },
                            );
                            append_opt_str(
                                b_polyphen.values(),
                                if polyphen_str.is_empty() {
                                    None
                                } else {
                                    Some(&polyphen_str)
                                },
                            );
                        } else {
                            b_sift.values().append_null();
                            b_polyphen.values().append_null();
                        }

                        // DOMAINS (List<Utf8> -- "&"-joined string per transcript)
                        if flags.everything {
                            let is_coding =
                                tc.cds_position.as_deref().is_some_and(|s| !s.is_empty());
                            if is_coding {
                                let domains_str = lookup_domains(
                                    tc.transcript_id.as_deref(),
                                    tc.protein_position.as_deref(),
                                    tc.amino_acids.as_deref(),
                                    ctx,
                                );
                                if !domains_str.is_empty() {
                                    b_domains.values().append_value(&domains_str);
                                } else {
                                    b_domains.values().append_null();
                                }
                            } else {
                                b_domains.values().append_null();
                            }
                        } else {
                            b_domains.values().append_null();
                        }

                        // miRNA
                        if flags.everything {
                            let ncrna = tx_opt.and_then(|tx| tx.ncrna_structure.as_deref());
                            let biotype_for_mirna = tc
                                .biotype_override
                                .as_deref()
                                .unwrap_or(tx_opt.map(|tx| tx.biotype.as_str()).unwrap_or(""));
                            let (cs, ce) = tc
                                .cdna_position
                                .as_deref()
                                .and_then(|p| {
                                    if let Some((a, b)) = p.split_once('-') {
                                        Some((a.parse::<usize>().ok()?, b.parse::<usize>().ok()?))
                                    } else {
                                        let v = p.parse::<usize>().ok()?;
                                        Some((v, v))
                                    }
                                })
                                .unwrap_or((0, 0));
                            let mirna_val = if cs > 0 {
                                mirna_structure_field(ncrna, biotype_for_mirna, Some(cs), Some(ce))
                            } else {
                                String::new()
                            };
                            append_opt_str(
                                b_mirna.values(),
                                if mirna_val.is_empty() {
                                    None
                                } else {
                                    Some(&mirna_val)
                                },
                            );
                        } else {
                            b_mirna.values().append_null();
                        }

                        // HGVS_OFFSET (Int64)
                        if flags.everything && hgvs_flags.hgvsc && tc.hgvsc.is_some() {
                            let tx_strand = tx_opt.map(|tx| tx.strand).unwrap_or(1);
                            let offset_val = row_variant
                                .as_ref()
                                .and_then(|v| v.hgvs_shift_for_strand(tx_strand))
                                .filter(|s| s.shift_length > 0)
                                .map(|s| {
                                    let signed = s.shift_length as i64;
                                    if tx_strand < 0 { -signed } else { signed }
                                });
                            match offset_val {
                                Some(v) => b_hgvs_offset.values().append_value(v),
                                None => b_hgvs_offset.values().append_null(),
                            }
                        } else {
                            b_hgvs_offset.values().append_null();
                        }
                    } // end for tc in row_assignments

                    // Close all outer lists
                    b_consequence.append(true);
                    b_impact.append(true);
                    b_symbol.append(true);
                    b_gene.append(true);
                    b_feature_type.append(true);
                    b_feature.append(true);
                    b_biotype.append(true);
                    b_exon.append(true);
                    b_intron.append(true);
                    b_hgvsc.append(true);
                    b_hgvsp.append(true);
                    b_cdna_position.append(true);
                    b_cds_position.append(true);
                    b_protein_position.append(true);
                    b_amino_acids.append(true);
                    b_codons.append(true);
                    b_distance.append(true);
                    b_strand.append(true);
                    b_flags.append(true);
                    b_symbol_source.append(true);
                    b_hgnc_id.append(true);
                    b_canonical.append(true);
                    b_mane.append(true);
                    b_mane_select.append(true);
                    b_mane_plus_clinical.append(true);
                    b_tsl.append(true);
                    b_appris.append(true);
                    b_ccds.append(true);
                    b_ensp.append(true);
                    b_swissprot.append(true);
                    b_trembl.append(true);
                    b_uniparc.append(true);
                    b_uniprot_isoform.append(true);
                    if include_refseq_fields {
                        b_refseq_match.append(true);
                        if include_source_field {
                            b_source.append(true);
                        }
                        b_refseq_offset.append(true);
                        b_given_ref.append(true);
                        b_used_ref.append(true);
                        b_bam_edit.append(true);
                    }
                    b_gene_pheno.append(true);
                    b_sift.append(true);
                    b_polyphen.append(true);
                    b_domains.append(true);
                    b_mirna.append(true);
                    b_hgvs_offset.append(true);
                } else {
                    // Cache-hit path: NULL lists for all transcript-level columns
                    b_consequence.append(false);
                    b_impact.append(false);
                    b_symbol.append(false);
                    b_gene.append(false);
                    b_feature_type.append(false);
                    b_feature.append(false);
                    b_biotype.append(false);
                    b_exon.append(false);
                    b_intron.append(false);
                    b_hgvsc.append(false);
                    b_hgvsp.append(false);
                    b_cdna_position.append(false);
                    b_cds_position.append(false);
                    b_protein_position.append(false);
                    b_amino_acids.append(false);
                    b_codons.append(false);
                    b_distance.append(false);
                    b_strand.append(false);
                    b_flags.append(false);
                    b_symbol_source.append(false);
                    b_hgnc_id.append(false);
                    b_canonical.append(false);
                    b_mane.append(false);
                    b_mane_select.append(false);
                    b_mane_plus_clinical.append(false);
                    b_tsl.append(false);
                    b_appris.append(false);
                    b_ccds.append(false);
                    b_ensp.append(false);
                    b_swissprot.append(false);
                    b_trembl.append(false);
                    b_uniparc.append(false);
                    b_uniprot_isoform.append(false);
                    if include_refseq_fields {
                        b_refseq_match.append(false);
                        if include_source_field {
                            b_source.append(false);
                        }
                        b_refseq_offset.append(false);
                        b_given_ref.append(false);
                        b_used_ref.append(false);
                        b_bam_edit.append(false);
                    }
                    b_gene_pheno.append(false);
                    b_sift.append(false);
                    b_polyphen.append(false);
                    b_domains.append(false);
                    b_mirna.append(false);
                    b_hgvs_offset.append(false);
                }

                // -- Frequency columns (29) --
                // AF columns: parse resolved frequency strings to Float32.
                for (i, af_val) in frequency_fields.af_values.iter().enumerate() {
                    if i < b_af.len() {
                        if af_val.is_empty() {
                            b_af[i].append_null();
                        } else {
                            match af_val.parse::<f32>() {
                                Ok(v) => b_af[i].append_value(v),
                                Err(_) => b_af[i].append_null(),
                            }
                        }
                    }
                }
                // Pad any missing AF columns (if frequency_fields has fewer than 27 entries).
                for i in frequency_fields.af_values.len()..b_af.len() {
                    b_af[i].append_null();
                }
                // MAX_AF
                if frequency_fields.max_af.is_empty() {
                    b_max_af.append_null();
                } else {
                    match frequency_fields.max_af.parse::<f32>() {
                        Ok(v) => b_max_af.append_value(v),
                        Err(_) => b_max_af.append_null(),
                    }
                }
                // MAX_AF_POPS
                append_opt_str(
                    &mut b_max_af_pops,
                    if frequency_fields.max_af_pops.is_empty() {
                        None
                    } else {
                        Some(&frequency_fields.max_af_pops)
                    },
                );

                // -- Variant-level columns (9) --
                // CLIN_SIG (List<Utf8>)
                if !variant_fields.clin_sig.is_empty() {
                    let vals = b_clin_sig.values();
                    for v in variant_fields.clin_sig.split(',') {
                        vals.append_value(v.trim());
                    }
                    b_clin_sig.append(true);
                } else {
                    b_clin_sig.append(false);
                }
                // SOMATIC
                append_opt_str(
                    &mut b_somatic,
                    if variant_fields.somatic.is_empty() {
                        None
                    } else {
                        Some(&variant_fields.somatic)
                    },
                );
                // PHENO
                append_opt_str(
                    &mut b_pheno,
                    if variant_fields.pheno.is_empty() {
                        None
                    } else {
                        Some(&variant_fields.pheno)
                    },
                );
                // PUBMED (List<Utf8>)
                if !variant_fields.pubmed.is_empty() {
                    let vals = b_pubmed.values();
                    for v in variant_fields.pubmed.split(',') {
                        vals.append_value(v.trim());
                    }
                    b_pubmed.append(true);
                } else {
                    b_pubmed.append(false);
                }
                // MOTIF_NAME, MOTIF_POS, HIGH_INF_POS, MOTIF_SCORE_CHANGE, TRANSCRIPTION_FACTORS
                // These are currently not populated (always NULL) — they require motif feature
                // consequence data that is not yet exposed in the per-transcript CSQ path.
                b_motif_name.append_null();
                b_motif_pos.append_null();
                b_high_inf_pos.append_null();
                b_motif_score_change.append_null();
                b_transcription_factors.append(false);

                // -- Cache-only columns (7) --
                // Read from the intermediate batch's cache_ columns.
                // clin_sig_allele (List<Utf8>, semicolon-separated)
                match schema
                    .index_of("cache_clin_sig_allele")
                    .ok()
                    .and_then(|idx| string_at(batch.column(idx).as_ref(), row))
                {
                    Some(v) if !v.is_empty() => {
                        let vals = b_clin_sig_allele.values();
                        for s in v.split(';') {
                            vals.append_value(s.trim());
                        }
                        b_clin_sig_allele.append(true);
                    }
                    _ => b_clin_sig_allele.append(false),
                }
                // clinical_impact
                match schema
                    .index_of("cache_clinical_impact")
                    .ok()
                    .and_then(|idx| string_at(batch.column(idx).as_ref(), row))
                {
                    Some(v) if !v.is_empty() => b_clinical_impact.append_value(&v),
                    _ => b_clinical_impact.append_null(),
                }
                // minor_allele
                match schema
                    .index_of("cache_minor_allele")
                    .ok()
                    .and_then(|idx| string_at(batch.column(idx).as_ref(), row))
                {
                    Some(v) if !v.is_empty() => b_minor_allele.append_value(&v),
                    _ => b_minor_allele.append_null(),
                }
                // minor_allele_freq (Float32)
                match schema
                    .index_of("cache_minor_allele_freq")
                    .ok()
                    .and_then(|idx| string_at(batch.column(idx).as_ref(), row))
                {
                    Some(v) if !v.is_empty() => match v.parse::<f32>() {
                        Ok(f) => b_minor_allele_freq.append_value(f),
                        Err(_) => b_minor_allele_freq.append_null(),
                    },
                    _ => b_minor_allele_freq.append_null(),
                }
                // clinvar_ids (List<Utf8>)
                match schema
                    .index_of("cache_clinvar_ids")
                    .ok()
                    .and_then(|idx| string_at(batch.column(idx).as_ref(), row))
                {
                    Some(v) if !v.is_empty() => {
                        let vals = b_clinvar_ids.values();
                        for s in v.split(',') {
                            vals.append_value(s.trim());
                        }
                        b_clinvar_ids.append(true);
                    }
                    _ => b_clinvar_ids.append(false),
                }
                // cosmic_ids (List<Utf8>)
                match schema
                    .index_of("cache_cosmic_ids")
                    .ok()
                    .and_then(|idx| string_at(batch.column(idx).as_ref(), row))
                {
                    Some(v) if !v.is_empty() => {
                        let vals = b_cosmic_ids.values();
                        for s in v.split(',') {
                            vals.append_value(s.trim());
                        }
                        b_cosmic_ids.append(true);
                    }
                    _ => b_cosmic_ids.append(false),
                }
                // dbsnp_ids (List<Utf8>)
                match schema
                    .index_of("cache_dbsnp_ids")
                    .ok()
                    .and_then(|idx| string_at(batch.column(idx).as_ref(), row))
                {
                    Some(v) if !v.is_empty() => {
                        let vals = b_dbsnp_ids.values();
                        for s in v.split(',') {
                            vals.append_value(s.trim());
                        }
                        b_dbsnp_ids.append(true);
                    }
                    _ => b_dbsnp_ids.append(false),
                }
            } // end if !skip_typed_cols
        } // end per-row loop

        // --- Build output columns ---
        let mut out_cols =
            Vec::with_capacity(self.vcf_field_count() + self.annotation_column_count());
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

        // Typed annotation columns are mode-dependent; skip building when not projected.
        if skip_typed_cols {
            for col_def in &self.annotation_column_defs {
                out_cols.push(new_null_array(&col_def.data_type, batch.num_rows()));
            }
        } else {
            // Transcript-level columns (42)
            out_cols.push(Arc::new(b_allele.finish()));
            out_cols.push(Arc::new(b_consequence.finish()));
            out_cols.push(Arc::new(b_impact.finish()));
            out_cols.push(Arc::new(b_symbol.finish()));
            out_cols.push(Arc::new(b_gene.finish()));
            out_cols.push(Arc::new(b_feature_type.finish()));
            out_cols.push(Arc::new(b_feature.finish()));
            out_cols.push(Arc::new(b_biotype.finish()));
            out_cols.push(Arc::new(b_exon.finish()));
            out_cols.push(Arc::new(b_intron.finish()));
            out_cols.push(Arc::new(b_hgvsc.finish()));
            out_cols.push(Arc::new(b_hgvsp.finish()));
            out_cols.push(Arc::new(b_cdna_position.finish()));
            out_cols.push(Arc::new(b_cds_position.finish()));
            out_cols.push(Arc::new(b_protein_position.finish()));
            out_cols.push(Arc::new(b_amino_acids.finish()));
            out_cols.push(Arc::new(b_codons.finish()));
            out_cols.push(Arc::new(b_existing_variation.finish()));
            out_cols.push(Arc::new(b_distance.finish()));
            out_cols.push(Arc::new(b_strand.finish()));
            out_cols.push(Arc::new(b_flags.finish()));
            out_cols.push(Arc::new(b_variant_class.finish()));
            out_cols.push(Arc::new(b_symbol_source.finish()));
            out_cols.push(Arc::new(b_hgnc_id.finish()));
            out_cols.push(Arc::new(b_canonical.finish()));
            out_cols.push(Arc::new(b_mane.finish()));
            out_cols.push(Arc::new(b_mane_select.finish()));
            out_cols.push(Arc::new(b_mane_plus_clinical.finish()));
            out_cols.push(Arc::new(b_tsl.finish()));
            out_cols.push(Arc::new(b_appris.finish()));
            out_cols.push(Arc::new(b_ccds.finish()));
            out_cols.push(Arc::new(b_ensp.finish()));
            out_cols.push(Arc::new(b_swissprot.finish()));
            out_cols.push(Arc::new(b_trembl.finish()));
            out_cols.push(Arc::new(b_uniparc.finish()));
            out_cols.push(Arc::new(b_uniprot_isoform.finish()));
            if include_refseq_fields {
                out_cols.push(Arc::new(b_refseq_match.finish()));
                if include_source_field {
                    out_cols.push(Arc::new(b_source.finish()));
                }
                out_cols.push(Arc::new(b_refseq_offset.finish()));
                out_cols.push(Arc::new(b_given_ref.finish()));
                out_cols.push(Arc::new(b_used_ref.finish()));
                out_cols.push(Arc::new(b_bam_edit.finish()));
            }
            out_cols.push(Arc::new(b_gene_pheno.finish()));
            out_cols.push(Arc::new(b_sift.finish()));
            out_cols.push(Arc::new(b_polyphen.finish()));
            out_cols.push(Arc::new(b_domains.finish()));
            out_cols.push(Arc::new(b_mirna.finish()));
            out_cols.push(Arc::new(b_hgvs_offset.finish()));
            // Frequency columns (29)
            for af_b in b_af.iter_mut() {
                out_cols.push(Arc::new(af_b.finish()));
            }
            out_cols.push(Arc::new(b_max_af.finish()));
            out_cols.push(Arc::new(b_max_af_pops.finish()));
            // Variant-level columns (9)
            out_cols.push(Arc::new(b_clin_sig.finish()));
            out_cols.push(Arc::new(b_somatic.finish()));
            out_cols.push(Arc::new(b_pheno.finish()));
            out_cols.push(Arc::new(b_pubmed.finish()));
            out_cols.push(Arc::new(b_motif_name.finish()));
            out_cols.push(Arc::new(b_motif_pos.finish()));
            out_cols.push(Arc::new(b_high_inf_pos.finish()));
            out_cols.push(Arc::new(b_motif_score_change.finish()));
            out_cols.push(Arc::new(b_transcription_factors.finish()));
            // Cache-only columns (7)
            out_cols.push(Arc::new(b_clin_sig_allele.finish()));
            out_cols.push(Arc::new(b_clinical_impact.finish()));
            out_cols.push(Arc::new(b_minor_allele.finish()));
            out_cols.push(Arc::new(b_minor_allele_freq.finish()));
            out_cols.push(Arc::new(b_clinvar_ids.finish()));
            out_cols.push(Arc::new(b_cosmic_ids.finish()));
            out_cols.push(Arc::new(b_dbsnp_ids.finish()));
        }

        debug_assert_eq!(
            out_cols.len(),
            self.schema.fields().len(),
            "annotate_vep(): output column builder order is out of sync with provider schema"
        );

        Ok(RecordBatch::try_new(self.schema.clone(), out_cols)?)
    }
}

/// Append an optional string value to a StringBuilder: non-empty Some → value, else → NULL.
fn append_opt_str(builder: &mut StringBuilder, val: Option<&str>) {
    match val {
        Some(v) if !v.is_empty() => builder.append_value(v),
        _ => builder.append_null(),
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

/// Reconstruct `FLAGS` string from promoted boolean columns when the ordered
/// transcript attributes are unavailable in `raw_object_json`.
fn flags_str_from_bools(cds_start_nf: bool, cds_end_nf: bool) -> Option<String> {
    match (cds_start_nf, cds_end_nf) {
        (true, true) => Some("cds_start_NF&cds_end_NF".to_string()),
        (true, false) => Some("cds_start_NF".to_string()),
        (false, true) => Some("cds_end_NF".to_string()),
        (false, false) => None,
    }
}

fn normalize_source_label(source: &str) -> Option<String> {
    match source {
        "" | "-" => None,
        "Ensembl" | "RefSeq" => Some(source.to_string()),
        value
            if matches!(
                value,
                "ensembl" | "ensembl_havana" | "havana" | "VEGA" | "vega"
            ) =>
        {
            Some("Ensembl".to_string())
        }
        value if matches!(value, "BestRefSeq" | "RefSeq" | "Gnomon") => Some("RefSeq".to_string()),
        other => Some(other.to_string()),
    }
}

fn json_unwrap_value(value: &Value) -> &Value {
    value.get("__value").unwrap_or(value)
}

fn json_extract_seq(value: &Value) -> Option<String> {
    let value = json_unwrap_value(value);
    value
        .get("seq")
        .and_then(Value::as_str)
        .or_else(|| {
            value
                .get("primary_seq")
                .map(json_unwrap_value)
                .and_then(|seq| seq.get("seq"))
                .and_then(Value::as_str)
        })
        .map(|seq| seq.to_ascii_uppercase())
}

fn synthesize_spliced_seq(
    five_prime_utr_seq: Option<&str>,
    translateable_seq: Option<&str>,
    three_prime_utr_seq: Option<&str>,
    cdna_coding_start: Option<usize>,
    cdna_coding_end: Option<usize>,
    cdna_seq: Option<&str>,
) -> Option<String> {
    let translateable_seq = translateable_seq.or_else(|| {
        let cdna_seq = cdna_seq?;
        let start = cdna_coding_start?.checked_sub(1)?;
        let end = cdna_coding_end?;
        if start >= end || end > cdna_seq.len() {
            return None;
        }
        Some(&cdna_seq[start..end])
    })?;

    let mut spliced_seq = String::with_capacity(
        five_prime_utr_seq.map_or(0, str::len)
            + translateable_seq.len()
            + three_prime_utr_seq.map_or(0, str::len),
    );
    if let Some(seq) = five_prime_utr_seq {
        spliced_seq.push_str(seq);
    }
    spliced_seq.push_str(translateable_seq);
    if let Some(seq) = three_prime_utr_seq {
        spliced_seq.push_str(seq);
    }
    (!spliced_seq.is_empty()).then_some(spliced_seq.to_ascii_uppercase())
}

fn push_unique_string(out: &mut Vec<String>, seen: &mut HashSet<String>, value: &str) {
    if seen.insert(value.to_string()) {
        out.push(value.to_string());
    }
}

fn parse_refseq_edit_attribute(attribute: &Value) -> Option<RefSeqEdit> {
    let value = attribute.get("value").and_then(Value::as_str)?;
    let parts: Vec<&str> = value.split_whitespace().collect();
    if !matches!(parts.len(), 2 | 3) {
        return None;
    }

    let start = parts[0].parse::<i64>().ok()?;
    let end = parts[1].parse::<i64>().ok()?;
    let replacement_len = (parts.len() == 3).then(|| parts[2].len());
    let same_len_substitution = replacement_len.is_some_and(|len| end - start + 1 == len as i64);
    let op_x_edit = attribute
        .get("description")
        .and_then(Value::as_str)
        .is_some_and(|description| description.contains("op=X"));

    Some(RefSeqEdit {
        start,
        end,
        replacement_len,
        skip_refseq_offset: same_len_substitution || op_x_edit,
    })
}

fn parse_raw_cdna_mapper_segments(vef_cache: Option<&Value>) -> Vec<TranscriptCdnaMapperSegment> {
    let Some(pairs) = vef_cache
        .and_then(|cache| cache.get("mapper"))
        .map(json_unwrap_value)
        .and_then(|mapper| mapper.get("exon_coord_mapper"))
        .map(json_unwrap_value)
        .and_then(|mapper| mapper.get("_pair_cdna"))
        .map(json_unwrap_value)
        .and_then(|pairs| pairs.get("CDNA"))
        .and_then(Value::as_array)
    else {
        return Vec::new();
    };

    let mut segments = Vec::with_capacity(pairs.len());
    for pair in pairs {
        let pair = json_unwrap_value(pair);
        let Some(from) = pair.get("from").map(json_unwrap_value) else {
            continue;
        };
        let Some(to) = pair.get("to").map(json_unwrap_value) else {
            continue;
        };
        let Some(ori) = pair
            .get("ori")
            .and_then(Value::as_i64)
            .and_then(|v| i8::try_from(v).ok())
        else {
            continue;
        };
        let Some(cdna_start) = from
            .get("start")
            .and_then(Value::as_i64)
            .and_then(|v| usize::try_from(v).ok())
        else {
            continue;
        };
        let Some(cdna_end) = from
            .get("end")
            .and_then(Value::as_i64)
            .and_then(|v| usize::try_from(v).ok())
        else {
            continue;
        };
        let Some(genomic_start) = to.get("start").and_then(Value::as_i64) else {
            continue;
        };
        let Some(genomic_end) = to.get("end").and_then(Value::as_i64) else {
            continue;
        };
        segments.push(TranscriptCdnaMapperSegment {
            genomic_start,
            genomic_end,
            cdna_start,
            cdna_end,
            ori,
        });
    }
    segments.sort_by_key(|segment| {
        (
            segment.genomic_start,
            segment.genomic_end,
            segment.cdna_start,
        )
    });
    // Ensembl TranscriptMapper can encode transcript-only insertions or other
    // complex gap semantics as multiple adjacent genomic pairs with cDNA jumps.
    // Our current fallback only replays simple monotonic pair mappings; if the
    // serialized mapper contains an internal cDNA discontinuity across
    // contiguous genomic bases, keep using the exon-based fallback until we
    // implement full Mapper gap semantics.
    // Traceability:
    // - Ensembl `Bio::EnsEMBL::Mapper` stores both Pair and Gap units
    // - VEP reuses the live TranscriptMapper via `genomic2cdna`
    //   https://github.com/Ensembl/ensembl-variation/blob/release/115/modules/Bio/EnsEMBL/Variation/BaseTranscriptVariation.pm#L478-L492
    segments
}

fn parse_transcript_raw_metadata(raw_object_json: &str) -> TranscriptRawMetadata {
    let Ok(root) = serde_json::from_str::<Value>(raw_object_json) else {
        return TranscriptRawMetadata::default();
    };
    let tx = json_unwrap_value(&root);
    let vef_cache = tx
        .get("_variation_effect_feature_cache")
        .map(json_unwrap_value);
    let display_xref_id = tx
        .get("display_xref")
        .map(json_unwrap_value)
        .and_then(|xref| xref.get("display_id"))
        .and_then(Value::as_str)
        .map(str::to_string);
    let source = tx
        .get("_source_cache")
        .and_then(Value::as_str)
        .and_then(normalize_source_label);
    let cdna_mapper_segments = parse_raw_cdna_mapper_segments(vef_cache);
    let gene_hgnc_id_native = tx
        .get("_gene_hgnc_id")
        .or_else(|| tx.get("gene_hgnc_id"))
        .and_then(Value::as_str)
        .map(str::to_string);
    let spliced_seq = tx.get("spliced_seq").and_then(json_extract_seq);
    let five_prime_utr_seq = vef_cache
        .and_then(|cache| cache.get("five_prime_utr"))
        .and_then(json_extract_seq);
    let three_prime_utr_seq = tx
        .get("three_prime_utr")
        .or_else(|| vef_cache.and_then(|cache| cache.get("three_prime_utr")))
        .and_then(json_extract_seq);
    let translateable_seq = tx
        .get("translateable_seq")
        .or_else(|| vef_cache.and_then(|cache| cache.get("translateable_seq")))
        .and_then(|value| match value {
            Value::String(seq) => Some(seq.to_ascii_uppercase()),
            _ => json_extract_seq(value),
        });

    let mut refseq_match_codes = Vec::new();
    let mut seen_refseq_match_codes = HashSet::new();
    let mut flags = Vec::new();
    let mut seen_flags = HashSet::new();
    let mut refseq_edits = Vec::new();
    let mut is_gencode_basic = false;
    let mut is_gencode_primary = false;

    if let Some(attributes) = tx.get("attributes").and_then(Value::as_array) {
        for attribute in attributes {
            let attribute = json_unwrap_value(attribute);
            let Some(code) = attribute.get("code").and_then(Value::as_str) else {
                continue;
            };
            match code {
                "gencode_basic" => is_gencode_basic = true,
                "gencode_primary" => is_gencode_primary = true,
                "cds_start_NF" | "cds_end_NF" => {
                    push_unique_string(&mut flags, &mut seen_flags, code)
                }
                _ if code.starts_with("rseq") => {
                    push_unique_string(&mut refseq_match_codes, &mut seen_refseq_match_codes, code);
                }
                _ if code.starts_with("_rna_edit") => {
                    if let Some(edit) = parse_refseq_edit_attribute(attribute) {
                        refseq_edits.push(edit);
                    }
                }
                _ => {}
            }
        }
    }
    refseq_edits.sort_by(|left, right| {
        left.start
            .cmp(&right.start)
            .then(left.end.cmp(&right.end))
            .then(left.replacement_len.cmp(&right.replacement_len))
    });

    TranscriptRawMetadata {
        display_xref_id,
        source,
        gene_hgnc_id_native,
        refseq_match: (!refseq_match_codes.is_empty()).then(|| refseq_match_codes.join("&")),
        refseq_edits,
        cdna_mapper_segments,
        spliced_seq,
        five_prime_utr_seq,
        three_prime_utr_seq,
        translateable_seq,
        flags_str: (!flags.is_empty()).then(|| flags.join("&")),
        is_gencode_basic,
        is_gencode_primary,
    }
}

fn is_ensembl_transcript(tx: &TranscriptFeature) -> bool {
    tx.source.as_deref() == Some("Ensembl") || tx.transcript_id.starts_with("ENST")
}

fn is_refseq_transcript(tx: &TranscriptFeature) -> bool {
    tx.source.as_deref() == Some("RefSeq")
        || matches!(
            tx.transcript_id.as_bytes().get(..2),
            Some(b"NM") | Some(b"NR") | Some(b"XM") | Some(b"XR")
        )
}

fn is_refseq_offset_transcript(tx: &TranscriptFeature) -> bool {
    tx.transcript_id.starts_with("NM_") || tx.transcript_id.starts_with("XM_")
}

fn parse_cdna_position_start(value: &str) -> Option<i64> {
    let value = value.trim();
    let mut chars = value.char_indices();
    let mut sign = 1_i64;
    let mut digit_start = 0;

    if let Some((_, '-')) = chars.clone().next() {
        sign = -1;
        digit_start = 1;
        chars.next();
    }

    let mut digit_end = digit_start;
    for (idx, ch) in value[digit_start..].char_indices() {
        if ch.is_ascii_digit() {
            digit_end = digit_start + idx + ch.len_utf8();
        } else {
            break;
        }
    }
    if digit_end == digit_start {
        return None;
    }

    value[digit_start..digit_end]
        .parse::<i64>()
        .ok()
        .map(|pos| pos * sign)
}

fn refseq_misalignment_offset(tx: &TranscriptFeature, cdna_start: i64) -> Option<i64> {
    if !is_refseq_offset_transcript(tx) || tx.refseq_edits.is_empty() {
        return None;
    }

    let mut offset = 0_i64;
    for edit in &tx.refseq_edits {
        if edit.end >= cdna_start {
            continue;
        }
        offset += refseq_edit_offset_delta(edit).unwrap_or(0);
    }

    (offset != 0).then_some(offset)
}

fn is_predicted_refseq_transcript(tx: &TranscriptFeature) -> bool {
    tx.transcript_id.starts_with("XM_") || tx.transcript_id.starts_with("XR_")
}

fn is_mitochondrial_chrom(chrom: &str) -> bool {
    matches!(
        chrom.strip_prefix("chr").unwrap_or(chrom),
        "M" | "MT" | "m" | "mt"
    )
}

fn is_default_refseq_transcript_id(tx: &TranscriptFeature) -> bool {
    let id = tx.transcript_id.as_str();
    let starts_with_refseq_accession = id.len() >= 4
        && id.as_bytes()[0].is_ascii_uppercase()
        && id.as_bytes()[1].is_ascii_uppercase()
        && id.as_bytes()[2] == b'_'
        && id.as_bytes()[3].is_ascii_digit();
    if starts_with_refseq_accession {
        return true;
    }

    if is_mitochondrial_chrom(&tx.chrom) {
        let is_mt_stable_id = (id.len() == 4 && id.chars().all(|ch| ch.is_ascii_digit()))
            || id
                .strip_prefix("rna-")
                .unwrap_or(id)
                .chars()
                .all(|ch| ch.is_ascii_uppercase() || ch.is_ascii_digit())
                && id.strip_prefix("rna-").unwrap_or(id).len() >= 3;
        if is_mt_stable_id {
            return true;
        }
    }

    tx.display_xref_id.as_deref().is_some_and(|display_id| {
        let is_refseq_accession = display_id.len() >= 4
            && display_id.as_bytes()[0].is_ascii_uppercase()
            && display_id.as_bytes()[1].is_ascii_uppercase()
            && display_id.as_bytes()[2] == b'_'
            && display_id.as_bytes()[3].is_ascii_digit();
        let is_mt_display_id =
            display_id.len() == 4 && display_id.chars().all(|ch| ch.is_ascii_digit());
        is_refseq_accession || is_mt_display_id
    })
}

fn passes_transcript_selection(
    tx: &TranscriptFeature,
    selection: TranscriptSelectionFlags,
) -> bool {
    if tx.transcript_id.is_empty() {
        return false;
    }

    if selection.gencode_basic && !tx.is_gencode_basic {
        return false;
    }
    if selection.gencode_primary && !tx.is_gencode_primary {
        return false;
    }
    if selection.exclude_predicted && is_predicted_refseq_transcript(tx) {
        return false;
    }

    match selection.source_mode {
        TranscriptSourceMode::Ensembl => is_ensembl_transcript(tx),
        TranscriptSourceMode::RefSeq => {
            is_refseq_transcript(tx)
                && (selection.all_refseq || is_default_refseq_transcript_id(tx))
        }
        TranscriptSourceMode::Merged => {
            if is_refseq_transcript(tx) {
                selection.all_refseq || is_default_refseq_transcript_id(tx)
            } else {
                is_ensembl_transcript(tx)
            }
        }
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

/// Read predictions directly into CompactPrediction without intermediate String allocations.
/// Reads `&str` from Arrow arrays and encodes amino acid/prediction as u8 in-place.
pub(crate) fn read_compact_predictions(col: &dyn Array, row: usize) -> Vec<CompactPrediction> {
    if col.is_null(row) {
        return Vec::new();
    }
    let Some(list_arr) = col.as_any().downcast_ref::<ListArray>() else {
        return Vec::new();
    };
    let offsets = list_arr.offsets();
    let start_off = offsets[row] as usize;
    let end_off = offsets[row + 1] as usize;
    if start_off == end_off {
        return Vec::new();
    }
    let values = list_arr.values();
    let struct_arr = values.as_struct();
    let Some(positions) = struct_arr.column_by_name("position") else {
        return Vec::new();
    };
    let Some(amino_acids) = struct_arr.column_by_name("amino_acid") else {
        return Vec::new();
    };
    let Some(predictions) = struct_arr.column_by_name("prediction") else {
        return Vec::new();
    };
    let Some(scores) = struct_arr.column_by_name("score") else {
        return Vec::new();
    };
    let pos_arr = positions.as_any().downcast_ref::<Int32Array>();
    let aa_arr = amino_acids.as_any().downcast_ref::<StringArray>();
    let aa_view = amino_acids.as_any().downcast_ref::<StringViewArray>();
    let pred_arr = predictions.as_any().downcast_ref::<StringArray>();
    let pred_view = predictions.as_any().downcast_ref::<StringViewArray>();
    let score_arr = scores.as_any().downcast_ref::<Float32Array>();

    let mut out = Vec::with_capacity(end_off - start_off);
    for i in start_off..end_off {
        let Some(pos) = pos_arr.and_then(|a| if a.is_null(i) { None } else { Some(a.value(i)) })
        else {
            continue;
        };
        // Read &str directly from Arrow buffer (zero-copy), encode to u8
        let aa_str = aa_arr
            .and_then(|a| if a.is_null(i) { None } else { Some(a.value(i)) })
            .or_else(|| aa_view.and_then(|a| if a.is_null(i) { None } else { Some(a.value(i)) }));
        let pred_str = pred_arr
            .and_then(|a| if a.is_null(i) { None } else { Some(a.value(i)) })
            .or_else(|| pred_view.and_then(|a| if a.is_null(i) { None } else { Some(a.value(i)) }));
        let score = score_arr.and_then(|a| if a.is_null(i) { None } else { Some(a.value(i)) });
        if let (Some(aa), Some(pred), Some(sc)) = (aa_str, pred_str, score) {
            if let Some(aa_idx) = CompactPrediction::encode_amino_acid(aa) {
                out.push(CompactPrediction {
                    position: pos,
                    amino_acid: aa_idx,
                    prediction: CompactPrediction::encode_prediction(pred),
                    score: sc,
                });
            }
        }
    }
    out
}

/// Read protein domain features from a `List<Struct<analysis, hseqname, start, end>>` column.
///
/// Traceability:
/// - VEP OutputFactory.pm DOMAINS output
///   <https://github.com/Ensembl/ensembl-vep/blob/release/115/modules/Bio/EnsEMBL/VEP/OutputFactory.pm#L1448-L1466>
fn read_protein_features(col: &dyn Array, row: usize) -> Vec<ProteinDomainFeature> {
    if col.is_null(row) {
        return Vec::new();
    }
    let Some(list_arr) = col.as_any().downcast_ref::<ListArray>() else {
        return Vec::new();
    };
    let offsets = list_arr.offsets();
    let start_off = offsets[row] as usize;
    let end_off = offsets[row + 1] as usize;
    if start_off == end_off {
        return Vec::new();
    }
    let values = list_arr.values();
    let struct_arr = values.as_struct();
    let analysis_col = struct_arr.column_by_name("analysis");
    let hseqname_col = struct_arr.column_by_name("hseqname");
    let Some(start_col) = struct_arr.column_by_name("start") else {
        return Vec::new();
    };
    let Some(end_col) = struct_arr.column_by_name("end") else {
        return Vec::new();
    };

    let start_arr = start_col.as_any().downcast_ref::<Int64Array>();
    let end_arr = end_col.as_any().downcast_ref::<Int64Array>();

    let mut out = Vec::with_capacity(end_off - start_off);
    for i in start_off..end_off {
        let s = start_arr.and_then(|a| if a.is_null(i) { None } else { Some(a.value(i)) });
        let e = end_arr.and_then(|a| if a.is_null(i) { None } else { Some(a.value(i)) });
        let (Some(s), Some(e)) = (s, e) else {
            continue;
        };

        let analysis = analysis_col.and_then(|c| {
            c.as_any()
                .downcast_ref::<StringArray>()
                .and_then(|a| {
                    if a.is_null(i) {
                        None
                    } else {
                        Some(a.value(i).to_string())
                    }
                })
                .or_else(|| {
                    c.as_any().downcast_ref::<StringViewArray>().and_then(|a| {
                        if a.is_null(i) {
                            None
                        } else {
                            Some(a.value(i).to_string())
                        }
                    })
                })
        });
        let hseqname = hseqname_col.and_then(|c| {
            c.as_any()
                .downcast_ref::<StringArray>()
                .and_then(|a| {
                    if a.is_null(i) {
                        None
                    } else {
                        Some(a.value(i).to_string())
                    }
                })
                .or_else(|| {
                    c.as_any().downcast_ref::<StringViewArray>().and_then(|a| {
                        if a.is_null(i) {
                            None
                        } else {
                            Some(a.value(i).to_string())
                        }
                    })
                })
        });

        out.push(ProteinDomainFeature {
            analysis,
            hseqname,
            start: s,
            end: e,
        });
    }
    out
}

/// Look up overlapping protein domains for a given transcript and protein position.
///
/// Formats matching domains as `analysis:hseqname` joined with `&`, with
/// spaces, semicolons, and equals signs replaced by underscores (matching VEP).
///
/// Traceability:
/// - VEP OutputFactory.pm DOMAINS output
///   <https://github.com/Ensembl/ensembl-vep/blob/release/115/modules/Bio/EnsEMBL/VEP/OutputFactory.pm#L1448-L1466>
/// - VEP BaseTranscriptVariation.pm get_overlapping_ProteinFeatures
///   <https://github.com/Ensembl/ensembl-variation/blob/release/115/modules/Bio/EnsEMBL/Variation/BaseTranscriptVariation.pm>
/// - VEP Mapper.pm map_insert: for insertions, translation_start > translation_end
///   <https://github.com/Ensembl/ensembl/blob/release/115/modules/Bio/EnsEMBL/Mapper.pm>
fn lookup_domains(
    transcript_id: Option<&str>,
    protein_position: Option<&str>,
    amino_acids: Option<&str>,
    ctx: &PreparedContext<'_>,
) -> String {
    let Some(tx_id) = transcript_id else {
        return String::new();
    };
    let Some(pos_str) = protein_position else {
        return String::new();
    };
    if pos_str.is_empty() {
        return String::new();
    }

    // Parse protein position: single int or range "start-end".
    let (prot_start, prot_end) = if let Some((a, b)) = pos_str.split_once('-') {
        let Ok(s) = a.parse::<i64>() else {
            return String::new();
        };
        let Ok(e) = b.parse::<i64>() else {
            return String::new();
        };
        // For insertions (amino_acids = "-/X"), VEP's Mapper.map_insert swaps
        // translation_start and translation_end so that start > end. This makes
        // the overlap check exclude features that touch only at the exact
        // insertion boundary. E.g., insertion at protein 408-409 becomes
        // tl_start=409, tl_end=408: overlap with [389-408] is 409<=408 → false.
        let is_insertion = amino_acids.is_some_and(|aa| aa.starts_with("-/"));
        if is_insertion {
            (e, s) // swap: start=409, end=408
        } else {
            (s, e)
        }
    } else if pos_str.contains('?') {
        return String::new();
    } else {
        let Ok(p) = pos_str.parse::<i64>() else {
            return String::new();
        };
        (p, p)
    };

    let Some(tl) = ctx.translation_by_tx.get(tx_id) else {
        return String::new();
    };

    let mut labels: Vec<String> = Vec::new();
    for pf in &tl.protein_features {
        // Check overlap: variant protein range [prot_start, prot_end] vs feature [pf.start, pf.end]
        if prot_start <= pf.end && prot_end >= pf.start {
            let parts: Vec<&str> = [pf.analysis.as_deref(), pf.hseqname.as_deref()]
                .iter()
                .filter_map(|p| *p)
                .collect();
            if parts.is_empty() {
                continue;
            }
            let mut label = parts.join(":");
            // Replace spaces, semicolons, and equals signs with underscores.
            label = label.replace(' ', "_").replace(';', "_").replace('=', "_");
            labels.push(label);
        }
    }
    labels.join("&")
}

/// Rebuild RefSeq/XM/XR CDS strings for transcript consequence evaluation.
///
/// Traceability:
/// - Ensembl Variation `BaseTranscriptVariationAllele::_get_peptide_alleles()`
///   consumes transcript-derived CDS/peptide sequence rather than a VCF-local
///   allele string
///   <https://github.com/Ensembl/ensembl-variation/blob/release/115/modules/Bio/EnsEMBL/Variation/BaseTranscriptVariationAllele.pm#L367-L509>
///
/// For BAM-edited RefSeq transcripts Ensembl derives translateable sequence
/// from the edited transcript object (`spliced_seq` + seq edits), not directly
/// from genomic exons. Prefer that edited CDS slice when present; otherwise
/// fall back to rebuilding from genomic FASTA for transcripts that still need
/// hydration.
fn edited_refseq_translation_cds_from_spliced_seq(tx: &TranscriptFeature) -> Option<String> {
    if tx.bam_edit_status.as_deref() != Some("ok") {
        return None;
    }
    let spliced_seq = tx.spliced_seq.as_deref()?;
    let start = tx.cdna_coding_start?.checked_sub(1)?;
    let end = tx.cdna_coding_end?;
    if start >= end || end > spliced_seq.len() {
        return None;
    }
    Some(spliced_seq[start..end].to_ascii_uppercase())
}

fn hydrate_refseq_translation_cds_from_reference<R>(
    reader: &mut fasta::io::indexed_reader::IndexedReader<R>,
    transcripts: &[TranscriptFeature],
    exons: &[ExonFeature],
    translations: &mut [TranslationFeature],
    translateable_seq_by_tx: &HashMap<String, String>,
    input_variant_intervals: &HashMap<String, Vec<(i64, i64)>>,
) -> Result<HashSet<String>>
where
    R: BufRead + Seek,
{
    let mut exons_by_tx: HashMap<&str, Vec<&ExonFeature>> = HashMap::new();
    for exon in exons {
        exons_by_tx
            .entry(exon.transcript_id.as_str())
            .or_default()
            .push(exon);
    }
    for tx_exons in exons_by_tx.values_mut() {
        tx_exons.sort_by_key(|exon| (exon.start, exon.end, exon.exon_number));
    }
    let translation_ids: HashSet<&str> = translations
        .iter()
        .map(|translation| translation.transcript_id.as_str())
        .collect();

    let mut hydrated_by_tx: HashMap<String, String> = HashMap::new();
    for tx in transcripts {
        if !translation_ids.contains(tx.transcript_id.as_str()) {
            continue;
        }
        if !is_refseq_transcript_for_hydration(tx) {
            continue;
        }
        // VEP keeps transcript-local sequence state on the Transcript object.
        // If we already have `_translateable_seq`, do not overwrite it with a
        // genomic reconstruction here.
        // https://github.com/Ensembl/ensembl/blob/release/115/modules/Bio/EnsEMBL/Transcript.pm
        // https://github.com/Ensembl/ensembl-variation/blob/release/115/modules/Bio/EnsEMBL/Variation/TranscriptVariationAllele.pm#L2470-L2481
        if translateable_seq_by_tx.contains_key(&tx.transcript_id) {
            continue;
        }
        let chrom = tx.chrom.strip_prefix("chr").unwrap_or(&tx.chrom);
        let overlaps_input = input_variant_intervals
            .get(chrom)
            .map(|intervals| interval_overlaps_any(intervals, tx.start, tx.end))
            .unwrap_or(false);
        if !overlaps_input {
            continue;
        }
        if let Some(edited_cds) = edited_refseq_translation_cds_from_spliced_seq(tx) {
            hydrated_by_tx.insert(tx.transcript_id.clone(), edited_cds);
            continue;
        }
        let (Some(cds_start), Some(cds_end)) = (tx.cds_start, tx.cds_end) else {
            continue;
        };
        let Some(tx_exons) = exons_by_tx.get(tx.transcript_id.as_str()) else {
            continue;
        };
        let mut genomic_cds = String::new();
        for exon in tx_exons {
            let seg_start = exon.start.max(cds_start);
            let seg_end = exon.end.min(cds_end);
            if seg_start > seg_end {
                continue;
            }
            let segment = read_reference_sequence(reader, chrom, seg_start, seg_end)?;
            if segment.len() != usize::try_from(seg_end - seg_start + 1).unwrap_or_default() {
                genomic_cds.clear();
                break;
            }
            genomic_cds.push_str(&segment);
        }
        if genomic_cds.is_empty() {
            continue;
        }
        let cds_sequence = if tx.strand >= 0 {
            genomic_cds
        } else {
            reverse_complement_dna(&genomic_cds).ok_or_else(|| {
                DataFusionError::Execution(format!(
                    "annotate_vep(): failed to reverse-complement CDS for transcript {}",
                    tx.transcript_id
                ))
            })?
        };
        hydrated_by_tx.insert(tx.transcript_id.clone(), cds_sequence);
    }

    let mut hydrated_transcript_ids = HashSet::new();
    for translation in translations.iter_mut() {
        let Some(cds_sequence) = hydrated_by_tx.get(&translation.transcript_id) else {
            continue;
        };
        translation.cds_sequence = Some(apply_cds_phase_padding(
            translation.cds_sequence.as_deref(),
            cds_sequence.clone(),
        ));
        hydrated_transcript_ids.insert(translation.transcript_id.clone());
    }

    Ok(hydrated_transcript_ids)
}

/// Hydrate `TranscriptFeature.cdna_seq` with the full spliced transcript cDNA
/// (all exon bases concatenated, strand-oriented) so that `three_prime_utr_seq()`
/// can extract the 3' UTR for stop-loss / frameshift extension distance.
///
/// Only populates coding transcripts that already have `cdna_coding_end` but lack
/// `spliced_seq` (which is the preferred source for UTR extraction).
///
/// Traceability:
/// - VEP `_three_prime_utr()` derives UTR from `$self->transcript->seq()` which
///   is the full spliced cDNA assembled from exons
///   <https://github.com/Ensembl/ensembl-variation/blob/release/115/modules/Bio/EnsEMBL/Variation/TranscriptVariationAllele.pm#L2412-L2418>
fn hydrate_transcript_cdna_from_reference<R>(
    reader: &mut fasta::io::indexed_reader::IndexedReader<R>,
    transcripts: &mut [TranscriptFeature],
    exons: &[ExonFeature],
    indel_intervals: &HashMap<String, Vec<(i64, i64)>>,
    all_intervals: &HashMap<String, Vec<(i64, i64)>>,
) -> Result<()>
where
    R: BufRead + Seek,
{
    let mut exons_by_tx: HashMap<&str, Vec<&ExonFeature>> = HashMap::new();
    for exon in exons {
        exons_by_tx
            .entry(exon.transcript_id.as_str())
            .or_default()
            .push(exon);
    }
    for tx_exons in exons_by_tx.values_mut() {
        tx_exons.sort_by_key(|exon| (exon.start, exon.end, exon.exon_number));
    }

    for tx in transcripts.iter_mut() {
        let Some(tx_exons) = exons_by_tx.get(tx.transcript_id.as_str()) else {
            continue;
        };

        let should_infer_implicit_refseq_deletions = tx.spliced_seq.is_some()
            && is_refseq_transcript_for_hydration(tx)
            && tx.bam_edit_status.as_deref() == Some("ok")
            && tx.refseq_edits.is_empty()
            && tx.cdna_mapper_segments.is_empty();

        // Prefer cache-native spliced_seq, which already matches the live
        // transcript object Ensembl would use. For cdna_seq, only skip when
        // it clearly contains explicit 3' UTR beyond any leading phase Ns;
        // merged cache rows often store CDS-like `_translateable_seq` here.
        if tx.spliced_seq.is_some() {
            if should_infer_implicit_refseq_deletions {
                let Some(genomic_cdna) =
                    read_spliced_transcript_cdna_from_reference(reader, tx, tx_exons)?
                else {
                    continue;
                };
                let inferred = infer_refseq_deletion_edits_from_sequences(
                    &genomic_cdna,
                    tx.spliced_seq.as_deref().unwrap_or_default(),
                );
                if !inferred.is_empty() {
                    tx.refseq_edits = inferred;
                }
            }
            continue;
        }
        let Some(coding_end) = tx.cdna_coding_end else {
            continue;
        };
        if cdna_seq_has_explicit_three_prime_utr(tx) {
            continue;
        }
        // Skip LoF biotype (VEP doesn't provide UTR for these).
        if tx.biotype.contains("LoF") {
            continue;
        }
        let chrom = tx.chrom.strip_prefix("chr").unwrap_or(&tx.chrom);
        let (Some(cds_start_g), Some(cds_end_g)) = (tx.cds_start, tx.cds_end) else {
            continue;
        };
        // Hydrate if: (a) any indel overlaps the CDS (potential frameshift), OR
        // (b) any variant (incl. SNV) overlaps the stop codon (potential stop_lost).
        let indel_overlaps_cds = indel_intervals
            .get(chrom)
            .map(|iv| interval_overlaps_any(iv, cds_start_g, cds_end_g))
            .unwrap_or(false);
        let (stop_start, stop_end) = if tx.strand >= 0 {
            (cds_end_g.saturating_sub(2), cds_end_g)
        } else {
            (cds_start_g, cds_start_g.saturating_add(2))
        };
        let any_overlaps_stop = all_intervals
            .get(chrom)
            .map(|iv| interval_overlaps_any(iv, stop_start, stop_end))
            .unwrap_or(false);
        if !indel_overlaps_cds && !any_overlaps_stop {
            continue;
        }
        // Only hydrate if exons extend past CDS (i.e., there IS a 3' UTR).
        // If total exonic length <= coding_end, there's no UTR to hydrate.
        let total_exonic: usize = tx_exons
            .iter()
            .map(|e| usize::try_from(e.end.saturating_sub(e.start).saturating_add(1)).unwrap_or(0))
            .sum();
        if total_exonic <= coding_end {
            continue;
        }
        let Some(cdna) = read_spliced_transcript_cdna_from_reference(reader, tx, tx_exons)? else {
            continue;
        };
        tx.cdna_seq = Some(cdna);
    }
    Ok(())
}

fn read_spliced_transcript_cdna_from_reference<R>(
    reader: &mut fasta::io::indexed_reader::IndexedReader<R>,
    tx: &TranscriptFeature,
    tx_exons: &[&ExonFeature],
) -> Result<Option<String>>
where
    R: BufRead + Seek,
{
    let chrom = tx.chrom.strip_prefix("chr").unwrap_or(&tx.chrom);

    // Read the entire transcript span in ONE FASTA query and extract exon
    // subsequences from it. This reduces repeated per-exon FASTA reads.
    // For very large transcripts (>500KB), fall back to per-exon reads.
    let tx_span_size = tx.end - tx.start + 1;
    let genomic_cdna = if tx_span_size > 500_000 {
        let mut genomic_cdna = String::new();
        for exon in tx_exons {
            let segment = read_reference_sequence(reader, chrom, exon.start, exon.end)?;
            let expected = usize::try_from(exon.end - exon.start + 1).unwrap_or_default();
            if segment.len() != expected {
                return Ok(None);
            }
            genomic_cdna.push_str(&segment);
        }
        genomic_cdna
    } else {
        let tx_span = read_reference_sequence(reader, chrom, tx.start, tx.end)?;
        let tx_span_len = usize::try_from(tx_span_size).unwrap_or_default();
        if tx_span.len() != tx_span_len {
            return Ok(None);
        }
        let mut genomic_cdna = String::new();
        for exon in tx_exons {
            let local_start = usize::try_from(exon.start - tx.start).unwrap_or(0);
            let local_end = usize::try_from(exon.end - tx.start + 1).unwrap_or(0);
            let Some(segment) = tx_span.get(local_start..local_end) else {
                return Ok(None);
            };
            genomic_cdna.push_str(segment);
        }
        genomic_cdna
    };
    if genomic_cdna.is_empty() {
        return Ok(None);
    }

    Ok(if tx.strand >= 0 {
        Some(genomic_cdna.to_ascii_uppercase())
    } else {
        reverse_complement_dna(&genomic_cdna)
    })
}

fn cdna_seq_has_explicit_three_prime_utr(tx: &TranscriptFeature) -> bool {
    let Some(seq) = tx.cdna_seq.as_deref() else {
        return false;
    };
    let Some(coding_end) = tx.cdna_coding_end else {
        return false;
    };
    let leading_n_count = seq.bytes().take_while(|&b| b == b'N' || b == b'n').count();
    coding_end < seq.len().saturating_sub(leading_n_count)
}

fn apply_translateable_seq_overrides(
    translations: &mut [TranslationFeature],
    translateable_seq_by_tx: &HashMap<String, String>,
    hydrated_transcript_ids: &HashSet<String>,
) {
    for translation in translations {
        if hydrated_transcript_ids.contains(&translation.transcript_id) {
            continue;
        }
        if let Some(seq) = translateable_seq_by_tx.get(&translation.transcript_id) {
            translation.cds_sequence = Some(seq.clone());
        }
    }
}

fn is_refseq_transcript_for_hydration(tx: &TranscriptFeature) -> bool {
    is_refseq_transcript(tx)
}

fn read_reference_sequence<R>(
    reader: &mut fasta::io::indexed_reader::IndexedReader<R>,
    chrom: &str,
    start: i64,
    end: i64,
) -> Result<String>
where
    R: BufRead + Seek,
{
    let region = Region::new(
        chrom,
        Position::try_from(usize::try_from(start).map_err(|e| {
            DataFusionError::Execution(format!(
                "annotate_vep(): invalid FASTA start {start} for {chrom}:{start}-{end}: {e}"
            ))
        })?)
        .map_err(|e| {
            DataFusionError::Execution(format!(
                "annotate_vep(): invalid FASTA start {start} for {chrom}:{start}-{end}: {e}"
            ))
        })?..=Position::try_from(usize::try_from(end).map_err(|e| {
            DataFusionError::Execution(format!(
                "annotate_vep(): invalid FASTA end {end} for {chrom}:{start}-{end}: {e}"
            ))
        })?)
        .map_err(|e| {
            DataFusionError::Execution(format!(
                "annotate_vep(): invalid FASTA end {end} for {chrom}:{start}-{end}: {e}"
            ))
        })?,
    );
    let record = reader.query(&region).map_err(|e| {
        DataFusionError::Execution(format!(
            "annotate_vep(): failed querying FASTA for {chrom}:{start}-{end}: {e}"
        ))
    })?;
    String::from_utf8(record.sequence().as_ref().to_vec()).map_err(|e| {
        DataFusionError::Execution(format!(
            "annotate_vep(): FASTA sequence for {chrom}:{start}-{end} is not valid UTF-8: {e}"
        ))
    })
}

fn collect_input_variant_intervals(
    batches: &[RecordBatch],
) -> Result<HashMap<String, Vec<(i64, i64)>>> {
    let mut intervals_by_chrom: HashMap<String, Vec<(i64, i64)>> = HashMap::new();
    for batch in batches {
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
            intervals_by_chrom
                .entry(chrom.strip_prefix("chr").unwrap_or(&chrom).to_string())
                .or_default()
                .push((start, end));
        }
    }

    for intervals in intervals_by_chrom.values_mut() {
        intervals.sort_unstable_by_key(|interval| interval.0);
        let mut merged = Vec::with_capacity(intervals.len());
        for (start, end) in intervals.drain(..) {
            if let Some((_, last_end)) = merged.last_mut() {
                if start <= *last_end {
                    *last_end = (*last_end).max(end);
                    continue;
                }
            }
            merged.push((start, end));
        }
        *intervals = merged;
    }

    Ok(intervals_by_chrom)
}

/// Collect intervals for indel-only variants (ref_len != alt_len).
/// Used to restrict cDNA hydration to transcripts that might have
/// frameshift or stop-loss consequences requiring UTR extension.
fn collect_indel_variant_intervals(
    batches: &[RecordBatch],
) -> Result<HashMap<String, Vec<(i64, i64)>>> {
    let mut intervals_by_chrom: HashMap<String, Vec<(i64, i64)>> = HashMap::new();
    for batch in batches {
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
        let ref_idx = schema.index_of("reference").ok();
        let alt_idx = schema.index_of("alternate").ok();
        for row in 0..batch.num_rows() {
            // Only include indels (different ref/alt lengths).
            let is_indel = match (ref_idx, alt_idx) {
                (Some(ri), Some(ai)) => {
                    let ref_len = string_at(batch.column(ri).as_ref(), row)
                        .map(|s| s.len())
                        .unwrap_or(0);
                    let alt_len = string_at(batch.column(ai).as_ref(), row)
                        .map(|s| s.len())
                        .unwrap_or(0);
                    ref_len != alt_len
                }
                _ => true, // If we can't determine, include conservatively.
            };
            if !is_indel {
                continue;
            }
            let Some(chrom) = string_at(batch.column(chrom_idx).as_ref(), row) else {
                continue;
            };
            let Some(start) = int64_at(batch.column(start_idx).as_ref(), row) else {
                continue;
            };
            let Some(end) = int64_at(batch.column(end_idx).as_ref(), row) else {
                continue;
            };
            intervals_by_chrom
                .entry(chrom.strip_prefix("chr").unwrap_or(&chrom).to_string())
                .or_default()
                .push((start, end));
        }
    }

    for intervals in intervals_by_chrom.values_mut() {
        intervals.sort_unstable_by_key(|interval| interval.0);
        let mut merged = Vec::with_capacity(intervals.len());
        for (start, end) in intervals.drain(..) {
            if let Some((_, last_end)) = merged.last_mut() {
                if start <= *last_end {
                    *last_end = (*last_end).max(end);
                    continue;
                }
            }
            merged.push((start, end));
        }
        *intervals = merged;
    }

    Ok(intervals_by_chrom)
}

fn interval_overlaps_any(intervals: &[(i64, i64)], start: i64, end: i64) -> bool {
    let idx = intervals.partition_point(|interval| interval.1 < start);
    idx < intervals.len() && intervals[idx].0 <= end
}

/// Parse cached TranscriptMapper exon-to-cDNA pairs from serialized transcript
/// `raw_object_json`.
///
/// Traceability:
/// - Ensembl VEP `AnnotationSource::Database::Transcript::prefetch_transcript_data()`
///   caches `mapper` on `_variation_effect_feature_cache`
///   <https://github.com/Ensembl/ensembl-vep/blob/release/115/modules/Bio/EnsEMBL/VEP/AnnotationSource/Database/Transcript.pm#L333-L352>
/// - Ensembl Variation `TranscriptVariationAllele::_get_cDNA_position()`
///   resolves transcript positions through TranscriptMapper `genomic2cdna`
///   <https://github.com/Ensembl/ensembl-variation/blob/release/115/modules/Bio/EnsEMBL/Variation/TranscriptVariationAllele.pm#L2683-L2765>
/// Read cdna_mapper_segments from a promoted List<Struct> parquet column.
fn cdna_mapper_segments_from_list_column(
    col: &dyn Array,
    row: usize,
) -> Vec<TranscriptCdnaMapperSegment> {
    use datafusion::arrow::array::{AsArray, Int8Array, Int64Array, ListArray, StructArray};

    let list_array = col.as_any().downcast_ref::<ListArray>();
    let Some(list_array) = list_array else {
        return Vec::new();
    };
    if list_array.is_null(row) {
        return Vec::new();
    }
    let start = list_array.value_offsets()[row] as usize;
    let end = list_array.value_offsets()[row + 1] as usize;
    if start == end {
        return Vec::new();
    }
    let values = list_array.values();
    let struct_array = values.as_any().downcast_ref::<StructArray>();
    let Some(struct_array) = struct_array else {
        return Vec::new();
    };
    let gs = struct_array
        .column_by_name("genomic_start")
        .and_then(|c| c.as_any().downcast_ref::<Int64Array>());
    let ge = struct_array
        .column_by_name("genomic_end")
        .and_then(|c| c.as_any().downcast_ref::<Int64Array>());
    let cs = struct_array
        .column_by_name("cdna_start")
        .and_then(|c| c.as_any().downcast_ref::<Int64Array>());
    let ce = struct_array
        .column_by_name("cdna_end")
        .and_then(|c| c.as_any().downcast_ref::<Int64Array>());
    let ori = struct_array
        .column_by_name("ori")
        .and_then(|c| c.as_any().downcast_ref::<Int8Array>());
    let (Some(gs), Some(ge), Some(cs), Some(ce), Some(ori)) = (gs, ge, cs, ce, ori) else {
        return Vec::new();
    };
    let mut segments = Vec::with_capacity(end - start);
    for i in start..end {
        segments.push(TranscriptCdnaMapperSegment {
            genomic_start: gs.value(i),
            genomic_end: ge.value(i),
            cdna_start: cs.value(i) as usize,
            cdna_end: ce.value(i) as usize,
            ori: ori.value(i),
        });
    }
    segments
}

fn apply_cds_phase_padding(existing_cds: Option<&str>, mut hydrated_cds: String) -> String {
    let leading_phase_padding = existing_cds
        .map(|seq| seq.bytes().take_while(|&b| b == b'N' || b == b'n').count())
        .unwrap_or(0);
    if leading_phase_padding == 0 {
        return hydrated_cds;
    }
    let mut padded = String::with_capacity(leading_phase_padding + hydrated_cds.len());
    padded.extend(std::iter::repeat_n('N', leading_phase_padding));
    padded.push_str(&hydrated_cds);
    hydrated_cds.clear();
    padded
}

fn reverse_complement_dna(seq: &str) -> Option<String> {
    let mut out = String::with_capacity(seq.len());
    for base in seq.as_bytes().iter().rev() {
        let comp = match base.to_ascii_uppercase() {
            b'A' => 'T',
            b'C' => 'G',
            b'G' => 'C',
            b'T' => 'A',
            b'N' => 'N',
            _ => return None,
        };
        out.push(comp);
    }
    Some(out)
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

pub(crate) fn int64_at(array: &dyn Array, row: usize) -> Option<i64> {
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

pub(crate) fn string_at(array: &dyn Array, row: usize) -> Option<String> {
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

// ---------------------------------------------------------------------------
// Streaming contig-by-contig ExecutionPlan
// ---------------------------------------------------------------------------

/// Bundled configuration for per-contig annotation (extracted from AnnotateProvider).
#[derive(Clone)]
struct ContigAnnotationConfig {
    vcf_table: String,
    options_json: Option<String>,
    cache_columns: Vec<String>,
    extended_probes: bool,
    translations_sift_table: Option<String>,
    flags: VepFlags,
    hgvs_flags: HgvsFlags,
    transcript_selection: TranscriptSelectionFlags,
    allowed_failed: i64,
    reference_fasta_path: Option<String>,
    upstream_distance: i64,
    downstream_distance: i64,
    input_buffer_size: usize,
    projection: Option<Vec<usize>>,
    annotation_column_count: usize,
    /// Maximum number of output rows (LIMIT pushdown).
    fetch_limit: Option<usize>,
    pick_flags: PickFlags,
    /// When true, use fjall KV store for variation lookup + SIFT instead of parquet.
    #[cfg(feature = "kv-cache")]
    use_fjall: bool,
    /// Shared fjall KV store handle (opened once, reused across contigs).
    #[cfg(feature = "kv-cache")]
    kv_store: Option<Arc<crate::kv_cache::VepKvStore>>,
    /// Shared fjall SIFT store (opened once, reused across contigs).
    #[cfg(feature = "kv-cache")]
    sift_kv_store: Option<crate::kv_cache::SiftKvStore>,
}

/// Leaf `ExecutionPlan` that processes contigs one at a time via a state-machine
/// stream, reclaiming memory after each contig completes.
struct ContigAnnotationExec {
    projected_schema: SchemaRef,
    full_schema: SchemaRef,
    contigs: Vec<String>,
    session: Arc<SessionContext>,
    cache: Arc<PartitionedParquetCache>,
    config: ContigAnnotationConfig,
    properties: PlanProperties,
}

impl ContigAnnotationExec {
    fn new(
        projected_schema: SchemaRef,
        full_schema: SchemaRef,
        contigs: Vec<String>,
        session: Arc<SessionContext>,
        cache: Arc<PartitionedParquetCache>,
        config: ContigAnnotationConfig,
    ) -> Self {
        let properties = PlanProperties::new(
            EquivalenceProperties::new(projected_schema.clone()),
            datafusion::physical_plan::Partitioning::UnknownPartitioning(1),
            EmissionType::Incremental,
            Boundedness::Bounded,
        );
        Self {
            projected_schema,
            full_schema,
            contigs,
            session,
            cache,
            config,
            properties,
        }
    }
}

impl Debug for ContigAnnotationExec {
    fn fmt(&self, f: &mut Formatter<'_>) -> std::fmt::Result {
        write!(
            f,
            "ContigAnnotationExec {{ contigs: {}, schema_fields: {} }}",
            self.contigs.len(),
            self.projected_schema.fields().len()
        )
    }
}

impl DisplayAs for ContigAnnotationExec {
    fn fmt_as(&self, _t: DisplayFormatType, f: &mut Formatter<'_>) -> std::fmt::Result {
        write!(f, "ContigAnnotationExec: contigs={}", self.contigs.len())
    }
}

impl ExecutionPlan for ContigAnnotationExec {
    fn name(&self) -> &str {
        "ContigAnnotationExec"
    }

    fn as_any(&self) -> &dyn Any {
        self
    }

    fn schema(&self) -> SchemaRef {
        self.projected_schema.clone()
    }

    fn properties(&self) -> &PlanProperties {
        &self.properties
    }

    fn children(&self) -> Vec<&Arc<dyn ExecutionPlan>> {
        vec![]
    }

    fn with_new_children(
        self: Arc<Self>,
        _children: Vec<Arc<dyn ExecutionPlan>>,
    ) -> Result<Arc<dyn ExecutionPlan>> {
        Ok(self)
    }

    fn execute(
        &self,
        _partition: usize,
        _context: Arc<TaskContext>,
    ) -> Result<SendableRecordBatchStream> {
        Ok(Box::pin(ContigAnnotationStream::new(
            self.projected_schema.clone(),
            self.full_schema.clone(),
            self.contigs.clone(),
            Arc::clone(&self.session),
            Arc::clone(&self.cache),
            self.config.clone(),
        )))
    }
}

// ---------------------------------------------------------------------------
// State-machine stream — window-based streaming annotation
// ---------------------------------------------------------------------------

/// Number of lookup batches to accumulate before hydrating + annotating.
/// Each window triggers a PreparedContext rebuild (~22ms).
/// With ~30 rows/batch: 1000 batches ≈ 30K variants per window.
const HYDRATION_WINDOW_SIZE: usize = 1000;
/// Ensembl VEP release/115 default `buffer_size`.
const VEP_INPUT_BUFFER_SIZE: usize = 5000;
/// Ensembl VEP transcript cache region size (`cache_region_size`).
const VEP_TRANSCRIPT_CACHE_REGION_SIZE_BP: i64 = 1_000_000;

#[derive(Debug, Clone, PartialEq, Eq, Hash)]
struct TranscriptCacheRegion {
    chrom: String,
    region_index: i64,
}

type PersistedBufferTranscripts =
    HashMap<TranscriptCacheRegion, HashMap<String, TranscriptFeature>>;

type FastaReader = fasta::io::indexed_reader::IndexedReader<fasta::io::BufReader<std::fs::File>>;

/// Everything needed to start streaming annotation for a contig.
/// Produced by `prepare_contig_context()`.
struct ContigReadyState {
    lookup_stream: SendableRecordBatchStream,
    coloc_sink: ColocatedSink,
    transcripts: Vec<TranscriptFeature>,
    translateable_seq_by_tx: HashMap<String, String>,
    exons: Vec<ExonFeature>,
    translations: Vec<TranslationFeature>,
    regulatory: Vec<RegulatoryFeature>,
    motifs: Vec<MotifFeature>,
    mirnas: Vec<MirnaFeature>,
    structural: Vec<StructuralFeature>,
    ephemeral_tables: Vec<String>,
    chrom: String,
    t_contig: Instant,
}

/// Mutable annotation state for window-based streaming within a single contig.
struct ContigAnnotationState {
    /// Dropped after exhaustion to reclaim BuildSide memory (COITrees, hash
    /// indices, concatenated VCF batch).
    lookup_stream: Option<SendableRecordBatchStream>,
    // Context data (transcripts/translations mutated by HGVS hydration).
    transcripts: Vec<TranscriptFeature>,
    translateable_seq_by_tx: HashMap<String, String>,
    exons: Vec<ExonFeature>,
    translations: Vec<TranslationFeature>,
    regulatory: Vec<RegulatoryFeature>,
    motifs: Vec<MotifFeature>,
    mirnas: Vec<MirnaFeature>,
    structural: Vec<StructuralFeature>,
    transcript_cache_regions: HashMap<String, Vec<TranscriptCacheRegion>>,
    persisted_buffer_transcripts: PersistedBufferTranscripts,
    // Colocated (built lazily from sink on first window).
    colocated_map: HashMap<ColocatedKey, ColocatedData>,
    colocated_map_built: bool,
    coloc_sink: ColocatedSink,
    // HGVS hydration tracking — already-hydrated transcripts are skipped.
    hydrated_cds_tx_ids: HashSet<String>,
    hgvs_reader: Option<FastaReader>,
    // SIFT state (same sliding-window pattern as before).
    sift_cache: SiftPolyphenCache,
    sift_direct: Option<SiftDirectReader>,
    loaded_sift_windows: HashSet<(String, i64)>,
    /// Fjall-backed SIFT store for lazy per-transcript lookups (used when use_fjall=true).
    #[cfg(feature = "kv-cache")]
    sift_kv: Option<crate::kv_cache::SiftKvStore>,
    // Annotation engine + provider.
    tmp_provider: AnnotateProvider,
    engine: TranscriptConsequenceEngine,
    // Window buffer.
    window_buffer: Vec<RecordBatch>,
    lookup_done: bool,
    // Cleanup + profiling.
    ephemeral_tables: Vec<String>,
    chrom: String,
    config: ContigAnnotationConfig,
    session: Arc<SessionContext>,
    cache: Arc<PartitionedParquetCache>,
    t_contig: Instant,
    contig_rows: usize,
}

type PrepareFuture = Pin<Box<dyn Future<Output = Result<Option<ContigReadyState>>> + Send>>;
type CleanupFuture = Pin<Box<dyn Future<Output = Result<()>> + Send>>;

enum StreamState {
    StartContig,
    /// Async setup: register tables, parallel context load + lookup stream creation.
    PreparingContig(PrepareFuture),
    /// Pull from lookup stream, accumulate windows, hydrate + annotate.
    AnnotatingContig(ContigAnnotationState),
    /// Yield annotated batches from the current window, then resume annotation.
    DrainingWindow {
        batches: VecDeque<RecordBatch>,
        annotation_state: ContigAnnotationState,
    },
    /// Deregister ephemeral tables after contig completes.
    CleaningUp(CleanupFuture),
    /// Deregister ephemeral tables after an error, then propagate the error.
    ErrorCleaningUp(CleanupFuture, DataFusionError),
    /// Final async cleanup (e.g., deregister global KV table) before Done.
    /// Unlike CleaningUp, transitions to Done instead of StartContig.
    FinalCleanup(CleanupFuture),
    Done,
}

/// Spawn an async future to deregister ephemeral tables.
fn make_cleanup_future(session: Arc<SessionContext>, tables: Vec<String>) -> CleanupFuture {
    Box::pin(async move {
        for tbl in &tables {
            crate::partitioned_cache::deregister_table(&session, tbl).await?;
        }
        Ok(())
    })
}

fn deregister_tables_sync(session: &SessionContext, tables: &[String]) {
    for table in tables {
        let _ = session.deregister_table(table);
    }
}

struct ContigAnnotationStream {
    projected_schema: SchemaRef,
    full_schema: SchemaRef,
    contigs: VecDeque<String>,
    session: Arc<SessionContext>,
    cache: Arc<PartitionedParquetCache>,
    config: ContigAnnotationConfig,
    state: StreamState,
    /// Rows emitted so far (for LIMIT pushdown).
    rows_emitted: usize,
}

impl ContigAnnotationStream {
    fn new(
        projected_schema: SchemaRef,
        full_schema: SchemaRef,
        contigs: Vec<String>,
        session: Arc<SessionContext>,
        cache: Arc<PartitionedParquetCache>,
        config: ContigAnnotationConfig,
    ) -> Self {
        Self {
            projected_schema,
            full_schema,
            contigs: VecDeque::from(contigs),
            session,
            cache,
            config,
            state: StreamState::StartContig,
            rows_emitted: 0,
        }
    }

    fn cleanup_registered_tables_on_drop(&mut self) {
        match &mut self.state {
            StreamState::AnnotatingContig(ann) => {
                deregister_tables_sync(&ann.session, &ann.ephemeral_tables);
                if ann.config.use_fjall {
                    let _ = ann.session.deregister_table("__vep_kv_variation");
                }
                ann.ephemeral_tables.clear();
            }
            StreamState::DrainingWindow {
                annotation_state: ann,
                ..
            } => {
                deregister_tables_sync(&ann.session, &ann.ephemeral_tables);
                if ann.config.use_fjall {
                    let _ = ann.session.deregister_table("__vep_kv_variation");
                }
                ann.ephemeral_tables.clear();
            }
            StreamState::CleaningUp(_) | StreamState::ErrorCleaningUp(_, _) => {
                if self.config.use_fjall {
                    let _ = self.session.deregister_table("__vep_kv_variation");
                }
            }
            StreamState::StartContig
            | StreamState::PreparingContig(_)
            | StreamState::FinalCleanup(_)
            | StreamState::Done => {
                if self.config.use_fjall {
                    let _ = self.session.deregister_table("__vep_kv_variation");
                }
            }
        }
    }
}

impl Drop for ContigAnnotationStream {
    fn drop(&mut self) {
        self.cleanup_registered_tables_on_drop();
    }
}

impl RecordBatchStream for ContigAnnotationStream {
    fn schema(&self) -> SchemaRef {
        self.projected_schema.clone()
    }
}

/// Hydrate a window of batches: compute variant intervals, hydrate new
/// transcripts (skip already-hydrated), apply translateable_seq overrides.
/// Mirrors the SIFT sliding-window pattern.
///
/// We intentionally do not gate this on `shift_hgvs`: start-codon indel
/// classification needs hydrated cDNA/spliced sequence even when callers only
/// request consequence terms and HGVSc shifting is disabled.
fn hydrate_window(
    transcripts: &mut [TranscriptFeature],
    exons: &[ExonFeature],
    translations: &mut [TranslationFeature],
    translateable_seq_by_tx: &HashMap<String, String>,
    hgvs_reader: &mut Option<FastaReader>,
    hydrated_cds_tx_ids: &mut HashSet<String>,
    window_batches: &[RecordBatch],
) -> Result<()> {
    let Some(reader) = hgvs_reader.as_mut() else {
        return Ok(());
    };

    let input_intervals = collect_input_variant_intervals(window_batches)?;
    let indel_intervals = collect_indel_variant_intervals(window_batches)?;

    // CDS hydration — skip transcripts already hydrated in previous windows.
    let newly_hydrated = hydrate_refseq_translation_cds_from_reference(
        reader,
        transcripts,
        exons,
        translations,
        translateable_seq_by_tx,
        &input_intervals,
    )?;
    // Merge newly hydrated IDs into cumulative set.
    let first_window = hydrated_cds_tx_ids.is_empty();
    hydrated_cds_tx_ids.extend(newly_hydrated.iter().cloned());

    // Apply translateable_seq overrides for newly hydrated transcripts.
    apply_translateable_seq_overrides(translations, translateable_seq_by_tx, &newly_hydrated);

    // cDNA hydration (already skips transcripts with existing spliced_seq).
    hydrate_transcript_cdna_from_reference(
        reader,
        transcripts,
        exons,
        &indel_intervals,
        &input_intervals,
    )?;

    if first_window && profiling_enabled() {
        eprintln!(
            "[VEP_PROFILE] first window HGVS hydration: {} CDS transcripts",
            hydrated_cds_tx_ids.len()
        );
    }
    Ok(())
}

fn split_batches_by_row_limit(batches: &[RecordBatch], row_limit: usize) -> Vec<Vec<RecordBatch>> {
    let mut out = Vec::new();
    let mut current = Vec::new();
    let mut current_rows = 0usize;

    for batch in batches {
        let mut offset = 0usize;
        while offset < batch.num_rows() {
            let remaining = row_limit.saturating_sub(current_rows);
            let take = remaining.min(batch.num_rows() - offset);
            current.push(batch.slice(offset, take));
            current_rows += take;
            offset += take;

            if current_rows == row_limit {
                out.push(std::mem::take(&mut current));
                current_rows = 0;
            }
        }
    }

    if !current.is_empty() {
        out.push(current);
    }

    out
}

fn buffer_variant_bounds(batches: &[RecordBatch]) -> Result<Option<(String, i64, i64)>> {
    let mut chrom: Option<String> = None;
    let mut min_start = i64::MAX;
    let mut max_end = i64::MIN;

    for batch in batches {
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

        for row in 0..batch.num_rows() {
            let Some(row_chrom) = string_at(batch.column(chrom_idx).as_ref(), row) else {
                continue;
            };
            let Some(start) = int64_at(batch.column(start_idx).as_ref(), row) else {
                continue;
            };
            let Some(end) = int64_at(batch.column(end_idx).as_ref(), row) else {
                continue;
            };
            if chrom.is_none() {
                chrom = Some(row_chrom);
            }
            min_start = min_start.min(start.min(end));
            max_end = max_end.max(start.max(end));
        }
    }

    Ok(chrom.map(|chrom| (chrom, min_start, max_end)))
}

fn normalized_chrom_name(chrom: &str) -> &str {
    chrom.strip_prefix("chr").unwrap_or(chrom)
}

fn cache_region_index(pos: i64) -> i64 {
    pos.saturating_sub(1) / VEP_TRANSCRIPT_CACHE_REGION_SIZE_BP
}

fn cache_regions_for_coords(
    chrom: &str,
    start: i64,
    end: i64,
    upstream_distance: i64,
    downstream_distance: i64,
) -> Vec<TranscriptCacheRegion> {
    let chrom = normalized_chrom_name(chrom).to_string();
    // Traceability:
    // - Ensembl VEP `AnnotationType::Transcript::up_down_size()` returns
    //   max(UPSTREAM_DISTANCE, DOWNSTREAM_DISTANCE) — a single scalar used
    //   symmetrically to bound the coarse cache-region fetch.
    //   <https://github.com/Ensembl/ensembl-vep/blob/release/115/modules/Bio/EnsEMBL/VEP/AnnotationType/Transcript.pm#L170-L177>
    // - Ensembl VEP `AnnotationSource::get_regions_from_coords()` applies it as
    //   `[start - up_down_size, end + up_down_size]`.
    //   <https://github.com/Ensembl/ensembl-vep/blob/release/115/modules/Bio/EnsEMBL/VEP/AnnotationSource.pm#L190>
    // Strand-aware asymmetric upstream/downstream gating happens later in
    // `upstream_downstream_term()`.
    let up_down_size = upstream_distance.max(downstream_distance);
    let query_start = start.min(end).saturating_sub(up_down_size);
    let query_end = start.max(end).saturating_add(up_down_size);
    let region_start = cache_region_index(query_start);
    let region_end = cache_region_index(query_end);

    (region_start..=region_end)
        .map(|region_index| TranscriptCacheRegion {
            chrom: chrom.clone(),
            region_index,
        })
        .collect()
}

fn transcript_cache_regions(transcript: &TranscriptFeature) -> Vec<TranscriptCacheRegion> {
    cache_regions_for_coords(&transcript.chrom, transcript.start, transcript.end, 0, 0)
}

fn collect_buffer_cache_regions(
    batches: &[RecordBatch],
    upstream_distance: i64,
    downstream_distance: i64,
) -> Result<HashSet<TranscriptCacheRegion>> {
    let mut regions = HashSet::new();

    for batch in batches {
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

            regions.extend(cache_regions_for_coords(
                &chrom,
                start,
                end,
                upstream_distance,
                downstream_distance,
            ));
        }
    }

    Ok(regions)
}

/// Retain only transcript objects that still belong to active Ensembl
/// transcript cache regions, mirroring `clean_cache()` pruning of cached
/// transcript objects between adjacent input buffers.
///
/// Traceability:
/// - Ensembl VEP `AnnotationSource::clean_cache()`
///   <https://github.com/Ensembl/ensembl-vep/blob/release/115/modules/Bio/EnsEMBL/VEP/AnnotationSource.pm#L392-L422>
/// - Ensembl VEP region cache population in `get_features_by_regions_uncached()`
///   <https://github.com/Ensembl/ensembl-vep/blob/release/115/modules/Bio/EnsEMBL/VEP/AnnotationSource.pm#L109-L139>
fn prune_persisted_buffer_transcripts(
    persisted_transcripts: &mut PersistedBufferTranscripts,
    active_regions: &HashSet<TranscriptCacheRegion>,
) {
    persisted_transcripts.retain(|region, _| active_regions.contains(region));
}

fn select_buffer_local_transcripts(
    transcripts: &[TranscriptFeature],
    chrom: &str,
    min_start: i64,
    max_end: i64,
    upstream_distance: i64,
    downstream_distance: i64,
) -> Vec<TranscriptFeature> {
    let chrom_norm = chrom.strip_prefix("chr").unwrap_or(chrom);
    // VEP's up_down_size = max(upstream, downstream), applied symmetrically.
    // See collect_buffer_cache_regions above for the full traceability.
    let up_down_size = upstream_distance.max(downstream_distance);
    let query_start = min_start.saturating_sub(up_down_size);
    let query_end = max_end.saturating_add(up_down_size);

    transcripts
        .iter()
        .filter(|tx| {
            let tx_chrom = tx.chrom.strip_prefix("chr").unwrap_or(&tx.chrom);
            tx_chrom == chrom_norm && tx.end >= query_start && tx.start <= query_end
        })
        .cloned()
        .collect()
}

fn build_buffer_local_transcripts(
    transcripts: &[TranscriptFeature],
    chrom: &str,
    min_start: i64,
    max_end: i64,
    upstream_distance: i64,
    downstream_distance: i64,
) -> Vec<TranscriptFeature> {
    let mut buffer_transcripts = select_buffer_local_transcripts(
        transcripts,
        chrom,
        min_start,
        max_end,
        upstream_distance,
        downstream_distance,
    );
    reset_buffer_local_hgnc_effective_values(&mut buffer_transcripts);
    apply_buffer_local_hgnc_propagation(&mut buffer_transcripts);
    buffer_transcripts
}

fn build_stateful_buffer_local_transcripts(
    transcripts: &[TranscriptFeature],
    transcript_cache_regions: &HashMap<String, Vec<TranscriptCacheRegion>>,
    persisted_transcripts: &mut PersistedBufferTranscripts,
    buffer_batches: &[RecordBatch],
    chrom: &str,
    min_start: i64,
    max_end: i64,
    upstream_distance: i64,
    downstream_distance: i64,
) -> Result<Vec<TranscriptFeature>> {
    let active_regions =
        collect_buffer_cache_regions(buffer_batches, upstream_distance, downstream_distance)?;
    prune_persisted_buffer_transcripts(persisted_transcripts, &active_regions);

    let mut buffer_transcripts = select_buffer_local_transcripts(
        transcripts,
        chrom,
        min_start,
        max_end,
        upstream_distance,
        downstream_distance,
    );

    // Reuse the same transcript objects across adjacent input buffers while
    // their 1 Mb cache regions remain active, matching Ensembl VEP's region
    // cache plus in-place `merge_features()` mutation behavior.
    //
    // Traceability:
    // - Ensembl VEP `AnnotationSource::get_features_by_regions_uncached()`
    //   caches feature objects by region
    //   <https://github.com/Ensembl/ensembl-vep/blob/release/115/modules/Bio/EnsEMBL/VEP/AnnotationSource.pm#L109-L139>
    // - Ensembl VEP `AnnotationType::Transcript::merge_features()`
    //   mutates transcript objects in place
    //   <https://github.com/Ensembl/ensembl-vep/blob/release/115/modules/Bio/EnsEMBL/VEP/AnnotationType/Transcript.pm#L246-L310>
    for tx in &mut buffer_transcripts {
        let persisted = active_regions.iter().find_map(|region| {
            persisted_transcripts
                .get(region)
                .and_then(|by_transcript| by_transcript.get(&tx.transcript_id))
        });
        if let Some(persisted) = persisted {
            *tx = persisted.clone();
        }
    }

    reset_persisted_hgnc_effective_values_outside_start_region(
        &mut buffer_transcripts,
        transcript_cache_regions,
        &active_regions,
    );
    apply_buffer_local_hgnc_propagation(&mut buffer_transcripts);

    for tx in &buffer_transcripts {
        if let Some(regions) = transcript_cache_regions.get(&tx.transcript_id) {
            for region in regions
                .iter()
                .filter(|region| active_regions.contains(*region))
            {
                persisted_transcripts
                    .entry(region.clone())
                    .or_default()
                    .insert(tx.transcript_id.clone(), tx.clone());
            }
        }
    }

    Ok(buffer_transcripts)
}

fn reset_buffer_local_hgnc_effective_values(transcripts: &mut [TranscriptFeature]) {
    for tx in transcripts {
        tx.gene_hgnc_id = tx.gene_hgnc_id_native.clone();
    }
}

fn reset_persisted_hgnc_effective_values_outside_start_region(
    transcripts: &mut [TranscriptFeature],
    transcript_cache_regions: &HashMap<String, Vec<TranscriptCacheRegion>>,
    active_regions: &HashSet<TranscriptCacheRegion>,
) {
    for tx in transcripts {
        let Some(regions) = transcript_cache_regions.get(&tx.transcript_id) else {
            continue;
        };
        let Some(start_region) = regions.first() else {
            continue;
        };
        let spans_multiple_regions = regions.len() > 1;
        if spans_multiple_regions && !active_regions.contains(start_region) {
            tx.gene_hgnc_id = tx.gene_hgnc_id_native.clone();
        }
    }
}

/// Apply the HGNC-relevant parts of Ensembl VEP `merge_features()` to a single
/// input buffer's transcript feature set.
///
/// Traceability:
/// - Ensembl VEP `AnnotationType::Transcript::merge_features()`
///   <https://github.com/Ensembl/ensembl-vep/blob/release/115/modules/Bio/EnsEMBL/VEP/AnnotationType/Transcript.pm#L246-L310>
fn apply_buffer_local_hgnc_propagation(transcripts: &mut [TranscriptFeature]) {
    #[derive(Default)]
    struct GeneFill {
        gene_symbol: Option<String>,
        gene_symbol_source: Option<String>,
        gene_hgnc_id_native: Option<String>,
    }

    let mut hgnc_by_symbol: HashMap<String, String> = HashMap::new();
    let mut gene_fill_by_stable_id: HashMap<String, GeneFill> = HashMap::new();

    for tx in transcripts.iter() {
        if let (Some(symbol), Some(hgnc_id)) =
            (tx.gene_symbol.as_deref(), tx.gene_hgnc_id_native.as_deref())
        {
            hgnc_by_symbol
                .entry(symbol.to_string())
                .or_insert_with(|| hgnc_id.to_string());
        }

        if let Some(gene_stable_id) = tx.gene_stable_id.as_deref() {
            let fill = gene_fill_by_stable_id
                .entry(gene_stable_id.to_string())
                .or_default();
            if fill.gene_symbol.is_none() {
                fill.gene_symbol = tx.gene_symbol.clone();
            }
            if fill.gene_symbol_source.is_none() {
                fill.gene_symbol_source = tx.gene_symbol_source.clone();
            }
            if fill.gene_hgnc_id_native.is_none() {
                fill.gene_hgnc_id_native = tx.gene_hgnc_id_native.clone();
            }
        }
    }

    for tx in transcripts.iter_mut() {
        tx.gene_hgnc_id = tx
            .gene_hgnc_id_native
            .clone()
            .or_else(|| tx.gene_hgnc_id.clone());

        if tx.gene_hgnc_id.is_none() {
            if let Some(hgnc_id) = tx
                .gene_symbol
                .as_deref()
                .and_then(|symbol| hgnc_by_symbol.get(symbol))
            {
                tx.gene_hgnc_id = Some(hgnc_id.clone());
            }
        }
    }

    for tx in transcripts.iter_mut() {
        let Some(gene_stable_id) = tx.gene_stable_id.as_deref() else {
            continue;
        };
        let Some(fill) = gene_fill_by_stable_id.get(gene_stable_id) else {
            continue;
        };

        if tx.gene_symbol.is_none() {
            tx.gene_symbol = fill.gene_symbol.clone();
        }
        if tx.gene_symbol_source.is_none() {
            tx.gene_symbol_source = fill.gene_symbol_source.clone();
        }
        if tx.gene_hgnc_id.is_none() {
            tx.gene_hgnc_id = fill.gene_hgnc_id_native.clone();
        }
    }
}

/// Annotate a window of batches: build PreparedContext, run SIFT loading,
/// annotate each batch, apply projection.
fn annotate_window(
    ann: &mut ContigAnnotationState,
    window_batches: &[RecordBatch],
    projection: Option<&[usize]>,
) -> Result<VecDeque<RecordBatch>> {
    let engine = &ann.engine;
    let csq_col_idx = ann.tmp_provider.vcf_field_count();
    let skip_csq = projection.is_some_and(|indices| !indices.contains(&csq_col_idx));
    let typed_cols_start = csq_col_idx + 2;
    let typed_cols_end = typed_cols_start + ann.tmp_provider.annotation_column_defs.len();
    let skip_typed_cols = projection.map_or(false, |indices| {
        !indices
            .iter()
            .any(|&i| i >= typed_cols_start && i < typed_cols_end)
    });
    let pick_requires_full_annotations = ann
        .config
        .pick_flags
        .requires_transcript_annotations(skip_csq, skip_typed_cols);
    let sift_enabled = ann.config.flags.everything;
    let mut out = VecDeque::with_capacity(window_batches.len());

    for buffer_batches in split_batches_by_row_limit(window_batches, ann.config.input_buffer_size) {
        let Some((chrom, min_start, max_end)) = buffer_variant_bounds(&buffer_batches)? else {
            continue;
        };
        let buffer_transcripts = build_stateful_buffer_local_transcripts(
            &ann.transcripts,
            &ann.transcript_cache_regions,
            &mut ann.persisted_buffer_transcripts,
            &buffer_batches,
            &chrom,
            min_start,
            max_end,
            ann.config.upstream_distance,
            ann.config.downstream_distance,
        )?;
        let buffer_tx_ids: HashSet<&str> = buffer_transcripts
            .iter()
            .map(|tx| tx.transcript_id.as_str())
            .collect();
        let buffer_exons: Vec<ExonFeature> = ann
            .exons
            .iter()
            .filter(|exon| buffer_tx_ids.contains(exon.transcript_id.as_str()))
            .cloned()
            .collect();
        let buffer_translations: Vec<TranslationFeature> = ann
            .translations
            .iter()
            .filter(|translation| buffer_tx_ids.contains(translation.transcript_id.as_str()))
            .cloned()
            .collect();
        let ctx = PreparedContext::new(
            &buffer_transcripts,
            &buffer_exons,
            &buffer_translations,
            &ann.regulatory,
            &ann.motifs,
            &ann.mirnas,
            &ann.structural,
        );

        for batch in &buffer_batches {
            // Lazy SIFT window loading (same pattern as before).
            if sift_enabled && ann.sift_direct.is_some() {
                let batch_needs_engine = if pick_requires_full_annotations {
                    true
                } else {
                    batch
                        .schema()
                        .index_of("cache_most_severe_consequence")
                        .ok()
                        .map_or(true, |idx| batch.column(idx).null_count() > 0)
                };
                if batch_needs_engine {
                    let schema = batch.schema();
                    if let (Ok(ci), Ok(si), Ok(ei)) = (
                        schema.index_of("chrom"),
                        schema.index_of("start"),
                        schema.index_of("end"),
                    ) {
                        let mut batch_chrom_bounds: HashMap<String, (i64, i64)> = HashMap::new();
                        for row in 0..batch.num_rows() {
                            if let (Some(c), Some(s), Some(e)) = (
                                string_at(batch.column(ci).as_ref(), row),
                                int64_at(batch.column(si).as_ref(), row),
                                int64_at(batch.column(ei).as_ref(), row),
                            ) {
                                let c_norm = c.strip_prefix("chr").unwrap_or(&c).to_string();
                                let entry = batch_chrom_bounds
                                    .entry(c_norm)
                                    .or_insert((i64::MAX, i64::MIN));
                                entry.0 = entry.0.min(s);
                                entry.1 = entry.1.max(e);
                            }
                        }
                        for (ch, (batch_min, batch_max)) in &batch_chrom_bounds {
                            let window_start = (batch_max / AnnotateProvider::SIFT_WINDOW_SIZE)
                                * AnnotateProvider::SIFT_WINDOW_SIZE;
                            let min_window_start = (batch_min / AnnotateProvider::SIFT_WINDOW_SIZE)
                                * AnnotateProvider::SIFT_WINDOW_SIZE;
                            let mut ws = min_window_start;
                            while ws <= window_start + AnnotateProvider::SIFT_WINDOW_SIZE {
                                let key = (ch.clone(), ws);
                                if !ann.loaded_sift_windows.contains(&key) {
                                    if let Some(ref direct) = ann.sift_direct {
                                        direct.load_window(
                                            ch,
                                            ws,
                                            ws + AnnotateProvider::SIFT_WINDOW_SIZE,
                                            &mut ann.sift_cache,
                                        )?;
                                    }
                                    ann.loaded_sift_windows.insert(key);
                                }
                                ws += AnnotateProvider::SIFT_WINDOW_SIZE;
                            }
                            ann.sift_cache.evict_before(*batch_min);
                        }
                    }
                }
            }

            #[cfg(not(feature = "kv-cache"))]
            let sift_kv: Option<()> = None;
            #[cfg(feature = "kv-cache")]
            let sift_kv = &ann.sift_kv;

            let annotated = ann.tmp_provider.annotate_batch_with_transcript_engine(
                batch,
                engine,
                &ctx,
                &ann.colocated_map,
                &mut ann.sift_cache,
                #[cfg(feature = "kv-cache")]
                sift_kv,
                #[cfg(not(feature = "kv-cache"))]
                &sift_kv,
                skip_csq,
                skip_typed_cols,
                &ann.config.flags,
                &ann.config.hgvs_flags,
                ann.config.transcript_selection,
                &ann.config.pick_flags,
                &mut ann.hgvs_reader,
            )?;

            if let Some(indices) = projection {
                out.push_back(annotated.project(indices)?);
            } else {
                out.push_back(annotated);
            }
        }
    }
    Ok(out)
}

impl Stream for ContigAnnotationStream {
    type Item = Result<RecordBatch>;

    fn poll_next(mut self: Pin<&mut Self>, cx: &mut TaskCtx<'_>) -> Poll<Option<Self::Item>> {
        loop {
            let rows_emitted = self.rows_emitted;
            let fetch_limit = self.config.fetch_limit;
            match &mut self.state {
                StreamState::Done => return Poll::Ready(None),

                StreamState::StartContig => {
                    // LIMIT pushdown: stop processing contigs if limit reached.
                    if fetch_limit.is_some_and(|limit| rows_emitted >= limit) {
                        self.state = StreamState::Done;
                        return Poll::Ready(None);
                    }
                    let Some(chrom) = self.contigs.pop_front() else {
                        // All contigs processed. Deregister the global KV
                        // variation table if fjall was used, via async cleanup
                        // future (safe on any Tokio runtime flavor).
                        #[cfg(feature = "kv-cache")]
                        if self.config.use_fjall {
                            let session = Arc::clone(&self.session);
                            let fut: CleanupFuture = Box::pin(async move {
                                crate::partitioned_cache::deregister_table(
                                    &session,
                                    "__vep_kv_variation",
                                )
                                .await
                                .ok();
                                Ok(())
                            });
                            // Transition to FinalCleanup which goes to Done
                            // (not StartContig, avoiding the infinite loop).
                            self.state = StreamState::FinalCleanup(fut);
                            continue;
                        }
                        self.state = StreamState::Done;
                        return Poll::Ready(None);
                    };
                    let session = Arc::clone(&self.session);
                    let cache = Arc::clone(&self.cache);
                    let config = self.config.clone();
                    let full_schema = self.full_schema.clone();

                    let fut: PrepareFuture = Box::pin(async move {
                        prepare_contig_context(session, cache, chrom, config, full_schema).await
                    });
                    self.state = StreamState::PreparingContig(fut);
                }

                StreamState::PreparingContig(fut) => match fut.as_mut().poll(cx) {
                    Poll::Pending => return Poll::Pending,
                    Poll::Ready(Err(e)) => {
                        self.state = StreamState::Done;
                        return Poll::Ready(Some(Err(e)));
                    }
                    Poll::Ready(Ok(None)) => {
                        // Contig skipped (no variation table).
                        self.state = StreamState::StartContig;
                    }
                    Poll::Ready(Ok(Some(ready))) => {
                        let projected_schema = self.projected_schema.clone();
                        let full_schema = self.full_schema.clone();
                        let session = Arc::clone(&self.session);
                        let cache = Arc::clone(&self.cache);
                        let config = self.config.clone();

                        // Derive VCF-only schema for AnnotateProvider.
                        let vcf_field_count = full_schema
                            .fields()
                            .len()
                            .saturating_sub(config.annotation_column_count);
                        let vcf_only_schema =
                            Schema::new(full_schema.fields()[..vcf_field_count].to_vec());

                        let tmp_provider = match AnnotateProvider::new(
                            Arc::clone(&session),
                            config.vcf_table.clone(),
                            String::new(),
                            AnnotationBackend::Parquet,
                            config.options_json.clone(),
                            vcf_only_schema,
                        ) {
                            Ok(provider) => provider,
                            Err(e) => {
                                self.state = StreamState::Done;
                                return Poll::Ready(Some(Err(e)));
                            }
                        };
                        let engine = TranscriptConsequenceEngine::new_with_hgvs_shift(
                            config.upstream_distance,
                            config.downstream_distance,
                            config.hgvs_flags.shift_hgvs,
                        );
                        let hgvs_reader = config.reference_fasta_path.as_deref().and_then(|path| {
                            fasta::io::indexed_reader::Builder::default()
                                .build_from_path(path)
                                .ok()
                        });
                        // SIFT source: when use_fjall, use SiftKvStore from fjall
                        // for lazy per-transcript lookups; otherwise use parquet
                        // SiftDirectReader.
                        #[cfg(feature = "kv-cache")]
                        let use_fjall_sift = config.use_fjall;
                        #[cfg(not(feature = "kv-cache"))]
                        let use_fjall_sift = false;

                        let sift_direct: Option<SiftDirectReader> = if config.flags.everything
                            && !use_fjall_sift
                        {
                            config
                                .translations_sift_table
                                .as_deref()
                                .and_then(|table| {
                                    let path = std::path::Path::new(table);
                                    if path.exists() {
                                        AnnotateProvider::build_sift_direct_reader(table)
                                    } else {
                                        None
                                    }
                                })
                                .or_else(|| {
                                    cache
                                        .context_path("translation_sift", &ready.chrom)
                                        .and_then(|p| {
                                            AnnotateProvider::build_sift_direct_reader(p.to_str()?)
                                        })
                                })
                        } else {
                            None
                        };

                        // Reuse the pre-opened SiftKvStore from config (opened
                        // once, shared across contigs) rather than re-opening
                        // the fjall DB every contig.
                        #[cfg(feature = "kv-cache")]
                        let sift_kv = if use_fjall_sift && config.flags.everything {
                            config.sift_kv_store.clone()
                        } else {
                            None
                        };

                        self.state = StreamState::AnnotatingContig(ContigAnnotationState {
                            lookup_stream: Some(ready.lookup_stream),
                            transcript_cache_regions: ready
                                .transcripts
                                .iter()
                                .map(|tx| (tx.transcript_id.clone(), transcript_cache_regions(tx)))
                                .collect(),
                            persisted_buffer_transcripts: HashMap::new(),
                            transcripts: ready.transcripts,
                            translateable_seq_by_tx: ready.translateable_seq_by_tx,
                            exons: ready.exons,
                            translations: ready.translations,
                            regulatory: ready.regulatory,
                            motifs: ready.motifs,
                            mirnas: ready.mirnas,
                            structural: ready.structural,
                            colocated_map: HashMap::new(),
                            colocated_map_built: false,
                            coloc_sink: ready.coloc_sink,
                            hydrated_cds_tx_ids: HashSet::new(),
                            hgvs_reader,
                            sift_cache: SiftPolyphenCache::new(),
                            sift_direct,
                            loaded_sift_windows: HashSet::new(),
                            #[cfg(feature = "kv-cache")]
                            sift_kv,
                            tmp_provider,
                            engine,
                            window_buffer: Vec::with_capacity(HYDRATION_WINDOW_SIZE),
                            lookup_done: false,
                            ephemeral_tables: ready.ephemeral_tables,
                            chrom: ready.chrom,
                            config,
                            session,
                            cache,
                            t_contig: ready.t_contig,
                            contig_rows: 0,
                        });
                    }
                },

                StreamState::AnnotatingContig(ann) => {
                    // Pull batches from lookup stream into window buffer.
                    //
                    // Both VariantLookupExec (parquet) and KvLookupExec (fjall)
                    // buffer matched rows internally and only emit after the
                    // input is exhausted, ensuring the colocated sink is fully
                    // populated when the first batch arrives here.
                    //
                    // LIMIT pushdown: once we have enough buffered rows to
                    // satisfy the limit, stop pulling from the lookup stream
                    // to avoid unnecessary annotation work.
                    let buffered_rows: usize = ann.window_buffer.iter().map(|b| b.num_rows()).sum();
                    let limit_buffered =
                        fetch_limit.is_some_and(|limit| rows_emitted + buffered_rows >= limit);
                    if !ann.lookup_done && !limit_buffered {
                        let stream = ann
                            .lookup_stream
                            .as_mut()
                            .expect("stream alive during lookup");
                        match stream.as_mut().poll_next(cx) {
                            Poll::Pending => return Poll::Pending,
                            Poll::Ready(Some(Ok(batch))) => {
                                ann.contig_rows += batch.num_rows();
                                ann.window_buffer.push(batch);
                                if ann.window_buffer.len() < HYDRATION_WINDOW_SIZE {
                                    continue; // Keep filling window.
                                }
                            }
                            Poll::Ready(Some(Err(e))) => {
                                let session = Arc::clone(&ann.session);
                                let tables = ann.ephemeral_tables.clone();
                                self.state = StreamState::ErrorCleaningUp(
                                    make_cleanup_future(session, tables),
                                    e,
                                );
                                continue;
                            }
                            Poll::Ready(None) => {
                                ann.lookup_done = true;
                                // Drop the lookup stream to reclaim BuildSide
                                // memory (COITrees, hash indices, concatenated
                                // VCF batch) — no longer needed.
                                ann.lookup_stream = None;
                                // NOTE: colocated sink is cleared AFTER the
                                // colocated map is built (below), not here.
                                // Clearing here would lose data for backends
                                // like fjall that emit batches immediately.
                                profile_end!(
                                    &format!("{}: 1. variation_lookup", ann.chrom),
                                    ann.t_contig,
                                    format!("{} rows", ann.contig_rows)
                                );
                            }
                        }
                    }

                    // Process window (full or final partial).
                    // Take ownership of state to work with it.
                    let StreamState::AnnotatingContig(mut ann) =
                        std::mem::replace(&mut self.state, StreamState::Done)
                    else {
                        unreachable!()
                    };

                    // LIMIT pushdown: skip remaining windows if limit reached.
                    let limit_reached = fetch_limit.is_some_and(|limit| rows_emitted >= limit);

                    if ann.window_buffer.is_empty() || limit_reached {
                        // No more data (or limit reached) — clean up.
                        // Drop heavy state eagerly before the async cleanup future runs.
                        profile_end!(
                            &format!("{}: TOTAL", ann.chrom),
                            ann.t_contig,
                            format!("{} rows", ann.contig_rows)
                        );
                        if profiling_enabled() {
                            eprintln!("[VEP_PROFILE] ------ contig {} END ------", ann.chrom);
                        }
                        // Eagerly reclaim per-contig memory.
                        ann.colocated_map = HashMap::new();
                        ann.transcripts = Vec::new();
                        ann.exons = Vec::new();
                        ann.translations = Vec::new();
                        ann.regulatory = Vec::new();
                        ann.motifs = Vec::new();
                        let fut = make_cleanup_future(
                            Arc::clone(&ann.session),
                            std::mem::take(&mut ann.ephemeral_tables),
                        );
                        self.state = StreamState::CleaningUp(fut);
                        continue;
                    }

                    // Build colocated map once (sink is fully populated now
                    // that the entire lookup stream has been drained).
                    if !ann.colocated_map_built {
                        if ann.config.flags.check_existing {
                            let mut guard = ann.coloc_sink.lock().unwrap();
                            ann.colocated_map = build_colocated_map_from_sink(&guard);
                            // Clear sink now that data has been copied to the map.
                            guard.clear();
                        }
                        ann.colocated_map_built = true;
                    }

                    // Take one window's worth of batches from the buffer.
                    let window_end = ann.window_buffer.len().min(HYDRATION_WINDOW_SIZE);
                    let window_batches: Vec<RecordBatch> =
                        ann.window_buffer.drain(..window_end).collect();

                    // Window-based HGVS hydration (like SIFT sliding window).
                    if let Err(e) = hydrate_window(
                        &mut ann.transcripts,
                        &ann.exons,
                        &mut ann.translations,
                        &ann.translateable_seq_by_tx,
                        &mut ann.hgvs_reader,
                        &mut ann.hydrated_cds_tx_ids,
                        &window_batches,
                    ) {
                        let fut = make_cleanup_future(
                            Arc::clone(&ann.session),
                            std::mem::take(&mut ann.ephemeral_tables),
                        );
                        self.state = StreamState::ErrorCleaningUp(fut, e);
                        continue;
                    }

                    // Annotate window.
                    let projection = ann.config.projection.clone();
                    match annotate_window(&mut ann, &window_batches, projection.as_deref()) {
                        Err(e) => {
                            let fut = make_cleanup_future(
                                Arc::clone(&ann.session),
                                std::mem::take(&mut ann.ephemeral_tables),
                            );
                            self.state = StreamState::ErrorCleaningUp(fut, e);
                            continue;
                        }
                        Ok(batches) => {
                            self.state = StreamState::DrainingWindow {
                                batches,
                                annotation_state: ann,
                            };
                        }
                    }
                }

                StreamState::DrainingWindow {
                    batches,
                    annotation_state: _,
                } => {
                    if let Some(batch) = batches.pop_front() {
                        // LIMIT pushdown: truncate or stop if we've reached the limit.
                        if let Some(limit) = fetch_limit {
                            let remaining = limit.saturating_sub(rows_emitted);
                            if remaining == 0 {
                                // Don't emit — fall through to the window-drained
                                // path below, which transitions to AnnotatingContig.
                                // The limit_reached check there triggers proper
                                // cleanup (deregister ephemeral tables).
                                continue;
                            }
                            if batch.num_rows() > remaining {
                                self.rows_emitted += remaining;
                                return Poll::Ready(Some(Ok(batch.slice(0, remaining))));
                            }
                        }
                        self.rows_emitted += batch.num_rows();
                        return Poll::Ready(Some(Ok(batch)));
                    }
                    // Window drained — back to AnnotatingContig for next window
                    // (or cleanup if window_buffer is empty).
                    let StreamState::DrainingWindow {
                        annotation_state: ann,
                        ..
                    } = std::mem::replace(&mut self.state, StreamState::Done)
                    else {
                        unreachable!()
                    };
                    self.state = StreamState::AnnotatingContig(ann);
                }

                StreamState::CleaningUp(fut) => match fut.as_mut().poll(cx) {
                    Poll::Pending => return Poll::Pending,
                    Poll::Ready(Err(e)) => {
                        self.state = StreamState::Done;
                        return Poll::Ready(Some(Err(e)));
                    }
                    Poll::Ready(Ok(())) => {
                        self.state = StreamState::StartContig;
                    }
                },

                StreamState::ErrorCleaningUp(fut, _) => match fut.as_mut().poll(cx) {
                    Poll::Pending => return Poll::Pending,
                    Poll::Ready(_) => {
                        // Cleanup done (ignore cleanup errors) — propagate original error.
                        let StreamState::ErrorCleaningUp(_, original_err) =
                            std::mem::replace(&mut self.state, StreamState::Done)
                        else {
                            unreachable!()
                        };
                        return Poll::Ready(Some(Err(original_err)));
                    }
                },

                StreamState::FinalCleanup(fut) => match fut.as_mut().poll(cx) {
                    Poll::Pending => return Poll::Pending,
                    Poll::Ready(_) => {
                        self.state = StreamState::Done;
                        return Poll::Ready(None);
                    }
                },
            }
        }
    }
}

// ---------------------------------------------------------------------------
// Per-contig setup: parallel context loading + lookup stream creation
// ---------------------------------------------------------------------------

/// Register ephemeral tables, load context data, and create the variation
/// lookup stream — all in preparation for window-based streaming annotation.
///
/// Context loading completes before the lookup stream is first polled; the
/// actual lookup I/O (build + probe) happens lazily on `poll_next`.
/// Returns `None` if the contig has no variation table (skip).
async fn prepare_contig_context(
    session: Arc<SessionContext>,
    cache: Arc<PartitionedParquetCache>,
    chrom: String,
    config: ContigAnnotationConfig,
    full_schema: SchemaRef,
) -> Result<Option<ContigReadyState>> {
    let t_contig = profile_start!();
    if profiling_enabled() {
        eprintln!("[VEP_PROFILE] ------ contig {chrom} START ------");
    }

    // Derive VCF-only schema.
    let vcf_field_count = full_schema
        .fields()
        .len()
        .saturating_sub(config.annotation_column_count);
    let vcf_only_schema = Schema::new(full_schema.fields()[..vcf_field_count].to_vec());

    // 1. Register ALL ephemeral tables upfront (variation + context).
    let mut ephemeral_tables: Vec<String> = Vec::new();

    // Variation table: either per-chrom parquet or global fjall KV store.
    #[cfg(feature = "kv-cache")]
    let use_fjall = config.use_fjall;
    #[cfg(not(feature = "kv-cache"))]
    let use_fjall = false;

    let var_table = if use_fjall {
        #[cfg(feature = "kv-cache")]
        {
            // Register the shared fjall KV store as a table (idempotent).
            let kv_table_name = "__vep_kv_variation".to_string();
            if !session.table_exist(&kv_table_name)? {
                let store = config
                    .kv_store
                    .as_ref()
                    .expect("kv_store must be set when use_fjall=true");
                let kv_provider = KvCacheTableProvider::from_store(Arc::clone(store));
                session.register_table(&kv_table_name, Arc::new(kv_provider))?;
            }
            // Don't add to ephemeral_tables — the global KV table persists
            // across contigs and is deregistered after the last contig.
            kv_table_name
        }
        #[cfg(not(feature = "kv-cache"))]
        {
            unreachable!("use_fjall requires kv-cache feature")
        }
    } else {
        let var_table =
            crate::partitioned_cache::register_chrom_parquet(&session, &cache, "variation", &chrom)
                .await?;
        let Some(var_table) = var_table else {
            if profiling_enabled() {
                eprintln!("[VEP_PROFILE] ------ contig {chrom} SKIP (no variation) ------");
            }
            return Ok(None);
        };
        ephemeral_tables.push(var_table.clone());
        var_table
    };

    let tx_table =
        crate::partitioned_cache::register_chrom_parquet(&session, &cache, "transcript", &chrom)
            .await?;
    if let Some(ref t) = tx_table {
        ephemeral_tables.push(t.clone());
    }
    let ex_table =
        crate::partitioned_cache::register_chrom_parquet(&session, &cache, "exon", &chrom).await?;
    if let Some(ref t) = ex_table {
        ephemeral_tables.push(t.clone());
    }
    let tl_table = crate::partitioned_cache::register_chrom_parquet(
        &session,
        &cache,
        "translation_core",
        &chrom,
    )
    .await?;
    if let Some(ref t) = tl_table {
        ephemeral_tables.push(t.clone());
    }
    let rg_table =
        crate::partitioned_cache::register_chrom_parquet(&session, &cache, "regulatory", &chrom)
            .await?;
    if let Some(ref t) = rg_table {
        ephemeral_tables.push(t.clone());
    }
    let mt_table =
        crate::partitioned_cache::register_chrom_parquet(&session, &cache, "motif", &chrom).await?;
    if let Some(ref t) = mt_table {
        ephemeral_tables.push(t.clone());
    }

    // 2. Create lookup stream + load context data.
    let worklist = MissWorklist::for_chrom(&chrom);

    // Lookup arm: build LookupProvider, create stream (cheap — build+probe
    // happens on first poll, NOT here).
    let coloc_sink: ColocatedSink = Arc::new(Mutex::new(HashMap::new()));
    let vcf_schema = session
        .table(&config.vcf_table)
        .await?
        .schema()
        .as_arrow()
        .clone();
    let cache_schema = session.table(&var_table).await?.schema().as_arrow().clone();
    let mut provider = LookupProvider::new(
        Arc::clone(&session),
        config.vcf_table.clone(),
        var_table.clone(),
        vcf_schema,
        cache_schema,
        config.cache_columns.clone(),
        config.extended_probes,
        config.allowed_failed,
        None, // reference_fasta_path is for HGVS hydration, not lookup
    )?;
    provider.set_vcf_filter(Some(col("chrom").eq(lit(&*chrom))));
    if config.flags.check_existing {
        provider.set_colocated_sink(Arc::clone(&coloc_sink));
    }
    let session_state = session.state();
    let plan = provider.scan(&session_state, None, &[], None).await?;
    let lookup_stream = plan.execute(0, session.task_ctx())?;

    // Context arm: load transcripts, exons, translations, etc.
    let t_ctx = profile_start!();
    let tmp_provider = AnnotateProvider::new(
        Arc::clone(&session),
        config.vcf_table.clone(),
        String::new(),
        AnnotationBackend::Parquet,
        config.options_json.clone(),
        vcf_only_schema,
    )?;

    let tx = if let Some(ref table) = tx_table {
        let (tx, seq) = tmp_provider.load_transcripts(table, &worklist).await?;
        let filtered: Vec<_> = tx
            .into_iter()
            .filter(|t| passes_transcript_selection(t, config.transcript_selection))
            .collect();
        (filtered, seq)
    } else {
        (Vec::new(), HashMap::new())
    };
    let (tx_vec, translateable_seq) = tx;
    let tx_ids: HashSet<String> = tx_vec.iter().map(|t| t.transcript_id.clone()).collect();

    let ex = if let Some(ref table) = ex_table {
        let raw = tmp_provider.load_exons(table, &worklist).await?;
        raw.into_iter()
            .filter(|e| tx_ids.contains(&e.transcript_id))
            .collect()
    } else {
        Vec::new()
    };
    let tl = if let Some(ref table) = tl_table {
        let raw = tmp_provider.load_translations(table, &worklist).await?;
        raw.into_iter()
            .filter(|t| tx_ids.contains(&t.transcript_id))
            .collect()
    } else {
        Vec::new()
    };
    let rg = if let Some(ref table) = rg_table {
        tmp_provider
            .load_regulatory_features(table, &worklist)
            .await?
    } else {
        Vec::new()
    };
    let mt = if let Some(ref table) = mt_table {
        tmp_provider.load_motif_features(table, &worklist).await?
    } else {
        Vec::new()
    };
    profile_end!(&format!("{chrom}: context_load"), t_ctx);

    Ok(Some(ContigReadyState {
        lookup_stream,
        coloc_sink,
        transcripts: tx_vec,
        translateable_seq_by_tx: translateable_seq,
        exons: ex,
        translations: tl,
        regulatory: rg,
        motifs: mt,
        // TODO: miRNA and structural features are not yet partitioned —
        // these are rare and handled by the monolithic path only.
        mirnas: Vec::new(),
        structural: Vec::new(),
        ephemeral_tables,
        chrom,
        t_contig,
    }))
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
        limit: Option<usize>,
    ) -> Result<Arc<dyn ExecutionPlan>> {
        let _store = build_store(self.backend, self.cache_source.clone());

        // Parse use_fjall option — when true, use fjall KV store for variation
        // lookup + SIFT while keeping context from partitioned parquet.
        #[cfg(feature = "kv-cache")]
        let use_fjall = self
            .options_json
            .as_deref()
            .and_then(|opts| Self::parse_json_bool_option(opts, "use_fjall"))
            .unwrap_or(false);
        #[cfg(not(feature = "kv-cache"))]
        let use_fjall = false;

        // Check for partitioned per-chromosome cache layout.
        // Opt-in/out via "partitioned": true/false in options_json.
        // Both parquet-only and fjall paths require partitioned context parquet.
        let partitioned_opt = self
            .options_json
            .as_deref()
            .and_then(|opts| Self::parse_json_bool_option(opts, "partitioned"));
        let partitioned_cache = if partitioned_opt != Some(false) {
            PartitionedParquetCache::detect(&self.cache_source)
        } else {
            None
        };
        // When explicitly requested or auto-detected, use partitioned path.
        if let Some(ref cache) = partitioned_cache {
            if profiling_enabled() {
                eprintln!(
                    "[VEP_PROFILE] detected partitioned cache: {} chroms in {}{}",
                    cache.available_chroms().len(),
                    self.cache_source,
                    if use_fjall {
                        " [fjall variation+sift]"
                    } else {
                        ""
                    },
                );
            }

            // Determine requested cache columns.
            // When using fjall, get schema from the KV store; otherwise from
            // a sample variation parquet file.
            #[cfg(feature = "kv-cache")]
            let kv_store_arc: Option<Arc<crate::kv_cache::VepKvStore>> = if use_fjall {
                let fjall_path = std::path::Path::new(&self.cache_source).join("variation.fjall");
                if !fjall_path.exists() {
                    return Err(DataFusionError::Execution(format!(
                        "annotate_vep(): use_fjall=true but no fjall store found at '{}'",
                        fjall_path.display()
                    )));
                }
                Some(Arc::new(crate::kv_cache::VepKvStore::open(&fjall_path)?))
            } else {
                None
            };

            let (available_cache_columns, sample_table_to_deregister) = if use_fjall {
                #[cfg(feature = "kv-cache")]
                {
                    let store = kv_store_arc.as_ref().unwrap();
                    let cols: HashSet<String> = store
                        .schema()
                        .fields()
                        .iter()
                        .map(|f| f.name().clone())
                        .collect();
                    (cols, None)
                }
                #[cfg(not(feature = "kv-cache"))]
                {
                    unreachable!("use_fjall requires kv-cache feature")
                }
            } else {
                let sample_chrom = &cache.available_chroms()[0];
                let sample_table = crate::partitioned_cache::register_chrom_parquet(
                    &self.session,
                    cache,
                    "variation",
                    sample_chrom,
                )
                .await?;
                let sample_table = sample_table.ok_or_else(|| {
                    DataFusionError::Execution(
                        "partitioned cache: no variation parquet for sample chrom".to_string(),
                    )
                })?;
                let cache_schema = self
                    .session
                    .table(&sample_table)
                    .await?
                    .schema()
                    .as_arrow()
                    .clone();
                let cols: HashSet<String> = cache_schema
                    .fields()
                    .iter()
                    .map(|f| f.name().clone())
                    .collect();
                (cols, Some(sample_table))
            };

            let mut preferred_columns = vec!["consequence_types", "most_severe_consequence"];
            for c in cache_lookup_column_names() {
                if !preferred_columns.contains(&c) {
                    preferred_columns.push(c);
                }
            }
            let requested_columns: Vec<&str> = preferred_columns
                .iter()
                .copied()
                .filter(|name| available_cache_columns.contains(*name))
                .collect();

            let extended_probes = self
                .options_json
                .as_deref()
                .and_then(|opts| Self::parse_json_bool_option(opts, "extended_probes"))
                .unwrap_or(true);

            let translations_sift_table = self
                .options_json
                .as_deref()
                .and_then(|opts| Self::parse_json_string_option(opts, "translations_sift_table"));

            // Deregister sample table (will be re-registered per contig).
            if let Some(ref tbl) = sample_table_to_deregister {
                crate::partitioned_cache::deregister_table(&self.session, tbl).await?;
            }

            return self
                .scan_with_transcript_engine_partitioned(
                    state,
                    projection,
                    &requested_columns,
                    extended_probes,
                    cache,
                    translations_sift_table.as_deref(),
                    #[cfg(feature = "kv-cache")]
                    kv_store_arc,
                    limit,
                )
                .await;
        }

        Err(DataFusionError::Execution(format!(
            "annotate_vep(): no partitioned cache detected at '{}'. \
             Expected a directory with a variation/ subdirectory containing per-chromosome parquet files.",
            self.cache_source
        )))
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::transcript_consequence::{
        CachedPredictions, FeatureType, ProteinDomainFeature, SiftPolyphenCache, TranslationFeature,
    };

    // ── format_appris ──────────────────────────────────────────────────

    #[test]
    fn test_format_appris_principal1() {
        assert_eq!(format_appris("principal1"), "P1");
    }

    #[test]
    fn test_format_appris_alternative2() {
        assert_eq!(format_appris("alternative2"), "A2");
    }

    #[test]
    fn test_format_appris_principal5() {
        assert_eq!(format_appris("principal5"), "P5");
    }

    #[test]
    fn test_format_appris_passthrough() {
        assert_eq!(format_appris("other"), "other");
    }

    #[test]
    fn test_pick_flags_parse_custom_order() {
        let flags = PickFlags::from_options_json(Some(
            "{\"flag_pick_allele_gene\":true,\"pick_order\":\"biotype,rank,mane_select,tsl,canonical,appris,ccds,length\"}",
        ))
        .expect("pick flags should parse");

        assert!(flags.flag_pick_allele_gene);
        assert_eq!(
            flags.pick_order,
            vec![
                PickCriterion::Biotype,
                PickCriterion::Rank,
                PickCriterion::ManeSelect,
                PickCriterion::Tsl,
                PickCriterion::Canonical,
                PickCriterion::Appris,
                PickCriterion::Ccds,
                PickCriterion::Length,
            ]
        );
    }

    #[test]
    fn test_pick_flags_reject_invalid_pick_order_criterion() {
        let err = PickFlags::from_options_json(Some(
            "{\"flag_pick_allele_gene\":true,\"pick_order\":\"mane_select,bogus\"}",
        ))
        .expect_err("invalid pick_order should fail")
        .to_string();

        assert!(err.contains("unsupported pick_order criterion 'bogus'"));
    }

    #[test]
    fn test_append_pick_flag_preserves_existing_flags_and_avoids_duplicates() {
        let mut flags = Some("cds_end_NF".to_string());
        append_pick_flag(&mut flags);
        append_pick_flag(&mut flags);
        assert_eq!(flags, Some("cds_end_NF&PICK".to_string()));

        let mut empty = None;
        append_pick_flag(&mut empty);
        assert_eq!(empty, Some("PICK".to_string()));
    }

    #[test]
    fn test_parse_appris_pick_rank_matches_release_115() {
        assert_eq!(parse_appris_pick_rank(Some("principal5")), 5);
        assert_eq!(parse_appris_pick_rank(Some("alternative2")), 12);
        assert_eq!(parse_appris_pick_rank(Some("other")), 100);
        assert_eq!(parse_appris_pick_rank(None), 100);
    }

    #[test]
    fn test_mark_flag_pick_allele_gene_requires_exact_source_match() {
        let mut tx_havana = make_tx("ENST00000001", Some("GENE1"), Some("HGNC"), None);
        tx_havana.gene_stable_id = Some("ENSG00000001".to_string());
        tx_havana.biotype = "protein_coding".to_string();
        tx_havana.source = Some("ensembl_havana".to_string());

        let mut tx_ensembl = make_tx("ENST00000002", Some("GENE1"), Some("HGNC"), None);
        tx_ensembl.gene_stable_id = Some("ENSG00000001".to_string());
        tx_ensembl.biotype = "protein_coding".to_string();
        tx_ensembl.source = Some("Ensembl".to_string());

        let transcripts = vec![tx_havana, tx_ensembl];
        let ctx = PreparedContext::new(&transcripts, &[], &[], &[], &[], &[], &[]);
        let mut assignments = vec![
            TranscriptConsequence {
                transcript_idx: Some(0),
                feature_type: FeatureType::Transcript,
                ..Default::default()
            },
            TranscriptConsequence {
                transcript_idx: Some(1),
                feature_type: FeatureType::Transcript,
                ..Default::default()
            },
        ];

        mark_flag_pick_allele_gene(
            &mut assignments,
            &ctx,
            &PickFlags {
                flag_pick_allele_gene: true,
                pick_order: vec![PickCriterion::Ensembl],
            },
        );

        assert_eq!(assignments[0].flags, None);
        assert_eq!(assignments[1].flags, Some("PICK".to_string()));
    }

    #[test]
    fn test_mark_flag_pick_allele_gene_prefers_cds_length_over_protein_length() {
        let mut tx_shorter_cds = make_tx("ENST00000011", Some("GENE1"), Some("HGNC"), None);
        tx_shorter_cds.gene_stable_id = Some("ENSG00000011".to_string());
        tx_shorter_cds.biotype = "protein_coding".to_string();

        let mut tx_longer_cds = make_tx("ENST00000012", Some("GENE1"), Some("HGNC"), None);
        tx_longer_cds.gene_stable_id = Some("ENSG00000011".to_string());
        tx_longer_cds.biotype = "protein_coding".to_string();

        let transcripts = vec![tx_shorter_cds, tx_longer_cds];
        let translations = vec![
            TranslationFeature {
                transcript_id: "ENST00000011".to_string(),
                cds_len: Some(90),
                protein_len: Some(100),
                translation_seq: None,
                cds_sequence: Some("A".repeat(90)),
                stable_id: None,
                version: None,
                protein_features: Vec::new(),
            },
            TranslationFeature {
                transcript_id: "ENST00000012".to_string(),
                cds_len: Some(99),
                protein_len: Some(20),
                translation_seq: None,
                cds_sequence: Some("A".repeat(99)),
                stable_id: None,
                version: None,
                protein_features: Vec::new(),
            },
        ];
        let ctx = PreparedContext::new(&transcripts, &[], &translations, &[], &[], &[], &[]);
        let mut assignments = vec![
            TranscriptConsequence {
                transcript_idx: Some(0),
                feature_type: FeatureType::Transcript,
                ..Default::default()
            },
            TranscriptConsequence {
                transcript_idx: Some(1),
                feature_type: FeatureType::Transcript,
                ..Default::default()
            },
        ];

        mark_flag_pick_allele_gene(
            &mut assignments,
            &ctx,
            &PickFlags {
                flag_pick_allele_gene: true,
                pick_order: vec![PickCriterion::Length],
            },
        );

        assert_eq!(assignments[0].flags, None);
        assert_eq!(assignments[1].flags, Some("PICK".to_string()));
    }

    #[test]
    fn test_mark_flag_pick_allele_gene_requires_gene_stable_id() {
        let mut tx_a = make_tx("ENST00000021", Some("GENE1"), Some("HGNC"), None);
        tx_a.biotype = "protein_coding".to_string();
        tx_a.is_canonical = false;

        let mut tx_b = make_tx("ENST00000022", Some("GENE1"), Some("HGNC"), None);
        tx_b.biotype = "protein_coding".to_string();
        tx_b.is_canonical = true;

        let transcripts = vec![tx_a, tx_b];
        let ctx = PreparedContext::new(&transcripts, &[], &[], &[], &[], &[], &[]);
        let mut assignments = vec![
            TranscriptConsequence {
                transcript_idx: Some(0),
                feature_type: FeatureType::Transcript,
                ..Default::default()
            },
            TranscriptConsequence {
                transcript_idx: Some(1),
                feature_type: FeatureType::Transcript,
                ..Default::default()
            },
        ];

        mark_flag_pick_allele_gene(
            &mut assignments,
            &ctx,
            &PickFlags {
                flag_pick_allele_gene: true,
                pick_order: vec![PickCriterion::Canonical],
            },
        );

        assert_eq!(assignments[0].flags, None);
        assert_eq!(assignments[1].flags, None);
    }

    // ── format_prediction ──────────────────────────────────────────────

    #[test]
    fn test_format_prediction_deleterious() {
        assert_eq!(format_prediction("deleterious", 0.01), "deleterious(0.01)");
    }

    #[test]
    fn test_format_prediction_probably_damaging() {
        assert_eq!(
            format_prediction("probably damaging", 0.999),
            "probably_damaging(0.999)"
        );
    }

    #[test]
    fn test_format_prediction_tolerated_low_confidence() {
        assert_eq!(
            format_prediction("tolerated - low confidence", 0.23),
            "tolerated_low_confidence(0.23)"
        );
    }

    // ── mirna_structure_field ──────────────────────────────────────────

    #[test]
    fn test_mirna_structure_field_preserves_distinct_stem_sides() {
        assert_eq!(
            mirna_structure_field(Some("(.)."), "miRNA", Some(1), Some(4)),
            "miRNA_loop&miRNA_stem&miRNA_stem"
        );
    }

    #[test]
    fn test_mirna_structure_field_open_and_close_stem_only() {
        assert_eq!(
            mirna_structure_field(Some("()"), "miRNA", Some(1), Some(2)),
            "miRNA_stem&miRNA_stem"
        );
    }

    #[test]
    fn test_mirna_structure_field_non_mirna_empty() {
        assert!(mirna_structure_field(Some("(.)."), "lncRNA", Some(1), Some(4)).is_empty());
    }

    // ── lookup_sift_polyphen ───────────────────────────────────────────

    fn make_sift_cache(entries: Vec<(&str, i32, &str, &str, f32, &str, f32)>) -> SiftPolyphenCache {
        let mut cache = SiftPolyphenCache::new();
        let mut by_tx: HashMap<String, CachedPredictions> = HashMap::new();
        for (tx_id, pos, aa, sift_pred, sift_score, pp_pred, pp_score) in entries {
            let preds = by_tx.entry(tx_id.to_string()).or_default();
            if let Some(aa_idx) = CompactPrediction::encode_amino_acid(aa) {
                preds.sift.push(CompactPrediction {
                    position: pos,
                    amino_acid: aa_idx,
                    prediction: CompactPrediction::encode_prediction(sift_pred),
                    score: sift_score,
                });
                preds.polyphen.push(CompactPrediction {
                    position: pos,
                    amino_acid: aa_idx,
                    prediction: CompactPrediction::encode_prediction(pp_pred),
                    score: pp_score,
                });
            }
        }
        for (tx_id, mut preds) in by_tx {
            preds.sort();
            cache.insert(tx_id, preds, i64::MAX);
        }
        cache
    }

    #[test]
    fn test_lookup_sift_polyphen_single_aa_match() {
        let mut cache = make_sift_cache(vec![(
            "ENST00000001",
            42,
            "I",
            "deleterious",
            0.01,
            "probably damaging",
            0.999,
        )]);
        let (sift, polyphen) = lookup_sift_polyphen(
            Some("ENST00000001"),
            Some("42"),
            Some("V/I"),
            &mut cache,
            &None,
        );
        assert_eq!(sift, "deleterious(0.01)");
        assert_eq!(polyphen, "probably_damaging(0.999)");
    }

    #[test]
    fn test_lookup_sift_polyphen_non_substitution_skipped() {
        let mut cache = make_sift_cache(vec![(
            "ENST00000001",
            42,
            "I",
            "deleterious",
            0.01,
            "benign",
            0.0,
        )]);
        // Multi-char alt amino acid — not a single substitution.
        let (sift, polyphen) = lookup_sift_polyphen(
            Some("ENST00000001"),
            Some("42"),
            Some("V/IL"),
            &mut cache,
            &None,
        );
        assert!(sift.is_empty());
        assert!(polyphen.is_empty());

        // Range position — indel, should be skipped.
        let (sift, polyphen) = lookup_sift_polyphen(
            Some("ENST00000001"),
            Some("42-43"),
            Some("V/I"),
            &mut cache,
            &None,
        );
        assert!(sift.is_empty());
        assert!(polyphen.is_empty());
    }

    #[test]
    fn test_lookup_sift_polyphen_missing_transcript() {
        let mut cache = SiftPolyphenCache::new();
        let (sift, polyphen) = lookup_sift_polyphen(
            Some("ENST_MISSING"),
            Some("42"),
            Some("V/I"),
            &mut cache,
            &None,
        );
        assert!(sift.is_empty());
        assert!(polyphen.is_empty());
    }

    // ── lookup_domains ─────────────────────────────────────────────────

    fn make_translation(
        tx_id: &str,
        protein_features: Vec<ProteinDomainFeature>,
    ) -> TranslationFeature {
        TranslationFeature {
            transcript_id: tx_id.to_string(),
            cds_len: None,
            protein_len: None,
            translation_seq: None,
            cds_sequence: None,
            translation_seq_canonical: None,
            cds_sequence_canonical: None,
            stable_id: None,
            version: None,
            protein_features,
        }
    }

    fn minimal_ctx(translations: &[TranslationFeature]) -> PreparedContext<'_> {
        PreparedContext::new(&[], &[], translations, &[], &[], &[], &[])
    }

    #[test]
    fn test_lookup_domains_single_overlap() {
        let translations = vec![make_translation(
            "ENST00000001",
            vec![ProteinDomainFeature {
                analysis: Some("Pfam".to_string()),
                hseqname: Some("PF00069".to_string()),
                start: 30,
                end: 100,
            }],
        )];
        let ctx = minimal_ctx(&translations);
        assert_eq!(
            lookup_domains(Some("ENST00000001"), Some("42"), None, &ctx),
            "Pfam:PF00069"
        );
    }

    #[test]
    fn test_lookup_domains_multiple_overlap() {
        let translations = vec![make_translation(
            "ENST00000001",
            vec![
                ProteinDomainFeature {
                    analysis: Some("Pfam".to_string()),
                    hseqname: Some("PF00069".to_string()),
                    start: 30,
                    end: 100,
                },
                ProteinDomainFeature {
                    analysis: Some("PROSITE profiles".to_string()),
                    hseqname: Some("PS50011".to_string()),
                    start: 40,
                    end: 50,
                },
            ],
        )];
        let ctx = minimal_ctx(&translations);
        assert_eq!(
            lookup_domains(Some("ENST00000001"), Some("42"), None, &ctx),
            "Pfam:PF00069&PROSITE_profiles:PS50011"
        );
    }

    #[test]
    fn test_lookup_domains_no_overlap() {
        let translations = vec![make_translation(
            "ENST00000001",
            vec![ProteinDomainFeature {
                analysis: Some("Pfam".to_string()),
                hseqname: Some("PF00069".to_string()),
                start: 100,
                end: 200,
            }],
        )];
        let ctx = minimal_ctx(&translations);
        assert!(lookup_domains(Some("ENST00000001"), Some("42"), None, &ctx).is_empty());
    }

    #[test]
    fn test_lookup_domains_spaces_replaced() {
        let translations = vec![make_translation(
            "ENST00000001",
            vec![ProteinDomainFeature {
                analysis: Some("Gene3D db".to_string()),
                hseqname: Some("1.10.510.10".to_string()),
                start: 1,
                end: 50,
            }],
        )];
        let ctx = minimal_ctx(&translations);
        assert_eq!(
            lookup_domains(Some("ENST00000001"), Some("10"), None, &ctx),
            "Gene3D_db:1.10.510.10"
        );
    }

    // ── HGVS_OFFSET sign ──────────────────────────────────────────────

    #[test]
    fn test_hgvs_offset_reverse_strand_is_negative() {
        // The HGVS_OFFSET logic in the CSQ builder negates shift_length for
        // reverse-strand transcripts (strand == -1).
        let shift_length: u32 = 3;
        let tx_strand: i8 = -1;
        let signed = shift_length as i64;
        let value = if tx_strand < 0 { -signed } else { signed };
        assert_eq!(value, -3);
    }

    #[test]
    fn test_hgvs_offset_forward_strand_is_positive() {
        let shift_length: u32 = 3;
        let tx_strand: i8 = 1;
        let signed = shift_length as i64;
        let value = if tx_strand < 0 { -signed } else { signed };
        assert_eq!(value, 3);
    }

    // ── buffer-local HGNC propagation ────────────────────────────────

    /// Minimal `TranscriptFeature` for HGNC propagation tests.
    fn make_tx(
        id: &str,
        gene_stable_id: Option<&str>,
        symbol: Option<&str>,
        symbol_source: Option<&str>,
        native_hgnc_id: Option<&str>,
    ) -> TranscriptFeature {
        TranscriptFeature {
            transcript_id: id.to_string(),
            chrom: "chr2".to_string(),
            start: 1,
            end: 100,
            strand: 1,
            biotype: "snoRNA".to_string(),
            cds_start: None,
            cds_end: None,
            cdna_coding_start: None,
            cdna_coding_end: None,
            cdna_mapper_segments: Vec::new(),
            mature_mirna_regions: Vec::new(),
            gene_stable_id: gene_stable_id.map(|s| s.to_string()),
            gene_symbol: symbol.map(|s| s.to_string()),
            gene_symbol_source: symbol_source.map(|s| s.to_string()),
            gene_hgnc_id_native: native_hgnc_id.map(|s| s.to_string()),
            gene_hgnc_id: native_hgnc_id.map(|s| s.to_string()),
            display_xref_id: None,
            source: None,
            refseq_match: None,
            refseq_edits: Vec::new(),
            is_gencode_basic: false,
            is_gencode_primary: false,
            bam_edit_status: None,
            has_non_polya_rna_edit: false,
            spliced_seq: None,
            five_prime_utr_seq: None,
            three_prime_utr_seq: None,
            translateable_seq: None,
            cdna_seq: None,
            version: None,
            cds_start_nf: false,
            cds_end_nf: false,
            flags_str: None,
            is_canonical: false,
            tsl: None,
            mane_select: None,
            mane_plus_clinical: None,
            translation_stable_id: None,
            gene_phenotype: false,
            ccds: None,
            swissprot: None,
            trembl: None,
            uniparc: None,
            uniprot_isoform: None,
            appris: None,
            ncrna_structure: None,
        }
    }

    fn make_selection_tx(id: &str, source: Option<&str>) -> TranscriptFeature {
        let mut tx = make_tx(id, None, None, None, None);
        tx.chrom = "1".to_string();
        tx.biotype = "protein_coding".to_string();
        tx.source = source.map(str::to_string);
        tx
    }

    #[test]
    fn edited_refseq_translation_cds_from_spliced_seq_prefers_edited_transcript_sequence() {
        let mut tx = make_tx("NM_EDIT.1", None, None, None, None);
        tx.bam_edit_status = Some("ok".to_string());
        tx.cdna_coding_start = Some(5);
        tx.cdna_coding_end = Some(13);
        tx.spliced_seq = Some("AAAATCGCACTTTGGG".to_string());

        assert_eq!(
            edited_refseq_translation_cds_from_spliced_seq(&tx).as_deref(),
            Some("TCGCACTTT")
        );
    }

    #[test]
    fn edited_refseq_translation_cds_from_spliced_seq_ignores_unedited_transcripts() {
        let mut tx = make_tx("NM_RAW.1", None, None, None, None);
        tx.cdna_coding_start = Some(5);
        tx.cdna_coding_end = Some(13);
        tx.spliced_seq = Some("AAAATCGCACTTTGGG".to_string());

        assert_eq!(edited_refseq_translation_cds_from_spliced_seq(&tx), None);
    }

    fn make_buffer_batch(chrom: &str, start: i64, end: i64) -> RecordBatch {
        let schema = Arc::new(Schema::new(vec![
            Field::new("chrom", DataType::Utf8, false),
            Field::new("start", DataType::Int64, false),
            Field::new("end", DataType::Int64, false),
        ]));
        RecordBatch::try_new(
            schema,
            vec![
                Arc::new(StringArray::from(vec![chrom])) as Arc<dyn Array>,
                Arc::new(Int64Array::from(vec![start])) as Arc<dyn Array>,
                Arc::new(Int64Array::from(vec![end])) as Arc<dyn Array>,
            ],
        )
        .unwrap()
    }

    #[test]
    fn test_transcript_selection_flags_reject_invalid_combinations() {
        let err =
            TranscriptSelectionFlags::from_options_json(Some("{\"refseq\":true,\"merged\":true}"))
                .expect_err("refseq+merged should be rejected")
                .to_string();
        assert!(err.contains("--refseq and --merged"));

        let err = TranscriptSelectionFlags::from_options_json(Some(
            "{\"refseq\":true,\"gencode_basic\":true}",
        ))
        .expect_err("refseq+gencode_basic should be rejected")
        .to_string();
        assert!(err.contains("--refseq and --gencode_basic"));

        let err = TranscriptSelectionFlags::from_options_json(Some(
            "{\"gencode_basic\":true,\"gencode_primary\":true}",
        ))
        .expect_err("gencode_basic+gencode_primary should be rejected")
        .to_string();
        assert!(err.contains("--gencode_basic and --gencode_primary"));

        let err = TranscriptSelectionFlags::from_options_json(Some("{\"all_refseq\":true}"))
            .expect_err("all_refseq without refseq/merged should be rejected")
            .to_string();
        assert!(err.contains("--all_refseq requires --refseq or --merged"));

        let err = TranscriptSelectionFlags::from_options_json(Some("{\"exclude_predicted\":true}"))
            .expect_err("exclude_predicted without refseq/merged should be rejected")
            .to_string();
        assert!(err.contains("--exclude_predicted requires --refseq or --merged"));
    }

    #[test]
    fn test_csq_placeholder_layout_matches_schema_width_for_all_modes() {
        for (everything, selection, expected_len) in [
            (false, TranscriptSelectionFlags::default(), 74usize),
            (
                false,
                TranscriptSelectionFlags {
                    source_mode: TranscriptSourceMode::RefSeq,
                    ..Default::default()
                },
                78,
            ),
            (
                false,
                TranscriptSelectionFlags {
                    source_mode: TranscriptSourceMode::Merged,
                    ..Default::default()
                },
                79,
            ),
            (true, TranscriptSelectionFlags::default(), 80),
            (
                true,
                TranscriptSelectionFlags {
                    source_mode: TranscriptSourceMode::RefSeq,
                    ..Default::default()
                },
                85,
            ),
            (
                true,
                TranscriptSelectionFlags {
                    source_mode: TranscriptSourceMode::Merged,
                    ..Default::default()
                },
                86,
            ),
        ] {
            let layout = CsqPlaceholderLayout::for_mode(everything, selection);
            assert_eq!(layout.fields.len(), expected_len);
        }
    }

    #[test]
    fn test_csq_placeholder_layout_aligns_refseq_and_merged_fields() {
        let mut frequency_fields = ColocatedFrequencyFields {
            af_values: vec![String::new(); AF_COLUMNS.len()],
            max_af: "0.9".to_string(),
            max_af_pops: "gnomADg_AFR".to_string(),
        };
        frequency_fields.af_values[0] = "0.1".to_string();
        frequency_fields.af_values[7] = "0.7".to_string();
        frequency_fields.af_values[17] = "0.17".to_string();
        let variant_fields = ColocatedVariantFields {
            existing_variation: "rs123".to_string(),
            clin_sig: "pathogenic".to_string(),
            somatic: "1".to_string(),
            pheno: "1".to_string(),
            pubmed: "12345".to_string(),
        };
        let entry = CsqPlaceholderEntry {
            allele: "G",
            consequence: "sequence_variant",
            impact: "MODIFIER",
            existing_variation: variant_fields.existing_variation.as_str(),
            variant_class: "SNV",
            frequency_fields: &frequency_fields,
            variant_fields: &variant_fields,
        };

        let refseq_selection = TranscriptSelectionFlags {
            source_mode: TranscriptSourceMode::RefSeq,
            ..Default::default()
        };
        let refseq_layout = CsqPlaceholderLayout::for_mode(false, refseq_selection);
        let mut refseq_row = String::new();
        refseq_layout.append_entry(&mut refseq_row, &entry);
        let refseq_values: Vec<&str> = refseq_row.split('|').collect();
        let refseq_fields = crate::golden_benchmark::csq_field_names_for_mode(false, true, false);
        assert_eq!(refseq_values.len(), refseq_fields.len());
        let refseq_index = |name: &str| {
            refseq_fields
                .iter()
                .position(|field| *field == name)
                .unwrap()
        };
        assert_eq!(refseq_values[refseq_index("REFSEQ_MATCH")], "");
        assert_eq!(refseq_values[refseq_index("REFSEQ_OFFSET")], "");
        assert_eq!(refseq_values[refseq_index("VARIANT_CLASS")], "SNV");
        assert_eq!(refseq_values[refseq_index("AF")], "0.1");
        assert_eq!(refseq_values[refseq_index("gnomADe_AFR")], "0.7");
        assert_eq!(refseq_values[refseq_index("gnomADg_AFR")], "0.17");
        assert_eq!(refseq_values[refseq_index("MAX_AF_POPS")], "gnomADg_AFR");

        let merged_selection = TranscriptSelectionFlags {
            source_mode: TranscriptSourceMode::Merged,
            ..Default::default()
        };
        let merged_layout = CsqPlaceholderLayout::for_mode(true, merged_selection);
        let mut merged_row = String::new();
        merged_layout.append_entry(&mut merged_row, &entry);
        let merged_values: Vec<&str> = merged_row.split('|').collect();
        let merged_fields = crate::golden_benchmark::csq_field_names_for_mode(true, false, true);
        assert_eq!(merged_values.len(), merged_fields.len());
        let merged_index = |name: &str| {
            merged_fields
                .iter()
                .position(|field| *field == name)
                .unwrap()
        };
        assert_eq!(merged_values[merged_index("SOURCE")], "");
        assert_eq!(merged_values[merged_index("REFSEQ_MATCH")], "");
        assert_eq!(merged_values[merged_index("VARIANT_CLASS")], "SNV");
        assert_eq!(merged_values[merged_index("CLIN_SIG")], "pathogenic");
        assert_eq!(merged_values[merged_index("PUBMED")], "12345");
        assert_eq!(merged_values[merged_index("MOTIF_NAME")], "");
        assert_eq!(merged_values[merged_index("TRANSCRIPTION_FACTORS")], "");
    }

    #[test]
    fn test_parse_transcript_raw_metadata_uses_direct_refseq_match_codes() {
        let raw = r#"{
          "__class":"Bio::EnsEMBL::Transcript",
          "__value":{
            "_source_cache":"RefSeq",
            "_gene_hgnc_id":"HGNC:5",
            "display_xref":{"display_id":"NM_000001"},
            "attributes":[
              {"__class":"Bio::EnsEMBL::Attribute","__value":{"code":"gencode_basic","value":"1"}},
              {"__class":"Bio::EnsEMBL::Attribute","__value":{"code":"gencode_primary","value":"1"}},
              {"__class":"Bio::EnsEMBL::Attribute","__value":{"code":"enst_refseq_compare","value":"ENST00000332831:cds_only,ENST00000619216:whole_transcript"}},
              {"__class":"Bio::EnsEMBL::Attribute","__value":{"code":"rseq_ens_match_cds","value":"1"}},
              {"__class":"Bio::EnsEMBL::Attribute","__value":{"code":"cds_start_NF","value":"1"}},
              {"__class":"Bio::EnsEMBL::Attribute","__value":{"code":"_rna_edit","value":"10 9 AAA"}},
              {"__class":"Bio::EnsEMBL::Attribute","__value":{"code":"_rna_edit","value":"20 20 G"}},
              {"__class":"Bio::EnsEMBL::Attribute","__value":{"code":"_rna_edit","value":"30 31"}},
              {"__class":"Bio::EnsEMBL::Attribute","__value":{"code":"_rna_edit","value":"40 40 T","description":"op=X"}}
            ]
          }
        }"#;

        let metadata = parse_transcript_raw_metadata(raw);
        assert_eq!(metadata.source.as_deref(), Some("RefSeq"));
        assert_eq!(metadata.display_xref_id.as_deref(), Some("NM_000001"));
        assert_eq!(metadata.gene_hgnc_id_native.as_deref(), Some("HGNC:5"));
        assert_eq!(metadata.refseq_match.as_deref(), Some("rseq_ens_match_cds"));
        assert_eq!(metadata.flags_str.as_deref(), Some("cds_start_NF"));
        assert!(metadata.is_gencode_basic);
        assert!(metadata.is_gencode_primary);
        assert_eq!(
            metadata.refseq_edits,
            vec![
                RefSeqEdit {
                    start: 10,
                    end: 9,
                    replacement_len: Some(3),
                    skip_refseq_offset: false,
                },
                RefSeqEdit {
                    start: 20,
                    end: 20,
                    replacement_len: Some(1),
                    skip_refseq_offset: true,
                },
                RefSeqEdit {
                    start: 30,
                    end: 31,
                    replacement_len: None,
                    skip_refseq_offset: false,
                },
                RefSeqEdit {
                    start: 40,
                    end: 40,
                    replacement_len: Some(1),
                    skip_refseq_offset: true,
                },
            ]
        );
    }

    #[test]
    fn test_refseq_misalignment_offset_matches_vep_rules() {
        let mut tx = make_selection_tx("NM_000001", Some("RefSeq"));
        tx.refseq_edits = vec![
            RefSeqEdit {
                start: 10,
                end: 9,
                replacement_len: Some(3),
                skip_refseq_offset: false,
            },
            RefSeqEdit {
                start: 20,
                end: 20,
                replacement_len: Some(1),
                skip_refseq_offset: true,
            },
            RefSeqEdit {
                start: 30,
                end: 31,
                replacement_len: None,
                skip_refseq_offset: false,
            },
            RefSeqEdit {
                start: 40,
                end: 40,
                replacement_len: Some(1),
                skip_refseq_offset: true,
            },
        ];

        assert_eq!(parse_cdna_position_start("35-36"), Some(35));
        assert_eq!(parse_cdna_position_start("35+2"), Some(35));
        assert_eq!(refseq_misalignment_offset(&tx, 35), Some(1));
        assert_eq!(refseq_misalignment_offset(&tx, 5), None);
        assert_eq!(refseq_misalignment_offset(&tx, 10), Some(3));

        tx.transcript_id = "NR_000001".to_string();
        assert_eq!(refseq_misalignment_offset(&tx, 35), None);
    }

    #[test]
    fn test_parse_transcript_raw_metadata_sorts_refseq_edits_by_cdna_position() {
        let raw = r#"{
          "__class":"Bio::EnsEMBL::Transcript",
          "__value":{
            "_source_cache":"BestRefSeq",
            "attributes":[
              {"__class":"Bio::EnsEMBL::Attribute","__value":{"code":"_rna_edit","value":"3723 3723 "}},
              {"__class":"Bio::EnsEMBL::Attribute","__value":{"code":"_rna_edit","value":"3228 3228 A"}},
              {"__class":"Bio::EnsEMBL::Attribute","__value":{"code":"_rna_edit","value":"1258 1258 "}}
            ]
          }
        }"#;

        let metadata = parse_transcript_raw_metadata(raw);

        assert_eq!(metadata.source.as_deref(), Some("RefSeq"));
        assert_eq!(
            metadata.refseq_edits,
            vec![
                RefSeqEdit {
                    start: 1258,
                    end: 1258,
                    replacement_len: None,
                    skip_refseq_offset: false,
                },
                RefSeqEdit {
                    start: 3228,
                    end: 3228,
                    replacement_len: Some(1),
                    skip_refseq_offset: true,
                },
                RefSeqEdit {
                    start: 3723,
                    end: 3723,
                    replacement_len: None,
                    skip_refseq_offset: false,
                },
            ]
        );
    }

    #[test]
    fn test_refseq_misalignment_offset_counts_same_coordinate_multibase_edit_as_full_insertion() {
        let mut tx = make_selection_tx("NM_001172437.2", Some("RefSeq"));
        tx.refseq_edits = vec![RefSeqEdit {
            start: 1447,
            end: 1447,
            replacement_len: Some(2),
            skip_refseq_offset: false,
        }];

        assert_eq!(refseq_misalignment_offset(&tx, 1447), None);
        assert_eq!(refseq_misalignment_offset(&tx, 1448), Some(2));
        assert_eq!(refseq_misalignment_offset(&tx, 2768), Some(2));
    }

    #[test]
    fn test_cdna_seq_has_explicit_three_prime_utr_ignores_phase_padding_only() {
        let mut tx = make_selection_tx("ENST_PHASE.1", Some("Ensembl"));
        tx.cdna_coding_end = Some(473);
        tx.cdna_seq = Some(format!("N{}", "A".repeat(473)));
        assert!(!cdna_seq_has_explicit_three_prime_utr(&tx));
    }

    #[test]
    fn test_cdna_seq_has_explicit_three_prime_utr_detects_real_utr() {
        let mut tx = make_selection_tx("ENST_UTR.1", Some("Ensembl"));
        tx.cdna_coding_end = Some(473);
        tx.cdna_seq = Some(format!("N{}TTAA", "A".repeat(473)));
        assert!(cdna_seq_has_explicit_three_prime_utr(&tx));
    }

    #[test]
    fn test_parse_transcript_raw_metadata_ignores_refseq_compare_without_direct_code() {
        let raw = r#"{
          "__class":"Bio::EnsEMBL::Transcript",
          "__value":{
            "attributes":[
              {"__class":"Bio::EnsEMBL::Attribute","__value":{"code":"enst_refseq_compare","value":"ENST00000332831:cds_only"}}
            ]
          }
        }"#;

        let metadata = parse_transcript_raw_metadata(raw);
        assert_eq!(metadata.refseq_match, None);
    }

    #[test]
    fn test_parse_transcript_raw_metadata_reads_gene_hgnc_id_fallback_key() {
        let raw = r#"{
          "__class":"Bio::EnsEMBL::Transcript",
          "__value":{
            "gene_hgnc_id":"HGNC:1100"
          }
        }"#;

        let metadata = parse_transcript_raw_metadata(raw);
        assert_eq!(metadata.gene_hgnc_id_native.as_deref(), Some("HGNC:1100"));
    }

    #[test]
    fn test_parse_transcript_raw_metadata_reads_nested_transcript_sequences() {
        let raw = r#"{
          "__class":"Bio::EnsEMBL::Transcript",
          "__value":{
            "_variation_effect_feature_cache":{
              "five_prime_utr":{
                "__class":"Bio::Seq",
                "__value":{"primary_seq":{"__class":"Bio::PrimarySeq","__value":{"seq":"aaaccc"}}}
              },
              "three_prime_utr":{
                "__class":"Bio::Seq",
                "__value":{"primary_seq":{"__class":"Bio::PrimarySeq","__value":{"seq":"gggttt"}}}
              },
              "translateable_seq":"atggcc"
            },
            "spliced_seq":{
              "__class":"Bio::Seq",
              "__value":{"primary_seq":{"__class":"Bio::PrimarySeq","__value":{"seq":"aaacccatggccgggttt"}}}
            }
          }
        }"#;

        let metadata = parse_transcript_raw_metadata(raw);
        assert_eq!(metadata.five_prime_utr_seq.as_deref(), Some("AAACCC"));
        assert_eq!(metadata.three_prime_utr_seq.as_deref(), Some("GGGTTT"));
        assert_eq!(metadata.translateable_seq.as_deref(), Some("ATGGCC"));
        assert_eq!(metadata.spliced_seq.as_deref(), Some("AAACCCATGGCCGGGTTT"));
    }

    #[test]
    fn test_parse_transcript_raw_metadata_reads_nested_cdna_mapper_segments() {
        let raw = r#"{
          "__class":"Bio::EnsEMBL::Transcript",
          "__value":{
            "_variation_effect_feature_cache":{
              "mapper":{
                "__class":"Bio::EnsEMBL::TranscriptMapper",
                "__value":{
                  "exon_coord_mapper":{
                    "__class":"Bio::EnsEMBL::Mapper",
                    "__value":{
                      "_pair_cdna":{
                        "CDNA":[
                          {
                            "from":{"__class":"Bio::EnsEMBL::Mapper::Unit","__value":{"start":1,"end":10,"id":"cdna"}},
                            "to":{"__class":"Bio::EnsEMBL::Mapper::Unit","__value":{"start":101,"end":110,"id":"genome"}},
                            "ori":1
                          },
                          {
                            "from":{"__class":"Bio::EnsEMBL::Mapper::Unit","__value":{"start":11,"end":20,"id":"cdna"}},
                            "to":{"__class":"Bio::EnsEMBL::Mapper::Unit","__value":{"start":201,"end":210,"id":"genome"}},
                            "ori":1
                          }
                        ]
                      }
                    }
                  }
                }
              }
            }
          }
        }"#;

        let metadata = parse_transcript_raw_metadata(raw);
        assert_eq!(
            metadata.cdna_mapper_segments,
            vec![
                TranscriptCdnaMapperSegment {
                    genomic_start: 101,
                    genomic_end: 110,
                    cdna_start: 1,
                    cdna_end: 10,
                    ori: 1,
                },
                TranscriptCdnaMapperSegment {
                    genomic_start: 201,
                    genomic_end: 210,
                    cdna_start: 11,
                    cdna_end: 20,
                    ori: 1,
                }
            ]
        );
    }

    #[test]
    fn test_parse_transcript_raw_metadata_preserves_gapped_cdna_mapper_segments() {
        let raw = r#"{
          "__class":"Bio::EnsEMBL::Transcript",
          "__value":{
            "_variation_effect_feature_cache":{
              "mapper":{
                "__class":"Bio::EnsEMBL::TranscriptMapper",
                "__value":{
                  "exon_coord_mapper":{
                    "__class":"Bio::EnsEMBL::Mapper",
                    "__value":{
                      "_pair_cdna":{
                        "CDNA":[
                          {
                            "from":{"__class":"Bio::EnsEMBL::Mapper::Unit","__value":{"start":1,"end":10,"id":"cdna"}},
                            "to":{"__class":"Bio::EnsEMBL::Mapper::Unit","__value":{"start":101,"end":110,"id":"genome"}},
                            "ori":1
                          },
                          {
                            "from":{"__class":"Bio::EnsEMBL::Mapper::Unit","__value":{"start":17,"end":20,"id":"cdna"}},
                            "to":{"__class":"Bio::EnsEMBL::Mapper::Unit","__value":{"start":111,"end":114,"id":"genome"}},
                            "ori":1
                          }
                        ]
                      }
                    }
                  }
                }
              }
            }
          }
        }"#;

        let metadata = parse_transcript_raw_metadata(raw);
        assert_eq!(
            metadata.cdna_mapper_segments,
            vec![
                TranscriptCdnaMapperSegment {
                    genomic_start: 101,
                    genomic_end: 110,
                    cdna_start: 1,
                    cdna_end: 10,
                    ori: 1,
                },
                TranscriptCdnaMapperSegment {
                    genomic_start: 111,
                    genomic_end: 114,
                    cdna_start: 17,
                    cdna_end: 20,
                    ori: 1,
                }
            ]
        );
    }

    #[test]
    fn test_synthesize_spliced_seq_rebuilds_full_transcript_from_utr_and_cds() {
        assert_eq!(
            synthesize_spliced_seq(
                Some("AAACCC"),
                Some("ATGGCC"),
                Some("GGGTTT"),
                Some(7),
                Some(12),
                None,
            )
            .as_deref(),
            Some("AAACCCATGGCCGGGTTT")
        );
    }

    #[test]
    fn test_passes_transcript_selection_matches_vep_refseq_filters() {
        let ensembl_tx = make_selection_tx("ENST00000311111", Some("Ensembl"));
        let nm_tx = make_selection_tx("NM_000001", Some("RefSeq"));
        let mut ccds_tx = make_selection_tx("CCDS1234.1", Some("RefSeq"));
        ccds_tx.display_xref_id = Some("CCDS1234".to_string());
        let xm_tx = make_selection_tx("XM_123456", Some("RefSeq"));
        let mut gencode_tx = make_selection_tx("ENST00000322222", Some("Ensembl"));
        gencode_tx.is_gencode_basic = true;
        gencode_tx.is_gencode_primary = true;

        let default_selection = TranscriptSelectionFlags::from_options_json(None).unwrap();
        assert!(passes_transcript_selection(&ensembl_tx, default_selection));
        assert!(!passes_transcript_selection(&nm_tx, default_selection));

        let refseq_selection =
            TranscriptSelectionFlags::from_options_json(Some("{\"refseq\":true}")).unwrap();
        assert!(!passes_transcript_selection(&ensembl_tx, refseq_selection));
        assert!(passes_transcript_selection(&nm_tx, refseq_selection));
        assert!(
            !passes_transcript_selection(&ccds_tx, refseq_selection),
            "CCDS rows should be excluded unless all_refseq is enabled"
        );

        let all_refseq_selection = TranscriptSelectionFlags::from_options_json(Some(
            "{\"merged\":true,\"all_refseq\":true}",
        ))
        .unwrap();
        assert!(passes_transcript_selection(
            &ensembl_tx,
            all_refseq_selection
        ));
        assert!(passes_transcript_selection(&nm_tx, all_refseq_selection));
        assert!(passes_transcript_selection(&ccds_tx, all_refseq_selection));

        let exclude_predicted_selection = TranscriptSelectionFlags::from_options_json(Some(
            "{\"merged\":true,\"all_refseq\":true,\"exclude_predicted\":true}",
        ))
        .unwrap();
        assert!(
            !passes_transcript_selection(&xm_tx, exclude_predicted_selection),
            "XM_/XR_ transcripts should be filtered by exclude_predicted"
        );

        let gencode_basic_selection =
            TranscriptSelectionFlags::from_options_json(Some("{\"gencode_basic\":true}")).unwrap();
        assert!(passes_transcript_selection(
            &gencode_tx,
            gencode_basic_selection
        ));
        assert!(!passes_transcript_selection(
            &ensembl_tx,
            gencode_basic_selection
        ));

        let gencode_primary_selection = TranscriptSelectionFlags::from_options_json(Some(
            "{\"merged\":true,\"gencode_primary\":true}",
        ))
        .unwrap();
        assert!(passes_transcript_selection(
            &gencode_tx,
            gencode_primary_selection
        ));
        assert!(!passes_transcript_selection(
            &nm_tx,
            gencode_primary_selection
        ));
    }

    #[test]
    fn test_buffer_local_hgnc_propagation_uses_native_symbol_donor() {
        let tx_donor = make_tx(
            "ENST00000919191",
            Some("ENSG00000182158"),
            Some("NBAS"),
            Some("HGNC"),
            Some("HGNC:15625"),
        );
        let mut tx_refseq = make_tx(
            "XR_007076390.1",
            Some("GENE:NBAS"),
            Some("NBAS"),
            Some("EntrezGene"),
            None,
        );
        tx_refseq.source = Some("RefSeq".to_string());

        let mut transcripts = vec![tx_donor, tx_refseq];
        apply_buffer_local_hgnc_propagation(&mut transcripts);
        assert_eq!(transcripts[1].gene_hgnc_id.as_deref(), Some("HGNC:15625"));
    }

    #[test]
    fn test_buffer_local_hgnc_propagation_ignores_non_native_effective_values() {
        let mut tx_promoted = make_tx(
            "ENST00000426186",
            Some("ENSG00000225475"),
            Some("ANAPC1P1"),
            Some("HGNC"),
            None,
        );
        tx_promoted.gene_hgnc_id = Some("HGNC:44150".to_string());
        let mut tx_refseq = make_tx(
            "NR_037931.2",
            Some("GENE:ANAPC1P1"),
            Some("ANAPC1P1"),
            Some("EntrezGene"),
            None,
        );
        tx_refseq.source = Some("RefSeq".to_string());

        let mut transcripts = vec![tx_promoted, tx_refseq];
        apply_buffer_local_hgnc_propagation(&mut transcripts);
        assert_eq!(
            transcripts[1].gene_hgnc_id, None,
            "cache-promoted HGNC IDs must not seed VEP-style propagation"
        );
    }

    #[test]
    fn test_buffer_local_hgnc_propagation_refills_same_gene_stable_id() {
        let tx_with_symbol = make_tx(
            "ENST00000111111",
            Some("ENSG00000123456"),
            Some("BRCA1"),
            Some("HGNC"),
            Some("HGNC:1100"),
        );
        let tx_missing = make_tx("ENST00000222222", Some("ENSG00000123456"), None, None, None);

        let mut transcripts = vec![tx_with_symbol, tx_missing];
        apply_buffer_local_hgnc_propagation(&mut transcripts);
        assert_eq!(transcripts[1].gene_symbol.as_deref(), Some("BRCA1"));
        assert_eq!(transcripts[1].gene_symbol_source.as_deref(), Some("HGNC"));
        assert_eq!(transcripts[1].gene_hgnc_id.as_deref(), Some("HGNC:1100"));
    }

    #[test]
    fn test_build_buffer_local_transcripts_scopes_to_expanded_range() {
        let mut tx_before = make_tx(
            "ENST00000000001",
            Some("ENSG_BEFORE"),
            Some("GENE1"),
            Some("HGNC"),
            Some("HGNC:1"),
        );
        tx_before.start = 10;
        tx_before.end = 20;
        let mut tx_inside = make_tx(
            "ENST00000000002",
            Some("ENSG_INSIDE"),
            Some("GENE2"),
            Some("HGNC"),
            Some("HGNC:2"),
        );
        tx_inside.start = 120;
        tx_inside.end = 180;
        let mut tx_after = make_tx(
            "ENST00000000003",
            Some("ENSG_AFTER"),
            Some("GENE3"),
            Some("HGNC"),
            Some("HGNC:3"),
        );
        tx_after.start = 400;
        tx_after.end = 450;

        let scoped = build_buffer_local_transcripts(
            &[tx_before, tx_inside.clone(), tx_after],
            "chr2",
            140,
            160,
            50,
            50,
        );

        assert_eq!(scoped.len(), 1);
        assert_eq!(scoped[0].transcript_id, tx_inside.transcript_id);
    }

    #[test]
    fn test_stateful_buffer_local_transcripts_carry_hgnc_across_adjacent_buffers() {
        let mut tx_donor = make_tx(
            "ENST_DONOR",
            Some("ENSG00000123456"),
            Some("LINC02778"),
            Some("HGNC"),
            Some("HGNC:54298"),
        );
        tx_donor.start = 100;
        tx_donor.end = 200;

        let mut tx_recipient = make_tx(
            "XR_RECIPIENT",
            Some("105378760"),
            Some("LINC02778"),
            Some("EntrezGene"),
            None,
        );
        tx_recipient.start = 100;
        tx_recipient.end = 700_000;

        let transcripts = vec![tx_donor, tx_recipient];
        let transcript_regions: HashMap<String, Vec<TranscriptCacheRegion>> = transcripts
            .iter()
            .map(|tx| (tx.transcript_id.clone(), transcript_cache_regions(tx)))
            .collect();
        let mut persisted_transcripts = HashMap::new();

        let first_buffer = vec![make_buffer_batch("chr2", 150, 150)];
        let first_scoped = build_stateful_buffer_local_transcripts(
            &transcripts,
            &transcript_regions,
            &mut persisted_transcripts,
            &first_buffer,
            "chr2",
            150,
            150,
            50,
            50,
        )
        .unwrap();
        let first_recipient = first_scoped
            .iter()
            .find(|tx| tx.transcript_id == "XR_RECIPIENT")
            .unwrap();
        assert_eq!(first_recipient.gene_hgnc_id.as_deref(), Some("HGNC:54298"));

        let second_buffer = vec![make_buffer_batch("chr2", 650_000, 650_000)];
        let second_scoped = build_stateful_buffer_local_transcripts(
            &transcripts,
            &transcript_regions,
            &mut persisted_transcripts,
            &second_buffer,
            "chr2",
            650_000,
            650_000,
            50,
            50,
        )
        .unwrap();
        assert_eq!(second_scoped.len(), 1);
        assert_eq!(second_scoped[0].transcript_id, "XR_RECIPIENT");
        assert_eq!(second_scoped[0].gene_hgnc_id.as_deref(), Some("HGNC:54298"));
    }

    #[test]
    fn test_stateful_buffer_local_transcripts_keep_hgnc_within_start_cache_region() {
        let mut tx_donor = make_tx(
            "ENST_DONOR",
            Some("ENSG00000175463"),
            Some("SUGCT"),
            Some("HGNC"),
            Some("HGNC:16001"),
        );
        tx_donor.chrom = "chr7".to_string();
        tx_donor.start = 40_134_044;
        tx_donor.end = 40_861_613;

        let mut tx_recipient = make_tx(
            "XR_007060157.1",
            Some("79690"),
            Some("SUGCT"),
            Some("EntrezGene"),
            None,
        );
        tx_recipient.chrom = "chr7".to_string();
        tx_recipient.start = 40_135_005;
        tx_recipient.end = 41_038_816;

        let transcripts = vec![tx_donor, tx_recipient];
        let transcript_regions: HashMap<String, Vec<TranscriptCacheRegion>> = transcripts
            .iter()
            .map(|tx| (tx.transcript_id.clone(), transcript_cache_regions(tx)))
            .collect();
        let mut persisted_transcripts = HashMap::new();

        let first_buffer = vec![make_buffer_batch("chr7", 40_500_000, 40_500_000)];
        let first_scoped = build_stateful_buffer_local_transcripts(
            &transcripts,
            &transcript_regions,
            &mut persisted_transcripts,
            &first_buffer,
            "chr7",
            40_500_000,
            40_500_000,
            50,
            50,
        )
        .unwrap();
        let first_recipient = first_scoped
            .iter()
            .find(|tx| tx.transcript_id == "XR_007060157.1")
            .unwrap();
        assert_eq!(first_recipient.gene_hgnc_id.as_deref(), Some("HGNC:16001"));

        let second_buffer = vec![make_buffer_batch("chr7", 40_986_831, 40_986_831)];
        let second_scoped = build_stateful_buffer_local_transcripts(
            &transcripts,
            &transcript_regions,
            &mut persisted_transcripts,
            &second_buffer,
            "chr7",
            40_986_831,
            40_986_831,
            50,
            50,
        )
        .unwrap();
        assert_eq!(second_scoped.len(), 1);
        assert_eq!(second_scoped[0].transcript_id, "XR_007060157.1");
        assert_eq!(second_scoped[0].gene_hgnc_id.as_deref(), Some("HGNC:16001"));
    }

    #[test]
    fn test_stateful_buffer_local_transcripts_prune_hgnc_across_cache_regions() {
        let mut tx_donor = make_tx(
            "ENST_DONOR",
            Some("ENSG00000123456"),
            Some("PDK1"),
            Some("HGNC"),
            Some("HGNC:8809"),
        );
        tx_donor.start = 100;
        tx_donor.end = 200;

        let mut tx_region0 = make_tx(
            "XR_REGION0",
            Some("5163"),
            Some("PDK1"),
            Some("EntrezGene"),
            None,
        );
        tx_region0.start = 100;
        tx_region0.end = 700_000;

        let mut tx_region1 = make_tx(
            "XR_REGION1",
            Some("5163"),
            Some("PDK1"),
            Some("EntrezGene"),
            None,
        );
        tx_region1.start = 1_050_000;
        tx_region1.end = 1_060_000;

        let transcripts = vec![tx_donor, tx_region0, tx_region1];
        let transcript_regions: HashMap<String, Vec<TranscriptCacheRegion>> = transcripts
            .iter()
            .map(|tx| (tx.transcript_id.clone(), transcript_cache_regions(tx)))
            .collect();
        let mut persisted_transcripts = HashMap::new();

        let first_buffer = vec![make_buffer_batch("chr2", 150, 150)];
        build_stateful_buffer_local_transcripts(
            &transcripts,
            &transcript_regions,
            &mut persisted_transcripts,
            &first_buffer,
            "chr2",
            150,
            150,
            50,
            50,
        )
        .unwrap();

        let second_buffer = vec![make_buffer_batch("chr2", 1_050_000, 1_050_000)];
        let second_scoped = build_stateful_buffer_local_transcripts(
            &transcripts,
            &transcript_regions,
            &mut persisted_transcripts,
            &second_buffer,
            "chr2",
            1_050_000,
            1_050_000,
            50,
            50,
        )
        .unwrap();
        assert_eq!(second_scoped.len(), 1);
        assert_eq!(second_scoped[0].transcript_id, "XR_REGION1");
        assert_eq!(second_scoped[0].gene_hgnc_id, None);
    }

    #[test]
    fn test_stateful_buffer_local_transcripts_do_not_carry_same_transcript_across_regions() {
        let mut tx_donor = make_tx(
            "ENST_DONOR",
            Some("ENSG00000181143"),
            Some("MUC16"),
            Some("HGNC"),
            Some("HGNC:15582"),
        );
        tx_donor.chrom = "chr19".to_string();
        tx_donor.start = 8_848_844;
        tx_donor.end = 9_010_390;

        let mut tx_recipient = make_tx(
            "NM_001414686.1",
            Some("94025"),
            Some("MUC16"),
            Some("EntrezGene"),
            None,
        );
        tx_recipient.chrom = "chr19".to_string();
        tx_recipient.start = 8_848_844;
        tx_recipient.end = 9_065_751;

        let transcripts = vec![tx_donor, tx_recipient];
        let transcript_regions: HashMap<String, Vec<TranscriptCacheRegion>> = transcripts
            .iter()
            .map(|tx| (tx.transcript_id.clone(), transcript_cache_regions(tx)))
            .collect();
        let mut persisted_transcripts = HashMap::new();

        let first_buffer = vec![make_buffer_batch("chr19", 8_900_000, 8_900_000)];
        let first_scoped = build_stateful_buffer_local_transcripts(
            &transcripts,
            &transcript_regions,
            &mut persisted_transcripts,
            &first_buffer,
            "chr19",
            8_900_000,
            8_900_000,
            50,
            50,
        )
        .unwrap();
        let first_recipient = first_scoped
            .iter()
            .find(|tx| tx.transcript_id == "NM_001414686.1")
            .unwrap();
        assert_eq!(first_recipient.gene_hgnc_id.as_deref(), Some("HGNC:15582"));

        let second_buffer = vec![make_buffer_batch("chr19", 9_058_432, 9_058_432)];
        let second_scoped = build_stateful_buffer_local_transcripts(
            &transcripts,
            &transcript_regions,
            &mut persisted_transcripts,
            &second_buffer,
            "chr19",
            9_058_432,
            9_058_432,
            50,
            50,
        )
        .unwrap();
        assert_eq!(second_scoped.len(), 1);
        assert_eq!(second_scoped[0].transcript_id, "NM_001414686.1");
        assert_eq!(second_scoped[0].gene_hgnc_id, None);
    }

    #[test]
    fn test_stateful_buffer_local_transcripts_clears_promoted_hgnc_before_later_muc16_buffer() {
        let mut tx_donor = make_tx(
            "ENST_DONOR",
            Some("ENSG00000181143"),
            Some("MUC16"),
            Some("HGNC"),
            Some("HGNC:15582"),
        );
        tx_donor.chrom = "chr19".to_string();
        tx_donor.start = 8_848_844;
        tx_donor.end = 9_010_390;

        let mut tx_recipient = make_tx(
            "NM_001414686.1",
            Some("94025"),
            Some("MUC16"),
            Some("EntrezGene"),
            None,
        );
        tx_recipient.chrom = "chr19".to_string();
        tx_recipient.start = 8_848_844;
        tx_recipient.end = 9_065_751;

        let transcripts = vec![tx_donor, tx_recipient];
        let transcript_regions: HashMap<String, Vec<TranscriptCacheRegion>> = transcripts
            .iter()
            .map(|tx| (tx.transcript_id.clone(), transcript_cache_regions(tx)))
            .collect();
        let mut persisted_transcripts = HashMap::new();

        let first_buffer = vec![
            make_buffer_batch("chr19", 8_900_000, 8_900_000),
            make_buffer_batch("chr19", 9_058_364, 9_058_364),
        ];
        let first_scoped = build_stateful_buffer_local_transcripts(
            &transcripts,
            &transcript_regions,
            &mut persisted_transcripts,
            &first_buffer,
            "chr19",
            8_900_000,
            9_058_364,
            50,
            50,
        )
        .unwrap();
        let first_recipient = first_scoped
            .iter()
            .find(|tx| tx.transcript_id == "NM_001414686.1")
            .unwrap();
        assert_eq!(first_recipient.gene_hgnc_id.as_deref(), Some("HGNC:15582"));

        let second_buffer = vec![make_buffer_batch("chr19", 9_058_432, 9_058_432)];
        let second_scoped = build_stateful_buffer_local_transcripts(
            &transcripts,
            &transcript_regions,
            &mut persisted_transcripts,
            &second_buffer,
            "chr19",
            9_058_432,
            9_058_432,
            50,
            50,
        )
        .unwrap();
        assert_eq!(second_scoped.len(), 1);
        assert_eq!(second_scoped[0].transcript_id, "NM_001414686.1");
        assert_eq!(second_scoped[0].gene_hgnc_id, None);
    }

    // ── csq_escape ────────────────────────────────────────────────────

    #[test]
    fn test_csq_escape_comma_becomes_ampersand() {
        // Multi-accession SWISSPROT values must not split CSQ entries (issue #93)
        assert_eq!(
            csq_escape("A0A0J9YXY3.52,P0DPF7.28"),
            "A0A0J9YXY3.52&P0DPF7.28"
        );
    }

    #[test]
    fn test_csq_escape_pipe_becomes_ampersand() {
        assert_eq!(csq_escape("a|b"), "a&b");
    }

    #[test]
    fn test_csq_escape_semicolon_percent_encoded() {
        assert_eq!(csq_escape("a;b"), "a%3Bb");
    }

    #[test]
    fn test_csq_escape_dash_becomes_empty() {
        assert_eq!(csq_escape("-"), "");
    }

    #[test]
    fn test_csq_escape_no_change() {
        let val = "Q9Y6K1.3";
        let escaped = csq_escape(val);
        assert_eq!(escaped, val);
        // Should be a borrowed Cow (no allocation)
        assert!(matches!(escaped, std::borrow::Cow::Borrowed(_)));
    }
}
