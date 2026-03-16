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
use std::io::{BufRead, Seek};
use std::sync::{Arc, Mutex};
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
use noodles_core::{Position, Region};
use noodles_fasta as fasta;
use serde_json::Value;

use crate::allele::{
    MatchedVariantAllele, vcf_to_vep_allele, vcf_to_vep_input_allele, vep_norm_end, vep_norm_start,
};
use crate::annotation_store::{AnnotationBackend, build_store};
#[cfg(feature = "kv-cache")]
use crate::kv_cache::KvCacheTableProvider;
use crate::lookup_provider::LookupProvider;
use crate::miss_worklist::{MissWorklist, collect_miss_worklist};
use crate::so_terms::{SoImpact, SoTerm, most_severe_term};
use crate::transcript_consequence::{
    CachedPredictions, CompactPrediction, ExonFeature, MirnaFeature, MotifFeature, PreparedContext,
    ProteinDomainFeature, RegulatoryFeature, SiftPolyphenCache, StructuralFeature, SvEventKind,
    SvFeatureKind, TranscriptCdnaMapperSegment, TranscriptConsequenceEngine, TranscriptFeature,
    TranslationFeature, VariantInput, is_vep_transcript,
};
use crate::variant_lookup_exec::{
    ColocatedCacheEntry, ColocatedKey, ColocatedSink, ColocatedSinkValue,
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
/// Output is sorted, `&`-joined unique SO terms.
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

    // Map variant cDNA positions to structure indices and collect SO terms.
    let mut has_stem = false;
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
            b'(' | b')' => has_stem = true,
            b'.' => has_loop = true,
            _ => {}
        }
    }

    // Sorted output: miRNA_loop < miRNA_stem (alphabetical).
    match (has_loop, has_stem) {
        (true, true) => "miRNA_loop&miRNA_stem".to_string(),
        (true, false) => "miRNA_loop".to_string(),
        (false, true) => "miRNA_stem".to_string(),
        (false, false) => String::new(),
    }
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
    cache: &SiftPolyphenCache,
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
                "gene_hgnc_id",
                "refseq_id",
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
        let mut refseq_ids: Vec<Option<String>> = Vec::new();

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
            let gene_hgnc_id_idx = schema.index_of("gene_hgnc_id").ok();
            let refseq_id_idx = schema.index_of("refseq_id").ok();
            let source_idx = schema.index_of("source").ok();
            let version_idx = schema.index_of("version").ok();
            let raw_json_idx = schema.index_of("raw_object_json").ok();
            let cds_start_nf_idx = schema.index_of("cds_start_nf").ok();
            let cds_end_nf_idx = schema.index_of("cds_end_nf").ok();
            let mirna_regions_idx = schema.index_of("mature_mirna_regions").ok();
            let cdna_seq_idx = schema.index_of("cdna_seq").ok();
            // Promoted columns (previously extracted from raw_object_json).
            let bam_edit_status_idx = schema.index_of("bam_edit_status").ok();
            let has_non_polya_rna_edit_idx = schema.index_of("has_non_polya_rna_edit").ok();
            let spliced_seq_idx = schema.index_of("spliced_seq").ok();
            let translateable_seq_idx = schema.index_of("translateable_seq").ok();
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
                // Read from promoted top-level parquet columns
                // (biodatageeks/datafusion-bio-formats#125, #126).
                if let Some(translateable_seq) =
                    translateable_seq_idx.and_then(|idx| string_at(batch.column(idx).as_ref(), row))
                {
                    translateable_seq_by_tx.insert(transcript_id.clone(), translateable_seq);
                }
                let cdna_mapper_segments = cdna_mapper_segments_idx
                    .map(|idx| {
                        cdna_mapper_segments_from_list_column(batch.column(idx).as_ref(), row)
                    })
                    .unwrap_or_default();
                let flags_str = flags_str_idx
                    .and_then(|idx| string_at(batch.column(idx).as_ref(), row))
                    .or_else(|| flags_str_from_bools(cds_start_nf, cds_end_nf));

                let gene_stable_id =
                    gene_stable_id_idx.and_then(|idx| string_at(batch.column(idx).as_ref(), row));
                let gene_symbol =
                    gene_symbol_idx.and_then(|idx| string_at(batch.column(idx).as_ref(), row));
                let gene_symbol_source = gene_symbol_source_idx
                    .and_then(|idx| string_at(batch.column(idx).as_ref(), row));
                let gene_hgnc_id =
                    gene_hgnc_id_idx.and_then(|idx| string_at(batch.column(idx).as_ref(), row));
                let refseq_id =
                    refseq_id_idx.and_then(|idx| string_at(batch.column(idx).as_ref(), row));
                let source = source_idx.and_then(|idx| string_at(batch.column(idx).as_ref(), row));
                let bam_edit_status =
                    bam_edit_status_idx.and_then(|idx| string_at(batch.column(idx).as_ref(), row));
                let has_non_polya_rna_edit = has_non_polya_rna_edit_idx
                    .and_then(|idx| bool_at(batch.column(idx).as_ref(), row))
                    .unwrap_or(false);
                let spliced_seq =
                    spliced_seq_idx.and_then(|idx| string_at(batch.column(idx).as_ref(), row));
                let cdna_seq =
                    cdna_seq_idx.and_then(|idx| string_at(batch.column(idx).as_ref(), row));
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
                    gene_hgnc_id,
                    source,
                    bam_edit_status,
                    has_non_polya_rna_edit,
                    spliced_seq,
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
                refseq_ids.push(refseq_id);
            }
        }

        backfill_missing_hgnc_ids(&mut out, &refseq_ids);
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
            let cds_seq_idx = schema
                .index_of("cds_sequence")
                .or_else(|_| schema.index_of("cds_seq"))
                .or_else(|_| schema.index_of("coding_sequence"))
                .ok();
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

    async fn scan_with_transcript_engine(
        &self,
        state: &dyn Session,
        projection: Option<&Vec<usize>>,
        requested_columns: &[&str],
        extended_probes: bool,
        cache_table: &str,
        transcripts_table: Option<&str>,
        exons_table: Option<&str>,
        translations_table: Option<&str>,
        regulatory_table: Option<&str>,
        motif_table: Option<&str>,
        mirna_table: Option<&str>,
        sv_table: Option<&str>,
    ) -> Result<Arc<dyn ExecutionPlan>> {
        let t_total = profile_start!();
        if profiling_enabled() {
            eprintln!("[VEP_PROFILE] ====== scan_with_transcript_engine START ======");
        }
        let flags = VepFlags::from_options_json(self.options_json.as_deref());
        let hgvs_flags = HgvsFlags::from_options_json(self.options_json.as_deref());

        // Build the lookup plan directly (bypassing SQL) so we can attach
        // a co-located data sink that piggybacks on the same cache scan.
        let coloc_sink: ColocatedSink = Arc::new(Mutex::new(HashMap::new()));

        let vcf_schema = self
            .session
            .table(&self.vcf_table)
            .await?
            .schema()
            .as_arrow()
            .clone();
        let cache_schema = self
            .session
            .table(cache_table)
            .await?
            .schema()
            .as_arrow()
            .clone();
        let cache_columns: Vec<String> = requested_columns.iter().map(|s| s.to_string()).collect();
        let allowed_failed = self
            .options_json
            .as_deref()
            .and_then(|opts| Self::parse_json_i64_option(opts, "failed"))
            .unwrap_or(0);
        let reference_fasta_path = self
            .options_json
            .as_deref()
            .and_then(|opts| Self::parse_json_string_option(opts, "reference_fasta_path"));
        Self::validate_hgvs_reference_fasta(hgvs_flags, reference_fasta_path.as_deref())?;
        let mut hgvs_reference_reader = if hgvs_flags.any() && hgvs_flags.shift_hgvs {
            reference_fasta_path
                .as_deref()
                .map(|path| {
                    fasta::io::indexed_reader::Builder::default()
                        .build_from_path(path)
                        .map_err(|e| {
                            DataFusionError::Execution(format!(
                                "failed to open indexed reference FASTA '{path}': {e}"
                            ))
                        })
                })
                .transpose()?
        } else {
            None
        };
        let mut provider = LookupProvider::new(
            Arc::clone(&self.session),
            self.vcf_table.clone(),
            cache_table.to_string(),
            vcf_schema,
            cache_schema,
            cache_columns,
            extended_probes,
            allowed_failed,
            None,
        )?;
        if flags.check_existing {
            provider.set_colocated_sink(Arc::clone(&coloc_sink));
        }

        let t_var = profile_start!();
        let plan = provider.scan(state, None, &[], None).await?;
        let task_ctx = self.session.task_ctx();
        let base_batches = datafusion::physical_plan::collect(plan, task_ctx).await?;
        let total_vcf_rows: usize = base_batches.iter().map(|b| b.num_rows()).sum();
        profile_end!(
            "1. variation_lookup (scan+collect)",
            t_var,
            format!(
                "{} VCF rows, {} batches",
                total_vcf_rows,
                base_batches.len()
            )
        );

        let t_intervals = profile_start!();
        let input_variant_intervals = collect_input_variant_intervals(&base_batches)?;
        let indel_variant_intervals = collect_indel_variant_intervals(&base_batches)?;
        profile_end!("2. collect_variant_intervals", t_intervals);

        // Build co-located variant aggregation from the piggybacked sink.
        let t_coloc = profile_start!();
        let colocated_map = if flags.check_existing {
            let guard = coloc_sink.lock().unwrap();
            build_colocated_map_from_sink(&guard)
        } else {
            HashMap::new()
        };
        profile_end!(
            "3. colocated_map_build",
            t_coloc,
            format!("{} entries", colocated_map.len())
        );

        let mut cache_hit_count = 0usize;
        let mut cache_miss_count = 0usize;
        let worklist =
            collect_miss_worklist(&base_batches, &mut cache_hit_count, &mut cache_miss_count)?;
        let has_cache_misses = !worklist.is_empty();
        if profiling_enabled() {
            let hit_rate = if total_vcf_rows > 0 {
                cache_hit_count as f64 / total_vcf_rows as f64 * 100.0
            } else {
                0.0
            };
            eprintln!(
                "[VEP_PROFILE] cache hits: {}, misses: {}, hit rate: {:.1}%",
                cache_hit_count, cache_miss_count, hit_rate
            );
        }

        let t_ctx_load = profile_start!();
        let (
            mut transcripts,
            translateable_seq_by_tx,
            exons,
            mut translations,
            regulatory,
            motifs,
            mirnas,
            structural,
        ) = if has_cache_misses {
            let merged = self
                .options_json
                .as_deref()
                .and_then(|opts| Self::parse_json_bool_option(opts, "merged"))
                .unwrap_or(false);

            let t_tx = profile_start!();
            let (tx, translateable_seq_by_tx) = if let Some(table) = transcripts_table {
                let (tx, translateable_seq_by_tx) = self.load_transcripts(table, &worklist).await?;
                (
                    tx.into_iter()
                        .filter(|t| is_vep_transcript(&t.transcript_id, merged))
                        .collect::<Vec<_>>(),
                    translateable_seq_by_tx,
                )
            } else {
                (Vec::new(), HashMap::new())
            };
            profile_end!(
                "4a. load_transcripts",
                t_tx,
                format!("{} transcripts", tx.len())
            );

            let tx_ids: HashSet<String> = tx.iter().map(|t| t.transcript_id.clone()).collect();

            let t_ex = profile_start!();
            let ex = if let Some(table) = exons_table {
                let raw = self.load_exons(table, &worklist).await?;
                raw.into_iter()
                    .filter(|e| tx_ids.contains(&e.transcript_id))
                    .collect()
            } else {
                Vec::new()
            };
            profile_end!("4b. load_exons", t_ex, format!("{} exons", ex.len()));

            let t_tl = profile_start!();
            let tl = if let Some(table) = translations_table {
                let raw = self.load_translations(table, &worklist).await?;
                raw.into_iter()
                    .filter(|t| tx_ids.contains(&t.transcript_id))
                    .collect()
            } else {
                Vec::new()
            };
            profile_end!(
                "4c. load_translations",
                t_tl,
                format!("{} translations", tl.len())
            );

            let t_rg = profile_start!();
            let rg = if let Some(table) = regulatory_table {
                self.load_regulatory_features(table, &worklist).await?
            } else {
                Vec::new()
            };
            profile_end!(
                "4d. load_regulatory",
                t_rg,
                format!("{} features", rg.len())
            );

            let t_mt = profile_start!();
            let mt = if let Some(table) = motif_table {
                self.load_motif_features(table, &worklist).await?
            } else {
                Vec::new()
            };
            profile_end!("4e. load_motif", t_mt);

            let t_mi = profile_start!();
            let mi = if let Some(table) = mirna_table {
                self.load_mirna_features(table, &worklist).await?
            } else {
                Vec::new()
            };
            profile_end!("4f. load_mirna", t_mi);

            let t_st = profile_start!();
            let st = if let Some(table) = sv_table {
                self.load_structural_features(table, &worklist).await?
            } else {
                Vec::new()
            };
            profile_end!("4g. load_structural", t_st);
            (tx, translateable_seq_by_tx, ex, tl, rg, mt, mi, st)
        } else {
            (
                Vec::new(),
                HashMap::new(),
                Vec::new(),
                Vec::new(),
                Vec::new(),
                Vec::new(),
                Vec::new(),
                Vec::new(),
            )
        };
        profile_end!("4. context_tables_total", t_ctx_load);

        let t_hydrate = profile_start!();
        let hydrated_transcript_ids = if let Some(reader) = hgvs_reference_reader.as_mut() {
            hydrate_refseq_translation_cds_from_reference(
                reader,
                &transcripts,
                &exons,
                &mut translations,
                &input_variant_intervals,
            )?
        } else {
            HashSet::new()
        };
        apply_translateable_seq_overrides(
            &mut translations,
            &translateable_seq_by_tx,
            &hydrated_transcript_ids,
        );
        profile_end!(
            "5a. hydrate_refseq_cds",
            t_hydrate,
            format!("{} hydrated", hydrated_transcript_ids.len())
        );

        let t_cdna = profile_start!();
        // Hydrate full spliced cDNA on transcripts that lack spliced_seq so
        // that three_prime_utr_seq() can find the UTR for stop-loss/frameshift
        // extension distance computation.
        if let Some(reader) = hgvs_reference_reader.as_mut() {
            hydrate_transcript_cdna_from_reference(
                reader,
                &mut transcripts,
                &exons,
                &indel_variant_intervals,
                &input_variant_intervals,
            )?;
        }
        profile_end!("5b. hydrate_transcript_cdna", t_cdna);

        let t_prep = profile_start!();
        let (upstream_distance, downstream_distance) = self.transcript_distance_config();
        let engine = TranscriptConsequenceEngine::new_with_hgvs_shift(
            upstream_distance,
            downstream_distance,
            hgvs_flags.shift_hgvs,
        );
        let ctx = PreparedContext::new(
            &transcripts,
            &exons,
            &translations,
            &regulatory,
            &motifs,
            &mirnas,
            &structural,
        );
        profile_end!(
            "6. prepared_context_build",
            t_prep,
            format!("{} tx_trees chroms", ctx.tx_trees.len())
        );

        // Build SIFT/PolyPhen prediction cache using lazy sliding windows.
        //
        // Instead of loading all translations upfront, windows are loaded
        // lazily as the annotation batch loop advances through genomic
        // positions. SQL queries only fire when a batch crosses a 5MB window
        // boundary. Entries whose genomic end falls behind the current batch
        // are evicted (safe because VCF is position-sorted).
        //
        // See biodatageeks/datafusion-bio-functions#38.
        let t_sift = profile_start!();
        let mut sift_cache = SiftPolyphenCache::new();
        let sift_table_name: Option<String> = self
            .options_json
            .as_deref()
            .and_then(|opts| Self::parse_json_string_option(opts, "translations_sift_table"))
            .or_else(|| translations_table.map(|s| s.to_string()));
        let sift_enabled = flags.everything && sift_table_name.is_some();
        let mut loaded_windows: HashSet<(String, i64)> = HashSet::new();
        // Try to build a direct parquet reader for sift, bypassing DataFusion SQL.
        let sift_direct: Option<SiftDirectReader> = if sift_enabled {
            let table_name = sift_table_name.as_deref().unwrap();
            // Resolve file path: try cache_source parent dir + table_name.parquet
            let cache_dir = std::path::Path::new(&self.cache_source)
                .parent()
                .map(|p| p.to_path_buf());
            cache_dir.and_then(|dir| {
                let path = dir.join(format!("{table_name}.parquet"));
                if path.exists() {
                    Self::build_sift_direct_reader(path.to_str()?)
                } else {
                    None
                }
            })
        } else {
            None
        };
        if profiling_enabled() && sift_direct.is_some() {
            eprintln!("[VEP_PROFILE] sift_direct_reader: enabled (bypassing DataFusion SQL)");
        }
        profile_end!("7. sift_polyphen_cache_init", t_sift);

        let t_annotate = profile_start!();
        let mut sift_load_ms = 0u128;
        let mut annotate_ms = 0u128;
        let mut annotated_batches = Vec::with_capacity(base_batches.len());
        for batch in &base_batches {
            // Lazily load SIFT/PolyPhen windows as the batch loop advances.
            if sift_enabled {
                let t_sift_win = std::time::Instant::now();
                let table = sift_table_name.as_deref().unwrap();
                let schema = batch.schema();
                if let (Ok(ci), Ok(si), Ok(ei)) = (
                    schema.index_of("chrom"),
                    schema.index_of("start"),
                    schema.index_of("end"),
                ) {
                    // Collect per-chrom min/max positions in this batch.
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

                    for (chrom, (batch_min, batch_max)) in &batch_chrom_bounds {
                        let window_start =
                            (batch_max / Self::SIFT_WINDOW_SIZE) * Self::SIFT_WINDOW_SIZE;

                        // Also load the window containing the minimum position
                        // (may differ from the max window).
                        let min_window_start =
                            (batch_min / Self::SIFT_WINDOW_SIZE) * Self::SIFT_WINDOW_SIZE;

                        // Load all windows from the min to the max + next window.
                        let mut ws = min_window_start;
                        while ws <= window_start + Self::SIFT_WINDOW_SIZE {
                            let key = (chrom.clone(), ws);
                            if !loaded_windows.contains(&key) {
                                if let Some(ref direct) = sift_direct {
                                    direct.load_window(
                                        chrom,
                                        ws,
                                        ws + Self::SIFT_WINDOW_SIZE,
                                        &mut sift_cache,
                                    )?;
                                } else {
                                    self.load_sift_window(
                                        table,
                                        chrom,
                                        ws,
                                        ws + Self::SIFT_WINDOW_SIZE,
                                        &mut sift_cache,
                                    )
                                    .await?;
                                }
                                loaded_windows.insert(key);
                            }
                            ws += Self::SIFT_WINDOW_SIZE;
                        }

                        // Evict entries whose genomic end is behind the batch
                        // minimum (position-sorted VCF guarantees no future
                        // batch will need them).
                        sift_cache.evict_before(*batch_min);
                    }
                }
                sift_load_ms += t_sift_win.elapsed().as_millis();
            }

            let t_ann = std::time::Instant::now();
            annotated_batches.push(self.annotate_batch_with_transcript_engine(
                batch,
                &engine,
                &ctx,
                &colocated_map,
                &sift_cache,
            )?);
            annotate_ms += t_ann.elapsed().as_millis();
        }
        if profiling_enabled() {
            eprintln!(
                "[VEP_PROFILE] 7a. sift_lazy_load_only...........................  {sift_load_ms}ms"
            );
            eprintln!(
                "[VEP_PROFILE] 7b. annotate_batches_only.........................  {annotate_ms}ms"
            );
        }
        profile_end!(
            "7+8. sift_lazy_load + annotate_batches",
            t_annotate,
            format!(
                "{} batches, {} total rows, {} sift windows loaded",
                annotated_batches.len(),
                annotated_batches
                    .iter()
                    .map(|b| b.num_rows())
                    .sum::<usize>(),
                loaded_windows.len()
            )
        );

        let t_project = profile_start!();
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
        let result = mem.scan(state, None, &[], None).await;
        profile_end!("9. projection + memtable", t_project);
        profile_end!(
            "TOTAL scan_with_transcript_engine",
            t_total,
            format!("{} VCF rows", total_vcf_rows)
        );
        if profiling_enabled() {
            eprintln!("[VEP_PROFILE] ====== scan_with_transcript_engine END ======");
        }
        result
    }

    fn annotate_batch_with_transcript_engine(
        &self,
        batch: &RecordBatch,
        engine: &TranscriptConsequenceEngine,
        ctx: &PreparedContext<'_>,
        colocated_map: &HashMap<ColocatedKey, ColocatedData>,
        sift_cache: &SiftPolyphenCache,
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
        let hgvs_flags = HgvsFlags::from_options_json(self.options_json.as_deref());
        let mut hgvs_reference_reader = if hgvs_flags.any() && hgvs_flags.shift_hgvs {
            self.options_json
                .as_deref()
                .and_then(|opts| Self::parse_json_string_option(opts, "reference_fasta_path"))
                .map(|path| {
                    fasta::io::indexed_reader::Builder::default()
                        .build_from_path(&path)
                        .map_err(|e| {
                            DataFusionError::Execution(format!(
                                "failed to open indexed reference FASTA '{path}': {e}"
                            ))
                        })
                })
                .transpose()?
        } else {
            None
        };

        let mut csq_builder = StringBuilder::with_capacity(batch.num_rows(), batch.num_rows() * 40);
        let mut most_builder =
            StringBuilder::with_capacity(batch.num_rows(), batch.num_rows() * 16);

        // Reusable buffers to avoid per-row/per-CSQ-entry String allocations.
        let mut csq_buf = String::with_capacity(4096);
        let mut terms_buf = String::with_capacity(128);

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

            // VEP skips star alleles entirely — no CSQ produced.
            if alt_allele == "*" {
                csq_builder.append_null();
                most_builder.append_null();
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
                            &flags,
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
            csq_buf.clear();
            if let Some(most_val) = &cached_most {
                use std::fmt::Write;
                // Cache hit — produce single CSQ entry with empty transcript fields.
                let csq_val = cached_csq.unwrap_or_default();
                let impact = SoTerm::from_str(most_val)
                    .map(|t| impact_label(t.impact()))
                    .unwrap_or_else(|| impact_label(SoImpact::Modifier));
                if flags.everything {
                    let _ = write!(
                        csq_buf,
                        "{vep_allele}|{csq_val}|{impact}|||||||||||||||{existing_var}||||\
                         {variant_class}|||||||||||||||||||||\
                         {batch3_suffix}|||||"
                    );
                } else {
                    let _ = write!(
                        csq_buf,
                        "{vep_allele}|{csq_val}|{impact}|||||||||||||||{existing_var}||||||||||||\
                         {variant_class}||||||||||||{batch3_suffix}"
                    );
                };
                most_str = most_val.clone();
            } else {
                use std::fmt::Write;
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

                // VEP skips star alleles entirely — no CSQ produced.
                if alt_allele == "*" {
                    csq_builder.append_null();
                    most_builder.append_null();
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

                // Build per-transcript CSQ entries into reusable buffer (already cleared above).
                for tc in &assignments {
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
                    let ccds = tx_opt.and_then(|tx| tx.ccds.as_deref()).unwrap_or("");
                    let swissprot = tx_opt.and_then(|tx| tx.swissprot.as_deref()).unwrap_or("");
                    let trembl_raw = tx_opt.and_then(|tx| tx.trembl.as_deref()).unwrap_or("");
                    let trembl = csq_escape(trembl_raw);
                    let uniparc = tx_opt.and_then(|tx| tx.uniparc.as_deref()).unwrap_or("");
                    let uniprot_isoform = tx_opt
                        .and_then(|tx| tx.uniprot_isoform.as_deref())
                        .unwrap_or("");

                    if flags.everything {
                        // HGVS_OFFSET: shift_length for the transcript's strand.
                        // VEP only emits HGVS_OFFSET when HGVSc was actually computed
                        // for this transcript variant allele. Non-transcript features
                        // (regulatory, intergenic) never get HGVS_OFFSET.
                        let hgvs_offset = if hgvs_flags.hgvsc && tc.hgvsc.is_some() {
                            let tx_strand = tx_opt.map(|tx| tx.strand).unwrap_or(1);
                            variant
                                .hgvs_shift_for_strand(tx_strand)
                                .filter(|s| s.shift_length > 0)
                                .map(|s| {
                                    // VEP emits negative HGVS_OFFSET for reverse-strand
                                    // transcripts (the shift direction is opposite).
                                    let signed = s.shift_length as i64;
                                    if tx_strand < 0 {
                                        (-signed).to_string()
                                    } else {
                                        signed.to_string()
                                    }
                                })
                                .unwrap_or_default()
                        } else {
                            String::new()
                        };
                        // MANE generic: VEP emits "MANE_Select" or "MANE_Plus_Clinical"
                        // depending on the transcript's MANE annotation.
                        // Traceability:
                        // - VEP OutputFactory.pm MANE output
                        //   https://github.com/Ensembl/ensembl-vep/blob/release/115/modules/Bio/EnsEMBL/VEP/OutputFactory.pm#L1548-L1560
                        let mane = if tx_opt.and_then(|tx| tx.mane_select.as_deref()).is_some() {
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
                        );
                        // DOMAINS: overlapping protein domain features.
                        // VEP gates DOMAINS on $pre->{coding} which requires
                        // a valid CDS coordinate mapping. Our cds_position is
                        // only set when the variant falls within the CDS region.
                        // Traceability:
                        // - VEP OutputFactory.pm line 1434: if($self->{domains} && $pre->{coding})
                        // - VEP BaseVariationFeatureOverlapAllele.pm _bvfo_preds lines 449-471
                        //   https://github.com/Ensembl/ensembl-variation/blob/release/115/modules/Bio/EnsEMBL/Variation/BaseVariationFeatureOverlapAllele.pm
                        let is_coding = tc.cds_position.as_deref().is_some_and(|s| !s.is_empty());
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
                                        Some((a.parse::<usize>().ok()?, b.parse::<usize>().ok()?))
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
                    } else {
                        // 74-field CSQ: 29 base + 12 Batch 1 + 33 Batch 3.
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
                if csq_buf.is_empty() {
                    let impact = impact_label(SoImpact::Modifier);
                    if flags.everything {
                        let _ = write!(
                            csq_buf,
                            "{vep_allele}|sequence_variant|{impact}|||||||||||||||{existing_var}||||\
                             {variant_class}|||||||||||||||||||||\
                             {batch3_suffix}|||||"
                        );
                    } else {
                        let _ = write!(
                            csq_buf,
                            "{vep_allele}|sequence_variant|{impact}|||||||||||||||{existing_var}||||||||||||\
                             {variant_class}||||||||||||{batch3_suffix}"
                        );
                    }
                }
            };

            csq_builder.append_value(&csq_buf);
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

/// Fill missing merged-mode `HGNC_ID` values for RefSeq rows using the paired
/// Ensembl transcript mapping in the merged cache, with a unique-gene-symbol
/// fallback when no explicit `refseq_id -> HGNC_ID` link exists.
///
/// Traceability:
/// - Ensembl VEP `OutputFactory::BaseTranscriptVariationAllele_to_output_hash()`
///   emits transcript `_gene_hgnc_id` directly
///   <https://github.com/Ensembl/ensembl-vep/blob/release/115/modules/Bio/EnsEMBL/VEP/OutputFactory.pm#L1474-L1484>
///
/// Our parquet cache materializes Ensembl and RefSeq transcript rows
/// separately. RefSeq rows often lack promoted `gene_hgnc_id`, while the paired
/// Ensembl row carries both `refseq_id` and `gene_hgnc_id`. This reconstructs
/// the VEP-visible merged transcript state without heuristics.
fn backfill_missing_hgnc_ids(transcripts: &mut [TranscriptFeature], refseq_ids: &[Option<String>]) {
    let mut refseq_gene_hgnc_by_id: HashMap<String, (Option<String>, String)> = HashMap::new();
    let mut symbol_to_hgnc_ids: HashMap<String, HashSet<String>> = HashMap::new();

    for (tx, refseq_id) in transcripts.iter().zip(refseq_ids.iter()) {
        let Some(hgnc_id) = tx.gene_hgnc_id.as_deref() else {
            continue;
        };
        if let Some(refseq_id) = refseq_id.as_deref() {
            refseq_gene_hgnc_by_id
                .entry(refseq_id.to_string())
                .or_insert_with(|| (tx.gene_symbol.clone(), hgnc_id.to_string()));
        }
        if let Some(symbol) = tx.gene_symbol.as_deref() {
            symbol_to_hgnc_ids
                .entry(symbol.to_string())
                .or_default()
                .insert(hgnc_id.to_string());
        }
    }

    let unique_hgnc_by_symbol: HashMap<String, String> = symbol_to_hgnc_ids
        .into_iter()
        .filter_map(|(symbol, ids)| {
            if ids.len() == 1 {
                ids.into_iter().next().map(|id| (symbol, id))
            } else {
                None
            }
        })
        .collect();

    for tx in transcripts.iter_mut() {
        if tx.gene_hgnc_id.is_some() {
            continue;
        }
        if let Some((mapped_symbol, hgnc_id)) =
            refseq_gene_hgnc_by_id.get(tx.transcript_id.as_str())
        {
            if mapped_symbol.as_deref() == tx.gene_symbol.as_deref() || tx.gene_symbol.is_none() {
                tx.gene_hgnc_id = Some(hgnc_id.clone());
                continue;
            }
        }
        if let Some(symbol) = tx.gene_symbol.as_deref() {
            if let Some(hgnc_id) = unique_hgnc_by_symbol.get(symbol) {
                tx.gene_hgnc_id = Some(hgnc_id.clone());
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
fn read_compact_predictions(col: &dyn Array, row: usize) -> Vec<CompactPrediction> {
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

/// Rebuild RefSeq/XM/XR CDS strings from the indexed genomic reference so the
/// transcript engine sees the same spliced reference sequence Ensembl Variation
/// uses through transcript objects during HGVS/coding consequence evaluation.
///
/// Traceability:
/// - Ensembl Variation `BaseTranscriptVariationAllele::_get_peptide_alleles()`
///   consumes transcript-derived CDS/peptide sequence rather than a VCF-local
///   allele string
///   <https://github.com/Ensembl/ensembl-variation/blob/release/115/modules/Bio/EnsEMBL/Variation/BaseTranscriptVariationAllele.pm#L367-L509>
///
/// In our merged parquet cache, a small subset of RefSeq transcript CDS strings
/// diverge from the genomic reference sequence that VEP's transcript objects
/// expose. Reconstructing the spliced CDS from the indexed FASTA restores the
/// same sequence basis without changing transcript/exon coordinates.
fn hydrate_refseq_translation_cds_from_reference<R>(
    reader: &mut fasta::io::indexed_reader::IndexedReader<R>,
    transcripts: &[TranscriptFeature],
    exons: &[ExonFeature],
    translations: &mut [TranslationFeature],
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
        let chrom = tx.chrom.strip_prefix("chr").unwrap_or(&tx.chrom);
        let overlaps_input = input_variant_intervals
            .get(chrom)
            .map(|intervals| interval_overlaps_any(intervals, tx.start, tx.end))
            .unwrap_or(false);
        if !overlaps_input {
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
        // Skip if we already have spliced_seq or cdna_seq.
        if tx.spliced_seq.is_some() {
            continue;
        }
        // Only need UTR for coding transcripts.
        let Some(coding_end) = tx.cdna_coding_end else {
            continue;
        };
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
        let Some(tx_exons) = exons_by_tx.get(tx.transcript_id.as_str()) else {
            continue;
        };
        // Only hydrate if exons extend past CDS (i.e., there IS a 3' UTR).
        // If total exonic length <= coding_end, there's no UTR to hydrate.
        let total_exonic: usize = tx_exons
            .iter()
            .map(|e| usize::try_from(e.end.saturating_sub(e.start).saturating_add(1)).unwrap_or(0))
            .sum();
        if total_exonic <= coding_end {
            continue;
        }
        // Read the entire transcript span in ONE FASTA query and extract
        // exon subsequences from it. This reduces ~8.7 FASTA reads per
        // transcript to just 1, cutting hydration I/O by ~8x.
        // For very large transcripts (>500KB), fall back to per-exon reads
        // to avoid excessive memory allocation.
        let tx_span_size = tx.end - tx.start + 1;
        if tx_span_size > 500_000 {
            // Per-exon fallback for large transcripts
            let mut genomic_cdna = String::new();
            let mut ok = true;
            for exon in tx_exons {
                let segment = read_reference_sequence(reader, chrom, exon.start, exon.end)?;
                let expected = usize::try_from(exon.end - exon.start + 1).unwrap_or_default();
                if segment.len() != expected {
                    ok = false;
                    break;
                }
                genomic_cdna.push_str(&segment);
            }
            if !ok || genomic_cdna.is_empty() {
                continue;
            }
            let cdna = if tx.strand >= 0 {
                genomic_cdna.to_ascii_uppercase()
            } else {
                let Some(rc) = reverse_complement_dna(&genomic_cdna) else {
                    continue;
                };
                rc
            };
            tx.cdna_seq = Some(cdna);
            continue;
        }
        let tx_span = read_reference_sequence(reader, chrom, tx.start, tx.end)?;
        let tx_span_len = usize::try_from(tx_span_size).unwrap_or_default();
        if tx_span.len() != tx_span_len {
            continue;
        }
        let mut genomic_cdna = String::new();
        let mut ok = true;
        for exon in tx_exons {
            let local_start = usize::try_from(exon.start - tx.start).unwrap_or(0);
            let local_end = usize::try_from(exon.end - tx.start + 1).unwrap_or(0);
            let Some(segment) = tx_span.get(local_start..local_end) else {
                ok = false;
                break;
            };
            genomic_cdna.push_str(segment);
        }
        if !ok || genomic_cdna.is_empty() {
            continue;
        }
        let cdna = if tx.strand >= 0 {
            genomic_cdna.to_ascii_uppercase()
        } else {
            let Some(rc) = reverse_complement_dna(&genomic_cdna) else {
                continue;
            };
            rc
        };
        tx.cdna_seq = Some(cdna);
    }
    Ok(())
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
    tx.source
        .as_deref()
        .and_then(normalize_source_label)
        .as_deref()
        == Some("RefSeq")
        || matches!(
            tx.transcript_id.as_bytes().get(..2),
            Some(b"NM") | Some(b"NR") | Some(b"XM") | Some(b"XR")
        )
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
                    &requested_columns,
                    extended_probes,
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
                if has_col(col_name) {
                    projected_exprs.push(format!(
                        "CAST(l.`cache_{col_name}` AS VARCHAR) AS `{col_name}`"
                    ));
                } else {
                    projected_exprs.push(format!("CAST(NULL AS VARCHAR) AS `{col_name}`"));
                }
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
    use crate::transcript_consequence::{
        CachedPredictions, ProteinDomainFeature, SiftPolyphenCache, TranslationFeature,
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
        let cache = make_sift_cache(vec![(
            "ENST00000001",
            42,
            "I",
            "deleterious",
            0.01,
            "probably damaging",
            0.999,
        )]);
        let (sift, polyphen) =
            lookup_sift_polyphen(Some("ENST00000001"), Some("42"), Some("V/I"), &cache);
        assert_eq!(sift, "deleterious(0.01)");
        assert_eq!(polyphen, "probably_damaging(0.999)");
    }

    #[test]
    fn test_lookup_sift_polyphen_non_substitution_skipped() {
        let cache = make_sift_cache(vec![(
            "ENST00000001",
            42,
            "I",
            "deleterious",
            0.01,
            "benign",
            0.0,
        )]);
        // Multi-char alt amino acid — not a single substitution.
        let (sift, polyphen) =
            lookup_sift_polyphen(Some("ENST00000001"), Some("42"), Some("V/IL"), &cache);
        assert!(sift.is_empty());
        assert!(polyphen.is_empty());

        // Range position — indel, should be skipped.
        let (sift, polyphen) =
            lookup_sift_polyphen(Some("ENST00000001"), Some("42-43"), Some("V/I"), &cache);
        assert!(sift.is_empty());
        assert!(polyphen.is_empty());
    }

    #[test]
    fn test_lookup_sift_polyphen_missing_transcript() {
        let cache = SiftPolyphenCache::new();
        let (sift, polyphen) =
            lookup_sift_polyphen(Some("ENST_MISSING"), Some("42"), Some("V/I"), &cache);
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
}
