//! Co-located variant metadata types plus allele-comparison and
//! genomic-shift helpers shared by the Lance variation-lookup path.
//!
//! Extracted from the former `variant_lookup_exec` module when the legacy
//! `VariantLookupExec` (Parquet/MemTable interval-join) lookup path was
//! removed. The Lance lookup (`lance_cache::lookup_exec`) and the annotation
//! engine (`annotate_provider`) consume these symbols.

use std::collections::{HashMap, HashSet};
use std::io::{BufRead, Seek};
use std::sync::{Arc, Mutex};

use datafusion::arrow::array::{
    Array, ArrayRef, LargeStringArray, RecordBatch, StringArray, StringViewArray,
};
use datafusion::common::{DataFusionError, Result};
use noodles_core::{Position, Region};
use noodles_fasta as fasta;

use crate::allele::{
    MatchedVariantAllele, VariantAlleleInput, allele_matches, get_matched_variant_alleles,
    vcf_to_vep_allele, vcf_to_vep_input_allele, vep_norm_end, vep_norm_start,
};

/// Zero-copy holder for a cache batch's AF columns, shared (by `Arc`) across all
/// co-located entries collected from that batch.
///
/// Holds one `ArrayRef` per `AF_COL_NAMES` entry (27); absent columns get an empty
/// placeholder. Cloning is an `Arc` ref-count bump — no string data is copied. A
/// row's AF value is read by index via [`AfColumns::value`], returning a borrowed
/// `&str` into the underlying Arrow buffer (`""` for null / absent / out-of-range).
#[derive(Debug, Clone)]
pub struct AfColumns(Arc<Vec<ArrayRef>>);

impl AfColumns {
    /// Wrap an already-collected set of column arrays (one per `AF_COL_NAMES`).
    pub fn new(columns: Vec<ArrayRef>) -> Self {
        AfColumns(Arc::new(columns))
    }

    /// Build from a record batch given per-column indices (`None` = column absent,
    /// stored as an empty placeholder so the index space stays aligned with
    /// `AF_COL_NAMES`). Each present column is an `Arc::clone` of the batch column —
    /// only these AF columns are retained, never the whole batch.
    pub fn from_batch(batch: &RecordBatch, af_indices: &[Option<usize>]) -> Self {
        let empty: ArrayRef = Arc::new(StringArray::from(Vec::<Option<&str>>::new()));
        let columns = af_indices
            .iter()
            .map(|idx| match idx {
                Some(i) => batch.column(*i).clone(),
                None => empty.clone(),
            })
            .collect();
        AfColumns(Arc::new(columns))
    }

    /// Number of AF columns (== `AF_COL_NAMES.len()` in production).
    pub fn len(&self) -> usize {
        self.0.len()
    }

    pub fn is_empty(&self) -> bool {
        self.0.is_empty()
    }

    /// True if both share the same underlying column allocation (no deep copy).
    pub fn ptr_eq(&self, other: &AfColumns) -> bool {
        Arc::ptr_eq(&self.0, &other.0)
    }

    /// Zero-copy raw AF value for `(col, row)`. Returns `""` for an absent column,
    /// out-of-range col/row, a null cell, or a non-string column — matching the
    /// former `Vec<String>` path's `unwrap_or_default()` (`String::new()`).
    pub fn value(&self, col: usize, row: usize) -> &str {
        let Some(array) = self.0.get(col) else {
            return "";
        };
        if row >= array.len() || array.is_null(row) {
            return "";
        }
        let array = array.as_ref();
        if let Some(a) = array.as_any().downcast_ref::<StringArray>() {
            a.value(row)
        } else if let Some(a) = array.as_any().downcast_ref::<StringViewArray>() {
            a.value(row)
        } else if let Some(a) = array.as_any().downcast_ref::<LargeStringArray>() {
            a.value(row)
        } else {
            ""
        }
    }

    /// Approximate retained heap size of the underlying column buffers (for probes).
    pub fn buffer_bytes(&self) -> usize {
        self.0
            .iter()
            .map(|a| a.get_array_memory_size())
            .sum::<usize>()
            + self.0.len() * std::mem::size_of::<ArrayRef>()
    }
}

/// A cache row's co-located metadata collected during streaming.
#[derive(Debug, Clone)]
pub struct ColocatedCacheEntry {
    pub variation_name: String,
    pub allele_string: String,
    pub matched_alleles: Vec<MatchedVariantAllele>,
    pub somatic: i64,
    pub pheno: i64,
    pub clin_sig: Option<String>,
    pub clin_sig_allele: Option<String>,
    pub pubmed: Option<String>,
    /// Zero-copy AF columns shared across all entries from the same cache batch.
    /// Indexed same as `AF_COL_NAMES`; read this row's value via [`Self::af_value`].
    pub af: AfColumns,
    /// This entry's row within `af`.
    pub af_row: u32,
}

impl ColocatedCacheEntry {
    /// Raw AF value for column `col` (`AF_COL_NAMES` order); `""` if null/absent.
    pub fn af_value(&self, col: usize) -> &str {
        self.af.value(col, self.af_row as usize)
    }

    /// Number of AF columns.
    pub fn af_len(&self) -> usize {
        self.af.len()
    }
}

/// AF column names collected into `ColocatedCacheEntry.af_values`.
/// Order must match `AF_COLUMNS` in `annotate_provider.rs`.
pub const AF_COL_NAMES: &[&str] = &[
    "AF",
    "AFR",
    "AMR",
    "EAS",
    "EUR",
    "SAS",
    "gnomADe",
    "gnomADe_AFR",
    "gnomADe_AMR",
    "gnomADe_ASJ",
    "gnomADe_EAS",
    "gnomADe_FIN",
    "gnomADe_MID",
    "gnomADe_NFE",
    "gnomADe_REMAINING",
    "gnomADe_SAS",
    "gnomADg",
    "gnomADg_AFR",
    "gnomADg_AMI",
    "gnomADg_AMR",
    "gnomADg_ASJ",
    "gnomADg_EAS",
    "gnomADg_FIN",
    "gnomADg_MID",
    "gnomADg_NFE",
    "gnomADg_REMAINING",
    "gnomADg_SAS",
];

pub(crate) fn output_allele_from_allele_string(allele_string: &str) -> Option<&str> {
    allele_string.split_once('/').map(|(_, alt)| alt)
}

#[derive(Clone, Copy)]
pub(crate) enum ShiftableIndelKind {
    Insertion,
    Deletion,
}

/// Traceability:
/// - Ensembl Variation `create_shift_hash()`
///   <https://github.com/Ensembl/ensembl-variation/blob/23c76f60b1592e4df86159cf5530bdc326120c3d/modules/Bio/EnsEMBL/Variation/TranscriptVariationAllele.pm#L365-L400>
///
/// Genomic shift state is only defined for simple insertion/deletion allele
/// strings. Substitutions and multi-ALT representations bypass this path.
pub(crate) fn parse_shiftable_indel(
    allele_string: &str,
) -> Option<(&str, &str, ShiftableIndelKind)> {
    let (ref_allele, alt_allele) = allele_string.split_once('/')?;
    if ref_allele == "-" && !alt_allele.is_empty() && alt_allele != "-" {
        return Some((ref_allele, alt_allele, ShiftableIndelKind::Insertion));
    }
    if alt_allele == "-" && !ref_allele.is_empty() && ref_allele != "-" {
        return Some((ref_allele, alt_allele, ShiftableIndelKind::Deletion));
    }
    None
}

/// Traceability:
/// - Ensembl Variation `create_shift_hash()`
///   <https://github.com/Ensembl/ensembl-variation/blob/23c76f60b1592e4df86159cf5530bdc326120c3d/modules/Bio/EnsEMBL/Variation/TranscriptVariationAllele.pm#L365-L400>
///
/// Shift-state construction queries the indexed genomic reference by absolute
/// chromosome coordinates before `_genomic_shift()` rotates repeat indels.
fn build_reference_region(chrom: &str, start: i64, end: i64) -> Result<Region> {
    let start = usize::try_from(start).map_err(|_| {
        DataFusionError::Execution(format!(
            "reference query start is negative or overflowed for {chrom}:{start}-{end}"
        ))
    })?;
    let end = usize::try_from(end).map_err(|_| {
        DataFusionError::Execution(format!(
            "reference query end is negative or overflowed for {chrom}:{start}-{end}"
        ))
    })?;
    let start = Position::try_from(start).map_err(|e| {
        DataFusionError::Execution(format!(
            "reference query start is invalid for {chrom}:{start}-{end}: {e}"
        ))
    })?;
    let end = Position::try_from(end).map_err(|e| {
        DataFusionError::Execution(format!(
            "reference query end is invalid for {chrom}:{start}-{end}: {e}"
        ))
    })?;
    Ok(Region::new(chrom, start..=end))
}

/// Traceability:
/// - Ensembl Variation `_genomic_shift()`
///   <https://github.com/Ensembl/ensembl-variation/blob/23c76f60b1592e4df86159cf5530bdc326120c3d/modules/Bio/EnsEMBL/Variation/TranscriptVariationAllele.pm#L411-L466>
///
/// The downstream flank used for VEP-style repeat shifting comes directly
/// from the indexed reference FASTA, not from local allele heuristics.
pub(crate) fn read_reference_sequence<R>(
    reader: &mut fasta::io::indexed_reader::IndexedReader<R>,
    chrom: &str,
    start: i64,
    end: i64,
) -> Result<String>
where
    R: BufRead + Seek,
{
    let region = build_reference_region(chrom, start, end)?;
    let record = reader.query(&region).map_err(|e| {
        DataFusionError::Execution(format!(
            "failed to query reference FASTA for {chrom}:{start}-{end}: {e}"
        ))
    })?;
    String::from_utf8(record.sequence().as_ref().to_vec()).map_err(|e| {
        DataFusionError::Execution(format!(
            "reference FASTA returned non-UTF8 sequence for {chrom}:{start}-{end}: {e}"
        ))
    })
}

/// Traceability:
/// - Ensembl Variation `perform_shift()`
///   <https://github.com/Ensembl/ensembl-variation/blob/23c76f60b1592e4df86159cf5530bdc326120c3d/modules/Bio/EnsEMBL/Variation/TranscriptVariationAllele.pm#L291-L351>
///
/// This is the source-equivalent positive-strand genomic branch of VEP's
/// indel shifting loop. It rotates the shifted sequence through the 3' flank
/// and advances genomic coordinates one base at a time until the next flank
/// base no longer matches.
pub(crate) fn perform_forward_genomic_shift(
    seq_to_check: &str,
    post_seq: &str,
    start: i64,
    end: i64,
) -> (usize, String, i64, i64) {
    let mut seq_to_check = seq_to_check.as_bytes().to_vec();
    let post_seq = post_seq.as_bytes();
    let indel_length = seq_to_check.len();
    let mut shift_length = 0usize;
    let mut start = start;
    let mut end = end;

    if indel_length == 0 || post_seq.len() < indel_length {
        return (
            0,
            String::from_utf8(seq_to_check).unwrap_or_default(),
            start,
            end,
        );
    }

    let loop_limiter = post_seq.len() - indel_length;
    for n in 0..=loop_limiter {
        let check_next = seq_to_check[0];
        if check_next != post_seq[n] {
            break;
        }

        shift_length += 1;
        seq_to_check.rotate_left(1);
        start += 1;
        end += 1;
    }

    (
        shift_length,
        String::from_utf8(seq_to_check).unwrap_or_default(),
        start,
        end,
    )
}

/// Traceability:
/// - Ensembl Variation `create_shift_hash()`
///   <https://github.com/Ensembl/ensembl-variation/blob/23c76f60b1592e4df86159cf5530bdc326120c3d/modules/Bio/EnsEMBL/Variation/TranscriptVariationAllele.pm#L365-L400>
/// - Ensembl Variation `_genomic_shift()`
///   <https://github.com/Ensembl/ensembl-variation/blob/23c76f60b1592e4df86159cf5530bdc326120c3d/modules/Bio/EnsEMBL/Variation/TranscriptVariationAllele.pm#L411-L466>
///
/// This materializes the VF-level genomic shift state VEP uses for colocated
/// matching: active compare space becomes the shifted indel representation,
/// while the original minimized representation is retained separately as
/// `unshifted_*`.
pub(crate) fn build_shifted_compare_state<R>(
    reader: &mut fasta::io::indexed_reader::IndexedReader<R>,
    chrom: &str,
    allele_string: &str,
    start: i64,
    end: i64,
) -> Result<Option<(String, i64, i64)>>
where
    R: BufRead + Seek,
{
    let Some((ref_allele, alt_allele, kind)) = parse_shiftable_indel(allele_string) else {
        return Ok(None);
    };

    let seq_to_check = match kind {
        ShiftableIndelKind::Insertion => alt_allele,
        ShiftableIndelKind::Deletion => ref_allele,
    };

    let flank_start = end + 1;
    if flank_start <= 0 {
        return Ok(None);
    }
    let flank_end = flank_start + 999;
    let post_seq = read_reference_sequence(reader, chrom, flank_start, flank_end)?;
    let (shift_length, shifted_seq, shifted_start, shifted_end) =
        perform_forward_genomic_shift(seq_to_check, &post_seq, start, end);

    if shift_length == 0 {
        return Ok(None);
    }

    let shifted_allele_string = match kind {
        ShiftableIndelKind::Insertion => format!("-/{shifted_seq}"),
        ShiftableIndelKind::Deletion => format!("{shifted_seq}/-"),
    };

    Ok(Some((shifted_allele_string, shifted_start, shifted_end)))
}

/// Two-pass allele matching for colocated variant collection.
///
/// This replicates VEP's `compare_existing()` logic: first match with the
/// active `VariationFeature` allele string, then (if unshifted state exists)
/// also match with unshifted alleles and merge the results.
///
/// Returns `None` when both passes produce zero matches (variant should be
/// skipped). Returns `Some(vec![])` for unknown-allele records that match on
/// exact coordinates.
pub(crate) fn compare_existing_variant_alleles(
    compare_allele_string: &str,
    compare_start: i64,
    compare_end: i64,
    unshifted_allele_string: Option<&str>,
    unshifted_start: Option<i64>,
    _unshifted_end: Option<i64>,
    existing_allele_string: &str,
    existing_start: i64,
    existing_end: i64,
) -> Option<Vec<MatchedVariantAllele>> {
    if !existing_allele_string.contains('/') {
        return (existing_start == compare_start && existing_end == compare_end)
            .then_some(Vec::new());
    }

    let mut matched_alleles = get_matched_variant_alleles(
        VariantAlleleInput {
            allele_string: compare_allele_string,
            pos: compare_start,
            strand: 1,
        },
        VariantAlleleInput {
            allele_string: existing_allele_string,
            pos: existing_start,
            strand: 1,
        },
    );

    if let (Some(unshifted_as), Some(unshifted_s)) = (unshifted_allele_string, unshifted_start) {
        let mut seen = matched_alleles.iter().cloned().collect::<HashSet<_>>();
        for matched in get_matched_variant_alleles(
            VariantAlleleInput {
                allele_string: unshifted_as,
                pos: unshifted_s,
                strand: 1,
            },
            VariantAlleleInput {
                allele_string: existing_allele_string,
                pos: existing_start,
                strand: 1,
            },
        ) {
            if seen.insert(matched.clone()) {
                matched_alleles.push(matched);
            }
        }
    }

    if matched_alleles.is_empty() {
        None
    } else {
        Some(matched_alleles)
    }
}

/// Shared sink key for co-located data collected during Lance variation lookup.
///
/// Key = VCF (chrom, input_start, input_end, input_allele_string), value =
/// cache entries attached to that specific parser/input allele.
pub type ColocatedKey = (String, i64, i64, String);

#[derive(Debug, Clone, Default)]
pub struct ColocatedSinkValue {
    pub entries: Vec<ColocatedCacheEntry>,
    /// Output allele component used by VEP's `add_colocated_variant_info()`
    /// fallback (`VariationFeature->{shifted_allele_string}`). This is derived
    /// from the active compare allele space, so IDs follow the same allele
    /// representation that produced `matched_alleles`.
    pub compare_output_allele: Option<String>,
    /// Original compare-space output allele component retained only when VEP
    /// defines genomic shift metadata. This is the only fallback used for
    /// frequency lookup, matching OutputFactory's `alt_orig_allele_string`.
    pub unshifted_output_allele: Option<String>,
}

/// Shared sink for co-located data collected during Lance variation lookup.
pub type ColocatedSink = Arc<Mutex<HashMap<ColocatedKey, ColocatedSinkValue>>>;

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn af_columns_zero_copy_value_handles_null_absent_and_empty() {
        use datafusion::arrow::array::{ArrayRef, StringArray};

        // col 0: ["T:0.93", null] ; col 1: ["", "A:0.1"]
        let c0: ArrayRef = Arc::new(StringArray::from(vec![Some("T:0.93"), None]));
        let c1: ArrayRef = Arc::new(StringArray::from(vec![Some(""), Some("A:0.1")]));
        let af = AfColumns::new(vec![c0, c1]);

        assert_eq!(af.len(), 2);
        // row 0
        assert_eq!(af.value(0, 0), "T:0.93");
        assert_eq!(af.value(1, 0), ""); // present but empty
        // row 1
        assert_eq!(af.value(0, 1), ""); // null -> ""
        assert_eq!(af.value(1, 1), "A:0.1");
        // absent column index -> ""
        assert_eq!(af.value(5, 0), "");
        // row out of range -> ""
        assert_eq!(af.value(0, 9), "");
    }

    #[test]
    fn output_allele_from_allele_string_extracts_alt_component() {
        assert_eq!(output_allele_from_allele_string("AA/-"), Some("-"));
        assert_eq!(output_allele_from_allele_string("-/AC"), Some("AC"));
        assert_eq!(output_allele_from_allele_string("AA"), None);
    }

    #[test]
    fn parse_shiftable_indel_only_accepts_simple_insertions_and_deletions() {
        assert!(matches!(
            parse_shiftable_indel("-/AC"),
            Some(("-", "AC", ShiftableIndelKind::Insertion))
        ));
        assert!(matches!(
            parse_shiftable_indel("AC/-"),
            Some(("AC", "-", ShiftableIndelKind::Deletion))
        ));
        assert!(parse_shiftable_indel("A/G").is_none());
        assert!(parse_shiftable_indel("A/G/T").is_none());
    }

    #[test]
    fn compare_existing_variant_alleles_matches_unknown_on_active_coords_only() {
        assert_eq!(
            compare_existing_variant_alleles(
                "-/TTTT",
                1735012,
                1735011,
                Some("-/TTTT"),
                Some(1735009),
                Some(1735008),
                "COSMIC_MUTATION",
                1735009,
                1735008,
            ),
            None
        );
    }

    #[test]
    fn compare_existing_variant_rejects_unknown_alleles_outside_both_coords() {
        assert_eq!(
            compare_existing_variant_alleles(
                "-/TTTT",
                1735012,
                1735011,
                Some("-/TTTT"),
                Some(1735009),
                Some(1735008),
                "COSMIC_MUTATION",
                1735007,
                1735006,
            ),
            None
        );
    }
}
