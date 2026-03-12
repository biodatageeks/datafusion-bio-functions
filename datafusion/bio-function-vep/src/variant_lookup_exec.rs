//! Self-contained left interval join ExecutionPlan for variant lookup.
//!
//! `VariantLookupExec` replaces the SQL-based join approach with a custom
//! physical operator that:
//! - Collects VCF rows (build side, small) into per-chromosome COITrees
//! - Streams cache rows (probe side, large), probing the trees for overlaps
//! - Applies `match_allele()` as a post-filter
//! - Emits unmatched VCF rows with NULL cache columns (LEFT JOIN semantics)

use std::any::Any;
use std::collections::{HashMap, HashSet};
use std::fmt::{Debug, Formatter};
use std::io::{BufRead, Seek};
use std::pin::Pin;
use std::sync::{Arc, Mutex};
use std::task::{Context, Poll};

use coitrees::{COITree, Interval, IntervalTree};
use datafusion::arrow::array::{
    Array, ArrayRef, Int8Array, Int16Array, Int32Array, Int64Array, LargeStringArray, RecordBatch,
    StringArray, StringViewArray, UInt32Array, UInt64Array,
};
use datafusion::arrow::datatypes::SchemaRef;
use datafusion::common::{DataFusionError, Result};
use datafusion::execution::{RecordBatchStream, SendableRecordBatchStream, TaskContext};
use datafusion::physical_expr::EquivalenceProperties;
use datafusion::physical_plan::execution_plan::{Boundedness, EmissionType};
use datafusion::physical_plan::{
    DisplayAs, DisplayFormatType, ExecutionPlan, ExecutionPlanProperties, PlanProperties,
};
use futures::{Stream, StreamExt};
use noodles_core::{Position, Region};
use noodles_fasta as fasta;

use crate::allele::{
    MatchedVariantAllele, VariantAlleleInput, allele_matches, get_matched_variant_alleles,
    vcf_to_vep_allele, vcf_to_vep_input_allele, vep_norm_end, vep_norm_start,
};
use crate::coordinate::CoordinateNormalizer;

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
    /// Raw AF column values (indexed same as `AF_COL_NAMES`).
    /// Each entry is the raw string like "T:0.9301" or "" if absent.
    pub af_values: Vec<String>,
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

fn read_optional_i64(
    batch: &RecordBatch,
    row: usize,
    column_idx: Option<usize>,
    name: &str,
) -> Option<i64> {
    column_idx.and_then(|idx| {
        let column = batch.column(idx);
        if column.is_null(row) {
            None
        } else {
            I64Accessor::new(column, name)
                .ok()
                .map(|accessor| accessor.value(row))
        }
    })
}

fn read_optional_string(
    batch: &RecordBatch,
    row: usize,
    column_idx: Option<usize>,
    name: &str,
) -> Option<String> {
    column_idx.and_then(|idx| {
        let column = batch.column(idx);
        if column.is_null(row) {
            None
        } else {
            StringAccessor::new(column, name)
                .ok()
                .map(|accessor| accessor.value(row).to_string())
        }
    })
}

fn push_unique_candidate(candidates: &mut Vec<usize>, seen: &mut HashSet<usize>, candidate: usize) {
    if seen.insert(candidate) {
        candidates.push(candidate);
    }
}

#[derive(Clone, Copy)]
enum ShiftableIndelKind {
    Insertion,
    Deletion,
}

fn parse_shiftable_indel(allele_string: &str) -> Option<(&str, &str, ShiftableIndelKind)> {
    let (ref_allele, alt_allele) = allele_string.split_once('/')?;
    if ref_allele == "-" && !alt_allele.is_empty() && alt_allele != "-" {
        return Some((ref_allele, alt_allele, ShiftableIndelKind::Insertion));
    }
    if alt_allele == "-" && !ref_allele.is_empty() && ref_allele != "-" {
        return Some((ref_allele, alt_allele, ShiftableIndelKind::Deletion));
    }
    None
}

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

fn read_reference_sequence<R>(
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
///   https://github.com/Ensembl/ensembl-variation/blob/23c76f60b1592e4df86159cf5530bdc326120c3d/modules/Bio/EnsEMBL/Variation/TranscriptVariationAllele.pm#L291-L351
///
/// This is the source-equivalent positive-strand genomic branch of VEP's
/// indel shifting loop. It rotates the shifted sequence through the 3' flank
/// and advances genomic coordinates one base at a time until the next flank
/// base no longer matches.
fn perform_forward_genomic_shift(
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
///   https://github.com/Ensembl/ensembl-variation/blob/23c76f60b1592e4df86159cf5530bdc326120c3d/modules/Bio/EnsEMBL/Variation/TranscriptVariationAllele.pm#L365-L400
/// - Ensembl Variation `_genomic_shift()`
///   https://github.com/Ensembl/ensembl-variation/blob/23c76f60b1592e4df86159cf5530bdc326120c3d/modules/Bio/EnsEMBL/Variation/TranscriptVariationAllele.pm#L411-L466
///
/// This materializes the VF-level genomic shift state VEP uses for colocated
/// matching: active compare space becomes the shifted indel representation,
/// while the original minimized representation is retained separately as
/// `unshifted_*`.
fn build_shifted_compare_state<R>(
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

/// Traceability:
/// - Ensembl VEP `InputBuffer::get_overlapping_vfs()`
///   https://github.com/Ensembl/ensembl-vep/blob/2beada0d57ca6234f467b14a6c60280f4a082717/modules/Bio/EnsEMBL/VEP/InputBuffer.pm#L311-L329
///
/// This mirrors VEP's prefilter by considering overlaps against both shifted
/// active coordinates and, when defined, unshifted minimized coordinates
/// before `compare_existing()` decides whether the existing variant is
/// allele-compatible.
fn collect_overlapping_candidates(
    build: &BuildSide,
    chrom: &str,
    start: i64,
    end: i64,
) -> Vec<usize> {
    let mut candidates = Vec::new();
    let mut seen = HashSet::new();
    let probe_start = start.min(end);
    let probe_end = start.max(end);

    if let Some(tree) = build.compare_trees.get(chrom) {
        tree.query(probe_start as i32, probe_end as i32, |interval| {
            push_unique_candidate(&mut candidates, &mut seen, *interval.metadata as usize);
        });
    }

    if let Some(tree) = build.unshifted_trees.get(chrom) {
        tree.query(probe_start as i32, probe_end as i32, |interval| {
            push_unique_candidate(&mut candidates, &mut seen, *interval.metadata as usize);
        });
    }

    candidates
}

/// Traceability:
/// - Ensembl VEP `AnnotationSource::Cache::VariationTabix::_annotate_pm()`
///   https://github.com/Ensembl/ensembl-vep/blob/2beada0d57ca6234f467b14a6c60280f4a082717/modules/Bio/EnsEMBL/VEP/AnnotationSource/Cache/VariationTabix.pm#L151-L189
///
/// The variation cache fetch path does not consider every overlapping existing
/// variant on the chromosome. Tabix is queried only for records whose START
/// coordinate falls inside the input VF window `start - 1 .. end + 1`, using
/// the active VF coordinates. We must replicate that prefilter after our
/// overlap-tree candidate collection, otherwise long existing variants that
/// begin before the VEP query window are incorrectly exposed to
/// `compare_existing()`.
fn existing_start_is_visible_to_input_row(input_row: &BuildRow, existing_start: i64) -> bool {
    let query_start = (input_row.compare_start - 1).min(input_row.compare_end + 1);
    let query_end = (input_row.compare_start - 1).max(input_row.compare_end + 1);
    existing_start >= query_start && existing_start <= query_end
}

/// Traceability:
/// - Ensembl VEP `compare_existing()`
///   https://github.com/Ensembl/ensembl-vep/blob/2beada0d57ca6234f467b14a6c60280f4a082717/modules/Bio/EnsEMBL/VEP/AnnotationType/Variation.pm#L146-L206
/// - Ensembl Variation `get_matched_variant_alleles()`
///   https://github.com/Ensembl/ensembl-variation/blob/23c76f60b1592e4df86159cf5530bdc326120c3d/modules/Bio/EnsEMBL/Variation/Utils/Sequence.pm#L1098-L1258
///
/// This ports VEP's existing-variant decision exactly for the offline cache
/// path: unknown-allele records match only on exact shifted coordinates, while
/// known-allele records are accepted only when shifted input alleles, or
/// explicitly defined unshifted input alleles, produce non-empty
/// `matched_alleles`.
fn compare_existing_variant(
    input_row: &BuildRow,
    existing_allele_string: &str,
    existing_start: i64,
    existing_end: i64,
) -> Option<Vec<MatchedVariantAllele>> {
    if !existing_allele_string.contains('/') {
        return (existing_start == input_row.compare_start
            && existing_end == input_row.compare_end)
            .then_some(Vec::new());
    }

    let mut matched_alleles = get_matched_variant_alleles(
        VariantAlleleInput {
            allele_string: &input_row.compare_allele_string,
            pos: input_row.compare_start,
            strand: 1,
        },
        VariantAlleleInput {
            allele_string: existing_allele_string,
            pos: existing_start,
            strand: 1,
        },
    );

    if let (Some(unshifted_allele_string), Some(unshifted_start)) = (
        input_row.unshifted_allele_string.as_deref(),
        input_row.unshifted_start,
    ) {
        let mut seen = matched_alleles.iter().cloned().collect::<HashSet<_>>();
        for matched in get_matched_variant_alleles(
            VariantAlleleInput {
                allele_string: unshifted_allele_string,
                pos: unshifted_start,
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

/// Shared sink key for co-located data collected during `VariantLookupExec`
/// streaming.
///
/// Key = VCF (chrom, input_start, input_end, input_allele_string), value =
/// cache entries attached to that specific parser/input allele.
pub type ColocatedKey = (String, i64, i64, String);

/// Shared sink for co-located data collected during `VariantLookupExec` streaming.
pub type ColocatedSink = Arc<Mutex<HashMap<ColocatedKey, Vec<ColocatedCacheEntry>>>>;

/// Physical execution plan for self-contained left interval join.
///
/// VCF (left/build side) is collected into COITrees keyed by chromosome.
/// Cache (right/probe side) is streamed, probing the trees for overlaps.
/// Unmatched VCF rows are emitted with NULL cache columns after the cache
/// stream is exhausted.
pub struct VariantLookupExec {
    vcf_input: Arc<dyn ExecutionPlan>,
    cache_input: Arc<dyn ExecutionPlan>,
    cache_columns: Vec<String>,
    vcf_has_chr: bool,
    coord_normalizer: CoordinateNormalizer,
    extended_probes: bool,
    allowed_failed: i64,
    reference_fasta_path: Option<String>,
    output_schema: SchemaRef,
    properties: PlanProperties,
    /// Optional sink for co-located data collected during probe phase.
    colocated_sink: Option<ColocatedSink>,
}

impl VariantLookupExec {
    #[allow(clippy::too_many_arguments)]
    pub fn new(
        vcf_input: Arc<dyn ExecutionPlan>,
        cache_input: Arc<dyn ExecutionPlan>,
        cache_columns: Vec<String>,
        vcf_has_chr: bool,
        coord_normalizer: CoordinateNormalizer,
        extended_probes: bool,
        allowed_failed: i64,
        reference_fasta_path: Option<String>,
        output_schema: SchemaRef,
    ) -> Self {
        let properties = PlanProperties::new(
            EquivalenceProperties::new(output_schema.clone()),
            vcf_input.output_partitioning().clone(),
            EmissionType::Incremental,
            Boundedness::Bounded,
        );

        Self {
            vcf_input,
            cache_input,
            cache_columns,
            vcf_has_chr,
            coord_normalizer,
            extended_probes,
            allowed_failed,
            reference_fasta_path,
            output_schema,
            properties,
            colocated_sink: None,
        }
    }

    /// Set the co-located data sink for piggybacked collection during probe.
    pub fn with_colocated_sink(mut self, sink: ColocatedSink) -> Self {
        self.colocated_sink = Some(sink);
        self
    }
}

impl Debug for VariantLookupExec {
    fn fmt(&self, f: &mut Formatter<'_>) -> std::fmt::Result {
        write!(
            f,
            "VariantLookupExec {{ cache_columns: {:?}, extended_probes: {}, allowed_failed: {} }}",
            self.cache_columns, self.extended_probes, self.allowed_failed
        )
    }
}

impl DisplayAs for VariantLookupExec {
    fn fmt_as(&self, _t: DisplayFormatType, f: &mut Formatter<'_>) -> std::fmt::Result {
        write!(
            f,
            "VariantLookupExec: columns={:?}, extended_probes={}, allowed_failed={}",
            self.cache_columns, self.extended_probes, self.allowed_failed
        )
    }
}

impl ExecutionPlan for VariantLookupExec {
    fn name(&self) -> &str {
        "VariantLookupExec"
    }

    fn as_any(&self) -> &dyn Any {
        self
    }

    fn schema(&self) -> SchemaRef {
        self.output_schema.clone()
    }

    fn properties(&self) -> &PlanProperties {
        &self.properties
    }

    fn children(&self) -> Vec<&Arc<dyn ExecutionPlan>> {
        vec![&self.vcf_input, &self.cache_input]
    }

    fn with_new_children(
        self: Arc<Self>,
        children: Vec<Arc<dyn ExecutionPlan>>,
    ) -> Result<Arc<dyn ExecutionPlan>> {
        assert_eq!(children.len(), 2);
        let mut exec = VariantLookupExec::new(
            children[0].clone(),
            children[1].clone(),
            self.cache_columns.clone(),
            self.vcf_has_chr,
            self.coord_normalizer.clone(),
            self.extended_probes,
            self.allowed_failed,
            self.reference_fasta_path.clone(),
            self.output_schema.clone(),
        );
        exec.colocated_sink = self.colocated_sink.clone();
        Ok(Arc::new(exec))
    }

    fn execute(
        &self,
        partition: usize,
        context: Arc<TaskContext>,
    ) -> Result<SendableRecordBatchStream> {
        let vcf_stream = self.vcf_input.execute(partition, context.clone())?;
        let cache_stream = self.cache_input.execute(partition, context)?;

        Ok(Box::pin(VariantLookupStream::new(
            vcf_stream,
            cache_stream,
            self.output_schema.clone(),
            self.cache_columns.clone(),
            self.vcf_has_chr,
            self.coord_normalizer.clone(),
            self.extended_probes,
            self.allowed_failed,
            self.reference_fasta_path.clone(),
            self.colocated_sink.clone(),
        )))
    }
}

// ---------------------------------------------------------------------------
// Stream implementation
// ---------------------------------------------------------------------------

/// Row data extracted from VCF for allele matching during probe.
struct BuildRow {
    vcf_ref: String,
    vcf_alt: String,
    /// Parser-level input allele string used for sink identity.
    input_allele_string: String,
    /// Parser-level input start used for sink identity.
    input_start: i64,
    /// Parser-level input end used for sink identity.
    input_end: i64,
    /// Minimized allele string used by VEP matched-alleles comparison.
    compare_allele_string: String,
    /// Minimized start used by VEP matched-alleles comparison.
    compare_start: i64,
    /// Minimized end used by exact unknown-allele shifted checks.
    compare_end: i64,
    /// Mirrors VEP's `unshifted_allele_string`, which is only defined when
    /// shifting logic has produced an original input representation to retain.
    unshifted_allele_string: Option<String>,
    /// VEP-normalized start (1-based). Stored for tree fallback coordinate check.
    vep_start: i64,
    /// VEP-normalized end (1-based, may be < start for insertions).
    vep_end: i64,
    /// Mirrors VEP's `unshifted_start`, only when defined upstream.
    unshifted_start: Option<i64>,
    /// Mirrors VEP's `unshifted_end`, only when defined upstream.
    unshifted_end: Option<i64>,
}

/// Materialized build side: VCF data + dual lookup indices.
///
/// Uses a two-level HashMap for fast O(1) exact lookups (primary),
/// with per-chromosome COITrees as fallback for shifted-coordinate indels.
struct BuildSide {
    /// All VCF batches concatenated.
    vcf_batch: RecordBatch,
    /// Per-row extracted data for allele matching.
    rows: Vec<BuildRow>,
    /// Primary: chrom → { (vep_norm_start, vep_norm_end) → [row_indices] }.
    /// Two-level to avoid per-row String allocation during probe — the outer
    /// lookup by `&str` is free.
    hash_index: HashMap<String, HashMap<(i64, i64), Vec<u32>>>,
    /// Fallback: per-chromosome interval trees for overlap queries.
    /// Intervals use raw 1-based coordinates (not VEP-normalized) so that
    /// shifted indels can be found via range overlap.
    trees: HashMap<String, COITree<u32, u32>>,
    /// Active minimized or shifted compare-space overlap trees used by VEP
    /// `get_overlapping_vfs()` before `compare_existing()` evaluates allele
    /// compatibility.
    compare_trees: HashMap<String, COITree<u32, u32>>,
    /// Unshifted input coordinates used by VEP's overlap prefilter before
    /// `compare_existing()` checks the unshifted allele string. Populated only
    /// for rows where upstream shift state exists.
    unshifted_trees: HashMap<String, COITree<u32, u32>>,
    /// Track which build rows have been matched.
    matched: Vec<bool>,
}

/// Cached column indices for the cache (probe) schema, resolved once on first batch.
struct CacheIndices {
    chrom: usize,
    start: usize,
    end: usize,
    allele_string: usize,
    failed: Option<usize>,
    /// Indices of requested cache columns in the cache batch schema.
    output_col_indices: Vec<usize>,
    /// Target data types for each output cache column (for casting).
    output_col_types: Vec<datafusion::arrow::datatypes::DataType>,
}

enum StreamState {
    /// Collecting VCF (build) batches.
    CollectBuild,
    /// Processing cache (probe) batches.
    ProcessProbe,
    /// Emitting unmatched build rows.
    EmitUnmatched,
    /// Done.
    Done,
}

/// Cached column indices for co-located metadata in the cache schema.
struct ColocIndices {
    variation_name: usize,
    somatic: Option<usize>,
    pheno: Option<usize>,
    clin_sig: Option<usize>,
    clin_sig_allele: Option<usize>,
    pubmed: Option<usize>,
    /// Column indices for AF fields (same order as `AF_COL_NAMES`).
    /// `None` if column not present in cache schema.
    af_indices: Vec<Option<usize>>,
}

struct VariantLookupStream {
    vcf_stream: Option<SendableRecordBatchStream>,
    cache_stream: Option<SendableRecordBatchStream>,
    schema: SchemaRef,
    cache_columns: Vec<String>,
    vcf_has_chr: bool,
    coord_normalizer: CoordinateNormalizer,
    /// When true, fall back to COITree on hash miss (for shifted indels).
    extended_probes: bool,
    /// Maximum allowed cache `failed` value.
    allowed_failed: i64,
    /// Optional indexed reference FASTA used to materialize Ensembl genomic
    /// shift state for colocated existing-variant matching.
    reference_fasta_path: Option<String>,
    state: StreamState,
    /// Accumulated VCF batches (before build side is materialized).
    vcf_batches: Vec<RecordBatch>,
    /// Materialized build side.
    build: Option<BuildSide>,
    /// Number of VCF columns in output schema.
    num_vcf_cols: usize,
    /// Cached cache-schema column indices (resolved on first probe batch).
    cache_indices: Option<CacheIndices>,
    /// Shared sink for co-located data.
    colocated_sink: Option<ColocatedSink>,
    /// Cached column indices for co-located collection.
    coloc_indices: Option<ColocIndices>,
}

impl VariantLookupStream {
    fn new(
        vcf_stream: SendableRecordBatchStream,
        cache_stream: SendableRecordBatchStream,
        schema: SchemaRef,
        cache_columns: Vec<String>,
        vcf_has_chr: bool,
        coord_normalizer: CoordinateNormalizer,
        extended_probes: bool,
        allowed_failed: i64,
        reference_fasta_path: Option<String>,
        colocated_sink: Option<ColocatedSink>,
    ) -> Self {
        let num_vcf_cols = schema.fields().len() - cache_columns.len();
        Self {
            vcf_stream: Some(vcf_stream),
            cache_stream: Some(cache_stream),
            schema,
            cache_columns,
            vcf_has_chr,
            coord_normalizer,
            extended_probes,
            allowed_failed,
            reference_fasta_path,
            state: StreamState::CollectBuild,
            vcf_batches: Vec::new(),
            build: None,
            num_vcf_cols,
            cache_indices: None,
            colocated_sink,
            coloc_indices: None,
        }
    }

    /// Resolve and cache column indices for the probe (cache) schema.
    fn resolve_cache_indices(&mut self, cache_schema: &SchemaRef) -> Result<()> {
        if self.cache_indices.is_some() {
            return Ok(());
        }

        let chrom = cache_schema.index_of("chrom")?;
        let start = cache_schema.index_of("start")?;
        let end = cache_schema.index_of("end")?;
        let allele_string = cache_schema.index_of("allele_string")?;
        let failed = cache_schema.index_of("failed").ok();

        let mut output_col_indices = Vec::with_capacity(self.cache_columns.len());
        let mut output_col_types = Vec::with_capacity(self.cache_columns.len());
        for (i, col_name) in self.cache_columns.iter().enumerate() {
            output_col_indices.push(cache_schema.index_of(col_name)?);
            output_col_types.push(self.schema.field(self.num_vcf_cols + i).data_type().clone());
        }

        self.cache_indices = Some(CacheIndices {
            chrom,
            start,
            end,
            allele_string,
            failed,
            output_col_indices,
            output_col_types,
        });
        Ok(())
    }

    /// Resolve co-located column indices from cache schema (best-effort).
    fn resolve_coloc_indices(&mut self, cache_schema: &SchemaRef) {
        if self.coloc_indices.is_some() || self.colocated_sink.is_none() {
            return;
        }
        let variation_name = match cache_schema.index_of("variation_name") {
            Ok(idx) => idx,
            Err(_) => return, // no variation_name column → nothing to collect
        };
        let af_indices: Vec<Option<usize>> = AF_COL_NAMES
            .iter()
            .map(|name| cache_schema.index_of(name).ok())
            .collect();
        self.coloc_indices = Some(ColocIndices {
            variation_name,
            somatic: cache_schema.index_of("somatic").ok(),
            pheno: cache_schema.index_of("phenotype_or_disease").ok(),
            clin_sig: cache_schema.index_of("clin_sig").ok(),
            clin_sig_allele: cache_schema.index_of("clin_sig_allele").ok(),
            pubmed: cache_schema.index_of("pubmed").ok(),
            af_indices,
        });
    }

    /// Traceability:
    /// - Ensembl VEP `InputBuffer::get_overlapping_vfs()`
    ///   https://github.com/Ensembl/ensembl-vep/blob/2beada0d57ca6234f467b14a6c60280f4a082717/modules/Bio/EnsEMBL/VEP/InputBuffer.pm#L311-L329
    /// - Ensembl VEP `compare_existing()`
    ///   https://github.com/Ensembl/ensembl-vep/blob/2beada0d57ca6234f467b14a6c60280f4a082717/modules/Bio/EnsEMBL/VEP/AnnotationType/Variation.pm#L168-L196
    ///
    /// The build side must retain both the shifted lookup coordinates and the
    /// upstream input state required for later shifted/unshifted matching.
    /// Materialize all VCF batches and build dual lookup indices (hash + COITree).
    fn materialize_build_side(&mut self) -> Result<()> {
        if self.vcf_batches.is_empty() {
            self.build = None;
            return Ok(());
        }

        let vcf_batch = if self.vcf_batches.len() == 1 {
            self.vcf_batches.pop().unwrap()
        } else {
            datafusion::arrow::compute::concat_batches(
                &self.vcf_batches[0].schema(),
                &self.vcf_batches,
            )
            .map_err(|e| DataFusionError::ArrowError(Box::new(e), None))?
        };
        self.vcf_batches.clear();

        let vcf_schema = vcf_batch.schema();
        let chrom_idx = vcf_schema.index_of("chrom")?;
        let start_idx = vcf_schema.index_of("start")?;
        let end_idx = vcf_schema.index_of("end")?;
        let ref_idx = vcf_schema.index_of("ref")?;
        let alt_idx = vcf_schema.index_of("alt")?;

        // Build side is small (~5-6M rows for WGS), Vec<String> is fine here.
        let chroms = as_string_values(vcf_batch.column(chrom_idx), "chrom")?;
        let starts = as_i64_values(vcf_batch.column(start_idx), "start")?;
        let ends = as_i64_values(vcf_batch.column(end_idx), "end")?;
        let refs = as_string_values(vcf_batch.column(ref_idx), "ref")?;
        let alts = as_string_values(vcf_batch.column(alt_idx), "alt")?;

        let num_rows = vcf_batch.num_rows();
        let mut reference_reader = self
            .reference_fasta_path
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
            .transpose()?;

        let mut rows = Vec::with_capacity(num_rows);
        let mut hash_index: HashMap<String, HashMap<(i64, i64), Vec<u32>>> = HashMap::new();
        let mut intervals_by_chrom: HashMap<String, Vec<Interval<u32>>> = HashMap::new();
        let mut compare_intervals_by_chrom: HashMap<String, Vec<Interval<u32>>> = HashMap::new();
        let mut unshifted_intervals_by_chrom: HashMap<String, Vec<Interval<u32>>> = HashMap::new();

        for i in 0..num_rows {
            let raw_chrom = &chroms[i];
            let norm_chrom = if self.vcf_has_chr {
                raw_chrom
                    .strip_prefix("chr")
                    .unwrap_or(raw_chrom)
                    .to_string()
            } else {
                raw_chrom.clone()
            };

            let start = starts[i];
            let end = ends[i];
            let vcf_ref = refs[i].clone();
            let vcf_alt = alts[i].clone();
            let one_based_start = if self.coord_normalizer.input_zero_based {
                start + 1
            } else {
                start
            };
            let one_based_end = end;
            let (input_ref, input_alt, input_start) =
                vcf_to_vep_input_allele(one_based_start, &vcf_ref, &vcf_alt);
            let input_allele_string = format!("{input_ref}/{input_alt}");
            let (compare_ref, compare_alt) = vcf_to_vep_allele(&vcf_ref, &vcf_alt);
            let compare_allele_string = format!("{compare_ref}/{compare_alt}");

            // VEP-normalized coordinates for the hash index (exact matching).
            let vep_start = vep_norm_start(one_based_start, &vcf_ref, &vcf_alt);
            let vep_end = vep_norm_end(one_based_start, &vcf_ref, &vcf_alt);

            // Hash index: chrom → (vep_norm_start, vep_norm_end) → [row_indices].
            hash_index
                .entry(norm_chrom.clone())
                .or_default()
                .entry((vep_start, vep_end))
                .or_default()
                .push(i as u32);

            // COITree fallback: use VEP-normalized coordinates (same space as
            // hash and cache). For insertions where vep_start > vep_end, the
            // interval is [vep_end, vep_start] to cover the insertion point.
            let tree_start = vep_start.min(vep_end);
            let tree_end = vep_start.max(vep_end);

            intervals_by_chrom
                .entry(norm_chrom.clone())
                .or_default()
                .push(Interval::new(tree_start as i32, tree_end as i32, i as u32));

            let mut active_compare_allele_string = compare_allele_string.clone();
            let mut active_compare_start = vep_start;
            let mut active_compare_end = vep_end;
            let mut unshifted_allele_string = None;
            let mut unshifted_start = None;
            let mut unshifted_end = None;

            if let Some(reader) = reference_reader.as_mut() {
                if let Some((shifted_allele_string, shifted_start, shifted_end)) =
                    build_shifted_compare_state(
                        reader,
                        &norm_chrom,
                        &compare_allele_string,
                        vep_start,
                        vep_end,
                    )?
                {
                    unshifted_allele_string = Some(compare_allele_string.clone());
                    unshifted_start = Some(vep_start);
                    unshifted_end = Some(vep_end);
                    active_compare_allele_string = shifted_allele_string;
                    active_compare_start = shifted_start;
                    active_compare_end = shifted_end;
                }
            }

            let compare_tree_start = active_compare_start.min(active_compare_end);
            let compare_tree_end = active_compare_start.max(active_compare_end);
            compare_intervals_by_chrom
                .entry(norm_chrom.clone())
                .or_default()
                .push(Interval::new(
                    compare_tree_start as i32,
                    compare_tree_end as i32,
                    i as u32,
                ));

            rows.push(BuildRow {
                vcf_ref,
                vcf_alt,
                input_allele_string,
                input_start,
                input_end: one_based_end,
                compare_allele_string: active_compare_allele_string,
                compare_start: active_compare_start,
                compare_end: active_compare_end,
                unshifted_allele_string,
                vep_start,
                vep_end,
                unshifted_start,
                unshifted_end,
            });
        }

        for (row_idx, row) in rows.iter().enumerate() {
            let Some(unshifted_start) = row.unshifted_start else {
                continue;
            };
            let Some(unshifted_end) = row.unshifted_end else {
                continue;
            };
            let chrom = if self.vcf_has_chr {
                chroms[row_idx]
                    .strip_prefix("chr")
                    .unwrap_or(&chroms[row_idx])
                    .to_string()
            } else {
                chroms[row_idx].clone()
            };
            let tree_start = unshifted_start.min(unshifted_end);
            let tree_end = unshifted_start.max(unshifted_end);
            unshifted_intervals_by_chrom
                .entry(chrom)
                .or_default()
                .push(Interval::new(
                    tree_start as i32,
                    tree_end as i32,
                    row_idx as u32,
                ));
        }

        let mut trees = HashMap::new();
        for (chrom, intervals) in intervals_by_chrom {
            trees.insert(chrom, COITree::new(&intervals));
        }

        let mut compare_trees = HashMap::new();
        for (chrom, intervals) in compare_intervals_by_chrom {
            compare_trees.insert(chrom, COITree::new(&intervals));
        }

        let mut unshifted_trees = HashMap::new();
        for (chrom, intervals) in unshifted_intervals_by_chrom {
            unshifted_trees.insert(chrom, COITree::new(&intervals));
        }

        let matched = vec![false; num_rows];

        self.build = Some(BuildSide {
            vcf_batch,
            rows,
            hash_index,
            trees,
            compare_trees,
            unshifted_trees,
            matched,
        });

        Ok(())
    }

    /// Process a single cache (probe) batch against the build-side indices.
    ///
    /// Strategy: try O(1) hash lookup first; on miss, fall back to COITree
    /// overlap query for shifted-coordinate indels.
    fn process_probe_batch(&mut self, cache_batch: &RecordBatch) -> Result<Option<RecordBatch>> {
        let cache_schema = cache_batch.schema();
        self.resolve_cache_indices(&cache_schema)?;
        self.resolve_coloc_indices(&cache_schema);
        let ci = self.cache_indices.as_ref().unwrap();

        let build = self.build.as_mut().unwrap();

        // Zero-copy accessors for the hot-loop columns.
        let cache_chroms = StringAccessor::new(cache_batch.column(ci.chrom), "chrom")?;
        let cache_starts = I64Accessor::new(cache_batch.column(ci.start), "start")?;
        let cache_ends = I64Accessor::new(cache_batch.column(ci.end), "end")?;
        let cache_alleles =
            StringAccessor::new(cache_batch.column(ci.allele_string), "allele_string")?;
        let cache_failed = ci
            .failed
            .map(|idx| I64Accessor::new(cache_batch.column(idx), "failed"))
            .transpose()?;

        let num_cache_rows = cache_batch.num_rows();
        let cache_zero_based = self.coord_normalizer.cache_zero_based;
        // Collect (vcf_row_idx, cache_row_idx) pairs for matched rows.
        let mut vcf_indices: Vec<u32> = Vec::new();
        let mut cache_indices_buf: Vec<u32> = Vec::new();

        for cache_row in 0..num_cache_rows {
            let cache_chrom = cache_chroms.value(cache_row);

            let cache_start_raw = cache_starts.value(cache_row);
            let cache_end_raw = cache_ends.value(cache_row);
            let failed = cache_failed
                .as_ref()
                .map(|accessor| accessor.value(cache_row))
                .unwrap_or(0);
            // Traceability:
            // - Ensembl VEP `filter_variation()`
            //   https://github.com/Ensembl/ensembl-vep/blob/2beada0d57ca6234f467b14a6c60280f4a082717/modules/Bio/EnsEMBL/VEP/AnnotationType/Variation.pm#L224-L227
            if failed > self.allowed_failed {
                continue;
            }
            let cache_start_1based = if cache_zero_based {
                cache_start_raw + 1
            } else {
                cache_start_raw
            };
            let cache_end_1based = cache_end_raw;

            let allele_str = cache_alleles.value(cache_row);

            // --- Primary: O(1) hash lookup by exact (chrom, start, end) ---
            // Outer lookup by &str avoids String allocation per row.
            let mut found_via_hash = false;
            if let Some(coord_map) = build.hash_index.get(cache_chrom) {
                if let Some(candidates) = coord_map.get(&(cache_start_1based, cache_end_1based)) {
                    for &vcf_idx in candidates {
                        let row = &build.rows[vcf_idx as usize];
                        if allele_matches(&row.vcf_ref, &row.vcf_alt, allele_str) {
                            vcf_indices.push(vcf_idx);
                            cache_indices_buf.push(cache_row as u32);
                            build.matched[vcf_idx as usize] = true;
                            found_via_hash = true;
                        }
                    }
                }
            }

            // --- Fallback: COITree overlap query for shifted indels ---
            // Only used when extended_probes=true. The tree uses VEP-normalized
            // coordinates, so after finding overlap candidates we verify that the
            // cache row's coordinates are compatible with the VCF variant's
            // VEP-normalized coordinates (either exact match or insertion-style
            // reversal where cache start/end == VCF vep_start/vep_end).
            if !found_via_hash && self.extended_probes {
                if let Some(tree) = build.trees.get(cache_chrom) {
                    let probe_start = cache_start_1based.min(cache_end_1based);
                    let probe_end = cache_start_1based.max(cache_end_1based);

                    tree.query(probe_start as i32, probe_end as i32, |interval| {
                        let vcf_idx = *interval.metadata;

                        // Skip VCF variants already matched via hash — tree
                        // fallback is only for variants the hash missed entirely.
                        // Without this guard, a cache row at a different position
                        // can pollute co-located data for already-matched variants.
                        if build.matched[vcf_idx as usize] {
                            return;
                        }

                        let row = &build.rows[vcf_idx as usize];

                        // Coordinate compatibility check:
                        // 1. Exact match (same as hash would find)
                        // 2. Swapped (insertion-style: cache start/end flipped)
                        // 3. Cache range contains VEP point (half-open vs closed)
                        let exact =
                            cache_start_1based == row.vep_start && cache_end_1based == row.vep_end;
                        let swapped =
                            cache_start_1based == row.vep_end && cache_end_1based == row.vep_start;
                        let contains =
                            cache_start_1based <= row.vep_start && cache_end_1based >= row.vep_end;
                        if !exact && !swapped && !contains {
                            return;
                        }

                        if allele_matches(&row.vcf_ref, &row.vcf_alt, allele_str) {
                            vcf_indices.push(vcf_idx);
                            cache_indices_buf.push(cache_row as u32);
                            build.matched[vcf_idx as usize] = true;
                        }
                    });
                }
            }
        }

        // --- Co-located data collection (piggybacked on same cache scan) ---
        // For every cache row with a non-empty variation_name, probe the VCF
        // COITree to find overlapping VCF positions and collect the entry.
        if let (Some(sink), Some(coloc_idx)) = (&self.colocated_sink, &self.coloc_indices) {
            let var_name_col = cache_batch.column(coloc_idx.variation_name);
            let var_names = StringAccessor::new(var_name_col, "variation_name")?;

            let mut local_buf: HashMap<ColocatedKey, Vec<ColocatedCacheEntry>> = HashMap::new();

            for cache_row in 0..num_cache_rows {
                if var_name_col.is_null(cache_row) {
                    continue;
                }
                let var_name = var_names.value(cache_row);
                if var_name.is_empty() {
                    continue;
                }

                let cache_chrom = cache_chroms.value(cache_row);
                let cache_start_raw = cache_starts.value(cache_row);
                let cache_end_raw = cache_ends.value(cache_row);
                let failed = cache_failed
                    .as_ref()
                    .map(|accessor| accessor.value(cache_row))
                    .unwrap_or(0);
                if failed > self.allowed_failed {
                    continue;
                }
                let cs1 = if cache_zero_based {
                    cache_start_raw + 1
                } else {
                    cache_start_raw
                };
                let ce1 = cache_end_raw;
                let allele_str = cache_alleles.value(cache_row);

                if !build.compare_trees.contains_key(cache_chrom)
                    && !build.unshifted_trees.contains_key(cache_chrom)
                {
                    continue;
                }
                for vcf_idx in collect_overlapping_candidates(build, cache_chrom, cs1, ce1) {
                    let vcf_row = &build.rows[vcf_idx];
                    if !existing_start_is_visible_to_input_row(vcf_row, cs1) {
                        continue;
                    }
                    let Some(matched_alleles) =
                        compare_existing_variant(vcf_row, allele_str, cs1, ce1)
                    else {
                        continue;
                    };

                    let key = (
                        cache_chrom.to_string(),
                        vcf_row.input_start,
                        vcf_row.input_end,
                        vcf_row.input_allele_string.clone(),
                    );
                    local_buf.entry(key).or_default().push(ColocatedCacheEntry {
                        variation_name: var_name.to_string(),
                        allele_string: allele_str.to_string(),
                        matched_alleles,
                        somatic: read_optional_i64(
                            cache_batch,
                            cache_row,
                            coloc_idx.somatic,
                            "somatic",
                        )
                        .unwrap_or(0),
                        pheno: read_optional_i64(
                            cache_batch,
                            cache_row,
                            coloc_idx.pheno,
                            "phenotype_or_disease",
                        )
                        .unwrap_or(0),
                        clin_sig: read_optional_string(
                            cache_batch,
                            cache_row,
                            coloc_idx.clin_sig,
                            "clin_sig",
                        ),
                        clin_sig_allele: read_optional_string(
                            cache_batch,
                            cache_row,
                            coloc_idx.clin_sig_allele,
                            "clin_sig_allele",
                        ),
                        pubmed: read_optional_string(
                            cache_batch,
                            cache_row,
                            coloc_idx.pubmed,
                            "pubmed",
                        ),
                        af_values: coloc_idx
                            .af_indices
                            .iter()
                            .map(|opt_idx| {
                                read_optional_string(cache_batch, cache_row, *opt_idx, "af")
                                    .unwrap_or_default()
                            })
                            .collect(),
                    });
                }
            }

            if !local_buf.is_empty() {
                let mut guard = sink.lock().unwrap();
                for (key, entries) in local_buf {
                    guard.entry(key).or_default().extend(entries);
                }
            }
        }

        if vcf_indices.is_empty() {
            return Ok(None);
        }

        // Build output batch: VCF columns via take + cache columns via take.
        let vcf_take = UInt32Array::from(vcf_indices);
        let cache_take = UInt32Array::from(cache_indices_buf);

        let mut output_columns: Vec<ArrayRef> = Vec::with_capacity(self.schema.fields().len());

        // VCF columns from build side.
        for col_idx in 0..self.num_vcf_cols {
            let taken =
                datafusion::arrow::compute::take(build.vcf_batch.column(col_idx), &vcf_take, None)
                    .map_err(|e| DataFusionError::ArrowError(Box::new(e), None))?;
            output_columns.push(taken);
        }

        // Cache columns — use pre-resolved indices and types.
        for (i, &cache_col_idx) in ci.output_col_indices.iter().enumerate() {
            let taken = datafusion::arrow::compute::take(
                cache_batch.column(cache_col_idx),
                &cache_take,
                None,
            )
            .map_err(|e| DataFusionError::ArrowError(Box::new(e), None))?;
            let target_type = &ci.output_col_types[i];
            if taken.data_type() != target_type {
                let casted = datafusion::arrow::compute::cast(&taken, target_type)
                    .map_err(|e| DataFusionError::ArrowError(Box::new(e), None))?;
                output_columns.push(casted);
            } else {
                output_columns.push(taken);
            }
        }

        RecordBatch::try_new(self.schema.clone(), output_columns)
            .map(Some)
            .map_err(|e| DataFusionError::ArrowError(Box::new(e), None))
    }

    /// Emit unmatched VCF rows with NULL cache columns.
    fn emit_unmatched(&mut self) -> Result<Option<RecordBatch>> {
        let build = match self.build.as_ref() {
            Some(b) => b,
            None => return Ok(None),
        };

        let unmatched_indices: Vec<u32> = build
            .matched
            .iter()
            .enumerate()
            .filter(|(_, m)| !**m)
            .map(|(i, _)| i as u32)
            .collect();

        if unmatched_indices.is_empty() {
            return Ok(None);
        }

        let take_indices = UInt32Array::from(unmatched_indices);
        let mut output_columns: Vec<ArrayRef> = Vec::with_capacity(self.schema.fields().len());

        for col_idx in 0..self.num_vcf_cols {
            let taken = datafusion::arrow::compute::take(
                build.vcf_batch.column(col_idx),
                &take_indices,
                None,
            )
            .map_err(|e| DataFusionError::ArrowError(Box::new(e), None))?;
            output_columns.push(taken);
        }

        let num_unmatched = take_indices.len();
        for i in 0..self.cache_columns.len() {
            let field = self.schema.field(self.num_vcf_cols + i);
            let null_array =
                datafusion::arrow::array::new_null_array(field.data_type(), num_unmatched);
            output_columns.push(null_array);
        }

        RecordBatch::try_new(self.schema.clone(), output_columns)
            .map(Some)
            .map_err(|e| DataFusionError::ArrowError(Box::new(e), None))
    }
}

impl Stream for VariantLookupStream {
    type Item = Result<RecordBatch>;

    fn poll_next(mut self: Pin<&mut Self>, cx: &mut Context<'_>) -> Poll<Option<Self::Item>> {
        loop {
            match self.state {
                StreamState::CollectBuild => {
                    let stream = self.vcf_stream.as_mut().unwrap();
                    match stream.poll_next_unpin(cx) {
                        Poll::Ready(Some(Ok(batch))) => {
                            if batch.num_rows() > 0 {
                                self.vcf_batches.push(batch);
                            }
                            continue;
                        }
                        Poll::Ready(Some(Err(e))) => return Poll::Ready(Some(Err(e))),
                        Poll::Ready(None) => {
                            self.vcf_stream = None;
                            if let Err(e) = self.materialize_build_side() {
                                return Poll::Ready(Some(Err(e)));
                            }
                            if self.build.is_none() {
                                self.state = StreamState::Done;
                                return Poll::Ready(None);
                            }
                            self.state = StreamState::ProcessProbe;
                            continue;
                        }
                        Poll::Pending => return Poll::Pending,
                    }
                }
                StreamState::ProcessProbe => {
                    let stream = self.cache_stream.as_mut().unwrap();
                    match stream.poll_next_unpin(cx) {
                        Poll::Ready(Some(Ok(batch))) => {
                            if batch.num_rows() == 0 {
                                continue;
                            }
                            match self.process_probe_batch(&batch) {
                                Ok(Some(output)) => return Poll::Ready(Some(Ok(output))),
                                Ok(None) => continue,
                                Err(e) => return Poll::Ready(Some(Err(e))),
                            }
                        }
                        Poll::Ready(Some(Err(e))) => return Poll::Ready(Some(Err(e))),
                        Poll::Ready(None) => {
                            self.cache_stream = None;
                            self.state = StreamState::EmitUnmatched;
                            continue;
                        }
                        Poll::Pending => return Poll::Pending,
                    }
                }
                StreamState::EmitUnmatched => {
                    self.state = StreamState::Done;
                    match self.emit_unmatched() {
                        Ok(Some(batch)) => return Poll::Ready(Some(Ok(batch))),
                        Ok(None) => return Poll::Ready(None),
                        Err(e) => return Poll::Ready(Some(Err(e))),
                    }
                }
                StreamState::Done => return Poll::Ready(None),
            }
        }
    }
}

impl RecordBatchStream for VariantLookupStream {
    fn schema(&self) -> SchemaRef {
        self.schema.clone()
    }
}

// ---------------------------------------------------------------------------
// Zero-copy accessors for Arrow string/integer columns
// ---------------------------------------------------------------------------

/// Zero-copy string accessor that wraps different Arrow string array types.
/// Returns `&str` without allocating, unlike `as_string_values()` which clones
/// every value into a `Vec<String>`.
enum StringAccessor<'a> {
    Utf8(&'a StringArray),
    Utf8View(&'a StringViewArray),
    LargeUtf8(&'a LargeStringArray),
}

impl<'a> StringAccessor<'a> {
    fn new(col: &'a ArrayRef, column_name: &str) -> Result<Self> {
        if let Some(arr) = col.as_any().downcast_ref::<StringArray>() {
            return Ok(Self::Utf8(arr));
        }
        if let Some(arr) = col.as_any().downcast_ref::<StringViewArray>() {
            return Ok(Self::Utf8View(arr));
        }
        if let Some(arr) = col.as_any().downcast_ref::<LargeStringArray>() {
            return Ok(Self::LargeUtf8(arr));
        }
        Err(DataFusionError::Execution(format!(
            "column '{column_name}' expected string array, got {:?}",
            col.data_type()
        )))
    }

    #[inline]
    fn value(&self, i: usize) -> &'a str {
        match self {
            Self::Utf8(arr) => {
                if arr.is_null(i) {
                    ""
                } else {
                    arr.value(i)
                }
            }
            Self::Utf8View(arr) => {
                if arr.is_null(i) {
                    ""
                } else {
                    arr.value(i)
                }
            }
            Self::LargeUtf8(arr) => {
                if arr.is_null(i) {
                    ""
                } else {
                    arr.value(i)
                }
            }
        }
    }
}

/// Zero-copy i64 accessor that wraps different Arrow integer array types.
enum I64Accessor<'a> {
    Int64(&'a Int64Array),
    Int32(&'a Int32Array),
    Int16(&'a Int16Array),
    Int8(&'a Int8Array),
    UInt64(&'a UInt64Array),
    UInt32(&'a UInt32Array),
}

impl<'a> I64Accessor<'a> {
    fn new(col: &'a ArrayRef, column_name: &str) -> Result<Self> {
        if let Some(arr) = col.as_any().downcast_ref::<Int64Array>() {
            return Ok(Self::Int64(arr));
        }
        if let Some(arr) = col.as_any().downcast_ref::<Int32Array>() {
            return Ok(Self::Int32(arr));
        }
        if let Some(arr) = col.as_any().downcast_ref::<Int16Array>() {
            return Ok(Self::Int16(arr));
        }
        if let Some(arr) = col.as_any().downcast_ref::<Int8Array>() {
            return Ok(Self::Int8(arr));
        }
        if let Some(arr) = col.as_any().downcast_ref::<UInt64Array>() {
            return Ok(Self::UInt64(arr));
        }
        if let Some(arr) = col.as_any().downcast_ref::<UInt32Array>() {
            return Ok(Self::UInt32(arr));
        }
        Err(DataFusionError::Execution(format!(
            "column '{column_name}' expected integer array, got {:?}",
            col.data_type()
        )))
    }

    #[inline]
    fn value(&self, i: usize) -> i64 {
        match self {
            Self::Int64(arr) => {
                if arr.is_null(i) {
                    0
                } else {
                    arr.value(i)
                }
            }
            Self::Int32(arr) => {
                if arr.is_null(i) {
                    0
                } else {
                    i64::from(arr.value(i))
                }
            }
            Self::Int16(arr) => {
                if arr.is_null(i) {
                    0
                } else {
                    i64::from(arr.value(i))
                }
            }
            Self::Int8(arr) => {
                if arr.is_null(i) {
                    0
                } else {
                    i64::from(arr.value(i))
                }
            }
            Self::UInt64(arr) => {
                if arr.is_null(i) {
                    0
                } else {
                    arr.value(i) as i64
                }
            }
            Self::UInt32(arr) => {
                if arr.is_null(i) {
                    0
                } else {
                    i64::from(arr.value(i))
                }
            }
        }
    }
}

// ---------------------------------------------------------------------------
// Helper functions for build-side column extraction (allocating, used once)
// ---------------------------------------------------------------------------

fn as_string_values(col: &ArrayRef, column_name: &str) -> Result<Vec<String>> {
    if let Some(arr) = col.as_any().downcast_ref::<StringArray>() {
        return Ok((0..arr.len())
            .map(|i| {
                if arr.is_null(i) {
                    String::new()
                } else {
                    arr.value(i).to_string()
                }
            })
            .collect());
    }
    if let Some(arr) = col.as_any().downcast_ref::<StringViewArray>() {
        return Ok((0..arr.len())
            .map(|i| {
                if arr.is_null(i) {
                    String::new()
                } else {
                    arr.value(i).to_string()
                }
            })
            .collect());
    }
    if let Some(arr) = col.as_any().downcast_ref::<LargeStringArray>() {
        return Ok((0..arr.len())
            .map(|i| {
                if arr.is_null(i) {
                    String::new()
                } else {
                    arr.value(i).to_string()
                }
            })
            .collect());
    }
    Err(DataFusionError::Execution(format!(
        "column '{column_name}' expected string array, got {:?}",
        col.data_type()
    )))
}

fn as_i64_values(col: &ArrayRef, column_name: &str) -> Result<Vec<i64>> {
    if let Some(arr) = col.as_any().downcast_ref::<Int64Array>() {
        return Ok((0..arr.len())
            .map(|i| if arr.is_null(i) { 0 } else { arr.value(i) })
            .collect());
    }
    if let Some(arr) = col.as_any().downcast_ref::<Int32Array>() {
        return Ok((0..arr.len())
            .map(|i| {
                if arr.is_null(i) {
                    0
                } else {
                    i64::from(arr.value(i))
                }
            })
            .collect());
    }
    if let Some(arr) = col.as_any().downcast_ref::<UInt64Array>() {
        return Ok((0..arr.len())
            .map(|i| {
                if arr.is_null(i) {
                    0
                } else {
                    arr.value(i) as i64
                }
            })
            .collect());
    }
    if let Some(arr) = col.as_any().downcast_ref::<UInt32Array>() {
        return Ok((0..arr.len())
            .map(|i| {
                if arr.is_null(i) {
                    0
                } else {
                    i64::from(arr.value(i))
                }
            })
            .collect());
    }
    Err(DataFusionError::Execution(format!(
        "column '{column_name}' expected integer array, got {:?}",
        col.data_type()
    )))
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::fs;
    use std::time::{SystemTime, UNIX_EPOCH};

    use datafusion::arrow::array::StringArray;
    use datafusion::arrow::datatypes::{DataType, Field, Schema};
    use datafusion::physical_plan::stream::RecordBatchStreamAdapter;

    fn empty_sendable_stream(schema: SchemaRef) -> SendableRecordBatchStream {
        let iter = futures::stream::iter(Vec::<Result<RecordBatch>>::new());
        Box::pin(RecordBatchStreamAdapter::new(schema, Box::pin(iter)))
    }

    fn write_reference_fasta(prefix: &str, chrom: &str, sequence: &str) -> String {
        let unique = SystemTime::now()
            .duration_since(UNIX_EPOCH)
            .unwrap()
            .as_nanos();
        let fasta_path = std::env::temp_dir().join(format!("{prefix}_{unique}.fa"));
        let header = format!(">{chrom}\n");
        fs::write(&fasta_path, format!("{header}{sequence}\n")).unwrap();
        fs::write(
            fasta_path.with_extension("fa.fai"),
            format!(
                "{chrom}\t{}\t{}\t{}\t{}\n",
                sequence.len(),
                header.len(),
                sequence.len(),
                sequence.len() + 1,
            ),
        )
        .unwrap();
        fasta_path.display().to_string()
    }

    #[test]
    fn compare_existing_variant_matches_shifted_or_unshifted_input() {
        let row = BuildRow {
            vcf_ref: "AAA".to_string(),
            vcf_alt: "A".to_string(),
            input_allele_string: "AA/-".to_string(),
            input_start: 101,
            input_end: 102,
            compare_allele_string: "AA/-".to_string(),
            compare_start: 101,
            compare_end: 102,
            unshifted_allele_string: Some("AAA/A".to_string()),
            vep_start: 101,
            vep_end: 102,
            unshifted_start: Some(100),
            unshifted_end: Some(102),
        };

        let matched = compare_existing_variant(&row, "AA/-", 101, 102).unwrap();
        assert_eq!(
            matched,
            vec![
                MatchedVariantAllele {
                    a_allele: "-".to_string(),
                    a_index: 0,
                    b_allele: "-".to_string(),
                    b_index: 0,
                },
                MatchedVariantAllele {
                    a_allele: "A".to_string(),
                    a_index: 0,
                    b_allele: "-".to_string(),
                    b_index: 0,
                },
            ]
        );
    }

    #[test]
    fn compare_existing_variant_allows_unknown_alleles_on_exact_shifted_coords_only() {
        let row = BuildRow {
            vcf_ref: "ACGT".to_string(),
            vcf_alt: "A".to_string(),
            input_allele_string: "CGT/-".to_string(),
            input_start: 101,
            input_end: 103,
            compare_allele_string: "CGT/-".to_string(),
            compare_start: 101,
            compare_end: 103,
            unshifted_allele_string: Some("ACGT/A".to_string()),
            vep_start: 101,
            vep_end: 103,
            unshifted_start: Some(100),
            unshifted_end: Some(103),
        };

        assert_eq!(
            compare_existing_variant(&row, "HGMD_MUTATION", 101, 103),
            Some(Vec::new())
        );
        assert_eq!(
            compare_existing_variant(&row, "HGMD_MUTATION", 100, 103),
            None
        );
    }

    #[test]
    fn compare_existing_variant_does_not_use_unshifted_matching_when_state_is_absent() {
        let row = BuildRow {
            vcf_ref: "AAA".to_string(),
            vcf_alt: "A".to_string(),
            input_allele_string: "AA/-".to_string(),
            input_start: 101,
            input_end: 102,
            compare_allele_string: "AA/-".to_string(),
            compare_start: 101,
            compare_end: 102,
            unshifted_allele_string: None,
            vep_start: 101,
            vep_end: 102,
            unshifted_start: None,
            unshifted_end: None,
        };

        assert_eq!(compare_existing_variant(&row, "AA/-", 100, 102), None);
    }

    #[test]
    fn compare_existing_variant_uses_compare_coords_for_unknown_insertions() {
        let row = BuildRow {
            vcf_ref: "TTA".to_string(),
            vcf_alt: "TATATATA".to_string(),
            input_allele_string: "TA/ATATATA".to_string(),
            input_start: 119247098,
            input_end: 119247099,
            compare_allele_string: "-/ATATA".to_string(),
            compare_start: 119247098,
            compare_end: 119247097,
            unshifted_allele_string: None,
            vep_start: 119247098,
            vep_end: 119247097,
            unshifted_start: None,
            unshifted_end: None,
        };

        assert_eq!(
            compare_existing_variant(&row, "HGMD_MUTATION", 119247098, 119247097),
            Some(Vec::new())
        );
        assert_eq!(
            compare_existing_variant(&row, "HGMD_MUTATION", 119247098, 119247099),
            None
        );
    }

    #[test]
    fn collect_overlapping_candidates_uses_active_compare_coordinates() {
        let schema = Arc::new(datafusion::arrow::datatypes::Schema::empty());
        let empty_batch = RecordBatch::new_empty(schema);
        let build = BuildSide {
            vcf_batch: empty_batch,
            rows: vec![BuildRow {
                vcf_ref: "CATACATATATATATATATATATATAT".to_string(),
                vcf_alt: "CATATATATATATAT".to_string(),
                input_allele_string: "ATACATATATATATATATATATATAT/ATATATATATATAT".to_string(),
                input_start: 62689176,
                input_end: 62689202,
                compare_allele_string: "ACATATATATATATATATATATAT/-".to_string(),
                compare_start: 62689177,
                compare_end: 62689188,
                unshifted_allele_string: None,
                vep_start: 62689177,
                vep_end: 62689188,
                unshifted_start: None,
                unshifted_end: None,
            }],
            hash_index: HashMap::new(),
            trees: HashMap::from([(
                "1".to_string(),
                COITree::new(&[Interval::new(62689177, 62689188, 0)]),
            )]),
            compare_trees: HashMap::from([(
                "1".to_string(),
                COITree::new(&[Interval::new(62689177, 62689188, 0)]),
            )]),
            unshifted_trees: HashMap::new(),
            matched: vec![false],
        };

        assert_eq!(
            collect_overlapping_candidates(&build, "1", 62689177, 62689188),
            vec![0]
        );
    }

    #[test]
    fn compare_existing_variant_uses_minimized_compare_allele_space_for_repeat_insertions() {
        let row = BuildRow {
            vcf_ref: "TTA".to_string(),
            vcf_alt: "TATATATA".to_string(),
            input_allele_string: "TA/ATATATA".to_string(),
            input_start: 119247098,
            input_end: 119247099,
            compare_allele_string: "-/ATATA".to_string(),
            compare_start: 119247098,
            compare_end: 119247097,
            unshifted_allele_string: None,
            vep_start: 119247098,
            vep_end: 119247097,
            unshifted_start: None,
            unshifted_end: None,
        };

        let matched =
            compare_existing_variant(&row, "-/A/ATA/ATATA/ATATATA", 119247098, 119247097).unwrap();

        assert_eq!(
            matched,
            vec![MatchedVariantAllele {
                a_allele: "ATATA".to_string(),
                a_index: 0,
                b_allele: "ATATA".to_string(),
                b_index: 2,
            }]
        );
    }

    #[test]
    fn materialize_build_side_tracks_compare_and_parser_spaces_separately() {
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
                Arc::new(Int64Array::from(vec![119247097])),
                Arc::new(Int64Array::from(vec![119247099])),
                Arc::new(StringArray::from(vec!["TTA"])),
                Arc::new(StringArray::from(vec!["TATATATA"])),
            ],
        )
        .unwrap();

        let mut stream = VariantLookupStream::new(
            empty_sendable_stream(schema.clone()),
            empty_sendable_stream(schema.clone()),
            schema,
            Vec::new(),
            false,
            CoordinateNormalizer::new(false, false),
            true,
            0,
            None,
            None,
        );
        stream.vcf_batches.push(batch);

        stream.materialize_build_side().unwrap();
        let build = stream.build.as_ref().unwrap();
        let row = &build.rows[0];

        assert_eq!(row.input_allele_string, "TA/ATATATA");
        assert_eq!(row.input_start, 119247098);
        assert_eq!(row.input_end, 119247099);
        assert_eq!(row.compare_allele_string, "-/ATATA");
        assert_eq!(row.compare_start, 119247098);
        assert_eq!(row.compare_end, 119247097);
        assert_eq!(row.vep_start, 119247098);
        assert_eq!(row.vep_end, 119247097);
    }

    #[test]
    fn materialize_build_side_materializes_shifted_and_unshifted_compare_state() {
        let fasta_path = write_reference_fasta("repeat_shift", "1", "ACACACACAC");
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
                Arc::new(Int64Array::from(vec![1])),
                Arc::new(Int64Array::from(vec![3])),
                Arc::new(StringArray::from(vec!["ACA"])),
                Arc::new(StringArray::from(vec!["A"])),
            ],
        )
        .unwrap();

        let mut stream = VariantLookupStream::new(
            empty_sendable_stream(schema.clone()),
            empty_sendable_stream(schema.clone()),
            schema,
            Vec::new(),
            false,
            CoordinateNormalizer::new(false, false),
            true,
            0,
            Some(fasta_path),
            None,
        );
        stream.vcf_batches.push(batch);

        stream.materialize_build_side().unwrap();
        let build = stream.build.as_ref().unwrap();
        let row = &build.rows[0];

        assert_eq!(row.compare_allele_string, "CA/-");
        assert_eq!(row.compare_start, 8);
        assert_eq!(row.compare_end, 9);
        assert_eq!(row.unshifted_allele_string.as_deref(), Some("CA/-"));
        assert_eq!(row.unshifted_start, Some(2));
        assert_eq!(row.unshifted_end, Some(3));
        assert_eq!(collect_overlapping_candidates(build, "1", 8, 9), vec![0]);
        assert_eq!(collect_overlapping_candidates(build, "1", 2, 3), vec![0]);
    }

    #[test]
    fn existing_start_visibility_matches_variation_tabix_query_window() {
        let deletion_row = BuildRow {
            vcf_ref: "CAACAACAAAAAA".to_string(),
            vcf_alt: "CAAAA".to_string(),
            input_allele_string: "AACAACAAAAAA/AAAA".to_string(),
            input_start: 27971600,
            input_end: 27971611,
            compare_allele_string: "CAACAAAA/-".to_string(),
            compare_start: 27971602,
            compare_end: 27971609,
            unshifted_allele_string: None,
            vep_start: 27971602,
            vep_end: 27971609,
            unshifted_start: None,
            unshifted_end: None,
        };

        assert!(existing_start_is_visible_to_input_row(
            &deletion_row,
            27971601
        ));
        assert!(existing_start_is_visible_to_input_row(
            &deletion_row,
            27971610
        ));
        assert!(!existing_start_is_visible_to_input_row(
            &deletion_row,
            27971600
        ));

        let insertion_row = BuildRow {
            vcf_ref: "A".to_string(),
            vcf_alt: "ATT".to_string(),
            input_allele_string: "-/TT".to_string(),
            input_start: 101,
            input_end: 100,
            compare_allele_string: "-/TT".to_string(),
            compare_start: 101,
            compare_end: 100,
            unshifted_allele_string: None,
            vep_start: 101,
            vep_end: 100,
            unshifted_start: None,
            unshifted_end: None,
        };

        assert!(existing_start_is_visible_to_input_row(&insertion_row, 100));
        assert!(existing_start_is_visible_to_input_row(&insertion_row, 101));
        assert!(!existing_start_is_visible_to_input_row(&insertion_row, 99));
    }
}
