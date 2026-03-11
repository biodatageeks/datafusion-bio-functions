//! Self-contained left interval join ExecutionPlan for variant lookup.
//!
//! `VariantLookupExec` replaces the SQL-based join approach with a custom
//! physical operator that:
//! - Collects VCF rows (build side, small) into per-chromosome COITrees
//! - Streams cache rows (probe side, large), probing the trees for overlaps
//! - Applies `match_allele()` as a post-filter
//! - Emits unmatched VCF rows with NULL cache columns (LEFT JOIN semantics)

use std::any::Any;
use std::collections::HashMap;
use std::fmt::{Debug, Formatter};
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

use crate::allele::{allele_matches, vep_norm_end, vep_norm_start};
use crate::coordinate::CoordinateNormalizer;

/// A cache row's co-located metadata collected during streaming.
#[derive(Debug, Clone)]
pub struct ColocatedCacheEntry {
    pub variation_name: String,
    pub allele_string: String,
    pub somatic: i64,
    pub pheno: i64,
    pub clin_sig_allele: Option<String>,
    pub pubmed: Option<String>,
}

/// Shared sink for co-located data collected during `VariantLookupExec` streaming.
/// Key = VCF (chrom, vep_start, vep_end), value = cache entries at that position.
pub type ColocatedSink = Arc<Mutex<HashMap<(String, i64, i64), Vec<ColocatedCacheEntry>>>>;

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
            "VariantLookupExec {{ cache_columns: {:?}, extended_probes: {} }}",
            self.cache_columns, self.extended_probes
        )
    }
}

impl DisplayAs for VariantLookupExec {
    fn fmt_as(&self, _t: DisplayFormatType, f: &mut Formatter<'_>) -> std::fmt::Result {
        write!(
            f,
            "VariantLookupExec: columns={:?}, extended_probes={}",
            self.cache_columns, self.extended_probes
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
    /// VEP-normalized start (1-based). Stored for tree fallback coordinate check.
    vep_start: i64,
    /// VEP-normalized end (1-based, may be < start for insertions).
    vep_end: i64,
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
    /// Track which build rows have been matched.
    matched: Vec<bool>,
}

/// Cached column indices for the cache (probe) schema, resolved once on first batch.
struct CacheIndices {
    chrom: usize,
    start: usize,
    end: usize,
    allele_string: usize,
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
    allele_string: usize,
    somatic: Option<usize>,
    pheno: Option<usize>,
    clin_sig_allele: Option<usize>,
    pubmed: Option<usize>,
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
        self.coloc_indices = Some(ColocIndices {
            variation_name,
            allele_string: cache_schema.index_of("allele_string").unwrap_or(usize::MAX),
            somatic: cache_schema.index_of("somatic").ok(),
            pheno: cache_schema.index_of("phenotype_or_disease").ok(),
            clin_sig_allele: cache_schema.index_of("clin_sig_allele").ok(),
            pubmed: cache_schema.index_of("pubmed").ok(),
        });
    }

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
        let mut rows = Vec::with_capacity(num_rows);
        let mut hash_index: HashMap<String, HashMap<(i64, i64), Vec<u32>>> = HashMap::new();
        let mut intervals_by_chrom: HashMap<String, Vec<Interval<u32>>> = HashMap::new();

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

            rows.push(BuildRow {
                vcf_ref,
                vcf_alt,
                vep_start,
                vep_end,
            });
        }

        let mut trees = HashMap::new();
        for (chrom, intervals) in intervals_by_chrom {
            trees.insert(chrom, COITree::new(&intervals));
        }

        let matched = vec![false; num_rows];

        self.build = Some(BuildSide {
            vcf_batch,
            rows,
            hash_index,
            trees,
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

        let num_cache_rows = cache_batch.num_rows();
        let cache_zero_based = self.coord_normalizer.cache_zero_based;
        // Collect (vcf_row_idx, cache_row_idx) pairs for matched rows.
        let mut vcf_indices: Vec<u32> = Vec::new();
        let mut cache_indices_buf: Vec<u32> = Vec::new();

        for cache_row in 0..num_cache_rows {
            let cache_chrom = cache_chroms.value(cache_row);

            let cache_start_raw = cache_starts.value(cache_row);
            let cache_end_raw = cache_ends.value(cache_row);
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
                        let exact = cache_start_1based == row.vep_start
                            && cache_end_1based == row.vep_end;
                        let swapped = cache_start_1based == row.vep_end
                            && cache_end_1based == row.vep_start;
                        let contains = cache_start_1based <= row.vep_start
                            && cache_end_1based >= row.vep_end;
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
        if let (Some(sink), Some(coloc_idx)) =
            (&self.colocated_sink, &self.coloc_indices)
        {
            let var_name_col = cache_batch.column(coloc_idx.variation_name);
            let var_names = StringAccessor::new(var_name_col, "variation_name")?;

            let mut local_buf: HashMap<(String, i64, i64), Vec<ColocatedCacheEntry>> =
                HashMap::new();

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
                let cs1 = if cache_zero_based {
                    cache_start_raw + 1
                } else {
                    cache_start_raw
                };
                let ce1 = cache_end_raw;

                let Some(tree) = build.trees.get(cache_chrom) else {
                    continue;
                };

                let probe_s = cs1.min(ce1);
                let probe_e = cs1.max(ce1);

                // Build entry lazily (only when tree has hits).
                let mut entry: Option<ColocatedCacheEntry> = None;

                tree.query(probe_s as i32, probe_e as i32, |interval| {
                    let vcf_idx = *interval.metadata as usize;
                    let vcf_row = &build.rows[vcf_idx];

                    // Co-located = exact position match (same as old SQL equi-join).
                    // Unlike the main allele-matched join, co-located collection
                    // requires cache (start,end) == VCF (vep_start,vep_end).
                    if cs1 != vcf_row.vep_start || ce1 != vcf_row.vep_end {
                        return;
                    }

                    let e = entry.get_or_insert_with(|| {
                        let allele_str_val = if coloc_idx.allele_string != usize::MAX {
                            StringAccessor::new(
                                cache_batch.column(coloc_idx.allele_string),
                                "allele_string",
                            )
                            .map(|a| a.value(cache_row).to_string())
                            .unwrap_or_default()
                        } else {
                            String::new()
                        };
                        ColocatedCacheEntry {
                            variation_name: var_name.to_string(),
                            allele_string: allele_str_val,
                            somatic: coloc_idx
                                .somatic
                                .and_then(|idx| {
                                    I64Accessor::new(cache_batch.column(idx), "somatic")
                                        .ok()
                                        .map(|a| a.value(cache_row))
                                })
                                .unwrap_or(0),
                            pheno: coloc_idx
                                .pheno
                                .and_then(|idx| {
                                    I64Accessor::new(cache_batch.column(idx), "pheno")
                                        .ok()
                                        .map(|a| a.value(cache_row))
                                })
                                .unwrap_or(0),
                            clin_sig_allele: coloc_idx.clin_sig_allele.and_then(|idx| {
                                let col = cache_batch.column(idx);
                                if col.is_null(cache_row) {
                                    return None;
                                }
                                StringAccessor::new(col, "clin_sig_allele")
                                    .ok()
                                    .map(|a| a.value(cache_row).to_string())
                            }),
                            pubmed: coloc_idx.pubmed.and_then(|idx| {
                                let col = cache_batch.column(idx);
                                if col.is_null(cache_row) {
                                    return None;
                                }
                                StringAccessor::new(col, "pubmed")
                                    .ok()
                                    .map(|a| a.value(cache_row).to_string())
                            }),
                        }
                    });

                    let key =
                        (cache_chrom.to_string(), vcf_row.vep_start, vcf_row.vep_end);
                    local_buf.entry(key).or_default().push(e.clone());
                });
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
