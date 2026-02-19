use std::pin::Pin;
use std::sync::Arc;
use std::task::{Context, Poll};

use ahash::AHashMap;
use datafusion::arrow::array::RecordBatch;
use datafusion::arrow::compute::concat_batches;
use datafusion::arrow::datatypes::SchemaRef;
use datafusion::common::Result;
use datafusion::execution::SendableRecordBatchStream;
use futures::{StreamExt, ready};

use crate::array_utils::get_join_col_arrays;

/// Per-contig group with indexed intervals: `(contig_name, sorted_intervals)`.
/// Each interval carries `(start, end, global_row_idx)`.
pub type IndexedGroups = Vec<(String, Vec<(i64, i64, usize)>)>;

/// Collects input batches into per-contig interval groups, sorting each group
/// by (start, end) once the input is exhausted.
///
/// Uses `AHashMap<String, usize>` (contig→index) + `Vec<Vec<(i64, i64)>>`
/// for O(1) amortized lookups instead of O(log n) BTreeMap string comparisons.
pub struct StreamCollector {
    input: SendableRecordBatchStream,
    columns: Arc<(String, String, String)>,
    contig_to_idx: AHashMap<String, usize>,
    groups: Vec<Vec<(i64, i64)>>,
    contig_names: Vec<String>,
    done: bool,
}

impl StreamCollector {
    pub fn new(input: SendableRecordBatchStream, columns: Arc<(String, String, String)>) -> Self {
        Self {
            input,
            columns,
            contig_to_idx: AHashMap::new(),
            groups: Vec::new(),
            contig_names: Vec::new(),
            done: false,
        }
    }

    /// Poll the input stream, collecting rows into grouped intervals.
    ///
    /// Returns `Poll::Ready(Ok(true))` when all input has been consumed and
    /// groups are sorted. Returns `Poll::Ready(Ok(false))` is never returned —
    /// the method either returns `Pending` or `Ready(Ok(true))`.
    pub fn poll_collect(&mut self, cx: &mut Context<'_>) -> Poll<Result<bool>> {
        if self.done {
            return Poll::Ready(Ok(true));
        }

        loop {
            let batch_opt = ready!(Pin::new(&mut self.input).poll_next_unpin(cx));

            match batch_opt {
                Some(Ok(batch)) => {
                    if batch.num_rows() == 0 {
                        continue;
                    }
                    let (contig_arr, start_arr, end_arr) = get_join_col_arrays(
                        &batch,
                        (&self.columns.0, &self.columns.1, &self.columns.2),
                    )?;
                    let start_resolved = start_arr.resolve_i64()?;
                    let end_resolved = end_arr.resolve_i64()?;
                    let starts = &*start_resolved;
                    let ends = &*end_resolved;

                    for i in 0..batch.num_rows() {
                        let contig = contig_arr.value(i);
                        let idx = match self.contig_to_idx.get(contig) {
                            Some(&idx) => idx,
                            None => {
                                let idx = self.groups.len();
                                self.groups.push(Vec::with_capacity(1024));
                                self.contig_names.push(contig.to_string());
                                self.contig_to_idx.insert(contig.to_string(), idx);
                                idx
                            }
                        };
                        self.groups[idx].push((starts[i], ends[i]));
                    }
                }
                Some(Err(e)) => {
                    self.done = true;
                    return Poll::Ready(Err(e));
                }
                None => {
                    self.done = true;
                    for intervals in &mut self.groups {
                        intervals.sort_unstable();
                    }
                    return Poll::Ready(Ok(true));
                }
            }
        }
    }

    /// Take ownership of the collected groups as a `Vec<(contig, intervals)>`,
    /// sorted by contig name. Used by merge, cluster, complement, subtract
    /// when iterating groups in deterministic order.
    pub fn take_groups(&mut self) -> Vec<(String, Vec<(i64, i64)>)> {
        let names = std::mem::take(&mut self.contig_names);
        let groups = std::mem::take(&mut self.groups);
        self.contig_to_idx.clear();

        let mut pairs: Vec<(String, Vec<(i64, i64)>)> = names.into_iter().zip(groups).collect();
        pairs.sort_unstable_by(|a, b| a.0.cmp(&b.0));
        pairs
    }

    /// Take ownership of the collected groups as an `AHashMap<String, Vec<(i64, i64)>>`.
    /// Used by callers that need O(1) key lookup (complement view_bounds, subtract right_groups).
    pub fn take_groups_as_map(&mut self) -> AHashMap<String, Vec<(i64, i64)>> {
        let names = std::mem::take(&mut self.contig_names);
        let groups = std::mem::take(&mut self.groups);
        self.contig_to_idx.clear();

        names.into_iter().zip(groups).collect()
    }
}

/// Like [`StreamCollector`] but stores full `RecordBatch` rows so that extra
/// columns beyond (contig, start, end) can be reconstructed via Arrow `take`.
///
/// Groups carry `(start, end, global_row_idx)` tuples.
pub struct FullBatchCollector {
    input: SendableRecordBatchStream,
    columns: Arc<(String, String, String)>,
    input_schema: SchemaRef,
    batches: Vec<RecordBatch>,
    contig_to_idx: AHashMap<String, usize>,
    groups: Vec<Vec<(i64, i64, usize)>>,
    contig_names: Vec<String>,
    row_offset: usize,
    done: bool,
    concatenated: Option<RecordBatch>,
}

impl FullBatchCollector {
    pub fn new(
        input: SendableRecordBatchStream,
        columns: Arc<(String, String, String)>,
        input_schema: SchemaRef,
    ) -> Self {
        Self {
            input,
            columns,
            input_schema,
            batches: Vec::new(),
            contig_to_idx: AHashMap::new(),
            groups: Vec::new(),
            contig_names: Vec::new(),
            row_offset: 0,
            done: false,
            concatenated: None,
        }
    }

    /// Poll the input stream, collecting full batches and building grouped
    /// interval tuples `(start, end, global_row_idx)`.
    pub fn poll_collect(&mut self, cx: &mut Context<'_>) -> Poll<Result<bool>> {
        if self.done {
            return Poll::Ready(Ok(true));
        }

        loop {
            let batch_opt = ready!(Pin::new(&mut self.input).poll_next_unpin(cx));

            match batch_opt {
                Some(Ok(batch)) => {
                    if batch.num_rows() == 0 {
                        continue;
                    }
                    let (contig_arr, start_arr, end_arr) = get_join_col_arrays(
                        &batch,
                        (&self.columns.0, &self.columns.1, &self.columns.2),
                    )?;
                    let start_resolved = start_arr.resolve_i64()?;
                    let end_resolved = end_arr.resolve_i64()?;
                    let starts = &*start_resolved;
                    let ends = &*end_resolved;

                    let base = self.row_offset;
                    for i in 0..batch.num_rows() {
                        let contig = contig_arr.value(i);
                        let idx = match self.contig_to_idx.get(contig) {
                            Some(&idx) => idx,
                            None => {
                                let idx = self.groups.len();
                                self.groups.push(Vec::with_capacity(1024));
                                self.contig_names.push(contig.to_string());
                                self.contig_to_idx.insert(contig.to_string(), idx);
                                idx
                            }
                        };
                        self.groups[idx].push((starts[i], ends[i], base + i));
                    }

                    self.row_offset += batch.num_rows();
                    self.batches.push(batch);
                }
                Some(Err(e)) => {
                    self.done = true;
                    return Poll::Ready(Err(e));
                }
                None => {
                    self.done = true;
                    // Sort each group by (start, end, row_idx)
                    for intervals in &mut self.groups {
                        intervals.sort_unstable();
                    }
                    // Concatenate all stored batches into one
                    if !self.batches.is_empty() {
                        let concatenated = concat_batches(&self.input_schema, &self.batches)?;
                        self.batches.clear();
                        self.concatenated = Some(concatenated);
                    }
                    return Poll::Ready(Ok(true));
                }
            }
        }
    }

    /// Take ownership of the collected groups sorted by contig name.
    pub fn take_groups(&mut self) -> IndexedGroups {
        let names = std::mem::take(&mut self.contig_names);
        let groups = std::mem::take(&mut self.groups);
        self.contig_to_idx.clear();

        let mut pairs: IndexedGroups = names.into_iter().zip(groups).collect();
        pairs.sort_unstable_by(|a, b| a.0.cmp(&b.0));
        pairs
    }

    /// Take the concatenated `RecordBatch` containing all input rows.
    pub fn take_concatenated(&mut self) -> Option<RecordBatch> {
        self.concatenated.take()
    }
}
