use std::collections::BTreeMap;
use std::pin::Pin;
use std::sync::Arc;
use std::task::{Context, Poll};

use datafusion::common::Result;
use datafusion::execution::SendableRecordBatchStream;
use futures::{StreamExt, ready};

use crate::array_utils::get_join_col_arrays;

/// Intervals grouped by contig, sorted by (start, end) within each group.
pub type GroupedIntervals = BTreeMap<String, Vec<(i64, i64)>>;

/// Default output batch size for collect-then-emit streams.
/// Larger than DataFusion's default 8192 to reduce flush overhead.
pub const DEFAULT_BATCH_SIZE: usize = 65_536;

/// Collects input batches into [`GroupedIntervals`], grouping by contig
/// and sorting each group by (start, end) once the input is exhausted.
pub struct StreamCollector {
    input: SendableRecordBatchStream,
    columns: Arc<(String, String, String)>,
    groups: GroupedIntervals,
    done: bool,
}

impl StreamCollector {
    pub fn new(input: SendableRecordBatchStream, columns: Arc<(String, String, String)>) -> Self {
        Self {
            input,
            columns,
            groups: BTreeMap::new(),
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
                        if let Some(vec) = self.groups.get_mut(contig) {
                            vec.push((starts[i], ends[i]));
                        } else {
                            self.groups
                                .insert(contig.to_string(), vec![(starts[i], ends[i])]);
                        }
                    }
                }
                Some(Err(e)) => {
                    self.done = true;
                    return Poll::Ready(Err(e));
                }
                None => {
                    self.done = true;
                    for intervals in self.groups.values_mut() {
                        intervals.sort_unstable();
                    }
                    return Poll::Ready(Ok(true));
                }
            }
        }
    }

    /// Take ownership of the collected groups. Returns an empty map if called
    /// before collection is complete.
    pub fn take_groups(&mut self) -> GroupedIntervals {
        std::mem::take(&mut self.groups)
    }
}
