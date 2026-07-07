//! Reassembles window-worker outputs into window-index order.
//!
//! Workers complete out of order (work-stealing). Each completed window is
//! recorded with its index; the drain releases batches only once every
//! lower-indexed window has been released, so output order is deterministic
//! regardless of completion order.
//
// TODO(parallel-window-driver): consumed by the parallel window driver added in
// a later phase; remove this `allow(dead_code)` once the driver wires it in.
#![allow(dead_code)]

use std::collections::HashMap;

use datafusion::arrow::record_batch::RecordBatch;

pub(crate) struct OrderedWindowDrain {
    next_emit: usize,
    pending: HashMap<usize, Vec<RecordBatch>>,
}

impl OrderedWindowDrain {
    pub(crate) fn new() -> Self {
        Self {
            next_emit: 0,
            pending: HashMap::new(),
        }
    }

    /// Record a completed window. Returns the batches now releasable in order
    /// (this window if it is `next_emit`, plus any contiguous buffered
    /// successors). Returns empty if this window must wait for an earlier one.
    pub(crate) fn complete(&mut self, index: usize, batches: Vec<RecordBatch>) -> Vec<RecordBatch> {
        self.pending.insert(index, batches);
        let mut released = Vec::new();
        while let Some(batches) = self.pending.remove(&self.next_emit) {
            released.extend(batches);
            self.next_emit += 1;
        }
        released
    }

    /// Number of completed-but-not-yet-released windows held in the buffer.
    pub(crate) fn buffered_windows(&self) -> usize {
        self.pending.len()
    }
}

impl Default for OrderedWindowDrain {
    fn default() -> Self {
        Self::new()
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use datafusion::arrow::array::Int64Array;
    use datafusion::arrow::datatypes::{DataType, Field, Schema};
    use std::sync::Arc;

    fn batch(tag: i64) -> RecordBatch {
        let schema = Arc::new(Schema::new(vec![Field::new("t", DataType::Int64, false)]));
        RecordBatch::try_new(schema, vec![Arc::new(Int64Array::from(vec![tag]))]).unwrap()
    }

    fn tags(batches: &[RecordBatch]) -> Vec<i64> {
        batches
            .iter()
            .map(|b| {
                b.column(0)
                    .as_any()
                    .downcast_ref::<Int64Array>()
                    .unwrap()
                    .value(0)
            })
            .collect()
    }

    #[test]
    fn releases_in_order_when_completed_in_order() {
        let mut drain = OrderedWindowDrain::new();
        assert_eq!(tags(&drain.complete(0, vec![batch(0)])), vec![0]);
        assert_eq!(tags(&drain.complete(1, vec![batch(1)])), vec![1]);
        assert_eq!(drain.buffered_windows(), 0);
    }

    #[test]
    fn buffers_out_of_order_then_releases_contiguously() {
        let mut drain = OrderedWindowDrain::new();
        // Window 2 finishes first -> nothing released, buffered.
        assert!(drain.complete(2, vec![batch(2)]).is_empty());
        // Window 1 finishes next -> still waiting on 0.
        assert!(drain.complete(1, vec![batch(1)]).is_empty());
        assert_eq!(drain.buffered_windows(), 2);
        // Window 0 finishes -> releases 0, 1, 2 in order.
        assert_eq!(tags(&drain.complete(0, vec![batch(0)])), vec![0, 1, 2]);
        assert_eq!(drain.buffered_windows(), 0);
    }
}
