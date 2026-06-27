//! Groups a contig's position-sorted input batches into indexed windows.
//!
//! A window is a contiguous run of whole input batches totalling at least
//! `target_variants` rows (no mid-batch slicing). The final window may be
//! smaller. Windows are the unit of work handed to the parallel window pool.
//
// TODO(parallel-window-driver): the planner is consumed by the parallel window
// driver added in a later phase; until then the items are crate-internal and
// unused outside tests. Remove this `allow(dead_code)` once the driver wires
// it in.
#![allow(dead_code)]

use datafusion::arrow::record_batch::RecordBatch;

/// A contiguous group of position-sorted input batches for one contig.
pub(crate) struct Window {
    pub index: usize,
    pub batches: Vec<RecordBatch>,
}

/// Accumulates input batches and emits a [`Window`] each time the pending
/// group reaches `target_variants` rows.
pub(crate) struct WindowPlanner {
    target_variants: usize,
    pending: Vec<RecordBatch>,
    pending_rows: usize,
    next_index: usize,
}

impl WindowPlanner {
    pub(crate) fn new(target_variants: usize) -> Self {
        Self {
            target_variants: target_variants.max(1),
            pending: Vec::new(),
            pending_rows: 0,
            next_index: 0,
        }
    }

    /// Push one input batch. Returns a [`Window`] when the pending group
    /// reaches `target_variants`. Empty batches are ignored.
    pub(crate) fn push(&mut self, batch: RecordBatch) -> Option<Window> {
        let rows = batch.num_rows();
        if rows == 0 {
            return None;
        }
        self.pending.push(batch);
        self.pending_rows += rows;
        if self.pending_rows >= self.target_variants {
            Some(self.flush())
        } else {
            None
        }
    }

    /// Emit any remaining pending batches as a final (possibly smaller) window.
    pub(crate) fn finish(&mut self) -> Option<Window> {
        if self.pending.is_empty() {
            None
        } else {
            Some(self.flush())
        }
    }

    fn flush(&mut self) -> Window {
        let index = self.next_index;
        self.next_index += 1;
        let batches = std::mem::take(&mut self.pending);
        self.pending_rows = 0;
        Window { index, batches }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use datafusion::arrow::array::Int64Array;
    use datafusion::arrow::datatypes::{DataType, Field, Schema};
    use std::sync::Arc;

    fn batch(rows: usize) -> RecordBatch {
        let schema = Arc::new(Schema::new(vec![Field::new("pos", DataType::Int64, false)]));
        let arr = Int64Array::from((0..rows as i64).collect::<Vec<_>>());
        RecordBatch::try_new(schema, vec![Arc::new(arr)]).unwrap()
    }

    #[test]
    fn groups_whole_batches_to_target_and_flushes_remainder() {
        let mut planner = WindowPlanner::new(5000);
        // 3000 pending -> below target, no window yet.
        assert!(planner.push(batch(3000)).is_none());
        // +3000 = 6000 >= 5000 -> emit window 0 with both batches.
        let w0 = planner.push(batch(3000)).expect("window 0");
        assert_eq!(w0.index, 0);
        assert_eq!(w0.batches.iter().map(|b| b.num_rows()).sum::<usize>(), 6000);
        // Remaining 3000 emitted by finish() as window 1.
        assert!(planner.push(batch(3000)).is_none());
        let w1 = planner.finish().expect("window 1");
        assert_eq!(w1.index, 1);
        assert_eq!(w1.batches.iter().map(|b| b.num_rows()).sum::<usize>(), 3000);
        // Nothing left.
        assert!(planner.finish().is_none());
    }

    #[test]
    fn ignores_empty_batches() {
        let mut planner = WindowPlanner::new(10);
        assert!(planner.push(batch(0)).is_none());
        assert!(planner.finish().is_none());
    }
}
