use std::collections::HashMap;
use std::sync::Arc;

use datafusion::arrow::array::{Int16Builder, Int32Builder, RecordBatch, StringBuilder};
use datafusion::arrow::datatypes::SchemaRef;
use datafusion::common::Result;

use crate::events::{ContigEvents, DenseContigDepth};

/// A single coverage block.
#[derive(Debug, Clone, PartialEq)]
pub struct CoverageBlock {
    pub contig: Arc<str>,
    pub pos_start: i32,
    pub pos_end: i32,
    pub coverage: i16,
}

/// Convert per-contig events (Vec) to coverage blocks.
///
/// Sorts the events in-place, merges same-position deltas, then walks with
/// cumulative sum. For coordinate-sorted BAM input the Vec is nearly sorted,
/// so `sort_unstable` (timsort) runs in ~O(n).
/// Zero-coverage gaps are NOT emitted.
pub fn events_to_coverage_blocks(contig: &str, events: &mut [(u32, i32)]) -> Vec<CoverageBlock> {
    if events.is_empty() {
        return Vec::new();
    }

    events.sort_unstable_by_key(|e| e.0);

    let contig_arc: Arc<str> = Arc::from(contig);
    let mut blocks = Vec::new();
    let mut cov: i32 = 0;
    let mut prev_cov: i32 = 0;
    let mut block_start: u32 = 0;

    let mut i = 0;
    while i < events.len() {
        let pos = events[i].0;
        // Merge all deltas at the same position
        while i < events.len() && events[i].0 == pos {
            cov += events[i].1;
            i += 1;
        }

        if prev_cov != 0 && cov != prev_cov {
            blocks.push(CoverageBlock {
                contig: Arc::clone(&contig_arc),
                pos_start: block_start as i32,
                pos_end: (pos - 1) as i32,
                coverage: prev_cov as i16,
            });
            if cov != 0 {
                block_start = pos;
            }
        } else if prev_cov == 0 && cov != 0 {
            block_start = pos;
        }

        prev_cov = cov;
    }

    blocks
}

/// Convert all contig events to a single RecordBatch of coverage blocks.
pub fn all_events_to_record_batch(
    contig_events: &mut HashMap<String, ContigEvents>,
    schema: &SchemaRef,
) -> datafusion::common::Result<RecordBatch> {
    let mut all_blocks = Vec::new();

    // Sort contigs for deterministic output
    let mut contigs: Vec<String> = contig_events.keys().cloned().collect();
    contigs.sort();

    for contig in &contigs {
        let events = &mut contig_events.get_mut(contig).unwrap().events;
        let blocks = events_to_coverage_blocks(contig, events);
        all_blocks.extend(blocks);
    }

    coverage_blocks_to_record_batch(&all_blocks, schema)
}

/// Convert a dense depth array to coverage blocks via prefix sum + RLE.
///
/// Walks the array linearly, maintaining a running coverage total.
/// The key optimization: most positions have `delta == 0` and are skipped
/// (~99.9% of positions for typical WES/WGS data).
/// Zero-coverage gaps are NOT emitted.
pub fn dense_depth_to_coverage_blocks(contig: &str, depth: &[i32]) -> Vec<CoverageBlock> {
    let contig_arc: Arc<str> = Arc::from(contig);
    let mut blocks = Vec::new();
    let mut cov: i32 = 0;
    let mut prev_cov: i32 = 0;
    let mut block_start: usize = 0;

    for (pos, &delta) in depth.iter().enumerate() {
        if delta == 0 {
            continue;
        }
        cov += delta;

        if prev_cov != 0 && cov != prev_cov {
            blocks.push(CoverageBlock {
                contig: Arc::clone(&contig_arc),
                pos_start: block_start as i32,
                pos_end: (pos - 1) as i32,
                coverage: prev_cov as i16,
            });
            if cov != 0 {
                block_start = pos;
            }
        } else if prev_cov == 0 && cov != 0 {
            block_start = pos;
        }

        prev_cov = cov;
    }

    blocks
}

/// Convert a dense depth slice to coverage blocks with a position offset.
///
/// Same RLE logic as `dense_depth_to_coverage_blocks` but positions are
/// shifted by `start_offset`, allowing callers to pass a sub-slice of
/// the full depth array (the touched region) instead of the entire contig.
pub fn dense_depth_to_coverage_blocks_bounded(
    contig: &str,
    depth_slice: &[i32],
    start_offset: usize,
) -> Vec<CoverageBlock> {
    let contig_arc: Arc<str> = Arc::from(contig);
    let mut blocks = Vec::new();
    let mut cov: i32 = 0;
    let mut prev_cov: i32 = 0;
    let mut block_start: usize = 0;

    for (i, &delta) in depth_slice.iter().enumerate() {
        if delta == 0 {
            continue;
        }
        let pos = start_offset + i;
        cov += delta;

        if prev_cov != 0 && cov != prev_cov {
            blocks.push(CoverageBlock {
                contig: Arc::clone(&contig_arc),
                pos_start: block_start as i32,
                pos_end: (pos - 1) as i32,
                coverage: prev_cov as i16,
            });
            if cov != 0 {
                block_start = pos;
            }
        } else if prev_cov == 0 && cov != 0 {
            block_start = pos;
        }

        prev_cov = cov;
    }

    blocks
}

/// Convert a single contig's dense depth to a RecordBatch.
///
/// Uses the tracked touched-range bounds from `DenseContigDepth` to scan
/// only the modified region instead of the full contig-length array.
pub fn dense_depth_to_record_batch(
    contig: &str,
    depth: &DenseContigDepth,
    schema: &SchemaRef,
) -> datafusion::common::Result<RecordBatch> {
    let blocks = dense_depth_to_coverage_blocks_bounded(
        contig,
        depth.touched_range(),
        depth.touched_start(),
    );
    coverage_blocks_to_record_batch(&blocks, schema)
}

/// Convert a single contig's dense depth to multiple RecordBatches, each at most `batch_size` rows.
pub fn dense_depth_to_record_batches(
    contig: &str,
    depth: &DenseContigDepth,
    schema: &SchemaRef,
    batch_size: usize,
) -> datafusion::common::Result<Vec<RecordBatch>> {
    let blocks = dense_depth_to_coverage_blocks_bounded(
        contig,
        depth.touched_range(),
        depth.touched_start(),
    );
    coverage_blocks_to_chunked_batches(&blocks, schema, batch_size)
}

/// Convert all contig events to multiple RecordBatches, each at most `batch_size` rows.
pub fn all_events_to_record_batches(
    contig_events: &mut HashMap<String, ContigEvents>,
    schema: &SchemaRef,
    batch_size: usize,
) -> datafusion::common::Result<Vec<RecordBatch>> {
    let mut all_blocks = Vec::new();

    let mut contigs: Vec<String> = contig_events.keys().cloned().collect();
    contigs.sort();

    for contig in &contigs {
        let events = &mut contig_events.get_mut(contig).unwrap().events;
        let blocks = events_to_coverage_blocks(contig, events);
        all_blocks.extend(blocks);
    }

    coverage_blocks_to_chunked_batches(&all_blocks, schema, batch_size)
}

/// Split coverage blocks into multiple RecordBatches of at most `batch_size` rows each.
pub fn coverage_blocks_to_chunked_batches(
    blocks: &[CoverageBlock],
    schema: &SchemaRef,
    batch_size: usize,
) -> datafusion::common::Result<Vec<RecordBatch>> {
    if blocks.is_empty() {
        return Ok(Vec::new());
    }
    let mut batches = Vec::with_capacity(blocks.len().div_ceil(batch_size));
    for chunk in blocks.chunks(batch_size) {
        let batch = coverage_blocks_to_record_batch(chunk, schema)?;
        batches.push(batch);
    }
    Ok(batches)
}

/// Convert a slice of coverage blocks to an Arrow RecordBatch.
///
/// Uses Arrow builders for single-pass construction without intermediate Vecs.
pub fn coverage_blocks_to_record_batch(
    blocks: &[CoverageBlock],
    schema: &SchemaRef,
) -> datafusion::common::Result<RecordBatch> {
    let n = blocks.len();
    let mut contig_builder = StringBuilder::with_capacity(n, n * 5);
    let mut start_builder = Int32Builder::with_capacity(n);
    let mut end_builder = Int32Builder::with_capacity(n);
    let mut cov_builder = Int16Builder::with_capacity(n);

    for b in blocks {
        contig_builder.append_value(b.contig.as_ref());
        start_builder.append_value(b.pos_start);
        end_builder.append_value(b.pos_end);
        cov_builder.append_value(b.coverage);
    }

    let batch = RecordBatch::try_new(
        schema.clone(),
        vec![
            Arc::new(contig_builder.finish()),
            Arc::new(start_builder.finish()),
            Arc::new(end_builder.finish()),
            Arc::new(cov_builder.finish()),
        ],
    )?;

    Ok(batch)
}

/// Lazy streaming emitter for per-base coverage output.
///
/// Holds the finalized dense depth delta array for a single contig and
/// generates one `RecordBatch` per `next_batch()` call, keeping memory
/// at O(batch_size) instead of materializing all positions at once.
pub struct PerBaseEmitter {
    contig: String,
    /// Delta array taken from `DenseContigDepth`. Prefix-summed lazily.
    depth: Vec<i32>,
    /// Next position index to emit.
    next_idx: usize,
    /// Exclusive end index.
    end_idx: usize,
    /// Running prefix sum carried across batches.
    cumulative_cov: i32,
}

impl PerBaseEmitter {
    /// Create a new emitter for the given contig.
    ///
    /// - `zero_based=true`: emits positions `[0..len)` (0-based)
    /// - `zero_based=false`: emits positions `[1..len]` (1-based), but the
    ///   underlying array index still starts at 1 (BAM 1-based convention).
    pub fn new(contig: String, depth_vec: Vec<i32>, zero_based: bool) -> Self {
        let len = depth_vec.len();
        let (next_idx, end_idx) = if zero_based {
            (0, len)
        } else {
            // 1-based: skip index 0, emit [1..len)
            (1, len)
        };
        Self {
            contig,
            depth: depth_vec,
            next_idx,
            end_idx,
            cumulative_cov: 0,
        }
    }

    /// Generate the next batch of up to `batch_size` per-base rows.
    ///
    /// Returns `None` when all positions have been emitted.
    pub fn next_batch(
        &mut self,
        schema: &SchemaRef,
        batch_size: usize,
    ) -> Option<Result<RecordBatch>> {
        if self.next_idx >= self.end_idx {
            return None;
        }

        let chunk_end = (self.next_idx + batch_size).min(self.end_idx);
        let n = chunk_end - self.next_idx;

        let mut contig_builder = StringBuilder::with_capacity(n, n * self.contig.len());
        let mut pos_builder = Int32Builder::with_capacity(n);
        let mut cov_builder = Int16Builder::with_capacity(n);

        for idx in self.next_idx..chunk_end {
            self.cumulative_cov += self.depth[idx];
            contig_builder.append_value(&self.contig);
            pos_builder.append_value(idx as i32);
            cov_builder.append_value(self.cumulative_cov as i16);
        }

        self.next_idx = chunk_end;

        let batch = RecordBatch::try_new(
            schema.clone(),
            vec![
                Arc::new(contig_builder.finish()),
                Arc::new(pos_builder.finish()),
                Arc::new(cov_builder.finish()),
            ],
        );
        Some(batch.map_err(Into::into))
    }

    /// Eagerly drain all remaining positions into batches.
    ///
    /// Used when flushing a previous contig's emitter during a multi-contig
    /// transition (rare edge case).
    pub fn flush_remaining(
        &mut self,
        schema: &SchemaRef,
        batch_size: usize,
    ) -> Result<Vec<RecordBatch>> {
        let mut batches = Vec::new();
        while let Some(result) = self.next_batch(schema, batch_size) {
            batches.push(result?);
        }
        Ok(batches)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use datafusion::arrow::array::{Int16Array, Int32Array, StringArray};

    #[test]
    fn test_single_read() {
        let mut events = vec![(100u32, 1i32), (110, -1)];

        let blocks = events_to_coverage_blocks("chr1", &mut events);
        assert_eq!(
            blocks,
            vec![CoverageBlock {
                contig: Arc::from("chr1"),
                pos_start: 100,
                pos_end: 109,
                coverage: 1,
            }]
        );
    }

    #[test]
    fn test_overlapping_reads() {
        // Read1: [0, 10), Read2: [5, 15)
        let mut events = vec![(0u32, 1i32), (10, -1), (5, 1), (15, -1)];

        let blocks = events_to_coverage_blocks("chr1", &mut events);
        assert_eq!(
            blocks,
            vec![
                CoverageBlock {
                    contig: Arc::from("chr1"),
                    pos_start: 0,
                    pos_end: 4,
                    coverage: 1,
                },
                CoverageBlock {
                    contig: Arc::from("chr1"),
                    pos_start: 5,
                    pos_end: 9,
                    coverage: 2,
                },
                CoverageBlock {
                    contig: Arc::from("chr1"),
                    pos_start: 10,
                    pos_end: 14,
                    coverage: 1,
                },
            ]
        );
    }

    #[test]
    fn test_gap_between_reads() {
        // Read1: [0, 5), Read2: [10, 15)
        let mut events = vec![(0u32, 1i32), (5, -1), (10, 1), (15, -1)];

        let blocks = events_to_coverage_blocks("chr1", &mut events);
        assert_eq!(
            blocks,
            vec![
                CoverageBlock {
                    contig: Arc::from("chr1"),
                    pos_start: 0,
                    pos_end: 4,
                    coverage: 1,
                },
                CoverageBlock {
                    contig: Arc::from("chr1"),
                    pos_start: 10,
                    pos_end: 14,
                    coverage: 1,
                },
            ]
        );
    }

    #[test]
    fn test_deletion_gap() {
        // 5M3D5M at position 0:
        // Events: (0, +1), (5, -1), (8, +1), (13, -1)
        let mut events = vec![(0u32, 1i32), (5, -1), (8, 1), (13, -1)];

        let blocks = events_to_coverage_blocks("chr1", &mut events);
        assert_eq!(
            blocks,
            vec![
                CoverageBlock {
                    contig: Arc::from("chr1"),
                    pos_start: 0,
                    pos_end: 4,
                    coverage: 1,
                },
                CoverageBlock {
                    contig: Arc::from("chr1"),
                    pos_start: 8,
                    pos_end: 12,
                    coverage: 1,
                },
            ]
        );
    }

    #[test]
    fn test_empty_input() {
        let mut events = Vec::new();
        let blocks = events_to_coverage_blocks("chr1", &mut events);
        assert!(blocks.is_empty());
    }

    #[test]
    fn test_record_batch_output() {
        let blocks = vec![CoverageBlock {
            contig: Arc::from("chr1"),
            pos_start: 0,
            pos_end: 9,
            coverage: 1,
        }];
        let schema = crate::schema::coverage_output_schema(true);
        let batch = coverage_blocks_to_record_batch(&blocks, &schema).unwrap();

        assert_eq!(batch.num_rows(), 1);
        assert_eq!(
            batch
                .column(0)
                .as_any()
                .downcast_ref::<StringArray>()
                .unwrap()
                .value(0),
            "chr1"
        );
        assert_eq!(
            batch
                .column(1)
                .as_any()
                .downcast_ref::<Int32Array>()
                .unwrap()
                .value(0),
            0
        );
        assert_eq!(
            batch
                .column(2)
                .as_any()
                .downcast_ref::<Int32Array>()
                .unwrap()
                .value(0),
            9
        );
        assert_eq!(
            batch
                .column(3)
                .as_any()
                .downcast_ref::<Int16Array>()
                .unwrap()
                .value(0),
            1
        );
    }

    // --- Tests for dense_depth_to_coverage_blocks ---

    #[test]
    fn test_dense_single_read() {
        // Simulate one 10bp read at position 100
        let mut depth = vec![0i32; 120];
        depth[100] = 1;
        depth[110] = -1;

        let blocks = dense_depth_to_coverage_blocks("chr1", &depth);
        assert_eq!(
            blocks,
            vec![CoverageBlock {
                contig: Arc::from("chr1"),
                pos_start: 100,
                pos_end: 109,
                coverage: 1,
            }]
        );
    }

    #[test]
    fn test_dense_overlapping_reads() {
        // Read1: [0, 10), Read2: [5, 15)
        let mut depth = vec![0i32; 20];
        depth[0] = 1;
        depth[5] += 1;
        depth[10] += -1;
        depth[15] = -1;

        let blocks = dense_depth_to_coverage_blocks("chr1", &depth);
        assert_eq!(
            blocks,
            vec![
                CoverageBlock {
                    contig: Arc::from("chr1"),
                    pos_start: 0,
                    pos_end: 4,
                    coverage: 1,
                },
                CoverageBlock {
                    contig: Arc::from("chr1"),
                    pos_start: 5,
                    pos_end: 9,
                    coverage: 2,
                },
                CoverageBlock {
                    contig: Arc::from("chr1"),
                    pos_start: 10,
                    pos_end: 14,
                    coverage: 1,
                },
            ]
        );
    }

    #[test]
    fn test_dense_gap_between_reads() {
        // Read1: [0, 5), Read2: [10, 15)
        let mut depth = vec![0i32; 20];
        depth[0] = 1;
        depth[5] = -1;
        depth[10] = 1;
        depth[15] = -1;

        let blocks = dense_depth_to_coverage_blocks("chr1", &depth);
        assert_eq!(
            blocks,
            vec![
                CoverageBlock {
                    contig: Arc::from("chr1"),
                    pos_start: 0,
                    pos_end: 4,
                    coverage: 1,
                },
                CoverageBlock {
                    contig: Arc::from("chr1"),
                    pos_start: 10,
                    pos_end: 14,
                    coverage: 1,
                },
            ]
        );
    }

    #[test]
    fn test_dense_empty() {
        let depth = vec![0i32; 100];
        let blocks = dense_depth_to_coverage_blocks("chr1", &depth);
        assert!(blocks.is_empty());
    }

    #[test]
    fn test_dense_matches_sparse() {
        // Verify dense and sparse produce identical results
        let mut events = vec![(0u32, 1i32), (5, 1), (10, -1), (15, -1)];
        let sparse_blocks = events_to_coverage_blocks("chr1", &mut events);

        let mut depth = vec![0i32; 20];
        depth[0] = 1;
        depth[5] = 1;
        depth[10] = -1;
        depth[15] = -1;
        let dense_blocks = dense_depth_to_coverage_blocks("chr1", &depth);

        assert_eq!(sparse_blocks, dense_blocks);
    }

    // --- Tests for PerBaseEmitter ---

    #[test]
    fn test_per_base_emitter_single_read() {
        // 10-position contig with a single 5bp read at position 2
        let mut depth = vec![0i32; 10];
        depth[2] = 1;
        depth[7] = -1;

        let schema = crate::schema::per_base_output_schema(true);
        let mut emitter = PerBaseEmitter::new("chr1".to_string(), depth, true);

        // Emit all in one big batch
        let batch = emitter.next_batch(&schema, 1000).unwrap().unwrap();
        assert!(emitter.next_batch(&schema, 1000).is_none());

        assert_eq!(batch.num_rows(), 10);

        let positions = batch
            .column(1)
            .as_any()
            .downcast_ref::<Int32Array>()
            .unwrap();
        let coverages = batch
            .column(2)
            .as_any()
            .downcast_ref::<Int16Array>()
            .unwrap();

        // Positions 0,1 → coverage 0; positions 2-6 → coverage 1; positions 7-9 → coverage 0
        assert_eq!(coverages.value(0), 0); // pos 0
        assert_eq!(coverages.value(1), 0); // pos 1
        assert_eq!(coverages.value(2), 1); // pos 2
        assert_eq!(coverages.value(6), 1); // pos 6
        assert_eq!(coverages.value(7), 0); // pos 7
        assert_eq!(coverages.value(9), 0); // pos 9

        // Verify positions are 0-based sequential
        for i in 0..10 {
            assert_eq!(positions.value(i), i as i32);
        }
    }

    #[test]
    fn test_per_base_emitter_chunking() {
        // 10-position contig with a read at position 3
        let mut depth = vec![0i32; 10];
        depth[3] = 1;
        depth[6] = -1;

        let schema = crate::schema::per_base_output_schema(true);
        let mut emitter = PerBaseEmitter::new("chr1".to_string(), depth, true);

        // Emit in chunks of 4
        let batch1 = emitter.next_batch(&schema, 4).unwrap().unwrap();
        assert_eq!(batch1.num_rows(), 4); // positions 0-3

        let batch2 = emitter.next_batch(&schema, 4).unwrap().unwrap();
        assert_eq!(batch2.num_rows(), 4); // positions 4-7

        let batch3 = emitter.next_batch(&schema, 4).unwrap().unwrap();
        assert_eq!(batch3.num_rows(), 2); // positions 8-9

        assert!(emitter.next_batch(&schema, 4).is_none());

        // Verify cumulative_cov carries across batches:
        // batch1: pos 0→0, 1→0, 2→0, 3→1
        let cov1 = batch1
            .column(2)
            .as_any()
            .downcast_ref::<Int16Array>()
            .unwrap();
        assert_eq!(cov1.value(3), 1);

        // batch2: pos 4→1, 5→1, 6→0, 7→0 (delta -1 at position 6)
        let cov2 = batch2
            .column(2)
            .as_any()
            .downcast_ref::<Int16Array>()
            .unwrap();
        assert_eq!(cov2.value(0), 1); // pos 4
        assert_eq!(cov2.value(1), 1); // pos 5
        assert_eq!(cov2.value(2), 0); // pos 6
        assert_eq!(cov2.value(3), 0); // pos 7
    }

    #[test]
    fn test_per_base_emitter_one_based() {
        // 1-based: emits positions [1..len)
        let mut depth = vec![0i32; 6]; // indices 0..5
        depth[1] = 1; // 1-based position 1
        depth[4] = -1; // 1-based position 4

        let schema = crate::schema::per_base_output_schema(false);
        let mut emitter = PerBaseEmitter::new("chr1".to_string(), depth, false);

        let batch = emitter.next_batch(&schema, 1000).unwrap().unwrap();
        assert!(emitter.next_batch(&schema, 1000).is_none());

        // Should emit 5 rows (positions 1,2,3,4,5)
        assert_eq!(batch.num_rows(), 5);

        let positions = batch
            .column(1)
            .as_any()
            .downcast_ref::<Int32Array>()
            .unwrap();
        let coverages = batch
            .column(2)
            .as_any()
            .downcast_ref::<Int16Array>()
            .unwrap();

        assert_eq!(positions.value(0), 1);
        assert_eq!(positions.value(4), 5);
        assert_eq!(coverages.value(0), 1); // pos 1: delta +1
        assert_eq!(coverages.value(1), 1); // pos 2
        assert_eq!(coverages.value(2), 1); // pos 3
        assert_eq!(coverages.value(3), 0); // pos 4: delta -1
        assert_eq!(coverages.value(4), 0); // pos 5
    }

    #[test]
    fn test_per_base_emitter_flush_remaining() {
        let mut depth = vec![0i32; 10];
        depth[2] = 1;
        depth[7] = -1;

        let schema = crate::schema::per_base_output_schema(true);
        let mut emitter = PerBaseEmitter::new("chr1".to_string(), depth, true);

        let batches = emitter.flush_remaining(&schema, 3).unwrap();
        // 10 positions / batch_size 3 → 4 batches (3+3+3+1)
        assert_eq!(batches.len(), 4);
        let total_rows: usize = batches.iter().map(|b| b.num_rows()).sum();
        assert_eq!(total_rows, 10);
    }
}
