use std::collections::HashMap;
use std::sync::Arc;

use datafusion::arrow::array::{Int16Builder, Int32Builder, RecordBatch, StringBuilder};
use datafusion::arrow::datatypes::SchemaRef;

use crate::events::{ContigEvents, DenseContigDepth};
use crate::schema::coverage_output_schema;

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
) -> datafusion::common::Result<RecordBatch> {
    let schema = coverage_output_schema();
    let mut all_blocks = Vec::new();

    // Sort contigs for deterministic output
    let mut contigs: Vec<String> = contig_events.keys().cloned().collect();
    contigs.sort();

    for contig in &contigs {
        let events = &mut contig_events.get_mut(contig).unwrap().events;
        let blocks = events_to_coverage_blocks(contig, events);
        all_blocks.extend(blocks);
    }

    coverage_blocks_to_record_batch(&all_blocks, &schema)
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
        let schema = coverage_output_schema();
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
}
