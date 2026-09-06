//! Genomic-region restriction for `annotate_vep()`: option parsing, run
//! planning over the input-buffer grid, and the row gate that keeps a range
//! cut exact on stateful (Merged/RefSeq) caches. Pure code; no session access.

use datafusion::common::{DataFusionError, Result};
use serde_json::Value;

/// One `regions` entry as given in `options_json`, 1-based closed.
#[derive(Clone, Debug, PartialEq, Eq)]
pub(crate) struct RegionSpec {
    pub chrom: String,
    pub start: Option<i64>,
    pub end: Option<i64>,
}

fn plan_err(msg: String) -> DataFusionError {
    DataFusionError::Plan(format!("annotate_vep(): regions {msg}"))
}

fn optional_coordinate(obj: &serde_json::Map<String, Value>, key: &str) -> Result<Option<i64>> {
    match obj.get(key) {
        None | Some(Value::Null) => Ok(None),
        Some(Value::Number(n)) => {
            let v = n
                .as_i64()
                .ok_or_else(|| plan_err(format!("'{key}' must be an integer")))?;
            if v < 1 {
                return Err(plan_err(format!("'{key}' must be >= 1, got {v}")));
            }
            Ok(Some(v))
        }
        Some(_) => Err(plan_err(format!("'{key}' must be an integer"))),
    }
}

/// Parse and validate the `regions` option. `Ok(None)` when absent.
pub(crate) fn parse_regions_option(options_json: Option<&str>) -> Result<Option<Vec<RegionSpec>>> {
    let Some(raw) = options_json else {
        return Ok(None);
    };
    let value: Value = serde_json::from_str(raw).map_err(|e| {
        DataFusionError::Plan(format!(
            "annotate_vep() options_json must be valid JSON: {e}"
        ))
    })?;
    let Some(entries) = value.get("regions") else {
        return Ok(None);
    };
    let entries = entries
        .as_array()
        .ok_or_else(|| plan_err("must be an array of {chrom, start?, end?} objects".into()))?;
    if entries.is_empty() {
        return Err(plan_err("must not be empty".into()));
    }
    let mut specs = Vec::with_capacity(entries.len());
    for entry in entries {
        let obj = entry
            .as_object()
            .ok_or_else(|| plan_err("entries must be objects".into()))?;
        let chrom = obj
            .get("chrom")
            .and_then(Value::as_str)
            .ok_or_else(|| plan_err("entries need a string 'chrom'".into()))?;
        if chrom.is_empty() || chrom.contains('`') {
            return Err(plan_err(format!("invalid chrom {chrom:?}")));
        }
        let start = optional_coordinate(obj, "start")?;
        let end = optional_coordinate(obj, "end")?;
        if let (Some(s), Some(e)) = (start, end) {
            if s > e {
                return Err(plan_err(format!("start {s} > end {e} on {chrom}")));
            }
        }
        specs.push(RegionSpec {
            chrom: chrom.to_string(),
            start,
            end,
        });
    }
    Ok(Some(specs))
}

/// A closed position interval on `start`; `None` sides are open.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub(crate) struct RunBounds {
    pub lo: Option<i64>,
    pub hi: Option<i64>,
}

impl RunBounds {
    pub const OPEN: RunBounds = RunBounds { lo: None, hi: None };

    pub fn is_open(&self) -> bool {
        self.lo.is_none() && self.hi.is_none()
    }

    pub fn contains(&self, pos: i64) -> bool {
        self.lo.is_none_or(|lo| pos >= lo) && self.hi.is_none_or(|hi| pos <= hi)
    }
}

/// Sort by lower bound and merge intervals that overlap or touch
/// (`hi + 1 >= next.lo`). Open sides absorb whatever they reach.
pub(crate) fn merge_bounds(mut bounds: Vec<RunBounds>) -> Vec<RunBounds> {
    bounds.sort_by_key(|b| (b.lo.unwrap_or(i64::MIN), b.hi.unwrap_or(i64::MAX)));
    let mut merged: Vec<RunBounds> = Vec::with_capacity(bounds.len());
    for b in bounds {
        match merged.last_mut() {
            Some(last)
                if last
                    .hi
                    .is_none_or(|hi| b.lo.is_none_or(|lo| lo <= hi.saturating_add(1))) =>
            {
                last.hi = match (last.hi, b.hi) {
                    (None, _) | (_, None) => None,
                    (Some(a), Some(c)) => Some(a.max(c)),
                };
            }
            _ => merged.push(b),
        }
    }
    merged
}

use std::collections::HashMap;

/// Per-contig merged intervals, keyed by the VCF's own contig spelling.
pub(crate) type ContigRuns = HashMap<String, Vec<RunBounds>>;

/// Map each spec onto the VCF contig it names (through the alias set, so
/// `1`/`chr1` and `M`/`MT`/`chrM`/`chrMT` match) and merge per contig. Specs
/// naming a contig the VCF does not have select nothing and are dropped.
pub(crate) fn resolve_regions(specs: &[RegionSpec], vcf_contigs: &[String]) -> ContigRuns {
    let mut by_alias: HashMap<String, &String> = HashMap::new();
    for contig in vcf_contigs {
        for alias in crate::cache::manifest::contig_alias_set(contig) {
            by_alias.entry(alias).or_insert(contig);
        }
    }
    let mut runs: HashMap<String, Vec<RunBounds>> = HashMap::new();
    for spec in specs {
        if let Some(contig) = by_alias.get(&spec.chrom) {
            runs.entry((*contig).clone()).or_default().push(RunBounds {
                lo: spec.start,
                hi: spec.end,
            });
        }
    }
    runs.into_iter()
        .map(|(k, v)| (k, merge_bounds(v)))
        .collect()
}

/// Keep only contigs that have a run, in the order given.
pub(crate) fn restrict_contigs(contigs: Vec<String>, runs: &ContigRuns) -> Vec<String> {
    contigs
        .into_iter()
        .filter(|c| runs.contains_key(c))
        .collect()
}

/// Map an interval onto whole input buffers. `boundary_positions[k]` is the
/// `start` of the first row of buffer `k`; the last entry (`i64::MAX`) is the
/// terminal boundary, so there are `len - 1` buffers. Returns `(bk, bk1)` with
/// `bk < bk1`: the first buffer whose start is at or below `lo` (buffer 0 when
/// `lo` precedes every row) through the last buffer whose start is at or below
/// `hi`, exclusive. Including buffer `bk` even when `lo` falls inside it is a
/// superset; the output trim removes the extra rows.
pub(crate) fn buffer_range_for_bounds(
    boundary_positions: &[i64],
    bounds: RunBounds,
) -> (usize, usize) {
    let b = boundary_positions.len().saturating_sub(1);
    if b == 0 {
        return (0, 0);
    }
    let bk = match bounds.lo {
        None => 0,
        // number of buffer starts <= lo, minus one; clamp at 0
        Some(lo) => boundary_positions[..b]
            .partition_point(|&p| p <= lo)
            .saturating_sub(1),
    };
    let bk1 = match bounds.hi {
        None => b,
        Some(hi) => boundary_positions[..b]
            .partition_point(|&p| p <= hi)
            .max(bk + 1)
            .min(b),
    };
    (bk, bk1)
}

/// One lookup pass over a contig: the original intervals it serves (for the
/// output trim) and, on the grid path, the whole-buffer range it must annotate.
#[derive(Clone, Debug, PartialEq, Eq)]
pub(crate) struct RunPlan {
    pub bounds: Vec<RunBounds>,
    pub buffers: Option<(usize, usize)>,
}

/// Plan the runs of one contig. Without a grid (`None`, the stateless Ensembl
/// path) every merged interval is its own run. With a grid, intervals map to
/// buffer ranges and ranges that overlap or touch merge into one run so no
/// buffer is warmed up and annotated twice.
pub(crate) fn plan_runs(bounds: &[RunBounds], boundary_positions: Option<&[i64]>) -> Vec<RunPlan> {
    let bounds = merge_bounds(bounds.to_vec());
    let Some(positions) = boundary_positions else {
        return bounds
            .into_iter()
            .map(|b| RunPlan {
                bounds: vec![b],
                buffers: None,
            })
            .collect();
    };
    if positions.len() < 2 {
        return Vec::new();
    }
    let mut runs: Vec<RunPlan> = Vec::new();
    for b in bounds {
        let (bk, bk1) = buffer_range_for_bounds(positions, b);
        match runs.last_mut() {
            Some(last) if last.buffers.is_some_and(|(_, end)| bk <= end) => {
                last.bounds.push(b);
                let (start, end) = last.buffers.unwrap_or((bk, bk1));
                last.buffers = Some((start, end.max(bk1)));
            }
            _ => runs.push(RunPlan {
                bounds: vec![b],
                buffers: Some((bk, bk1)),
            }),
        }
    }
    runs
}

use datafusion::arrow::array::{BooleanArray, Int64Array};
use datafusion::arrow::compute::{cast, filter_record_batch};
use datafusion::arrow::datatypes::DataType;
use datafusion::arrow::record_batch::RecordBatch;

/// Rank gate for one run: drops the position-tie rows that belong to the
/// previous buffer, routes warm-up rows and emit rows, and stops at the emit
/// end. Same semantics as the sharded worker's inline gate.
#[derive(Debug)]
pub(crate) struct RunGate {
    to_skip: usize,
    next_rank: usize,
    emit_start_row: usize,
    emit_end_row: usize,
    warm_up_start_row: usize,
    reached_emit: bool,
}

pub(crate) struct GateOutput {
    pub warm_up: Option<RecordBatch>,
    pub emit: Option<RecordBatch>,
    /// True on the batch that first reaches `emit_start_row` (or when there
    /// is no warm-up region, on the first non-empty batch). The caller replays
    /// the collected warm-up rows before pushing `emit`.
    pub reached_emit: bool,
    /// True once every row up to `emit_end_row` has been seen.
    pub done: bool,
}

impl RunGate {
    pub fn new(
        skip_leading_rows: usize,
        warm_up_start_row: usize,
        emit_start_row: usize,
        emit_end_row: usize,
    ) -> Self {
        Self {
            to_skip: skip_leading_rows,
            next_rank: warm_up_start_row,
            emit_start_row,
            emit_end_row,
            warm_up_start_row,
            reached_emit: false,
        }
    }

    pub fn needs_warm_up(&self) -> bool {
        self.emit_start_row > self.warm_up_start_row
    }

    pub fn feed(&mut self, mut batch: RecordBatch) -> GateOutput {
        let mut out = GateOutput {
            warm_up: None,
            emit: None,
            reached_emit: false,
            done: false,
        };
        if self.to_skip > 0 {
            let drop = self.to_skip.min(batch.num_rows());
            batch = batch.slice(drop, batch.num_rows() - drop);
            self.to_skip -= drop;
        }
        let n = batch.num_rows();
        if n == 0 {
            return out;
        }
        let batch_start = self.next_rank;
        let batch_end = batch_start + n;
        self.next_rank = batch_end;
        let warm_end = self.emit_start_row.saturating_sub(batch_start).min(n);
        let emit_to = self.emit_end_row.saturating_sub(batch_start).min(n);
        if warm_end > 0 {
            out.warm_up = Some(batch.slice(0, warm_end));
        }
        if !self.reached_emit && batch_end >= self.emit_start_row {
            self.reached_emit = true;
            out.reached_emit = true;
        }
        if emit_to > warm_end {
            out.emit = Some(batch.slice(warm_end, emit_to - warm_end));
        }
        out.done = batch_end >= self.emit_end_row;
        out
    }
}

/// Keep the rows whose `start` lies inside at least one interval. Used after
/// annotation (before projection) because stateful runs annotate whole
/// buffers and indexed reads are overlap-based.
pub(crate) fn filter_batch_to_bounds(
    batch: &RecordBatch,
    start_idx: usize,
    bounds: &[RunBounds],
) -> Result<RecordBatch> {
    if bounds.iter().any(RunBounds::is_open) || batch.num_rows() == 0 {
        return Ok(batch.clone());
    }
    let starts = cast(batch.column(start_idx), &DataType::Int64)?;
    let starts = starts
        .as_any()
        .downcast_ref::<Int64Array>()
        .ok_or_else(|| DataFusionError::Internal("start column did not cast to Int64".into()))?;
    let mask: BooleanArray = starts
        .iter()
        .map(|v| v.map(|pos| bounds.iter().any(|b| b.contains(pos))))
        .collect();
    Ok(filter_record_batch(batch, &mask)?)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn parse_absent_and_present() {
        assert_eq!(parse_regions_option(None).unwrap(), None);
        assert_eq!(
            parse_regions_option(Some(r#"{"everything":true}"#)).unwrap(),
            None
        );
        let specs = parse_regions_option(Some(
            r#"{"regions":[{"chrom":"chr1","start":10,"end":20},{"chrom":"chr2"},{"chrom":"chrX","start":5},{"chrom":"chrY","end":9}]}"#,
        ))
        .unwrap()
        .unwrap();
        assert_eq!(specs.len(), 4);
        assert_eq!(
            specs[0],
            RegionSpec {
                chrom: "chr1".into(),
                start: Some(10),
                end: Some(20)
            }
        );
        assert_eq!(
            specs[1],
            RegionSpec {
                chrom: "chr2".into(),
                start: None,
                end: None
            }
        );
        assert_eq!(
            specs[2],
            RegionSpec {
                chrom: "chrX".into(),
                start: Some(5),
                end: None
            }
        );
        assert_eq!(
            specs[3],
            RegionSpec {
                chrom: "chrY".into(),
                start: None,
                end: Some(9)
            }
        );
    }

    #[test]
    fn parse_rejects_bad_shapes() {
        for bad in [
            r#"{"regions":"chr1"}"#,
            r#"{"regions":[]}"#,
            r#"{"regions":[{"start":1}]}"#,
            r#"{"regions":[{"chrom":""}]}"#,
            r#"{"regions":[{"chrom":"chr1","start":0}]}"#,
            r#"{"regions":[{"chrom":"chr1","start":"5"}]}"#,
            r#"{"regions":[{"chrom":"chr1","start":9,"end":8}]}"#,
            r#"{"regions":[{"chrom":"chr`1"}]}"#,
        ] {
            let err = parse_regions_option(Some(bad)).unwrap_err().to_string();
            assert!(err.contains("regions"), "{bad}: {err}");
        }
    }

    #[test]
    fn bounds_contains_and_open() {
        assert!(RunBounds::OPEN.is_open());
        assert!(RunBounds::OPEN.contains(1));
        let b = RunBounds {
            lo: Some(10),
            hi: Some(20),
        };
        assert!(b.contains(10) && b.contains(20));
        assert!(!b.contains(9) && !b.contains(21));
        assert!(
            RunBounds {
                lo: None,
                hi: Some(5)
            }
            .contains(1)
        );
        assert!(
            !RunBounds {
                lo: Some(5),
                hi: None
            }
            .contains(4)
        );
    }

    #[test]
    fn merge_bounds_sorts_and_merges_overlapping_and_adjacent() {
        let merged = merge_bounds(vec![
            RunBounds {
                lo: Some(30),
                hi: Some(40),
            },
            RunBounds {
                lo: Some(10),
                hi: Some(20),
            },
            RunBounds {
                lo: Some(21),
                hi: Some(25),
            }, // adjacent to 10-20
            RunBounds {
                lo: Some(35),
                hi: Some(50),
            }, // overlaps 30-40
        ]);
        assert_eq!(
            merged,
            vec![
                RunBounds {
                    lo: Some(10),
                    hi: Some(25)
                },
                RunBounds {
                    lo: Some(30),
                    hi: Some(50)
                }
            ]
        );
    }

    #[test]
    fn merge_bounds_open_sides_absorb() {
        let merged = merge_bounds(vec![
            RunBounds {
                lo: Some(100),
                hi: Some(200),
            },
            RunBounds {
                lo: None,
                hi: Some(50),
            },
            RunBounds {
                lo: Some(150),
                hi: None,
            },
        ]);
        assert_eq!(
            merged,
            vec![
                RunBounds {
                    lo: None,
                    hi: Some(50)
                },
                RunBounds {
                    lo: Some(100),
                    hi: None
                }
            ]
        );
        assert_eq!(
            merge_bounds(vec![
                RunBounds::OPEN,
                RunBounds {
                    lo: Some(1),
                    hi: Some(2)
                }
            ]),
            vec![RunBounds::OPEN]
        );
    }

    fn s(chrom: &str, start: Option<i64>, end: Option<i64>) -> RegionSpec {
        RegionSpec {
            chrom: chrom.into(),
            start,
            end,
        }
    }

    #[test]
    fn resolve_matches_aliases_and_keeps_vcf_spelling() {
        let vcf = vec!["chr1".to_string(), "chr2".to_string(), "chrM".to_string()];
        let runs = resolve_regions(
            &[
                s("1", Some(10), Some(20)),
                s("chr1", Some(15), Some(30)),
                s("MT", None, None),
                s("chr9", None, None),
            ],
            &vcf,
        );
        assert_eq!(runs.len(), 2, "chr9 is not in the VCF: {runs:?}");
        assert_eq!(
            runs["chr1"],
            vec![RunBounds {
                lo: Some(10),
                hi: Some(30)
            }]
        );
        assert_eq!(runs["chrM"], vec![RunBounds::OPEN]);
    }

    #[test]
    fn restrict_contigs_keeps_vcf_order_and_drops_unlisted() {
        let mut runs = ContigRuns::new();
        runs.insert("chr3".into(), vec![RunBounds::OPEN]);
        runs.insert("chr1".into(), vec![RunBounds::OPEN]);
        let kept = restrict_contigs(vec!["chr1".into(), "chr2".into(), "chr3".into()], &runs);
        assert_eq!(kept, vec!["chr1".to_string(), "chr3".to_string()]);
        assert!(restrict_contigs(vec!["chr2".into()], &runs).is_empty());
    }

    #[test]
    fn buffer_range_maps_to_whole_buffers() {
        // Buffers: [0..100), [100..200), [200..300), [300..]
        let pos = [10, 100, 200, 300, i64::MAX];
        assert_eq!(
            buffer_range_for_bounds(
                &pos,
                RunBounds {
                    lo: Some(150),
                    hi: Some(250)
                }
            ),
            (1, 3)
        );
        // exactly on a boundary
        assert_eq!(
            buffer_range_for_bounds(
                &pos,
                RunBounds {
                    lo: Some(100),
                    hi: Some(200)
                }
            ),
            (1, 3)
        );
        // below the first row: clamp to buffer 0
        assert_eq!(
            buffer_range_for_bounds(
                &pos,
                RunBounds {
                    lo: Some(1),
                    hi: Some(5)
                }
            ),
            (0, 1)
        );
        // beyond the last boundary: through the last buffer
        assert_eq!(
            buffer_range_for_bounds(
                &pos,
                RunBounds {
                    lo: Some(350),
                    hi: Some(900)
                }
            ),
            (3, 4)
        );
        // open sides
        assert_eq!(
            buffer_range_for_bounds(
                &pos,
                RunBounds {
                    lo: None,
                    hi: Some(150)
                }
            ),
            (0, 2)
        );
        assert_eq!(
            buffer_range_for_bounds(
                &pos,
                RunBounds {
                    lo: Some(250),
                    hi: None
                }
            ),
            (2, 4)
        );
        assert_eq!(buffer_range_for_bounds(&pos, RunBounds::OPEN), (0, 4));
    }

    #[test]
    fn buffer_range_on_position_tie_includes_the_earlier_buffer() {
        // A boundary at pos 100 with rows_before_pos > 0 means rows at 100 sit
        // in BOTH buffer 0 and buffer 1; a run starting at 100 maps to buffer 1
        // and the gate's skip_leading_rows drops the earlier buffer's ties.
        let pos = [10, 100, 200, i64::MAX];
        assert_eq!(
            buffer_range_for_bounds(
                &pos,
                RunBounds {
                    lo: Some(100),
                    hi: Some(100)
                }
            ),
            (1, 2)
        );
        // lo strictly inside buffer 0 that ends with a tie at 100: buffer 0.
        assert_eq!(
            buffer_range_for_bounds(
                &pos,
                RunBounds {
                    lo: Some(99),
                    hi: Some(100)
                }
            ),
            (0, 2)
        );
    }

    #[test]
    fn plan_runs_without_grid_is_one_run_per_interval() {
        let b = [
            RunBounds {
                lo: Some(1),
                hi: Some(5),
            },
            RunBounds {
                lo: Some(50),
                hi: None,
            },
        ];
        assert_eq!(
            plan_runs(&b, None),
            vec![
                RunPlan {
                    bounds: vec![b[0]],
                    buffers: None
                },
                RunPlan {
                    bounds: vec![b[1]],
                    buffers: None
                },
            ]
        );
    }

    #[test]
    fn plan_runs_with_grid_merges_touching_buffer_ranges() {
        let pos = [10, 100, 200, 300, 400, i64::MAX];
        let b = [
            RunBounds {
                lo: Some(20),
                hi: Some(30),
            }, // buffer 0
            RunBounds {
                lo: Some(120),
                hi: Some(130),
            }, // buffer 1 -> touches (0,1)+(1,2)
            RunBounds {
                lo: Some(350),
                hi: Some(360),
            }, // buffer 3
        ];
        assert_eq!(
            plan_runs(&b, Some(&pos)),
            vec![
                RunPlan {
                    bounds: vec![b[0], b[1]],
                    buffers: Some((0, 2))
                },
                RunPlan {
                    bounds: vec![b[2]],
                    buffers: Some((3, 4))
                },
            ]
        );
    }

    #[test]
    fn plan_runs_with_empty_grid_yields_no_runs() {
        assert!(plan_runs(&[RunBounds::OPEN], Some(&[i64::MAX])).is_empty());
        assert!(plan_runs(&[RunBounds::OPEN], Some(&[])).is_empty());
    }

    use datafusion::arrow::array::{Array, Int64Array, StringArray};
    use datafusion::arrow::datatypes::{DataType, Field, Schema};
    use std::sync::Arc;

    fn batch(starts: &[i64]) -> RecordBatch {
        let schema = Arc::new(Schema::new(vec![
            Field::new("chrom", DataType::Utf8, false),
            Field::new("start", DataType::Int64, false),
        ]));
        RecordBatch::try_new(
            schema,
            vec![
                Arc::new(StringArray::from(vec!["chr1"; starts.len()])) as Arc<dyn Array>,
                Arc::new(Int64Array::from(starts.to_vec())) as Arc<dyn Array>,
            ],
        )
        .unwrap()
    }

    fn starts(b: &RecordBatch) -> Vec<i64> {
        b.column(1)
            .as_any()
            .downcast_ref::<Int64Array>()
            .unwrap()
            .values()
            .to_vec()
    }

    #[test]
    fn gate_splits_skip_warmup_emit_and_stops() {
        // ranks after skipping 1 tie row: 10,11,12 | 13,14,15 | 16,17
        let mut gate = RunGate::new(1, 10, 13, 16);
        assert!(gate.needs_warm_up());
        let out = gate.feed(batch(&[99, 100, 101])); // 99 skipped; ranks 10,11 warm-up
        assert_eq!(starts(out.warm_up.as_ref().unwrap()), vec![100, 101]);
        assert!(out.emit.is_none() && !out.reached_emit && !out.done);
        let out = gate.feed(batch(&[102, 103, 104])); // rank 12 warm-up, 13,14 emit
        assert_eq!(starts(out.warm_up.as_ref().unwrap()), vec![102]);
        assert_eq!(starts(out.emit.as_ref().unwrap()), vec![103, 104]);
        assert!(out.reached_emit && !out.done);
        let out = gate.feed(batch(&[105, 106, 107])); // 15 emit; 16,17 discarded
        assert!(out.warm_up.is_none());
        assert_eq!(starts(out.emit.as_ref().unwrap()), vec![105]);
        assert!(out.done);
    }

    #[test]
    fn gate_without_warm_up_emits_from_the_first_row() {
        let mut gate = RunGate::new(0, 0, 0, 2);
        assert!(!gate.needs_warm_up());
        let out = gate.feed(batch(&[1, 2, 3]));
        assert!(out.warm_up.is_none());
        assert_eq!(starts(out.emit.as_ref().unwrap()), vec![1, 2]);
        assert!(out.reached_emit && out.done);
    }

    #[test]
    fn gate_reports_reached_emit_once_and_handles_empty_batches() {
        let mut gate = RunGate::new(0, 0, 2, usize::MAX);
        let out = gate.feed(batch(&[]));
        assert!(out.warm_up.is_none() && out.emit.is_none() && !out.reached_emit);
        let out = gate.feed(batch(&[1, 2, 3]));
        assert!(out.reached_emit);
        let out = gate.feed(batch(&[4]));
        assert!(
            !out.reached_emit,
            "reached_emit fires only on the crossing batch"
        );
        assert_eq!(starts(out.emit.as_ref().unwrap()), vec![4]);
        assert!(!out.done, "open upper rank never stops");
    }

    #[test]
    fn filter_batch_keeps_rows_inside_any_interval() {
        let b = batch(&[5, 10, 15, 20, 25, 30]);
        let out = filter_batch_to_bounds(
            &b,
            1,
            &[
                RunBounds {
                    lo: Some(10),
                    hi: Some(15),
                },
                RunBounds {
                    lo: Some(30),
                    hi: None,
                },
            ],
        )
        .unwrap();
        assert_eq!(starts(&out), vec![10, 15, 30]);
        let same = filter_batch_to_bounds(&b, 1, &[RunBounds::OPEN]).unwrap();
        assert_eq!(same.num_rows(), 6);
    }

    #[test]
    fn filter_batch_accepts_uint32_start() {
        use datafusion::arrow::array::UInt32Array;
        let schema = Arc::new(Schema::new(vec![Field::new(
            "start",
            DataType::UInt32,
            false,
        )]));
        let b = RecordBatch::try_new(
            schema,
            vec![Arc::new(UInt32Array::from(vec![1u32, 7, 9])) as Arc<dyn Array>],
        )
        .unwrap();
        let out = filter_batch_to_bounds(
            &b,
            0,
            &[RunBounds {
                lo: Some(7),
                hi: Some(9),
            }],
        )
        .unwrap();
        assert_eq!(out.num_rows(), 2);
    }
}
