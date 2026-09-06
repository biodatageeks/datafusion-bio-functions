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
}
