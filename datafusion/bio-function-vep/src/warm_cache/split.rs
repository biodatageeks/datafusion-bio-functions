use std::collections::{BTreeMap, BTreeSet};
use std::fmt;

#[derive(Debug, Clone, Copy)]
pub struct FrequencyFields<'a> {
    pub minor_allele_freq: Option<f64>,
    pub af: Option<&'a str>,
    pub gnomadg: Option<&'a str>,
    pub gnomade: Option<&'a str>,
}

#[derive(Debug, Clone, Copy)]
pub struct PositionFrequency {
    pub start: i64,
    pub max_af: f64,
}

#[derive(Debug, Clone, PartialEq)]
pub enum WarmSplitError {
    UnsortedPosition { previous: i64, current: i64 },
}

impl fmt::Display for WarmSplitError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::UnsortedPosition { previous, current } => write!(
                f,
                "variation rows must be sorted by start for streaming warm selection: previous={previous}, current={current}"
            ),
        }
    }
}

impl std::error::Error for WarmSplitError {}

pub type Result<T> = std::result::Result<T, WarmSplitError>;

#[derive(Debug)]
pub struct WarmPositionCollector {
    af_threshold: f64,
    position_radius: i64,
    current_start: Option<i64>,
    current_max_af: f64,
    warm_positions: BTreeSet<i64>,
}

impl WarmPositionCollector {
    pub fn new(af_threshold: f64, position_radius: i64) -> Self {
        Self {
            af_threshold,
            position_radius: position_radius.max(0),
            current_start: None,
            current_max_af: 0.0,
            warm_positions: BTreeSet::new(),
        }
    }

    pub fn push(&mut self, start: i64, max_af: f64) -> Result<()> {
        let max_af = if max_af.is_finite() { max_af } else { 0.0 };

        match self.current_start {
            None => {
                self.current_start = Some(start);
                self.current_max_af = max_af;
            }
            Some(current) if start == current => {
                self.current_max_af = self.current_max_af.max(max_af);
            }
            Some(current) if start > current => {
                self.flush_current();
                self.current_start = Some(start);
                self.current_max_af = max_af;
            }
            Some(previous) => {
                return Err(WarmSplitError::UnsortedPosition {
                    previous,
                    current: start,
                });
            }
        }

        Ok(())
    }

    pub fn finish(mut self) -> BTreeSet<i64> {
        self.flush_current();
        self.warm_positions
    }

    fn flush_current(&mut self) {
        let Some(start) = self.current_start.take() else {
            return;
        };

        if self.current_max_af < self.af_threshold {
            self.current_max_af = 0.0;
            return;
        }

        for offset in -self.position_radius..=self.position_radius {
            self.warm_positions.insert(start.saturating_add(offset));
        }
        self.current_max_af = 0.0;
    }
}

pub fn max_af_from_pairs(value: Option<&str>) -> f64 {
    let Some(value) = value else {
        return 0.0;
    };

    value
        .split(',')
        .filter_map(|part| {
            let raw = part
                .rsplit_once(':')
                .map(|(_, freq)| freq)
                .unwrap_or(part)
                .trim();
            raw.parse::<f64>().ok().filter(|v| v.is_finite())
        })
        .fold(0.0, f64::max)
}

pub fn max_global_af(fields: &FrequencyFields<'_>) -> f64 {
    [
        fields
            .minor_allele_freq
            .filter(|value| value.is_finite())
            .unwrap_or(0.0),
        max_af_from_pairs(fields.af),
        max_af_from_pairs(fields.gnomadg),
        max_af_from_pairs(fields.gnomade),
    ]
    .into_iter()
    .fold(0.0, f64::max)
}

pub fn select_warm_positions<I>(rows: I, af_threshold: f64, position_radius: i64) -> BTreeSet<i64>
where
    I: IntoIterator<Item = PositionFrequency>,
{
    let mut max_af_by_position: BTreeMap<i64, f64> = BTreeMap::new();
    for row in rows {
        let max_af = if row.max_af.is_finite() {
            row.max_af
        } else {
            0.0
        };
        max_af_by_position
            .entry(row.start)
            .and_modify(|existing| *existing = existing.max(max_af))
            .or_insert(max_af);
    }

    let source_positions: BTreeSet<i64> = max_af_by_position.keys().copied().collect();
    let radius = position_radius.max(0);
    let mut warm_positions = BTreeSet::new();

    for (start, max_af) in max_af_by_position {
        if max_af < af_threshold {
            continue;
        }

        for offset in -radius..=radius {
            let candidate = start.saturating_add(offset);
            if source_positions.contains(&candidate) {
                warm_positions.insert(candidate);
            }
        }
    }

    warm_positions
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn max_af_parses_pair_lists() {
        assert_eq!(max_af_from_pairs(Some("A:0.1,C:2.6e-05")), 0.1);
        assert_eq!(max_af_from_pairs(Some("A:0,C:0")), 0.0);
        assert_eq!(max_af_from_pairs(None), 0.0);
    }

    #[test]
    fn max_af_parses_scalar_strings() {
        assert_eq!(max_af_from_pairs(Some("0.0125")), 0.0125);
        assert_eq!(max_af_from_pairs(Some("2.6e-05")), 2.6e-05);
    }

    #[test]
    fn max_global_af_uses_largest_frequency_column() {
        let af = max_global_af(&FrequencyFields {
            minor_allele_freq: Some(0.004),
            af: Some("A:0.001,C:0.003"),
            gnomadg: Some("G:0.018"),
            gnomade: Some("T:0.007"),
        });

        assert_eq!(af, 0.018);
    }

    #[test]
    fn warm_position_selection_keeps_complete_positions() {
        let rows = [
            PositionFrequency {
                start: 99,
                max_af: 0.0,
            },
            PositionFrequency {
                start: 100,
                max_af: 0.0,
            },
            PositionFrequency {
                start: 100,
                max_af: 0.02,
            },
            PositionFrequency {
                start: 101,
                max_af: 0.0,
            },
            PositionFrequency {
                start: 103,
                max_af: 0.0,
            },
        ];

        let selected = select_warm_positions(rows, 0.01, 1);

        assert!(selected.contains(&99));
        assert!(selected.contains(&100));
        assert!(selected.contains(&101));
        assert!(!selected.contains(&103));
    }

    #[test]
    fn warm_position_selection_uses_position_max_af() {
        let rows = [
            PositionFrequency {
                start: 10,
                max_af: 0.0,
            },
            PositionFrequency {
                start: 10,
                max_af: 0.015,
            },
            PositionFrequency {
                start: 11,
                max_af: 0.0,
            },
        ];

        let selected = select_warm_positions(rows, 0.01, 0);

        assert_eq!(selected.iter().copied().collect::<Vec<_>>(), vec![10]);
    }

    #[test]
    fn streaming_collector_collapses_duplicate_sorted_positions() {
        let mut collector = WarmPositionCollector::new(0.01, 1);
        collector.push(100, 0.0).unwrap();
        collector.push(100, 0.02).unwrap();
        collector.push(101, 0.0).unwrap();

        let selected = collector.finish();

        assert_eq!(
            selected.iter().copied().collect::<Vec<_>>(),
            vec![99, 100, 101]
        );
    }

    #[test]
    fn streaming_collector_rejects_unsorted_positions() {
        let mut collector = WarmPositionCollector::new(0.01, 1);
        collector.push(100, 0.02).unwrap();

        let err = collector.push(99, 0.0).unwrap_err();

        assert_eq!(
            err,
            WarmSplitError::UnsortedPosition {
                previous: 100,
                current: 99
            }
        );
    }
}
