use std::any::Any;

use crate::{QcModule, TidyRow};

/// Max phred value tracked (0..=93 covers all realistic phred scores).
const QUAL_MAX: usize = 94;

/// FastQC "Per Base Sequence Quality": per-position phred distribution,
/// summarized as mean + median + quartiles + 10th/90th percentiles.
#[derive(Debug, Default)]
pub struct PerBaseQuality {
    /// position -> histogram of phred values (0..QUAL_MAX)
    hist: Vec<[u64; QUAL_MAX]>,
}

impl PerBaseQuality {
    pub fn new() -> Self {
        Self::default()
    }

    fn ensure_len(&mut self, len: usize) {
        if self.hist.len() < len {
            self.hist.resize(len, [0u64; QUAL_MAX]);
        }
    }
}

/// FastQC "nth value" percentile over an integer histogram.
fn percentile(hist: &[u64; QUAL_MAX], total: u64, p: f64) -> f64 {
    if total == 0 {
        return 0.0;
    }
    let rank = (p * total as f64).ceil().max(1.0) as u64;
    let mut cum = 0u64;
    for (q, &c) in hist.iter().enumerate() {
        cum += c;
        if cum >= rank {
            return q as f64;
        }
    }
    (QUAL_MAX - 1) as f64
}

impl QcModule for PerBaseQuality {
    fn name(&self) -> &'static str {
        "per_base_quality"
    }

    fn update(&mut self, _name: &[u8], _seq: &[u8], qual: &[u8]) {
        self.ensure_len(qual.len());
        for (i, &q) in qual.iter().enumerate() {
            let phred = (q.saturating_sub(33)) as usize;
            let phred = phred.min(QUAL_MAX - 1);
            self.hist[i][phred] += 1;
        }
    }

    fn merge(&mut self, other: &dyn QcModule) {
        let o = other
            .as_any()
            .downcast_ref::<PerBaseQuality>()
            .expect("merge type mismatch");
        self.ensure_len(o.hist.len());
        for (sh, oh) in self.hist.iter_mut().zip(o.hist.iter()) {
            for (s, &o) in sh.iter_mut().zip(oh.iter()) {
                *s += o;
            }
        }
    }

    fn finalize(&self, out: &mut Vec<TidyRow>) {
        let m = "per_base_quality";
        let mut worst_median = f64::INFINITY;
        let mut worst_q1 = f64::INFINITY;
        for (i, h) in self.hist.iter().enumerate() {
            let total: u64 = h.iter().sum();
            if total == 0 {
                continue;
            }
            let pos = (i + 1) as i32;
            let sum: u64 = h.iter().enumerate().map(|(q, &c)| q as u64 * c).sum();
            let mean = sum as f64 / total as f64;
            let median = percentile(h, total, 0.50);
            let q1 = percentile(h, total, 0.25);
            let q3 = percentile(h, total, 0.75);
            let p10 = percentile(h, total, 0.10);
            let p90 = percentile(h, total, 0.90);
            worst_median = worst_median.min(median);
            worst_q1 = worst_q1.min(q1);
            for (metric, v) in [
                ("mean", mean),
                ("median", median),
                ("q1", q1),
                ("q3", q3),
                ("p10", p10),
                ("p90", p90),
            ] {
                out.push(TidyRow {
                    module: m,
                    label: None,
                    position: Some(pos),
                    metric: metric.to_string(),
                    value: Some(v),
                    value_str: None,
                });
            }
        }
        let status = if worst_q1 < 5.0 || worst_median < 20.0 {
            "FAIL"
        } else if worst_q1 < 10.0 || worst_median < 25.0 {
            "WARN"
        } else {
            "PASS"
        };
        out.push(TidyRow::status(m, status));
    }

    fn as_any(&self) -> &dyn Any {
        self
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn per_base_quality_mean_per_position() {
        let mut m = PerBaseQuality::new();
        // '!' = phred 0, 'I' = phred 40, '5' = phred 20
        m.update(b"", b"AA", b"!I"); // pos1 -> 0, pos2 -> 40
        m.update(b"", b"AA", b"I5"); // pos1 -> 40, pos2 -> 20
        let mut rows = Vec::new();
        m.finalize(&mut rows);
        let mean_at = |pos: i32| {
            rows.iter()
                .find(|r| r.position == Some(pos) && r.metric == "mean")
                .and_then(|r| r.value)
                .unwrap()
        };
        assert!((mean_at(1) - 20.0).abs() < 1e-9);
        assert!((mean_at(2) - 30.0).abs() < 1e-9);
    }

    #[test]
    fn per_base_quality_merge_variable_lengths() {
        // Merging accumulators built from different-length reads must extend
        // the shorter histogram (ensure_len), not panic or drop positions.
        let mut a = PerBaseQuality::new();
        a.update(b"", b"AAA", b"III"); // 3 positions, phred 40
        let mut b = PerBaseQuality::new();
        b.update(b"", b"AAAAA", b"!!!!!"); // 5 positions, phred 0
        a.merge(&b);
        let mut rows = Vec::new();
        a.finalize(&mut rows);
        let mean_at = |pos: i32| {
            rows.iter()
                .find(|r| r.position == Some(pos) && r.metric == "mean")
                .and_then(|r| r.value)
        };
        // pos1: (40 + 0) / 2 = 20; pos5 only from the longer read = 0
        assert!((mean_at(1).unwrap() - 20.0).abs() < 1e-9);
        assert!((mean_at(5).unwrap() - 0.0).abs() < 1e-9);
    }
}
