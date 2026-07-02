use std::any::Any;

use crate::{QcModule, TidyRow};

/// FastQC "Per Base Sequence Content": per-position percentage of G/A/T/C
/// (N and other bases are ignored). Emitted with metrics G, A, T, C to match
/// FastQC's column order.
#[derive(Debug, Default)]
pub struct PerBaseContent {
    /// position -> [g, a, t, c]
    counts: Vec<[u64; 4]>,
}

impl PerBaseContent {
    pub fn new() -> Self {
        Self::default()
    }

    fn ensure_len(&mut self, len: usize) {
        if self.counts.len() < len {
            self.counts.resize(len, [0u64; 4]);
        }
    }
}

impl QcModule for PerBaseContent {
    fn name(&self) -> &'static str {
        "per_base_content"
    }

    fn update(&mut self, seq: &[u8], _qual: &[u8]) {
        self.ensure_len(seq.len());
        for (i, &b) in seq.iter().enumerate() {
            match b {
                b'G' | b'g' => self.counts[i][0] += 1,
                b'A' | b'a' => self.counts[i][1] += 1,
                b'T' | b't' => self.counts[i][2] += 1,
                b'C' | b'c' => self.counts[i][3] += 1,
                _ => {},
            }
        }
    }

    fn merge(&mut self, other: &dyn QcModule) {
        let o = other
            .as_any()
            .downcast_ref::<PerBaseContent>()
            .expect("merge type mismatch");
        self.ensure_len(o.counts.len());
        for (s, o) in self.counts.iter_mut().zip(o.counts.iter()) {
            for k in 0..4 {
                s[k] += o[k];
            }
        }
    }

    fn finalize(&self, out: &mut Vec<TidyRow>) {
        let m = "per_base_content";
        let mut worst_diff = 0f64;
        for (i, c) in self.counts.iter().enumerate() {
            let total = c[0] + c[1] + c[2] + c[3];
            if total == 0 {
                continue;
            }
            let pos = (i + 1) as i32;
            let pct: Vec<f64> = c.iter().map(|&x| x as f64 / total as f64 * 100.0).collect();
            // metrics in FastQC column order: G, A, T, C
            for (metric, v) in [("G", pct[0]), ("A", pct[1]), ("T", pct[2]), ("C", pct[3])] {
                out.push(TidyRow {
                    module: m,
                    label: None,
                    position: Some(pos),
                    metric: metric.to_string(),
                    value: Some(v),
                    value_str: None,
                });
            }
            // FastQC flags |%A-%T| or |%G-%C| deviations.
            worst_diff = worst_diff.max((pct[1] - pct[2]).abs()).max((pct[0] - pct[3]).abs());
        }
        let status = if worst_diff > 20.0 {
            "FAIL"
        } else if worst_diff > 10.0 {
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
    fn per_base_content_percentages() {
        let mut m = PerBaseContent::new();
        m.update(b"GA", b"II");
        m.update(b"GT", b"II");
        let mut rows = Vec::new();
        m.finalize(&mut rows);
        let at = |p: i32, metric: &str| {
            rows.iter()
                .find(|r| r.position == Some(p) && r.metric == metric)
                .and_then(|r| r.value)
                .unwrap()
        };
        // pos1: both G -> G=100%
        assert!((at(1, "G") - 100.0).abs() < 1e-9);
        // pos2: A and T -> 50% each
        assert!((at(2, "A") - 50.0).abs() < 1e-9);
        assert!((at(2, "T") - 50.0).abs() < 1e-9);
    }
}
