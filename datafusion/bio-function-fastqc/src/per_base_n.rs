use std::any::Any;

use crate::{QcModule, TidyRow};

/// FastQC "Per Base N Content": per-position percentage of N bases (the output
/// header says "N-Count" but the values are percentages).
#[derive(Debug, Default)]
pub struct PerBaseN {
    /// position -> [n_count, non_n_count]
    counts: Vec<[u64; 2]>,
}

impl PerBaseN {
    pub fn new() -> Self {
        Self::default()
    }

    fn ensure_len(&mut self, len: usize) {
        if self.counts.len() < len {
            self.counts.resize(len, [0u64; 2]);
        }
    }
}

impl QcModule for PerBaseN {
    fn name(&self) -> &'static str {
        "per_base_n"
    }

    fn update(&mut self, seq: &[u8], _qual: &[u8]) {
        self.ensure_len(seq.len());
        for (i, &b) in seq.iter().enumerate() {
            if b == b'N' || b == b'n' {
                self.counts[i][0] += 1;
            } else {
                self.counts[i][1] += 1;
            }
        }
    }

    fn merge(&mut self, other: &dyn QcModule) {
        let o = other
            .as_any()
            .downcast_ref::<PerBaseN>()
            .expect("merge type mismatch");
        self.ensure_len(o.counts.len());
        for (s, o) in self.counts.iter_mut().zip(o.counts.iter()) {
            s[0] += o[0];
            s[1] += o[1];
        }
    }

    fn finalize(&self, out: &mut Vec<TidyRow>) {
        let m = "per_base_n";
        let mut worst = 0f64;
        for (i, c) in self.counts.iter().enumerate() {
            let total = c[0] + c[1];
            if total == 0 {
                continue;
            }
            let pct = c[0] as f64 / total as f64 * 100.0;
            worst = worst.max(pct);
            out.push(TidyRow {
                module: m,
                label: None,
                position: Some((i + 1) as i32),
                metric: "pct".to_string(),
                value: Some(pct),
                value_str: None,
            });
        }
        // FastQC: warn if any position exceeds 5% N, error above 20%.
        let status = if worst > 20.0 {
            "FAIL"
        } else if worst > 5.0 {
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
    fn per_base_n_percentage_per_position() {
        let mut m = PerBaseN::new();
        m.update(b"NAA", b"III"); // pos1 N
        m.update(b"AAA", b"III");
        let mut rows = Vec::new();
        m.finalize(&mut rows);
        let at = |p: i32| {
            rows.iter()
                .find(|r| r.position == Some(p) && r.metric == "pct")
                .and_then(|r| r.value)
                .unwrap()
        };
        assert!((at(1) - 50.0).abs() < 1e-9); // 1 of 2 reads has N at pos1
        assert!((at(2) - 0.0).abs() < 1e-9);
    }
}
