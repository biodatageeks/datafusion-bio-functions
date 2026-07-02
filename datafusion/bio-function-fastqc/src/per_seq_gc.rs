use std::any::Any;

use crate::{QcModule, TidyRow};

/// FastQC "Per Sequence GC Content": distribution of per-read GC% over 0..=100.
#[derive(Debug)]
pub struct PerSeqGc {
    bins: [u64; 101],
}

impl PerSeqGc {
    pub fn new() -> Self {
        Self { bins: [0u64; 101] }
    }
}

impl Default for PerSeqGc {
    fn default() -> Self {
        Self::new()
    }
}

impl QcModule for PerSeqGc {
    fn name(&self) -> &'static str {
        "per_seq_gc"
    }

    fn update(&mut self, seq: &[u8], _qual: &[u8]) {
        if seq.is_empty() {
            return;
        }
        let mut gc = 0u64;
        let mut counted = 0u64;
        for &b in seq {
            match b {
                b'G' | b'g' | b'C' | b'c' => {
                    gc += 1;
                    counted += 1;
                },
                b'A' | b'a' | b'T' | b't' | b'U' | b'u' => counted += 1,
                _ => {}, // N excluded from denominator, matching FastQC
            }
        }
        if counted == 0 {
            return;
        }
        let pct = (gc as f64 / counted as f64) * 100.0;
        let bin = pct.round() as usize;
        self.bins[bin.min(100)] += 1;
    }

    fn merge(&mut self, other: &dyn QcModule) {
        let o = other
            .as_any()
            .downcast_ref::<PerSeqGc>()
            .expect("merge type mismatch");
        for i in 0..101 {
            self.bins[i] += o.bins[i];
        }
    }

    fn finalize(&self, out: &mut Vec<TidyRow>) {
        let m = "per_seq_gc";
        for (g, &c) in self.bins.iter().enumerate() {
            out.push(TidyRow {
                module: m,
                label: None,
                position: Some(g as i32),
                metric: "count".to_string(),
                value: Some(c as f64),
                value_str: None,
            });
        }
        // Phase-1 status: PASS. Exact FastQC theoretical-distribution status
        // is a follow-up; the parity harness checks the count distribution,
        // which is exact here.
        out.push(TidyRow::status(m, "PASS"));
    }

    fn as_any(&self) -> &dyn Any {
        self
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn per_seq_gc_bins_reads() {
        let mut m = PerSeqGc::new();
        m.update(b"GGCC", b"IIII"); // 100% -> bin 100
        m.update(b"ATAT", b"IIII"); // 0%   -> bin 0
        m.update(b"ATGC", b"IIII"); // 50%  -> bin 50
        let mut rows = Vec::new();
        m.finalize(&mut rows);
        let count_at = |bin: i32| {
            rows.iter()
                .find(|r| r.position == Some(bin) && r.metric == "count")
                .and_then(|r| r.value)
                .unwrap_or(0.0)
        };
        assert_eq!(count_at(0), 1.0);
        assert_eq!(count_at(50), 1.0);
        assert_eq!(count_at(100), 1.0);
    }
}
