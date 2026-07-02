use std::any::Any;

use crate::{QcModule, TidyRow};

const QUAL_MAX: usize = 94;

/// FastQC "Per Sequence Quality Scores": distribution of each read's mean
/// quality. FastQC computes the mean of the raw ASCII quality bytes with
/// integer division, then subtracts the Phred+33 offset — equivalently
/// floor(sum(phred)/len).
#[derive(Debug)]
pub struct PerSeqQuality {
    hist: [u64; QUAL_MAX],
}

impl PerSeqQuality {
    pub fn new() -> Self {
        Self {
            hist: [0u64; QUAL_MAX],
        }
    }
}

impl Default for PerSeqQuality {
    fn default() -> Self {
        Self::new()
    }
}

impl QcModule for PerSeqQuality {
    fn name(&self) -> &'static str {
        "per_seq_quality"
    }

    fn update(&mut self, _name: &[u8], _seq: &[u8], qual: &[u8]) {
        if qual.is_empty() {
            return;
        }
        let sum: u64 = qual.iter().map(|&q| q as u64).sum();
        let mean_ascii = sum / qual.len() as u64; // integer division, like FastQC
        let phred = mean_ascii.saturating_sub(33) as usize;
        self.hist[phred.min(QUAL_MAX - 1)] += 1;
    }

    fn merge(&mut self, other: &dyn QcModule) {
        let o = other
            .as_any()
            .downcast_ref::<PerSeqQuality>()
            .expect("merge type mismatch");
        for (s, &v) in self.hist.iter_mut().zip(o.hist.iter()) {
            *s += v;
        }
    }

    fn finalize(&self, out: &mut Vec<TidyRow>) {
        let m = "per_seq_quality";
        let first = self.hist.iter().position(|&c| c > 0);
        let last = self.hist.iter().rposition(|&c| c > 0);
        let (Some(first), Some(last)) = (first, last) else {
            out.push(TidyRow::status(m, "PASS"));
            return;
        };
        // Emit the contiguous observed range (FastQC prints zeros within it).
        let mut mode_phred = first;
        for phred in first..=last {
            let c = self.hist[phred];
            if c > self.hist[mode_phred] {
                mode_phred = phred;
            }
            out.push(TidyRow {
                module: m,
                label: None,
                position: Some(phred as i32),
                metric: "count".to_string(),
                value: Some(c as f64),
                value_str: None,
            });
        }
        // FastQC: warn if the most frequent mean quality < 27, error if < 20.
        let status = if mode_phred < 20 {
            "FAIL"
        } else if mode_phred < 27 {
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
    fn per_seq_quality_bins_mean() {
        let mut m = PerSeqQuality::new();
        m.update(b"", b"AAAA", b"IIII"); // 'I'=73 -> mean 73 -> phred 40
        m.update(b"", b"AAAA", b"!!!!"); // '!'=33 -> mean 33 -> phred 0
        let mut rows = Vec::new();
        m.finalize(&mut rows);
        let count_at = |p: i32| {
            rows.iter()
                .find(|r| r.position == Some(p) && r.metric == "count")
                .and_then(|r| r.value)
                .unwrap_or(0.0)
        };
        assert_eq!(count_at(40), 1.0);
        assert_eq!(count_at(0), 1.0);
        // contiguous range 0..=40 emitted
        assert_eq!(count_at(20), 0.0);
    }
}
