use std::collections::HashMap;

use crate::{QcModule, TidyRow};

/// FastQC "Sequence Length Distribution": histogram of read lengths.
#[derive(Debug, Default)]
pub struct SeqLength {
    lengths: HashMap<usize, u64>,
}

impl SeqLength {
    pub fn new() -> Self {
        Self::default()
    }
}

impl QcModule for SeqLength {
    fn as_any(&self) -> &dyn std::any::Any {
        self
    }

    fn name(&self) -> &'static str {
        "seq_length"
    }

    fn update(&mut self, _name: &[u8], seq: &[u8], _qual: &[u8]) {
        *self.lengths.entry(seq.len()).or_insert(0) += 1;
    }

    fn merge(&mut self, other: &dyn QcModule) {
        let o = other
            .as_any()
            .downcast_ref::<SeqLength>()
            .expect("merge type mismatch");
        for (&len, &c) in &o.lengths {
            *self.lengths.entry(len).or_insert(0) += c;
        }
    }

    fn finalize(&self, out: &mut Vec<TidyRow>) {
        let m = "seq_length";
        let mut lens: Vec<usize> = self.lengths.keys().copied().collect();
        lens.sort_unstable();
        for len in &lens {
            out.push(TidyRow {
                module: m,
                label: None,
                position: Some(*len as i32),
                metric: "count".to_string(),
                value: Some(self.lengths[len] as f64),
                value_str: None,
            });
        }
        // FastQC: warn if reads are not all the same length; error if any read
        // has length 0.
        let status = if lens.contains(&0) {
            "FAIL"
        } else if lens.len() > 1 {
            "WARN"
        } else {
            "PASS"
        };
        out.push(TidyRow::status(m, status));
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn seq_length_histogram_and_status() {
        let mut m = SeqLength::new();
        m.update(b"", b"ACGT", b"IIII");
        m.update(b"", b"ACGT", b"IIII");
        m.update(b"", b"ACGTACGT", b"IIIIIIII");
        let mut rows = Vec::new();
        m.finalize(&mut rows);
        let at = |l: i32| {
            rows.iter()
                .find(|r| r.position == Some(l) && r.metric == "count")
                .and_then(|r| r.value)
                .unwrap_or(0.0)
        };
        assert_eq!(at(4), 2.0);
        assert_eq!(at(8), 1.0);
        // multiple lengths -> WARN
        assert!(
            rows.iter()
                .any(|r| r.metric == "status" && r.value_str.as_deref() == Some("WARN"))
        );
    }
}
