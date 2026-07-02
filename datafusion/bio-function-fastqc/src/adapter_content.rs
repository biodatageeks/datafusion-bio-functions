use std::any::Any;

use crate::{QcModule, TidyRow};

/// FastQC's built-in adapter k-mers (Configuration/adapter_list.txt). All are
/// 12 bp, so `LONGEST_ADAPTER` is 12.
const ADAPTERS: [(&str, &[u8]); 6] = [
    ("Illumina Universal Adapter", b"AGATCGGAAGAG"),
    ("Illumina Small RNA 3' Adapter", b"TGGAATTCTCGG"),
    ("Illumina Small RNA 5' Adapter", b"GATCGTCGGACT"),
    ("Nextera Transposase Sequence", b"CTGTCTCTTATA"),
    ("PolyA", b"AAAAAAAAAAAA"),
    ("PolyG", b"GGGGGGGGGGGG"),
];
const LONGEST_ADAPTER: usize = 12;

/// FastQC "Adapter Content": per-position cumulative fraction of reads in which
/// an adapter has appeared at or before that position.
///
/// FastQC increments each adapter from its first-found index to
/// `longestSequence - longestAdapter`. We instead accumulate a histogram of the
/// first-found index per adapter and take the cumulative sum at finalize — an
/// order-independent, parallel-merge-safe formulation that is identical to
/// FastQC for uniform read lengths (and uses the final longest length as the
/// upper bound for mixed lengths).
#[derive(Debug)]
pub struct AdapterContent {
    total: u64,
    longest_seq: usize,
    /// per adapter: histogram of the first index at which it was found
    first_index: Vec<Vec<u64>>,
}

impl AdapterContent {
    pub fn new() -> Self {
        Self {
            total: 0,
            longest_seq: 0,
            first_index: vec![Vec::new(); ADAPTERS.len()],
        }
    }
}

impl Default for AdapterContent {
    fn default() -> Self {
        Self::new()
    }
}

impl QcModule for AdapterContent {
    fn name(&self) -> &'static str {
        "adapter_content"
    }

    fn update(&mut self, _name: &[u8], seq: &[u8], _qual: &[u8]) {
        self.total += 1;
        self.longest_seq = self.longest_seq.max(seq.len());
        for (a, (_, adapter)) in ADAPTERS.iter().enumerate() {
            if seq.len() < adapter.len() {
                continue;
            }
            if let Some(idx) = seq.windows(adapter.len()).position(|w| w == *adapter) {
                let h = &mut self.first_index[a];
                if h.len() <= idx {
                    h.resize(idx + 1, 0);
                }
                h[idx] += 1;
            }
        }
    }

    fn merge(&mut self, other: &dyn QcModule) {
        let o = other
            .as_any()
            .downcast_ref::<AdapterContent>()
            .expect("merge type mismatch");
        self.total += o.total;
        self.longest_seq = self.longest_seq.max(o.longest_seq);
        for (a, oh) in o.first_index.iter().enumerate() {
            let h = &mut self.first_index[a];
            if h.len() < oh.len() {
                h.resize(oh.len(), 0);
            }
            for (s, &v) in h.iter_mut().zip(oh.iter()) {
                *s += v;
            }
        }
    }

    fn finalize(&self, out: &mut Vec<TidyRow>) {
        let m = "adapter_content";
        if self.longest_seq < LONGEST_ADAPTER || self.total == 0 {
            out.push(TidyRow::status(m, "PASS"));
            return;
        }
        let max_pos0 = self.longest_seq - LONGEST_ADAPTER; // 0-based, inclusive
        let mut worst = 0f64;
        for (a, (name, _)) in ADAPTERS.iter().enumerate() {
            let hist = &self.first_index[a];
            let mut running = 0u64;
            for pos in 0..=max_pos0 {
                running += hist.get(pos).copied().unwrap_or(0);
                let pct = running as f64 / self.total as f64 * 100.0;
                worst = worst.max(pct);
                out.push(TidyRow {
                    module: m,
                    label: Some((*name).to_string()),
                    position: Some((pos + 1) as i32),
                    metric: "pct".to_string(),
                    value: Some(pct),
                    value_str: None,
                });
            }
        }
        // FastQC: warn if any position exceeds 5% adapter, error above 10%.
        let status = if worst > 10.0 {
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
    fn adapter_content_cumulative_from_found_index() {
        let mut m = AdapterContent::new();
        // 2 reads with the Illumina Universal Adapter at index 4, 2 without.
        let with = b"ACGTAGATCGGAAGAGTTTTTTTTTTTTTTTT"; // adapter at index 4, len 32
        let without = b"ACGTACGTACGTACGTACGTACGTACGTACGT";
        m.update(b"", with, &[b'I'; 32]);
        m.update(b"", with, &[b'I'; 32]);
        m.update(b"", without, &[b'I'; 32]);
        m.update(b"", without, &[b'I'; 32]);
        let mut rows = Vec::new();
        m.finalize(&mut rows);
        let at = |pos: i32| {
            rows.iter()
                .find(|r| {
                    r.label.as_deref() == Some("Illumina Universal Adapter")
                        && r.position == Some(pos)
                })
                .and_then(|r| r.value)
                .unwrap()
        };
        // longest=32, longest_adapter=12 -> positions 1..=21. Found at index 4
        // -> cumulative 50% (2/4) from 1-based position 5 onward; 0 before.
        assert!((at(4) - 0.0).abs() < 1e-9);
        assert!((at(5) - 50.0).abs() < 1e-9);
        assert!((at(21) - 50.0).abs() < 1e-9);
    }
}
