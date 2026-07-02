use std::any::Any;
use std::collections::HashMap;

use crate::{QcModule, TidyRow};

/// FastQC caps distinct-sequence tracking to bound memory.
const MAX_TRACKED: usize = 100_000;
/// Only the first N bases are used as the dedup key (FastQC uses 50).
const KEY_PREFIX: usize = 50;

/// FastQC "Sequence Duplication Levels".
#[derive(Debug, Default)]
pub struct DuplicationLevels {
    counts: HashMap<Vec<u8>, u64>,
    overflow_obs: u64,
}

impl DuplicationLevels {
    pub fn new() -> Self {
        Self::default()
    }

    fn key(seq: &[u8]) -> Vec<u8> {
        let n = seq.len().min(KEY_PREFIX);
        seq[..n].to_ascii_uppercase()
    }

    fn add_count(&mut self, key: Vec<u8>, n: u64) {
        if let Some(c) = self.counts.get_mut(&key) {
            *c += n;
        } else if self.counts.len() < MAX_TRACKED {
            self.counts.insert(key, n);
        } else {
            self.overflow_obs += n;
        }
    }

    fn level_bin(count: u64) -> &'static str {
        match count {
            0 => unreachable!(),
            1 => "1",
            2 => "2",
            3 => "3",
            4 => "4",
            5 => "5",
            6 => "6",
            7 => "7",
            8 => "8",
            9 => "9",
            10..=49 => ">10",
            50..=99 => ">50",
            100..=499 => ">100",
            500..=999 => ">500",
            1000..=4999 => ">1k",
            5000..=9999 => ">5k",
            _ => ">10k+",
        }
    }
}

const BINS: [&str; 16] = [
    "1", "2", "3", "4", "5", "6", "7", "8", "9", ">10", ">50", ">100", ">500", ">1k", ">5k",
    ">10k+",
];

impl QcModule for DuplicationLevels {
    fn name(&self) -> &'static str {
        "dup_levels"
    }

    fn update(&mut self, seq: &[u8], _qual: &[u8]) {
        let key = Self::key(seq);
        self.add_count(key, 1);
    }

    fn merge(&mut self, other: &dyn QcModule) {
        let o = other
            .as_any()
            .downcast_ref::<DuplicationLevels>()
            .expect("merge type mismatch");
        for (k, &c) in &o.counts {
            self.add_count(k.clone(), c);
        }
        self.overflow_obs += o.overflow_obs;
    }

    fn finalize(&self, out: &mut Vec<TidyRow>) {
        let m = "dup_levels";
        let distinct = self.counts.len() as u64;
        let tracked_obs: u64 = self.counts.values().sum();
        let total_obs = tracked_obs + self.overflow_obs;

        let mut level_counts: HashMap<&'static str, u64> = HashMap::new();
        for &c in self.counts.values() {
            *level_counts.entry(Self::level_bin(c)).or_insert(0) += 1;
        }
        for bin in BINS {
            let n = *level_counts.get(bin).unwrap_or(&0);
            let pct = if distinct > 0 {
                n as f64 / distinct as f64 * 100.0
            } else {
                0.0
            };
            out.push(TidyRow {
                module: m,
                label: Some(bin.to_string()),
                position: None,
                metric: "pct".to_string(),
                value: Some(pct),
                value_str: None,
            });
        }

        let pct_dup = if total_obs > 0 {
            (total_obs - distinct.min(total_obs)) as f64 / total_obs as f64 * 100.0
        } else {
            0.0
        };
        out.push(TidyRow::num(m, "pct_dup", pct_dup));

        let status = if pct_dup > 50.0 {
            "FAIL"
        } else if pct_dup > 20.0 {
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
    fn dup_levels_and_merge() {
        let mut a = DuplicationLevels::new();
        a.update(b"AAAA", b"IIII");
        a.update(b"AAAA", b"IIII");
        let mut b = DuplicationLevels::new();
        b.update(b"AAAA", b"IIII");
        b.update(b"CCCC", b"IIII");
        a.merge(&b);
        // AAAA x3, CCCC x1 -> 2 distinct, 4 observations.
        let mut rows = Vec::new();
        a.finalize(&mut rows);
        let pct_dup = rows
            .iter()
            .find(|r| r.metric == "pct_dup")
            .and_then(|r| r.value)
            .unwrap();
        assert!((pct_dup - 50.0).abs() < 1e-9); // (4-2)/4*100
    }
}
