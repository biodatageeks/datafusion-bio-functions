use std::any::Any;
use std::collections::HashMap;

use crate::{QcModule, TidyRow};

/// FastQC truncates a read to its first 50bp before using it as a key (same as
/// the duplication module — FastQC shares one tracker between them).
const KEY_PREFIX: usize = 50;
/// A sequence is "overrepresented" above this % of total (limits.txt: warn 0.1).
const WARN_PCT: f64 = 0.1;
/// error threshold (limits.txt: overrepresented error 1).
const FAIL_PCT: f64 = 1.0;
/// FastQC stops tracking new distinct sequences after this many uniques
/// (`OverRepresentedSeqs.OBSERVATION_CUTOFF`) to bound memory. Any sequence
/// above the 0.1% overrepresentation threshold is frequent enough to be seen
/// well before the freeze, so the reported set is unaffected.
const OBSERVATION_CUTOFF: usize = 100_000;

/// FastQC "Overrepresented sequences": lists 50bp-truncated sequences making up
/// more than 0.1% of the library, with count, percentage and a "Possible
/// Source". Contaminant/adapter source matching is not yet implemented, so the
/// source is always "No Hit" (matches libraries with no known contaminants).
#[derive(Debug, Default)]
pub struct OverrepresentedSeqs {
    seqs: HashMap<Vec<u8>, u64>,
    count: u64,
    /// once `seqs` reaches OBSERVATION_CUTOFF we stop adding new keys
    frozen: bool,
}

impl OverrepresentedSeqs {
    pub fn new() -> Self {
        Self::default()
    }

    fn key(seq: &[u8]) -> &[u8] {
        &seq[..seq.len().min(KEY_PREFIX)]
    }
}

impl QcModule for OverrepresentedSeqs {
    fn name(&self) -> &'static str {
        "overrepresented"
    }

    fn update(&mut self, _name: &[u8], seq: &[u8], _qual: &[u8]) {
        // Mirrors FastQC OverRepresentedSeqs.processSequence: count every read,
        // but stop adding NEW distinct sequences once frozen.
        self.count += 1;
        let key = Self::key(seq);
        if let Some(c) = self.seqs.get_mut(key) {
            *c += 1;
        } else if !self.frozen {
            self.seqs.insert(key.to_vec(), 1);
            if self.seqs.len() == OBSERVATION_CUTOFF {
                self.frozen = true;
            }
        }
    }

    fn merge(&mut self, other: &dyn QcModule) {
        let o = other
            .as_any()
            .downcast_ref::<OverrepresentedSeqs>()
            .expect("merge type mismatch");
        // Each partition table is bounded to OBSERVATION_CUTOFF keys, so the
        // union is bounded to CUTOFF * n_partitions.
        self.count += o.count;
        self.frozen |= o.frozen;
        for (k, &c) in &o.seqs {
            *self.seqs.entry(k.clone()).or_insert(0) += c;
        }
    }

    fn finalize(&self, out: &mut Vec<TidyRow>) {
        let m = "overrepresented";
        if self.count == 0 {
            out.push(TidyRow::status(m, "PASS"));
            return;
        }
        // Keep sequences above the warn threshold, sorted by count desc then
        // sequence (FastQC sorts by count desc; the tiebreak keeps us
        // deterministic across partition merges).
        let mut keepers: Vec<(&Vec<u8>, u64, f64)> = self
            .seqs
            .iter()
            .map(|(s, &c)| (s, c, c as f64 / self.count as f64 * 100.0))
            .filter(|&(_, _, pct)| pct > WARN_PCT)
            .collect();
        keepers.sort_by(|a, b| b.1.cmp(&a.1).then_with(|| a.0.cmp(b.0)));

        let mut worst = 0f64;
        for (seq, count, pct) in &keepers {
            worst = worst.max(*pct);
            let label = String::from_utf8_lossy(seq).into_owned();
            out.push(TidyRow {
                module: m,
                label: Some(label.clone()),
                position: None,
                metric: "count".to_string(),
                value: Some(*count as f64),
                value_str: None,
            });
            out.push(TidyRow {
                module: m,
                label: Some(label.clone()),
                position: None,
                metric: "pct".to_string(),
                value: Some(*pct),
                value_str: None,
            });
            out.push(TidyRow {
                module: m,
                label: Some(label),
                position: None,
                metric: "source".to_string(),
                value: None,
                value_str: Some("No Hit".to_string()),
            });
        }
        let status = if worst > FAIL_PCT {
            "FAIL"
        } else if worst > WARN_PCT {
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
    fn freezes_at_100k_unique_but_counts_existing() {
        let mut m = OverrepresentedSeqs::new();
        for i in 0..120_000u32 {
            m.update(b"", format!("{i:0>50}").as_bytes(), b"");
        }
        assert_eq!(m.seqs.len(), 100_000); // frozen at the cutoff
        assert_eq!(m.count, 120_000); // total still counted
        // A sequence already in the table keeps accumulating after the freeze.
        let existing = format!("{:0>50}", 0u32);
        m.update(b"", existing.as_bytes(), b"");
        assert_eq!(m.seqs.len(), 100_000);
        assert_eq!(*m.seqs.get(existing.as_bytes()).unwrap(), 2);
    }

    #[test]
    fn overrepresented_lists_frequent_sequences() {
        let mut m = OverrepresentedSeqs::new();
        for _ in 0..10 {
            m.update(b"", b"AAAACCCC", b"IIIIIIII"); // 10/11 ~= 90.9%
        }
        m.update(b"", b"GGGGTTTT", b"IIIIIIII"); // 1/11 ~= 9.1%
        let mut rows = Vec::new();
        m.finalize(&mut rows);
        let count_of = |label: &str| {
            rows.iter()
                .find(|r| r.label.as_deref() == Some(label) && r.metric == "count")
                .and_then(|r| r.value)
        };
        assert_eq!(count_of("AAAACCCC"), Some(10.0));
        assert_eq!(count_of("GGGGTTTT"), Some(1.0));
        // > 1% -> FAIL
        assert!(
            rows.iter()
                .any(|r| r.metric == "status" && r.value_str.as_deref() == Some("FAIL"))
        );
    }
}
