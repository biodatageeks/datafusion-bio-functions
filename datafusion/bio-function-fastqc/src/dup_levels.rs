use std::any::Any;
use std::collections::HashMap;

use crate::{QcModule, TidyRow};

/// FastQC uses the first 50bp of a read as its duplication key, but ONLY for
/// reads longer than 75bp; reads <= 75bp are keyed on the whole sequence
/// (FastQC 0.12.1 `SequenceDeDuplication`/`OverRepresentedSeqs`:
/// `if (seq.length() > 75) seq = seq.substring(0, 50)`). Truncating shorter
/// reads to 50bp would merge distinct 51-75bp reads that share a prefix.
const KEY_PREFIX: usize = 50;
/// Reads longer than this are truncated to `KEY_PREFIX`; shorter reads are kept whole.
const TRUNCATE_ABOVE: usize = 75;

/// FastQC "Sequence Duplication Levels" (matches FastQC 0.12.1 exactly for
/// libraries with <= 100k distinct sequences, i.e. the range where FastQC does
/// not freeze its observation table and applies no count correction).
///
/// We do NOT freeze at FastQC's 100k-unique cap: keeping every key makes the
/// parallel merge exact and deterministic. `getCorrectedCount` is implemented
/// faithfully but is the identity while `count_at_unique_limit == count`, which
/// always holds here — so beyond 100k distinct sequences we are more accurate
/// than (and diverge from) FastQC's sampled estimate.
/// FastQC stops tracking NEW distinct sequences after this many uniques, to
/// bound memory (`OverRepresentedSeqs.OBSERVATION_CUTOFF`).
const OBSERVATION_CUTOFF: usize = 100_000;

#[derive(Debug, Default)]
pub struct DuplicationLevels {
    /// truncated sequence -> number of observations
    seqs: HashMap<Vec<u8>, u64>,
    /// total sequences seen
    count: u64,
    /// once `seqs` reaches OBSERVATION_CUTOFF we stop adding new keys
    frozen: bool,
    /// total observations at the moment we froze (feeds `corrected_count`)
    count_at_unique_limit: u64,
}

/// FastQC's 16 duplication-level bins.
const LABELS: [&str; 16] = [
    "1", "2", "3", "4", "5", "6", "7", "8", "9", ">10", ">50", ">100", ">500", ">1k", ">5k",
    ">10k+",
];

/// FastQC's dupSlot mapping: slot = f(dupLevel-1).
fn dup_slot(dup_level: u64) -> usize {
    let s = dup_level as i64 - 1;
    if !(0..=9999).contains(&s) {
        15
    } else if s > 4999 {
        14
    } else if s > 999 {
        13
    } else if s > 499 {
        12
    } else if s > 99 {
        11
    } else if s > 49 {
        10
    } else if s > 9 {
        9
    } else {
        s as usize
    }
}

/// FastQC's corrected-count estimator. Returns the observed count unchanged
/// whenever the full library was tracked (`count_at_unique_limit == total`).
fn corrected_count(count_at_limit: u64, total: u64, dup_level: u64, num_obs: u64) -> f64 {
    if count_at_limit == total {
        return num_obs as f64;
    }
    if total - num_obs < count_at_limit {
        return num_obs as f64;
    }
    let mut p_not_seeing = 1.0f64;
    let limit_of_caring = 1.0 - (num_obs as f64 / (num_obs as f64 + 0.01));
    for i in 0..count_at_limit {
        p_not_seeing *= ((total - i) as f64 - dup_level as f64) / (total - i) as f64;
        if p_not_seeing < limit_of_caring {
            p_not_seeing = 0.0;
            break;
        }
    }
    let p_seeing = 1.0 - p_not_seeing;
    num_obs as f64 / p_seeing
}

impl DuplicationLevels {
    pub fn new() -> Self {
        Self::default()
    }

    fn key(seq: &[u8]) -> &[u8] {
        if seq.len() > TRUNCATE_ABOVE {
            &seq[..KEY_PREFIX]
        } else {
            seq
        }
    }
}

impl QcModule for DuplicationLevels {
    fn name(&self) -> &'static str {
        "dup_levels"
    }

    fn update(&mut self, _name: &[u8], seq: &[u8], _qual: &[u8]) {
        // Mirrors FastQC OverRepresentedSeqs.processSequence: count every read,
        // but stop adding NEW distinct sequences once frozen.
        self.count += 1;
        let key = Self::key(seq);
        if let Some(c) = self.seqs.get_mut(key) {
            *c += 1;
            if !self.frozen {
                self.count_at_unique_limit = self.count;
            }
        } else if !self.frozen {
            self.seqs.insert(key.to_vec(), 1);
            self.count_at_unique_limit = self.count;
            if self.seqs.len() == OBSERVATION_CUTOFF {
                self.frozen = true;
            }
        }
    }

    fn merge(&mut self, other: &dyn QcModule) {
        let o = other
            .as_any()
            .downcast_ref::<DuplicationLevels>()
            .expect("merge type mismatch");
        // Union the per-partition tables (each already bounded to
        // OBSERVATION_CUTOFF keys, so the merged table is bounded to
        // CUTOFF * n_partitions). A single-partition scan reproduces FastQC's
        // in-order freeze exactly; multi-partition is a bounded approximation.
        self.count += o.count;
        self.count_at_unique_limit += o.count_at_unique_limit;
        self.frozen |= o.frozen;
        for (k, &c) in &o.seqs {
            *self.seqs.entry(k.clone()).or_insert(0) += c;
        }
    }

    fn finalize(&self, out: &mut Vec<TidyRow>) {
        let m = "dup_levels";
        // Observations counted while the table was still growing. Equals the
        // total whenever we never froze (< 100k distinct), so small libraries
        // are unchanged and exact; past the freeze this drives FastQC's
        // corrected-count estimator.
        let count_at_limit = self.count_at_unique_limit;

        // Collate: duplication level -> number of distinct sequences at it.
        let mut collated: HashMap<u64, u64> = HashMap::new();
        for &c in self.seqs.values() {
            *collated.entry(c).or_insert(0) += 1;
        }

        // Correct + accumulate into raw / dedup totals and the 16 bins,
        // iterating dup levels in sorted order for deterministic float sums.
        let mut levels: Vec<u64> = collated.keys().copied().collect();
        levels.sort_unstable();
        let mut total_percentages = [0f64; 16];
        let mut dedup_total = 0f64;
        let mut raw_total = 0f64;
        for dup_level in levels {
            let num_obs = collated[&dup_level];
            let corrected = corrected_count(count_at_limit, self.count, dup_level, num_obs);
            dedup_total += corrected;
            raw_total += corrected * dup_level as f64;
            total_percentages[dup_slot(dup_level)] += corrected * dup_level as f64;
        }

        for (i, label) in LABELS.iter().enumerate() {
            let pct = if raw_total > 0.0 {
                total_percentages[i] / raw_total * 100.0
            } else {
                0.0
            };
            out.push(TidyRow {
                module: m,
                label: Some(label.to_string()),
                position: None,
                metric: "pct".to_string(),
                value: Some(pct),
                value_str: None,
            });
        }

        // FastQC's headline number: % of sequences remaining after dedup.
        let percent_different = if raw_total == 0.0 {
            100.0
        } else {
            dedup_total / raw_total * 100.0
        };
        out.push(TidyRow::num(m, "total_dedup_pct", percent_different));

        // FastQC status thresholds (Configuration/limits.txt): warn < 70,
        // error < 50, on the % remaining after deduplication.
        let status = if percent_different < 50.0 {
            "FAIL"
        } else if percent_different < 70.0 {
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

    fn pct_at(rows: &[TidyRow], label: &str) -> f64 {
        rows.iter()
            .find(|r| r.label.as_deref() == Some(label) && r.metric == "pct")
            .and_then(|r| r.value)
            .unwrap()
    }

    fn scalar(rows: &[TidyRow], metric: &str) -> f64 {
        rows.iter()
            .find(|r| r.metric == metric)
            .and_then(|r| r.value)
            .unwrap()
    }

    #[test]
    fn freezes_at_100k_unique_but_keeps_counting() {
        // FastQC stops tracking NEW sequences after 100k unique, but keeps
        // counting total observations. This bounds memory on diverse libraries.
        let mut m = DuplicationLevels::new();
        for i in 0..120_000u32 {
            m.update(b"", format!("{i:0>50}").as_bytes(), b"");
        }
        assert_eq!(m.seqs.len(), 100_000); // frozen at the observation cutoff
        assert_eq!(m.count, 120_000); // total keeps growing past the freeze
        assert_eq!(m.count_at_unique_limit, 100_000); // count when we froze
    }

    #[test]
    fn freeze_correction_matches_fastqc() {
        // 100k unique singletons then 50k copies of one already-tracked seq.
        // FastQC 0.12.1 reports exactly 66.667% remaining after dedup and puts
        // the duplicated sequence in the ">10k+" bin — we must reproduce that
        // via the frozen-table correction (verified against a real golden run).
        let mut m = DuplicationLevels::new();
        for i in 0..100_000u32 {
            m.update(b"", format!("{i:0>50}").as_bytes(), b"");
        }
        let dup = format!("{:0>50}", 0u32);
        for _ in 0..50_000 {
            m.update(b"", dup.as_bytes(), b"");
        }
        let mut rows = Vec::new();
        m.finalize(&mut rows);
        assert!((scalar(&rows, "total_dedup_pct") - 66.6667).abs() < 1e-3);
        assert!((pct_at(&rows, "1") - 66.666).abs() < 1e-2);
        assert!((pct_at(&rows, ">10k+") - 33.334).abs() < 1e-2);
    }

    #[test]
    fn dup_levels_percentages_and_dedup() {
        // AAAA x3, CCCC x1 -> collated {3:1, 1:1}; rawTotal = 3+1 = 4.
        let mut a = DuplicationLevels::new();
        a.update(b"", b"AAAA", b"IIII");
        a.update(b"", b"AAAA", b"IIII");
        let mut b = DuplicationLevels::new();
        b.update(b"", b"AAAA", b"IIII");
        b.update(b"", b"CCCC", b"IIII");
        a.merge(&b);
        let mut rows = Vec::new();
        a.finalize(&mut rows);
        // bin "3" gets the AAAA reads: 3/4*100 = 75; bin "1" gets CCCC: 1/4*100 = 25.
        assert!((pct_at(&rows, "3") - 75.0).abs() < 1e-9);
        assert!((pct_at(&rows, "1") - 25.0).abs() < 1e-9);
        // 2 distinct / 4 total -> 50% remaining after dedup.
        assert!((scalar(&rows, "total_dedup_pct") - 50.0).abs() < 1e-9);
    }

    #[test]
    fn dup_levels_all_unique_is_100pct_dedup() {
        let mut m = DuplicationLevels::new();
        for s in [b"AAAA".as_slice(), b"CCCC".as_slice(), b"GGGG".as_slice()] {
            m.update(b"", s, b"IIII");
        }
        let mut rows = Vec::new();
        m.finalize(&mut rows);
        assert!((pct_at(&rows, "1") - 100.0).abs() < 1e-9);
        assert!((scalar(&rows, "total_dedup_pct") - 100.0).abs() < 1e-9);
    }

    #[test]
    fn dup_levels_keeps_whole_read_up_to_75bp() {
        // FastQC only truncates reads > 75bp. Two 60bp reads that share the first
        // 50bp but differ afterwards are DISTINCT (keyed on the whole read), not
        // duplicates. (Regression guard for the 51-75bp truncation bug.)
        let read_a = vec![b'A'; 60];
        let mut read_b = vec![b'A'; 60];
        read_b[55] = b'C';
        let qual = [b'I'; 60];
        let mut m = DuplicationLevels::new();
        m.update(b"", &read_a, &qual);
        m.update(b"", &read_b, &qual);
        let mut rows = Vec::new();
        m.finalize(&mut rows);
        // 2 distinct singletons -> bin "1" = 100%, 100% remaining after dedup.
        assert!((pct_at(&rows, "1") - 100.0).abs() < 1e-9);
        assert!((scalar(&rows, "total_dedup_pct") - 100.0).abs() < 1e-9);
    }

    #[test]
    fn dup_levels_truncates_reads_over_75bp_to_50bp() {
        // Reads > 75bp are keyed on their first 50bp, so two 80bp reads sharing the
        // first 50bp but differing afterwards collapse to one duplicated key.
        let read_a = vec![b'A'; 80];
        let mut read_b = vec![b'A'; 80];
        read_b[60] = b'C';
        let qual = [b'I'; 80];
        let mut m = DuplicationLevels::new();
        m.update(b"", &read_a, &qual);
        m.update(b"", &read_b, &qual);
        let mut rows = Vec::new();
        m.finalize(&mut rows);
        // 1 distinct seen twice -> bin "2" = 100%, 50% remaining after dedup.
        assert!((pct_at(&rows, "2") - 100.0).abs() < 1e-9);
        assert!((scalar(&rows, "total_dedup_pct") - 50.0).abs() < 1e-9);
    }
}
