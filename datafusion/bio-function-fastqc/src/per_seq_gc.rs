use std::any::Any;
use std::collections::HashMap;

use crate::{QcModule, TidyRow};

/// FastQC "Per Sequence GC Content".
///
/// Matches FastQC exactly: raw counts are accumulated per `(read_length,
/// gc_count)` (integer, so the cross-partition merge stays associative), then
/// `finalize` interpolates each read's real GC% linearly across the two nearest
/// integer 0..=100 bins — FastQC's smoothing that yields fractional counts.
#[derive(Debug, Default)]
pub struct PerSeqGc {
    /// read length -> counts indexed by gc_count (0..=length)
    by_len: HashMap<usize, Vec<u64>>,
}

impl PerSeqGc {
    pub fn new() -> Self {
        Self::default()
    }

    /// FastQC truncates each read before the GC calculation to bound the number
    /// of cached length-models: reads >1000 bp are cut to a multiple of 1000,
    /// reads >100 bp to a multiple of 100 (so a 101 bp read uses its first
    /// 100 bp). This is specific to Per Sequence GC Content.
    fn truncated_len(len: usize) -> usize {
        if len > 1000 {
            (len / 1000) * 1000
        } else if len > 100 {
            (len / 100) * 100
        } else {
            len
        }
    }
}

impl QcModule for PerSeqGc {
    fn name(&self) -> &'static str {
        "per_seq_gc"
    }

    fn update(&mut self, _name: &[u8], seq: &[u8], _qual: &[u8]) {
        // FastQC truncates the read first, then counts G/C over the truncated
        // portion and models against the truncated length.
        let len = Self::truncated_len(seq.len());
        if len == 0 {
            return;
        }
        let gc = seq[..len]
            .iter()
            .filter(|&&b| matches!(b, b'G' | b'g' | b'C' | b'c'))
            .count();
        let counts = self
            .by_len
            .entry(len)
            .or_insert_with(|| vec![0u64; len + 1]);
        counts[gc] += 1;
    }

    fn merge(&mut self, other: &dyn QcModule) {
        let o = other
            .as_any()
            .downcast_ref::<PerSeqGc>()
            .expect("merge type mismatch");
        for (len, counts) in &o.by_len {
            let e = self
                .by_len
                .entry(*len)
                .or_insert_with(|| vec![0u64; len + 1]);
            for (i, &c) in counts.iter().enumerate() {
                e[i] += c;
            }
        }
    }

    fn finalize(&self, out: &mut Vec<TidyRow>) {
        let m = "per_seq_gc";
        // Replicate FastQC's GCModel exactly: for a read of length L, a GC count
        // `pos` claims the percentage bins round((pos-0.5)*100/L)..=round((pos+
        // 0.5)*100/L), each weighted 1/claimingCounts[p], where claimingCounts[p]
        // is how many gc-counts (0..=L) map onto bin p. Lengths are processed in
        // sorted order so the float accumulation is deterministic regardless of
        // how partitions merged.
        let mut bins = [0f64; 101];
        let mut lens: Vec<usize> = self.by_len.keys().copied().collect();
        lens.sort_unstable();
        for len in lens {
            if len == 0 {
                continue;
            }
            let lenf = len as f64;
            // Java Math.round == round-half-up; percentages are non-negative so
            // f64::round (half away from zero) is equivalent here.
            let bounds = |pos: usize| -> (usize, usize) {
                let mut low = pos as f64 - 0.5;
                let mut high = pos as f64 + 0.5;
                if low < 0.0 {
                    low = 0.0;
                }
                if high < 0.0 {
                    high = 0.0;
                }
                if high > lenf {
                    high = lenf;
                }
                if low > lenf {
                    low = lenf;
                }
                let lp = (low * 100.0 / lenf).round() as usize;
                let hp = (high * 100.0 / lenf).round() as usize;
                (lp, hp)
            };
            // Pass 1: claiming counts.
            let mut claiming = [0i64; 101];
            for pos in 0..=len {
                let (lp, hp) = bounds(pos);
                for c in claiming.iter_mut().take(hp + 1).skip(lp) {
                    *c += 1;
                }
            }
            // Pass 2: distribute each read's weight over its claimed bins.
            let counts = &self.by_len[&len];
            for (gc, &c) in counts.iter().enumerate() {
                if c == 0 {
                    continue;
                }
                let (lp, hp) = bounds(gc);
                for p in lp..=hp {
                    bins[p] += c as f64 / claiming[p] as f64;
                }
            }
        }
        for (g, &c) in bins.iter().enumerate() {
            out.push(TidyRow {
                module: m,
                label: None,
                position: Some(g as i32),
                metric: "count".to_string(),
                value: Some(c),
                value_str: None,
            });
        }
        // Phase-1 status: PASS. Exact FastQC theoretical-distribution status
        // (warn/fail on deviation from the modelled normal) is a follow-up.
        out.push(TidyRow::status(m, "PASS"));
    }

    fn as_any(&self) -> &dyn Any {
        self
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn counts(rows: &[TidyRow]) -> Vec<&TidyRow> {
        rows.iter().filter(|r| r.metric == "count").collect()
    }

    #[test]
    fn per_seq_gc_emits_full_axis_with_peak_near_true_gc() {
        // 50%-GC reads (len 8) -> mass concentrates near bin 50.
        let mut m = PerSeqGc::new();
        for _ in 0..10 {
            m.update(b"", b"ATGCATGC", b"IIIIIIII");
        }
        let mut rows = Vec::new();
        m.finalize(&mut rows);
        let c = counts(&rows);
        assert_eq!(c.len(), 101); // bins 0..=100
        assert!(rows.iter().any(|r| r.metric == "status"));
        let peak = c
            .iter()
            .max_by(|a, b| a.value.partial_cmp(&b.value).unwrap())
            .unwrap();
        assert!((40..=60).contains(&peak.position.unwrap()));
    }

    #[test]
    fn per_seq_gc_truncates_reads_over_100bp() {
        // A 101bp read is truncated to its first 100bp: a G/C at position 101
        // is dropped, so it must match a 100bp read with the same first 100bp.
        let mut long_seq = vec![b'A'; 100];
        long_seq.push(b'G'); // 101st base, must be ignored
        let mut a = PerSeqGc::new();
        a.update(b"", &long_seq, &[b'I'; 101]);
        let mut b = PerSeqGc::new();
        b.update(b"", &[b'A'; 100], &[b'I'; 100]);
        let (mut ra, mut rb) = (Vec::new(), Vec::new());
        a.finalize(&mut ra);
        b.finalize(&mut rb);
        assert_eq!(ra, rb);
    }

    #[test]
    fn per_seq_gc_merge_equals_combined() {
        // merge(a, b) must equal accumulating everything into one set, and
        // finalize must be deterministic -> underpins partition invariance.
        let g1: &[&[u8]] = &[b"ATGCATGC", b"GGGGCCCC"];
        let g2: &[&[u8]] = &[b"ATATATAT", b"ATGCGCAT"];
        let (mut a, mut b, mut ab) = (PerSeqGc::new(), PerSeqGc::new(), PerSeqGc::new());
        for s in g1 {
            a.update(b"", s, b"IIIIIIII");
            ab.update(b"", s, b"IIIIIIII");
        }
        for s in g2 {
            b.update(b"", s, b"IIIIIIII");
            ab.update(b"", s, b"IIIIIIII");
        }
        a.merge(&b);
        let (mut r_merged, mut r_combined) = (Vec::new(), Vec::new());
        a.finalize(&mut r_merged);
        ab.finalize(&mut r_combined);
        assert_eq!(r_merged, r_combined);
    }
}
