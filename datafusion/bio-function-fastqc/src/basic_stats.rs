use std::any::Any;

use crate::{QcModule, TidyRow};

/// FastQC "Basic Statistics": read count, length range, total bases, GC%.
#[derive(Debug, Default)]
pub struct BasicStats {
    n_seq: u64,
    total_bases: u64,
    gc_bases: u64,
    /// A+T+G+C only — FastQC's %GC denominator excludes N (and any other base).
    atgc_bases: u64,
    n_bases: u64,
    min_len: Option<u64>,
    max_len: u64,
}

impl BasicStats {
    pub fn new() -> Self {
        Self::default()
    }
}

impl QcModule for BasicStats {
    fn name(&self) -> &'static str {
        "basic_stats"
    }

    fn update(&mut self, _name: &[u8], seq: &[u8], _qual: &[u8]) {
        self.n_seq += 1;
        let len = seq.len() as u64;
        self.total_bases += len;
        self.max_len = self.max_len.max(len);
        self.min_len = Some(self.min_len.map_or(len, |m| m.min(len)));
        for &b in seq {
            match b {
                b'G' | b'g' | b'C' | b'c' => {
                    self.gc_bases += 1;
                    self.atgc_bases += 1;
                }
                b'A' | b'a' | b'T' | b't' => self.atgc_bases += 1,
                b'N' | b'n' => self.n_bases += 1,
                _ => {}
            }
        }
    }

    fn merge(&mut self, other: &dyn QcModule) {
        let o = other
            .as_any()
            .downcast_ref::<BasicStats>()
            .expect("merge type mismatch");
        self.n_seq += o.n_seq;
        self.total_bases += o.total_bases;
        self.gc_bases += o.gc_bases;
        self.atgc_bases += o.atgc_bases;
        self.n_bases += o.n_bases;
        self.max_len = self.max_len.max(o.max_len);
        self.min_len = match (self.min_len, o.min_len) {
            (Some(a), Some(b)) => Some(a.min(b)),
            (a, b) => a.or(b),
        };
    }

    fn finalize(&self, out: &mut Vec<TidyRow>) {
        let m = "basic_stats";
        out.push(TidyRow::num(m, "n_seq", self.n_seq as f64));
        out.push(TidyRow::num(m, "total_bases", self.total_bases as f64));
        out.push(TidyRow::num(m, "min_len", self.min_len.unwrap_or(0) as f64));
        out.push(TidyRow::num(m, "max_len", self.max_len as f64));
        // FastQC computes %GC as (G+C)*100/(A+T+G+C) with integer division and
        // N excluded. We emit the full-precision value over A/T/G/C (FastQC
        // displays its floor, e.g. 47.23 -> "47").
        let gc_pct = if self.atgc_bases > 0 {
            self.gc_bases as f64 / self.atgc_bases as f64 * 100.0
        } else {
            0.0
        };
        out.push(TidyRow::num(m, "gc_pct", gc_pct));
        // FastQC Basic Statistics never warns/fails.
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
    fn basic_stats_counts_and_gc() {
        let mut m = BasicStats::new();
        // 2 reads: "ACGT" (gc 2/4), "GGGGCA" (gc 5/6)
        m.update(b"", b"ACGT", b"IIII");
        m.update(b"", b"GGGGCA", b"IIIIII");
        let mut rows = Vec::new();
        m.finalize(&mut rows);
        let get = |metric: &str| {
            rows.iter()
                .find(|r| r.metric == metric)
                .and_then(|r| r.value)
                .unwrap()
        };
        assert_eq!(get("n_seq"), 2.0);
        assert_eq!(get("min_len"), 4.0);
        assert_eq!(get("max_len"), 6.0);
        assert_eq!(get("total_bases"), 10.0);
        assert!((get("gc_pct") - 70.0).abs() < 1e-9); // (2+5)/10*100
        assert!(
            rows.iter()
                .any(|r| r.metric == "status" && r.value_str.as_deref() == Some("PASS"))
        );
    }

    #[test]
    fn basic_stats_n_bases_and_length_range() {
        let mut m = BasicStats::new();
        m.update(b"", b"NNNN", b"IIII"); // 0 GC, len 4, all N
        m.update(b"", b"GCGCGC", b"IIIIII"); // 6 GC, len 6
        let mut rows = Vec::new();
        m.finalize(&mut rows);
        let get = |k: &str| {
            rows.iter()
                .find(|r| r.metric == k)
                .and_then(|r| r.value)
                .unwrap()
        };
        assert_eq!(get("min_len"), 4.0);
        assert_eq!(get("max_len"), 6.0);
        assert_eq!(get("total_bases"), 10.0); // full length, incl N
        // FastQC excludes N from %GC: (6 G/C) / (6 A+T+G+C) * 100 = 100%.
        assert!((get("gc_pct") - 100.0).abs() < 1e-9);
    }
}
