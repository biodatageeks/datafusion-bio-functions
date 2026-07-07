use std::any::Any;
use std::collections::HashMap;

use crate::{QcModule, TidyRow};

/// FastQC uses the first 50bp of a read as its key, but ONLY for reads longer
/// than 75bp; reads <= 75bp are keyed on the whole sequence (same rule and
/// tracker as the duplication module — FastQC 0.12.1
/// `OverRepresentedSeqs.processSequence`: `if (seq.length() > 75) seq =
/// seq.substring(0, 50)`).
const KEY_PREFIX: usize = 50;
/// Reads longer than this are truncated to `KEY_PREFIX`; shorter reads are kept whole.
const TRUNCATE_ABOVE: usize = 75;
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
/// Source". The source is resolved by running each overrepresented sequence
/// through the built-in contaminant list (a faithful port of FastQC's
/// `ContaminentFinder.findContaminantHit`); a matching contaminant is reported
/// as `name (percentID% over Nbp)`, otherwise "No Hit".
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
        if seq.len() > TRUNCATE_ABOVE {
            &seq[..KEY_PREFIX]
        } else {
            seq
        }
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
            // FastQC OverRepresentedSeqs.contaminantHit(): run the observed
            // (<=50bp) sequence through the contaminant finder; null -> "No Hit",
            // otherwise the hit's toString() = "name (percentID% over Nbp)".
            let source = contaminant::find_contaminant_hit(seq)
                .map(|h| h.to_string())
                .unwrap_or_else(|| "No Hit".to_string());
            out.push(TidyRow {
                module: m,
                label: Some(label),
                position: None,
                metric: "source".to_string(),
                value: None,
                value_str: Some(source),
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

/// Faithful port of FastQC 0.12.1's contaminant matcher
/// (`uk.ac.babraham.FastQC.Sequence.Contaminant`). Used to fill the
/// "Possible Source" column of the overrepresented-sequences table.
pub(crate) mod contaminant {
    use std::fmt;
    use std::sync::OnceLock;

    /// A single named contaminant with its uppercased forward sequence and its
    /// reverse-complement, mirroring FastQC's `Contaminant`.
    struct Contaminant {
        name: String,
        forward: Vec<u8>,
        reverse: Vec<u8>,
    }

    /// FastQC's `ContaminantHit`. `direction` is not retained because it never
    /// affects `toString()` (the only thing we emit).
    pub struct ContaminantHit {
        name: String,
        length: usize,
        percent_id: i64,
    }

    impl fmt::Display for ContaminantHit {
        // ContaminantHit.toString(): name + " (" + percentID + "% over " + length + "bp)"
        fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
            write!(
                f,
                "{} ({}% over {}bp)",
                self.name, self.percent_id, self.length
            )
        }
    }

    impl Contaminant {
        /// Build a contaminant; returns `None` if the sequence contains any
        /// non-ACGT character. FastQC throws `IllegalArgumentException` here; we
        /// skip such entries at load time instead of panicking.
        fn new(name: &str, sequence: &str) -> Option<Contaminant> {
            let forward: Vec<u8> = sequence.to_ascii_uppercase().into_bytes();
            let mut reverse = vec![0u8; forward.len()];
            let n = forward.len();
            for (c, &b) in forward.iter().enumerate() {
                // Java switch: G<->C, A<->T; reverse[(len-1)-c] = complement.
                let comp = match b {
                    b'G' => b'C',
                    b'A' => b'T',
                    b'T' => b'A',
                    b'C' => b'G',
                    _ => return None,
                };
                reverse[(n - 1) - c] = comp;
            }
            Some(Contaminant {
                name: name.to_string(),
                forward,
                reverse,
            })
        }

        /// Port of `Contaminant.findMatch(query)`. `query` must already be
        /// uppercased (the caller does it once).
        fn find_match(&self, query: &[u8]) -> Option<ContaminantHit> {
            let qlen = query.len();

            // Special case for queries 8..20bp: an exact substring is a 100% hit.
            if (8..20).contains(&qlen) {
                if contains(&self.forward, query) {
                    return Some(ContaminantHit {
                        name: self.name.clone(),
                        length: qlen,
                        percent_id: 100,
                    });
                }
                if contains(&self.reverse, query) {
                    return Some(ContaminantHit {
                        name: self.name.clone(),
                        length: qlen,
                        percent_id: 100,
                    });
                }
            }

            // General case: slide the contaminant across the query allowing one
            // mismatch per run and requiring a >20bp match. Offsets may be
            // negative, so all offset arithmetic is done in i64 (as the Java does
            // with int, bounds-checking i+offset against [0, cb.len)).
            let flen = self.forward.len() as i64;
            let start_off = -(flen - 20);
            let end_off = qlen as i64 - 20;

            let mut best: Option<ContaminantHit> = None;

            let mut offset = start_off;
            while offset < end_off {
                if let Some(hit) = self.inner_find_match(&self.forward, query, offset)
                    && best.as_ref().is_none_or(|b| hit.length > b.length)
                {
                    best = Some(hit);
                }
                offset += 1;
            }

            let mut offset = start_off;
            while offset < end_off {
                if let Some(hit) = self.inner_find_match(&self.reverse, query, offset)
                    && best.as_ref().is_none_or(|b| hit.length > b.length)
                {
                    best = Some(hit);
                }
                offset += 1;
            }

            best
        }

        /// Port of the private `Contaminant.findMatch(ca, cb, offset, direction)`.
        /// `ca` is the contaminant (forward/reverse), `cb` is the query.
        fn inner_find_match(&self, ca: &[u8], cb: &[u8], offset: i64) -> Option<ContaminantHit> {
            let mut best: Option<ContaminantHit> = None;
            let mut mismatch_count: i64 = 0;
            let mut start: i64 = 0;
            let mut end: i64 = 0;
            let cblen = cb.len() as i64;

            for i in 0..ca.len() as i64 {
                let j = i + offset;
                if j < 0 {
                    start = i + 1;
                    continue;
                }
                if j >= cblen {
                    break;
                }

                if ca[i as usize] == cb[j as usize] {
                    end = i;
                } else {
                    mismatch_count += 1;
                    if mismatch_count > 1 {
                        // End of this run; record it if it's long enough.
                        let run = 1 + (end - start);
                        if run > 20 {
                            let id = ((run - (mismatch_count - 1)) * 100) / run;
                            consider(&mut best, &self.name, run, id);
                        }
                        start = i + 1;
                        end = i + 1;
                        mismatch_count = 0;
                    }
                }
            }

            // See if we ended with a match.
            let run = 1 + (end - start);
            if run > 20 {
                let id = ((run - mismatch_count) * 100) / run;
                consider(&mut best, &self.name, run, id);
            }

            best
        }
    }

    /// Record a candidate hit if it's longer than the current best (tie-break on
    /// higher percentID), matching the Java bookkeeping.
    fn consider(best: &mut Option<ContaminantHit>, name: &str, run: i64, id: i64) {
        let replace = match best {
            None => true,
            Some(b) => (b.length as i64) < run || (b.length as i64 == run && b.percent_id < id),
        };
        if replace {
            *best = Some(ContaminantHit {
                name: name.to_string(),
                length: run as usize,
                percent_id: id,
            });
        }
    }

    /// Byte-slice substring test (Java `String.contains`).
    fn contains(haystack: &[u8], needle: &[u8]) -> bool {
        if needle.is_empty() {
            return true;
        }
        if needle.len() > haystack.len() {
            return false;
        }
        haystack.windows(needle.len()).any(|w| w == needle)
    }

    /// The built-in FastQC 0.12.1 contaminant list, embedded at compile time
    /// (a verbatim copy of `Configuration/contaminant_list.txt`).
    const CONTAMINANT_LIST: &str = include_str!("contaminant_list.txt");

    fn contaminants() -> &'static [Contaminant] {
        static LIST: OnceLock<Vec<Contaminant>> = OnceLock::new();
        LIST.get_or_init(|| {
            let mut v = Vec::new();
            // Port of ContaminentFinder.makeContaminantList().
            for line in CONTAMINANT_LIST.lines() {
                if line.starts_with('#') {
                    continue; // Skip comments
                }
                if line.trim().is_empty() {
                    continue; // Skip blank lines
                }
                // Java: line.split("\\t+") -> exactly 2 fields (name, sequence).
                // Splitting on single tabs and dropping the empties produced by
                // runs of tabs reproduces the \t+ collapse for this file (no line
                // has a leading tab, so no leading-empty field is dropped).
                let sections: Vec<&str> = line.split('\t').filter(|s| !s.is_empty()).collect();
                if sections.len() != 2 {
                    continue; // Malformed line, ignore
                }
                if let Some(c) = Contaminant::new(sections[0], sections[1]) {
                    v.push(c);
                }
            }
            v
        })
    }

    /// Port of `ContaminentFinder.findContaminantHit(sequence)`: the best hit over
    /// all contaminants, keeping the longest (tie-break on higher percentID).
    pub fn find_contaminant_hit(sequence: &[u8]) -> Option<ContaminantHit> {
        let query = sequence.to_ascii_uppercase();
        let mut best: Option<ContaminantHit> = None;
        for c in contaminants() {
            let this = match c.find_match(&query) {
                Some(h) => h,
                None => continue,
            };
            let replace = match &best {
                None => true,
                Some(b) => {
                    this.length > b.length
                        || (this.length == b.length && this.percent_id > b.percent_id)
                }
            };
            if replace {
                best = Some(this);
            }
        }
        best
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

    #[test]
    fn keys_whole_read_up_to_75bp_then_truncates() {
        // FastQC only truncates reads > 75bp to their first 50bp. Two 60bp reads
        // sharing the first 50bp but differing afterwards are DISTINCT keys...
        let mut short = OverrepresentedSeqs::new();
        let a60 = vec![b'A'; 60];
        let mut b60 = vec![b'A'; 60];
        b60[55] = b'C';
        short.update(b"", &a60, b"");
        short.update(b"", &b60, b"");
        assert_eq!(short.seqs.len(), 2, "51-75bp reads must be keyed whole");

        // ...but two 80bp reads sharing the first 50bp collapse to one key.
        let mut long = OverrepresentedSeqs::new();
        let a80 = vec![b'A'; 80];
        let mut b80 = vec![b'A'; 80];
        b80[60] = b'C';
        long.update(b"", &a80, b"");
        long.update(b"", &b80, b"");
        assert_eq!(long.seqs.len(), 1, ">75bp reads must be truncated to 50bp");
    }

    // "Illumina Single End Adapter 1" from the embedded contaminant_list.txt.
    const KNOWN_ADAPTER: &[u8] = b"GATCGGAAGAGCTCGTATGCCGTCTTCTGCTTG"; // 33bp

    #[test]
    fn contaminant_exact_match_is_100pct() {
        let hit = contaminant::find_contaminant_hit(KNOWN_ADAPTER)
            .expect("exact contaminant should match");
        assert_eq!(
            hit.to_string(),
            "Illumina Single End Adapter 1 (100% over 33bp)"
        );
    }

    #[test]
    fn contaminant_no_hit_returns_none() {
        // Non-ACGT bytes can never equal an (ACGT-only) contaminant base, so no
        // run of >20 matches can form: guaranteed "No Hit".
        let seq = b"NNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNN"; // 40bp
        assert!(contaminant::find_contaminant_hit(seq).is_none());
    }

    #[test]
    fn contaminant_one_mismatch_still_hits() {
        // Introduce a single substitution into the 33bp adapter: one mismatch is
        // tolerated within a run, so it still hits at id = (33-1)*100/33 = 96%.
        let mut seq = KNOWN_ADAPTER.to_vec();
        seq[16] = if seq[16] == b'A' { b'C' } else { b'A' };
        let hit = contaminant::find_contaminant_hit(&seq)
            .expect(">=20bp near-match with one mismatch should hit");
        assert_eq!(
            hit.to_string(),
            "Illumina Single End Adapter 1 (96% over 33bp)"
        );
    }
}
