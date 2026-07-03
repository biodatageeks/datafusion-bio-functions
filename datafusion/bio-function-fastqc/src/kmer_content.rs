use std::any::Any;
use std::collections::HashMap;

use crate::{QcModule, TidyRow};

/// FastQC default Kmer length (Configuration allows 2..=10; MIN==MAX==7).
const KMER_SIZE: usize = 7;
/// FastQC extrapolates the reported Count from the 2% sample by this factor.
const COUNT_DISPLAY_FACTOR: u64 = 5;

#[derive(Debug, Clone, Default)]
struct Kmer {
    count: u64,
    /// per start-position occurrence counts (index = 0-based start of the k-mer)
    positions: Vec<u64>,
}

impl Kmer {
    fn incr(&mut self, pos: usize) {
        self.count += 1;
        if self.positions.len() <= pos {
            self.positions.resize(pos + 1, 0);
        }
        self.positions[pos] += 1;
    }
}

/// FastQC "Kmer Content": faithfully replicates FastQC 0.12.1 `KmerContent`:
/// 2% read sampling (every 50th read), 500 bp truncation, per-position obs/exp
/// with a binomial (×4^k Bonferroni) p-value, keeping k-mers whose best position
/// clears `p < 0.01 && obs/exp > 5`, reported top-20 by max obs/exp.
///
/// The every-50th sampling is file-order dependent, so exact FastQC parity holds
/// only for a single-partition scan (plain `.fastq`/`.gz`).
#[derive(Debug, Default)]
pub struct KmerContent {
    skip_count: u64,
    longest_sequence: usize,
    /// total k-mers observed at each start position (INCLUDING N-containing ones)
    total_at_position: Vec<u64>,
    kmers: HashMap<[u8; KMER_SIZE], Kmer>,
}

impl KmerContent {
    pub fn new() -> Self {
        Self::default()
    }

    #[cfg(test)]
    fn total_kmer_count(&self) -> u64 {
        self.total_at_position.iter().sum()
    }

    #[cfg(test)]
    fn kmer_count(&self, kmer: &[u8]) -> u64 {
        <[u8; KMER_SIZE]>::try_from(kmer)
            .ok()
            .and_then(|k| self.kmers.get(&k))
            .map(|e| e.count)
            .unwrap_or(0)
    }
}

impl QcModule for KmerContent {
    fn name(&self) -> &'static str {
        "kmer_content"
    }

    fn update(&mut self, _name: &[u8], seq: &[u8], _qual: &[u8]) {
        // 2% sampling: only every 50th read is processed (FastQC skipCount % 50).
        self.skip_count += 1;
        if !self.skip_count.is_multiple_of(50) {
            return;
        }
        // Truncate reads longer than 500 bp (FastQC memory guard).
        let seq = if seq.len() > 500 { &seq[..500] } else { seq };
        if seq.len() > self.longest_sequence {
            self.longest_sequence = seq.len();
        }
        if seq.len() < KMER_SIZE {
            return;
        }
        for i in 0..=(seq.len() - KMER_SIZE) {
            // addKmerCount runs BEFORE the N check, so totals include N k-mers.
            if self.total_at_position.len() <= i {
                self.total_at_position.resize(i + 1, 0);
            }
            self.total_at_position[i] += 1;

            let window = &seq[i..i + KMER_SIZE];
            if window.contains(&b'N') {
                continue;
            }
            let key: [u8; KMER_SIZE] = window.try_into().unwrap();
            self.kmers.entry(key).or_default().incr(i);
        }
    }

    fn merge(&mut self, other: &dyn QcModule) {
        let o = other
            .as_any()
            .downcast_ref::<KmerContent>()
            .expect("merge type mismatch");
        self.skip_count += o.skip_count;
        self.longest_sequence = self.longest_sequence.max(o.longest_sequence);
        if self.total_at_position.len() < o.total_at_position.len() {
            self.total_at_position.resize(o.total_at_position.len(), 0);
        }
        for (i, &v) in o.total_at_position.iter().enumerate() {
            self.total_at_position[i] += v;
        }
        for (key, ok) in &o.kmers {
            let e = self.kmers.entry(*key).or_default();
            e.count += ok.count;
            if e.positions.len() < ok.positions.len() {
                e.positions.resize(ok.positions.len(), 0);
            }
            for (i, &v) in ok.positions.iter().enumerate() {
                e.positions[i] += v;
            }
        }
    }

    fn finalize(&self, out: &mut Vec<TidyRow>) {
        let m = "kmer_content";
        if self.kmers.is_empty() || self.longest_sequence < KMER_SIZE {
            out.push(TidyRow::status(m, "PASS"));
            return;
        }
        let groups_len = self.longest_sequence - KMER_SIZE + 1;
        let total_kmer_count: u64 = self.total_at_position.iter().sum();
        let bonferroni = 4f64.powi(KMER_SIZE as i32);

        let mut enriched: Vec<Enriched> = Vec::new();

        for (key, k) in &self.kmers {
            let expected_proportion = k.count as f64 / total_kmer_count as f64;
            let mut max_obs_exp = 0f64;
            let mut max_pos = 0usize;
            let mut lowest_p = 1f64;
            for g in 0..groups_len {
                let total_group_count = self.total_at_position.get(g).copied().unwrap_or(0);
                let total_group_hits = k.positions.get(g).copied().unwrap_or(0);
                let predicted = expected_proportion * total_group_count as f64;
                let obs_exp = if predicted > 0.0 {
                    total_group_hits as f64 / predicted
                } else {
                    0.0
                };
                if obs_exp > max_obs_exp {
                    max_obs_exp = obs_exp;
                    max_pos = g + 1;
                }
                if total_group_hits as f64 > predicted {
                    // FastQC: (1 - binomialCDF(hits)) * 4^k, computed via the
                    // upper tail to preserve precision on tiny p-values.
                    let p = binom_sf(total_group_count, expected_proportion, total_group_hits)
                        * bonferroni;
                    if p < 0.01 && obs_exp > 5.0 && p < lowest_p {
                        lowest_p = p;
                    }
                }
            }
            if lowest_p < 0.01 {
                enriched.push(Enriched {
                    seq: *key,
                    count: k.count,
                    lowest_p,
                    max_obs_exp,
                    max_pos: max_pos.max(1),
                });
            }
        }

        // Sort descending by max obs/exp (FastQC Kmer.compareTo), keep top 20.
        let enriched = rank_enriched(enriched);

        // Status from the top k-mer's p-value: -log10(p) vs kmer warn=2 / error=5.
        let status = match enriched.first() {
            Some(top) => {
                let s = -top.lowest_p.log10();
                if s > 5.0 {
                    "FAIL"
                } else if s > 2.0 {
                    "WARN"
                } else {
                    "PASS"
                }
            }
            None => "PASS",
        };

        for e in &enriched {
            let label = String::from_utf8_lossy(&e.seq).into_owned();
            let row = |metric: &str, value: f64| TidyRow {
                module: m,
                label: Some(label.clone()),
                position: None,
                metric: metric.to_string(),
                value: Some(value),
                value_str: None,
            };
            // Count matches FastQC's reported (extrapolated) value = raw × 5.
            out.push(row("count", (e.count * COUNT_DISPLAY_FACTOR) as f64));
            out.push(row("pvalue", e.lowest_p));
            out.push(row("obs_exp_max", e.max_obs_exp));
            out.push(row("max_position", e.max_pos as f64));
        }
        out.push(TidyRow::status(m, status));
    }

    fn as_any(&self) -> &dyn Any {
        self
    }
}

struct Enriched {
    seq: [u8; KMER_SIZE],
    count: u64,
    lowest_p: f64,
    max_obs_exp: f64,
    max_pos: usize, // 1-based
}

/// Rank enriched k-mers descending by max obs/exp (FastQC `Kmer.compareTo`) and
/// keep the top 20.
fn rank_enriched(mut enriched: Vec<Enriched>) -> Vec<Enriched> {
    // Break obs/exp ties by sequence ascending so the reported top-20 is
    // deterministic run-to-run (on real data many k-mers tie at the max obs/exp
    // and HashMap iteration order would otherwise decide which 20 are kept).
    enriched.sort_by(|a, b| {
        b.max_obs_exp
            .partial_cmp(&a.max_obs_exp)
            .unwrap()
            .then_with(|| a.seq.cmp(&b.seq))
    });
    enriched.truncate(20);
    enriched
}

/// P(X <= k) for X ~ Binomial(n, p), summing exact PMF terms in log space.
/// Mirrors Apache-commons `BinomialDistribution.cumulativeProbability`. Retained
/// as the lower-tail reference; the enrichment p-value uses [`binom_sf`] instead
/// (see its doc for why the upper tail is computed directly).
#[cfg(test)]
fn binom_cdf(n: u64, p: f64, k: u64) -> f64 {
    if k >= n {
        return 1.0;
    }
    if p <= 0.0 {
        return 1.0;
    }
    let ln_p = p.ln();
    let ln_q = (1.0 - p).ln();
    let mut cum = 0.0f64;
    for i in 0..=k {
        let ln_choose =
            ln_gamma((n + 1) as f64) - ln_gamma((i + 1) as f64) - ln_gamma((n - i + 1) as f64);
        cum += (ln_choose + (i as f64) * ln_p + ((n - i) as f64) * ln_q).exp();
    }
    cum.min(1.0)
}

/// P(X > k) = P(X >= k+1) for X ~ Binomial(n, p), summed directly from the
/// upper tail. FastQC computes `1 - BinomialDistribution.cumulativeProbability`,
/// but forming `1 - CDF` loses all precision when the tail is tiny (the CDF is
/// stored as ~1 and the subtraction keeps only absolute precision). Summing the
/// upper-tail PMF terms keeps full relative precision, matching FastQC's value.
fn binom_sf(n: u64, p: f64, k: u64) -> f64 {
    if k >= n {
        return 0.0;
    }
    if p <= 0.0 {
        return 0.0;
    }
    let ln_p = p.ln();
    let ln_q = (1.0 - p).ln();
    let mean = n as f64 * p;
    let mut sum = 0.0f64;
    for i in (k + 1)..=n {
        let ln_choose =
            ln_gamma((n + 1) as f64) - ln_gamma((i + 1) as f64) - ln_gamma((n - i + 1) as f64);
        let term = (ln_choose + (i as f64) * ln_p + ((n - i) as f64) * ln_q).exp();
        sum += term;
        // Past the mode the PMF decays monotonically; stop once terms no longer
        // affect the sum (bounds cost for large n while staying exact).
        if i as f64 > mean && term < sum * 1e-17 {
            break;
        }
    }
    sum.min(1.0)
}

/// Lanczos approximation for ln Gamma.
fn ln_gamma(x: f64) -> f64 {
    const G: [f64; 8] = [
        676.5203681218851,
        -1259.1392167224028,
        771.323428777653,
        -176.615029162140,
        12.507343278687,
        -0.13857109526572,
        9.984369578019572e-6,
        1.5056327351493116e-7,
    ];
    if x < 0.5 {
        std::f64::consts::PI / ((std::f64::consts::PI * x).sin() * ln_gamma(1.0 - x).exp())
    } else {
        let x = x - 1.0;
        let mut a = 0.9999999999998099;
        let t = x + 7.5;
        for (i, &g) in G.iter().enumerate() {
            a += g / (x + (i as f64) + 1.0);
        }
        0.5 * (2.0 * std::f64::consts::PI).ln() + (x + 0.5) * t.ln() - t + a.ln()
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn samples_every_50th_read_and_skips_n_kmers() {
        // Only the 50th read is processed (2% sampling). Use k=7 default.
        let mut m = KmerContent::new();
        let seq = b"AAAAAAACCCCCCC"; // len 14, no N -> 8 windows
        for _ in 0..49 {
            m.update(b"x", b"GGGGGGG", b"IIIIIII");
        }
        m.update(b"x", seq, &[b'I'; 14]);
        assert_eq!(m.total_kmer_count(), (seq.len() - 7 + 1) as u64);
        assert_eq!(m.kmer_count(b"AAAAAAA"), 1);

        // N-containing kmers are counted in the total but not stored per-kmer.
        let mut m2 = KmerContent::new();
        for _ in 0..50 {
            m2.update(b"x", b"NAAAAAAA", &[b'I'; 8]); // len 8 -> windows "NAAAAAA","AAAAAAA"
        }
        assert_eq!(m2.total_kmer_count(), 2); // both windows counted in totals
        assert_eq!(m2.kmer_count(b"AAAAAAA"), 1); // only the non-N window stored
        assert_eq!(m2.kmer_count(b"NAAAAAA"), 0); // N window not stored
    }

    #[test]
    fn rank_enriched_breaks_obs_exp_ties_by_sequence_ascending() {
        // Two k-mers with identical max obs/exp must come out in a stable,
        // sequence-ascending order regardless of HashMap iteration order, so
        // the reported top-20 is reproducible run-to-run.
        let mk = |seq: [u8; KMER_SIZE]| Enriched {
            seq,
            count: 1,
            lowest_p: 0.0,
            max_obs_exp: 100.0,
            max_pos: 1,
        };
        // Provided in DESCENDING sequence order; expect ASCENDING out.
        let ranked = rank_enriched(vec![mk(*b"CCCCCCC"), mk(*b"AAAAAAA")]);
        let order: Vec<[u8; KMER_SIZE]> = ranked.iter().map(|e| e.seq).collect();
        assert_eq!(order, vec![*b"AAAAAAA", *b"CCCCCCC"]);
    }

    #[test]
    fn binom_cdf_matches_known_values() {
        // Binomial(10, 0.5): P(X<=5) = 0.623046875 exactly.
        assert!((binom_cdf(10, 0.5, 5) - 0.623_046_875).abs() < 1e-9);
        // P(X<=10) = 1.0
        assert!((binom_cdf(10, 0.5, 10) - 1.0).abs() < 1e-12);
    }

    #[test]
    fn binom_sf_matches_known_values_and_is_precise_in_the_tail() {
        // P(X>5) for Binom(10,0.5) = 1 - 0.623046875 = 0.376953125 exactly.
        assert!((binom_sf(10, 0.5, 5) - 0.376_953_125).abs() < 1e-12);
        // P(X>10) = 0.
        assert_eq!(binom_sf(10, 0.5, 10), 0.0);
        // Extreme tail: the direct sum keeps full relative precision where
        // 1 - binom_cdf loses it. P(X>10) for Binom(40, 1e-3) computed exactly.
        let sf = binom_sf(40, 1e-3, 10);
        assert!(
            sf > 1e-30 && sf < 1e-20,
            "tail should be tiny and positive: {sf}"
        );
    }
}
