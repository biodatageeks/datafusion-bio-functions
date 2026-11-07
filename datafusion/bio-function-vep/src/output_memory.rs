//! Opt-in formatter ownership estimates, separate from performance trace gates.
//! Counts simultaneously active formatting work, not allocator RSS or engine
//! builders. CSQ accounting is an upper bound if callers retain Arrow aliases.
use std::sync::{Mutex, OnceLock};

#[derive(Clone, Copy, Debug)]
pub(crate) enum Kind {
    Csq = 0,
    Lines = 1,
    Body = 2,
}

#[derive(Default)]
struct Totals {
    current: [usize; 3],
    peak: [usize; 3],
    combined_peak: usize,
    batches: usize,
    rows: usize,
}

static ENABLED: OnceLock<bool> = OnceLock::new();
static TOTALS: Mutex<Totals> = Mutex::new(Totals {
    current: [0; 3],
    peak: [0; 3],
    combined_peak: 0,
    batches: 0,
    rows: 0,
});

pub(crate) fn enabled() -> bool {
    *ENABLED.get_or_init(|| {
        crate::pipeline_trace::enabled_from_env_value(
            std::env::var("VEP_OUTPUT_MEMORY").ok().as_deref(),
        )
    })
}

/// Tracks the lifetime of an owned formatting allocation. No allocation,
/// locking, or traversal of batches occurs when diagnostics are disabled.
#[derive(Debug)]
pub(crate) struct Hold {
    kind: Kind,
    bytes: usize,
    rows: usize,
}

impl Hold {
    pub(crate) fn new(kind: Kind, bytes: usize, rows: usize) -> Option<Self> {
        if !enabled() {
            return None;
        }
        let mut totals = TOTALS.lock().unwrap_or_else(|e| e.into_inner());
        let index = kind as usize;
        totals.current[index] += bytes;
        totals.peak[index] = totals.peak[index].max(totals.current[index]);
        if matches!(kind, Kind::Csq) {
            totals.batches += 1;
            totals.rows += rows;
        }
        let combined: usize = totals.current.iter().sum();
        if combined > totals.combined_peak {
            totals.combined_peak = combined;
            eprintln!(
                "[VEP_OUTPUT_MEMORY] scope=formatter csq_upper_bytes={} lines_capacity_bytes={} body_capacity_bytes={} active_batches={} active_rows={} concurrent_upper_peak_bytes={} csq_upper_peak_bytes={} lines_peak_bytes={} body_peak_bytes={}",
                totals.current[0],
                totals.current[1],
                totals.current[2],
                totals.batches,
                totals.rows,
                totals.combined_peak,
                totals.peak[0],
                totals.peak[1],
                totals.peak[2]
            );
        }
        Some(Self { kind, bytes, rows })
    }
}

impl Drop for Hold {
    fn drop(&mut self) {
        let mut totals = TOTALS.lock().unwrap_or_else(|e| e.into_inner());
        totals.current[self.kind as usize] -= self.bytes;
        if matches!(self.kind, Kind::Csq) {
            totals.batches -= 1;
            totals.rows -= self.rows;
        }
    }
}
