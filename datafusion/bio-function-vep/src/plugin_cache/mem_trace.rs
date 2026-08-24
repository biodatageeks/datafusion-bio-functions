//! Per-reservation memory tracing for the plugin-cache build.
//!
//! DataFusion's pool-exhaustion error names the *consumer* ("HashJoinInput"),
//! and the tier-join plan contains two `HashJoinExec`s -- the semi-join inside
//! `plugin_variation_probe` and the tier LEFT JOIN itself. Both register under
//! that same name, so the error alone cannot say which one consumed the pool,
//! and the answer decides whether a bigger budget is a fix or a distraction.
//!
//! This wraps any [`MemoryPool`] and tracks each *reservation* separately
//! (keyed by consumer identity, not name), reporting each one's peak. Opt-in
//! via `VEPYR_TRACE_MEMORY=1`; unset it costs nothing because the pool is only
//! wrapped when the variable is set.

use std::collections::HashMap;
use std::sync::Mutex;
use std::sync::atomic::{AtomicUsize, Ordering};

use datafusion::common::Result;
use datafusion::execution::memory_pool::{MemoryConsumer, MemoryPool, MemoryReservation};
use log::info;

pub const TRACE_ENV: &str = "VEPYR_TRACE_MEMORY";

pub fn tracing_enabled() -> bool {
    std::env::var(TRACE_ENV).is_ok_and(|v| !v.is_empty() && v != "0")
}

#[derive(Debug, Default)]
struct Tracked {
    /// Stable label for this reservation, e.g. "HashJoinInput#2".
    label: String,
    peak: usize,
    /// Highest 256 MiB step already logged, so a growing reservation reports
    /// progress without one line per batch.
    logged_steps: usize,
}

/// Wraps a pool and reports the peak of every individual reservation.
#[derive(Debug)]
pub struct TracingPool {
    inner: std::sync::Arc<dyn MemoryPool>,
    seen: Mutex<HashMap<usize, Tracked>>,
    next_id: AtomicUsize,
}

const STEP: usize = 256 * 1024 * 1024;

impl TracingPool {
    pub fn new(inner: std::sync::Arc<dyn MemoryPool>) -> Self {
        Self {
            inner,
            seen: Mutex::new(HashMap::new()),
            next_id: AtomicUsize::new(0),
        }
    }

    /// Identity of the reservation's consumer. Two consumers sharing a name
    /// are distinct objects, and the reservation owns its consumer, so the
    /// address is stable for as long as we need to attribute growth to it.
    fn key(reservation: &MemoryReservation) -> usize {
        std::ptr::from_ref(reservation.consumer()) as usize
    }

    fn record(&self, reservation: &MemoryReservation, additional: usize) {
        let key = Self::key(reservation);
        let size = reservation.size() + additional;
        let mut seen = self.seen.lock().expect("memory trace poisoned");
        let id = self.next_id.load(Ordering::Relaxed);
        let entry = seen.entry(key).or_insert_with(|| {
            self.next_id.fetch_add(1, Ordering::Relaxed);
            Tracked {
                label: format!("{}#{id}", reservation.consumer().name()),
                ..Tracked::default()
            }
        });
        entry.peak = entry.peak.max(size);
        let step = size / STEP;
        if step > entry.logged_steps {
            entry.logged_steps = step;
            info!(
                "mem_trace: {} grew to {:.2} GB (pool reserved {:.2} GB)",
                entry.label,
                size as f64 / 1024.0 / 1024.0 / 1024.0,
                self.inner.reserved() as f64 / 1024.0 / 1024.0 / 1024.0,
            );
        }
    }

    /// Log every reservation's peak, largest first. Call once the build is
    /// done (or has failed) -- this is the line that says which operator
    /// actually consumed the pool.
    pub fn report(&self) {
        let seen = self.seen.lock().expect("memory trace poisoned");
        let mut rows: Vec<_> = seen.values().map(|t| (t.label.clone(), t.peak)).collect();
        rows.sort_by_key(|(_, peak)| std::cmp::Reverse(*peak));
        info!("mem_trace: peak per reservation ({} total):", rows.len());
        for (label, peak) in rows {
            info!(
                "mem_trace:   {label:<24} {:.2} GB",
                peak as f64 / 1024.0 / 1024.0 / 1024.0
            );
        }
    }
}

impl MemoryPool for TracingPool {
    fn register(&self, consumer: &MemoryConsumer) {
        self.inner.register(consumer);
    }

    fn unregister(&self, consumer: &MemoryConsumer) {
        self.inner.unregister(consumer);
    }

    fn grow(&self, reservation: &MemoryReservation, additional: usize) {
        self.record(reservation, additional);
        self.inner.grow(reservation, additional);
    }

    fn shrink(&self, reservation: &MemoryReservation, shrink: usize) {
        self.inner.shrink(reservation, shrink);
    }

    fn try_grow(&self, reservation: &MemoryReservation, additional: usize) -> Result<()> {
        self.record(reservation, additional);
        self.inner.try_grow(reservation, additional)
    }

    fn reserved(&self) -> usize {
        self.inner.reserved()
    }
}
