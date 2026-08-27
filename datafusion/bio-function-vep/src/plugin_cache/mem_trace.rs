//! Per-reservation memory tracing for the plugin-cache build.
//!
//! DataFusion's pool-exhaustion error names the *consumer* ("HashJoinInput"),
//! and the build has two join stages -- the semi-join inside
//! `plugin_variation_probe` and the tier LEFT JOIN itself. Both register under
//! that same name, so stage-scoped attribution is needed before retrying.
//!
//! This wraps any [`MemoryPool`] and tracks each *reservation* separately
//! (keyed by consumer identity, not name), reporting each one's peak. Opt-in
//! via `VEPYR_TRACE_MEMORY=1`. The wrapper is always present in the adaptive
//! join context so it can attribute rejected reservations, but detailed
//! tracking and logging are disabled when the variable is unset.

use std::collections::{HashMap, HashSet};
use std::sync::Mutex;
use std::sync::atomic::{AtomicUsize, Ordering};

use datafusion::common::Result;
use datafusion::execution::memory_pool::{
    MemoryConsumer, MemoryLimit, MemoryPool, MemoryReservation,
};
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
    failed_consumers: Mutex<HashSet<String>>,
    next_id: AtomicUsize,
    trace: bool,
}

const STEP: usize = 256 * 1024 * 1024;

impl TracingPool {
    pub fn new(inner: std::sync::Arc<dyn MemoryPool>) -> Self {
        Self {
            inner,
            seen: Mutex::new(HashMap::new()),
            failed_consumers: Mutex::new(HashSet::new()),
            next_id: AtomicUsize::new(0),
            trace: tracing_enabled(),
        }
    }

    /// Remove any failure attributed to the preceding stage/query.
    pub fn clear_failures(&self) {
        self.failed_consumers
            .lock()
            .expect("memory failure trace poisoned")
            .clear();
    }

    fn consumer_failed(&self, prefix: &str) -> bool {
        self.failed_consumers
            .lock()
            .expect("memory failure trace poisoned")
            .iter()
            .any(|name| name.starts_with(prefix))
    }

    /// True when a reservation rejected during the current stage belonged to
    /// a hash join build. The error itself is also checked at the retry site,
    /// so an unrelated concurrent rejection cannot cause a false retry.
    pub fn hash_join_failed(&self) -> bool {
        self.consumer_failed("HashJoinInput")
    }

    /// True when the final query's external sorter was the consumer whose
    /// reservation failed while an unspillable hash build remained live.
    pub fn external_sorter_failed(&self) -> bool {
        self.consumer_failed("ExternalSorter")
    }

    /// Identity of the reservation's consumer. Two consumers sharing a name
    /// are distinct objects, and the reservation owns its consumer, so the
    /// address is stable for as long as we need to attribute growth to it.
    fn key(reservation: &MemoryReservation) -> usize {
        std::ptr::from_ref(reservation.consumer()) as usize
    }

    fn record(&self, reservation: &MemoryReservation, additional: usize) {
        if !self.trace {
            return;
        }
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
        if !self.trace {
            return;
        }
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
        let result = self.inner.try_grow(reservation, additional);
        if result.is_err() {
            self.failed_consumers
                .lock()
                .expect("memory failure trace poisoned")
                .insert(reservation.consumer().name().to_string());
        }
        result
    }

    fn reserved(&self) -> usize {
        self.inner.reserved()
    }

    fn memory_limit(&self) -> MemoryLimit {
        self.inner.memory_limit()
    }
}

/// Report how a batch vector will be *charged* versus what it actually holds.
///
/// `HashJoinExec` reserves per batch with `get_record_batch_memory_size`, which
/// de-duplicates buffers only within one batch. So a vector whose batches share
/// buffers is charged once per sharing batch. This logs both numbers plus the
/// max batch size, which separates the two ways that happens: oversized parents
/// (DataFusion's `BatchSplitStream` re-slices them into zero-copy children) and
/// batches that already share on arrival.
///
/// Gated on [`tracing_enabled`]; costs nothing unset.
pub fn describe_batches(label: &str, batches: &[datafusion::arrow::record_batch::RecordBatch]) {
    use datafusion::arrow::array::Array;
    use std::collections::HashSet;

    if !tracing_enabled() {
        return;
    }
    let mut charged = 0usize;
    let mut global: HashSet<usize> = HashSet::new();
    let mut distinct = 0usize;
    let mut max_rows = 0usize;
    let rows: usize = batches.iter().map(|b| b.num_rows()).sum();

    for batch in batches {
        max_rows = max_rows.max(batch.num_rows());
        let mut local: HashSet<usize> = HashSet::new();
        for column in batch.columns() {
            for buffer in column.to_data().buffers() {
                let key = buffer.as_ptr() as usize;
                if local.insert(key) {
                    charged += buffer.capacity();
                }
                if global.insert(key) {
                    distinct += buffer.capacity();
                }
            }
        }
    }
    let gb = |n: usize| n as f64 / 1024.0 / 1024.0 / 1024.0;
    info!(
        "mem_trace: {label}: {} batches, {rows} rows, max {max_rows} rows/batch -- \
         charged {:.2} GB vs {:.2} GB distinct ({:.1}x)",
        batches.len(),
        gb(charged),
        gb(distinct),
        charged as f64 / distinct.max(1) as f64,
    );
}
