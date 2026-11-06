//! Execution-wide cache admission. A lease remains charged while any reader
//! owns the entry, even after the entry leaves a residency map.
use std::sync::Arc;
use std::sync::atomic::{AtomicUsize, Ordering};

use datafusion::common::{DataFusionError, Result};

#[derive(Debug)]
pub(crate) struct CacheBudget {
    limit: usize,
    used: AtomicUsize,
    peak: AtomicUsize,
}

impl CacheBudget {
    pub(crate) fn new(limit: usize) -> Arc<Self> {
        Arc::new(Self {
            limit,
            used: AtomicUsize::new(0),
            peak: AtomicUsize::new(0),
        })
    }

    pub(crate) fn from_env(name: &str, default_mb: usize) -> Result<Arc<Self>> {
        let mb = match std::env::var(name) {
            Ok(value) => value.parse::<usize>().map_err(|_| {
                DataFusionError::Execution(format!(
                    "{name} must be a nonnegative integer number of MiB, got {value:?}"
                ))
            })?,
            Err(std::env::VarError::NotPresent) => default_mb,
            Err(error) => return Err(DataFusionError::Execution(format!("{name}: {error}"))),
        };
        let bytes = mb
            .checked_mul(1024 * 1024)
            .ok_or_else(|| DataFusionError::Execution(format!("{name} is too large")))?;
        Ok(Self::new(bytes))
    }

    pub(crate) fn admit(self: &Arc<Self>, bytes: usize) -> Option<BudgetLease> {
        if self.limit == 0 {
            return None;
        }
        let mut used = self.used.load(Ordering::Relaxed);
        loop {
            let next = used.checked_add(bytes)?;
            if next > self.limit {
                return None;
            }
            match self
                .used
                .compare_exchange_weak(used, next, Ordering::Relaxed, Ordering::Relaxed)
            {
                Ok(_) => {
                    self.peak.fetch_max(next, Ordering::Relaxed);
                    return Some(BudgetLease {
                        budget: Arc::clone(self),
                        bytes,
                    });
                }
                Err(observed) => used = observed,
            }
        }
    }

    pub(crate) fn limit(&self) -> usize {
        self.limit
    }
    pub(crate) fn used(&self) -> usize {
        self.used.load(Ordering::Relaxed)
    }
    pub(crate) fn peak(&self) -> usize {
        self.peak.load(Ordering::Relaxed)
    }
}

#[derive(Debug)]
pub(crate) struct BudgetLease {
    budget: Arc<CacheBudget>,
    bytes: usize,
}

impl Drop for BudgetLease {
    fn drop(&mut self) {
        self.budget.used.fetch_sub(self.bytes, Ordering::Relaxed);
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn pinned_payload_stays_charged_until_last_reader_drops() {
        let budget = CacheBudget::new(64);
        let resident = Arc::new(budget.admit(64).unwrap());
        let reader = Arc::clone(&resident);
        drop(resident);
        assert_eq!(budget.used(), 64);
        assert!(budget.admit(1).is_none());
        drop(reader);
        assert_eq!(budget.used(), 0);
        assert!(budget.admit(64).is_some());
        assert!(CacheBudget::new(0).admit(0).is_none());
        assert!(budget.admit(65).is_none());
    }

    #[test]
    fn concurrent_admission_cannot_overcommit() {
        let budget = CacheBudget::new(32);
        std::thread::scope(|scope| {
            let workers: Vec<_> = (0..16).map(|_| scope.spawn(|| budget.admit(8))).collect();
            let leases: Vec<_> = workers
                .into_iter()
                .filter_map(|w| w.join().unwrap())
                .collect();
            assert_eq!(leases.len(), 4);
            assert_eq!(budget.used(), 32);
        });
        assert_eq!(budget.used(), 0);
        assert_eq!(budget.peak(), 32);
    }
}
