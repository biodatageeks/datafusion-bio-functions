//! Bounded decoded Parquet pages. No prediction value is quantized: cached
//! Arrow binary buffers retain the original versioned SIFT/PolyPhen blobs.
use super::SinglePathParquetSiftLookup;
use crate::cache_budget::{BudgetLease, CacheBudget};
use crate::parquet_cache::page_dir::IoCounters;
use datafusion::arrow::array::RecordBatch;
use datafusion::common::{DataFusionError, Result};
use std::collections::HashMap;
use std::sync::atomic::{AtomicU64, Ordering};
use std::sync::{Arc, Mutex};
use tokio::sync::OnceCell;

const MAX_ENTRIES: usize = 4096;
const METADATA_BYTES: usize = 512 * 1024;
const FILL_BYTES: usize = 16 * 1024 * 1024;
type PageRange = (u64, u64);

pub(crate) struct PageCacheBudget {
    pub(crate) payload: Arc<CacheBudget>,
    scratch: Arc<CacheBudget>,
}

impl PageCacheBudget {
    pub(crate) fn from_env() -> Result<Arc<Self>> {
        Ok(Arc::new(Self {
            payload: CacheBudget::from_env("VEP_SIFT_CACHE_MB", 256)?,
            scratch: CacheBudget::new(2 * FILL_BYTES),
        }))
    }

    #[cfg(test)]
    pub(crate) fn new(bytes: usize) -> Arc<Self> {
        Arc::new(Self {
            payload: CacheBudget::new(bytes),
            scratch: CacheBudget::new(2 * FILL_BYTES),
        })
    }
}

pub(super) struct CachedPage {
    pub(super) batch: RecordBatch,
    pub(super) keys_sorted: bool,
    _lease: BudgetLease,
}

type FillResult = std::result::Result<Option<Arc<CachedPage>>, Arc<str>>;
struct Slot {
    value: OnceCell<FillResult>,
    _metadata: BudgetLease,
}

#[derive(Default)]
struct State {
    entries: HashMap<PageRange, (Arc<Slot>, u64)>,
    clock: u64,
}

#[derive(Default)]
struct Counters {
    hits: AtomicU64,
    misses: AtomicU64,
    loads: AtomicU64,
    evictions: AtomicU64,
    fallback: AtomicU64,
    rows: AtomicU64,
    decoded_bytes: AtomicU64,
    physical_bytes: AtomicU64,
    read_ops: AtomicU64,
}

pub(crate) struct DecodedPageCache {
    budget: Arc<PageCacheBudget>,
    state: Mutex<State>,
    counters: Counters,
    // The table's retained capacity remains charged even when entries leave.
    _metadata: Option<BudgetLease>,
}

impl DecodedPageCache {
    pub(crate) fn new(budget: Arc<PageCacheBudget>) -> Self {
        let metadata = budget.payload.admit(METADATA_BYTES);
        let entries = if metadata.is_some() {
            HashMap::with_capacity(MAX_ENTRIES)
        } else {
            HashMap::new()
        };
        Self {
            budget,
            state: Mutex::new(State { entries, clock: 0 }),
            counters: Counters::default(),
            _metadata: metadata,
        }
    }

    fn lock(&self) -> Result<std::sync::MutexGuard<'_, State>> {
        self.state.lock().map_err(|error| {
            DataFusionError::Execution(format!("SIFT page cache lock poisoned: {error}"))
        })
    }

    fn evict_one(&self) -> Result<bool> {
        let removed = {
            let mut state = self.lock()?;
            let oldest = state
                .entries
                .iter()
                .filter(|(_, (slot, _))| slot.value.get().is_some())
                .min_by_key(|(_, (_, age))| *age)
                .map(|(&range, _)| range);
            oldest.and_then(|range| state.entries.remove(&range))
        };
        if removed.is_some() {
            self.counters.evictions.fetch_add(1, Ordering::Relaxed);
        }
        // Dropping the map's Arc does not release a pinned reader's payload.
        Ok(removed.is_some())
    }

    fn admit(&self, bytes: usize) -> Result<Option<BudgetLease>> {
        loop {
            if let Some(lease) = self.budget.payload.admit(bytes) {
                return Ok(Some(lease));
            }
            if !self.evict_one()? {
                return Ok(None);
            }
        }
    }

    pub(super) async fn get(
        &self,
        lookup: &SinglePathParquetSiftLookup,
        range: PageRange,
    ) -> Result<Option<Arc<CachedPage>>> {
        if self._metadata.is_none() {
            return Ok(None);
        }
        let slot = loop {
            let mut state = self.lock()?;
            state.clock = state.clock.wrapping_add(1);
            let clock = state.clock;
            if let Some((slot, age)) = state.entries.get_mut(&range) {
                *age = clock;
                self.counters.hits.fetch_add(1, Ordering::Relaxed);
                break Arc::clone(slot);
            }
            if state.entries.len() >= MAX_ENTRIES {
                drop(state);
                if !self.evict_one()? {
                    return Ok(None);
                }
                continue;
            }
            // Nonblocking metadata admission. Never wait while holding this lock.
            let Some(metadata) = self.budget.payload.admit(512) else {
                drop(state);
                if !self.evict_one()? {
                    return Ok(None);
                }
                continue;
            };
            let slot = Arc::new(Slot {
                value: OnceCell::new(),
                _metadata: metadata,
            });
            state.entries.insert(range, (Arc::clone(&slot), clock));
            self.counters.misses.fetch_add(1, Ordering::Relaxed);
            break slot;
        };
        // One initialization for a granule, with decompression outside the
        // metadata lock. Errors are shared with concurrent waiters, not misses.
        let result = slot
            .value
            .get_or_init(|| async {
                self.fill(lookup, range)
                    .await
                    .map_err(|error| Arc::<str>::from(error.to_string()))
            })
            .await
            .clone();
        if !matches!(result, Ok(Some(_))) {
            let mut state = self.lock()?;
            if state
                .entries
                .get(&range)
                .is_some_and(|(current, _)| Arc::ptr_eq(current, &slot))
            {
                state.entries.remove(&range);
            }
        }
        result.map_err(|error| DataFusionError::Execution(error.to_string()))
    }

    async fn fill(
        &self,
        lookup: &SinglePathParquetSiftLookup,
        range: PageRange,
    ) -> Result<Option<Arc<CachedPage>>> {
        let Some((resident_bound, scratch_bound)) = lookup.decoded_page_bound(range) else {
            self.counters.fallback.fetch_add(1, Ordering::Relaxed);
            return Ok(None);
        };
        if scratch_bound > FILL_BYTES || resident_bound > self.budget.payload.limit() {
            self.counters.fallback.fetch_add(1, Ordering::Relaxed);
            return Ok(None);
        }
        // A fixed 16 MiB scratch permit also caps bulk fills at two across all
        // active/prefetched contigs of the annotation execution.
        let Some(_scratch) = self.budget.scratch.admit(FILL_BYTES) else {
            self.counters.fallback.fetch_add(1, Ordering::Relaxed);
            return Ok(None);
        };
        let Some(mut lease) = self.admit(resident_bound)? else {
            self.counters.fallback.fetch_add(1, Ordering::Relaxed);
            return Ok(None);
        };
        let counters = IoCounters::new();
        let batch = lookup.read_page_range(range, &counters).await?;
        let bytes = batch.get_array_memory_size() + std::mem::size_of::<CachedPage>() + 128;
        if bytes > resident_bound {
            return Err(DataFusionError::Execution(format!(
                "decoded SIFT page exceeded its metadata bound: {bytes} > {resident_bound}"
            )));
        }
        lease.shrink(bytes);
        let (physical_bytes, read_ops) = counters.snapshot();
        self.counters.loads.fetch_add(1, Ordering::Relaxed);
        self.counters
            .rows
            .fetch_add(batch.num_rows() as u64, Ordering::Relaxed);
        self.counters
            .decoded_bytes
            .fetch_add(bytes as u64, Ordering::Relaxed);
        self.counters
            .physical_bytes
            .fetch_add(physical_bytes, Ordering::Relaxed);
        self.counters
            .read_ops
            .fetch_add(read_ops, Ordering::Relaxed);
        let keys_sorted = batch
            .column_by_name("key")
            .and_then(|array| {
                array
                    .as_any()
                    .downcast_ref::<datafusion::arrow::array::UInt64Array>()
            })
            .is_some_and(|keys| keys.values().windows(2).all(|pair| pair[0] <= pair[1]));
        Ok(Some(Arc::new(CachedPage {
            batch,
            keys_sorted,
            _lease: lease,
        })))
    }
}

impl Drop for DecodedPageCache {
    fn drop(&mut self) {
        if std::env::var_os("VEP_ENGINE_PROFILE").is_some_and(|value| value == "1") {
            let c = &self.counters;
            eprintln!(
                "[VEP_SIFT_PAGE_CACHE] hits={} misses={} page_loads={} decoded_rows={} decoded_bytes={} physical_bytes={} read_ops={} evictions={} admission_fallbacks={} execution_payload_peak_bytes={} execution_scratch_peak_bytes={}",
                c.hits.load(Ordering::Relaxed),
                c.misses.load(Ordering::Relaxed),
                c.loads.load(Ordering::Relaxed),
                c.rows.load(Ordering::Relaxed),
                c.decoded_bytes.load(Ordering::Relaxed),
                c.physical_bytes.load(Ordering::Relaxed),
                c.read_ops.load(Ordering::Relaxed),
                c.evictions.load(Ordering::Relaxed),
                c.fallback.load(Ordering::Relaxed),
                self.budget.payload.peak(),
                self.budget.scratch.peak()
            );
        }
    }
}

#[cfg(test)]
mod tests {
    use super::super::write_sift_parquet;
    use super::*;
    use datafusion::arrow::array::{Array, BinaryArray, UInt64Array};
    use datafusion::arrow::datatypes::{DataType, Field, Schema};

    fn fixture(keys: &[u64], payload: &[u8]) -> tempfile::NamedTempFile {
        let schema = Arc::new(Schema::new(vec![
            Field::new("key", DataType::UInt64, false),
            Field::new("sift", DataType::Binary, true),
            Field::new("poly", DataType::Binary, true),
        ]));
        let batch = RecordBatch::try_new(
            schema,
            vec![
                Arc::new(UInt64Array::from(keys.to_vec())),
                Arc::new(BinaryArray::from_vec(vec![payload; keys.len()])),
                Arc::new(BinaryArray::from(vec![None::<&[u8]>; keys.len()])),
            ],
        )
        .unwrap();
        let file = tempfile::NamedTempFile::new().unwrap();
        write_sift_parquet(file.path(), &[batch]).unwrap();
        file
    }

    async fn lookup(
        file: &tempfile::NamedTempFile,
        budget: Arc<PageCacheBudget>,
    ) -> SinglePathParquetSiftLookup {
        SinglePathParquetSiftLookup::open(
            file.path(),
            vec!["key".into(), "sift".into(), "poly".into()],
        )
        .await
        .unwrap()
        .with_page_cache(budget)
    }

    fn assert_same(a: &RecordBatch, b: &RecordBatch) {
        assert_eq!(a.schema(), b.schema());
        assert_eq!(a.num_rows(), b.num_rows());
        for (a, b) in a.columns().iter().zip(b.columns()) {
            assert_eq!(a.to_data(), b.to_data());
        }
    }

    #[tokio::test]
    async fn cached_selection_preserves_unsorted_page_order_and_duplicate_keys() {
        let file = fixture(&[3, 1, 2, 2, 0], b"payload");
        let lookup = lookup(&file, PageCacheBudget::new(4 * 1024 * 1024)).await;
        let (actual, present) = lookup.take_keys(&[2, 3, 2]).await.unwrap();
        let (expected, old_present) = lookup.take_keys_uncached(&[2, 3, 2]).await.unwrap();
        assert_same(&actual, &expected);
        assert_eq!(present, old_present);
        assert_eq!(actual.num_rows(), 3);
        assert_eq!(
            lookup
                .decoded_cache
                .as_ref()
                .unwrap()
                .counters
                .loads
                .load(Ordering::Relaxed),
            1
        );
    }

    #[tokio::test]
    async fn new_positions_on_resident_page_require_no_new_decode() {
        let file = fixture(
            &(0..1000)
                .map(|i| (7u64 << 32) | (i * 2))
                .collect::<Vec<_>>(),
            b"unchanged",
        );
        let lookup = lookup(&file, PageCacheBudget::new(4 * 1024 * 1024)).await;
        let first = (7u64 << 32) | 100;
        let last = (7u64 << 32) | 200;
        let (a, _) = lookup.take_keys(&[first]).await.unwrap();
        let cache = lookup.decoded_cache.as_ref().unwrap();
        let loads = cache.counters.loads.load(Ordering::Relaxed);
        assert!(loads > 0);
        let (b, _) = lookup.take_keys(&[last, last + 1]).await.unwrap();
        assert_eq!(cache.counters.loads.load(Ordering::Relaxed), loads);
        assert_eq!(a.num_rows(), 1);
        assert_eq!(b.num_rows(), 1); // the odd key is a genuine negative result
        let (old, _) = lookup.take_keys_uncached(&[last, last + 1]).await.unwrap();
        assert_same(&b, &old);
    }

    #[tokio::test]
    async fn concurrent_readers_share_one_fill_and_preserve_pinned_accounting() {
        let file = fixture(&(0..100).collect::<Vec<_>>(), b"payload");
        let budget = PageCacheBudget::new(4 * 1024 * 1024);
        let lookup = lookup(&file, Arc::clone(&budget)).await;
        let results = futures::future::join_all((0..16).map(|i| {
            let lookup = &lookup;
            async move { lookup.take_keys(&[i]).await.unwrap() }
        }))
        .await;
        assert!(results.iter().all(|(b, _)| b.num_rows() == 1));
        let cache = lookup.decoded_cache.as_ref().unwrap();
        assert_eq!(cache.counters.loads.load(Ordering::Relaxed), 1);
        let range = lookup.page_dir.resolve_page_ranges(&[1])[0];
        let pinned = cache.get(&lookup, range).await.unwrap().unwrap();
        assert!(cache.evict_one().unwrap());
        assert!(budget.payload.used() > METADATA_BYTES);
        drop(pinned);
        assert_eq!(budget.payload.used(), METADATA_BYTES);
        lookup.take_keys(&[2]).await.unwrap();
        assert_eq!(cache.counters.loads.load(Ordering::Relaxed), 2);
        assert!(budget.payload.peak() <= budget.payload.limit());
        drop(lookup);
        assert_eq!(budget.payload.used(), 0);
    }

    #[tokio::test]
    async fn zero_tiny_and_busy_budgets_fall_back_without_negative_caching() {
        let file = fixture(&[0, 2, u64::MAX], b"payload");
        for bytes in [0, 1, METADATA_BYTES] {
            let lookup = lookup(&file, PageCacheBudget::new(bytes)).await;
            let (actual, present) = lookup.take_keys(&[0, 1, u64::MAX]).await.unwrap();
            let (expected, _) = lookup.take_keys_uncached(&[0, 1, u64::MAX]).await.unwrap();
            assert_same(&actual, &expected);
            assert_eq!(present, vec![0, u64::MAX]);
            assert_eq!(
                lookup
                    .decoded_cache
                    .as_ref()
                    .unwrap()
                    .counters
                    .loads
                    .load(Ordering::Relaxed),
                0
            );
        }
        let budget = PageCacheBudget::new(4 * 1024 * 1024);
        let lookup = lookup(&file, Arc::clone(&budget)).await;
        let a = budget.scratch.admit(FILL_BYTES).unwrap();
        let b = budget.scratch.admit(FILL_BYTES).unwrap();
        assert!(budget.scratch.admit(1).is_none());
        lookup.take_keys(&[2]).await.unwrap();
        assert_eq!(
            lookup
                .decoded_cache
                .as_ref()
                .unwrap()
                .counters
                .loads
                .load(Ordering::Relaxed),
            0
        );
        drop((a, b));
        lookup.take_keys(&[2]).await.unwrap();
        assert_eq!(
            lookup
                .decoded_cache
                .as_ref()
                .unwrap()
                .counters
                .loads
                .load(Ordering::Relaxed),
            1
        );
    }

    #[tokio::test]
    async fn shard_identity_keeps_equal_uids_separate_and_v1_bits_unchanged() {
        use crate::cache_common::{
            SiftBlobVersion, deserialize_position_predictions_versioned, serialize_position_entries,
        };
        use crate::transcript_consequence::CompactPrediction;
        let payload = serialize_position_entries(&[CompactPrediction {
            position: 7,
            amino_acid: CompactPrediction::encode_amino_acid("Z").unwrap(),
            prediction: 8,
            score: 0.12345679,
        }]);
        let key = ((u32::MAX as u64) << 32) | 7;
        let first = fixture(&[key], &payload);
        let second = fixture(&[key], b"different");
        let budget = PageCacheBudget::new(4 * 1024 * 1024);
        let a = lookup(&first, Arc::clone(&budget)).await;
        let b = lookup(&second, budget).await;
        let (actual, _) = a.take_keys(&[key]).await.unwrap();
        let bytes = actual
            .column_by_name("sift")
            .unwrap()
            .as_any()
            .downcast_ref::<BinaryArray>()
            .unwrap()
            .value(0);
        assert_eq!(bytes, payload);
        let predictions = deserialize_position_predictions_versioned(
            7,
            bytes,
            &[],
            SiftBlobVersion::V1Interleaved,
        )
        .unwrap();
        assert_eq!(
            predictions.sift[0].amino_acid,
            CompactPrediction::encode_amino_acid("Z").unwrap()
        );
        assert_eq!(predictions.sift[0].prediction, 8);
        assert_eq!(predictions.sift[0].score.to_bits(), 0.12345679f32.to_bits());
        let (other, _) = b.take_keys(&[key]).await.unwrap();
        assert_ne!(actual.column(1).to_data(), other.column(1).to_data());
    }

    #[tokio::test]
    async fn both_v2_predictor_grids_remain_bit_exact_in_cached_pages() {
        use crate::cache_common::{
            POLY_CODEC, SIFT_CODEC, deserialize_position_entries_v2, serialize_position_entries_v2,
        };
        use crate::transcript_consequence::CompactPrediction;
        for (codec, score) in [
            (SIFT_CODEC, 4.0f32 / 100.0),
            (POLY_CODEC, 791.0f32 / 1000.0),
        ] {
            let entries: Vec<_> = (0..26)
                .map(|aa| CompactPrediction {
                    position: 9,
                    amino_acid: aa,
                    prediction: aa % 9,
                    score,
                })
                .collect();
            let payload = serialize_position_entries_v2(&entries, codec).unwrap();
            let file = fixture(&[9], &payload);
            let lookup = lookup(&file, PageCacheBudget::new(4 * 1024 * 1024)).await;
            let (actual, _) = lookup.take_keys(&[9]).await.unwrap();
            let bytes = actual
                .column_by_name("sift")
                .unwrap()
                .as_any()
                .downcast_ref::<BinaryArray>()
                .unwrap()
                .value(0);
            assert_eq!(bytes, payload);
            let decoded = deserialize_position_entries_v2(9, bytes, codec).unwrap();
            assert_eq!(decoded.len(), entries.len());
            for (before, after) in entries.iter().zip(decoded) {
                assert_eq!(before.amino_acid, after.amino_acid);
                assert_eq!(before.prediction, after.prediction);
                assert_eq!(before.score.to_bits(), after.score.to_bits());
            }
        }
    }

    #[tokio::test]
    async fn bulk_fill_errors_propagate_and_do_not_become_cached_absences() {
        let file = fixture(&[0, 1, 2], b"payload");
        let original = std::fs::read(file.path()).unwrap();
        let lookup = lookup(&file, PageCacheBudget::new(4 * 1024 * 1024)).await;
        std::fs::write(file.path(), b"broken").unwrap();
        let results = futures::future::join_all((0..8).map(|_| lookup.take_keys(&[1]))).await;
        assert!(results.iter().all(Result::is_err));
        std::fs::write(file.path(), original).unwrap();
        assert_eq!(lookup.take_keys(&[1]).await.unwrap().0.num_rows(), 1);
    }

    #[tokio::test]
    async fn oversized_page_uses_point_fallback_without_bulk_materialization() {
        let payload = vec![b'x'; 5 * 1024 * 1024];
        let file = fixture(&[1], &payload);
        let lookup = lookup(&file, PageCacheBudget::new(64 * 1024 * 1024)).await;
        let (actual, _) = lookup.take_keys(&[1]).await.unwrap();
        assert_eq!(
            actual
                .column_by_name("sift")
                .unwrap()
                .as_any()
                .downcast_ref::<BinaryArray>()
                .unwrap()
                .value(0),
            payload
        );
        let cache = lookup.decoded_cache.as_ref().unwrap();
        assert_eq!(cache.counters.loads.load(Ordering::Relaxed), 0);
        assert!(cache.counters.fallback.load(Ordering::Relaxed) > 0);
    }
}
