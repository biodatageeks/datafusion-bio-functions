# Fjall 3.1.0 Optimization Proposals

Findings from reviewing fjall 3.0.3–3.1.0 changes and analyzing sorted WGS VCF access patterns
against a 1.2B-entry variation cache. These proposals complement the existing decisions in `design.md`
without changing the spec requirements.

**fjall version upgrade path:** 3.0.2 → 3.1.0 (latest as of 2026-03-14)

**Changelog 3.0.2 → 3.1.0:**
- **3.0.3** (2025-02-27): Correctness fix for database recovery after 3.0.1; byteview 0.10.1
- **3.0.4** (2025-03-05): Simplified sealed memtable flush after recovery (#257); FIFO docs
- **3.1.0** (2025-03-07): **Compaction filters** (`fjall::compaction::filter` module); MSRV reduced to 1.90

---

## Optimization 1: Update fjall to 3.1.0

**Change:** `fjall = { version = "3.1.0", optional = true }` in `bio-function-vep/Cargo.toml`

**Why:** Picks up the 3.0.3 recovery correctness fix (critical bugfix) and enables compaction filters (Optimization 2). No breaking API changes — `KeyspaceCreateOptions`, `Keyspace`, `DatabaseBuilder`, and `Leveled` APIs are identical to 3.0.2.

**Impact:** Trivial change, very low risk. Prerequisite for Optimization 2.

---

## Optimization 2: Compaction Filters for Post-Ingest Optimization

**Problem addressed:** The design identifies two pain points around post-ingest compaction:
- Decision 5: `start_ingestion()` creates a new run requiring follow-up compaction
- Risk table: `major_compact()` hardcodes a 64 MB target that can silently undo the configured 256 MiB `table_target_size`

**Key finding:** `Keyspace::major_compact()` is **not listed in the public API** for fjall 3.1.0. The design's concern about its 64 MB hardcoded target may be moot — but this also means there is no public "compact everything" escape hatch.

**New fjall 3.1.0 API:**

```rust
// Compaction filter module: fjall::compaction::filter
trait CompactionFilter: Send {
    fn filter_item(&mut self, item: ItemAccessor, ctx: Context) -> Result<Verdict>;
    fn finish(&mut self) {}  // cleanup hook
}

trait Factory: Send + Sync + RefUnwindSafe {
    fn name(&self) -> &str;
    fn make_filter(&self, ctx: &Context) -> Box<dyn CompactionFilter>;
}

enum Verdict {
    Keep,               // retain unchanged
    Remove,             // delete + tombstone
    RemoveWeak,         // delete without tombstone (may resurrect old versions)
    ReplaceValue(Slice), // rewrite value in-place
    Destroy,            // delete without tombstone (only for never-updated keys)
}

struct Context {
    is_last_level: bool,  // whether compacting into final level
}

// Registration: on DatabaseBuilder, not KeyspaceCreateOptions
Database::builder(path)
    .with_compaction_filter_factories(Arc::new(|keyspace_name| {
        // Return filter factory per keyspace, or None
    }))
```

**Proposed usage for VEP cache:**

| Use Case | Verdict | When | Impact |
|---|---|---|---|
| Duplicate position dedup | `Destroy` older duplicate | Incremental re-ingests produce overlapping runs | Eliminates stale entries without full rebuild |
| Tombstone cleanup | `Destroy` | `ctx.is_last_level == true` | Reclaims space at final level |
| Value recompression | `ReplaceValue(new_compressed)` | Dictionary evolution | Allows dictionary upgrades without full rebuild |
| Lazy format migration | `ReplaceValue(migrated)` | FORMAT_V0 → future versions | Zero-downtime format upgrades during background compaction |

**Critical benefit:** Compaction filters run during **normal leveled compaction** which respects `with_table_target_size(256 MiB)`. This eliminates the need for `major_compact()` entirely — the post-ingest path becomes:

```
start_ingestion() → finish() → natural leveled compaction (with filters) → steady state
```

No 64 MB target-size risk. No manual compaction step. The ingested run naturally merges into the LSM tree via the configured leveled strategy.

**Impact assessment:**
- Eliminates post-ingest `major_compact()` step (~5-15 min for 74-80 GB DB)
- Enables incremental cache updates (re-ingest updated chromosomes) without full rebuild
- Zero overhead when filter returns `Verdict::Keep` (the common case during reads)
- **Estimated build time savings:** ~5-15 min (replaces `major_compact()` with natural compaction)

---

## Optimization 3: Range-Scan Merge-Join for Sorted VCF

**Problem addressed:** The current design treats each VCF variant as an independent `get()` call. For 6M lookups against 1.2B entries, each `get()` independently traverses bloom filters (6 levels × hash computation), does binary search in pinned index blocks, and reads/probes data blocks. This ignores the fact that both the VCF input and fjall's key ordering are sorted by `(chrom, position)`.

**Proposed approach:** Replace per-variant `get()` with a merge-join between the sorted VCF stream and a fjall `range()` iterator:

```rust
// Current: 6M independent get() calls
for variant in vcf_batch {
    store.get(encode_key(variant.chrom, variant.start))?;
}

// Proposed: single range iterator, merge-join with VCF cursor
let first_key = encode_key(chrom, batch_first_start);
let last_key = encode_key(chrom, batch_last_start);
let mut kv_iter = keyspace.range(first_key..=last_key);
let mut vcf_idx = 0;

while let Some(kv_entry) = kv_iter.next() {
    let kv_pos = decode_position(&kv_entry.key);

    // Advance VCF cursor past positions not in cache (novel variants)
    while vcf_idx < batch.len() && batch[vcf_idx].start < kv_pos {
        emit_null_row(vcf_idx);  // miss detected structurally
        vcf_idx += 1;
    }

    if vcf_idx < batch.len() && batch[vcf_idx].start == kv_pos {
        let entry = decompress(&kv_entry.value);
        match_alleles(&batch[vcf_idx], &entry);
        vcf_idx += 1;
    }
}
// Remaining VCF variants past iterator → novel variants
```

**Operation-by-operation comparison (6M lookups):**

| Operation | Per-variant `get()` | Range iterator merge-join |
|---|---|---|
| Bloom filter checks | 6M × ~6 levels = 36M checks | **0** — iterator walks blocks directly |
| Index block lookups | 6M binary searches | ~4.5M block transitions only |
| Data block reads | 6M individual lookups (with cache) | Sequential block walk, each block read once |
| Within-block search | 6M hash probes | Iterator cursor advance — O(1) per entry |

**Interaction with extended coordinate probes:** With a range iterator already positioned near the target, checking ±3 nearby positions for indel normalization is a cursor seek, not 3 separate `get()` calls with independent bloom/index traversals.

**Hybrid strategy for sparse regions:**

```rust
let density = variants_in_region / estimated_blocks_in_region;
if density > 0.5 {
    range_scan_merge_join(...)   // dense: WGS, most chromosomes
} else {
    individual_point_lookups(...) // sparse: centromeres, panels with gaps
}
```

**Impact assessment:**

| Component | Savings | Mechanism |
|---|---|---|
| Bloom filter bypass | ~0.5-1s | 36M hash computations eliminated |
| Reduced index traversals | ~0.5-1s | Block-boundary transitions only |
| Sequential data block I/O | ~1-2s | OS prefetch optimal for sequential walk |
| Extended probe optimization | ~0.3-0.5s | Nearby probes = cursor seeks |
| **Total** | **~2-4s** | On 6M-variant WGS |

- **Complexity:** Medium-high — requires rethinking `KvLookupExec` from per-row `get()` to streaming merge-join
- **Risk:** If fjall's `range()` iterator has per-entry overhead beyond raw block walking (e.g., per-entry merge across LSM levels), savings may be lower. Needs benchmarking.
- **Composability:** Excellent — sequential I/O pattern is exactly what Decisions 10-11 (prefetch) and OS readahead optimize for

---

## Optimization 4: Reduced `max_cached_files` for Sequential Access

**Problem addressed:** The design proposes `max_cached_files(512)` as headroom above ~289-300 L6 SSTs. But sorted VCF access is strictly sequential through chromosome-ordered SSTs — earlier chromosomes' SSTs are never revisited.

**Proposed change:**

```rust
Database::builder(path)
    .max_cached_files(Some(128))  // 2 chroms × ~13 SSTs × ~3 levels + headroom
```

**Rationale:** With 256 MiB SSTs and ~80 GB DB:
- ~310 total SSTs at L6
- Each chromosome spans ~14 SSTs (80 GB / 22 chromosomes / 256 MiB)
- Active working set: current + next chromosome (prefetch) ≈ ~30 SSTs
- 128 provides 4x headroom above working set

**Impact assessment:**
- ~30-50 MB less kernel memory for file descriptor metadata
- Slightly faster LRU management (smaller pool)
- **Risk:** If user does random-access queries across chromosomes, this causes file re-opens. Safe for the primary sequential WGS annotation workload.

---

## Optimization 5: `worker_threads(1)` for Read-Only Annotation

**Problem addressed:** fjall defaults to `min(4, num_cores)` background worker threads for compaction and flush. During read-only annotation (the primary workload), none of these fire — the threads are parked idle.

**What the threads do:**

| Worker task | Triggered by | During annotation? |
|---|---|---|
| Memtable flush | `insert()` filling memtable | No — no writes |
| Compaction | Flush producing new L0 SSTs | No — no flushes |
| Blob GC | Compaction producing stale blobs | No — no compaction |

**Cost of idle threads:** Each thread allocates ~8 MB stack (macOS/Linux default), consumes a kernel thread descriptor and Mach port (macOS), and causes minor scheduler noise.

**Proposed change:**

```rust
// Read-only open (VepKvStore::open / open_with_cache_size):
Database::builder(path)
    .worker_threads(1)  // minimum allowed by fjall
    .cache_size(512 * 1024 * 1024)
    .open()

// Write open (VepKvStore::create — unchanged):
Database::builder(path)
    .worker_threads(4)  // need compaction parallelism
    .manual_journal_persist(true)
    .open()
```

This separation already exists in the code — `open()` is for reads, `create()` is for writes. The change is adding `.worker_threads(1)` to `open()`.

**Impact assessment:**

| Metric | Default (4 threads) | Proposed (1 thread) | Savings |
|---|---|---|---|
| Stack memory | ~32 MB | ~8 MB | ~24 MB |
| Kernel threads | 4 | 1 | 3 fewer |
| Context switches | ~100-200/s (idle wakeups) | ~25-50/s | 75% fewer |

- **Risk:** Very low — worker thread is needed only for startup maintenance (stale journal cleanup). After that it sits idle for the entire read-only annotation run.
- **When NOT to do this:** If write-back caching, background maintenance tasks, or compaction-during-read are ever added to the read path.

---

## Optimization 6: `size_of()` for Adaptive Decompression Buffer Sizing

**Problem addressed:** Current `decompress_into_buffer_with_retry` (`kv_store.rs:66-105`) uses a 16x heuristic for initial buffer size and retries with exponential growth when the destination is too small. This wastes allocations.

**fjall API:** `Keyspace::size_of(key) -> Result<Option<u32>>` returns the compressed value size without reading the full value (reads only the index, which is already pinned).

**Proposed change:**

```rust
// Before get(), pre-size the buffer using known compressed size
if let Some(compressed_size) = self.data.size_of(&key_buf).map_err(fjall_err)? {
    let estimated = (compressed_size as usize).saturating_mul(compression_ratio_hint);
    if buf.capacity() < estimated {
        buf.reserve(estimated - buf.capacity());
    }
}
let raw = self.data.get(&key_buf).map_err(fjall_err)?;
```

**Impact assessment:**
- Eliminates retry loops in decompression for the common case
- `size_of()` reads only pinned index blocks — effectively free
- ~0.1-0.3s saved over 6M lookups (fewer allocations, no retries)
- Also useful for KV separation decision: sample `size_of()` to build value-size histogram without reading full values

---

## Optimization 7: `journal_compression(Lz4)` for Fallback Ingest

**Problem addressed:** For ingest code paths that still use the journal (e.g., `batch_insert_raw()` fallback, metadata writes), journal I/O is uncompressed by default.

**Proposed change:**

```rust
// Write open (VepKvStore::create):
Database::builder(path)
    .journal_compression(CompressionType::Lz4)
    .manual_journal_persist(true)
    .open()
```

**Impact assessment:**
- ~30-50% journal I/O reduction during batch insert fallback path
- Only relevant when NOT using `start_ingestion()` (which bypasses journals entirely)
- Zero impact on read path — journals are not read during annotation
- LZ4 compression is ~3 GB/s, negligible CPU overhead

---

## Performance Model: 6M Lookups Against 1.2B Cache

Scaled from the design's 5M/1.15B model to the target workload.

### DB Geometry at 1.2B Entries

| Parameter | Value |
|---|---|
| Cache entries | 1.2B positions |
| Key size | 10B |
| Total data blocks (8 KiB) | ~9.5M |
| On-disk DB size | ~80 GB |
| Bloom filters (L0-L5, pinned) | ~210 MB |
| Index blocks (all levels, pinned) | ~170 MB |
| Pinned metadata total | ~380 MB |
| Cold-start pinned load (NVMe 3 GB/s) | ~127 ms |

### 6M Variant WGS Lookup Profile

| Parameter | Value |
|---|---|
| Total lookups | 6M |
| Hits (95%) | 5.7M |
| Misses (5%) | 300K |
| Distinct data blocks touched | ~4.5M (~47% of 9.5M) |
| Data I/O volume | ~36 GB |

### End-to-End Estimates

| Configuration | Time | vs. Default |
|---|---|---|
| Default fjall (current code, 3.0.2) | **~40-50s** | baseline |
| Tuned (pinned + hash + 8K blocks) | **~24-26s** | 1.7-2x |
| Tuned + prefetch (Decisions 10-11) | **~19-22s** | 2-2.5x |
| Tuned + prefetch + range-scan (Opt 3) | **~15-18s** | 2.5-3x |
| LMDB/heed (mmap, estimated) | **~8-12s** | 4-5x |

### Detailed Breakdown (Tuned Config, 6M Lookups)

| Component | Calculation | Time | % of Total |
|---|---|---|---|
| DB open + pinned load | 380 MB / 3 GB/s | ~0.13s | <1% |
| Novel variant rejection | 300K × ~200 ns (pinned bloom) | ~0.06s | <1% |
| Data block I/O | 4.5M blocks × 8 KiB = 36 GB / 3 GB/s | **~12s** | **50%** |
| In-block lookups | 5.7M × ~0.7 µs (hash probe) | **~4s** | 16% |
| Decompression | 5.7M × ~1 µs (zstd dict) | **~6s** | 25% |
| Allele matching | 5.7M × ~0.3 µs | ~2s | 8% |
| `expect_hits` penalty | ~6K extra reads × 50 µs | ~0.3s | 1% |
| **Total (tuned)** | | **~24-26s** | |

After tuning, the two dominant costs are:
- **Data block I/O (~12s, 50%)** — irreducible I/O floor; target of prefetch optimizations
- **Decompression (~6s, 25%)** — CPU-bound; target of projection-aware deserialization and dictionary tuning

Prefetch overlaps these two, saving ~4-5s. Range-scan merge-join additionally eliminates per-lookup bloom/index overhead, saving another ~2-4s.

### Storage Device Sensitivity

| Device | Sequential Read | Data Block I/O (4.5M × 8 KiB) | Total (tuned) |
|---|---|---|---|
| NVMe SSD | ~3 GB/s | ~12s | ~24-26s |
| SATA SSD | ~500 MB/s | ~72s | ~84-86s |
| HDD | ~150 MB/s | ~240s | ~252s |

Data block I/O dominates so heavily after tuning that storage throughput is the primary performance determinant. NVMe is assumed throughout; SATA users should consider LMDB/heed (mmap) which benefits more from OS page cache management.

---

## Summary: Impact Assessment Matrix

| # | Optimization | Requires | Est. Savings | Complexity | Risk |
|---|---|---|---|---|---|
| 1 | Update to fjall 3.1.0 | Version bump | Recovery bugfix + enables #2 | Trivial | Very low |
| 2 | Compaction filters (eliminate `major_compact()`) | fjall 3.1.0 | ~5-15 min build time | Medium | Low |
| 3 | Range-scan merge-join for sorted VCF | fjall 3.0.2+ | ~2-4s per 6M WGS | High | Medium |
| 4 | Lower `max_cached_files` to ~128 | fjall 3.0.2+ | ~30-50 MB memory | Trivial | Low |
| 5 | `worker_threads(1)` for read-only | fjall 3.0.2+ | 3 threads, ~24 MB | Trivial | Very low |
| 6 | `size_of()` for buffer pre-sizing | fjall 3.0.2+ | ~0.1-0.3s | Low | Very low |
| 7 | `journal_compression(Lz4)` for fallback | fjall 3.0.2+ | ~30-50% journal I/O | Trivial | Very low |

**Highest-impact recommendations:**
1. **#1 + #2** — Update fjall, use compaction filters to solve the `major_compact()` problem (build-time savings)
2. **#3** — Range-scan merge-join (most significant read-path optimization not in the current design)
3. **#5** — Worker threads (free memory savings, zero risk)

---

## Changes to Existing Design Decisions

These optimizations do not change the spec requirements. They affect design decisions as follows:

| Decision | Change | Reason |
|---|---|---|
| Decision 5 (start_ingestion) | Post-ingest path changes: natural leveled compaction replaces `major_compact()` | Compaction filters (Opt 2) run during natural compaction, which respects `table_target_size(256 MiB)` |
| Risk: `major_compact()` 64 MB target | **Resolved** — `major_compact()` is not in fjall 3.1.0 public API; natural compaction with configured target is the correct path | Version upgrade eliminates the risk entirely |
| Open question: optimal look-ahead buffer size | Range-scan merge-join (Opt 3) subsumes look-ahead for the sequential case | Look-ahead (Decision 11) still valuable for the `get()` code path (sparse regions, random access) |
| Database-level settings | Add `.worker_threads(1)` to read-only open; reduce `max_cached_files` to 128 | Opts 4-5, no impact on functionality |

### Resolved Open Questions

- **`major_compact()` behavior:** Not in fjall 3.1.0 public API. Use natural leveled compaction with `table_target_size(256 MiB)` + compaction filters instead.
