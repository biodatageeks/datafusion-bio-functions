# Streaming + Selective Variation BTree Index — Design Spec

**Date:** 2026-06-19
**Status:** IMPLEMENTED — streaming is the only production path; the eager
materialized index is removed (kept solely as a `#[cfg(test)]` parity oracle).
The `VEP_STREAMING_VARIATION_INDEX` flag and eager fallback described in §9 were
removed after the t1–t8 parity gate passed (per user decision). Byte-identical
on chr1 t1–t8.
**Component:** `datafusion/bio-function-vep` — Lance variation lookup runtime
**Related memory:** `vep-lookup-streaming-and-startup-barrier`, `vep-t8-rss-root-cause`

## 1. Problem

The per-contig variation lookup builds its position index **eagerly and fully**
before the first probe. `SinglePathLanceVariationLookup::open`
(`variation_runtime.rs:27`) calls `load_start_index_from_lance_btree`, which
`load_u32_btree_index` (`row_index.rs:246`) implements by reading the **entire**
BTree leaf (`page_data.lance`) into a `Vec<(u32,u64)>`, sorting it, and building
two parallel `Vec`s held resident for the contig's lifetime.

Measured cost (chr1, merged cache, `VEP_PROFILE=1 VEP_LANCE_PROFILE=1`):

```
lookup_wait               0.948 s   (worker blocked on first lookup batch)
└─ open (state build)     0.704 s   row_ids=88,153,966  unique_positions=83,219,731
   ├─ dataset_open        0.000 s
   └─ index_load          0.704 s
      ├─ segments         0.000 s
      ├─ pairs_read       0.395 s   decode 88.15M (pos,row_id) pairs from page_data
      ├─ sort             0.045 s   (already globally sorted → near-linear)
      └─ build vecs+valid ~0.264 s
resident index RAM        ~1.06 GB  (88.15M×4B + 88.15M×8B)
transient build RAM       ~1.41 GB  (88.15M × 16B pairs Vec) → ~2.5 GB peak
```

This single structure is **both** the startup stall (0.70 s of 0.95 s
`lookup_wait`) **and** the ~2.6 GB/partition RSS documented in
`vep-t8-rss-root-cause`. The VCF being annotated has ~323 K rows — we build an
83 M-position index to answer ~323 K queries.

## 2. Goals / Non-Goals

**Goals**
- Eliminate the eager full-materialization. Resident index RAM target: **< 5 MB
  per partition** (down from ~1.06 GB).
- Remove the upfront build barrier from `lookup_wait`: first output emits after
  reading only the pages covering the first VCF batch, not the whole index.
- When VCF partitions are contiguous position ranges, each partition reads only
  its band of `page_data` (≈ 1/N of 682 MB), in parallel.
- **Byte-identical** annotation output (0 CSQ mismatches on the chr1 e2e gate).

**Non-Goals**
- Reducing *total* leaf bytes read across all partitions (inherent for uniform
  full-coverage queries — confirmed acceptable in design review).
- Changing the cache on-disk format (no rebuild; read-only change).
- Changing `take_rows` (main-data fetch) — already per-match, untouched.
- Persisting a prebuilt mmap sidecar (a later, orthogonal optimization).

## 3. Background: Lance BTree physical layout (verified, lance-index 7.0.0)

A scalar BTree index is two files under `<dataset>/_indices/<uuid>/`:

| File | Size (chr1) | Contents |
|---|---|---|
| `page_lookup.lance` | 192 KB | sub-index: `BTreeMap<min_value → PageRecord{max, page_number}>` (`btree.rs:617`); one record per page |
| `page_data.lance` | 682 MB | leaf: `(values=position u32, ids=row_id u64)`, globally sorted by position, in pages of `DEFAULT_BTREE_BATCH_SIZE = 4096` rows (`btree.rs:74`) |

- `page_data` is **globally sorted by position** (confirmed: `sort` step on 88 M
  rows is 0.045 s → input already ordered).
- `IndexReader::read_record_batch(page_number, batch_size)` (`scalar.rs:195`)
  reads exactly **one page**; Lance navigates `page_lookup` → page, then reads
  only that page (`btree.rs:924-928`). Selective reads are first-class.
- `BTreeLookup::pages_in_range((lo, hi))` (`btree.rs:697`) returns the page
  numbers overlapping a value range.

The current code bypasses `page_lookup` entirely and bulk-reads `page_data`.

## 4. Design

Three orthogonal mechanisms compose. They are **independent** — clarified in
review:

| Mechanism | Cuts | Does NOT cut |
|---|---|---|
| **Selective seek** (page_lookup → first page) | IO/partition → 1/N, no upfront barrier | resident RAM |
| **Lazy buffered advance** (5-page window, discard consumed) | resident RAM → MB | total IO |
| **page_lookup re-read** (back-step floor) | — (correctness) | — |

### 4.1 Selective seek

The VCF positions themselves drive the cursor — no explicit `[w_start, w_end)`
plumbing needed. On the first `resolve()` call, binary-search the in-memory page
directory for the batch's minimum query position to get the starting page
number, then `read_record_batch(page)` from there. As query positions advance,
pages are pulled forward. A partition whose VCF covers a contiguous sub-range
thus touches only its band; round-robin partitions degrade gracefully to lazy
(RAM win retained, IO win lost).

### 4.2 Lazy buffered advance with a fixed 5-page window

The cursor holds a sliding ring buffer of **5 pages**: the page containing the
cursor (P) and ±2 neighbours `{P-2, P-1, P, P+1, P+2}`. On crossing into the
next page: evict the trailing page, prefetch the new leading page.

```
page = 4096 rows ≈ 11.6 kb genomic (chr1 density 0.354 rows/bp)
behind (−2) ≈ 23 kb  ≫  real back-step ~32 bp (build_probe_starts del-shift)
ahead  (+2) ≈ 23 kb  →  read-ahead depth (latency hiding)
resident: 5 × 4096 × (4B+8B) = 240 KB / partition
```

The backward bound is from `build_probe_starts` (`lookup_exec.rs:2375-2400`):
deletion-repeat shift is `del_len.min(32)` upward and `saturating_sub(1)`
downward — i.e. the only backward motion across sorted batches is ~tens of bp,
≪ one page. Two trailing pages is ~700× margin.

### 4.3 page_lookup re-read floor (correctness, not perf)

If a query position ever falls before the oldest retained page (pathological:
unsorted VCF, oversized SV), binary-search `page_lookup` → that page →
`read_record_batch` (one ~32 KB page). Log it (so the window can be tuned). The
window is the fast path; this guarantees correctness regardless of window size.

## 5. Invariants (byte-identical output)

1. **Per-position row_id order.** `load_u32_btree_index` sorts pairs by
   `(position, row_id)`, so within a matched position the materialized resolve
   emits row_ids **ascending**. The streaming cursor MUST sort the row_ids it
   collects for each matched position ascending before emitting. (Page-stored
   order within equal positions is not guaranteed row_id-ascending.)
2. **Match semantics identical** to
   `PositionRowIdIndex::resolve_sorted_positions_from_cursor`
   (`row_index.rs:97`): for each query position, emit all row_ids whose position
   equals it; advance monotonically.
3. **Forward-monotonic cursor** with bounded backward tolerance (5-page window);
   any out-of-window back-step resolved via 4.3, never silently dropped.
4. **Edge clamp:** at the band's first/last page fewer than 5 pages exist; clamp
   without error.

## 6. Interfaces

### 6.1 New types (`row_index.rs`)

```rust
/// In-memory page directory loaded from `page_lookup.lance`.
/// Maps a u32 position to the page_number whose [min,max] may contain it.
pub struct PositionPageDirectory {
    /// Per-page minimum position, ascending. page_mins[i] is page i's min.
    page_mins: Vec<u32>,
    /// Total pages and rows-per-page, for offset/edge math.
    num_pages: u32,
    batch_size: u64,
    num_rows: usize,
}
impl PositionPageDirectory {
    pub async fn load(dataset: &lance::Dataset) -> Result<Self>;
    /// First page that may contain `position` (partition_point on page_mins).
    pub fn first_page_for(&self, position: u32) -> u32;
    pub fn num_pages(&self) -> u32;
}
```

### 6.2 New type (`variation_runtime.rs`)

```rust
/// Per-partition forward streaming cursor over page_data with a 5-page window.
pub struct StreamingPositionCursor {
    reader: Arc<dyn lance_index::scalar::IndexReader>, // page_data reader (shared)
    directory: Arc<PositionPageDirectory>,             // shared
    /// Ring buffer of up to WINDOW_PAGES decoded pages, ascending page_number.
    window: VecDeque<DecodedPage>,
    /// Index into the *front* page of `window` of the next unconsumed row.
    front_off: usize,
    seeded: bool,
}
struct DecodedPage { page_number: u32, positions: Vec<u32>, row_ids: Vec<u64> }

pub const WINDOW_PAGES: usize = 5; // cursor page ± 2
```

### 6.3 Changed signatures

```rust
// variation_runtime.rs — open no longer materializes the index.
impl SinglePathLanceVariationLookup {
    pub async fn open(path: &Path, projection: Vec<String>) -> Result<Self>;
    // fields: dataset, projection, page_data_reader: Arc<dyn IndexReader>,
    //         directory: Arc<PositionPageDirectory>   (was: index: PositionRowIdIndex)
    pub fn new_cursor(&self) -> StreamingPositionCursor;
    pub async fn resolve_and_take(
        &self,
        sorted_unique_starts: &[u32],
        cursor: &mut StreamingPositionCursor,    // was: &mut usize
    ) -> Result<TakenVariationRows>;
}
```

```rust
// lookup_exec.rs — per-partition cursor type changes.
// line 358:
lance_cursors: HashMap<String, StreamingPositionCursor>,  // was HashMap<String, usize>
// line 1704 obtains/creates the cursor via lookup.new_cursor() on first use.
```

`ResolvedRowIds`, `TakenVariationRows`, and `take_rows` are unchanged.

## 7. Risks / open questions (resolved in Task 1 spike)

- **R1 — page_lookup schema.** Exact column names/types in `page_lookup.lance`
  are not yet confirmed from disk. Spike: open it, log schema + sample rows.
- **R2 — partition locality.** Whether VCF partitions are contiguous position
  ranges (selective win) or round-robin (lazy-only). Spike: log per-partition
  min/max start on chr1 forks=4. Design degrades gracefully either way.
- **R3 — async reads under `block_in_place`.** `resolve_and_take` is already
  called inside `block_in_place(block_on(...))` (`lookup_exec.rs:1706`); the
  cursor's per-page `read_record_batch` awaits must run in that same context.
  Spike: confirm a page read works from the cursor in a `block_in_place` test.

## 8. Success criteria / gates

- **G1 Parity:** chr1 merged e2e — `CSQ count mismatch: 0`, all 86 shared CSQ
  fields 100% (the gate already used this session).
- **G2 RAM:** resident index RAM < 5 MB/partition (from ~1.06 GB); peak process
  RSS drop measured on chr1 at forks ∈ {0, 4, 8}.
- **G3 Latency:** `lookup_wait` reduced; `open` build cost ≈ 0; first-output time
  drops. Reported via the existing `[VEP_PROFILE]` line.
- **G4 No regressions:** `cargo test -p datafusion-bio-function-vep` green;
  existing `variation_runtime.rs` tests pass unchanged.
- **G5 End-to-end wall time:** total annotation wall-clock reduced by removing
  the serial pre-output build barrier. The eager build is per-contig and
  serialized before that contig's first output, so its cost is on the critical
  path; streaming overlaps it with annotation. Expected saving ≈ the eager
  `open` build time that is *not* otherwise hidden — order 0.7 s/contig on chr1
  (more when context load is shorter than the build, less when it fully hides
  the build). Measure whole-run wall time eager (`=0`) vs streaming (`=1`) on
  chr1 single-contig (forks=0, so no cross-contig overlap masks it) and on a
  multi-contig run (e.g. chr21+chr22) where per-contig barriers accumulate.
  Target: streaming wall time ≤ eager, by ≈ Σ(per-contig unhidden build time).

## 8a. Results (2026-06-19, chr1 merged, t1–t8)

Benchmark: chr1, `--threads 1..8`, eager (`=0`) vs streaming (`=1`), full e2e
(bcftools norm + annotate + VEP diff). Wall = `/usr/bin/time` real; RSS = peak
process resident.

| threads | streaming wall | eager wall | streaming lookup_wait | eager lookup_wait | parity |
|---|---|---|---|---|---|
| 1 | 169.4 s | 166.5 s | **0.244 s** | **0.972 s** | PASS (both) |
| 2 | 150.5 s | — | 0.000 s | — | PASS |
| 3 | 150.1 s | — | 0.000 s | — | PASS |
| 4 | 146.7 s | 151.2 s | 0.000 s | 0.000 s | PASS |
| 5 | 148.2 s | — | 0.000 s | — | PASS |
| 6 | 148.1 s | — | 0.000 s | — | PASS |
| 7 | 150.3 s | — | 0.000 s | — | PASS |
| 8 | 149.0 s | 149.7 s | 0.000 s | 0.000 s | PASS |

- **G1 Parity: PASS at every t1–t8** (0 CSQ mismatch, all 86 fields 100%). The
  hard gate — streaming is byte-identical.
- **G3 lookup_wait:** at t=1 (serial, single partition) the build-wait drops
  **0.972 s → 0.244 s**. At t≥2 *both* read 0.000 s — the within-contig parallel
  pipeline already overlaps the eager build behind context load
  (`prepare_total ≈ 1.07 s`), so the eager barrier was already hidden there. The
  wall-time barrier therefore only surfaces at t=1 on a single contig (and on
  the first partition of each contig in multi-contig runs).
- **G5 wall time:** within run-to-run noise (~147–169 s, harness-dominated by
  norm+diff). No regression — streaming ≤ eager at t4/t8, ≈ eager at t1.
- **G2 RAM:** whole-process peak RSS is too noisy (7.5–10.5 GB, ±2 GB
  run-to-run) to isolate the ~1 GB variation-index delta. **Important negative
  finding:** peak RSS (~10 GB) is *not* driven by the variation index — it is
  dominated by SIFT/AF/context/output — so eliminating the index's ~1 GB
  resident allocation does not move *peak* RSS measurably here. The resident
  reduction (1.06 GB → ~240 KB/partition) is real and proven structurally + by
  unit tests (the 88 M-pair `Vec`s are never allocated), but is not the t8 peak
  driver in the current OnceCell-shared configuration.

**Verdict:** byte-identical, no regression, removes the eager index build (the
t=1 startup barrier) and ~1 GB of resident allocation. It does **not** by itself
cut steady-state wall time at t≥2 (the eager build was already hidden) or the t8
*peak* RSS (driven by other components). Default flip is therefore a judgment
call (see §9).

## 9. Rollout

- Internal runtime change; public `annotate_vep` API unchanged → no vepyr API
  edit, rebuild only.
- Gate behind `VEP_STREAMING_VARIATION_INDEX` env (default ON after G1–G4 pass;
  fallback to the eager `load_start_index_from_lance_btree` path retained for
  one release as `=0`).
- Keep the eager path code until G1–G4 verified on ≥2 contigs (chr1 + chr22).

## 10. Task 1 spike findings (2026-06-19 — all risks cleared)

Confirmed against a real indexed dataset (20,002 rows → 5 pages) via
`spike_dump_btree_index_structure` and a crate-wide code read.

**R1 — `page_lookup.lance` schema (CONFIRMED):** columns
`min: UInt32 (nullable)`, `max: UInt32 (nullable)`, `null_count: UInt32`,
`page_idx: UInt32`. One record per page (5 records for 5 pages). The page size
is in the schema **metadata**: `batch_size = "4096"`, `range_partitioned =
"false"`. → `PositionPageDirectory::load` reads `min` + `page_idx`, and SHOULD
read `batch_size` from `lookup.schema()` metadata rather than hardcoding
`DEFAULT_BTREE_BATCH_SIZE` (fallback to 4096 if absent).

**`page_data.lance` schema (CONFIRMED):** columns `values: UInt32`,
`ids: UInt64` (= `BTREE_VALUES_COLUMN`/`BTREE_IDS_COLUMN`).
`read_record_batch(0, 4096)` returns exactly **4096 rows = one page** →
per-page selective reads work as designed. Cursor decodes `values`/`ids` by
name.

**R3 — async read under `block_in_place` (CONFIRMED):** a per-page
`read_record_batch` driven inside `block_in_place(Handle::block_on(..))`
returns 4096 rows — matches the call context at `lookup_exec.rs:1706`. No
runtime-nesting issue.

**R2 — partition locality (RESOLVED: contiguous ranges).** `KvLookupExec`
inherits its VCF input's partitioning (`lookup_exec.rs:155,285`); **no
`RoundRobinBatch`/`RepartitionExec` exists anywhere in the crate** (grep clean).
The lookup runs as *"N spawned, position-ordered partitions"*
(`annotate_provider.rs:5224`) each writing *"their own position-ordered VCF"*
shard (`:9372`), and `window_planner` groups *"a contiguous run of
position-sorted batches"*. Position-ordered shard concatenation requires each
partition to own a contiguous, disjoint, ascending position range → **selective
1/N reads are achievable** (not just lazy). Round-robin would break sorted
shard output, which the e2e gate validates, so range-partitioning holds by
construction.

**Spike-driven plan adjustments:**
1. Task 2 `PositionPageDirectory::load`: column names `min`/`page_idx` are
   correct as written; additionally read `batch_size` from the page_lookup
   schema metadata (key `"batch_size"`), default 4096.
2. `min` is nullable in the schema; positions are non-null in practice, but the
   loader must downcast tolerantly (use `value(i)` only for non-null rows).
3. `load_btree_segments` was made `pub(crate)` (used by Task 2/3); retained.
