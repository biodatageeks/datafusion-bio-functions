# Fjall Path Optimization Analysis

## 1. Projection-Aware Deserialization (Late Column Skipping)

**Current state:** Every position entry is fully decompressed, then `PositionEntryReader::new()` walks ALL columns to compute `col_offsets` (`position_entry.rs:306-317`), even when only 3-5 of 78 columns are requested.

**Problem:** The `new()` constructor iterates through all `num_cols` column headers to build `col_offsets`, calling `column_packed_size()` for each — including string columns which read 4 bytes of `total_string_bytes`. This is ~78 column boundary calculations per lookup, when the user typically requests 5-10.

**Proposed fix:** Lazy column offset computation — only compute offsets up to the max requested column index, or compute on-demand.

| Approach | Effort | Impact |
|----------|--------|--------|
| Skip-to-column: stop offset walk at max needed col | Low | ~10-15% faster deserialization for typical 5-col queries |
| Two-tier layout: allele table + column directory (pre-computed offsets in header) | Medium | O(1) column access, eliminates sequential walk entirely |

**Impact:** For 5M WGS variants, saving ~70 unnecessary column walks per lookup = ~350M column boundary skips avoided.

---

## 2. Column-Group Split Layout (FORMAT_V1)

**Current state:** Each position entry stores ALL 78 columns in one zstd-compressed blob. Every point lookup decompresses everything, even for a 1-column query.

**Problem:** With 78 VEP columns, a typical `variation_name + consequence_type + impact` query still decompresses the full blob (~100-500 bytes after zstd, ~200-2000 bytes uncompressed). For the colocated path (`cache_exec.rs:652-801`), BOTH the primary match AND the colocated collection deserialize the same entry twice with different column sets.

**Proposed fix — Column-group split layout (FORMAT_V1):**

| Group | Columns | Use case |
|-------|---------|----------|
| **Core** (key value) | allele_string, end, failed | Always needed — allele matching + filtering |
| **Annotation** | consequence_type, impact, variation_name, clin_sig, ... | Primary lookup output |
| **Frequency** | AF, gnomADe, gnomADg, ... | Colocated / frequency queries |
| **Extended** | All remaining low-use columns | Rare |

Store as separate fjall keys: `[chrom_code][start][group_id]`. Each group independently compressed. The "Core" group (~20 bytes) is always fetched; annotation/frequency groups only on match.

| Metric | Current (all-in-one) | Column-group split |
|--------|---------------------|-------------------|
| Bytes decompressed per miss | 100-500 (full entry, realized bloom rejects most) | 20-40 (core only) |
| Bytes decompressed per hit | 100-500 | 20-40 (core) + 50-200 (matched group) |
| Fjall point lookups per hit | 1 | 2 (core + needed group) |
| Compression ratio | Best (all columns share dict) | Slightly worse per group |

**Impact:** ~30-50% less decompression work on hits (most bytes are in unused frequency/extended columns). Trade-off: extra point lookup per hit, but fjall's block cache makes second lookup nearly free for adjacent keys.

---

## 3. Key Buffer Allocation on Hot Path

**Current state** (`kv_store.rs:246-256`):
```rust
pub fn get_position_entry(&self, chrom_code: u16, start: i64) -> Result<Option<fjall::UserValue>> {
    let mut key_buf = Vec::with_capacity(10);  // ALLOCATION per lookup
    encode_position_key_buf(chrom_code, start, &mut key_buf);
    self.data.get(&key_buf)
}
```

**Problem:** `Vec::with_capacity(10)` allocates 10 bytes on the heap for every single lookup. For 5M variants with ~2-4 probes each = 10-20M heap allocations.

**Fix:** Use a stack-allocated `[u8; 10]` buffer instead:
```rust
let mut key_buf = [0u8; 10];
key_buf[0..2].copy_from_slice(&chrom_code.to_be_bytes());
key_buf[2..10].copy_from_slice(&start.to_be_bytes());
self.data.get(&key_buf)
```

**Impact:** Eliminates ~10-20M heap allocations per WGS annotation. Low effort, measurable throughput gain.

---

## 4. Batch Prefetch / Multi-Get

**Current state:** Each VCF row triggers 1-4 sequential `store.get_position_entry_fast()` calls (`cache_exec.rs:588`). For sorted VCF input, consecutive rows often hit adjacent positions on the same chromosome.

**Problem:** Sequential point lookups can't exploit OS read-ahead or fjall's block cache warming as effectively as batch access.

**Proposed fix:** Buffer N rows (e.g., 64-128), collect all probe keys, sort them, then iterate lookups in key order. This maximizes block cache reuse since fjall stores data in sorted LSM order.

Alternatively, fjall doesn't expose a native `multi_get()`, but you could:
- Pre-sort the probe keys within each batch before looking them up
- Group probes by chrom_code so all lookups for the same chromosome are adjacent

**Impact:** Estimated 15-25% reduction in I/O for cold-cache scenarios. Low impact when block cache is warm.

---

## 5. Context Tables Still in Parquet — Migration Candidates

**Current state:** Only the variation cache is in fjall. These tables are still loaded from Parquet via SQL:

| Table | Size (chr1) | Access Pattern | Fjall Suitability |
|-------|-------------|----------------|-------------------|
| **Transcripts** | 32 MB | Interval overlap by chrom+start/end | Medium — range-scan needed |
| **Exons** | 7 MB | By transcript_id after transcript load | Good — point lookup by transcript_id |
| **Translations** | 229 MB | By transcript_id (biggest bottleneck) | **Excellent** — point lookup by transcript_id |
| **Regulatory** | Small | Interval overlap | Low — already uses interval filter SQL |
| **SIFT/PolyPhen** | In translations | Windowed by position | Already has `SiftKvStore` |

### High-value migration target: Translations

Current bottleneck: 229MB for chr1, loaded via `SELECT * FROM translations WHERE chrom IN (...)` — even for a single cache miss.

**Proposed KV layout for translations:**
- Key: `[transcript_id bytes]` (variable-length string key)
- Value: zstd-compressed struct with `cds_len, protein_len, translation_seq, cds_sequence, protein_features`
- Access: After transcript interval join identifies relevant transcript_ids, do point lookups

**Impact:** Eliminates the 15.1s context-table loading bottleneck (from profiling data). Cold-start annotation for 1000 variants drops from ~40s to ~25s.

### Secondary target: Exons

- Key: `[transcript_id bytes]`
- Value: list of `(start, end, strand)` tuples
- Loaded after transcript interval join identifies relevant transcript_ids

---

## 6. Colocated Collection Redundancy

**Current state** (`cache_exec.rs:652-801`): For each probe position, the colocated collection path re-iterates ALL alleles at that position, re-reads `failed`, `end`, `variation_name`, etc. — many of which were already read in the primary match loop above (`cache_exec.rs:607-627`).

**Problem:** Double deserialization of the same `PositionEntryReader` fields. For colocated-heavy queries, this nearly doubles the per-entry work.

**Fix:** Merge the primary match loop and colocated collection into a single pass. During the primary allele iteration, collect colocated data for non-matched alleles simultaneously.

**Impact:** ~20-30% less deserialization work when colocated sink is active. Medium effort (refactor the two loops into one).

---

## 7. String Allocations in Colocated Path

**Current state** (`cache_exec.rs:659-663`):
```rust
let chrom_norm = chrom.to_string();           // alloc per row
let input_allele_string = format!("{input_ref}/{input_alt}");   // alloc per row
let compare_allele_string = format!("{compare_ref}/{compare_alt}"); // alloc per row
```

These are inside the `for row in 0..num_rows` loop. For 5M rows, that's 15M+ string allocations.

**Fix:** Pre-allocate and reuse `String` buffers:
```rust
let mut allele_buf = String::with_capacity(32);
// In loop:
allele_buf.clear();
write!(allele_buf, "{input_ref}/{input_alt}").unwrap();
```

**Impact:** Reduces GC pressure and allocation overhead. Low effort.

---

## 8. Sorted Ingestion API (`start_ingestion()`)

**Current state** (`loader.rs`): Uses `batch_insert_raw()` via `OwnedWriteBatch`, which routes through memtable -> flush -> compaction.

**Problem:** The design spec (Decision 5) explicitly recommends `Keyspace::start_ingestion()` for sorted data, estimating 3-6x faster build time. This is not yet implemented.

**Fix:** Pre-sort entries by key within each batch, then use the bulk ingestion API:
```rust
let mut ingestion = store.data.start_ingestion()?;
for (key, value) in sorted_entries {
    ingestion.write(&key, &value)?;
}
ingestion.finish()?;
```

**Impact:** Cache build time drops from ~hours to ~tens of minutes for 1.17B variants. High value for initial setup.

---

## 9. Per-Batch HashMap Grouping in Loader

**Current state** (`loader.rs:450-458`): Each batch creates a new `HashMap<(String, i64), Vec<usize>>` for position grouping, allocating String keys per row.

**Problem:** If Parquet data is pre-sorted by (chrom, start), the HashMap is unnecessary — consecutive rows share the same key.

**Fix:** Detect sorted input and use a streaming group-by (track previous key, flush on change). Falls back to HashMap for unsorted input.

**Impact:** Eliminates O(batch_size) String allocations per batch during loading. Medium effort.

---

## Prioritized Summary

| # | Optimization | Impact | Effort | Category |
|---|-------------|--------|--------|----------|
| 3 | Stack-allocated key buffer | High throughput | Very low | Hot path |
| 7 | Reuse string buffers in colocated | Medium throughput | Low | Hot path |
| 1 | Lazy column offset computation | Medium throughput | Low | Deserialization |
| 6 | Merge primary + colocated loops | Medium throughput | Medium | Hot path |
| 5 | Migrate translations to fjall | **Very high** (eliminates 15s bottleneck) | Medium-high | Architecture |
| 8 | Sorted ingestion API | High (3-6x build speed) | Medium | Ingestion |
| 4 | Batch/sorted prefetch | Medium I/O | Medium | Hot path |
| 2 | Column-group split layout (V1) | High (30-50% less decompress) | High | Layout redesign |
| 9 | Streaming group-by in loader | Medium (build only) | Medium | Ingestion |

Items 3, 7, 1 are quick wins. Item 5 (translations -> fjall) is the single highest-impact architectural change. Item 2 (column-group layout) is the most impactful layout change but requires a FORMAT_V1 migration.
