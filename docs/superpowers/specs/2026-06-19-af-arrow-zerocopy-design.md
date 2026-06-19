# Arrow Zero-Copy Colocated AF — Design Spec

**Date:** 2026-06-19
**Status:** Draft (pre-implementation)
**Component:** `datafusion/bio-function-vep` — colocated existing-variant data
**Related memory:** `vep-peak-rss-allocation-churn`, `vep-lookup-streaming-and-startup-barrier`

## 1. Problem

Deep-dive (dhat, chr1 t2) found the t8 peak RSS (~8–12 GB) is **allocation
churn**, not a leak or cache:

```
dhat Total:     109 GB in 466,269,798 allocations  (~1,440 allocs/variant, 323K variants)
dhat At t-gmax: 4.15 GB live in 12,227,480 blocks   (peak live heap, t2)
dhat At t-end:  293 KB                               (no leak — everything frees)
```

mimalloc RSS (~6 GB t2) vs 4.15 GB live ⇒ ~30% is fragmentation, ~70% genuine
live working set. Both are driven by the per-variant allocation storm.

The single largest identified contributor: **colocated existing-variant data
stores allele-frequency values as `Vec<String>` (27 strings/variant), built
TWICE**:

1. Lookup collect → `ColocatedCacheEntry.af_values: Vec<String>` (`colocated.rs:35`)
   — reads 27 Arrow `Utf8` AF columns, `.to_string()` each value.
2. `build_colocated_map_from_sink` (`annotate_provider.rs:2095`, clone at
   `:2130`) → `ColocatedEntry.af_values: Vec<String>` — a second 27/variant.

≈ 54 AF-string allocations/variant × 323K ≈ **17M allocations for AF alone**,
plus the same double-copy for `variation_name`, `allele_string`, `clin_sig`,
`clin_sig_allele`, `pubmed`. The consumers (`annotate_provider.rs:1649`,
`:1962`) only ever read these **as `&str`** — owned `String`s are never needed.

## 2. Goals / Non-Goals

**Goals**
- Eliminate per-value `String` allocation for colocated AF (and the other
  colocated string fields) by referencing the Arrow arrays the lookup already
  produces, accessed by row index (zero-copy `&str`).
- Eliminate the second copy in `build_colocated_map_from_sink` (→ `Arc::clone`).
- **Byte-identical** annotation output (the bytes are unchanged, just viewed).
- Reduce dhat `Total` allocations ~10× and `t-gmax` live heap below 4.15 GB
  (t2); reduce t8 peak RSS.

**Non-Goals**
- `matched_alleles` (`Vec<MatchedVariantAllele>`, computed not from Arrow) —
  separate minor fix (SmallVec/arena), out of scope here.
- The CSQ output-formatting churn (`per_column.join("&")`, `annotate_provider.rs:2056`)
  — secondary; addressed only as an optional buffer-reuse step (Task 6).
- Changing the lookup/streaming index, cache format, or the colocated matching
  semantics.

## 3. Design

**Core idea:** the lookup's `take_rows` → `unbundle_af_columns` already returns a
`RecordBatch` whose 27 AF columns are Arrow `StringArray`s (3 allocations per
column: offsets + values buffers, `Arc`-shared) containing exactly the matched
rows. Instead of copying each variant's row into a `Vec<String>`, retain the
`Arc`'d columns and store the variant's row index.

```rust
/// The matched-rows AF columns from one take_rows batch, shared by all variants
/// collected from that batch. Project to ONLY the columns colocated needs.
#[derive(Clone)]
struct AfColumns(Arc<Vec<StringArray>>); // len == AF_COLUMNS.len()

struct ColocatedEntry {
    af: AfColumns,   // Arc::clone — refcount bump, no data copy
    af_row: u32,     // this variant's row within `af`
    // variation_name / allele_string / clin_sig / clin_sig_allele / pubmed:
    // likewise (StringArray ref + row) rather than owned String.
    ...
}

impl ColocatedEntry {
    fn af_value(&self, col: usize) -> &str {
        let a = &self.af.0[col];
        if a.is_null(self.af_row as usize) { "" } else { a.value(self.af_row as usize) }
    }
}
```

Memory: the colocated map retains the take_rows AF buffers (matched rows only,
columnar) instead of millions of scattered `String`s — same logical bytes,
packed, with ~27 allocations/batch instead of ~54/variant.

## 4. Invariants (byte-identical)

1. `af_value(col)` returns the exact bytes the current `af_values[col]` holds,
   including `""` for null/absent (today `String::new()`).
2. Field order/index mapping (`AF_COLUMNS` / `AF_COL_NAMES`) unchanged.
3. The `is_present` / "all empty" checks (`annotate_provider.rs:1671`) preserved.
4. Colocated matching (`compare_existing`, allele matching) unchanged — only the
   AF/string *storage* changes, not the matching logic.

## 5. Interfaces — confirmed by Task 0 spike (2026-06-19)

**Exact sites (all `datafusion/bio-function-vep/src/`):**

- `ColocatedCacheEntry` struct: `colocated.rs:24-36`; field `af_values: Vec<String>` at `:35`.
- `ColocatedEntry` (annotate side): `annotate_provider.rs:1599-1610`; `af_values` at `:1609`.
- `ColocatedSinkValue`: `colocated.rs:329-341`; `ColocatedSink` (Arc<Mutex<HashMap<Key,Value>>>) `:344`.
- **Collect site (the per-variant `.to_string()` loop):** `lookup_exec.rs:1133-1140`
  builds `Vec<String>` via `batch_string_value(batch, *idx, row).unwrap_or_default()`;
  pushes `ColocatedCacheEntry` at `:1141-1151`. Inside the per-row loop `for &row in rows`
  (`:1079`), under the `if let (Some(buf), Some(prepared), Some(ci))` guard (`:1093`).
  `WarmColocIndices.af_indices: Vec<Option<usize>>` (`lookup_exec.rs:388`); resolved
  once per batch by `resolve_batch_coloc_indices` (`:1986`, `coloc_indices` at `:1072`).
- **Second copy:** `build_colocated_map_from_sink` `annotate_provider.rs:2095-2144`;
  deep clone `af_values: ce.af_values.clone()` at `:2130`. Source owned (drained via
  `std::mem::take` in `drain_colocated_sink` `:2146`). Accumulates into a target map
  via `merge_colocated_delta` (`:2161`).
- **AF consumer:** the ONLY `ColocatedEntry.af_values` read is the
  `ColocatedFrequencyFields` builder `annotate_provider.rs:1962` (`idx >= af_values.len()`)
  + `:1966` (`let raw = &af_values[idx]`, then `raw.is_empty()`/`raw.split(',')`). Needs
  `&str` only. (The `:1649/:1671` sites read the separate post-aggregation
  `ColocatedVariantFields`/`ColocatedFrequencyFields`, NOT `ColocatedEntry` — unaffected.)
- **Probe to update:** `colocated_map_bytes` `annotate_provider.rs:8762-8765` (the
  `for af in &e.af_values { bytes += af.capacity() }` block).
- **`batch_string_value`** (`lookup_exec.rs:1929`) downcasts to `StringArray` /
  `StringViewArray` / `LargeStringArray`. Post-`unbundle_af_columns` (`af_bundle.rs:195`)
  AF columns are concrete `StringArray` (Utf8), **null-free, empty-string-as-absent**
  (`af_bundle.rs:123-141` appends `""` for absent). `AF_COL_NAMES` / `AF_COLUMNS` = 27 each.
- **No `take_rows` exists** — collect reads the full `&RecordBatch` indexed by `row_usize`.

**Design decisions locked from the spike:**

1. **`AfColumns(Arc<Vec<ArrayRef>>)`, not `Arc<Vec<StringArray>>`.** Storing `ArrayRef`
   + a read-time downcast lets the accessor handle all three string-array types exactly
   like `batch_string_value` (byte-identical) and avoids any per-batch conversion. Absent
   columns (`af_indices[i] == None`) get an empty `StringArray` placeholder; the accessor
   returns `""` for null / `row >= len` / absent / non-string. 27 entries, indexed by
   `AF_COL_NAMES` order.
2. **AF field only in this slice.** `variation_name`/`allele_string`/`clin_sig`/
   `clin_sig_allele`/`pubmed` stay `String`/`Option<String>` for now — `variation_name`
   is a dedup HashMap key (`:2107/:2120`) and the rest are `Option<String>`, adding risk
   for ~1/3 of the win. Deferred to the §7 follow-up, gated on G3.

## 6. Risks

- **R1 — retained batch over-retention (ELEVATED by spike).** Spec §3 assumed
  `take_rows` yields matched-rows-only arrays. It does NOT — there is no `take_rows`;
  the collect site reads the full `&RecordBatch`. So `Arc`-retaining a column keeps
  ALL of that batch's rows alive, not just the matched ones, for as long as any entry
  from that batch lives in the (accumulating) colocated map. Mitigation taken: build
  `AfColumns` from only the 27 AF column `ArrayRef`s (one `Arc::clone` each — the whole
  batch / other columns are NOT retained), and entries from the same batch+call SHARE
  one `AfColumns` Arc. **Decision: accept full-batch AF-column retention** (it is the
  only approach that meets the primary G2 allocation-count goal: ~2 allocs/batch vs
  ~54 String allocs/variant), and **measure peak in Task 5 (G3)**. Fallback if G3
  regresses: per-call `arrow::compute::take` to gather matched rows into compact
  columns (bounds retention, but ~27 take ops/call ≈ no count win for the common
  small-match case — hence not the default).
- **R2 — exact collect/sink sites unknown.** The precise code building
  `ColocatedCacheEntry.af_values` from Arrow and the sink drain path are not yet
  read. Task 0 spike locates them.
- **R3 — null vs empty.** Today absent → `String::new()` (`""`). Arrow `Utf8`
  may represent absent as null OR empty-string; `af_value` must map both to `""`.
- **R4 — lifetime across windows.** Colocated map outlives individual take_rows
  batches; the `Arc` refs keep arrays alive correctly (intended), but confirm no
  `'a` borrow leaks into the map (must be owned `Arc`, not `&`).

## 7. Success criteria / gates

- **G1 Parity:** chr1 t8 e2e — 0 CSQ mismatch, 86/86 fields 100%.
- **G2 Allocations:** dhat `Total` blocks ↓ ≥5× (target ~10×); `t-gmax` live heap
  < 4.15 GB (t2). Measured via the `dhat-heap` vepyr feature.
- **G3 Peak RSS:** chr1 t8 peak RSS materially below current ~8–12 GB.
- **G4 No regressions:** `cargo test -p datafusion-bio-function-vep` green;
  clippy clean.

## 7b. Measured results (2026-06-19, chr1 merged, post-implementation)

Commits: `39dcc68` (AfColumns, pre-swap) → `f0dbb1b` (zero-copy swap wired in).

- **G1 Parity (t8):** ✅ PASS. ALL 86 shared CSQ fields 100%, 0 mismatch over
  4,737,090 rows (every AF column, MAX_AF, MAX_AF_POPS, CLIN_SIG, …). Byte-identical.
- **G4 Tests/clippy:** ✅ 756 lib tests pass; clippy `-D warnings` clean.
- **G3 Peak RSS (t8):** ~7.9 GB (in-proc getrusage) / 8.3 GB (`/usr/bin/time -l`),
  vs prior ~8–12 GB. **Modest / net-neutral.** Root cause: this change trades many
  small live `String`s for fewer, larger Arc-retained full-batch AF columns (R1).
  Churn/CPU ↓ but live working set ~flat → peak barely moves. (The `colocated_MB`
  probe over-counts shared Arcs — not a leak.)
- **G2 Allocations (dhat, t2):** NOT completed. dhat heap instrumentation slows the
  pipeline ~700× (≈78 var/s vs ≈50k), so a full chr1 t2 run is ~70 min and the
  Total/t-gmax summary only prints on profiler drop (no partial read). Killed at ~28%.
  **Predicted (not measured):** Total ↓ by the AF portion only — ~54 allocs/variant of
  ~1,440/variant ≈ **~4%** of the 466M baseline; t-gmax ~flat (consistent with G3).
- **Regression sweep (t1–t8 A/B, before `39dcc68` vs after `f0dbb1b`, probes dormant):**
  ✅ No regression. Throughput after ≈ before within <1% at every thread count
  (t1 12.6k, t2 28.5k, t4 44.8k, t8 63.8k var/s); scaling intact (~5.1× t1→t8).

**Conclusion.** The refactor is correct, byte-identical, and throughput-neutral, and it
removes a real (if small, ~4%) source of allocation churn. But AF zero-copy is a
**churn/CPU lever applied to a peak-RSS-bound problem** — so it does not move the
headline RSS number, and the R1 full-batch retention slightly offsets the churn win on
peak. **Keep the commit** (no downside, clean foundation). Next lever for peak/throughput
is the **engine CSQ output path** (per-variant `Vec<String>`/`join` → reused byte buffer /
byte-chunk output), which is both the dominant churn source (~the rest of the 1,440/variant)
and a live-set reducer. If peak RSS specifically is the target, the deferred
`arrow::compute::take` matched-rows gather would *reverse* R1 and cut retention.

## 8. Rollout

Internal change, public API unchanged → no vepyr edit, rebuild only. No on-disk
format change. Land behind the parity gate; no feature flag needed (byte-
identical). Keep the `[VEP_RSS]` + dhat probes for before/after measurement.
