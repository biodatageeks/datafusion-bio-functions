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

## 5. Interfaces (subject to Task 0 spike confirmation)

- `colocated.rs`: `ColocatedCacheEntry.af_values: Vec<String>` →
  `af: AfColumns, af_row: u32` (+ accessor). Same for the other string fields.
- `colocated.rs`: `ColocatedSinkValue` carries `Arc`'d arrays + row index.
- `annotate_provider.rs`: `ColocatedEntry` mirrors the above; `build_colocated_
  map_from_sink` becomes `Arc::clone` not deep clone; consumers call
  `entry.af_value(idx)` instead of `&entry.af_values[idx]`.
- Lookup collect site (variant_lookup_exec.rs / lance_cache/lookup_exec.rs):
  store array refs + row index instead of `.to_string()`.

## 6. Risks

- **R1 — retained batch over-retention.** Holding `Arc<StringArray>` keeps that
  column alive; ensure only the needed AF (+ string) columns are retained, not
  the whole take_rows batch. Mitigation: build `AfColumns` by projecting/cloning
  only the required columns (cheap — column `Arc` clone).
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

## 8. Rollout

Internal change, public API unchanged → no vepyr edit, rebuild only. No on-disk
format change. Land behind the parity gate; no feature flag needed (byte-
identical). Keep the `[VEP_RSS]` + dhat probes for before/after measurement.
