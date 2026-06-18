# Transcript-Clone Elimination Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Eliminate the unconditional per-buffer `TranscriptFeature` clone in the within-contig parallel annotation path, replacing it with a copy-on-write overlay so unmodified transcripts are served as shared references into `base_transcripts` and only genuinely-mutated transcripts are cloned.

**Architecture:** `BufferAnnotationContext.transcripts` changes from `Vec<TranscriptFeature>` (always cloned) to `Vec<CowTranscript>`, where `CowTranscript` is `Borrowed(usize)` (an index into `shared.base_transcripts`) or `Owned(Box<TranscriptFeature>)` — the same borrowed-or-owned shape as `std::borrow::Cow`, but storing an index instead of a reference (which is what keeps `BufferAnnotationContext` free of a borrow lifetime and lets the stored struct outlive the build call). The buffer build starts every transcript `Borrowed` and converts to `Owned` (via `to_mut`, which clones from base) only on a write that actually changes a field — using compare-before-write so today's idempotent assignments (e.g. resetting `gene_hgnc_id` to a value it already holds) never trigger a clone. `prepared_context_from_buffer` resolves each `CowTranscript` to a `&TranscriptFeature`, so `PreparedContext` (`Vec<&'a TranscriptFeature>`) and the entire `transcript_consequence.rs` engine are **unchanged**. `Borrowed` transcripts also share the base instance's geometry `OnceLock` (built once, not rebuilt per clone).

**Tech Stack:** Rust 2024, DataFusion 52.1, Arrow 57. Single crate `datafusion-bio-function-vep`. Bench harness `examples/bench_annotate_vcf.rs` (mimalloc global allocator — required for the perf numbers). Backend: lance variation cache.

**Why this approach (measured):** Instrumenting the real build path on chr1 200k `--everything` showed **75.9% of buffer transcripts are byte-identical to their `base_transcripts` counterpart after all mutations** (HGNC donation, `gene_symbol`/`gene_symbol_source` fill, sequence overrides, persisted reuse), so ~76% of clones are eliminated; the result is stable across threads=1 (75.87%) and threads=8 (75.94%). Source problem doc: `docs/superpowers/plans/2026-06-18-transcript-clone-elimination.md`.

**Risk note:** This is the same code region as the recently-fixed partition-boundary HGNC parity bug. Output MUST stay byte-identical at threads ∈ {1,2,4,8}. Because `CowTranscript::Owned` is a *verbatim* clone of the base transcript that is then mutated exactly as today, parity is preserved by construction: the only thing the compare-before-write logic can get wrong is taking an unnecessary clone (a perf miss, never a wrong byte). The mandatory verification gate (Task 9) is the byte-identical `diff` from the problem doc.

**TDD adaptation:** The externally-observable behavior is *unchanged* (byte-identical), so the primary "test" is the existing unit-test suite plus the byte-identical parity diff — not new assertions on new behavior. Where a task introduces a new pure helper (`CowTranscript` methods, compare-before-write donor application), it gets a focused unit test written first. Integration tasks are verified by `cargo test` + the parity gate.

---

## File Structure

All changes are in one file: `datafusion/bio-function-vep/src/annotate_provider.rs`. No changes to `transcript_consequence.rs` (the engine), to `PreparedContext`, or to the output formatter — that is the central property of this design.

Functions/types touched (current line numbers, will drift as edits land):

| Symbol | Current line | Change |
|---|---|---|
| `CowTranscript` (new enum) | — | new type + impl near `BufferAnnotationContext` (~9754) |
| `BufferAnnotationContext.transcripts` | 9755 | type `Vec<TranscriptFeature>` → `Vec<CowTranscript>` |
| `apply_hgnc_donors` | 10947 | operate on `&mut [CowTranscript]` + base, compare-before-write |
| `collect_hgnc_donors` | 10911 | unchanged (read-only over `&TranscriptFeature`) |
| `apply_buffer_local_hgnc_propagation` | 10942 | operate on `&mut [CowTranscript]` + base |
| `reset_buffer_local_hgnc_effective_values` | 10868 | `&mut [CowTranscript]` + base, compare-before-write |
| `reset_persisted_hgnc_effective_values_outside_start_region` | 10874 | `&mut [CowTranscript]` + base, compare-before-write |
| `build_buffer_local_transcripts` | 10760 | build `Vec<CowTranscript>` internally, return `Vec<TranscriptFeature>` (test path, signature unchanged) |
| `build_stateful_buffer_local_transcripts_cow` (new) | — | the real builder, returns `Vec<CowTranscript>` |
| `build_stateful_buffer_local_transcripts` | 10781 | thin owned wrapper over `_cow` (signature unchanged, keeps all existing tests green) |
| `apply_partition_transcript_overrides` | 10302 | `&mut [CowTranscript]` + base |
| `prepare_buffer_annotation_context` | 10356 | call `_cow`, pass base to overrides |
| `prepared_context_from_buffer` | 10446 | resolve `CowTranscript` → `&TranscriptFeature` |

`select_buffer_local_transcripts` (10730) is **left unchanged** (still returns `Vec<TranscriptFeature>`); it is used by the separate hydration path (`hydrate_worker_window`, 10497), which is out of scope (see "Out of scope").

---

## Task 1: Introduce the `CowTranscript` type with unit tests

**Files:**
- Modify: `datafusion/bio-function-vep/src/annotate_provider.rs` (add type + impl just above `struct BufferAnnotationContext` at ~9754)
- Test: same file, `#[cfg(test)] mod tests` (existing test module)

- [ ] **Step 1: Write the failing test**

Add to the existing test module in `annotate_provider.rs`. Reuse the existing test transcript constructor used by nearby tests (search the test module for how `TranscriptFeature` test instances are built — e.g. a `test_transcript(...)` helper or struct literal at ~14635/14654; use whichever the surrounding tests use). The test below assumes a helper `tf(id, hgnc)` that builds a `TranscriptFeature` with `transcript_id = id` and `gene_hgnc_id = gene_hgnc_id_native = hgnc`. If no such helper exists, add a minimal one in the test module mirroring the struct-literal at line 14654.

```rust
#[test]
fn cow_transcript_borrowed_resolves_to_base_owned_overrides() {
    let base = vec![tf("ENST1", Some("HGNC:1")), tf("ENST2", Some("HGNC:2"))];

    // Borrowed resolves to the base entry by index.
    let borrowed0 = CowTranscript::Borrowed(0);
    assert_eq!(borrowed0.as_ref(&base).transcript_id, "ENST1");
    assert_eq!(borrowed0.transcript_id(&base), "ENST1");

    // to_mut converts Borrowed -> Owned by cloning the base entry, then mutates
    // the owned copy without touching base.
    let mut bt = CowTranscript::Borrowed(1);
    bt.to_mut(&base).gene_hgnc_id = Some("HGNC:999".to_string());
    assert_eq!(bt.as_ref(&base).gene_hgnc_id.as_deref(), Some("HGNC:999"));
    assert!(matches!(bt, CowTranscript::Owned(_)));
    assert_eq!(base[1].gene_hgnc_id.as_deref(), Some("HGNC:2")); // base untouched

    // into_owned materializes a Borrowed into a clone equal to base.
    let owned = CowTranscript::Borrowed(0).into_owned(&base);
    assert_eq!(owned, base[0]);
}
```

- [ ] **Step 2: Run test to verify it fails**

Run: `cargo test -p datafusion-bio-function-vep cow_transcript_borrowed_resolves -- --nocapture`
Expected: FAIL to **compile** with "cannot find type `CowTranscript`".

- [ ] **Step 3: Write minimal implementation**

Insert immediately above `struct BufferAnnotationContext` (~line 9754):

```rust
/// Copy-on-write transcript: the same borrowed-or-owned shape as
/// `std::borrow::Cow`, but the `Borrowed` arm stores an *index* into
/// `SharedContigAnnotationContext::base_transcripts` instead of a reference, so
/// the value owns no lifetime and can be stored in `BufferAnnotationContext`.
/// `Borrowed` = no clone, geometry cache shared with the base instance;
/// `Owned` = a clone mutated by buffer-local HGNC propagation / partition
/// overrides. Replaces the former unconditional `Vec<TranscriptFeature>` clone
/// of every in-range transcript.
enum CowTranscript {
    Borrowed(usize),
    Owned(Box<TranscriptFeature>),
}

impl CowTranscript {
    /// Resolve to a borrowed transcript, indexing into `base` for `Borrowed`.
    /// Like `Cow::as_ref`, but resolves the index against `base` rather than
    /// dereferencing a stored reference.
    fn as_ref<'a>(&'a self, base: &'a [TranscriptFeature]) -> &'a TranscriptFeature {
        match self {
            CowTranscript::Borrowed(idx) => &base[*idx],
            CowTranscript::Owned(tx) => tx,
        }
    }

    fn transcript_id<'a>(&'a self, base: &'a [TranscriptFeature]) -> &'a str {
        self.as_ref(base).transcript_id.as_str()
    }

    /// Get a mutable transcript, converting `Borrowed` -> `Owned` (one clone
    /// from `base`) on first call. Idempotent for an already-`Owned` value.
    /// Mirrors `Cow::to_mut`.
    fn to_mut(&mut self, base: &[TranscriptFeature]) -> &mut TranscriptFeature {
        if let CowTranscript::Borrowed(idx) = self {
            *self = CowTranscript::Owned(Box::new(base[*idx].clone()));
        }
        match self {
            CowTranscript::Owned(tx) => tx,
            CowTranscript::Borrowed(_) => unreachable!("converted to Owned above"),
        }
    }

    /// Materialize into an owned `TranscriptFeature` (clones a `Borrowed`).
    /// Mirrors `Cow::into_owned`. Used only by the owned-returning test wrapper
    /// `build_*_local_transcripts`.
    fn into_owned(self, base: &[TranscriptFeature]) -> TranscriptFeature {
        match self {
            CowTranscript::Borrowed(idx) => base[idx].clone(),
            CowTranscript::Owned(tx) => *tx,
        }
    }
}
```

- [ ] **Step 4: Run test to verify it passes**

Run: `cargo test -p datafusion-bio-function-vep cow_transcript_borrowed_resolves -- --nocapture`
Expected: PASS.

- [ ] **Step 5: Commit**

```bash
git add datafusion/bio-function-vep/src/annotate_provider.rs
git commit -m "feat(vep): add CowTranscript copy-on-write transcript representation"
```

---

## Task 2: Convert the HGNC mutators to compare-before-write over `&mut [CowTranscript]`

These four functions currently take `&mut [TranscriptFeature]` and write unconditionally. Convert them to `&mut [CowTranscript]` + `base`, computing the desired value and only calling `to_mut` when it differs from the current value. `collect_hgnc_donors` is **unchanged** (it already reads `&TranscriptFeature`).

**Files:**
- Modify: `datafusion/bio-function-vep/src/annotate_provider.rs` (lines 10868-10987 region)
- Test: same file, test module

- [ ] **Step 1: Write the failing test**

```rust
#[test]
fn apply_hgnc_donors_only_clones_changed_transcripts() {
    // tf_full(id, gene_stable_id, gene_symbol, gene_symbol_source, hgnc_native)
    // builds a TranscriptFeature with gene_hgnc_id == hgnc_native.
    let base = vec![
        tf_full("ENST_A", Some("ENSG1"), Some("SYM1"), Some("HGNC"), Some("HGNC:1")),
        // ENST_B has no native HGNC but shares SYM1 -> should receive HGNC:1.
        tf_full("ENST_B", Some("ENSG1"), Some("SYM1"), Some("HGNC"), None),
    ];
    let mut buf = vec![CowTranscript::Borrowed(0), CowTranscript::Borrowed(1)];

    let donors: Vec<&TranscriptFeature> = base.iter().collect();
    let (hgnc_by_symbol, gene_fill) = collect_hgnc_donors(donors);
    apply_hgnc_donors(&mut buf, &base, &hgnc_by_symbol, &gene_fill);

    // ENST_A unchanged -> still Borrowed (no clone).
    assert!(matches!(buf[0], CowTranscript::Borrowed(_)));
    // ENST_B received the donated HGNC id -> Owned.
    assert!(matches!(buf[1], CowTranscript::Owned(_)));
    assert_eq!(buf[1].as_ref(&base).gene_hgnc_id.as_deref(), Some("HGNC:1"));
}
```

If `tf_full` does not exist, add it to the test module as a thin wrapper over the existing struct-literal pattern (line ~14654), setting `gene_stable_id`, `gene_symbol`, `gene_symbol_source`, `gene_hgnc_id_native`, and `gene_hgnc_id = gene_hgnc_id_native.clone()`.

- [ ] **Step 2: Run test to verify it fails**

Run: `cargo test -p datafusion-bio-function-vep apply_hgnc_donors_only_clones -- --nocapture`
Expected: FAIL to compile — `apply_hgnc_donors` still takes `&mut [TranscriptFeature]`.

- [ ] **Step 3: Write the implementation**

Replace `apply_hgnc_donors` (10947-10987) with this compare-before-write version. The two-pass structure and exact fill semantics are preserved; only the write is now conditional.

```rust
fn apply_hgnc_donors(
    transcripts: &mut [CowTranscript],
    base: &[TranscriptFeature],
    hgnc_by_symbol: &HashMap<String, String>,
    gene_fill_by_stable_id: &HashMap<String, HgncGeneFill>,
) {
    // Pass 1: effective HGNC id = native, else current, else symbol-donated.
    for bt in transcripts.iter_mut() {
        let cur = bt.as_ref(base);
        let mut desired = cur
            .gene_hgnc_id_native
            .clone()
            .or_else(|| cur.gene_hgnc_id.clone());
        if desired.is_none() {
            if let Some(hgnc_id) = cur
                .gene_symbol
                .as_deref()
                .and_then(|symbol| hgnc_by_symbol.get(symbol))
            {
                desired = Some(hgnc_id.clone());
            }
        }
        if desired != cur.gene_hgnc_id {
            bt.to_mut(base).gene_hgnc_id = desired;
        }
    }

    // Pass 2: gene-stable-id fill of symbol / source / hgnc when still None.
    for bt in transcripts.iter_mut() {
        let cur = bt.as_ref(base);
        let Some(gene_stable_id) = cur.gene_stable_id.as_deref() else {
            continue;
        };
        let Some(fill) = gene_fill_by_stable_id.get(gene_stable_id) else {
            continue;
        };
        let new_symbol = if cur.gene_symbol.is_none() {
            fill.gene_symbol.clone()
        } else {
            cur.gene_symbol.clone()
        };
        let new_source = if cur.gene_symbol_source.is_none() {
            fill.gene_symbol_source.clone()
        } else {
            cur.gene_symbol_source.clone()
        };
        let new_hgnc = if cur.gene_hgnc_id.is_none() {
            fill.gene_hgnc_id_native.clone()
        } else {
            cur.gene_hgnc_id.clone()
        };
        if new_symbol != cur.gene_symbol
            || new_source != cur.gene_symbol_source
            || new_hgnc != cur.gene_hgnc_id
        {
            let tx = bt.to_mut(base);
            tx.gene_symbol = new_symbol;
            tx.gene_symbol_source = new_source;
            tx.gene_hgnc_id = new_hgnc;
        }
    }
}
```

Replace `reset_buffer_local_hgnc_effective_values` (10868-10872):

```rust
fn reset_buffer_local_hgnc_effective_values(
    transcripts: &mut [CowTranscript],
    base: &[TranscriptFeature],
) {
    for bt in transcripts {
        let cur = bt.as_ref(base);
        let desired = cur.gene_hgnc_id_native.clone();
        if cur.gene_hgnc_id != desired {
            bt.to_mut(base).gene_hgnc_id = desired;
        }
    }
}
```

Replace `reset_persisted_hgnc_effective_values_outside_start_region` (10874-10891):

```rust
fn reset_persisted_hgnc_effective_values_outside_start_region(
    transcripts: &mut [CowTranscript],
    base: &[TranscriptFeature],
    transcript_cache_regions: &HashMap<String, Vec<TranscriptCacheRegion>>,
    active_regions: &HashSet<TranscriptCacheRegion>,
) {
    for bt in transcripts {
        let cur = bt.as_ref(base);
        let Some(regions) = transcript_cache_regions.get(&cur.transcript_id) else {
            continue;
        };
        let Some(start_region) = regions.first() else {
            continue;
        };
        let spans_multiple_regions = regions.len() > 1;
        if spans_multiple_regions && !active_regions.contains(start_region) {
            let desired = cur.gene_hgnc_id_native.clone();
            if cur.gene_hgnc_id != desired {
                bt.to_mut(base).gene_hgnc_id = desired;
            }
        }
    }
}
```

Replace `apply_buffer_local_hgnc_propagation` (10942-10945):

```rust
fn apply_buffer_local_hgnc_propagation(
    transcripts: &mut [CowTranscript],
    base: &[TranscriptFeature],
) {
    let (hgnc_by_symbol, gene_fill_by_stable_id) =
        collect_hgnc_donors(transcripts.iter().map(|bt| bt.as_ref(base)));
    apply_hgnc_donors(transcripts, base, &hgnc_by_symbol, &gene_fill_by_stable_id);
}
```

> Note: `collect_hgnc_donors` (10911-10940) is unchanged — its `impl IntoIterator<Item = &'a TranscriptFeature>` bound is satisfied by `transcripts.iter().map(|bt| bt.as_ref(base))`.

This step does NOT yet update the callers (`build_buffer_local_transcripts`, `build_stateful_buffer_local_transcripts`) — they are converted in Tasks 3-4. The crate will not fully compile until Task 4; run only the targeted unit test, which compiles the test module against the new signatures.

- [ ] **Step 4: Run the targeted test**

Run: `cargo test -p datafusion-bio-function-vep apply_hgnc_donors_only_clones -- --nocapture`
Expected: the crate fails to compile at the *callers* of the four changed functions (Tasks 3-4 fix them). If you want this task to compile standalone, temporarily proceed to Task 3+4 before running the suite. Otherwise verify by inspection that the test module + the four functions type-check, and commit; the green bar lands at Task 4 Step 4.

- [ ] **Step 5: Commit**

```bash
git add datafusion/bio-function-vep/src/annotate_provider.rs
git commit -m "refactor(vep): HGNC mutators compare-before-write over CowTranscript"
```

---

## Task 3: Convert `build_buffer_local_transcripts` (test path) to build via `CowTranscript`

This non-stateful builder has **no production caller** (only tests). Keep its public signature returning `Vec<TranscriptFeature>` so its tests are unchanged; internally build `Vec<CowTranscript>` over the passed slice and materialize at the end. This re-greens the test path and exercises the Task-2 functions.

**Files:**
- Modify: `datafusion/bio-function-vep/src/annotate_provider.rs:10760-10779`

- [ ] **Step 1: Replace the function body**

```rust
fn build_buffer_local_transcripts(
    transcripts: &[TranscriptFeature],
    chrom: &str,
    min_start: i64,
    max_end: i64,
    upstream_distance: i64,
    downstream_distance: i64,
) -> Vec<TranscriptFeature> {
    let mut buffer = select_buffer_local_indices(
        transcripts,
        chrom,
        min_start,
        max_end,
        upstream_distance,
        downstream_distance,
    );
    reset_buffer_local_hgnc_effective_values(&mut buffer, transcripts);
    apply_buffer_local_hgnc_propagation(&mut buffer, transcripts);
    buffer
        .into_iter()
        .map(|bt| bt.into_owned(transcripts))
        .collect()
}
```

- [ ] **Step 2: Add the index-selecting helper next to `select_buffer_local_transcripts` (~10730)**

`select_buffer_local_transcripts` stays as-is (owned, hydration path). Add a sibling that returns `Vec<CowTranscript>` with the **identical** filter predicate (keep the two predicates literally identical — a divergence here is a parity bug):

```rust
/// Same range filter as `select_buffer_local_transcripts`, but returns shared
/// indices (`CowTranscript::Borrowed`) into `transcripts` instead of clones.
fn select_buffer_local_indices(
    transcripts: &[TranscriptFeature],
    chrom: &str,
    min_start: i64,
    max_end: i64,
    upstream_distance: i64,
    downstream_distance: i64,
) -> Vec<CowTranscript> {
    let chrom_norm = chrom.strip_prefix("chr").unwrap_or(chrom);
    let up_down_size = upstream_distance.max(downstream_distance);
    let query_start = min_start.saturating_sub(up_down_size);
    let query_end = max_end.saturating_add(up_down_size);

    transcripts
        .iter()
        .enumerate()
        .filter(|(_, tx)| {
            let tx_chrom = tx.chrom.strip_prefix("chr").unwrap_or(&tx.chrom);
            tx_chrom == chrom_norm && tx.end >= query_start && tx.start <= query_end
        })
        .map(|(idx, _)| CowTranscript::Borrowed(idx))
        .collect()
}
```

- [ ] **Step 3: Run the build_buffer_local_transcripts tests**

Run: `cargo test -p datafusion-bio-function-vep build_buffer_local -- --nocapture`
Expected: the existing `build_buffer_local_transcripts` test at ~15539 (and any others matching) PASS — byte-identical owned output to before.

- [ ] **Step 4: Commit**

```bash
git add datafusion/bio-function-vep/src/annotate_provider.rs
git commit -m "refactor(vep): build_buffer_local_transcripts via CowTranscript (owned output unchanged)"
```

---

## Task 4: Add `build_stateful_buffer_local_transcripts_cow` + keep the owned wrapper

Split the stateful builder: the real logic moves to a `_cow` function returning `Vec<CowTranscript>`; the existing public name becomes a thin owned wrapper so **all ~16 existing test call sites stay unchanged** (15591, 15648, 15711, 15777, 15835, 15854, 15903, 15922, 15973, 15992, 16051, 16065, 16114, 16133, 16185, 16204).

**Files:**
- Modify: `datafusion/bio-function-vep/src/annotate_provider.rs:10781-10866`

- [ ] **Step 1: Replace the function with `_cow` + wrapper**

```rust
/// Owned-returning wrapper preserving the original signature for the unit
/// tests. Production code calls `_cow` directly (see prepare_buffer_annotation_context).
fn build_stateful_buffer_local_transcripts(
    transcripts: &[TranscriptFeature],
    transcript_cache_regions: &HashMap<String, Vec<TranscriptCacheRegion>>,
    persisted_transcripts: &mut PersistedCowTranscripts,
    buffer_batches: &[RecordBatch],
    chrom: &str,
    min_start: i64,
    max_end: i64,
    upstream_distance: i64,
    downstream_distance: i64,
) -> Result<Vec<TranscriptFeature>> {
    let buffer = build_stateful_buffer_local_transcripts_cow(
        transcripts,
        transcript_cache_regions,
        persisted_transcripts,
        buffer_batches,
        chrom,
        min_start,
        max_end,
        upstream_distance,
        downstream_distance,
    )?;
    Ok(buffer.into_iter().map(|bt| bt.into_owned(transcripts)).collect())
}

fn build_stateful_buffer_local_transcripts_cow(
    transcripts: &[TranscriptFeature],
    transcript_cache_regions: &HashMap<String, Vec<TranscriptCacheRegion>>,
    persisted_transcripts: &mut PersistedCowTranscripts,
    buffer_batches: &[RecordBatch],
    chrom: &str,
    min_start: i64,
    max_end: i64,
    upstream_distance: i64,
    downstream_distance: i64,
) -> Result<Vec<CowTranscript>> {
    let active_regions =
        collect_buffer_cache_regions(buffer_batches, upstream_distance, downstream_distance)?;
    prune_persisted_buffer_transcripts(persisted_transcripts, &active_regions);

    let mut buffer = select_buffer_local_indices(
        transcripts,
        chrom,
        min_start,
        max_end,
        upstream_distance,
        downstream_distance,
    );

    // Cross-buffer reuse: restore a persisted (already-mutated) copy. This is a
    // genuine divergence from base, so it becomes Owned.
    for bt in &mut buffer {
        let id = bt.transcript_id(transcripts);
        let persisted = active_regions.iter().find_map(|region| {
            persisted_transcripts
                .get(region)
                .and_then(|by_transcript| by_transcript.get(id))
        });
        if let Some(persisted) = persisted {
            *bt = CowTranscript::Owned(Box::new(persisted.clone()));
        }
    }

    reset_persisted_hgnc_effective_values_outside_start_region(
        &mut buffer,
        transcripts,
        transcript_cache_regions,
        &active_regions,
    );

    // Donors are the FULL active cache-region transcript set (partition-invariant).
    let region_donors: Vec<&TranscriptFeature> = transcripts
        .iter()
        .filter(|tx| {
            transcript_cache_regions
                .get(&tx.transcript_id)
                .is_some_and(|regions| regions.iter().any(|r| active_regions.contains(r)))
        })
        .collect();
    let (hgnc_by_symbol, gene_fill_by_stable_id) = collect_hgnc_donors(region_donors);
    apply_hgnc_donors(
        &mut buffer,
        transcripts,
        &hgnc_by_symbol,
        &gene_fill_by_stable_id,
    );

    // Re-persist only the transcripts that actually diverged from base (Owned).
    // Borrowed transcripts equal base, so restoring them next buffer is a no-op;
    // skipping them is behavior-preserving and avoids a clone into the persist map.
    for bt in &buffer {
        let CowTranscript::Owned(tx) = bt else {
            continue;
        };
        if let Some(regions) = transcript_cache_regions.get(&tx.transcript_id) {
            for region in regions.iter().filter(|region| active_regions.contains(*region)) {
                persisted_transcripts
                    .entry(region.clone())
                    .or_default()
                    .insert(tx.transcript_id.clone(), (**tx).clone());
            }
        }
    }

    Ok(buffer)
}
```

> Behavior note (covered by the Task-9 parity gate): re-persisting only `Owned` transcripts is the one subtle change. It is behavior-preserving because a `Borrowed` transcript is byte-equal to `base`, and the next buffer's `select_buffer_local_indices` produces that same `Borrowed` (effective `gene_hgnc_id == native`), so restoring a persisted-but-identical copy would yield the same value. If the parity gate ever fails and bisects to here, fall back to persisting all transcripts: replace the `let ... else { continue }` with materializing every `bt` via `bt.as_ref(transcripts)` and persist `bt.as_ref(transcripts).clone()`.

- [ ] **Step 2: Run the full stateful-builder test set**

Run: `cargo test -p datafusion-bio-function-vep build_stateful -- --nocapture`
Expected: all stateful-builder tests PASS (owned wrapper returns identical `Vec<TranscriptFeature>`).

- [ ] **Step 3: Run the whole crate test suite**

Run: `cargo test -p datafusion-bio-function-vep 2>&1 | tail -20`
Expected: all tests PASS (the production path still uses the owned wrapper at this point, so nothing observable changed yet).

- [ ] **Step 4: Commit**

```bash
git add datafusion/bio-function-vep/src/annotate_provider.rs
git commit -m "refactor(vep): build_stateful_buffer_local_transcripts_cow returns CowTranscript"
```

---

## Task 5: Switch `BufferAnnotationContext` + production build to `CowTranscript`

Now make the production path actually carry `Vec<CowTranscript>` end-to-end: change the struct field, call `_cow`, and update `apply_partition_transcript_overrides`, the `buffer_tx_ids` derivation, and `prepared_context_from_buffer`.

**Files:**
- Modify: `datafusion/bio-function-vep/src/annotate_provider.rs:9754-9759` (struct), `10302-10312` (overrides), `10356-10409` (prepare), `10446-10482` (prepared_context_from_buffer)

- [ ] **Step 1: Change the struct field type (9755)**

```rust
struct BufferAnnotationContext {
    transcripts: Vec<CowTranscript>,
    exon_indices: Vec<usize>,
    translation_indices: Vec<usize>,
    translation_overrides: HashMap<String, TranslationPartitionState>,
}
```

- [ ] **Step 2: Convert `apply_partition_transcript_overrides` (10302-10312)**

```rust
fn apply_partition_transcript_overrides(
    transcripts: &mut [CowTranscript],
    base: &[TranscriptFeature],
    overrides: &HashMap<String, TranscriptPartitionState>,
) {
    for bt in transcripts {
        let Some(override_state) = overrides.get(bt.transcript_id(base)) else {
            continue;
        };
        override_state.apply_to(bt.to_mut(base));
    }
}
```

- [ ] **Step 3: Update `prepare_buffer_annotation_context` (10368-10409)**

Change the builder call to `_cow`, pass `base` into the overrides + tx-id derivation:

```rust
    let mut transcripts = build_stateful_buffer_local_transcripts_cow(
        &shared.base_transcripts,
        &shared.transcript_cache_regions,
        &mut worker.persisted_buffer_transcripts,
        buffer_batches,
        chrom,
        min_start,
        max_end,
        config.upstream_distance,
        config.downstream_distance,
    )?;
    apply_partition_transcript_overrides(
        &mut transcripts,
        &shared.base_transcripts,
        &worker.transcript_overrides,
    );
    let buffer_tx_ids: HashSet<&str> = transcripts
        .iter()
        .map(|bt| bt.transcript_id(&shared.base_transcripts))
        .collect();
```

The rest of `prepare_buffer_annotation_context` (exon/translation indices, struct construction at 10404-10409) is unchanged — `transcripts` is now `Vec<CowTranscript>` and slots directly into the struct.

> Borrow check: `buffer_tx_ids` borrows `&shared.base_transcripts` immutably for the lifetime of `BufferAnnotationContext`'s construction; `transcripts` (the `Vec<CowTranscript>`) owns indices, not borrows, so moving it into the struct at 10404 does not conflict.

- [ ] **Step 4: Update `prepared_context_from_buffer` (10451)**

```rust
    let transcript_refs = buffer_context
        .transcripts
        .iter()
        .map(|bt| bt.as_ref(&shared.base_transcripts))
        .collect();
```

The `exon_refs` / `translation_refs` / `PreparedContext::new_from_refs` call below are unchanged. `transcript_refs` is `Vec<&'a TranscriptFeature>` exactly as before (lifetime `'a` from `shared`), so `PreparedContext<'a>` is byte-for-byte the same shape.

- [ ] **Step 5: Build + run the full suite**

Run: `cargo build -p datafusion-bio-function-vep 2>&1 | tail -20`
Expected: compiles clean.

Run: `cargo test -p datafusion-bio-function-vep 2>&1 | tail -20`
Expected: all tests PASS.

- [ ] **Step 6: Commit**

```bash
git add datafusion/bio-function-vep/src/annotate_provider.rs
git commit -m "feat(vep): production buffer annotation carries CowTranscript (CoW transcripts)"
```

---

## Task 6: Lint + format

**Files:**
- Modify: `datafusion/bio-function-vep/src/annotate_provider.rs` (as needed)

- [ ] **Step 1: Format**

Run: `cargo fmt`
Expected: no diff or only whitespace.

- [ ] **Step 2: Clippy (warnings as errors)**

Run: `cargo clippy -p datafusion-bio-function-vep --all-targets --features lance-cache -- -D warnings 2>&1 | tail -30`
Expected: clean. Likely fixups: a `match` that clippy wants as `if let`, or an unused `mut`. Apply the minimal idiomatic fix clippy names.

- [ ] **Step 3: Commit**

```bash
git add -A
git commit -m "style(vep): fmt + clippy for CowTranscript CoW path"
```

---

## Task 7: Build the release bench

**Files:** none (build only)

- [ ] **Step 1: Build the instrumented bench**

Run:
```bash
cargo build --release --features lance-cache --example bench_annotate_vcf 2>&1 | tail -5
```
Expected: `Finished release` with no errors. (The bench sets `mimalloc` as `#[global_allocator]` — required for the perf numbers.)

---

## Task 8: Reproduce the indexed input if missing

**Files:** none (data prep)

- [ ] **Step 1: Ensure the indexed 200k input exists**

Run:
```bash
test -f /tmp/chr1_200k.vcf.gz && echo "input present" || \
bcftools view /Users/mwiewior/workspace/data_vepyr/HG002_GRCh38_1_22_v4.2.1_benchmark.vcf.gz chr1 \
 | awk 'BEGIN{c=0} /^#/{print;next} {if(c<200000){print;c++}else{exit}}' \
 | bgzip > /tmp/chr1_200k.vcf.gz && tabix -p vcf -f /tmp/chr1_200k.vcf.gz
```
Expected: `input present` (it currently exists) or a freshly bgzipped+tabixed file.

---

## Task 9: MANDATORY byte-identical parity gate + scaling check

This is the acceptance gate. Output MUST be byte-identical at every thread count AND identical to the pre-change serial baseline.

**Files:** none (verification)

- [ ] **Step 1: Capture the pre-change serial baseline (from before this branch's work)**

If a known-good serial output exists, reuse it. Otherwise generate from the parent commit of this plan's work. With the current tree (post-change), produce the four thread outputs:

```bash
BIN=./target/release/examples/bench_annotate_vcf
CACHE=/Users/mwiewior/workspace/data_vepyr/115_GRCh38_merged
FASTA=/Users/mwiewior/workspace/data_vepyr/Homo_sapiens.GRCh38.dna.primary_assembly.fa
for T in 1 2 4 8; do
  $BIN --input /tmp/chr1_200k.vcf.gz --cache $CACHE \
    --output /tmp/c${T}.vcf --backend lance --everything \
    --reference-fasta $FASTA --threads $T --no-progress
done
```
Expected: four files `/tmp/c{1,2,4,8}.vcf`, each run exits 0.

- [ ] **Step 2: GATE — all thread counts byte-identical**

```bash
for T in 2 4 8; do
  if diff <(grep -v '^#' /tmp/c1.vcf) <(grep -v '^#' /tmp/c${T}.vcf) >/dev/null; then
    echo "t${T}: IDENTICAL to t1"
  else
    echo "t${T}: DIFFERS — GATE FAILED"; diff <(grep -v '^#' /tmp/c1.vcf) <(grep -v '^#' /tmp/c${T}.vcf) | head
  fi
done
```
Expected: `IDENTICAL to t1` for t2, t4, t8. Any DIFFERS = stop and debug (see superpowers:systematic-debugging); do NOT proceed.

- [ ] **Step 3: GATE — no serial regression vs the pre-change baseline**

Compare `/tmp/c1.vcf` against the saved pre-change serial output (e.g. `/tmp/c1_baseline.vcf` produced from the commit before Task 1). If no saved baseline exists, generate one by `git stash`-ing this branch's changes, building, and running threads=1 to `/tmp/c1_baseline.vcf`, then restoring.

```bash
diff <(grep -v '^#' /tmp/c1_baseline.vcf) <(grep -v '^#' /tmp/c1.vcf) >/dev/null \
  && echo "serial output unchanged" || { echo "SERIAL REGRESSION — GATE FAILED"; }
```
Expected: `serial output unchanged`.

- [ ] **Step 4: Scaling + engine-CPU-sum check**

Run the timed comparison with profiling, threads=8, and confirm the engine-CPU-sum inflation shrank toward the ~1/3 clock floor and wall-clock improved below the 4.94s baseline:

```bash
VEP_PROFILE=1 VEP_ENGINE_PROFILE=1 ./target/release/examples/bench_annotate_vcf \
  --input /tmp/chr1_200k.vcf.gz --cache $CACHE --output /tmp/c8_prof.vcf \
  --backend lance --everything --reference-fasta $FASTA --threads 8 --no-progress 2>&1 \
  | grep -iE 'engine|wall|tx_query|transcript_|total' | tail -40
```
Expected (success criteria from the problem doc):
- `t1` wall-clock no slower than the 16.0s baseline (ideally faster from less allocation).
- `t8` wall-clock below 4.94s.
- engine-CPU-sum inflation at t8 drops from +37% toward ~+12% (the clock-only floor).

Record the actual numbers in the commit message. If wall-clock did not improve despite parity holding, the win was smaller than the 76% clone-elimination predicted — capture the profile and note it; this does not fail the gate (parity is the gate), but it informs whether the follow-up (shared geometry cache) is worth pursuing.

- [ ] **Step 5: Commit the result**

```bash
git add -A
git commit -m "perf(vep): eliminate ~76% of per-buffer transcript clones via CoW CowTranscript

Parity: byte-identical at threads {1,2,4,8}; serial output unchanged.
Scaling: t8 <wall>s (was 4.94s); engine-CPU-sum inflation <X>% (was +37%)."
```

---

## Out of scope (follow-ups, not this plan)

- **Hydration path** (`hydrate_worker_window` at 10497 → `select_buffer_local_transcripts`): also clones the window transcript set every window for before/after snapshotting. Separate optimization; only active when a reference FASTA is provided and CDS hydration fires (`hydrate>0`). Left untouched here — `select_buffer_local_transcripts` keeps its owned signature for it.
- **Shared geometry cache for Owned transcripts:** `CowTranscript::Owned` clones reset the `GeometryCache OnceLock`, so the ~24% mutated transcripts still rebuild geometry. Decoupling `GeometryCache` into a shared side-cache (keyed by `transcript_id`) would share geometry for those too, but requires threading a cache parameter through the leaf functions that currently reach the cache via `&TranscriptFeature`. Only pursue if Task 9 Step 4 shows the residual geometry rebuild is still a measurable cost.
- **Output channel buffer RSS** (`PARALLEL_ANNOTATE_OUT_CAP` × N partitions): the other axis of the t8 RSS growth identified in the problem doc; orthogonal to transcript clones.
- **The ~1/3 boost-vs-all-core clock scaling:** physics, not fixable in software.

## Results (executed 2026-06-18, commits 7304470 / f731ba9 / e7ff0a0)

**Mandatory parity gate — PASS:**
- Gate A: `t1 ≡ t2 ≡ t4 ≡ t8` byte-identical data rows (200,000 each).
- Gate B: CoW `t1` data byte-identical to the pre-refactor (7304470) serial output.
- The only full-file `cmp` difference is non-deterministic `##FORMAT` header passthrough ordering (present even between my own t1/t8 runs; sorted `##FORMAT` sets identical) — unrelated to this change.

**Behavior-neutral on the unit suite:** 923 pass / 35 pre-existing failures, identical set before and after (the 35 fail at clean commit 570961d too — unrelated to transcript clones; not addressed by this plan).

**Clone elimination (measured pre-implementation, instrumented build path):** 75.9% of buffer transcripts stay `Borrowed` (clone eliminated); stable t1 75.87% / t8 75.94%.

**Performance (A/B, same machine, min of 2 reps):**

| metric | baseline (7304470) | CoW (e7ff0a0) |
|---|---|---|
| t1 wall | 17.09s | 16.90s |
| t8 wall | 6.06s | **5.61s (~7.4% faster)** |
| `evaluate_prepared` CPU-sum t1→t8 | 4.783 → 5.721s (+19.6%) | 4.746 → 5.610s (+18.2%) |
| `tx_query_total` CPU-sum t1→t8 | 4.695 → 5.624s (+19.8%) | 4.662 → 5.512s (+18.2%) |

CoW lowers the absolute engine-CPU-sum at both t1 and t8 and trims the t8 inflation. The win is real but modest (smaller than the problem doc's optimistic +37%→+12% framing): consistent with the structural analysis that the parallel partitions are mostly disjoint (limiting cross-worker cache dedup) and that a large part of the t8 inflation is boost-vs-all-core clock scaling (unfixable in software). The dominant recovered cost is within-worker geometry-`OnceLock` reuse across sliding buffers + reduced allocation.

**Follow-up signal:** only ~24% of buffer transcripts become `Owned` (and still rebuild geometry), so the heavier "shared geometry side-cache" follow-up would target a minority — pursue only if a geometry-heavy workload shows it dominating.

## Self-Review

- **Spec coverage:** problem doc §3 "eliminate the per-buffer transcript clone → immutable shared transcripts" → Tasks 1-5. §3 "engine's per-transcript HGNC read consults the side-table" → superseded by the CoW design (engine unchanged; HGNC stays on the resolved `&TranscriptFeature`), which the measurement justified as lower-risk. §4 verification gate → Task 9 (verbatim repro + diff). §5 success criteria → Task 9 Step 4.
- **Type consistency:** `CowTranscript::{get(base), to_mut(base), transcript_id(base), into_owned(base)}` used identically in Tasks 1-5. `build_stateful_buffer_local_transcripts` keeps its owned signature (wrapper); `_cow` is the new `Vec<CowTranscript>` variant. `select_buffer_local_indices` mirrors `select_buffer_local_transcripts`'s predicate exactly.
- **Placeholder scan:** none — every step shows the full function body it introduces.
- **Known divergence flagged:** re-persisting only `Owned` transcripts (Task 4) — documented, with a concrete fallback if the parity gate bisects to it.
