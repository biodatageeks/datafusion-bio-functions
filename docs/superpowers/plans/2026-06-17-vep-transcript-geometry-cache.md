# VEP Transcript cDNA-Geometry Cache Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Compute each transcript's sorted cDNA coordinate table + coding bounds **once per
transcript** (cached on `PreparedContext`) instead of rebuilding it per (variant × transcript),
eliminating redundancy **C** (and free-riders **#2** and half of **D**) from
`docs/superpowers/specs/2026-06-17-vep-everything-redundancy-analysis.md`.

**Architecture:** `PreparedContext` already owns one `&TranscriptFeature` per `tx_idx` and the exon
slices grouped by transcript. We add a parallel `Vec<TranscriptGeometry>` indexed by `tx_idx`, built
**eagerly** when the context is constructed (once per contig, reused across every batch / variant /
ALT). The hot path fetches `ctx.tx_geometry(tx_idx)` and threads `Option<&TranscriptGeometry>` into
the small set of leaf builders (`transcript_cdna_coords`, `coding_cdna_bounds`,
`genomic_to_cdna_index_for_transcript`, `raw_cdna_position_from_genomic`,
`unshifted_cdna_bounds_for_hgvs_shift`). When the option is `Some`, they read the cached table;
when `None` (every existing direct-call unit test) they compute exactly as today. This makes the
change **behavior-preserving by construction** — the e2e parity gate (0 mismatches) is the proof.

**Tech Stack:** Rust 2024, single crate `datafusion/bio-function-vep`. No new deps. The cached value
type `TranscriptCdnaCoord` (`transcript_consequence.rs:7301-7307`) is already `Copy`.

**Why eager, not lazy:** for chr-scale annotation nearly every transcript is overlapped by ≥1
variant, so eager ≈ VEP's lazy `{_variation_effect_feature_cache}`. Eager needs no interior
mutability / locks and is trivially `Sync` across the per-contig partition threads. (Lazy
`Vec<OnceLock<…>>` is a documented follow-up if we ever annotate sparse variant sets against a full
transcriptome — see "Out of scope".)

**Out of scope (separate follow-ups):**
- The `add_intron_splice_terms` exon re-sort (`transcript_consequence.rs:3131-3132`, inside
  `splice_terms` ~1.65 s) — entangled with exon-**number** ordering; do it after this lands.
- Redundancy **E** (AF/colocated per-row string clones) and **B** (HGVSc `within_feature` pre-gate).
- Converting `use_cdna_mapper_for_general_coords` to always-on by pre-populating
  `cdna_mapper_segments` for all transcripts — considered and rejected here (changes RefSeq gating
  semantics; higher correctness risk than threading an additive optional param).

---

## File Structure

| File | Responsibility | Change |
|---|---|---|
| `datafusion/bio-function-vep/src/transcript_consequence.rs` | engine, `PreparedContext`, cDNA-geometry leaf builders | Add `TranscriptGeometry` + `build_transcript_geometry`; add `tx_geometry` field + accessor to `PreparedContext`; add optional `geometry` param to leaf builders; thread from engine hot loop |
| `datafusion/bio-function-vep/src/hgvs.rs` | HGVSc formatting (`format_hgvsc`, `coding_cdna_bounds`) | Add optional `geometry` param to `coding_cdna_bounds`, `format_hgvsc`/`format_hgvsc_profiled`/`format_hgvsc_inner`; use cached coding bounds |

No new files. All changes are additive optional parameters + one new struct, keeping every existing
direct caller compiling via `None`.

---

## Invariants (must hold after every task)

1. `cargo build -p datafusion-bio-function-vep` is green.
2. `cargo test -p datafusion-bio-function-vep` is green (existing unit tests unchanged in behavior).
3. The geometry table is a **pure function of `(tx, tx_exons)`** — variant- and allele-invariant.
   Caching must not change any output byte. Proven by Task 5's parity gate.

---

### Task 1: `TranscriptGeometry` value + builder

**Files:**
- Modify: `datafusion/bio-function-vep/src/transcript_consequence.rs` (add after `transcript_cdna_coords`, i.e. after line `7561`)
- Test: same file, `#[cfg(test)] mod tests` (append near the existing `genomic_to_cdna_index_*` tests ~`12020`)

- [ ] **Step 1: Write the failing test**

Append to the test module in `transcript_consequence.rs`:

```rust
#[test]
fn transcript_geometry_matches_direct_builders() {
    // Two-exon + strand transcript; reuse the same fixture style as
    // `genomic_to_cdna_index_multi_exon`.
    let tx = test_transcript_plus_strand_two_exon();
    let exons = test_exons_plus_strand_two_exon();
    let exon_refs: Vec<&ExonFeature> = exons.iter().collect();

    let geom = build_transcript_geometry(&tx, &exon_refs);

    assert_eq!(geom.cdna_coords, transcript_cdna_coords(&tx, &exon_refs));
    assert_eq!(geom.coding_bounds, crate::hgvs::coding_cdna_bounds(&tx, &exon_refs, None));
}
```

> If `test_transcript_plus_strand_two_exon` / `test_exons_plus_strand_two_exon` helpers don't yet
> exist, reuse the inline fixture already built in `genomic_to_cdna_index_multi_exon`
> (`transcript_consequence.rs:12006`) by copying its `TranscriptFeature` + `Vec<ExonFeature>`
> construction into the test body instead of calling helpers.

- [ ] **Step 2: Run test to verify it fails**

Run: `cargo test -p datafusion-bio-function-vep transcript_geometry_matches_direct_builders`
Expected: FAIL — `cannot find function 'build_transcript_geometry'` and
`coding_cdna_bounds` takes 2 args not 3 (the 3-arg form arrives in Task 3; for Task 1 call the
current 2-arg form `crate::hgvs::coding_cdna_bounds(&tx, &exon_refs)` — update the test in Task 3).

- [ ] **Step 3: Write the struct + builder**

Insert after line `7561` (right after `transcript_cdna_coords`'s closing brace):

```rust
/// Per-transcript cDNA geometry, computed once and reused across every variant
/// and ALT allele overlapping the transcript. Mirrors Ensembl's
/// `Transcript->{_variation_effect_feature_cache}` (built by the first
/// overlapping variant, reused by all the rest). Stored on `PreparedContext`,
/// keyed by `tx_idx`.
#[derive(Debug, Clone, Default, PartialEq, Eq)]
pub(crate) struct TranscriptGeometry {
    /// Genomic-sorted exon→cDNA segment table (the expensive build). `None`
    /// when the transcript has no usable exon set in this context slice.
    pub cdna_coords: Option<Vec<TranscriptCdnaCoord>>,
    /// (start_codon, stop_codon) in cDNA space, or `None` for non-coding.
    pub coding_bounds: Option<(usize, usize)>,
}

/// Build the per-transcript geometry once. Pure function of `(tx, tx_exons)`.
pub(crate) fn build_transcript_geometry(
    tx: &TranscriptFeature,
    tx_exons: &[&ExonFeature],
) -> TranscriptGeometry {
    TranscriptGeometry {
        cdna_coords: transcript_cdna_coords(tx, tx_exons),
        coding_bounds: crate::hgvs::coding_cdna_bounds(tx, tx_exons, None),
    }
}
```

> For Task 1 only, temporarily call `crate::hgvs::coding_cdna_bounds(tx, tx_exons)` (2-arg). Task 3
> adds the 3rd `geometry` arg and you update this call to pass `None`.

- [ ] **Step 4: Run test to verify it passes**

Run: `cargo test -p datafusion-bio-function-vep transcript_geometry_matches_direct_builders`
Expected: PASS

- [ ] **Step 5: Commit**

```bash
git add datafusion/bio-function-vep/src/transcript_consequence.rs
git commit -m "feat(vep): add TranscriptGeometry value + builder (redundancy C, phase 1)"
```

---

### Task 2: Eager `tx_geometry` cache on `PreparedContext` + accessor

**Files:**
- Modify: `datafusion/bio-function-vep/src/transcript_consequence.rs:811-824` (struct), `:1000-1011` (constructor return), `:900-1012` (`impl`)
- Test: same file, test module

- [ ] **Step 1: Write the failing test**

```rust
#[test]
fn prepared_context_caches_geometry_per_transcript() {
    let tx = test_transcript_plus_strand_two_exon();
    let exons = test_exons_plus_strand_two_exon();
    let ctx = PreparedContext::new(
        std::slice::from_ref(&tx),
        &exons,
        &[], &[], &[], &[], &[],
    );
    // transcripts are sorted by id; find our tx's index.
    let idx = ctx
        .transcripts
        .iter()
        .position(|t| t.transcript_id == tx.transcript_id)
        .unwrap();
    let exon_refs: Vec<&ExonFeature> = exons.iter().collect();
    let expected = build_transcript_geometry(&tx, &exon_refs);
    assert_eq!(ctx.tx_geometry(idx), &expected);
}
```

- [ ] **Step 2: Run test to verify it fails**

Run: `cargo test -p datafusion-bio-function-vep prepared_context_caches_geometry_per_transcript`
Expected: FAIL — no method `tx_geometry`, no field `tx_geometry`.

- [ ] **Step 3: Add the field, eager build, and accessor**

Add the field to the struct at `transcript_consequence.rs:817` (right after `transcripts`):

```rust
    /// Per-transcript cDNA geometry, aligned 1:1 with `transcripts` by index.
    /// Built once at context construction; read by `tx_geometry(tx_idx)`.
    tx_geometry: Vec<TranscriptGeometry>,
```

In `new_from_refs`, after the transcripts are sorted (`transcript_consequence.rs:953`, the
`tx_refs.sort_unstable_by(...)` line) and `exons_by_tx_genomic` is built (`:941-943`), insert:

```rust
        // Precompute per-transcript cDNA geometry once (redundancy C). Genomic
        // exon order is used; the builders re-sort internally so either exon set
        // is acceptable. Empty slice → None coords (consumers fall back).
        let tx_geometry: Vec<TranscriptGeometry> = tx_refs
            .iter()
            .map(|tx| {
                let tx_exons = exons_by_tx_genomic
                    .get(tx.transcript_id.as_str())
                    .map(Vec::as_slice)
                    .unwrap_or(&[]);
                build_transcript_geometry(tx, tx_exons)
            })
            .collect();
```

Add `tx_geometry,` to the `Self { … }` literal at `transcript_consequence.rs:1004` (next to
`transcripts: tx_refs,`).

Add the accessor inside `impl<'a> PreparedContext<'a>` (after `new_from_refs`, before line `1012`):

```rust
    /// Cached per-transcript cDNA geometry for `tx_idx` (index into `transcripts`).
    #[inline]
    pub(crate) fn tx_geometry(&self, tx_idx: usize) -> &TranscriptGeometry {
        &self.tx_geometry[tx_idx]
    }
```

- [ ] **Step 4: Run test to verify it passes**

Run: `cargo test -p datafusion-bio-function-vep prepared_context_caches_geometry_per_transcript`
Expected: PASS. Also run `cargo build -p datafusion-bio-function-vep` — expected: green (the new
field is private; all `PreparedContext { … }` literals are inside this constructor).

- [ ] **Step 5: Commit**

```bash
git add datafusion/bio-function-vep/src/transcript_consequence.rs
git commit -m "feat(vep): build per-transcript geometry once on PreparedContext"
```

---

### Task 3: Optional `geometry` param on the leaf builders (cache-or-compute)

Goal: each leaf builder gains a trailing `geometry: Option<&TranscriptGeometry>`. `Some` → use the
cached table; `None` → identical to today. Bottom-up so each compiles before the next.

**Files:**
- Modify: `datafusion/bio-function-vep/src/hgvs.rs:1446-1458` (`coding_cdna_bounds`) and its 3 call
  sites (`:1184`, `:1514`, `:1618`)
- Modify: `datafusion/bio-function-vep/src/transcript_consequence.rs:7487-7500`
  (`genomic_to_cdna_index_for_transcript`), `:7702-7736` (`raw_cdna_position_from_genomic`),
  `:7748-7804` (`unshifted_cdna_bounds_for_hgvs_shift`)
- Test: `transcript_consequence.rs` test module

- [ ] **Step 1: Write the failing equality tests**

```rust
#[test]
fn coding_cdna_bounds_cache_equals_compute() {
    let tx = test_transcript_plus_strand_two_exon();
    let exons = test_exons_plus_strand_two_exon();
    let exon_refs: Vec<&ExonFeature> = exons.iter().collect();
    let geom = build_transcript_geometry(&tx, &exon_refs);
    let computed = crate::hgvs::coding_cdna_bounds(&tx, &exon_refs, None);
    let cached = crate::hgvs::coding_cdna_bounds(&tx, &exon_refs, Some(&geom));
    assert_eq!(computed, cached);
}

#[test]
fn raw_cdna_position_cache_equals_compute() {
    let tx = test_transcript_plus_strand_two_exon();
    let exons = test_exons_plus_strand_two_exon();
    let exon_refs: Vec<&ExonFeature> = exons.iter().collect();
    let geom = build_transcript_geometry(&tx, &exon_refs);
    for pos in [exons[0].start, exons[0].end, exons[1].start, exons[1].end] {
        assert_eq!(
            raw_cdna_position_from_genomic(&tx, &exon_refs, pos, None),
            raw_cdna_position_from_genomic(&tx, &exon_refs, pos, Some(&geom)),
        );
    }
}
```

- [ ] **Step 2: Run to verify they fail**

Run: `cargo test -p datafusion-bio-function-vep cache_equals_compute`
Expected: FAIL — arity mismatch (`coding_cdna_bounds`/`raw_cdna_position_from_genomic` take fewer
args).

- [ ] **Step 3a: `coding_cdna_bounds` (hgvs.rs:1446)**

Replace the signature + add the early cache hit:

```rust
fn coding_cdna_bounds(
    tx: &TranscriptFeature,
    tx_exons: &[&ExonFeature],
    geometry: Option<&crate::transcript_consequence::TranscriptGeometry>,
) -> Option<(usize, usize)> {
    if let Some(bounds) = geometry.and_then(|g| g.coding_bounds) {
        return Some(bounds);
    }
    if let (Some(start_codon), Some(stop_codon)) = (tx.cdna_coding_start, tx.cdna_coding_end) {
        return Some((start_codon, stop_codon));
    }
    // ... unchanged body from :1450 onward ...
```

Update its 3 call sites:
- `hgvs.rs:1184` `|| coding_cdna_bounds(tx, tx_exons).is_some()` → `|| coding_cdna_bounds(tx, tx_exons, None).is_some()`
- `hgvs.rs:1514` `coding_cdna_bounds(tx, tx_exons)` → `coding_cdna_bounds(tx, tx_exons, geometry)` (see Step 3d — `shift_to_hgvs_coding_coordinates` gains the param)
- `hgvs.rs:1618` `|| coding_cdna_bounds(tx, tx_exons).is_some()` → `|| coding_cdna_bounds(tx, tx_exons, None).is_some()`
- `transcript_consequence.rs:7561` (Task 1 builder) → `crate::hgvs::coding_cdna_bounds(tx, tx_exons, None)`
- Task 1 test's direct call → add `, None`.

> Make `coding_cdna_bounds` and `shift_to_hgvs_coding_coordinates` visible to the test/engine:
> they are private in `hgvs.rs`; mark `pub(crate)` if a cross-module caller needs them (the engine
> calls `format_hgvsc`, not these, so `pub(crate)` is only needed for the Task-1/Task-3 tests that
> reference `crate::hgvs::coding_cdna_bounds`). Add `pub(crate)` to `coding_cdna_bounds`.

- [ ] **Step 3b: `genomic_to_cdna_index_for_transcript` (transcript_consequence.rs:7487)**

```rust
pub(crate) fn genomic_to_cdna_index_for_transcript(
    tx: &TranscriptFeature,
    tx_exons: &[&ExonFeature],
    pos: i64,
    geometry: Option<&TranscriptGeometry>,
) -> Option<usize> {
    if use_cdna_mapper_for_general_coords(tx) {
        return tx
            .cdna_mapper_segments
            .iter()
            .find_map(|segment| mapper_segment_cdna_index(segment, pos));
    }
    if let Some(coords) = geometry.and_then(|g| g.cdna_coords.as_deref()) {
        // Binary-searchable cached table replaces exon_segments rebuild + linear scan.
        return cdna_index_from_sorted_coords(tx.strand, pos, coords);
    }
    genomic_to_cdna_index(tx_exons, tx.strand, pos)
}
```

Add a helper next to `genomic_to_cdna_index` (`:7283`) that maps a genomic pos through an already
genomic-sorted `&[TranscriptCdnaCoord]` (binary search on `start`):

```rust
/// Map a genomic position to its 1-based cDNA index using a prebuilt,
/// genomic-start-sorted coord table. Equivalent to `genomic_to_cdna_index`
/// but O(log n) and zero-alloc.
fn cdna_index_from_sorted_coords(
    strand: i8,
    pos: i64,
    coords: &[TranscriptCdnaCoord],
) -> Option<usize> {
    let i = coords.partition_point(|c| c.end < pos);
    let seg = coords.get(i)?;
    if pos < seg.start || pos > seg.end {
        return None;
    }
    let local = if strand >= 0 {
        usize::try_from(pos - seg.start).ok()?
    } else {
        usize::try_from(seg.end - pos).ok()?
    };
    Some(seg.cdna_start + local)
}
```

> Add a unit test `cdna_index_from_sorted_coords_matches_linear` asserting it equals
> `genomic_to_cdna_index(exon_refs, strand, pos)` for every base across the two-exon fixture, both
> strands. (Write it in this step, run, confirm PASS — this is the load-bearing equivalence.)

Update `genomic_to_cdna_index_for_transcript`'s callers to pass geometry where available, else
`None`: `coding_cdna_bounds` body (`hgvs.rs:1455-1456`) passes `geometry`; the `hgvsp` helpers
(`genomic_to_cdna_index_for_hgvsp`, `:7336`; `genomic_to_cds_index_for_hgvsp`, `:7339`) pass `None`
for now (out of scope); `unshifted_cdna_bounds_for_hgvs_shift` (Step 3c) passes its `geometry`.

- [ ] **Step 3c: `raw_cdna_position_from_genomic` (:7702) & `unshifted_cdna_bounds_for_hgvs_shift` (:7748)**

Add `geometry: Option<&TranscriptGeometry>` as the trailing param to both. In
`raw_cdna_position_from_genomic`, replace the fallthrough at `:7734-7735`:

```rust
    let coords = match geometry.and_then(|g| g.cdna_coords.as_deref()) {
        Some(coords) => return raw_cdna_position_from_sorted_coords(
            tx.strand, genomic_pos, coords.iter().copied(),
        ),
        None => transcript_cdna_coords(tx, tx_exons)?,
    };
    raw_cdna_position_from_sorted_coords(tx.strand, genomic_pos, coords)
```

In `unshifted_cdna_bounds_for_hgvs_shift`, replace `:7756` `let coords = transcript_cdna_coords(tx, tx_exons)?;`
with a cache-or-compute that yields a `Vec<TranscriptCdnaCoord>` (the `(None,None)` branch at
`:7780-7793` iterates `coords`, so materialize once):

```rust
    let coords = match geometry.and_then(|g| g.cdna_coords.as_deref()) {
        Some(coords) => coords.to_vec(),
        None => transcript_cdna_coords(tx, tx_exons)?,
    };
```

and pass `geometry` into the four `genomic_to_cdna_index_for_transcript(tx, tx_exons, …)` calls in
its body (`:7760, :7761, :7798, :7799`).

- [ ] **Step 3d: ripple through the two `raw_cdna_position_from_genomic` callers inside hgvs.rs**

`hgvs_cdna_position_from_genomic` (`hgvs.rs:1460`) calls `raw_cdna_position_from_genomic` at
`:1465` and the `native_refseq_*` helpers. Add `geometry: Option<&TranscriptGeometry>` to
`hgvs_cdna_position_from_genomic` and `shift_to_hgvs_coding_coordinates` (`:1508`), threading it to
`raw_cdna_position_from_genomic(:1465, …, geometry)` and `coding_cdna_bounds(:1514, …, geometry)`.
Pass `None` to the `native_refseq_*` helpers (out of scope; they keep computing).

- [ ] **Step 4: Run the equality tests**

Run: `cargo test -p datafusion-bio-function-vep cache_equals_compute cdna_index_from_sorted_coords`
Expected: PASS. Then `cargo test -p datafusion-bio-function-vep` — expected: all green (all
production call sites currently pass `None`, so behavior is byte-identical).

- [ ] **Step 5: Commit**

```bash
git add datafusion/bio-function-vep/src/transcript_consequence.rs datafusion/bio-function-vep/src/hgvs.rs
git commit -m "feat(vep): thread optional cached geometry through cDNA leaf builders"
```

---

### Task 4: Wire the cache into the engine hot path

Goal: the engine fetches `ctx.tx_geometry(tx_idx)` once per (variant × transcript) and passes
`Some(&geometry)` into the cDNA-position + HGVSc path. After this, the per-call rebuild is gone for
the common (non-mapper) transcripts; the cache is hit instead.

**Files:**
- Modify: `datafusion/bio-function-vep/src/transcript_consequence.rs:1157-1220` (engine block),
  `:7820` (`compute_cdna_position` gains `geometry`)
- Modify: `datafusion/bio-function-vep/src/hgvs.rs:171-263` (`format_hgvsc`, `format_hgvsc_profiled`,
  `format_hgvsc_inner` gain `geometry`)
- Test: an end-to-end-ish engine test asserting cached == uncached CSQ for a fixture batch

- [ ] **Step 1: Write the failing test**

Add to the test module — reuse the existing engine-level fixture pattern from
`compute_cdna_position_uses_transcript_mapper_segments` (`:12691`) but drive the public engine:

```rust
#[test]
fn engine_output_identical_with_geometry_cache() {
    // Build a small transcripts+exons fixture and a handful of variants that
    // hit exonic, intronic, and boundary positions.
    let (transcripts, exons, translations) = small_engine_fixture();
    let variants = small_engine_variants();
    let engine = TranscriptConsequenceEngine::new(5000, 5000);

    let ctx = PreparedContext::new(
        &transcripts, &exons, &translations, &[], &[], &[], &[],
    );
    // Golden = current behavior captured before wiring (snapshot the Vec<TranscriptConsequence>
    // hgvsc/cdna_position fields). After wiring, the engine must reproduce them exactly.
    for v in &variants {
        let got = engine.evaluate_variant_prepared(v, &ctx);
        insta_like_assert_stable(v, &got); // compare hgvsc/cdna_position to the recorded golden
    }
}
```

> Practical form: capture the golden by running this test once **before** Step 3 wiring (it will
> read from `None`-path) and store the expected `(hgvsc, cdna_position)` tuples inline in the test.
> Then apply Step 3 and confirm the same tuples. If a separate snapshot helper doesn't exist, encode
> the expected values as a literal `Vec<(&str, Option<&str>)>` in the test.

- [ ] **Step 2: Run to verify it fails / captures golden**

Run: `cargo test -p datafusion-bio-function-vep engine_output_identical_with_geometry_cache`
Expected (pre-wiring): PASS against the `None` path — this **records the golden**. Commit the golden
literals so Step 4 re-runs against them.

- [ ] **Step 3: Thread geometry from the engine**

In the per-transcript closure, after `let tx = ctx.transcripts[tx_idx];` (`:1158`) add:

```rust
                let tx_geometry = ctx.tx_geometry(tx_idx);
```

At `:1220` change `compute_cdna_position(variant, tx, tx_exons)` →
`compute_cdna_position(variant, tx, tx_exons, Some(tx_geometry))`, and add the trailing
`geometry: Option<&TranscriptGeometry>` param to `compute_cdna_position` (`:7820`), threading it to
its internal `raw_cdna_position_from_genomic` / `genomic_to_cdna_index_for_transcript` calls (pass
`geometry` through; the mapper fast-paths ignore it).

At `:1368` / `:1355` change `format_hgvsc(tx, tx_exons_for_hgvsc, …)` /
`format_hgvsc_profiled(...)` to pass `Some(tx_geometry)` as a new trailing arg, and add
`geometry: Option<&TranscriptGeometry>` to `format_hgvsc` (`:171`), `format_hgvsc_profiled` (`:197`),
and `format_hgvsc_inner` (`:224`), threading it into `hgvs_cdna_position_from_genomic`,
`shift_to_hgvs_coding_coordinates`, and `coding_cdna_bounds` inside the inner formatter.

> Note: `tx_exons_for_hgvsc` is the **genomic**-sorted set; the cached `cdna_coords` was built from
> the genomic set too, so they are consistent. Leave the `exon_str`/`intron_str` calls (`:1218-1219`)
> on `tx_exons` (exon-number set) untouched — they need exon numbering, not the coord table.

- [ ] **Step 4: Run the golden test + full suite**

Run: `cargo test -p datafusion-bio-function-vep engine_output_identical_with_geometry_cache`
Expected: PASS (same golden tuples, now served from cache).
Run: `cargo test -p datafusion-bio-function-vep` — expected: all green.

- [ ] **Step 5: Commit**

```bash
git add datafusion/bio-function-vep/src/transcript_consequence.rs datafusion/bio-function-vep/src/hgvs.rs
git commit -m "perf(vep): serve cDNA geometry from per-transcript cache in engine hot path"
```

---

### Task 5: Verification gate — parity + profile delta (NO success claim without evidence)

**Files:** none (validation only).

- [ ] **Step 1: Full unit + clippy**

Run: `cargo test -p datafusion-bio-function-vep && cargo clippy -p datafusion-bio-function-vep -- -D warnings`
Expected: all green, no warnings.

- [ ] **Step 2: e2e parity gate WITHOUT `--skip-compare`**

This change touches HGVSc/cDNA-position correctness; the gate must stay **0 mismatches** vs the
Ensembl VEP reference (the same gate cited in the redundancy analysis: 4,737,090 rows, all 86 CSQ
fields). Run the project's e2e parity harness on chr1 **without** `--skip-compare`.
Expected: `0 mismatches`. If any field diverges, the cache is being read where the compute path
would have branched differently (likely a mapper/genomic-order mismatch) — STOP and diff, do not
proceed.

- [ ] **Step 3: Re-profile with the `vep-perf-profiling` skill**

Invoke the `vep-perf-profiling` skill (single-thread, all levels) on the chr1 workload. Read
`evaluate_prepared` → `transcript_hgvsc` → `coord_map` and the `fallback` line.
Expected direction: `coord_map` and `fallback` drop (cached table replaces rebuild + linear scan);
total engine drops toward the analysis's conservative **~2–4 s**. Record the before/after numbers in
the commit message / handoff. Do **not** claim the win until this shows the drop — trust the
profile, not the prediction.

- [ ] **Step 4: Record outcome**

Append measured deltas to `docs/superpowers/specs/2026-06-17-vep-everything-redundancy-analysis.md`
under redundancy C ("Measured after geometry cache: coord_map X→Y s, fallback A→B s, engine N→M s").

- [ ] **Step 5: Commit**

```bash
git add docs/superpowers/specs/2026-06-17-vep-everything-redundancy-analysis.md
git commit -m "docs(vep): record measured impact of per-transcript geometry cache"
```

---

## Self-Review

- **Spec coverage:** redundancy C (per-(variant×transcript) coord-table + coding-bounds rebuild) →
  Tasks 1–4; #2 (per-ALT) and half of D (`coding_cdna_bounds` recompute in fallback) are subsumed
  because the cache is keyed per transcript and reused across rows/ALTs and across the
  fast/fallback HGVSc paths. The `add_intron_splice_terms` exon re-sort and redundancies B/E are
  explicitly **out of scope** (noted up top), not silently dropped.
- **Placeholders:** none — every code step shows real signatures/bodies; the two "golden capture"
  steps describe the exact mechanism (record literals from the `None` path, re-assert after wiring).
- **Type consistency:** `TranscriptGeometry { cdna_coords: Option<Vec<TranscriptCdnaCoord>>,
  coding_bounds: Option<(usize, usize)> }` and `geometry: Option<&TranscriptGeometry>` are used
  identically in every task; `coding_cdna_bounds` is 3-arg everywhere after Task 3;
  `cdna_index_from_sorted_coords` returns `Option<usize>` matching `genomic_to_cdna_index`.
- **Risk:** additive optional params + `None` at every untouched call site ⇒ behavior-preserving
  until the engine opts in (Task 4); the parity gate (Task 5) is the hard backstop.
