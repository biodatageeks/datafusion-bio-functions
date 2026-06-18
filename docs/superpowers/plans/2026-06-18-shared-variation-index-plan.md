# Shared Variation-Lookup Index — Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Build each contig's lance variation lookup (`SinglePathLanceVariationLookup` = dataset + `PositionRowIdIndex`) once and share it across all per-partition lookup streams via an `Arc`, eliminating the ~2.6 GB/partition fixed RSS cost (t8 ~23 GB → ~6–7 GB), byte-identical.

**Architecture:** The variation lookup currently lives in per-partition `KvLookupStream` state (`lance_chroms`), built lazily on first probe — so N partitions build N copies of the 88M-entry index. We hoist it to a **shared `Arc<tokio::sync::OnceCell<Arc<SinglePathLanceVariationLookup>>>`** owned by the `KvLookupExec` and cloned (`Arc`) into every stream. The first stream to probe builds it via `OnceCell::get_or_try_init` (single-flight — concurrent fan-out callers await the one build, no race, one ~2.6 GB transient); all streams then share the one `Arc<lookup>`. The per-partition walk `cursor` (`lance_cursors`) stays per-stream. `resolve_and_take(&self, …, cursor)` and `take_rows` are unchanged. This faithfully implements the spec's "build once per contig, share among workers" (design spec `docs/superpowers/specs/2026-06-18-shared-variation-index-design.md`); we use a shared single-flight `OnceCell` rather than eager-build-in-`prepare_contig_context` because it reuses the existing path/projection build logic verbatim (lowest parity risk) while still guaranteeing exactly one build.

**Tech Stack:** Rust 2024, DataFusion 52.1, Arrow 57, `lance`/`lance-file`/`lance-index`, `tokio` (`tokio::sync::OnceCell`). Crate `datafusion-bio-function-vep`, feature `lance-cache`. Bench: `examples/bench_annotate_vcf.rs` (mimalloc global allocator).

**TDD adaptation:** Externally-observable behavior is unchanged (byte-identical), so the primary verification is the e2e byte-identical parity gate + the peak-RSS measurement, plus "unit-suite failure set unchanged vs baseline" (the branch has 35 pre-existing unrelated failures — see `transcript-clone-elimination.md`). No new unit assertions on new behavior; the change is a placement refactor whose correctness is proven by parity.

---

## File Structure

One file changed plus verification: `datafusion/bio-function-vep/src/kv_cache/cache_exec.rs`.

| Symbol | Current | Change |
|---|---|---|
| `KvLookupExec` (struct, `:76`) | per-backend fields | add `lance_lookup_cell: Arc<OnceCell<Arc<SinglePathLanceVariationLookup>>>` |
| `KvLookupExec::new_lance` (`:260`) | builds exec | init `lance_lookup_cell: Arc::new(OnceCell::new())` |
| `KvLookupExec::with_new_children` (lance arm, `:416`) | rebuilds exec | carry `lance_lookup_cell: Arc::clone(&self.lance_lookup_cell)` forward |
| `KvLookupExec::execute` (stream init, `:1903`) | `lance_chroms: HashMap::new()` | seed stream with `Arc::clone(&self.lance_lookup_cell)` |
| `KvLookupStream` (struct, `:497`) | `lance_chroms: HashMap<String, Option<…>>` | replace with `lance_lookup_cell: Arc<OnceCell<Arc<…>>>`; keep `lance_cursors` |
| `ensure_cold_parquet_lookup` (lance arm, `:2589`) | `lance_chroms.insert(open(...))` | `get_or_try_init` into the shared cell |
| resolve site (`:3642`) | `self.lance_chroms.get(&chrom)` | `self.lance_lookup_cell.get()` |

`SinglePathLanceVariationLookup` (`lance_cache/variation_runtime.rs`) is unchanged. `lance_cursors` stays per-stream.

---

## Task 1: Add the shared `OnceCell` to `KvLookupExec`

**Files:** Modify `datafusion/bio-function-vep/src/kv_cache/cache_exec.rs` (struct `:76`, `new_lance` `:260`, `with_new_children` lance arm `:416`).

- [ ] **Step 1: Confirm imports**

At the top of `cache_exec.rs`, ensure these are in scope (add if missing, behind the existing `#[cfg(feature = "lance-cache")]` import grouping where the other lance types are imported):

```rust
#[cfg(feature = "lance-cache")]
use tokio::sync::OnceCell;
// `SinglePathLanceVariationLookup` and `Arc` are already imported (used by lance_chroms today).
```

Run: `grep -n "use tokio::sync::OnceCell\|SinglePathLanceVariationLookup\|use std::sync::Arc" datafusion/bio-function-vep/src/kv_cache/cache_exec.rs | head`
Expected: `Arc` and `SinglePathLanceVariationLookup` already present; add the `OnceCell` use if not shown.

- [ ] **Step 2: Add the field to `KvLookupExec`**

In `struct KvLookupExec { … }` (`:76`), after the `target_partitions: usize,` field, add:

```rust
    /// Per-contig lance variation lookup (dataset + PositionRowIdIndex), built
    /// once and shared across all per-partition streams via the Arc. Single-flight
    /// via OnceCell so simultaneous fan-out probes build it exactly once.
    #[cfg(feature = "lance-cache")]
    lance_lookup_cell: Arc<OnceCell<Arc<SinglePathLanceVariationLookup>>>,
```

- [ ] **Step 3: Initialize it in `new_lance`**

In `KvLookupExec::new_lance` (`:260`), in the returned `Ok(Self { … })` literal, add:

```rust
            #[cfg(feature = "lance-cache")]
            lance_lookup_cell: Arc::new(OnceCell::new()),
```

(Other `new_*` constructors — `new`, `new_indexed_parquet` — do not have this field under `#[cfg(feature = "lance-cache")]`; if the struct field is unconditionally compiled, they must init it too. It is gated `#[cfg(feature = "lance-cache")]`, so only lance-cache builds see the field and only those constructors need it. If a non-lance constructor is compiled under `lance-cache`, add the same init line there — the compiler will name the exact constructor.)

- [ ] **Step 4: Carry it forward in `with_new_children`**

In `with_new_children`, in the `VariationLookupStorage::Lance { cache_root } => KvLookupExec::new_lance(…)?` arm (`:416`), the new exec is rebuilt via `new_lance`, which initializes a **fresh empty** cell — that would discard a build from a prior optimization pass. To preserve sharing across re-planning, after the `let mut exec = …` block (where `exec` is reassigned for all arms, near `:430`), add:

```rust
        #[cfg(feature = "lance-cache")]
        if let VariationLookupStorage::Lance { .. } = &self.variation_storage {
            exec.lance_lookup_cell = Arc::clone(&self.lance_lookup_cell);
        }
```

Place this alongside the existing post-construction `exec = exec.with_*(…)` calls (e.g. right after `exec = exec.with_target_partitions(self.target_partitions);`). This carries the same cell (and any already-built `Arc<lookup>`) into the re-planned exec.

- [ ] **Step 5: Build (lance-cache) to verify it compiles**

Run: `cargo build -p datafusion-bio-function-vep --features lance-cache 2>&1 | grep -E "error|Finished" | head`
Expected: `Finished` (the field is added but not yet used by the stream — may warn "field never read"; acceptable until Task 3).

- [ ] **Step 6: Commit**

```bash
git add datafusion/bio-function-vep/src/kv_cache/cache_exec.rs
git commit -m "feat(vep): add shared lance variation-lookup OnceCell to KvLookupExec"
```

---

## Task 2: Thread the shared cell into `KvLookupStream`, drop `lance_chroms`

**Files:** Modify `cache_exec.rs` — `KvLookupStream` struct (`:497`), stream init in `execute` (`:1903`).

- [ ] **Step 1: Replace the field in `KvLookupStream`**

In `struct KvLookupStream { … }` (`:497`), replace:

```rust
    #[cfg(feature = "lance-cache")]
    lance_chroms: HashMap<String, Option<SinglePathLanceVariationLookup>>,
    #[cfg(feature = "lance-cache")]
    lance_cursors: HashMap<String, usize>,
```

with:

```rust
    /// Shared (Arc-cloned from the exec) per-contig variation lookup, built once.
    #[cfg(feature = "lance-cache")]
    lance_lookup_cell: Arc<OnceCell<Arc<SinglePathLanceVariationLookup>>>,
    /// Per-partition forward walk cursor into the shared index (keyed by chrom).
    #[cfg(feature = "lance-cache")]
    lance_cursors: HashMap<String, usize>,
```

- [ ] **Step 2: Seed it in `execute`**

In `KvLookupExec::execute` where the `KvLookupStream { … }` literal is built (init at `:1903`), replace:

```rust
            #[cfg(feature = "lance-cache")]
            lance_chroms: HashMap::new(),
            #[cfg(feature = "lance-cache")]
            lance_cursors: HashMap::new(),
```

with:

```rust
            #[cfg(feature = "lance-cache")]
            lance_lookup_cell: Arc::clone(&self.lance_lookup_cell),
            #[cfg(feature = "lance-cache")]
            lance_cursors: HashMap::new(),
```

- [ ] **Step 3: Build (expect errors at the two use sites)**

Run: `cargo build -p datafusion-bio-function-vep --features lance-cache 2>&1 | grep -E "error\[|lance_chroms" | head`
Expected: errors at `ensure_cold_parquet_lookup` (`~2589`) and the resolve site (`~3642`) referencing the removed `lance_chroms`. These are fixed in Task 3. (Compile-unit with Task 3.)

- [ ] **Step 4: Commit (with Task 3 — single compile unit)**

Do not commit separately; proceed to Task 3 and commit together once it builds.

---

## Task 3: Build once via `get_or_try_init`; read from the shared cell

**Files:** Modify `cache_exec.rs` — `ensure_cold_parquet_lookup` lance arm (`:2589`), resolve site (`:3642`).

- [ ] **Step 1: Replace the lazy per-stream build with single-flight shared build**

In `ensure_cold_parquet_lookup`, replace the entire lance branch body (the `if !self.lance_chroms.contains_key(chrom) { … }` block, `:2589`–`:2632`) with a build into the shared cell. Keep the **identical** path + projection derivation:

```rust
        // Build the shared per-contig lookup exactly once (single-flight across
        // all partition streams). Concurrent fan-out probes await the one build.
        if self.lance_lookup_cell.get().is_none() {
            let cache = PartitionedLanceCache::detect(cache_root.to_string_lossy().as_ref())
                .ok_or_else(|| {
                    DataFusionError::Execution(format!(
                        "lance variation lookup but no variation.lance manifest under {}",
                        cache_root.display()
                    ))
                })?;
            let path = lance_variation_path_for_chrom(&cache, chrom).ok_or_else(|| {
                DataFusionError::Execution(format!(
                    "lance variation lookup but no Lance dataset for {chrom} in {}",
                    cache_root.display()
                ))
            })?;
            let projection = cold_parquet_projection_columns(cache_columns, collect_colocated);
            let cell = Arc::clone(&self.lance_lookup_cell);
            let open_started = Instant::now();
            let lookup_arc = tokio::task::block_in_place(|| {
                tokio::runtime::Handle::current().block_on(async {
                    cell.get_or_try_init(|| async {
                        SinglePathLanceVariationLookup::open(&path, projection)
                            .await
                            .map(Arc::new)
                    })
                    .await
                    .cloned()
                })
            })?;
            let open_elapsed = open_started.elapsed();
            if self.profile_enabled {
                self.profile.position_index_load += open_elapsed;
                self.profile.position_index_loaded += 1;
                self.profile.position_index_rows += lookup_arc.row_ids_len() as u64;
            }
            if self.profile_detailed || std::env::var_os("VEP_LANCE_PROFILE").is_some() {
                eprintln!(
                    "[vep-lance-profile] open chrom={} path(shared) projected_cols(shared) collect_colocated={} row_ids={} unique_positions={} open_s={:.3}",
                    chrom,
                    collect_colocated,
                    lookup_arc.row_ids_len(),
                    lookup_arc.unique_positions(),
                    open_elapsed.as_secs_f64(),
                );
            }
        }
        self.lance_cursors.entry(chrom.to_string()).or_insert(0);
```

Notes for the implementer:
- `path` and `projection` are **moved into** the async closure; that's why the profiling line above no longer prints `path`/`projected_cols` (they're consumed). If you want them in the log, bind `let projected_cols = projection.len(); let path_disp = path.display().to_string();` before the closure and use those — purely cosmetic.
- `get_or_try_init` runs the closure **once** even under concurrent calls on the shared cell; the other callers await and return the same `&Arc<…>`. `.cloned()` yields an owned `Arc` (cheap refcount bump) so the borrow of the cell doesn't outlive the `block_on`.
- The build's `open_elapsed`/profile counters now reflect the single build for the first stream and ~0 (cache hit, no closure run) for the rest — expected.

- [ ] **Step 2: Update the resolve site**

At the resolve call site (`:3642`), replace:

```rust
                    let lookup = self
                        .lance_chroms
                        .get(&chrom)
                        .and_then(Option::as_ref)
                        .ok_or_else(|| {
                            DataFusionError::Execution(format!(
                                "lance variation lookup not open for {chrom}"
                            ))
                        })?;
```

with:

```rust
                    let lookup = self
                        .lance_lookup_cell
                        .get()
                        .ok_or_else(|| {
                            DataFusionError::Execution(format!(
                                "lance variation lookup not built for {chrom}"
                            ))
                        })?;
```

`lookup` is now `&Arc<SinglePathLanceVariationLookup>`; the following `lookup.resolve_and_take(&starts, cursor)` call works unchanged (`Arc` derefs to the lookup, `resolve_and_take` takes `&self`).

- [ ] **Step 3: Build clean**

Run: `cargo build -p datafusion-bio-function-vep --features lance-cache 2>&1 | grep -E "error|warning: unused|field never read|Finished" | head -20`
Expected: `Finished`, no errors. If "unused import `HashMap`" appears and `HashMap` is no longer used elsewhere, leave it (it is used by many other fields); if clippy flags it later, Task 5 handles it.

- [ ] **Step 4: Run the full lib suite; compare failure set to baseline**

```bash
cargo test -p datafusion-bio-function-vep --features lance-cache --lib 2>&1 \
  | grep -E "^test .* FAILED" | grep -vE "test result:" | sed 's/ \.\.\. FAILED//' | sort -u > /tmp/sv_fail.txt
echo "failures: $(wc -l < /tmp/sv_fail.txt) (baseline 35)"
```
Expected: 35 failures, the same pre-existing set (capture a baseline with `git stash` if `/tmp/b.txt` from prior sessions is gone: stash this change, run the same command into `/tmp/b.txt`, pop). `comm -13 /tmp/b.txt /tmp/sv_fail.txt` must be empty (no new failures).

- [ ] **Step 5: Commit (Tasks 2+3)**

```bash
git add datafusion/bio-function-vep/src/kv_cache/cache_exec.rs
git commit -m "feat(vep): share per-contig lance variation lookup across partition streams

Replace per-stream lance_chroms (N copies of the 88M-entry PositionRowIdIndex)
with a shared Arc<OnceCell<Arc<SinglePathLanceVariationLookup>>> on KvLookupExec,
cloned into each KvLookupStream. get_or_try_init single-flights the build so it
runs exactly once even under simultaneous fan-out. Per-partition lance_cursors
unchanged. resolve_and_take / take_rows unchanged."
```

---

## Task 4: fmt + clippy

**Files:** `cache_exec.rs` as needed.

- [ ] **Step 1: Format**

Run: `cargo fmt`
Expected: no diff or whitespace only.

- [ ] **Step 2: Clippy (-D warnings)**

Run: `cargo clippy -p datafusion-bio-function-vep --features lance-cache --lib --example bench_annotate_vcf -- -D warnings 2>&1 | grep -E "error|warning|Finished" | head -20`
Expected: clean. Likely fixups: an unused import (`SinglePathLanceVariationLookup` is still used by the field type, so it stays), or a needless `.clone()`. Apply the minimal idiomatic fix clippy names.

- [ ] **Step 3: Commit**

```bash
git add -A
git commit -m "style(vep): fmt + clippy for shared variation-lookup cell"
```

---

## Task 5: Build the release bench + ensure input

**Files:** none (build/data).

- [ ] **Step 1: Build the bench**

Run: `cargo build --release --features lance-cache --example bench_annotate_vcf 2>&1 | tail -1`
Expected: `Finished release`.

- [ ] **Step 2: Ensure indexed input**

Run:
```bash
test -f /tmp/chr1_200k.vcf.gz && echo present || \
bcftools view /Users/mwiewior/workspace/data_vepyr/HG002_GRCh38_1_22_v4.2.1_benchmark.vcf.gz chr1 \
 | awk 'BEGIN{c=0} /^#/{print;next} {if(c<200000){print;c++}else{exit}}' \
 | bgzip > /tmp/chr1_200k.vcf.gz && tabix -p vcf -f /tmp/chr1_200k.vcf.gz
```
Expected: `present` or a fresh bgzipped+tabixed file.

---

## Task 6: MANDATORY parity gate — 0 mismatches at threads {1,2,4,8} + vs pre-change serial

**Files:** none (verification).

- [ ] **Step 1: Capture the pre-change baseline serial output**

Build the pre-change binary (the commit just before Task 1) and run threads=1 to `/tmp/sv_c1_baseline.vcf`:
```bash
PRE=$(git rev-parse HEAD~5)   # adjust: the commit before "add shared lance ... OnceCell"; verify with git log --oneline
git stash push -- datafusion/bio-function-vep/src/kv_cache/cache_exec.rs 2>/dev/null || true
git checkout "$PRE" -- datafusion/bio-function-vep/src/kv_cache/cache_exec.rs
cargo build --release --features lance-cache --example bench_annotate_vcf 2>&1 | tail -1
CACHE=/Users/mwiewior/workspace/data_vepyr/115_GRCh38_merged
FASTA=/Users/mwiewior/workspace/data_vepyr/Homo_sapiens.GRCh38.dna.primary_assembly.fa
./target/release/examples/bench_annotate_vcf --input /tmp/chr1_200k.vcf.gz --cache $CACHE \
  --output /tmp/sv_c1_baseline.vcf --backend lance --everything --reference-fasta $FASTA \
  --threads 1 --no-progress
git checkout HEAD -- datafusion/bio-function-vep/src/kv_cache/cache_exec.rs
git stash pop 2>/dev/null || true
cargo build --release --features lance-cache --example bench_annotate_vcf 2>&1 | tail -1
```
Expected: `/tmp/sv_c1_baseline.vcf` with 200000 data rows; working tree back to the shared-index version.

- [ ] **Step 2: Produce the four thread outputs (post-change)**

```bash
BIN=./target/release/examples/bench_annotate_vcf
CACHE=/Users/mwiewior/workspace/data_vepyr/115_GRCh38_merged
FASTA=/Users/mwiewior/workspace/data_vepyr/Homo_sapiens.GRCh38.dna.primary_assembly.fa
for T in 1 2 4 8; do
  $BIN --input /tmp/chr1_200k.vcf.gz --cache $CACHE --output /tmp/sv_c${T}.vcf \
    --backend lance --everything --reference-fasta $FASTA --threads $T --no-progress
  echo "t${T}: $(grep -vc '^#' /tmp/sv_c${T}.vcf) rows"
done
```
Expected: 200000 rows each.

- [ ] **Step 3: GATE — 0 mismatches across thread counts**

```bash
for T in 2 4 8; do
  if diff <(grep -v '^#' /tmp/sv_c1.vcf) <(grep -v '^#' /tmp/sv_c${T}.vcf) >/dev/null; then
    echo "t${T}: 0 mismatches vs t1"
  else echo "t${T}: MISMATCH — GATE FAILED"; diff <(grep -v '^#' /tmp/sv_c1.vcf) <(grep -v '^#' /tmp/sv_c${T}.vcf) | head; fi
done
```
Expected: `0 mismatches` for t2, t4, t8.

- [ ] **Step 4: GATE — 0 mismatches vs pre-change serial**

```bash
diff <(grep -v '^#' /tmp/sv_c1_baseline.vcf) <(grep -v '^#' /tmp/sv_c1.vcf) >/dev/null \
  && echo "serial output unchanged (0 mismatches)" || { echo "SERIAL REGRESSION — GATE FAILED"; diff <(grep -v '^#' /tmp/sv_c1_baseline.vcf) <(grep -v '^#' /tmp/sv_c1.vcf) | head; }
```
Expected: `serial output unchanged (0 mismatches)`. (Header `##FORMAT` line ordering is non-deterministic passthrough — ignore; the gate is data rows only.)

Any mismatch → STOP, debug with superpowers:systematic-debugging. Do not proceed.

---

## Task 7: Benchmark — RSS + timing, threads {1,2,4,8}, A/B vs pre-change

**Files:** none (measurement).

- [ ] **Step 1: Stage both binaries**

```bash
cp ./target/release/examples/bench_annotate_vcf /tmp/sv_bench_after
PRE=$(git rev-parse HEAD~6)   # the commit before Task 1; verify via git log --oneline
git checkout "$PRE" -- datafusion/bio-function-vep/src/kv_cache/cache_exec.rs
cargo build --release --features lance-cache --example bench_annotate_vcf 2>&1 | tail -1
cp ./target/release/examples/bench_annotate_vcf /tmp/sv_bench_before
git checkout HEAD -- datafusion/bio-function-vep/src/kv_cache/cache_exec.rs
cargo build --release --features lance-cache --example bench_annotate_vcf 2>&1 | tail -1
```
Expected: `/tmp/sv_bench_before` and `/tmp/sv_bench_after` exist; working tree restored to shared-index version.

- [ ] **Step 2: Measure peak RSS + wall time, both binaries, t1/t2/t4/t8 (min of 2)**

```bash
CACHE=/Users/mwiewior/workspace/data_vepyr/115_GRCh38_merged
FASTA=/Users/mwiewior/workspace/data_vepyr/Homo_sapiens.GRCh38.dna.primary_assembly.fa
measure() { # $1=binary $2=threads -> "wall_s rss_gb" (min over 2 reps)
  local bw=99999 br=1e18
  for r in 1 2; do
    /usr/bin/time -l "$1" --input /tmp/chr1_200k.vcf.gz --cache $CACHE --output /tmp/sv_m.vcf \
      --backend lance --everything --reference-fasta $FASTA --threads "$2" --no-progress >/tmp/sv_m.log 2>&1
    local w=$(grep -E '^[[:space:]]*real|elapsed' /tmp/sv_m.log | head -1 | awk '{print $1}')
    [ -z "$w" ] && w=$(grep real /tmp/sv_m.log | awk '{print $2}')
    local b=$(grep "maximum resident set size" /tmp/sv_m.log | awk '{print $1}')
    awk -v a="$w" -v c="$bw" 'BEGIN{exit !(a<c)}' && bw=$w
    awk -v a="$b" -v c="$br" 'BEGIN{exit !(a<c)}' && br=$b
  done
  awk -v w="$bw" -v r="$br" 'BEGIN{printf "%6ss  %5.2f GB", w, r/1073741824}'
}
printf "%-8s %-22s %-22s\n" threads BEFORE AFTER
for T in 1 2 4 8; do
  printf "%-8s %-22s %-22s\n" "$T" "$(measure /tmp/sv_bench_before $T)" "$(measure /tmp/sv_bench_after $T)"
done
rm -f /tmp/sv_m.vcf /tmp/sv_m.log
```
Expected (success criteria):
- **t8 peak RSS: ~23 GB (before) → ~6–7 GB (after).** Primary win.
- t1 RSS roughly unchanged (~5 GB).
- Wall time t1–t8: unchanged or slightly better (one BTree decode at startup instead of N).

Note: `/usr/bin/time` on macOS prints wall as the line `<real> real` and RSS as `<bytes>  maximum resident set size`; the helper extracts both. If your `time` formatting differs, capture RSS with `/usr/bin/time -l` (always present) and wall via the bench's own emitted timing or `date`-bracketing.

- [ ] **Step 3: Record results + commit a results note**

Append the measured BEFORE/AFTER table to the design spec under a "Results" heading, then:
```bash
git add docs/superpowers/specs/2026-06-18-shared-variation-index-design.md
git commit -m "docs(vep): shared variation-index measured results (t8 RSS <before> -> <after>)"
rm -f /tmp/sv_bench_before /tmp/sv_bench_after
```

---

## Self-Review

- **Spec coverage:** §4 architecture (per-contig shared lookup) → Tasks 1–3. §5 component changes (`KvLookupExec` field, stream field, `ensure_cold_parquet_lookup`, resolve site) → Tasks 1–3 (implemented via shared `OnceCell` rather than eager-in-`prepare_contig_context`; architecture note documents the equivalence). §7 concurrency safety (read-only index, per-stream cursor, lance concurrent reads) → preserved: `lance_cursors` untouched, `resolve_and_take(&self)` untouched. §9 gates → Tasks 6 (parity) + 7 (RSS/timing). Non-goals (backend removal, region-scoping, SIFT) → untouched.
- **Placeholder scan:** none — every code step shows the full replacement. `HEAD~N` offsets in Tasks 6–7 are annotated "verify via git log" because exact depth depends on how many commits land (fmt/clippy may add one); the implementer resolves the actual pre-Task-1 commit hash.
- **Type consistency:** `lance_lookup_cell: Arc<OnceCell<Arc<SinglePathLanceVariationLookup>>>` used identically in `KvLookupExec` (Task 1), `KvLookupStream` (Task 2), and both build/read sites (Task 3). `get_or_try_init` returns `Result<&Arc<…>>`; `.cloned()` → owned `Arc`. Resolve yields `&Arc<…>`, deref-compatible with `resolve_and_take(&self, …)`.
- **Known risk:** `with_new_children` re-planning must carry the cell (Task 1 Step 4) or a re-plan would silently rebuild per-stream — covered, and the parity/RSS gates would catch a regression.
