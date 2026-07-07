# Arrow Zero-Copy Colocated AF — Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: superpowers:subagent-driven-development or superpowers:executing-plans. Steps use `- [ ]`.

**Goal:** Replace colocated existing-variant `af_values: Vec<String>` (and the other colocated string fields) with zero-copy references into the Arrow `StringArray`s the lookup already produces, eliminating ~17M+ per-variant string allocations, byte-identically.

**Architecture:** The lookup `take_rows` → `unbundle_af_columns` yields 27 `Utf8` columns of matched rows. Store `Arc`'d column refs + a row index per colocated entry; expose `&str` accessors. Removes both the lookup-collect copy and the `build_colocated_map_from_sink` clone.

**Tech Stack:** Rust 2024, Arrow 58 (`StringArray`, `ArrayRef`), DataFusion 53.

## Global Constraints

- **Byte-identical** output: chr1 t8 e2e 0 CSQ mismatch, 86/86 fields 100%.
- Absent/null AF → `""` (matches today's `String::new()`).
- AF index order = `AF_COLUMNS` / `AF_COL_NAMES`, unchanged.
- Per task: `cargo test -p datafusion-bio-function-vep` + `cargo clippy -p datafusion-bio-function-vep -- -D warnings` green.
- Measure with `dhat-heap` (vepyr feature) + `[VEP_RSS]` probes (already present, uncommitted) before/after.
- No on-disk format change; public API unchanged (rebuild vepyr only).

---

## File Structure

- `colocated.rs` — `ColocatedCacheEntry`, `ColocatedSinkValue`: replace `af_values: Vec<String>` (+ `variation_name`/`allele_string`/`clin_sig`/`clin_sig_allele`/`pubmed`) with Arrow-backed `(AfColumns, row)` + accessors. Add `AfColumns`.
- `annotate_provider.rs` — `ColocatedEntry` mirror; `build_colocated_map_from_sink` (`:2095`) → `Arc::clone`; consumers (`:1649`, `:1962`, `:1671`) → accessor calls.
- Lookup collect site (Task 0 locates: `variant_lookup_exec.rs` and/or `lance_cache/lookup_exec.rs`) — store array refs + row instead of `.to_string()`.

---

## Task 0: Spike — locate exact collect/sink/consume sites

**Files:** read-only; write findings into the spec §5.

- [ ] **Step 1:** Find where `ColocatedCacheEntry.af_values` is built from Arrow:
  `grep -rn "ColocatedCacheEntry {" datafusion/bio-function-vep/src/` and read the constructor (which AF columns, how nulls handled, what batch/arrays are in scope).
- [ ] **Step 2:** Read `ColocatedSinkValue` (`colocated.rs:330`) and the sink drain in `build_colocated_map_from_sink` (`annotate_provider.rs:2095-2160`) — confirm where the second copy / `.clone()` happens and what owns the arrays.
- [ ] **Step 3:** Read all consumers of `af_values` (`annotate_provider.rs:1649, 1671, 1962`) and the other string fields — confirm they only need `&str`.
- [ ] **Step 4:** Confirm the take_rows batch columns are `StringArray` after `unbundle_af_columns`, and whether nulls vs `""` are used. Write exact file:line + types into spec §5/§6 (R3). Decide: retain projected columns vs the whole batch (R1).
- [ ] **Step 5:** Commit findings (docs only).

---

## Task 1: `AfColumns` + accessors on `ColocatedCacheEntry`

**Files:** `colocated.rs` (+ test there).

**Interfaces produced:** `AfColumns(Arc<Vec<StringArray>>)`; `ColocatedCacheEntry { af: AfColumns, af_row: u32, ... }` with `fn af_value(&self, col: usize) -> &str` and `fn af_all_empty(&self) -> bool`.

- [ ] **Step 1: Failing test** (in `colocated.rs` tests): build an `AfColumns` from two `StringArray`s (one with a null, one `""`), wrap in a `ColocatedCacheEntry`, assert `af_value(0)` returns the value, `af_value(null_col)` returns `""`, and `af_all_empty()` correct.
- [ ] **Step 2:** Run → fails (type/methods absent).
- [ ] **Step 3:** Implement `AfColumns` + the field swap + accessors:
```rust
#[derive(Clone)]
pub struct AfColumns(pub Arc<Vec<StringArray>>);
impl ColocatedCacheEntry {
    pub fn af_value(&self, col: usize) -> &str {
        match self.af.0.get(col) {
            Some(a) if !a.is_null(self.af_row as usize) => a.value(self.af_row as usize),
            _ => "",
        }
    }
    pub fn af_len(&self) -> usize { self.af.0.len() }
    pub fn af_all_empty(&self) -> bool {
        (0..self.af.0.len()).all(|c| self.af_value(c).is_empty())
    }
}
```
- [ ] **Step 4:** Run → passes.
- [ ] **Step 5:** Commit.

---

## Task 2: Build colocated entries from Arrow refs at the collect site

**Files:** the collect site from Task 0; `colocated.rs`.

**Interfaces consumed:** `AfColumns` (Task 1).

- [ ] **Step 1: Failing test:** a collect-path unit test (or extend an existing colocated-collection test) that feeds a take_rows-shaped batch and asserts the resulting entry's `af_value(idx)` equals the input array value — and that NO `.to_string()` of AF happens (assert via the value being a borrow: compare `as_ptr` of `af_value` into the array's data buffer if feasible, else just value-equality).
- [ ] **Step 2:** Run → fails.
- [ ] **Step 3:** At the collect site, project the needed AF (+ string) columns into `AfColumns` once per batch (`Arc::clone` the columns), and store `(af.clone(), row_idx)` per matched variant instead of `.to_string()` loops. Apply the same `(StringArray, row)` swap to `variation_name`/`allele_string`/`clin_sig`/`clin_sig_allele`/`pubmed`.
- [ ] **Step 4:** Run → passes.
- [ ] **Step 5:** Commit.

---

## Task 3: `ColocatedEntry` (annotate side) + drop the second copy

**Files:** `annotate_provider.rs`.

- [ ] **Step 1: Failing test:** unit test of `build_colocated_map_from_sink` (`:2095`) asserting the produced `ColocatedEntry.af_value(idx)` matches the sink entry, and that the conversion does not deep-copy AF (the entry's `af` Arc points at the same allocation as the sink's — assert `Arc::ptr_eq`).
- [ ] **Step 2:** Run → fails.
- [ ] **Step 3:** Change `ColocatedEntry` to `{ af: AfColumns, af_row: u32, ... }` with the same accessors; make `build_colocated_map_from_sink` move/`Arc::clone` the refs (the `:2130` `.clone()` becomes a ref bump). Keep `ColocatedFrequencyFields` output path (`:2056`) reading via `af_value`.
- [ ] **Step 4:** Run → passes.
- [ ] **Step 5:** Commit.

---

## Task 4: Switch all consumers to the accessor

**Files:** `annotate_provider.rs` (`:1649`, `:1671`, `:1962`, and any others from Task 0).

- [ ] **Step 1: Failing test:** an annotation-level test that exercises a variant with mixed populated/empty AF and asserts the CSQ frequency fields are byte-identical to a captured expected (use an existing golden-style test if present).
- [ ] **Step 2:** Run → fails to compile (field removed) or mismatches.
- [ ] **Step 3:** Replace `&entry.af_values[idx]` / `af_values.get(idx).map(String::as_str)` / `af_values.iter().all(String::is_empty)` with `entry.af_value(idx)` / `entry.af_all_empty()`. Remove the now-dead `colocated_map_bytes` `af_values` references (or update the probe).
- [ ] **Step 4:** Run → passes; `cargo build` clean.
- [ ] **Step 5:** Commit.

---

## Task 5: e2e parity + allocation/RSS measurement gate

**Files:** results into spec §7.

- [ ] **Step 1:** Rebuild vepyr (`uv run maturin develop --release`).
- [ ] **Step 2:** chr1 t8 parity (`VEP_PROFILE=1 ... --threads 8`): assert 0 CSQ mismatch, 86/86 fields (G1).
- [ ] **Step 3:** dhat: rebuild `--features dhat-heap`, run chr1 t2, record `Total` blocks + `t-gmax` live (G2) vs baseline (466M / 4.15 GB).
- [ ] **Step 4:** Peak RSS: `/usr/bin/time -l` chr1 t8 before/after (G3).
- [ ] **Step 5:** Record G1–G4 in spec §7; update memory `vep-peak-rss-allocation-churn`. Commit.

---

## Self-Review

Spec §1 chain → Tasks 1–4 (both copies + consumers); §3 design → Tasks 1–3; §4 invariants → tests in Tasks 1/4; §6 risks → R1 Task 0/2 (project columns), R2 Task 0, R3 Task 1 (null→""), R4 Task 3 (`Arc` owned); §7 gates → Task 5. `matched_alleles` + CSQ-output buffer reuse explicitly deferred (spec §2 non-goals) — note as follow-up if G2/G3 fall short.

**Follow-up (if peak still high after AF fix):** apply the same Arrow-ref pattern to `variation_name`/`allele_string`/etc (already in Task 2/3 scope), then `matched_alleles` (SmallVec/arena) and CSQ-output scratch-buffer reuse (`:2056`).
