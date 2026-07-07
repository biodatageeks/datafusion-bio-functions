# CSQ Output-Assembly Churn Reduction — Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Cut the ~60% of per-variant allocations in the in-repo CSQ output assembly (`annotate_provider.rs`), byte-identically, as staged + individually-measured steps.

**Architecture:** Four independent byte-identical refactors, each gated by a fast chr22 output byte-diff and a dhat re-measure: (1) `push_unique_value(&str)`, (2) eliminate per-transcript throwaway temporaries feeding `write!(csq_buf,…)`, (3) de-duplicate the typed-column pass's recomputation, (4) reuse colocated `frequency_fields`/`variant_fields` scratch across rows. Measure after each; stop when ROI flattens.

**Tech Stack:** Rust 2024, Arrow 58, DataFusion 53; vepyr (sibling, path-dep) for e2e runs; dhat-heap for churn measurement.

**Spec:** `docs/superpowers/specs/2026-06-19-csq-output-churn-design.md`.

## Global Constraints

- **Byte-identical** output: each step's chr22 output VCF must be `diff`-identical to the Task-0 baseline; final chr1 t8 e2e = 0 CSQ mismatch, 86/86 fields 100%.
- Per task: `cargo test -p datafusion-bio-function-vep` + `cargo clippy -p datafusion-bio-function-vep -- -D warnings` green.
- No change to CSQ field semantics/order, the output RecordBatch schema, the public API, or the on-disk format. Rebuild vepyr only.
- Cross-crate VCF byte serialization (`vcf_sink`/sibling `datafusion_bio_format_vcf::serializer`) is OUT of scope (spec §1/§2).
- Reused scratch buffers MUST be `.clear()`-ed before each use (no stale-data leak).
- **Stop-when-flat rule:** after each step's dhat re-measure, if the `Total` blocks reduction is negligible, stop and report rather than proceeding on faith (spec §4.3).
- Baseline (chr22, single-thread, current HEAD `f0dbb1b`): dhat `Total` = **76,070,619 blocks**, `t-gmax` = **1.75 GB**; 50,861 variants.
- Repo root: `/Users/mwiewior/research/git/datafusion-bio-functions`. vepyr: `/Users/mwiewior/research/git/vepyr`. vepyr scripts: `/Users/mwiewior/research/git/vepyr/e2e-testing/scripts`.

---

## File Structure

- `datafusion/bio-function-vep/src/annotate_provider.rs` — all four steps. `push_unique_value` (`:1749`), CSQ-text loop (`:5930–6252`), typed-column loop (`:6315–6900`), `ColocatedData::frequency_fields` (`:1955–2089`) / `variant_fields` (`:1858–1942`), per-batch scratch declarations (`:5563–5572`), per-batch context for new scratch.
- `datafusion/bio-function-vep/src/annotate_provider.rs` `#[cfg(test)] mod tests` (`:~12262`) — characterization unit test for `push_unique_value`.
- `/tmp/csq_step_gate.sh` (created Task 0) — reusable chr22 byte-diff + optional dhat gate.

> **Note on code completeness:** Step 1 includes verbatim before/after code. Steps 2–4 act on many sites inside a 12k-line function; each lists the EXACT `file:line` sites and the transformation rule, and instructs the implementer to read those lines first, then apply the rule. The chr22 byte-diff (Task 0 harness) is the hard correctness gate for every step.

---

## Task 0: Baseline capture + measurement harness

**Files:**
- Create: `/tmp/csq_step_gate.sh`
- Create (artifacts): `/tmp/csq_baseline_chr22.vcf`

**Interfaces:**
- Produces: `csq_step_gate.sh <label>` — annotates chr22 (non-dhat release), diffs output VCF vs baseline, prints `BYTE_IDENTICAL` / `BYTE_DIFF`.

- [ ] **Step 1: Ensure a clean non-dhat release of current HEAD.** The working tree may have a dhat-instrumented `.so` installed from prior measurement.

```bash
cd /Users/mwiewior/research/git/vepyr
env -u VIRTUAL_ENV -u CONDA_PREFIX uv run maturin develop --release > /tmp/csq_build_base.log 2>&1
tail -1 /tmp/csq_build_base.log   # expect: 🛠 Installed vepyr-0.1.0
```

- [ ] **Step 2: Capture the chr22 baseline output VCF.**

```bash
cd /Users/mwiewior/research/git/vepyr/e2e-testing/scripts
env -u VEP_PROFILE -u VIRTUAL_ENV -u CONDA_PREFIX uv run python run_annotation_fast.py \
  chr22 --cache merged --forks 0 --threads 1 --force --backend lance --skip-compare > /tmp/csq_base_run.out 2>/dev/null
# locate the produced output VCF (name pattern vepyr_lance_chr22_merged.vcf); copy it
OUT=$(grep -oE '[^ ]*vepyr_lance_chr22_merged\.vcf' /tmp/csq_base_run.out | head -1)
echo "output=$OUT"; cp "$OUT" /tmp/csq_baseline_chr22.vcf
wc -l /tmp/csq_baseline_chr22.vcf
```
Expected: a non-empty VCF; record its path in `$OUT` for reuse.

- [ ] **Step 3: Write the reusable gate script.**

```bash
cat > /tmp/csq_step_gate.sh <<'EOF'
#!/bin/zsh
# Usage: csq_step_gate.sh <label>   (run AFTER rebuilding vepyr release)
set -u
S=/Users/mwiewior/research/git/vepyr/e2e-testing/scripts
cd $S
env -u VEP_PROFILE -u VIRTUAL_ENV -u CONDA_PREFIX uv run python run_annotation_fast.py \
  chr22 --cache merged --forks 0 --threads 1 --force --backend lance --skip-compare > /tmp/csq_${1}_run.out 2>/dev/null
OUT=$(grep -oE '[^ ]*vepyr_lance_chr22_merged\.vcf' /tmp/csq_${1}_run.out | head -1)
if diff -q "$OUT" /tmp/csq_baseline_chr22.vcf >/dev/null; then
  echo "BYTE_IDENTICAL ($1)"
else
  echo "BYTE_DIFF ($1) -- FAIL"; diff "$OUT" /tmp/csq_baseline_chr22.vcf | head -20
fi
EOF
chmod +x /tmp/csq_step_gate.sh
```

- [ ] **Step 4: Sanity-check the gate against the unchanged build (must be identical).**

```bash
cd /Users/mwiewior/research/git/vepyr && env -u VIRTUAL_ENV -u CONDA_PREFIX uv run maturin develop --release >/dev/null 2>&1
zsh /tmp/csq_step_gate.sh sanity
```
Expected: `BYTE_IDENTICAL (sanity)`.

- [ ] **Step 5: Commit the harness reference (docs only — script lives in /tmp, record it in the plan dir).**

```bash
cd /Users/mwiewior/research/git/datafusion-bio-functions
# (no repo files change in Task 0; nothing to commit unless you choose to vendor the script)
echo "Task 0 baseline captured: /tmp/csq_baseline_chr22.vcf + /tmp/csq_step_gate.sh"
```

---

## Task 1: `push_unique_value(&str)` — allocate only on insert

**Files:**
- Modify: `datafusion/bio-function-vep/src/annotate_provider.rs:1749-1754` (fn) + call sites `:1897`, `:2027`, `:2045`
- Test: `annotate_provider.rs` `#[cfg(test)] mod tests`

**Interfaces:**
- Produces: `fn push_unique_value(values: &mut Vec<String>, value: &str)` (was `value: impl Into<String>`).

- [ ] **Step 1: Write a characterization test** (pins dedup + first-seen order; passes before and after — guards the refactor). Add to the tests module:

```rust
#[test]
fn push_unique_value_dedups_and_preserves_first_seen_order() {
    let mut v: Vec<String> = Vec::new();
    push_unique_value(&mut v, "a");
    push_unique_value(&mut v, "b");
    push_unique_value(&mut v, "a"); // duplicate -> ignored
    push_unique_value(&mut v, "c");
    assert_eq!(v, vec!["a".to_string(), "b".to_string(), "c".to_string()]);
}
```

- [ ] **Step 2: Run it against the current impl (should already pass — it's a characterization test).**

Run: `cargo test -p datafusion-bio-function-vep --lib -- push_unique_value_dedups`
Expected: PASS (with the current `impl Into<String>` signature, called with `"a"` etc.).
If it does NOT compile because the current signature differs, that confirms the call-shape; proceed to Step 3.

- [ ] **Step 3: Read the current fn + call sites, then refactor.** Read `annotate_provider.rs:1749-1754` and the three call sites at `:1897`, `:2027`, `:2045`. Replace the fn with:

```rust
fn push_unique_value(values: &mut Vec<String>, value: &str) {
    if !values.iter().any(|existing| existing == value) {
        values.push(value.to_string());
    }
}
```
Then at each call site pass a `&str` and DROP the `.clone()`/`.to_string()`:
- `:1897` `push_unique_value(&mut …, value.clone())` → `push_unique_value(&mut …, &value)` (or `value.as_str()`).
- `:2027` `push_unique_value(&mut per_column[idx], chosen.clone())` → `push_unique_value(&mut per_column[idx], &chosen)`.
- `:2045` `push_unique_value(&mut …, pop_name.to_string())` → `push_unique_value(&mut …, pop_name)` (if `pop_name` is already `&str`) or `&pop_name`.

Read the exact current argument expressions before editing — pass whatever yields `&str` without allocating.

- [ ] **Step 4: Run tests + clippy.**

Run: `cargo test -p datafusion-bio-function-vep --lib -- push_unique_value_dedups` → PASS
Run: `cargo clippy -p datafusion-bio-function-vep --lib -- -D warnings` → clean

- [ ] **Step 5: Rebuild vepyr + byte-diff gate.**

```bash
cd /Users/mwiewior/research/git/vepyr && env -u VIRTUAL_ENV -u CONDA_PREFIX uv run maturin develop --release >/dev/null 2>&1
zsh /tmp/csq_step_gate.sh step1
```
Expected: `BYTE_IDENTICAL (step1)`.

- [ ] **Step 6: Commit.**

```bash
cd /Users/mwiewior/research/git/datafusion-bio-functions
git add datafusion/bio-function-vep/src/annotate_provider.rs
git commit -m "perf(vep): push_unique_value takes &str, allocates only on insert

Drops the per-call String alloc-then-maybe-drop (~0.87M dhat blocks). Byte-
identical (chr22 diff). Characterization test for dedup/order.

Co-Authored-By: Claude Opus 4.8 (1M context) <noreply@anthropic.com>
Claude-Session: https://claude.ai/code/session_01AC49866apkEsgUTSmYDseE"
```

---

## Task 2: Eliminate per-transcript throwaway temporaries (CSQ-text loop)

**Files:**
- Modify: `datafusion/bio-function-vep/src/annotate_provider.rs:5930-6252` (per-transcript CSQ-text loop) + add scratch near `:5563-5572`.

**Interfaces:**
- Consumes: the reused `csq_buf` (`:5563`), `terms_buf` (`:5932`).
- Produces: no signature changes; same `csq_builder.append_value(&csq_buf)` output.

**Transformation rule:** Each temporary below is built as an owned `String` solely to be interpolated into a `write!(csq_buf, …)`. For each, write the underlying value directly into `csq_buf` (or pass a `&str`/borrowed `Cow`) instead of allocating. Integers: `write!`-format the integer (drop `.to_string()`). Uppercased text (`bam_edit`): uppercase into a single reused `String` scratch declared near `:5563` and `.clear()`-ed per transcript. `Cow`-returning escapes: keep borrowed; do not `.into_owned()`. For functions that currently return owned `String` (sift/polyphen/domains/mirna/hgvsp), only convert the hot+safe ones to write-into-buffer or `&str`/`Cow`; leave the rest if conversion risks lifetime/byte changes (R4) — the goal is removing allocs, not a rewrite.

**Exact sites to address (read each before editing):** `:5973` distance, `:5981` hgvsp, `:6001` bam_edit, `:6014` refseq_offset, `:6026` tsl_str, `:6043`/`:6045` swissprot/trembl (Cow), `:6071` hgvs_offset, `:6098` appris_str, `:6104` sift_str/polyphen_str, `:6132` domains, `:6150` mirna_str.

- [ ] **Step 1: Read the loop `:5930-6252`** and list, for each site above, the current temporary and the exact `write!(csq_buf, …)` it feeds. Confirm `csq_buf` is `.clear()`-ed per row (`:5705`) — the per-transcript writes append within a row's CSQ list separated as today.

- [ ] **Step 2: Add reused scratch** near `:5563` if needed for the uppercase case, e.g.:

```rust
let mut upper_buf = String::with_capacity(32);
```
(Only add scratch you actually use; `.clear()` before each use.)

- [ ] **Step 3: Apply the transformation rule** to the integer/`&str`/`Cow` sites first (lowest risk): `distance`, `refseq_offset`, `tsl_str`, `hgvs_offset` (integers → direct `write!`), `swissprot`/`trembl` (keep `Cow` borrowed). Then `bam_edit` via `upper_buf`. Then evaluate sift/polyphen/domains/mirna/hgvsp case-by-case.

- [ ] **Step 4: Build + unit tests + clippy.**

Run: `cargo test -p datafusion-bio-function-vep --lib` → all pass
Run: `cargo clippy -p datafusion-bio-function-vep --lib -- -D warnings` → clean

- [ ] **Step 5: Rebuild vepyr + byte-diff gate (THE correctness gate for this step).**

```bash
cd /Users/mwiewior/research/git/vepyr && env -u VIRTUAL_ENV -u CONDA_PREFIX uv run maturin develop --release >/dev/null 2>&1
zsh /tmp/csq_step_gate.sh step2
```
Expected: `BYTE_IDENTICAL (step2)`. If `BYTE_DIFF`, revert the offending site (the diff output names the field) — a temporary likely changed formatting (e.g. integer width, escaping).

- [ ] **Step 6: dhat re-measure + decision.**

```bash
cd /Users/mwiewior/research/git/vepyr && env -u VIRTUAL_ENV -u CONDA_PREFIX uv run maturin develop --release --features dhat-heap >/dev/null 2>&1
cd /Users/mwiewior/research/git/vepyr/e2e-testing/scripts && rm -f dhat-heap.json
env -u VIRTUAL_ENV -u CONDA_PREFIX uv run python run_annotation_fast.py chr22 --cache merged --forks 0 --threads 1 --force --backend lance --skip-compare >/dev/null 2>/tmp/csq_dhat_step2.err
grep -E "dhat: Total|t-gmax" /tmp/csq_dhat_step2.err
python3 /tmp/dhat_attribute.py dhat-heap.json | head -12
```
Record `Total` blocks vs 76.07M baseline. **Decision:** if the reduction is negligible, stop here and report; else continue.

- [ ] **Step 7: Rebuild non-dhat release (so later steps' byte-diff runs are fast) + commit.**

```bash
cd /Users/mwiewior/research/git/vepyr && env -u VIRTUAL_ENV -u CONDA_PREFIX uv run maturin develop --release >/dev/null 2>&1
cd /Users/mwiewior/research/git/datafusion-bio-functions
git add datafusion/bio-function-vep/src/annotate_provider.rs
git commit -m "perf(vep): write CSQ-text fragments directly, drop per-transcript temporaries

Byte-identical (chr22 diff); dhat Total <RECORD> blocks (was 76.07M).

Co-Authored-By: Claude Opus 4.8 (1M context) <noreply@anthropic.com>
Claude-Session: https://claude.ai/code/session_01AC49866apkEsgUTSmYDseE"
```

---

## Task 3: De-duplicate the typed-column pass

**Files:**
- Modify: `datafusion/bio-function-vep/src/annotate_provider.rs:6315-6900` (typed-column loop) + per-transcript scratch.

**Interfaces:**
- Consumes: per-transcript values already computed in the CSQ-text loop (Task 2).
- Produces: same typed-column builder outputs (`b_consequence` etc.).

**Transformation rule:** The typed-column loop is a 2nd pass over `sorted_indices` that recomputes values the CSQ-text loop already produced — `bam_edit` (`:6334`), `refseq_offset` (`:6322`), and `terms.iter().collect::<Vec<&str>>().join("&")` (`:6341`, duplicating `terms_buf` from `:5932`). Cache these per-transcript during the CSQ-text loop into a reusable `Vec` of small structs (cleared per row), then read them here instead of recomputing. Do NOT merge the two passes (byte-identity risk) — only remove the redundant recomputation.

- [ ] **Step 1: Read both loops** (`:5930-6252` and `:6315-6900`) and identify the exact recomputed values + their types. Define a per-row scratch, e.g.:

```rust
struct TxScratch { terms_joined: String, bam_edit: String, refseq_offset: String }
let mut tx_scratch: Vec<TxScratch> = Vec::new(); // declared near :5563, cleared per row
```
(Adjust fields to exactly what is recomputed; reuse Strings via `.clear()` where possible.)

- [ ] **Step 2: In the CSQ-text loop**, after computing `terms`/`bam_edit`/`refseq_offset`, store them into `tx_scratch[transcript_position]` (push or index by the same iteration order used by the typed-column loop).

- [ ] **Step 3: In the typed-column loop**, replace the recomputation at `:6322`/`:6334`/`:6341` with reads from `tx_scratch`. For `:6341`, append the cached `terms_joined` instead of re-`join`ing.

- [ ] **Step 4: Build + unit tests + clippy.**

Run: `cargo test -p datafusion-bio-function-vep --lib` → pass
Run: `cargo clippy -p datafusion-bio-function-vep --lib -- -D warnings` → clean

- [ ] **Step 5: byte-diff gate.**

```bash
cd /Users/mwiewior/research/git/vepyr && env -u VIRTUAL_ENV -u CONDA_PREFIX uv run maturin develop --release >/dev/null 2>&1
zsh /tmp/csq_step_gate.sh step3
```
Expected: `BYTE_IDENTICAL (step3)`. If `BYTE_DIFF`, the cache indexing or the two passes' iteration order diverged — verify both iterate `sorted_indices` identically.

- [ ] **Step 6: dhat re-measure + decision** (same commands as Task 2 Step 6, output `/tmp/csq_dhat_step3.err`). Record vs baseline; decide continue/stop.

- [ ] **Step 7: Rebuild non-dhat release + commit** (same shape as Task 2 Step 7; message: `perf(vep): cache per-transcript values, drop typed-column recomputation`).

---

## Task 4: Reuse colocated scratch across rows

**Files:**
- Modify: `datafusion/bio-function-vep/src/annotate_provider.rs:1858-2089` (`variant_fields`/`frequency_fields`) + their per-row call site `:5664-5681` + per-batch scratch.

**Interfaces:**
- Consumes: `ColocatedData` (unchanged).
- Produces: same `ColocatedFrequencyFields`/`ColocatedVariantFields` *values*; only allocation lifetime changes (scratch reused across rows).

**Transformation rule:** `frequency_fields` (`:1955`) allocates per row: `per_column: Vec<Vec<String>>` (`:1961`) and, per (entry × AF column), a `HashMap<String,String>` (`:1987`) + `HashSet<String>` (`:1988`). Hoist these into reusable scratch owned by the per-batch loop context, passed in by `&mut`, and `.clear()`-ed per row/column instead of re-allocated. Same for `variant_fields`'s per-entry `HashMap` (`:1885`). Pool the output `ColocatedFrequencyFields`/`ColocatedVariantFields` (reuse + clear `af_values`/`max_af`/`max_af_pops`) rather than constructing fresh per row.

- [ ] **Step 1: Read `:1858-2089`** and the call site `:5664-5681`. Decide the scratch carrier: either add fields to an existing per-batch context struct or pass `&mut` scratch params into `frequency_fields`/`variant_fields`. Keep the functions' return *values* identical.

- [ ] **Step 2: Introduce reusable scratch** near `:5563` (per-batch), e.g.:

```rust
let mut per_column_scratch: Vec<Vec<String>> = vec![Vec::new(); AF_COLUMNS.len()];
let mut freq_map_scratch: HashMap<String, String> = HashMap::new();
let mut allele_set_scratch: HashSet<String> = HashSet::new();
```
Change `frequency_fields`/`variant_fields` to take `&mut` these and `.clear()` each at the top of every row/column iteration (clear the inner `Vec`s in `per_column_scratch`, not reallocate the outer).

- [ ] **Step 3: Apply** — replace the per-row/per-column `vec!`/`HashMap::new()`/`HashSet::new()` allocations (`:1961`, `:1987`, `:1988`, `:1885`) with `.clear()` on the passed-in scratch. Preserve the `.join("&")` outputs exactly.

- [ ] **Step 4: Build + unit tests + clippy.** The existing colocated tests (`colocated_repeat_indel_id_fallback_does_not_drive_frequency_fields`, and the AF tests) must still pass.

Run: `cargo test -p datafusion-bio-function-vep --lib` → pass
Run: `cargo clippy -p datafusion-bio-function-vep --lib -- -D warnings` → clean

- [ ] **Step 5: byte-diff gate.**

```bash
cd /Users/mwiewior/research/git/vepyr && env -u VIRTUAL_ENV -u CONDA_PREFIX uv run maturin develop --release >/dev/null 2>&1
zsh /tmp/csq_step_gate.sh step4
```
Expected: `BYTE_IDENTICAL (step4)`. If `BYTE_DIFF`, a scratch buffer wasn't cleared (stale data) — audit every `.clear()`.

- [ ] **Step 6: dhat re-measure + decision** (output `/tmp/csq_dhat_step4.err`). Record vs baseline.

- [ ] **Step 7: Rebuild non-dhat release + commit** (message: `perf(vep): reuse colocated frequency/variant scratch across rows`).

---

## Task 5: Final gate — parity + throughput + RSS + write-up

**Files:**
- Modify (docs): `docs/superpowers/specs/2026-06-19-csq-output-churn-design.md` §8 (record G1–G4); memory note `vep-peak-rss-allocation-churn`.

- [ ] **Step 1: Full chr1 t8 parity** (drop `--skip-compare`).

```bash
cd /Users/mwiewior/research/git/vepyr/e2e-testing/scripts
env -u VIRTUAL_ENV -u CONDA_PREFIX uv run python run_annotation_fast.py chr1 --cache merged --forks 0 --threads 8 --force --backend lance > /tmp/csq_parity.out 2>/dev/null
grep -E "ALL .* CSQ fields match|match at 100" /tmp/csq_parity.out
```
Expected: `ALL 86 shared CSQ fields match at 100%!` (G1).

- [ ] **Step 2: t1–t8 throughput A/B** vs `f0dbb1b` (reuse the prior sweep method: stash any uncommitted instrumentation, build each commit non-dhat, run t1/t2/t4/t8 ×2, compare `Done: … variants/s`). Expected: no regression vs the `f0dbb1b` numbers (t1 12.6k, t2 28.5k, t4 44.8k, t8 63.8k var/s); ideally faster (G3).

- [ ] **Step 3: Peak RSS** at t8.

```bash
cd /Users/mwiewior/research/git/vepyr/e2e-testing/scripts
VEP_PROFILE=1 /usr/bin/time -l env -u VIRTUAL_ENV -u CONDA_PREFIX uv run python run_annotation_fast.py chr1 --cache merged --forks 0 --threads 8 --force --backend lance --skip-compare >/dev/null 2>/tmp/csq_rss.err
grep -iE "maximum resident set size|peak_rss" /tmp/csq_rss.err | tail
```

- [ ] **Step 4: dhat final** (chr22) total + attribution vs the 76.07M baseline (G2). Then rebuild non-dhat release to leave a clean build installed.

- [ ] **Step 5: Record G1–G4** in spec §8 + update memory `vep-peak-rss-allocation-churn` with the measured per-step block reductions and the final picture. Commit (docs only).

---

## Self-Review

Spec coverage: spec §3 steps 1–4 → Tasks 1–4; §4 protocol → Task 0 harness + each task's byte-diff/dhat steps + Task 5 final gate; §5 invariants → byte-diff gate + push_unique_value characterization test; §6 risks → R1/R2 (byte-diff), R3 (per-step decision in Task 2/3/4 Step 6), R4 (Task 2 case-by-case rule); §8 gates → Task 5. Placeholder scan: Step 1 has verbatim code; Steps 2–4 deliberately instruct read-then-transform (documented in File Structure note) because the sites span a 12k-line function and verbatim reproduction would risk fabrication — the byte-diff is the safety net. Type consistency: `push_unique_value(&mut Vec<String>, &str)` used consistently; scratch types named per task.

**Honest caveat:** Steps 2–4 give exact sites + rules, not full before/after code (the code isn't read verbatim in this plan). The executing subagent reads each site first. This mirrors the AF plan's spike approach and is the reason every step has a hard byte-diff gate. If a step's dhat delta is negligible (Step 6 decision), stop — do not complete all four on faith.
