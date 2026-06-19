# CSQ Output-Assembly Churn Reduction — Design Spec

**Date:** 2026-06-19
**Status:** Draft (pre-implementation)
**Component:** `datafusion/bio-function-vep` — per-variant CSQ output assembly
**Related memory:** `vep-peak-rss-allocation-churn`, `vep-annotation-bottlenecks`
**Predecessor:** `2026-06-19-af-arrow-zerocopy-design.md` (AF zero-copy — done; taught
the lesson this spec is built on: measure per-step, count ≠ peak ≠ wall).

## 1. Problem & evidence

dhat on chr22 (50,861 variants, single-thread, current HEAD `f0dbb1b`):
**76,070,619 allocation blocks** (~1,496 allocs/variant — matches the chr1 ~1,440
baseline). Per-call-stack attribution of those blocks:

| category | blocks | % of all allocs |
|---|---|---|
| `transcript_consequence` / `annotate_batch` CSQ-text materialize | 26.6M | **35.0%** |
| CSQ / VCF output assembly (`AnnotateProvider::annotate_batch`, `push_unique_value`, vcf_sink) | 19.0M | **25.0%** |
| transcript_overlap | 7.8M | 10.2% |
| lance_lookup | 6.3M | 8.3% |
| hgvsc | 4.3M | 5.7% |
| colocated_match (incl. the AF area already done) | 3.2M | 4.2% |
| other (few but byte-huge buffer allocs, 206 GB bytes) | 8.8M | 11.6% |

Independently, the time profile (single-thread chr1, engine-bound ~93%,
`lookup_wait` 0.30s) shows the same target dominating wall-time:
`transcript_output_materialize` 4.02s + `csq_format` 2.41s + `colocated_fields`
1.87s ≈ half of 14.5s engine time.

Two independent signals (allocation count + wall-time) converge: **the per-variant
CSQ output assembly is ~60% of allocations and ~half of engine time.** Contrast the
just-completed AF zero-copy work, which was only ~4% of allocs (`colocated_match`)
and moved nothing measurable — this target is ~15× larger and evidence-backed.

### Confirmed allocation sites (mapper, all `bio-function-vep/src/`)

- `annotate_provider.rs:1749` `push_unique_value(values: &mut Vec<String>, value: impl Into<String>)`
  — `value.into()` allocates a `String` on **every** call, including duplicates that
  are then dropped (~0.87M blocks). Call sites `:1897`, `:2027`, `:2045`.
- `annotate_provider.rs:5930–6252` — per-transcript CSQ-text loop. The output buffer
  `csq_buf` (`:5563`, 4096-cap) is already reused via `write!`; the churn is ~12
  throwaway owned-String temporaries built **per transcript** to feed each `write!`:
  `distance.map(|d| d.to_string())` (`:5973`), `format_hgvsp_output` (`:5981`),
  `bam_edit` uppercase (`:6001`), `refseq_offset.to_string()` (`:6014`),
  `tsl_str.to_string()` (`:6026`), `hgvs_offset.to_string()` (`:6071`),
  `appris_str` (`:6098`), `sift_str`/`polyphen_str` (`:6104`), `domains` (`:6132`),
  `mirna_str` (`:6150`).
- `annotate_provider.rs:6315–6900` — typed-column loop, a 2nd pass over the same
  transcripts that **recomputes** `bam_edit` (`:6334`), `refseq_offset` (`:6322`),
  and a fresh `terms.iter().collect::<Vec<&str>>().join("&")` (`:6341`) already built
  into `terms_buf` at `:5932`.
- `annotate_provider.rs:1955–2089` `ColocatedData::frequency_fields` (called per row
  `:5666`): `per_column: Vec<Vec<String>>` allocated per row (`:1961`); a
  `HashMap<String,String>` (`:1987`) + `HashSet<String>` (`:1988`) allocated **per
  (entry × AF column)**; many `.to_string()` + per-column `.join("&")` (`:2072`).
  `ColocatedData::variant_fields` (`:1858`) similar.

### Scope decision (confirmed with stakeholder)

**In-repo only, staged + measured.** The cross-crate VCF byte-serialization
(`vcf_sink::format_vcf_body_chunk` → sibling `datafusion_bio_format_vcf::serializer::
batch_to_vcf_lines` → per-batch `Vec<u8>`; the 435/232/134 MB byte sites) is the
peak-*bytes* lever but requires a sibling `datafusion-bio-formats` change + version
coordination. **Out of scope here**; recorded as a follow-up.

## 2. Goals / Non-Goals

**Goals**
- Reduce per-variant allocation churn in the in-repo CSQ assembly, **byte-identical**.
- Headline success = chr1 t8 throughput (variants/s) no worse, ideally better;
  secondary = dhat `Total` blocks ↓ (vs 76.07M chr22 baseline) and peak RSS
  (benefits indirectly via reduced fragmentation).
- Measure after each step; **stop when the measured ROI flattens** (the AF lesson).

**Non-Goals**
- Cross-crate VCF byte serialization / sibling `datafusion-bio-formats` changes.
- Changing CSQ field semantics, ordering, or the output RecordBatch schema/contract.
- The lance_lookup (8.3%) / transcript_overlap (10.2%) / hgvsc (5.7%) categories —
  not output-assembly; separate efforts if pursued.
- Merging the two transcript passes structurally (too risky for byte-identity — we
  *cache* across them instead).

## 3. Design — staged byte-identical steps (ordered by ROI ÷ risk)

Each step is an independent TDD commit gated by a fast byte-diff + a dhat re-measure.

**Step 1 — `push_unique_value(&str)`.** Change to `fn push_unique_value(values:
&mut Vec<String>, value: &str)`: linear dedup-scan on `&str` first, allocate
(`value.to_string()`) only on the insert path. Update call sites `:1897`/`:2027`/
`:2045` to pass `&str` (drop the `.clone()`/`.to_string()`). ~0.87M blocks; trivial.

**Step 2 — Per-transcript temporary elimination (CSQ-text loop `:5930–6252`).**
For each temporary feeding a `write!(csq_buf, …)`, write the underlying value
directly instead of materializing an owned `String`: integers via `write!`
formatting (no `.to_string()`); `&str` borrows where the source outlives the write;
keep `Cow` borrowed (don't force `.into_owned()`); uppercase (`bam_edit`) into a
small reused scratch `String` (`.clear()` per transcript) rather than a fresh alloc.
Functions returning owned `String` (sift/polyphen/domains/mirna/hgvsp) are evaluated
case-by-case: prefer a `write!`-into-buffer or `&str`/`Cow` return; only refactor the
ones that are hot and safe. The bulk of the 35%.

**Step 3 — De-duplicate the typed-column pass (`:6315+`).** Cache the per-transcript
values already computed in Step-2's loop (`bam_edit`, `refseq_offset`, the joined
`terms`) into a per-transcript scratch vector (`.clear()`-reused per row), and read
them in the typed-column pass instead of recomputing. Passes stay structurally
separate (byte-identity safety); only the redundant recomputation is removed.

**Step 4 — Reuse colocated scratch across rows (`frequency_fields`/`variant_fields`
`:1955+`).** Hoist `per_column: Vec<Vec<String>>` and the per-(entry×column)
`HashMap`/`HashSet` into reusable scratch owned by the per-batch context,
`.clear()`-ed per row/column instead of re-allocated. Pool the
`ColocatedFrequencyFields`/`ColocatedVariantFields` output structs (reuse + clear)
rather than constructing per row. The colocated 25%.

The map's owned-`String` interfaces are preserved; only allocation *lifetime* changes
(reuse instead of per-row birth/death). Byte output is unchanged.

## 4. Measurement protocol (per step)

1. **Byte-identical self-check (fast, the per-step gate):** annotate chr22
   (`--skip-compare`, non-dhat release) before vs after the step; `diff` the two
   output VCFs — MUST be identical. ~6s/run, no 23 GB reference needed.
2. **Churn delta:** dhat chr22 after the step (or each pair) → record `Total` blocks
   and `t-gmax` vs the 76.07M / 1.75 GB baseline.
3. **Decision:** if a step's measured block reduction is negligible, stop / re-scope
   rather than proceed to the next step on faith.

**Final gate (after the last step taken):**
- chr1 t8 e2e parity: 0 CSQ mismatch, 86/86 fields 100% (byte-identical).
- t1–t8 throughput A/B (before = `f0dbb1b`, after) — no regression, ideally faster.
- peak RSS (getrusage + `/usr/bin/time -l`) at t8.
- `cargo test -p datafusion-bio-function-vep` + clippy `-D warnings` green.

## 5. Invariants (byte-identical)

1. Every CSQ field's bytes, field order, and the `csq`/`most_severe_consequence`/
   typed-column outputs are unchanged — only allocation strategy changes.
2. Reused scratch buffers are `.clear()`-ed before each use; no stale data leaks.
3. `push_unique_value` dedup semantics (first-seen order, equality) unchanged.
4. Colocated map output values identical; only their allocation lifetime changes.

## 6. Risks

- **R1 — stale scratch.** Reused buffers retaining prior-row data. Mitigation:
  `.clear()` discipline; the per-step chr22 byte-diff catches any leak immediately.
- **R2 — cross-pass cache drift.** Step 3's cached temporaries diverging from the
  recompute they replace. Mitigation: byte-diff; keep the two passes separate.
- **R3 — diminishing returns.** Mitigation: measure-after-each; the protocol's
  decision point (§4.3) is the explicit stop signal.
- **R4 — `Cow`/borrow lifetime friction.** Some temporaries borrow data that doesn't
  outlive the `write!`. Mitigation: those keep a (reused-scratch) owned form; not
  every temporary must become zero-alloc — only the hot, safe ones.

## 7. Rollout

Internal change, public API + on-disk format + RecordBatch schema unchanged →
rebuild vepyr only, no sibling-crate coordination. Land each step behind its byte-diff
gate; full chr1 parity + throughput/RSS A/B at the end. Keep the `[VEP_RSS]` probes
for the before/after measurement.

## 8. Success criteria / gates (summary)

- **G1 Parity:** chr1 t8 — 0 CSQ mismatch, 86/86 fields 100%.
- **G2 Churn:** dhat chr22 `Total` blocks materially below 76.07M (target: remove a
  large fraction of the ~45M output-assembly blocks; exact target set by per-step
  measurement, not promised up front).
- **G3 Throughput:** t1–t8 A/B no regression vs `f0dbb1b`; ideally faster (engine-bound,
  so fewer allocs/variant should show).
- **G4 No regressions:** tests + clippy green.
