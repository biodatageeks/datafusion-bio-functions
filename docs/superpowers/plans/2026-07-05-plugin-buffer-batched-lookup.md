# Plugin Buffer-Batched Lookup + CSQ Integration — Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Replace the plugin cache's whole-file in-memory lookup with a buffer-batched, page-scoped PageDir take (reusing the variation lookup's primitives without touching it), then wire plugin CSQ fields into VEP output — validated against the AlphaMissense chr22 golden.

**Architecture:** Per annotation buffer, each plugin does **one** PageDir take reading only the candidate pages for that buffer's positions (`take_buffer`), producing a compact `PluginBufferSlice`. The transcript engine then does synchronous per-transcript `probe`s (filter by `allele_string` + per-transcript discriminator) into that slice. Built on the shared `parquet_cache::page_dir` primitives; `SinglePathParquetVariationLookup` is untouched.

**Tech Stack:** Rust (edition 2024), DataFusion 53 / Arrow 58, `parquet` crate, existing `parquet_cache::page_dir`, the `plugin_cache` module (Tasks 1–9 already committed).

**Spec:** `docs/superpowers/specs/2026-07-05-custom-vep-plugin-caches-design.md` §5.

**Starting point:** `feat/plugin-cache-alphamissense` at the committed option-B checkpoint (in-memory `PluginLookup` green, 16 tests). The chr22 cache is built at `/tmp/plugin_cache/plugin/alphamissense/chr22.parquet` (1,481,489 rows, has `protein_variant`). Golden at `/Users/mwiewior/research/git/vepyr/sandbox/HG002_chr22_everything_hgvs_merged_am.vcf`; input `HG002_chr22.vcf`.

## Global Constraints

- DataFusion 53.0.0 / Arrow 58.0.0 / Rust edition 2024 — do not bump.
- **Do not modify** `src/parquet_cache/variation_lookup.rs` (`SinglePathParquetVariationLookup`) — parity-critical. Reuse only `parquet_cache::page_dir::{PageDir, CoalescingAsyncReader, IoCounters, selection_from_ranges, selection_from_offsets}`.
- Reads are **once per buffer per plugin** (`take_buffer`), page-scoped — never load a whole shard. Per-transcript `probe` is a sync in-memory filter (no Rayon/async pool).
- Lookup key: `(start:u32, allele_string, <match discriminator values…>)`. Miss → empty (the per-transcript gate).
- Parity is the gate: plugin CSQ fields byte-identical to VEP on chr22; w1-vs-w4 body byte-identical; the width-alignment test (`annotate_provider.rs:14961`) stays green.
- `cargo fmt` + `cargo clippy -- -D warnings` pass before every commit. Commit author is the repo's configured identity (now `Marek Wiewiorka`).

## File Structure

- `src/plugin_cache/lookup.rs` — **rewrite** in-memory map → `PluginLookup` (PageDir + `take_buffer`) and `PluginBufferSlice` (+ `PluginScalar`, `probe`).
- `src/plugin_cache/registry.rs` — **rewrite** `probe_all` → `take_buffer_all` returning per-plugin slices; keep `csq_fields`, `EngineAttrs`.
- `src/plugin_cache/build.rs`, `write.rs`, `source_manifest.rs`, `cache_manifest.rs`, `normalize.rs`, `provider.rs` — unchanged (option-B data side already done).
- `src/annotate_provider.rs` — **modify** (Task I1): per-buffer `take_buffer_all`, per-transcript probe, append CSQ fields (3 emission paths + placeholder layout).
- `src/golden_benchmark.rs`, `src/vcf_sink.rs` — **modify** (Task I1): append plugin field names to the CSQ header/format list.

---

## Task R1: `PluginLookup` — PageDir open + `take_buffer`

Rewrite the in-memory `PluginLookup` into the batched, page-scoped form.

**Files:**
- Rewrite: `src/plugin_cache/lookup.rs`

**Interfaces:**
- Consumes: `parquet_cache::page_dir::{PageDir, CoalescingAsyncReader, IoCounters, selection_from_ranges, selection_from_offsets}`.
- Produces: `struct PluginLookup { path, meta, page_dir, start_leaf, projection: Vec<String> }`; `async fn open(shard: &Path, match_columns: Vec<String>, value_columns: Vec<String>) -> Result<PluginLookup>`; `async fn take_buffer(&self, sorted_unique_starts: &[u32]) -> Result<RecordBatch>` returning columns `[start, allele_string, <match cols…>, <value cols…>]`; `fn projection(&self) -> &[String]`.

- [ ] **Step 1: Write the failing test** in `lookup.rs` (build a real tiered shard with `PluginShardWriter`, then `take_buffer`):

```rust
#[cfg(test)]
mod tests {
    use super::*;
    use crate::plugin_cache::source_manifest::{MatchColumn, ValueColumn, ValueType};
    use crate::plugin_cache::write::{PluginShardWriter, plugin_output_schema};
    use datafusion::arrow::array::{Float32Array, Int8Array, StringArray, UInt32Array};
    use datafusion::arrow::record_batch::RecordBatch;
    use std::sync::Arc;

    fn write_shard(path: &std::path::Path) {
        // three positions; 100 has two protein_variants (multi-isoform)
        let matches = vec![MatchColumn { column: "protein_variant".into(), engine_attr: "amino_acid_change".into() }];
        let vals = vec![ValueColumn { column: "am".into(), csq_field: "am".into(), ty: ValueType::Float32 }];
        let schema = plugin_output_schema(&matches, &vals);
        let batch = RecordBatch::try_new(schema.clone(), vec![
            Arc::new(StringArray::from(vec!["22","22","22"])),
            Arc::new(UInt32Array::from(vec![100u32,100,200])),
            Arc::new(UInt32Array::from(vec![100u32,100,200])),
            Arc::new(StringArray::from(vec!["A/G","A/G","C/T"])),
            Arc::new(StringArray::from(vec!["R12G","R78G","S9L"])),
            Arc::new(Float32Array::from(vec![0.0392f32,0.0427,0.5])),
            Arc::new(Int8Array::from(vec![1i8,1,1])),
        ]).unwrap();
        let mut w = PluginShardWriter::create(path, schema).unwrap();
        w.write(&batch).unwrap();
        w.finish().unwrap();
    }

    #[tokio::test(flavor = "multi_thread")]
    async fn take_buffer_reads_candidate_rows() {
        let dir = tempfile::tempdir().unwrap();
        let path = dir.path().join("chr22.parquet");
        write_shard(&path);
        let lk = PluginLookup::open(&path, vec!["protein_variant".into()], vec!["am".into()]).await.unwrap();
        // buffer requests only position 100 → both isoform rows returned, 200 excluded
        let batch = lk.take_buffer(&[100]).await.unwrap();
        assert_eq!(batch.num_rows(), 2);
        let names: Vec<_> = batch.schema().fields().iter().map(|f| f.name().clone()).collect();
        assert_eq!(names, vec!["start","allele_string","protein_variant","am"]);
        // empty request → empty batch
        assert_eq!(lk.take_buffer(&[]).await.unwrap().num_rows(), 0);
    }
}
```

- [ ] **Step 2: Run test to verify it fails**

Run: `cargo test -p datafusion-bio-function-vep plugin_cache::lookup::tests::take_buffer_reads_candidate_rows`
Expected: FAIL — `open`/`take_buffer` signatures changed.

- [ ] **Step 3: Write the implementation.** Replace the whole body of `lookup.rs` above the tests. Model `take_buffer` on `SinglePathParquetVariationLookup::{open, exact_offsets, take_payload}` (read `src/parquet_cache/variation_lookup.rs:56-240` first), generalized to the plugin projection:

```rust
//! Buffer-batched, page-scoped plugin lookup. Per buffer, one PageDir take reads
//! only the candidate pages for the buffer's positions (never the whole shard);
//! per-transcript matching is a sync filter over the resulting PluginBufferSlice.
use std::collections::HashMap;
use std::path::{Path, PathBuf};

use datafusion::arrow::array::{Array, Float32Array, Int32Array, LargeStringArray, StringArray, StringViewArray, UInt32Array};
use datafusion::arrow::record_batch::RecordBatch;
use datafusion::common::{DataFusionError, Result};
use futures::TryStreamExt;
use parquet::arrow::ProjectionMask;
use parquet::arrow::arrow_reader::{ArrowReaderMetadata, ArrowReaderOptions};
use parquet::arrow::async_reader::ParquetRecordBatchStreamBuilder;
use std::collections::HashSet;

use crate::parquet_cache::page_dir::{
    CoalescingAsyncReader, IoCounters, PageDir, selection_from_offsets, selection_from_ranges,
};

const COALESCE_GAP_BYTES: u64 = 512 * 1024;

#[derive(Debug, Clone, PartialEq)]
pub enum PluginScalar { Str(String), F32(f32), I32(i32), Null }

pub struct PluginLookup {
    path: PathBuf,
    meta: ArrowReaderMetadata,
    page_dir: PageDir,
    start_leaf: usize,
    projection: Vec<String>, // allele_string, <match…>, <value…> (start added at take)
}

impl PluginLookup {
    pub async fn open(shard: &Path, match_columns: Vec<String>, value_columns: Vec<String>) -> Result<Self> {
        let file = std::fs::File::open(shard).map_err(|e| DataFusionError::Execution(format!("open plugin shard '{}': {e}", shard.display())))?;
        #[allow(deprecated)]
        let meta = ArrowReaderMetadata::load(&file, ArrowReaderOptions::new().with_page_index(true))
            .map_err(|e| DataFusionError::Execution(format!("load plugin metadata: {e}")))?;
        let start_leaf = meta.parquet_schema().columns().iter().position(|c| c.name() == "start")
            .ok_or_else(|| DataFusionError::Execution("plugin shard has no 'start' column".into()))?;
        let page_dir = PageDir::build(&meta, start_leaf)?;
        let mut projection = vec!["allele_string".to_string()];
        projection.extend(match_columns);
        projection.extend(value_columns);
        Ok(Self { path: shard.to_path_buf(), meta, page_dir, start_leaf, projection })
    }

    pub fn projection(&self) -> &[String] { &self.projection }

    async fn open_async(&self) -> Result<tokio::fs::File> {
        tokio::fs::File::open(&self.path).await.map_err(|e| DataFusionError::Execution(format!("reopen plugin shard '{}': {e}", self.path.display())))
    }

    pub async fn take_buffer(&self, sorted_unique_starts: &[u32]) -> Result<RecordBatch> {
        let arrow_schema = self.meta.schema();
        // projected payload schema: start + projection columns, in a stable order
        let mut cols: Vec<String> = vec!["start".to_string()];
        cols.extend(self.projection.iter().cloned());
        let root_indices: Vec<usize> = cols.iter().filter_map(|n| arrow_schema.index_of(n).ok()).collect();
        if sorted_unique_starts.is_empty() {
            let mask = ProjectionMask::roots(self.meta.parquet_schema(), root_indices);
            // build an empty stream just to get the projected schema
            let builder = ParquetRecordBatchStreamBuilder::new_with_metadata(
                CoalescingAsyncReader::new(self.open_async().await?, IoCounters::new(), COALESCE_GAP_BYTES),
                self.meta.clone()).with_projection(mask);
            return Ok(RecordBatch::new_empty(builder.schema().clone()));
        }
        let counters = IoCounters::new();
        let probe_set: HashSet<u32> = sorted_unique_starts.iter().copied().collect();
        let probes64: Vec<u64> = sorted_unique_starts.iter().map(|&s| s as u64).collect();
        let ranges = self.page_dir.resolve_ranges(&probes64);

        // phase 2: start-only read → exact offsets
        let start_mask = ProjectionMask::leaves(self.meta.parquet_schema(), [self.start_leaf]);
        let mut sstream = ParquetRecordBatchStreamBuilder::new_with_metadata(
            CoalescingAsyncReader::new(self.open_async().await?, counters.clone(), COALESCE_GAP_BYTES), self.meta.clone())
            .with_projection(start_mask).with_row_selection(selection_from_ranges(&ranges)).with_batch_size(8192)
            .build().map_err(|e| DataFusionError::Execution(format!("build start stream: {e}")))?;
        let mut offsets = Vec::new();
        let mut cur = ranges.iter().flat_map(|&(a,b)| a..b);
        while let Some(b) = sstream.try_next().await.map_err(|e| DataFusionError::Execution(format!("read start batch: {e}")))? {
            let sa = b.column(0).as_any().downcast_ref::<UInt32Array>().ok_or_else(|| DataFusionError::Execution("start not UInt32".into()))?;
            for &v in sa.values() {
                let off = cur.next().ok_or_else(|| DataFusionError::Execution("range cursor underflow".into()))?;
                if probe_set.contains(&v) { offsets.push(off); }
            }
        }

        // phase 3: projected payload take at exact offsets
        let mask = ProjectionMask::roots(self.meta.parquet_schema(), root_indices);
        let builder = ParquetRecordBatchStreamBuilder::new_with_metadata(
            CoalescingAsyncReader::new(self.open_async().await?, counters.clone(), COALESCE_GAP_BYTES), self.meta.clone())
            .with_projection(mask).with_row_selection(selection_from_offsets(&offsets)).with_batch_size(8192);
        let proj_schema = builder.schema().clone();
        let mut pstream = builder.build().map_err(|e| DataFusionError::Execution(format!("build payload stream: {e}")))?;
        let mut taken = Vec::new();
        while let Some(b) = pstream.try_next().await.map_err(|e| DataFusionError::Execution(format!("read payload: {e}")))? { taken.push(b); }
        if taken.is_empty() { return Ok(RecordBatch::new_empty(proj_schema)); }
        datafusion::arrow::compute::concat_batches(&proj_schema, &taken).map_err(|e| DataFusionError::Execution(format!("concat payload: {e}")))
    }
}
```

(Keep the `string_value` / `decode_scalar` helpers from the current file — Task R2 uses them.)

- [ ] **Step 4: Run test to verify it passes**

Run: `cargo test -p datafusion-bio-function-vep plugin_cache::lookup::tests::take_buffer_reads_candidate_rows`
Expected: PASS (2 rows for position 100; empty for `[]`).

- [ ] **Step 5: Commit** (registry won't compile yet — that's Task R3; commit lookup + a temporary `#[allow(dead_code)]` on helpers if clippy complains, or fold R1+R2+R3 into one commit at the end of R3). Prefer: do R1→R3 then commit once.

```bash
cargo fmt && git add datafusion/bio-function-vep/src/plugin_cache/lookup.rs
# (deferred commit — see R3)
```

---

## Task R2: `PluginBufferSlice` — per-buffer working set + `probe`

**Files:**
- Modify: `src/plugin_cache/lookup.rs` (add `PluginBufferSlice`)

**Interfaces:**
- Consumes: the `take_buffer` batch (`[start, allele_string, <match…>, <value…>]`).
- Produces: `struct PluginBufferSlice { rows: HashMap<(u32,String,Vec<Option<String>>), Vec<PluginScalar>>, n_match: usize, n_values: usize }`; `fn from_batch(batch: &RecordBatch, n_match: usize, n_values: usize) -> Result<PluginBufferSlice>`; `fn probe(&self, start: u32, allele_string: &str, match_values: &[Option<String>]) -> Option<Vec<PluginScalar>>`.

- [ ] **Step 1: Write the failing test** in `lookup.rs`:

```rust
#[tokio::test(flavor = "multi_thread")]
async fn slice_probe_gates_on_discriminator() {
    let dir = tempfile::tempdir().unwrap();
    let path = dir.path().join("chr22.parquet");
    write_shard(&path); // from Task R1 test module
    let lk = PluginLookup::open(&path, vec!["protein_variant".into()], vec!["am".into()]).await.unwrap();
    let batch = lk.take_buffer(&[100]).await.unwrap();
    let slice = PluginBufferSlice::from_batch(&batch, 1, 1).unwrap();
    // isoform-specific hits
    assert!(matches!(slice.probe(100,"A/G",&[Some("R78G".into())]).unwrap()[0], PluginScalar::F32(v) if (v-0.0427).abs()<1e-6));
    assert!(matches!(slice.probe(100,"A/G",&[Some("R12G".into())]).unwrap()[0], PluginScalar::F32(v) if (v-0.0392).abs()<1e-6));
    // gate: no aa-change (non-missense) → miss
    assert!(slice.probe(100,"A/G",&[None]).is_none());
    // wrong allele → miss
    assert!(slice.probe(100,"C/T",&[Some("R12G".into())]).is_none());
}
```

- [ ] **Step 2: Run to verify it fails.** Run: `cargo test -p datafusion-bio-function-vep plugin_cache::lookup::tests::slice_probe_gates_on_discriminator` → FAIL.

- [ ] **Step 3: Implement `PluginBufferSlice`** in `lookup.rs` (batch columns: 0=start, 1=allele_string, 2..2+n_match=match cols, 2+n_match..=value cols):

```rust
pub struct PluginBufferSlice {
    rows: HashMap<(u32, String, Vec<Option<String>>), Vec<PluginScalar>>,
}

impl PluginBufferSlice {
    pub fn from_batch(batch: &RecordBatch, n_match: usize, n_values: usize) -> Result<Self> {
        let mut rows = HashMap::new();
        if batch.num_rows() == 0 { return Ok(Self { rows }); }
        let start = batch.column(0).as_any().downcast_ref::<UInt32Array>()
            .ok_or_else(|| DataFusionError::Execution("start not UInt32".into()))?;
        for r in 0..batch.num_rows() {
            let allele = string_value(batch.column(1).as_ref(), r)?
                .ok_or_else(|| DataFusionError::Execution("null allele_string".into()))?;
            let mut mv = Vec::with_capacity(n_match);
            for c in 0..n_match { mv.push(string_value(batch.column(2 + c).as_ref(), r)?); }
            let mut vals = Vec::with_capacity(n_values);
            for c in 0..n_values { vals.push(decode_scalar(batch.column(2 + n_match + c).as_ref(), r)?); }
            rows.insert((start.value(r), allele, mv), vals);
        }
        Ok(Self { rows })
    }

    pub fn probe(&self, start: u32, allele_string: &str, match_values: &[Option<String>]) -> Option<Vec<PluginScalar>> {
        self.rows.get(&(start, allele_string.to_string(), match_values.to_vec())).cloned()
    }
}
```

- [ ] **Step 4: Run to verify it passes.** Same command → PASS.

- [ ] **Step 5: Commit** — deferred to R3.

---

## Task R3: `PluginRegistry` — `take_buffer_all` + per-slice probe

**Files:**
- Rewrite: `src/plugin_cache/registry.rs`

**Interfaces:**
- Consumes: `PluginLookup`, `PluginBufferSlice`, `PluginScalar` (R1/R2), `EngineAttrs`, `discover_plugins`, `canonical_chrom_label`.
- Produces: `PluginRegistry::open(cache_root, chrom)` (unchanged surface, now stores `n_match`/`n_values` per plugin); `async fn take_buffer_all(&self, sorted_unique_starts: &[u32]) -> Result<BufferSlices>`; `struct BufferSlices` with `fn probe_all(&self, start: u32, allele_string: &str, attrs: &EngineAttrs) -> Vec<PluginScalar>`; `fn csq_fields(&self) -> Vec<String>`; `fn is_empty(&self)`.

- [ ] **Step 1: Write the failing test** — build the chr22-style shard + manifest (reuse the current registry test fixture), then:

```rust
let reg = PluginRegistry::open(cache_root, "22").await.unwrap();
let slices = reg.take_buffer_all(&[100]).await.unwrap();
let attrs = EngineAttrs { amino_acid_change: Some("R78G".into()) };
let hit = slices.probe_all(100, "A/G", &attrs);
assert!(matches!(hit[0], PluginScalar::F32(v) if (v-0.0427).abs()<1e-6));
// non-missense (no aa-change) → Null (gate)
let none = slices.probe_all(100, "A/G", &EngineAttrs::default());
assert_eq!(none, vec![PluginScalar::Null]);
```

- [ ] **Step 2: Run to verify it fails.** Run: `cargo test -p datafusion-bio-function-vep plugin_cache::registry` → FAIL.

- [ ] **Step 3: Implement.** `PluginEntry` keeps `name, csq_fields, match_engine_attrs, n_match, n_values, lookup: Option<PluginLookup>`. `take_buffer_all` calls `lookup.take_buffer(starts)` for each plugin (→ `PluginBufferSlice::from_batch`), collecting `Vec<(entry_meta, Option<PluginBufferSlice>)>` into `BufferSlices`. `BufferSlices::probe_all` mirrors the old `probe_all`: per plugin, build `match_values` from `attrs.get(engine_attr)`, `slice.probe(...)`, else `Null`-pad `csq_fields.len()`. Keep `EngineAttrs` and `csq_fields` as-is.

- [ ] **Step 4: Run to verify it passes.** Same command → PASS. Then the whole module: `cargo test -p datafusion-bio-function-vep plugin_cache` → all green.

- [ ] **Step 5: Commit R1+R2+R3 together**

```bash
cargo fmt && cargo clippy -p datafusion-bio-function-vep -- -D warnings
git add datafusion/bio-function-vep/src/plugin_cache/lookup.rs datafusion/bio-function-vep/src/plugin_cache/registry.rs
git commit -m "feat(plugin-cache): buffer-batched page-scoped lookup (take_buffer + BufferSlice)"
```

---

## Task I1: Wire plugin CSQ fields into annotation output

The 6-site integration (see spec §5.2 and the map in `docs/superpowers/plans/2026-07-05-alphamissense-plugin-prototype.md`). Do it behind an `Option<PluginRegistry>` so a no-plugin run is byte-identical to today.

**Files:**
- Modify: `src/annotate_provider.rs` — per-buffer `take_buffer_all`; per-transcript `probe_all`; append fields in the 3 emission paths (`~5627`, `~6147-6153`, `~6170`) + `CsqPlaceholderLayout` (`~1489-1584`); open the registry in `prepare_contig_context` (`~12120`) onto `SharedContigAnnotationContext` (`~8905`); thread to `annotate_batch_with_transcript_engine` (`~5234`, call at `~11222`).
- Modify: `src/golden_benchmark.rs` (`csq_field_names_for_mode_with_pick`, `~661`) and `src/vcf_sink.rs` (`csq_header_description`, `~637`) — append `registry.csq_fields()` (concatenated at the join sites; the names are runtime `String`s).

**Interfaces:**
- Consumes: `PluginRegistry::{open, take_buffer_all, csq_fields}`, `BufferSlices::probe_all`, `EngineAttrs`, `PluginScalar`.
- Produces: plugin CSQ columns in the emitted `##INFO CSQ Format` header and each transcript's CSQ string.

- [ ] **Step 1: Write the failing integration test** (`src/annotate_provider.rs` tests or `tests/`): a 1-variant chr22 VCF at a known AlphaMissense missense position with ≥2 transcripts (one missense, one non-coding); a pre-built plugin shard on that `(start, allele, protein_variant)`; assert the emitted CSQ has `am_pathogenicity`/`am_class` **only** on the matching-discriminator transcript line and empty on the other; and a variant with no shard entry emits empty. Keep the fixture tiny (reuse `PluginShardWriter`).

- [ ] **Step 2: Run to confirm it fails** (fields absent / mis-gated). Run: `cargo test -p datafusion-bio-function-vep plugin_csq_injection -- --nocapture` → FAIL.

- [ ] **Step 3: Implement the wiring.**
  1. **Field names:** append `registry.csq_fields()` after the existing trailing fields at the `golden_benchmark.rs` join site (used by both header and layout), and in `vcf_sink.rs::csq_header_description`.
  2. **Per buffer:** at the top of `annotate_batch_with_transcript_engine`, if `Some(registry)`, collect the batch's unique `start`s (from the input columns) and `let slices = block_on(registry.take_buffer_all(&starts))?` — in the same block_on-valid seam variation's cold probe uses.
  3. **Per transcript:** where the CSQ line is built, compute `EngineAttrs { amino_acid_change: aa_change_from(amino_acids, protein_position) }` (format `{ref}{pos}{alt}` from `Amino_acids` `"V/A"` + `Protein_position`; `None` if absent), then `let pv = slices.probe_all(start, input_allele_string, &attrs)`; format each `PluginScalar` (F32 → VEP number format; Str → verbatim; Null/`Null` → empty) and append in field order across the 3 paths + `CsqPlaceholderLayout::append_entry`.
  4. **Contig setup:** in `prepare_contig_context`, `let plugin_registry = match &config.cache_root { Some(root) => Some(PluginRegistry::open(root, &chrom).await?), None => None };` (only when non-empty), store on `SharedContigAnnotationContext`, thread to the worker call.

- [ ] **Step 4: Run tests.** Integration test → PASS. Full suite incl. the width-alignment test: `cargo test -p datafusion-bio-function-vep` → PASS (no-plugin paths unchanged).

- [ ] **Step 5: Commit**

```bash
cargo fmt && cargo clippy -p datafusion-bio-function-vep -- -D warnings
git add datafusion/bio-function-vep/src/annotate_provider.rs datafusion/bio-function-vep/src/golden_benchmark.rs datafusion/bio-function-vep/src/vcf_sink.rs
git commit -m "feat(plugin-cache): inject plugin CSQ fields per transcript (buffer-batched probe)"
```

---

## Task V1: AlphaMissense chr22 golden parity gate

**Files:**
- Create: `datafusion/bio-function-vep/scripts/parity_alphamissense_chr22.sh` (+ a comparison harness or reuse the golden-VEP comparator).

**Interfaces:** Consumes the chr22 plugin cache (`/tmp/plugin_cache`) + the wiring (I1) + the golden VCF.

- [ ] **Step 1: Annotate** `HG002_chr22.vcf` with vepyr, AlphaMissense enabled (`--cache_root /tmp/plugin_cache` or the vepyr flag that points at the plugin cache), emitting `am_pathogenicity`/`am_class`.

- [ ] **Step 2: Compare** per `(chrom,pos,ref,alt,Feature)` against `HG002_chr22_everything_hgvs_merged_am.vcf`, restricted to the two `am_*` fields. Expected: **1,912 populated lines match; the 3,905 empty lines stay empty** (the gate). Record any mismatch like the existing summary reports.

- [ ] **Step 3: w1-vs-w4 body byte-identity** on the plugin-enabled annotation. Expected: identical bodies.

- [ ] **Step 4: Commit** the parity script + a short results note.

```bash
git add datafusion/bio-function-vep/scripts/parity_alphamissense_chr22.sh
git commit -m "test(plugin-cache): AlphaMissense chr22 golden parity + worker byte-identity"
```

---

## Self-Review Notes

- **Spec coverage:** §5.1 components (R1 `PluginLookup`/`take_buffer`, R2 `PluginBufferSlice`, R3 `take_buffer_all`); §5.2 flow (I1); §5.3 generic per-variant/per-transcript (R2/R3 empty-discriminator path — covered by the existing `probe(...,&[])` per-variant tests); §5.4 concurrency/memory (R1 page-scoped, I1 block_on seam); §3.4 discriminator (R2 gate test); §8 parity (V1).
- **Untouched:** `variation_lookup.rs` (Global Constraint); build/data-side modules (option-B already done).
- **Executor unknowns flagged inline:** the exact vepyr flag to point at the plugin cache (V1 Step 1), the `Amino_acids`→`amino_acid_change` formatter for multi-residue edge cases (I1 Step 3 — SNV missense is single-residue), and the precise block_on seam in `annotate_batch_with_transcript_engine` (I1 Step 2, mirror `lookup_exec.rs:1676`).
- **Prototype pragmatism:** V1 uses the already-built `/tmp/plugin_cache` chr22 cache; a production build wires the full-genome AlphaMissense_hg38 via the driver.
