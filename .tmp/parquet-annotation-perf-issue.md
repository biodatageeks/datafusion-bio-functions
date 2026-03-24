# Parquet `annotate_vep`: avoid full-chrom context scans and wide materialization on cache misses

## Summary

The parquet-backed `annotate_vep()` path is still dominated by one-time parquet scanning/materialization work in the annotation phase, not by transcript consequence evaluation itself.

Current behavior on parquet:

- `lookup_variants()` scans the variation parquet and materializes lookup output.
- If **any** row is missing `cache_most_severe_consequence`, `annotate_vep()` loads **all** transcripts/exons/translations/regulatory/motif rows for every chromosome present in the input.
- Those context loaders use `SELECT *` and fully materialize the returned batches into Rust structs.
- The actual transcript consequence engine is comparatively cheap once context is in memory.

This issue is scoped to the **annotation-side parquet path**. It does **not** change the variation parquet lookup algorithm itself.

## Measured behavior

Local profile run:

```bash
VEP_PROFILE=1 cargo run -p datafusion-bio-function-vep --release --example profile_annotation -- \
  vep-benchmark/data/HG002_chr1.vcf.gz \
  /Users/mwiewior/research/data/vep/chr1-vep \
  1000
```

Observed timings:

- total: ~40.3s
- `variation_lookup (scan+collect)`: ~24.9s
- `context_tables_total`: ~15.1s
- `PreparedContext::new()`: ~26.6ms
- `sift_lazy_load + annotate_batches`: ~91.7ms

The files involved for chr1 are roughly:

- variation parquet: 2.8GB
- translation parquet: 229MB
- transcript parquet: 32MB
- exon parquet: 7.0MB
- regulatory parquet: 6.8MB

There is also a whole-chr1 benchmark run showing end-to-end runtime around 69s on another local path, and a slower real-world whole-chr1 runtime of about 120s has been observed as well. The exact number varies by machine/config, but the stage split above makes the current bottleneck ordering clear.

## Main bottlenecks in the annotation path

### 1. Any cache miss triggers whole-chromosome context loading

Code path:

- `datafusion/bio-function-vep/src/annotate_provider.rs`
- `scan_with_transcript_engine()`
- current miss gate around:
  - `has_cache_misses`
  - `load_transcripts`
  - `load_exons`
  - `load_translations`
  - `load_regulatory_features`
  - `load_motif_features`
  - `load_mirna_features`
  - `load_structural_features`

Today we do:

- inspect lookup output
- if any miss exists, collect `vcf_chroms`
- load all context rows for those chromosomes via `WHERE chrom IN (...)`

That means a small miss set can still load:

- all transcripts on chr1
- all exons on chr1
- all translations on chr1
- all regulatory rows on chr1

### 2. Context loaders use `SELECT *`

Examples:

- `load_transcripts()` uses `SELECT * FROM {table}{filter}`
- `load_exons()` uses `SELECT * FROM {table}{filter}`
- `load_translations()` uses `SELECT * FROM {table}{filter}`

This is especially expensive for translations. In default mode we do not need every column, but we still request every column and decode everything that survives the row filter.

### 3. We load exons/translations by chromosome, not by the actual transcript set we will evaluate

The engine only needs exons/translations for transcripts that overlap cache-miss variants. Today we first load all transcripts for the chromosome, but then we still load all exons and all translations for the same chromosome rather than semi-joining on transcript IDs.

### 4. Extra materialization/copying around the annotation step

The pipeline currently:

- collects lookup output into `base_batches`
- builds `annotated_batches`
- wraps those in a `MemTable`
- scans the `MemTable` again

That creates avoidable copies and memory pressure.

## Secondary observation

The current miss-rate accounting is not reliable for profiling because it uses `.any()` while mutating hit/miss counters, so the counters stop after the first batch that contains a miss. The boolean is fine; the reported totals are not.

This is not the core perf problem, but it makes the profile output misleading.

## Proposed code changes

### A. Replace chrom-only miss handling with a concrete miss worklist

Add a small helper struct in `annotate_provider.rs`, e.g.:

```rust
struct MissWorkItem {
    chrom: String,
    start: i64,
    end: i64,
}

struct MissWorklist {
    items: Vec<MissWorkItem>,
    by_chrom: HashMap<String, Vec<(i64, i64)>>,
}
```

Add a helper like:

```rust
fn collect_cache_miss_worklist(base_batches: &[RecordBatch]) -> Result<MissWorklist>
```

Behavior:

- scan all `base_batches` once
- collect rows where `cache_most_severe_consequence` is null
- normalize chromosome names
- merge nearby/overlapping intervals per chromosome to keep predicate lists small

This replaces the current `has_cache_misses` + `vcf_chroms` logic.

### B. Load transcripts/regulatory/motif/miRNA/SV by miss intervals, not whole chromosomes

Replace `chrom_filter_clause()` with interval-aware filtering support.

Proposed helper:

```rust
fn interval_filter_clauses(worklist: &MissWorklist) -> Vec<String>
```

For each chromosome, emit one or more clauses like:

```sql
(chrom = '1' AND start <= {max_end} AND "end" >= {min_start})
```

Then query:

```sql
SELECT <projected cols>
FROM `{table}`
WHERE <interval disjunction>
```

Apply this to:

- `load_transcripts`
- `load_regulatory_features`
- `load_motif_features`
- `load_mirna_features`
- `load_structural_features`

Important: do not use chromosome-only predicates once a miss worklist exists.

### C. Load exons and translations by `transcript_id`

After loading transcripts, derive:

```rust
let transcript_ids: HashSet<String> = tx.iter().map(|t| t.transcript_id.clone()).collect();
```

Then introduce:

```rust
async fn load_exons_for_transcript_ids(...)
async fn load_translations_for_transcript_ids(...)
```

These should query by `transcript_id` instead of chromosome:

```sql
SELECT <projected cols>
FROM `{table}`
WHERE transcript_id IN (...)
```

If the `IN (...)` list is too large, batch it into chunks or register a small in-memory table of transcript IDs and do a join.

This is the highest-value annotation-side change after miss scoping.

### D. Make projections mode-aware instead of `SELECT *`

Introduce explicit column lists.

Suggested minimum sets:

#### transcripts: default

- `transcript_id`
- `chrom`
- `start`
- `end`
- `strand`
- `biotype`
- `cds_start`
- `cds_end`
- `cdna_coding_start`
- `cdna_coding_end`
- `gene_stable_id`
- `gene_symbol`
- `gene_symbol_source`
- `gene_hgnc_id`
- `source`
- `version`
- `cds_start_nf`
- `cds_end_nf`
- `flags_str`
- `cdna_mapper_segments`
- batch-1 metadata currently emitted in CSQ

Only load these when required:

- `spliced_seq`
- `cdna_seq`
- `translateable_seq`
- `bam_edit_status`
- `has_non_polya_rna_edit`
- `mature_mirna_regions`

#### translations: default

- `transcript_id`
- `cds_len`
- `protein_len`
- `translation_seq`
- `cds_sequence`
- `stable_id`
- `version`

Load `protein_features` only when the output path needs `DOMAINS`.

Do not load:

- `sift_predictions`
- `polyphen_predictions`

in the main translation loader; keep those in the existing windowed loader.

#### regulatory / motif / miRNA / SV

Project only the columns actually read in the corresponding loader.

### E. Remove the final `MemTable` round-trip

Current code:

- builds `projected_batches`
- creates `MemTable`
- calls `mem.scan(...)`

Instead:

- return a `MemoryExec` directly if we still need batch materialization
- or, preferably, introduce a small custom execution plan / stream that yields annotated batches directly

This is lower impact than A-D but still worthwhile once the bigger scans are reduced.

### F. Fix miss-rate profiling output while touching this path

Do not use `.any()` while accumulating counters.

Instead:

- compute the full miss worklist once
- derive `has_cache_misses = !worklist.items.is_empty()`
- log real hit/miss counts from the full scan

## Suggested implementation order

1. Introduce `MissWorklist` and replace the current miss gate.
2. Make transcript/regulatory/motif/miRNA/SV loaders interval-aware.
3. Load exons/translations by transcript IDs.
4. Replace `SELECT *` with explicit projection lists.
5. Remove the final `MemTable` scan.
6. Fix profiling counters.

## Expected impact

On the profiled parquet annotation path, the realistic target for this issue is roughly:

- save most of the current `context_tables_total` cost
- reduce translation load time substantially
- cut end-to-end runtime by about 10-20s on current local chr1 runs

Relative speedup depends on workload:

- sample / localized miss workloads: larger gain
- full chr1 workloads with broadly distributed misses: smaller gain

But in both cases, this should be the highest-value work item **within the annotation path** before revisiting the variation parquet lookup algorithm itself.

## Files expected to change

- `datafusion/bio-function-vep/src/annotate_provider.rs`
- possibly `datafusion/bio-function-vep/src/transcript_consequence.rs` only if loader output contracts need small adjustments
- maybe a new helper module if the miss-worklist / projection logic grows too large

## Acceptance criteria

- `annotate_vep()` remains behaviorally identical on existing tests.
- Default parquet annotation no longer loads whole-chromosome exons/translations for a tiny miss set.
- Translation loader no longer uses `SELECT *`.
- Profiling output reports correct miss counts.
- Re-profile chr1 parquet runs and record before/after numbers for:
  - `variation_lookup`
  - `context_tables_total`
  - `load_translations`
  - total runtime
