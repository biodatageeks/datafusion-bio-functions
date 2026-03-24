# Plan: Multithreaded VEP annotation with partitioned parquet cache

> Saved for future reference. This is the full parallel design to implement after
> the contig-by-contig single-threaded refactor is complete and benchmarked.

## Context

The current `scan_with_transcript_engine()` collects ALL VCF rows into memory, loads context tables (transcripts, exons, translations, etc.) for ALL chromosomes at once, builds a single `PreparedContext`, and annotates sequentially. For WGS VCFs (~5M+ variants), this means:
- Huge upfront memory for all context across all chromosomes
- Single-threaded annotation loop
- No benefit from the new per-chromosome parquet cache layout

The new parquet cache (`115_GRCh38_vep/`) is already partitioned per chromosome:
```
variation/chr1.parquet   (882 row groups, 88M rows for chr1)
transcript/chr1.parquet
exon/chr1.parquet
translation_core/chr1.parquet
translation_sift/chr1.parquet
regulatory/chr1.parquet
motif/chr1.parquet
```

**Goal**: Process VCF contig-by-contig, loading only that contig's context, then scanning variation for that contig with multiple DataFusion target partitions in parallel.

## Architecture

### Separate codepath, easy A/B benchmarking

Keep the existing `scan_with_transcript_engine()` completely untouched as the baseline. Add a new `scan_with_transcript_engine_partitioned()` alongside it with its own self-contained code (accepting some duplication for isolation).

**Runtime dispatch** in `AnnotateProvider::scan()`:
- Detect partitioned layout: `cache_source` is a directory containing `variation/` subdir with `chrN.parquet` files
- If partitioned AND `options_json` doesn't have `"partitioned": false` → new path
- Otherwise → existing monolithic path (unchanged baseline)
- Can force either path via `"partitioned": true/false` in options_json for benchmarking

This means both codepaths coexist, same binary, switchable per query. No risk to existing behavior.

### Processing model

```
VCF (position-sorted)
  → split into contig groups (batch-level)
  → for each contig:
      1. Register per-chrom parquet files:
         - variation/chrN.parquet  (with target_partitions = N for parallel scan)
         - transcript/chrN.parquet
         - exon/chrN.parquet
         - translation_core/chrN.parquet
         - translation_sift/chrN.parquet
         - regulatory/chrN.parquet
         - motif/chrN.parquet
      2. Load shared context ONCE: transcripts, exons, translations → PreparedContext for this contig
      3. Scan variation with N DataFusion partitions (row group parallelism)
      4. Each partition streams: variation batch → probe COITree → annotate → emit
         - All partitions share the contig's PreparedContext (read-only Arc)
         - SIFT/PolyPhen: per-partition cache (no lock contention)
      5. Collect annotated batches, proceed to next contig
  → concatenate all contig results
```

### Key changes

#### 1. `PartitionedParquetCache` struct (new file: `partitioned_cache.rs`)

```rust
struct PartitionedParquetCache {
    base_dir: PathBuf,
    /// Available chromosomes detected from variation/*.parquet
    available_chroms: Vec<String>,
}

impl PartitionedParquetCache {
    fn detect(path: &str) -> Option<Self>;  // checks dir layout
    fn variation_path(&self, chrom: &str) -> PathBuf;
    fn transcript_path(&self, chrom: &str) -> PathBuf;
    fn exon_path(&self, chrom: &str) -> PathBuf;
    // etc.
}
```

#### 2. VCF contig splitting (in `annotate_provider.rs`)

```rust
fn split_batches_by_contig(batches: &[RecordBatch]) -> Vec<(String, Vec<RecordBatch>)>
```

Since VCF is position-sorted, this walks batches and splits on contig transitions. Batches that span a contig boundary get sliced. Returns ordered `(chrom, batches)` pairs.

#### 3. `scan_with_transcript_engine_partitioned()` — new method on `AnnotateProvider`

The main contig-by-contig loop:

```rust
async fn scan_with_transcript_engine_partitioned(&self, ...) -> Result<Arc<dyn ExecutionPlan>> {
    // 1. Collect VCF (small — build side for COITree)
    let base_batches = collect(vcf_plan, task_ctx).await?;
    let contig_groups = split_batches_by_contig(&base_batches);

    let mut all_annotated = Vec::new();

    for (chrom, vcf_batches) in contig_groups {
        // 2. Register per-chrom parquet files as ephemeral tables
        let var_table = register_chrom_parquet(session, cache, "variation", &chrom).await?;
        let tx_table = register_chrom_parquet(session, cache, "transcript", &chrom).await?;
        let ex_table = register_chrom_parquet(session, cache, "exon", &chrom).await?;
        // ... translation_core, translation_sift, regulatory, motif

        // 3. Build contig context (shared across all partitions)
        let transcripts = load_transcripts(tx_table, &chrom).await?;
        let exons = load_exons(ex_table, &chrom).await?;
        let translations = load_translations(tl_table, &chrom).await?;
        let ctx = Arc::new(PreparedContext::new(&transcripts, &exons, ...));
        let coitree = Arc::new(build_coitree_for_contig(&vcf_batches));
        let engine = Arc::new(TranscriptConsequenceEngine::new(...));

        // 4. Create variation scan plan (N partitions from row groups)
        let var_scan = session.table(&var_table).await?.create_physical_plan().await?;

        // 5. Wrap with ContigAnnotateExec — fuses scan + lookup + annotation
        //    Each partition streams: variation batch → probe COITree → annotate → emit
        let annotate_exec = ContigAnnotateExec::new(
            var_scan,       // N partitions
            coitree,        // Arc-shared
            ctx,            // Arc-shared
            engine,         // Arc-shared
            colocated_map,  // Arc-shared
            sift_config,    // per-partition clone
        );

        // 6. Collect annotated results from all partitions
        let contig_batches = collect(Arc::new(annotate_exec), task_ctx).await?;
        all_annotated.extend(contig_batches);

        // 7. Deregister ephemeral tables to free memory
        deregister_chrom_tables(session, &chrom).await?;
    }

    MemTable::try_new(schema, vec![all_annotated])?.scan(state, None, &[], None).await
}
```

#### 4. Streaming annotation fused with variation scan partitions (`ContigAnnotateExec`)

Each DataFusion physical partition of the variation scan streams directly into annotation. No intermediate collect step for lookup results.

**`ContigAnnotateExec`** — new custom `ExecutionPlan`:
- Wraps the variation parquet scan plan (inherits its N partitions from row groups)
- Holds `Arc` shared state: COITree (VCF build side), `PreparedContext`, `TranscriptConsequenceEngine`, `ColocatedMap`
- `execute(partition)` returns a stream that:
  1. Pulls batches from the variation scan's partition
  2. Probes the COITree for VCF-cache matches (reuses `VariantLookupExec` join logic)
  3. Runs `annotate_batch_with_transcript_engine` on matched results
  4. Emits annotated batches
- Each partition gets its own `SiftPolyphenCache` (no lock contention)

```
variation/chr1.parquet (882 row groups)
  → DataFusion splits into N partitions (target_partitions config)
  → Partition 0: [rg 0..110]  ─stream─→ probe COITree ─→ annotate ─→ emit
  → Partition 1: [rg 111..220] ─stream─→ probe COITree ─→ annotate ─→ emit
  → ...
  → Partition 7: [rg 771..881] ─stream─→ probe COITree ─→ annotate ─→ emit
     All share: Arc<PreparedContext>, Arc<COITree>, Arc<Engine>
```

The VCF COITree build side is small and immutable — `Arc`-shared across all partitions. The `PreparedContext` (transcripts, exons, interval trees) is read-only after construction — also `Arc`-shared.

For SIFT/PolyPhen: each partition maintains its own `SiftPolyphenCache` and `SiftDirectReader` (or per-partition `SiftKvStore` handle). No cross-partition coordination needed since each partition covers a disjoint genomic range within the chromosome.

#### 5. Configuration

New option in `options_json`:
```json
{
  "target_partitions": 8,     // DataFusion target partitions per contig scan
  "partitioned_cache": true   // auto-detected, but can be forced
}
```

Also configurable via `SessionConfig::with_target_partitions()`.

## Files to modify

| File | Change |
|------|--------|
| `annotate_provider.rs` | Add `scan_with_transcript_engine_partitioned()`, `split_batches_by_contig()`, contig loop; dispatch from `scan()` when partitioned layout detected |
| `contig_annotate_exec.rs` | **New** — `ContigAnnotateExec` ExecutionPlan: fuses variation scan + COITree probe + annotation. N output partitions matching the variation scan. Each partition streams through annotation with `Arc`-shared context |
| `partitioned_cache.rs` | **New** — `PartitionedParquetCache` detection and path resolution |
| `annotate_table_function.rs` | No change — cache_source path passed through as-is |
| `annotation_store.rs` | No change — `AnnotationBackend::Parquet` handles both layouts |
| `lookup_provider.rs` | Extract COITree-building + probe logic into reusable functions (shared with `ContigAnnotateExec`) |
| `variant_lookup_exec.rs` | Extract COITree build + match logic into standalone fns callable from `ContigAnnotateExec` |
| `config.rs` | Add `target_partitions` config option |
| `lib.rs` | Register new modules |
| `miss_worklist.rs` | No change needed (already filters by chrom) |

## Duplication strategy

The following will be duplicated between the two codepaths (acceptable for benchmarking isolation):
- COITree building (extracted to shared helper, but each path calls it independently)
- Context loading (partitioned path loads from per-chrom files; monolithic loads from full table with WHERE clause)
- Batch annotation (both call `annotate_batch_with_transcript_engine`, but partitioned path wires it into `ContigAnnotateExec` streaming)

Shared without duplication:
- `annotate_batch_with_transcript_engine()` — same function, called from both paths
- `TranscriptConsequenceEngine`, `PreparedContext` — same data structures
- `SiftPolyphenCache` — same cache logic
- Allele matching, coordinate normalization — unchanged

## Prerequisites

Before implementing this plan, the contig-by-contig single-threaded refactor must be completed first:
1. Update the current parquet path to process contig-by-contig with lazy context loading
2. Benchmark the contig-by-contig approach vs monolithic baseline
3. Then layer on the multithreaded `ContigAnnotateExec` streaming partitions

## Verification

1. Unit test: `split_batches_by_contig` with multi-chrom batches and boundary splits
2. Integration test: `annotate_vep()` with partitioned parquet cache path, verify identical output to monolithic path
3. Benchmark: WGS VCF with `VEP_PROFILE=1`, compare timing of partitioned vs monolithic
4. Check memory: peak RSS should be ~1/24th of context tables (one chrom at a time)
