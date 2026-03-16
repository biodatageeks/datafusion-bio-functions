# Implementation Tasks

## Status Key
- [x] Completed and merged
- [ ] Not started

## 0. Phase 0: Fix Profiling & MissWorklist Foundation — DONE

- [x] 0.1 Fix miss-rate profiling counters (replace `.any()` short-circuit with full iteration)
- [x] 0.2 Implement MissWorklist (collect_miss_worklist, interval coalescing, chrom normalization)
- [x] 0.3 Implement interval filter SQL builder (interval_filter_sql with fallback for >50 intervals)

## 1. Phase 1: Interval-Aware Context Loading — DONE

- [x] 1.1 Update all 7 loaders from `chroms: &HashSet<String>` to `worklist: &MissWorklist`
- [x] 1.2 Column projection on all loaders (especially load_translations: skip 132 MB sift/polyphen)
- [x] 1.3 Interval predicates for regulatory/motif/mirna/structural

## 2. Phase 2: Rust-Side Transcript ID Filtering — DONE

- [x] 2.1 Filter exons by transcript_ids HashSet in Rust
- [x] 2.2 Filter translations by transcript_ids HashSet in Rust

## 3. Phase 3: Sift Prediction Optimization — DONE

- [x] 3.1 Compact sift predictions: HashMap<(i32,String),(String,f32)> → sorted Vec<CompactPrediction> with u8 encoding
- [x] 3.2 Zero-copy Arrow parsing: read &str directly from Arrow buffers, encode to u8 in-place
- [x] 3.3 Direct parquet-rs sift reader with cached ArrowReaderMetadata (bypass DataFusion SQL)
- [x] 3.4 Split profiling: separate sift_lazy_load_only from annotate_batches_only

## 4. Phase 4: VCF Output & Cleanup — DONE

- [x] 4.1 Add --output=<path> VCF sink to profile_annotation
- [x] 4.2 Proper VCF output with INFO, FORMAT, sample columns (using bio.vcf.field.field_type metadata)
- [x] 4.3 VCF round-trip tests with noodles-vcf
- [x] 4.4 Remove dead code (resolve_parquet_path, ProteinPrediction, read_protein_predictions, chrom_filter_clause, decode_amino_acid)
- [x] 4.5 CSQ string buffer reuse (write! instead of Vec<String> + join)
- [x] 4.6 Address PR review feedback (two rounds from Claude and Codex)

## 5. Phase 5: Benchmarking & Validation — DONE

- [x] 5.1 Profile with --everything: 107.7s → 79.2s (-27%)
- [x] 5.2 Golden benchmark: 80/80 CSQ fields at 100% accuracy, 0 mismatches
- [x] 5.3 All tests pass (587 tests), clippy clean, fmt clean

## Measured Results (chr1 --everything, 319K variants)

| Optimization | Sift load | Context load | Total |
|---|---|---|---|
| Baseline (master) | 48.0s | 13.6s | 107.7s |
| + Column projection | 48.0s | 0.5s | ~94s |
| + Compact predictions | 35.8s | 0.5s | ~92s |
| + Zero-copy Arrow parsing | 26.1s | 0.5s | ~82s |
| + Direct parquet-rs reader | 23.7s | 0.5s | 79.2s |

---

## 6. Phase 6: Sift Loading Gate for Fjall Backend — TODO

**Critical for fjall integration.** Without this, fjall gives 2.8x speedup (79s → 28s). With this, fjall gives 13-16x speedup (79s → 5-6s).

### 6.1 Gate sift window loading on cache-miss batches
- [ ] 6.1.1 In the annotation batch loop, before sift window loading, check if the batch has any rows where `cache_most_severe_consequence` is null
- [ ] 6.1.2 Skip sift window loading entirely for batches where ALL rows are cache hits
- [ ] 6.1.3 Preserve the `loaded_windows` HashSet and eviction logic for miss-containing batches
- [ ] 6.1.4 Verify correctness: parquet backend (0% hits) behavior is unchanged
- [ ] 6.1.5 Verify correctness: fjall backend (95%+ hits) produces identical output

### 6.2 Context loading optimization for high hit rates
- [ ] 6.2.1 When MissWorklist is empty (100% cache hits), skip ALL context loading (transcripts, exons, translations, regulatory, etc.)
- [ ] 6.2.2 When MissWorklist has few intervals, context loading already benefits from interval predicates (Phase 1)

### 6.3 Benchmark with fjall backend
- [ ] 6.3.1 Convert chr1-vep variation parquet → fjall using `load_cache` example
- [ ] 6.3.2 Profile annotation with fjall backend, without sift gate (expect ~28s)
- [ ] 6.3.3 Profile annotation with fjall backend + sift gate (expect ~5-6s)
- [ ] 6.3.4 Golden benchmark correctness with fjall backend
- [ ] 6.3.5 Compare: parquet 79s → fjall 28s → fjall+gate 5-6s

### 6.4 Estimated impact

| Scenario | variation_lookup | sift_load | annotate | Total |
|---|---|---|---|---|
| Parquet (0% hits) | 39s | 24s | 8s | **79s** |
| Fjall (95% hits), current sift | 2.5s | 24s | 0.5s | **~28s** |
| **Fjall (95% hits) + sift gate** | **2.5s** | **~2s** | **0.5s** | **~5-6s** |
| **Fjall (98% hits) + sift gate** | **2.5s** | **~1.5s** | **0.2s** | **~4-5s** |
