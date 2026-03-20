# Optimize Annotation Context Loading for Parquet-Backed VEP

## Why

The parquet-backed `annotate_vep()` annotation path has two major bottlenecks in context loading. Profiling of chr1 `--everything` annotation (319K variants from HG002, 107.7s total) shows:

| Stage | Time | % |
|-------|------|---|
| variation_lookup | 33.9s | 31% |
| **sift window loading** (~51 windows on translation parquet) | **~38-43s** | **~37%** |
| annotation engine (HGVS + consequence) | ~15-20s | ~17% |
| **load_translations** (`SELECT *`) | **13.1s** | **12%** |
| load_transcripts/exons/regulatory | 0.5s | <1% |
| other | 2.5s | 2% |

The two translation-related stages account for **~51-56s (48-52%)** of total runtime:

1. **`load_translations` uses `SELECT *`** which reads 132 MB of sift/polyphen prediction data that the main loader immediately discards. These columns are loaded separately by the windowed `load_sift_window()` path. The actual core translation data is only ~5 MB.

2. **Sift window loading queries a bloated parquet file 51 times.** The translation parquet has 90 row groups of 256 rows (instead of 4-8 RGs of 3-6K rows), creating 90.8 MB of page/header overhead and a footer with 4,320 column-chunk stat entries (90 RGs × 48 physical cols). Each `session.sql()` window query must process this entire footer, costing ~800-900ms per query.

3. **Any cache miss triggers full-chromosome loading** -- a single miss on chr1 loads ALL 47.8K transcripts, 355.8K exons, 22.8K translations, and 35.9K regulatory features for the entire chromosome.

4. **Exons and translations are filtered by chromosome, not by transcript ID.**

Additionally, the profiling counters are broken (`.any()` short-circuits hit/miss counting), and the final annotated batches take an unnecessary round-trip through `MemTable`.

This proposal covers the code-side optimizations. It is complemented by parquet layout improvements in [datafusion-bio-formats#131](https://github.com/biodatageeks/datafusion-bio-formats/issues/131):
- Translation split into `translation_core` (sorted by `transcript_id`, 4-8 RGs) and `translation_sift` (sorted by `chrom, start`, 4-8 RGs, 6 columns)
- Transcript table: 4-8 RGs sorted by `(chrom, start)`
- Exon table: 8-16 RGs sorted by `(transcript_id, start)`
- Regulatory table: 4 RGs sorted by `(chrom, start)`
- `sorting_columns` declared in parquet footers

**Note:** metadata/source columns (species, assembly, cache_version, source_*) are negligible (0.13 MB in translation) and do not need to be dropped for performance.

## What Changes

### A. Miss Worklist: Replace chrom-only miss gate with concrete miss intervals

Replace the current `has_cache_misses` boolean + `vcf_chroms` set with a `MissWorklist` that collects the exact genomic intervals of cache-miss variants. This scopes all downstream context loading to the regions that actually need annotation.

**Current flow:**
```
base_batches → .any(has_miss?) → if true: load ALL context for vcf_chroms
```

**New flow:**
```
base_batches → collect_miss_worklist() → if non-empty: load context for miss intervals only
```

The worklist merges nearby/overlapping miss positions per chromosome into coalesced intervals (with configurable padding, default 1 Mb) to keep predicate lists compact.

### B. Interval-aware context loading for position-queried tables

Replace `chrom_filter_clause()` with interval-aware predicates for tables queried by genomic position:

**Affected loaders:**
- `load_transcripts` -- `WHERE chrom = '1' AND start <= {max_end} AND "end" >= {min_start}`
- `load_regulatory_features` -- same pattern
- `load_motif_features` -- same pattern
- `load_mirna_features` -- same pattern
- `load_structural_features` -- same pattern

For multiple intervals on the same chromosome, emit a disjunction:
```sql
WHERE (chrom = '1' AND start <= 1000500 AND "end" >= 999500)
   OR (chrom = '1' AND start <= 5000500 AND "end" >= 4999500)
```

With the revised parquet layouts (multiple RGs sorted by `chrom, start`), DataFusion prunes row groups via min/max statistics on `start`, skipping RGs that cannot overlap any miss interval.

### C. Filter exons and translations by transcript ID in Rust

After loading transcripts (which now returns only those overlapping miss intervals), derive the set of needed transcript IDs:

```rust
let transcript_ids: HashSet<String> = transcripts.iter().map(|t| t.stable_id.clone()).collect();
```

Then load exons and translations with the same interval/chrom filter + column projection (changes D), and filter in Rust:

```rust
let exons = self.load_exons_projected(table, &worklist).await?;
let exons: Vec<_> = exons.into_iter()
    .filter(|e| transcript_ids.contains(&e.transcript_id))
    .collect();
```

This is simpler and faster than SQL `IN (...)` because:
- The query runs once per chromosome — exon load is already ~125ms, translation ~0.3s (with projection)
- For distributed misses, the transcript set can be 40K+ IDs — a SQL `IN (...)` clause that large would hurt query planning more than it helps
- A Rust `HashSet` filter over the loaded results costs microseconds
- No need for semi-join fallback logic or temporary tables

### D. Mode-aware column projection

Replace `SELECT *` with explicit column lists in all loaders. The required columns depend on annotation mode.

**Transcript default projection:**
```
stable_id, version, chrom, start, end, strand, biotype, source, is_canonical,
gene_stable_id, gene_symbol, gene_symbol_source, gene_hgnc_id, refseq_id,
cds_start, cds_end, cdna_coding_start, cdna_coding_end,
translation_stable_id, exon_count,
cds_start_nf, cds_end_nf, flags_str, cdna_mapper_segments,
mane_select, mane_plus_clinical, ccds, tsl, appris,
swissprot, trembl, uniparc, uniprot_isoform, gene_phenotype
```

**Columns loaded only when needed:**
- `spliced_seq`, `cdna_seq`, `translateable_seq` -- only when transcript engine needs sequence data
- `codon_table`, `peptide_seq` -- only for protein consequence calculation
- `bam_edit_status`, `has_non_polya_rna_edit` -- only when BAM edits are relevant
- `mature_mirna_regions` -- only for miRNA annotation
- `ncrna_structure` -- only for ncRNA annotation

**Translation core default projection:**
```
transcript_id, stable_id, version, cds_len, protein_len,
translation_seq, cds_sequence
```

**Load `protein_features` only when DOMAINS output is requested.**

Do **not** load `sift_predictions` / `polyphen_predictions` in the main translation loader -- these remain in the existing windowed `load_sift_window()` path which queries `translation_sift` by position.

**Exon projection:**
```
transcript_id, exon_number, start, end, strand, stable_id, version,
phase, end_phase, is_current, is_constitutive, gene_stable_id
```

**Regulatory / motif / miRNA / structural:** project only columns actually read by the corresponding feature struct constructor.

### E. Reduce sift window loading cost via layout change

The current `load_sift_window()` is called 51 times per chromosome (5 Mb windows covering chr1's 249 Mb). Benchmarking shows:

- **Direct parquet-rs read** (cached metadata, RG selection, column projection): **~486ms avg** per window
- **DataFusion `session.sql()`** (same query): **~529ms avg** per window
- **DataFusion overhead: only ~43ms per query** — not the ~800ms previously estimated

The per-window cost is dominated by **parquet I/O and decompression of the nested sift/polyphen `list<struct>` arrays**, not by DataFusion planning. Reading 1-2 RGs of 256 rows each containing ~7,600 prediction items per translation is inherently expensive (~160ms for 256 rows, ~960ms for 1,792 rows).

Bulk-loading all sift/polyphen data is not viable: the nested arrays decompress to **~70 GB** for chr1 (22.8K translations × ~7,600 prediction items each). The windowed approach with eviction must be preserved for memory safety and correctness (removing it previously caused data mismatches).

**The primary optimization is therefore a layout change, not a code change.** The current translation parquet has 90 RGs of 256 rows with 1 MB alignment padding between each RG (89 MB total padding). Splitting into a dedicated `translation_sift` file with 4-8 larger RGs eliminates:
- 89 MB of inter-RG padding (the file shrinks from ~229 MB to ~135 MB)
- Per-RG page header overhead across fewer, larger groups
- Redundant column data (only 6 columns instead of 48)

With larger RGs, each window reads fewer, bigger batches — which is more I/O efficient for the nested list arrays. Expected per-window improvement: ~486ms → ~200-300ms, saving ~10-15s total across 51 windows.

A secondary code optimization is to **cache the `ArrowReaderMetadata`** and reuse it across all 51 window queries, saving the ~43ms per-query DataFusion overhead (51 × 43ms = ~2.2s). This is modest but free:

```rust
// ONCE before annotation loop:
let file = std::fs::File::open(&translation_path)?;
let arrow_metadata = ArrowReaderMetadata::load(&file, Default::default())?;

// PER WINDOW (replaces session.sql()):
let file = std::fs::File::open(&translation_path)?;
let reader = ParquetRecordBatchReaderBuilder::new_with_metadata(file, arrow_metadata.clone())
    .with_projection(projection.clone())
    .with_row_groups(matching_rgs)
    .build()?;
```

Window boundaries, `SIFT_WINDOW_SIZE`, `loaded_windows` tracking, `evict_before()` logic, cache population, and memory profile all remain unchanged.

### F. Replace MemTable round-trip with MemoryExec

Replace the current:
```rust
let mem = MemTable::try_new(projected_schema, vec![projected_batches])?;
mem.scan(state, None, &[], None).await
```

With a direct `MemoryExec`:
```rust
let exec = MemoryExec::try_new(&[projected_batches], projected_schema, None)?;
Ok(Arc::new(exec))
```

This avoids the table registration overhead and unnecessary extra scan.

### G. Fix miss-rate profiling counters

Replace the `.any()` pattern with a full scan that computes both the boolean and accurate counters:

```rust
let mut cache_hit_count = 0usize;
let mut cache_miss_count = 0usize;
for batch in &base_batches {
    // count all hits/misses across all batches
}
let has_cache_misses = cache_miss_count > 0;
```

## Expected Impact

Measured on chr1 `--everything` annotation (319K variants, HG002, 107.7s total).

Benchmarking confirmed that DataFusion `session.sql()` overhead is only ~43ms per sift window query — the ~500ms per-window cost is dominated by parquet I/O and decompression of nested `list<struct>` arrays, not by query planning.

| Change | Type | Savings | Mechanism |
|--------|------|---------|-----------|
| **D. Column projection on `load_translations`** | Code | **~12.8s** | Skip 132 MB sift/polyphen in main loader (`SELECT cols` not `SELECT *`) |
| **E. Cached metadata sift reader** | Code | **~2.2s** | Eliminate 51× DataFusion planning overhead (51 × 43ms) |
| B. Interval predicates + C. transcript_id joins | Code | ~0.4s | Scope context loading to miss regions |
| A. Miss worklist | Code | ~0s (enabler) | Foundation for B, C, and accurate profiling |
| F. MemoryExec + G. Profiling fix | Code | ~0.03s | Cleanup |
| **E. Sift layout change** (split `translation_sift`, fix RG size) | Layout | **~10-15s** | Fewer, larger RGs → more efficient I/O for nested arrays; eliminate 89 MB inter-RG padding |
| **Total code-only savings** | | **~15.4s** | |
| **Total code + layout savings** | | **~25-30s** | |

| Scenario | Runtime | Speedup |
|----------|---------|---------|
| **Current baseline** | **107.7s** | — |
| Code changes only (D + E + B/C + rest) | ~92s | **-15%** |
| **Code + layout changes** | **~78-83s** | **-23-28%** |

The single highest-impact code-only change is **column projection on `load_translations` (D)** — skip 132 MB of sift/polyphen data in the main loader for ~12.8s savings. This can ship immediately.

The highest-impact layout change is **splitting and re-sizing the translation parquet** ([datafusion-bio-formats#131](https://github.com/biodatageeks/datafusion-bio-formats/issues/131)) — fewer, larger RGs with no 1 MB inter-RG padding improves sift window I/O by ~10-15s.

**Note on the sift bottleneck:** The remaining ~25s of sift window loading (after all optimizations) is dominated by decompressing and materializing nested `list<struct<i32, string, string, f32>>` arrays — each window reads ~400 translations with ~7,600 prediction items each (~1.2M struct items). This is a fundamental I/O + decompression cost that cannot be eliminated by query planning or layout optimizations alone. Further improvement would require changing the sift/polyphen data encoding (e.g., pre-computed lookup tables, position-keyed KV store similar to the variation cache).

## Affected Code

- `datafusion/bio-function-vep/src/annotate_provider.rs` -- main changes (miss worklist, interval loading, projections, MemoryExec)
- `datafusion/bio-function-vep/src/annotate_table_function.rs` -- update test table schemas if projections change loader contracts
- Possibly a new `datafusion/bio-function-vep/src/miss_worklist.rs` helper module

**No breaking changes** -- `annotate_vep()` output is behaviorally identical; only internal loading strategy changes.

## Relationships

- **Complemented by** [datafusion-bio-formats#131](https://github.com/biodatageeks/datafusion-bio-formats/issues/131) for parquet layout improvements (translation split, multi-RG, sort keys). Code changes can ship independently; layout changes multiply their impact.
- **Complements** `add-fjall-variation-lookup` which optimizes the variation lookup step (~33.9s); this proposal optimizes the context loading + sift steps (~51-56s)
- The `load_sift_window()` query logic is **not changed** -- it already uses interval predicates and column projection. It benefits from the layout change (smaller footer, fewer RGs) without code changes.

## Critical Interaction: Sift Loading Gate for Fjall Backend

**Measured after implementation:** Switching to the fjall KV backend alone reduces variation_lookup from 39s to ~2.5s, but sift loading remains at ~24s because the current code loads sift windows for ALL batch positions regardless of cache hit/miss status. With fjall providing `cache_most_severe_consequence` and 95%+ hit rate, most variants are cache hits that skip the transcript engine entirely — but sift windows are still loaded for their position range.

**Estimated impact with fjall:**

| Scenario | Time | vs 79s baseline |
|----------|------|-----------------|
| Fjall only (current sift code) | ~28s | -64% |
| **Fjall + miss-gated sift loading** | **~5-6s** | **-92 to -94%** |

**Required change:** Gate sift window loading on batches that contain cache misses:

```rust
// In the annotation batch loop, before loading sift windows:
if sift_enabled {
    // Check if this batch has any cache misses that need the transcript engine.
    let has_miss = batch_has_cache_miss(batch, cached_most_idx);
    if !has_miss {
        // All rows are cache hits — skip sift loading for this batch.
        // Cache-hit path uses pre-computed CSQ, no SIFT/PolyPhen needed.
    } else {
        // Load sift windows for this batch's position range (existing code).
    }
}
```

With 95% cache hit rate, ~95% of batches contain only cache hits → only ~5% of batches trigger sift window loading → ~5 windows loaded instead of 51 → sift cost drops from 24s to ~2s.

This is a ~5-line change that transforms the fjall backend from a 2.8x speedup to a 13-16x speedup. It should be implemented alongside or immediately after the fjall backend integration.
