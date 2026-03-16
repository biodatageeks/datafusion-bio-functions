# Optimize Annotation Context Loading -- Design

## Context

The `scan_with_transcript_engine()` method in `annotate_provider.rs` is the core annotation pipeline. After the variation lookup resolves cache hits, it loads genomic context tables to run the transcript consequence engine on cache misses. This context loading is the second-largest bottleneck (~15s on chr1) after the variation lookup itself (~25s).

### Stakeholders

- VEP annotation pipeline users (faster end-to-end annotation)
- `bio-function-vep` crate maintainers
- polars-bio downstream consumers

## Goals / Non-Goals

**Goals:**
- Reduce `context_tables_total` from ~15s to <2s for typical workloads (localized miss sets)
- Reduce `context_tables_total` from ~15s to <5s even for broadly distributed miss workloads
- Load only the genomic context actually needed for the specific cache-miss variants
- Eliminate unnecessary column materialization via mode-aware projections
- Fix broken profiling counters so optimization impact can be measured accurately
- Remove unnecessary MemTable round-trip at the end of the pipeline

**Non-Goals:**
- Changing the variation lookup algorithm (covered by `add-fjall-variation-lookup`)
- Changing the transcript consequence engine itself
- Changing the `load_sift_window()` path (already optimized with interval predicates + column projection)
- Supporting non-parquet context table backends
- Changing the annotation output schema or semantics

---

## Alternatives Considered

| Alternative | Description | Verdict |
|---|---|---|
| Pre-materialized context per variant | Build a denormalized table joining all context for each known variant position | Enormous storage (~TB scale), impractical for cache |
| R-tree / interval index on context tables | Build spatial index over transcripts/regulatory features | Over-engineered; parquet RG min/max statistics with proper sort order achieve similar pruning |
| Lazy context loading per batch | Load context on-demand as each batch is annotated | More complex streaming architecture; eager loading with miss worklist is simpler and sufficient |
| **Miss worklist + interval predicates + transcript_id joins** | Scope loading to miss regions, then join by transcript ID | **Chosen** -- straightforward, leverages existing parquet/DataFusion predicate pushdown |

---

## Detailed Design

### 1. MissWorklist

New helper struct (either inline in `annotate_provider.rs` or in a new `miss_worklist.rs` module):

```rust
/// A genomic interval representing one or more coalesced cache-miss positions.
struct MissInterval {
    chrom: String,
    start: i64,
    end: i64,
}

/// Collected and coalesced cache-miss positions from lookup output.
struct MissWorklist {
    intervals: Vec<MissInterval>,
    /// Intervals grouped by normalized chromosome name for efficient lookup.
    by_chrom: HashMap<String, Vec<(i64, i64)>>,
    /// Total number of individual miss positions before coalescing.
    miss_count: usize,
    /// Total number of cache hits.
    hit_count: usize,
}
```

**Construction** (`collect_miss_worklist`):

1. Iterate **all** `base_batches` (no `.any()` short-circuit).
2. For each row where `cache_most_severe_consequence` is null, record `(chrom, start, end)`.
3. Normalize chromosome names (strip/add "chr" prefix for matching).
4. Sort positions per chromosome and merge intervals within a configurable padding distance (default: 1 Mb). The padding accounts for the fact that a variant at position P may need transcripts/regulatory features that start up to ~1 Mb away.
5. Record accurate hit/miss counts from the full scan.

**Coalescing rationale:**

Without coalescing, a VCF with 1000 scattered misses on chr1 would generate 1000 separate `OR` clauses. Coalescing with 1 Mb padding merges nearby misses into larger windows. For a typical WGS sample where misses cluster in regions with novel variants, this reduces to ~10-50 windows per chromosome while still excluding large swaths of the chromosome.

The padding value is conservative: the longest human transcript (RBFOX1) spans ~1.7 Mb, but >99% of transcripts are <500 Kb. A 1 Mb pad ensures we capture overlapping transcripts for virtually all miss positions.

### 2. Interval-Aware SQL Predicate Builder

Replace `chrom_filter_clause()` with:

```rust
fn interval_filter_sql(worklist: &MissWorklist) -> String
```

For each chromosome's coalesced intervals, generate:

```sql
(chrom = '1' AND start <= {interval_end} AND "end" >= {interval_start})
```

Multiple intervals become OR-ed clauses. Both "chr" and bare chromosome names are included (same as current `chrom_filter_clause` expansion logic).

**Predicate pushdown interaction with new parquet layouts:**

| Table | Sort key | RG count | Pushdown effect |
|---|---|---|---|
| transcript | `(chrom, start)` | 4-8 | Prunes RGs where `min(start) > interval_end` |
| regulatory | `(chrom, start)` | 4 | Same |
| motif | `(chrom, start)` | varies | Same |
| structural | `(chrom, start)` | varies | Same |

Note: The `"end" >= {interval_start}` predicate cannot benefit from start-sorted RG pruning (see Limitations section).

### 3. Rust-Side Transcript ID Filtering for Exons and Translations

After `load_transcripts()` returns the interval-filtered transcript set:

```rust
let transcript_ids: HashSet<String> = transcripts.iter()
    .map(|t| t.stable_id.clone())
    .collect();
```

Load exons and translations using the same interval/chrom filter + column projection (section 4), then filter in Rust:

```rust
// Load with interval predicates + column projection (already fast: ~125ms exons, ~300ms translations)
let mut exons = self.load_exons_projected(table, &worklist).await?;
exons.retain(|e| transcript_ids.contains(&e.transcript_id));

let mut translations = self.load_translations_projected(table, &worklist).await?;
translations.retain(|t| transcript_ids.contains(&t.transcript_id));
```

**Why not SQL `IN (...)`:**

- The query runs **once per chromosome** — exon/translation loads are already fast with column projection
- For distributed misses (typical WGS), the transcript set can be **40K+ IDs** on chr1. A SQL `IN (...)` clause that large hurts query planning more than it helps. Benchmarking showed DataFusion overhead is ~43ms per query for simple predicates; a 40K-element `IN` clause would be significantly worse.
- A Rust `HashSet::contains()` filter over loaded results costs microseconds for any set size
- Eliminates all complexity around semi-join fallbacks, temporary tables, and threshold tuning

### 4. Column Projection Constants

Define projection column lists as constants (or configuration) rather than inline strings:

```rust
/// Default transcript columns needed by the annotation engine.
const TRANSCRIPT_DEFAULT_COLS: &[&str] = &[
    "stable_id", "version", "chrom", "start", "end", "strand",
    "biotype", "source", "is_canonical",
    "gene_stable_id", "gene_symbol", "gene_symbol_source", "gene_hgnc_id", "refseq_id",
    "cds_start", "cds_end", "cdna_coding_start", "cdna_coding_end",
    "translation_stable_id", "exon_count",
    "cds_start_nf", "cds_end_nf", "flags_str", "cdna_mapper_segments",
    "mane_select", "mane_plus_clinical", "ccds", "tsl", "appris",
    "swissprot", "trembl", "uniparc", "uniprot_isoform", "gene_phenotype",
];

/// Additional transcript columns for sequence-dependent annotation.
const TRANSCRIPT_SEQUENCE_COLS: &[&str] = &[
    "spliced_seq", "cdna_seq", "translateable_seq",
    "codon_table", "peptide_seq",
];

/// Translation core columns (never includes sift/polyphen).
const TRANSLATION_CORE_COLS: &[&str] = &[
    "transcript_id", "stable_id", "version",
    "cds_len", "protein_len", "translation_seq", "cds_sequence",
];

/// Additional translation column for DOMAINS output.
const TRANSLATION_DOMAINS_COL: &str = "protein_features";

/// Exon columns.
const EXON_COLS: &[&str] = &[
    "transcript_id", "exon_number", "start", "end", "strand",
    "stable_id", "version", "phase", "end_phase",
    "is_current", "is_constitutive", "gene_stable_id",
];
```

Each loader builds the SELECT list by combining the required base columns with mode-dependent extras. If a column doesn't exist in the table schema (checked via DataFusion's schema introspection, which the loaders already do for `chrom`), it's silently skipped.

### 5. Cached Metadata Sift Window Reader

#### Problem

`load_sift_window()` is called 51 times per chromosome (5 Mb windows × 249 Mb chr1). Each call runs `session.sql()`.

Benchmarking (direct parquet-rs vs DataFusion `session.sql()` on the same windows) shows:
- **Direct parquet-rs** (cached metadata, RG selection, projection): **~486ms avg** per window
- **DataFusion `session.sql()`**: **~529ms avg** per window
- **DataFusion overhead: only ~43ms per query**

The per-window cost is dominated by **parquet I/O and decompression of nested `list<struct>` arrays**, not by query planning. Each window reads 1-5 RGs containing ~256-1,792 rows, where each translation has ~7,600 sift prediction items in a `list<struct<i32, string, string, f32>>` column.

Bulk-loading all sift/polyphen data is not viable: the nested arrays decompress to **~70 GB** for chr1 (22.8K translations × ~7,600 items each). The windowed approach with eviction must be preserved — removing it previously caused data mismatches due to missed transcripts.

#### Solution (two-part)

**Part 1 (code change, ~2.2s savings):** Cache `ArrowReaderMetadata` and use direct parquet-rs reads to eliminate the ~43ms per-query DataFusion overhead:

```rust
use parquet::arrow::arrow_reader::{ArrowReaderMetadata, ParquetRecordBatchReaderBuilder};
use parquet::arrow::ProjectionMask;

// ONCE before annotation loop: read footer, pre-compute RG ranges
let file = std::fs::File::open(&translation_path)?;
let arrow_metadata = ArrowReaderMetadata::load(&file, Default::default())?;
let parquet_schema = arrow_metadata.metadata().file_metadata().schema_descr_ptr();
let projection = ProjectionMask::roots(&parquet_schema, [end_idx, tid_idx, sift_idx, poly_idx]);

let rg_ranges: Vec<(i64, i64)> = /* extract start min / end max from each RG's stats */;

// PER WINDOW (replaces session.sql()):
let matching_rgs: Vec<usize> = rg_ranges.iter().enumerate()
    .filter(|(_, (s, e))| *s <= win_end && *e >= win_start)
    .map(|(i, _)| i)
    .collect();

let file = std::fs::File::open(&translation_path)?;  // just an fd
let reader = ParquetRecordBatchReaderBuilder::new_with_metadata(file, arrow_metadata.clone())
    .with_projection(projection.clone())
    .with_row_groups(matching_rgs)
    .build()?;

for batch in reader {
    // Same sift/polyphen parsing + cache insertion as current code
}
```

**Part 2 (layout change, ~10-15s savings):** The dominant cost is I/O + decompression of nested sift/polyphen arrays. Splitting the translation parquet ([datafusion-bio-formats#131](https://github.com/biodatageeks/datafusion-bio-formats/issues/131)) into a dedicated `translation_sift` file with 4-8 larger RGs:
- Eliminates 89 MB of 1 MB inter-RG alignment padding (detected in the current file: exactly 1,048,594 bytes between every consecutive RG)
- Reduces file size from ~229 MB to ~135 MB (data + padding → data only)
- Larger RGs mean fewer file seeks and more sequential reads for the nested arrays
- Expected per-window improvement: ~486ms → ~200-300ms

**What stays the same:** Window boundaries, `SIFT_WINDOW_SIZE`, `loaded_windows` tracking, `evict_before()` eviction, cache population logic, memory profile (~20 MB per window).

**Resolving the translation file path:**

The current code receives the translation table name (e.g., `"115_GRCh38_translation_1_vep"`) as a registered DataFusion table. To get the file path for direct parquet-rs access, either:
- Pass the file path through `AnnotateExec` alongside the table name (preferred — already available from the options JSON)
- Query DataFusion's catalog to resolve the registered parquet table's file path

**Remaining irreducible cost:**

After both optimizations, ~25s of sift loading remains for chr1 `--everything`. This is the cost of decompressing and materializing ~155M nested struct items across 51 windows. Further improvement would require a fundamentally different sift/polyphen data representation (e.g., position-keyed KV store, pre-computed lookup tables, or delta-encoded prediction arrays).

### 6. MemoryExec Replacement

Replace:
```rust
let mem = MemTable::try_new(projected_schema, vec![projected_batches])?;
mem.scan(state, None, &[], None).await
```

With:
```rust
let exec = MemoryExec::try_new(&[projected_batches], projected_schema, None)?;
Ok(Arc::new(exec))
```

`MemoryExec` is the execution plan that `MemTable` internally creates during `.scan()`. Using it directly skips the table abstraction layer and avoids the implicit schema validation pass.

### 6. Profiling Fix

Replace:
```rust
let has_cache_misses = base_batches.iter().any(|batch| { ... });
```

With:
```rust
let worklist = collect_miss_worklist(&base_batches)?;
let has_cache_misses = !worklist.intervals.is_empty();
// worklist.hit_count and worklist.miss_count are accurate
```

The profiling output now reports correct totals since `collect_miss_worklist` iterates all batches unconditionally.

---

## Execution Flow (After Changes)

```
scan_with_transcript_engine()
  │
  ├── 1. Collect base_batches from lookup output
  ├── 2. Build colocated map
  ├── 3. collect_miss_worklist(base_batches)        ← NEW
  │       → MissWorklist { intervals, by_chrom, hit_count, miss_count }
  │
  ├── if worklist.is_empty():
  │     → skip all context loading (same as today)
  │
  ├── 4a. load_transcripts(table, &worklist)         ← interval predicates
  │       → Vec<TranscriptFeature>
  │       → derive transcript_ids: HashSet<String>
  │
  ├── 4b. load_exons(table, &worklist) + filter by transcript_ids  ← projected + Rust filter
  ├── 4c. load_translations(table, &worklist) + filter by transcript_ids  ← projected + Rust filter
  ├── 4d. load_regulatory(table, &worklist)          ← interval predicates
  ├── 4e. load_motif(table, &worklist)               ← interval predicates
  ├── 4f. load_mirna(table, &worklist)               ← interval predicates
  ├── 4g. load_structural(table, &worklist)          ← interval predicates
  │
  ├── 5-8. PreparedContext, annotate batches, sift lazy load (unchanged)
  │
  └── 9. MemoryExec (not MemTable)                   ← simplified
```

---

## Limitations

1. **Interval overlap pruning is one-sided.** For start-sorted parquet, the predicate `start <= max_end` prunes "future" row groups, but `"end" >= min_start` cannot prune any RGs because `end` values within a start-sorted RG are uncorrelated with the RG boundaries. This is inherent to 1D sorting for interval data.

2. **Coalescing padding may over-fetch.** The 1 Mb default padding means isolated misses within 1 Mb of each other merge into one window. For pathological cases (evenly spaced misses every 500 Kb across chr1), coalescing would merge all misses into one giant window covering the entire chromosome. In practice, VCF misses cluster, making this unlikely.

3. **Schema evolution.** Column projection lists must be kept in sync with the parquet schemas produced by `datafusion-bio-formats`. If new columns are added to the parquet files, the projection constants need updating. Mitigated by silently skipping missing columns.

---

## Fjall Backend Interaction: Sift Loading Gate

### Problem

With the fjall KV backend, `cache_most_severe_consequence` is populated for known variants (~95-98% of WGS input). Cache-hit rows skip the transcript engine entirely, using pre-computed CSQ. However, sift window loading in the annotation batch loop is triggered by **batch position**, not by cache-miss status:

```rust
// Current code: loads sift for EVERY batch regardless of hits/misses
for batch in &base_batches {
    if sift_enabled {
        // Load sift windows for this batch's position range
        // → runs even when all rows in the batch are cache hits
    }
    annotated_batches.push(annotate_batch(..., &sift_cache)?);
}
```

With fjall + 95% hit rate, the annotation engine runs for only ~16K variants (5% of 321K), but sift windows are loaded for all 51 position windows because all 10,840 batches pass through the loop.

### Solution

Gate sift window loading on whether the batch contains any cache misses:

```rust
for batch in &base_batches {
    if sift_enabled {
        let has_miss = cached_most_idx.map_or(true, |idx| {
            (0..batch.num_rows()).any(|row| batch.column(idx).is_null(row))
        });
        if has_miss {
            // Load sift windows (existing code)
        }
    }
    annotated_batches.push(annotate_batch(..., &sift_cache)?);
}
```

With 95% hit rate, ~95% of batches contain only hits → only ~540 batches (5%) trigger sift loading → ~5 windows loaded instead of 51 → sift cost drops from 24s to ~2s.

### Estimated impact

| Scenario | variation_lookup | sift_load | annotate | Total |
|---|---|---|---|---|
| Parquet (current, 0% hits) | 39s | 24s | 8s | **79s** |
| Fjall, current sift code | 2.5s | 24s | 0.5s | **~28s** |
| **Fjall + miss-gated sift** | **2.5s** | **~2s** | **0.5s** | **~6s** |

---

## Expected Impact

### Measured baseline: chr1 `--everything` annotation (319K variants, HG002, 107.7s)

| Stage | Time | % |
|---|---|---|
| variation_lookup | 33.9s | 31% |
| **sift window loading** (~51 windows on translation parquet) | **~38-43s** | **~37%** |
| annotation engine (HGVS + consequence) | ~15-20s | ~17% |
| **load_translations** (`SELECT *`, reads 132 MB sift/polyphen it discards) | **13.1s** | **12%** |
| load_transcripts/exons/regulatory | 0.5s | <1% |
| other (colocated map, hydrate cdna, prepared context) | 2.5s | 2% |

### Translation parquet: corrected understanding

The 228.7 MB translation file breaks down as:
- **Sift/polyphen nested list data**: 131.9 MB (57.7%) -- `list<struct<position, amino_acid, prediction, score>>`
- **Page/header overhead from 90 tiny RGs**: 90.8 MB (39.7%)
- Core data + sequences: 5.8 MB (2.5%)
- Metadata/source columns: 0.13 MB (negligible)

The earlier analysis incorrectly attributed sift/polyphen physical leaf columns to metadata columns due to logical-to-physical index mismatch with nested types. Metadata columns are NOT a performance concern.

### Per-window sift benchmark results

Benchmarked direct parquet-rs vs DataFusion `session.sql()` on 10 windows (5 Mb each, spread across chr1):

| Approach | Avg | Min | Max |
|---|---|---|---|
| Direct parquet-rs (cached metadata) | 486ms | 162ms | 958ms |
| DataFusion session.sql() | 529ms | 164ms | 1019ms |
| **DataFusion overhead** | **~43ms** | | |

Per-window cost scales with row count (162ms for 256 rows, 958ms for 1,792 rows) and is dominated by decompression of nested `list<struct>` sift/polyphen arrays, NOT by query planning.

### Code-only improvements (no layout changes)

| Change | Savings | Mechanism |
|---|---|---|
| **Column projection on `load_translations`** | **~12.8s** | Skip 132 MB of sift/polyphen columns (read ~5 MB instead of 228 MB) |
| Cached metadata sift reader | ~2.2s | Eliminate 51 × 43ms DataFusion overhead |
| Miss worklist + interval predicates | ~0.4s | Scope transcript/regulatory loading to miss regions |
| MemoryExec, profiling fix | ~0.03s | Minor |
| **Code-only total** | **~15.4s** | **107.7s → ~92s (-14%)** |

**Column projection is the highest-value immediate fix** — it can ship today with zero layout dependency.

### Layout improvements (requires datafusion-bio-formats#131)

| Change | Savings | Mechanism |
|---|---|---|
| Split `translation_sift` + fix RG size | **~10-15s** | Fewer, larger RGs → more efficient I/O for nested arrays; eliminate 89 MB inter-RG alignment padding |
| Multi-RG transcript/regulatory tables | ~0.2s | Enable interval predicate RG pruning |
| Re-sort exon by `(transcript_id, start)` | ~0.05s | Enable transcript_id RG pruning |
| **Layout total** | **~10-15s** | |

### Combined improvements

| Scenario | Runtime | Speedup |
|---|---|---|
| **Current baseline** | **107.7s** | — |
| Code changes only | ~92s | -15% |
| **Code + layout combined** | **~78-83s** | **-23-28%** |

### Irreducible sift cost

After all optimizations, ~25s of sift window loading remains for chr1 `--everything`. This is the cost of decompressing ~155M nested struct items (`list<struct<i32, string, string, f32>>`) across 51 windows. Further improvement requires changing the sift/polyphen data representation (e.g., position-keyed KV store, pre-computed lookup tables).

### Context loading impact for different miss workloads

After code + layout changes:

| Workload | context_tables before | context_tables after |
|---|---|---|
| Localized (100 variants in 5 Mb region) | ~13.6s | ~0.2s |
| Moderate (5000 variants, clustered) | ~13.6s | ~0.5s |
| Distributed (319K variants, full chr1) | ~13.6s | ~0.3s (projection alone handles most of it) |

Note: context_tables savings are dominated by column projection in all cases because the 132 MB sift/polyphen data is skipped entirely once `SELECT *` is replaced with explicit column lists.
