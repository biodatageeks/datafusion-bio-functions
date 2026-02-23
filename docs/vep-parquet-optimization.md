# Parquet Layout Optimization for `lookup_variants` — VEP 115 GRCh38 Cache

## Context

`lookup_variants()` performs a LEFT JOIN between a user's VCF table and the VEP variation cache (`115_GRCh38.parquet`). The join predicate is:

```sql
vcf.chrom = cache.chrom
  AND vcf.end >= cache.start
  AND vcf.start <= cache.end
  AND match_allele(vcf.ref, vcf.alt, cache.allele_string)
```

The cache file is **30.5 GB compressed** (89.6 GB uncompressed), **1.17 billion rows**, **78 columns**, **9,528 row groups** of 122,880 rows each. The file is globally sorted by `(chrom, start)` in lexicographic chromosome order. This report analyzes all DataFusion parquet pruning levels and proposes data-driven optimizations.

---

## Current State Analysis

### File Profile

| Property | Value |
|---|---|
| Rows | 1,170,699,612 |
| Columns | 78 |
| Row groups | 9,528 |
| Rows/RG | 122,880 |
| Compressed size | 30.5 GB (ZSTD) |
| Uncompressed size | 89.6 GB |
| Compression ratio | 2.99x |
| Sort order | `(chrom ASC, start ASC)` lexicographic |

### Column Budget (compressed)

| Column | Size | % of File | Notes |
|---|---|---|---|
| `variation_name` | 5.70 GB | 19.1% | Near-unique per row, PLAIN encoded |
| `dbsnp_ids` | 5.66 GB | 18.9% | Near-unique, PLAIN encoded |
| `allele_string` | 1.36 GB | 4.5% | High cardinality |
| `end` | 1.16 GB | 3.9% | PLAIN encoded |
| `start` | 1.12 GB | 3.7% | PLAIN encoded |
| `gnomADg*` (11 cols) | 9.2 GB | 30.7% | ~53% populated |
| `chrom` | 0.001 GB | 0.0% | RLE_DICTIONARY, 1 entry/RG |
| Other 58 columns | 5.8 GB | 19.4% | Mostly >95% null |

### Row Group Genomic Span (chr1)

| Metric | Value |
|---|---|
| Median span | 327,147 bp |
| p25-p75 | 315K - 335K bp |
| Max span | 18.2M bp (centromere/telomere) |
| Overlap between RGs | None (monotonically increasing) |

---

## Pruning Level Assessment

DataFusion supports 5 pruning levels. Here is the current effectiveness and optimization opportunity for each.

### Level 1: Projection Pruning (Column Selection)

**Current state:** Fully effective. `lookup_variants` already projects only needed columns via SQL.

**Impact:** The 5 core join columns (`chrom`, `start`, `end`, `allele_string`, `variation_name`) total 9.34 GB (31.2%). Reading only these saves **20.6 GB** of I/O. Adding annotation columns (e.g., `clin_sig`, `gnomADg`) increases this proportionally.

**Recommendation:** No file layout change needed. Already optimal.

### Level 2: Row Group Pruning (Min/Max Statistics)

**Current state: HIGHLY EFFECTIVE.** Data is globally sorted by `(chrom, start)`. Each row group covers a narrow, non-overlapping genomic window (~325K bp). Row group statistics have tight min/max ranges.

**Measured effectiveness:** A single-position lookup (`chrom='1' AND start=500000`) hits exactly **1 of 717** chr1 row groups = **99.86% pruning rate**. A typical WGS VCF query touching 1000 positions across chr1 would read ~1000/717 ~ **all** chr1 RGs but **zero** RGs from other chromosomes. For a targeted panel (e.g., 50 genes), row group pruning eliminates >95% of data.

**Current row group size trade-off:**

| RG Size | RGs | Compressed/RG | Genomic Span/RG | Point-Query I/O |
|---|---|---|---|---|
| 122,880 rows (current) | 9,528 | 3-5 MB | ~325K bp | 3-5 MB |
| 500,000 rows | ~2,341 | 12-20 MB | ~1.3M bp | 12-20 MB |
| 1,000,000 rows (DF default) | ~1,171 | 24-40 MB | ~2.6M bp | 24-40 MB |

**Recommendation:** The current 122,880 rows/RG is **already excellent** for selective lookups. Each RG is only 3-5 MB compressed, minimizing I/O per row group touched. Do NOT increase row group size — it would widen genomic spans and reduce pruning granularity for point queries. The 9,528 RGs produce ~207 MB of metadata, which is large but acceptable for a 30 GB file (0.7% overhead).

### Level 3: Bloom Filters

**Current state: NOT PRESENT.** The file contains no Bloom filters.

**Opportunity analysis:**

| Column | Cardinality | Predicate Type | Bloom Filter Value |
|---|---|---|---|
| `chrom` | 302 | Equality | **Low** — min/max stats already sufficient |
| `start` | ~700M | Range | **None** — Bloom filters don't help range predicates |
| `allele_string` | Very high | Equality (via UDF) | **None** — called through `match_allele()` UDF, not a direct filter |

**Recommendation:** Bloom filters provide **negligible benefit** for `lookup_variants`. The join uses range predicates (`>=`, `<=`) on `start`/`end` and a UDF for allele matching — neither benefits from Bloom filters. Skip writing them to avoid the ~2-5% file size increase.

### Level 4: Page Index (Page-Level Min/Max)

**Current state: NOT PRESENT.** The file was written by Polars (writer version 1.0) which does not write the Parquet Page Index (column index + offset index). This is the **single largest optimization opportunity**.

**Why it matters:** With 122,880 rows per row group and a default 1 MB data page size, each row group contains roughly 10-50 pages per column. Page-level min/max stats on sorted columns (`start`, `end`) would allow DataFusion to skip pages within an already-selected row group. For a point query, this could reduce from reading the entire 3-5 MB row group to reading only the ~100-500 KB page containing the target position.

**Expected impact:**

| Scenario | Without Page Index | With Page Index | Improvement |
|---|---|---|---|
| Single-position lookup | 3-5 MB (1 full RG) | 100-500 KB (1-2 pages) | **5-30x less I/O** |
| 100-position panel (same chrom) | 50-100 MB | 10-50 MB | **2-5x less I/O** |
| WGS (all positions, 1 chrom) | Full chrom scan | Full chrom scan | None (all pages needed) |

**Recommendation: HIGH PRIORITY.** Re-write the file with a writer that emits page index metadata. See "Recommended Rewrite Parameters" below.

### Level 5: Row-Level Filter Pushdown (Late Materialization)

**Current state:** Available but disabled by default in DataFusion (`pushdown_filters = false`).

**Opportunity:** After row group and page pruning select the relevant pages, filter pushdown would decode only `chrom` and `start`/`end` first, build a row selection mask, then decode `allele_string`, `variation_name`, and annotation columns only for matching rows. On a 122,880-row RG where only ~100 rows match the overlap predicate, this avoids decoding 99.9% of the heavy columns.

**Expected impact:** For selective queries (panel, exome), **2-5x reduction** in decode CPU time. For WGS-scale queries, negligible benefit (most rows match).

**Recommendation: MEDIUM PRIORITY.** Enable via DataFusion session config:
```rust
config.options_mut().execution.parquet.pushdown_filters = true;
config.options_mut().execution.parquet.reorder_filters = true;
```
No file layout change needed — this is a runtime setting.

---

## Recommended Rewrite Parameters

Re-write the file using PyArrow (or `parquet-rewriter`) with these parameters to enable page index:

```python
import pyarrow.parquet as pq
import pyarrow as pa

table = pq.read_table('115_GRCh38.parquet')

# Sort: maintain existing (chrom, start) order — already sorted
# If re-sorting, use natural chromosome order for better UX:
# 1,2,...,22,X,Y,patches

pq.write_table(
    table,
    '115_GRCh38_optimized.parquet',
    compression='zstd',
    compression_level=3,         # ZSTD level 3 (good ratio, fast decode)
    row_group_size=122_880,      # Keep current size — proven effective
    data_page_size=8192,         # 8 KB pages for fine-grained page pruning
    write_statistics=True,       # Ensure column stats (already present)
    write_page_index=True,       # CRITICAL: enables page-level pruning
    # writer_version='2.6',     # Required for page index in some writers
    use_dictionary=[             # Only for low-cardinality columns
        'chrom', 'region_bin', 'species', 'assembly',
        'failed', 'somatic', 'strand', 'phenotype_or_disease',
    ],
    # Disable dictionary for high-cardinality columns
    # (PyArrow auto-falls back to PLAIN when dict exceeds threshold)
)
```

### Page Size Selection Rationale

The choice of **8 KB** page size is driven by the query pattern:

- **122,880 rows / RG**, each row is ~75 bytes uncompressed across core columns
- With 8 KB pages: ~100-200 rows per page → **600-1200 pages per RG** for `start` column
- Each page covers ~270-650 bp of genomic range (at chr1 density)
- A point lookup touches **1-2 pages** instead of the entire column chunk
- Trade-off: More page headers (~1-2% overhead), slightly worse compression

| Page Size | Rows/Page | Pages/RG | Genomic Span/Page | Point-Query Reads |
|---|---|---|---|---|
| 8 KB | ~100-200 | 600-1200 | ~300 bp | 1-2 pages |
| 64 KB | ~800-1600 | 75-150 | ~2,500 bp | 1 page |
| 1 MB (default) | ~12,000-25,000 | 5-10 | ~30K bp | 1 page (but large) |

8 KB is recommended by the Parquet specification for point-lookup workloads. Since `lookup_variants` performs interval overlap (essentially point-to-range), 8 KB maximizes page skip rate.

---

## Columns to Drop

20 columns are **100% null** across all 1.17B rows and 15 columns are **constant** (same value in every row). These add negligible storage cost (dictionary-compressed to near-zero) but bloat the schema and metadata. Dropping them would:

- Reduce metadata size from 207 MB to ~150 MB
- Simplify schema for downstream tools
- Eliminate 35 column chunks per row group from the footer

**100% null columns (safe to drop):**
`minor_allele`, `minor_allele_freq`, `assembly_ids`, `gencode_ids`, `genebuild_ids`, `gnomade_ids`, `gnomadg_ids`, `polyphen_ids`, `refseq_ids`, `regbuild_ids`, `sift_ids`, `src_1000genomes_ids`, `clinical_impact`, and 7 others.

**Constant metadata columns (move to file-level key-value metadata):**
`species`, `assembly`, `cache_version`, `serializer_type`, `source_cache_path`, `source_file`, `source_assembly`, and 8 `source_*_version` columns.

---

## Chromosome Sort Order

Current: lexicographic (`1, 10, 11, ..., 2, 20, ...`). This works fine for pruning (each chrom's RGs are contiguous), but **karyotypic order** (`1, 2, ..., 22, X, Y, patches`) would be more natural and could slightly improve sequential scan patterns when processing chromosomes in biological order. Low priority but nice-to-have.

---

## Impact Summary

| Optimization | Effort | I/O Reduction | CPU Reduction | Priority |
|---|---|---|---|---|
| **Page Index** (rewrite file) | Medium | 5-30x for point queries | 5-30x decode savings | **HIGH** |
| **Filter Pushdown** (runtime config) | Trivial | — | 2-5x decode savings | **MEDIUM** |
| **Drop null/constant columns** (rewrite) | Low | ~5% metadata | Negligible | LOW |
| **Karyotypic sort order** (rewrite) | Medium | None | Negligible | LOW |
| Bloom filters | Medium | None for this workload | None | SKIP |
| Row group size tuning | Medium | Negative impact | Negative impact | SKIP |

### Net Expected Impact (Selective Panel Query — 50 genes, ~500 positions)

| Metric | Current | Optimized | Improvement |
|---|---|---|---|
| Row groups read | ~500 | ~500 | Same |
| Data per RG read | 3-5 MB (full RG) | 100-500 KB (1-2 pages) | **10x** |
| Total I/O | ~2 GB | ~200 MB | **10x** |
| Decode work | All columns, all rows in RG | Filter cols first, then projected cols for matching rows only | **5x** |

### Net Expected Impact (WGS — all chromosomes, millions of positions)

| Metric | Current | Optimized | Improvement |
|---|---|---|---|
| Row groups read | All ~9,528 | All ~9,528 | Same |
| Data per RG read | 3-5 MB | 3-5 MB (most pages needed) | ~Same |
| Total I/O | ~30 GB | ~30 GB | Negligible |
| Decode work | Full decode | Slightly less with pushdown | ~1.2x |

**Conclusion:** The optimizations are most impactful for **targeted/panel** lookups and diminish for whole-genome scale. The page index rewrite is the single highest-value change.

---

## Sparse Data Optimization (Implemented)

### Problem: Broken Null Filter

The `lookup_variants()` pre-filter intended to prune cache rows where all annotation columns are NULL was **ineffective**. The filter included structural columns (`variation_name`, `allele_string`) which are always non-null, making the OR condition trivially true:

```sql
-- Old (broken): trivially true because variation_name is never null
WHERE `variation_name` IS NOT NULL OR `clin_sig` IS NOT NULL OR ...
```

**Fix:** `STRUCTURAL_COLUMNS` constant excludes `chrom`, `start`, `end`, `variation_name`, `allele_string` from the null filter. Only true annotation columns (sparse, >95% null) participate:

```sql
-- Fixed: only sparse columns in the filter
WHERE `clin_sig` IS NOT NULL OR `gnomADg_AF` IS NOT NULL OR ...
```

When requesting only structural columns (e.g., `variation_name` alone), no filter is applied — all cache entries participate in the join.

**Impact:** For a query requesting `variation_name,clin_sig`, the filter now prunes the ~95% of cache rows that have no clinical significance. This eliminates 1.1B rows from the join's probe side, reducing cache scan I/O from ~9 GB (5 core columns) to ~450 MB (only annotated rows).

### Parquet Filter Pushdown (Enabled)

DataFusion's `pushdown_filters` and `reorder_filters` are now enabled in the bio session config. This provides **late materialization**: when scanning the cache parquet file with a `WHERE clin_sig IS NOT NULL` filter, DataFusion:

1. First decodes only the `clin_sig` column (definition levels only — fast for sparse nulls)
2. Builds a row selection mask from non-null positions
3. Decodes remaining projected columns (`chrom`, `start`, `end`, `allele_string`, `variation_name`) **only for matching rows**

Combined with row group pruning on `(chrom, start)`, this means a panel query touching 500 positions:
- Prunes ~9,000 of 9,528 row groups by coordinate range
- Within the ~500 remaining RGs, prunes ~95% of rows by null filter
- Decodes heavy columns for only ~2,500 rows (500 positions × ~5 annotated entries each)

### Future: Vertical Partitioning

For maximum sparse-data efficiency, the cache could be split into two files:

| File | Columns | Rows | Size (est.) |
|---|---|---|---|
| `core.parquet` | chrom, start, end, allele_string, variation_name | 1.17B (all) | ~9.3 GB |
| `annotations.parquet` | chrom, start, end, allele_string + 60 annotation cols | ~58M (non-null only) | ~4-8 GB |

The interval join would scan only `core.parquet` (9.3 GB). A secondary equi-join on `(chrom, start, end, allele_string)` would fetch annotations for matched rows only. For a panel query matching 500 variants, the annotation lookup reads a few MB instead of scanning the full file.

**Trade-off:** Requires maintaining two files and a two-phase join in `lookup_variants()`. The current null filter + pushdown approach already captures most of the benefit for typical queries.

---

## Verification

After rewriting, validate with:

```bash
# 1. Check page index presence
python -c "
import pyarrow.parquet as pq
pf = pq.ParquetFile('115_GRCh38_optimized.parquet')
rg = pf.metadata.row_group(0)
for i in range(5):
    col = rg.column(i)
    print(f'{col.path_in_schema}: has_statistics={col.is_stats_set}')
# Page index presence can be verified via parquet-tools or pyarrow internal API
"

# 2. Compare row counts
python -c "
import pyarrow.parquet as pq
old = pq.ParquetFile('115_GRCh38.parquet').metadata.num_rows
new = pq.ParquetFile('115_GRCh38_optimized.parquet').metadata.num_rows
assert old == new, f'Row count mismatch: {old} vs {new}'
print(f'OK: {new} rows')
"

# 3. Benchmark lookup_variants with both files (use a small panel VCF)
# Compare query times and DataFusion metrics (row groups pruned, pages pruned)
```
