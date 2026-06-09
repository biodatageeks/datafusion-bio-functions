# Lance Variation Cache Layout Design

## Goal

Design Lance-backed variation-cache artifacts for the VEP `everything` path, using both:

- Lance 2.1 as the stable baseline.
- Lance 2.2 as an enhanced benchmark target for features that are not available or not equally effective in 2.1.

The design must optimize two access patterns:

- Warm cache: full scan over the columns needed by the `everything` path.
- Cold cache: batched point lookups, with a target batch size of 2,000 positions.

The logical annotation output must remain compatible with the current `indexed_parquet` path. Physical schema changes are allowed, including struct packing, as long as runtime reconstruction exposes the same logical fields.

## Inputs Reviewed

The current implementation discovers requested variation cache columns in `annotate_provider.rs`:

- It starts with `consequence_types` and `most_severe_consequence`.
- It appends `cache_lookup_column_names()`.
- It filters that list to columns present in the variation artifact.
- Warm and cold readers then add runtime key/match columns such as `position_key`, `variant_keys`, `allele_string`, `end`, and `failed`.

For the profiled chr1 cache, `consequence_types` and `most_severe_consequence` are not present in `chr1_warm.parquet` or `chr1_cold.parquet`, so the current runtime does not read them from variation storage.

Profiled files:

```text
/Users/mwiewior/workspace/data_vepyr/115_GRCh38_merged/variation/chr1_warm.parquet
/Users/mwiewior/workspace/data_vepyr/115_GRCh38_merged/variation/chr1_cold.parquet
```

Relevant current code paths:

- `datafusion/bio-function-vep/src/annotate_provider.rs`
- `datafusion/bio-function-vep/src/warm_cache/reader.rs`
- `datafusion/bio-function-vep/src/warm_cache/chrom_cache.rs`
- `datafusion/bio-function-vep/src/warm_cache/cold_parquet.rs`
- `datafusion/bio-function-vep/src/kv_cache/cache_exec.rs`
- `datafusion/bio-function-vep/src/variant_lookup_exec.rs`

Lance format references:

- Lance encoding strategy: <https://lance.org/format/file/encoding/>
- Lance versioning: <https://lance.org/format/file/versioning/>
- Lance table format: <https://lance.org/format/table/>
- LanceDB scalar indexes: <https://docs.lancedb.com/indexing/scalar-index>

## Chr1 Profile Summary

```text
chr1_warm.parquet
  rows:       3,628,123
  size:       338 MiB
  row groups: 8
  row-group rows: mostly 500,000

chr1_cold.parquet
  rows:       84,525,843
  size:       3.3 GiB
  row groups: 10,319
  row-group rows: mostly 8,192
```

Position-key structure:

```text
warm unique position_key values: 2,655,552
warm rows per position:          1.3662
warm multi-row positions:        21.0425%
warm max rows for one position:  50

cold unique position_key values: 80,564,179
cold rows per position:          1.0492
cold multi-row positions:        4.4287%
cold max rows for one position:  37
```

Warm `variant_keys`:

```text
rows:        3,628,123
avg length:  1.635 keys/list
max length:  99 keys/list
```

Whole-cache metadata checks across non-empty variation Parquet files:

```text
minor_allele       all null across 1,170,699,612 rows
minor_allele_freq  all null across 1,170,699,612 rows

clinical_impact    not globally null, but 99.999864% null
strand             not globally null, present only in other_cold.parquet
```

All-null ID/helper columns across the same metadata pass:

```text
assembly_ids
gencode_ids
genebuild_ids
gnomade_ids
gnomadg_ids
polyphen_ids
refseq_ids
regbuild_ids
sift_ids
src_1000genomes_ids
```

## Design Decision

Use one logical Lance variation dataset with tier-aware physical fragments:

```text
variation.lance/
  tier = warm fragments
  tier = cold fragments
```

The logical schema is common. Warm and cold fragments may use different field metadata, struct packing, fragment sizing, and index choices.

If the Lance writer/runtime cannot enforce tier-specific encoding reliably inside one dataset, keep the same logical design but materialize two datasets:

```text
variation_warm.lance
variation_cold.lance
```

That fallback preserves the key design properties while avoiding a forced compromise between scan and point-lookup layouts.

## Non-Goals

- Do not change VEP annotation semantics.
- Do not remove logical fields solely because they are null in chr1.
- Do not remove logical fields that are all-null across the cache; use constant-null encoding instead.
- Do not depend on Lance 2.2-only features for the stable production baseline.
- Do not replace SIFT/PolyPhen translation storage in this design.
- Do not remove the existing sidecar variant Bloom or position index until Lance indexes benchmark faster or simpler for the same workload.

## Access Paths

### Warm Full Scan

Warm fragments are optimized for scanning the full `everything` variation projection:

```text
position_key
variant_keys
allele_string
end
failed
variation_name
somatic
phenotype_or_disease
clin_sig
clin_sig_allele
clinical_impact
pubmed
AF and all population frequency columns
clinvar_ids
cosmic_ids
dbsnp_ids
minor_allele
minor_allele_freq
```

Scan-time priorities:

- Minimize decompression CPU for dense fields.
- Keep fields that are read together near each other in the physical plan.
- Preserve enough columnar independence to avoid pulling unused source metadata.

### Cold Batched Point Lookup

Cold fragments are optimized for batched equality probes on `position_key`:

```text
batch size target: 2,000 position_key values
primary filter:    position_key IN (...)
secondary match:   allele_string, end, failed
payload:           everything colocated fields
```

Lookup-time priorities:

- Keep `position_key` top-level and independently indexed.
- Avoid packing `position_key` into row-major payload structs.
- Preserve sorted physical order by `position_key`.
- Keep match-critical payload close enough to reduce I/O per candidate row.

## Fragment Layout

### Warm Fragments

Warm fragments should be written by chromosome and tier, sorted by:

```text
position_key, start, end, allele_string
```

Target fragment sizing:

```text
initial target: 500,000 rows per fragment
benchmark range: 250,000 to 1,000,000 rows
```

Rationale:

- Current warm Parquet uses roughly 500,000-row row groups.
- Warm full scans amortize initialization and benefit from larger contiguous reads.
- Very small fragments would increase metadata overhead without helping full scans.

### Cold Fragments

Cold fragments should be written by chromosome and tier, sorted by:

```text
position_key, start, end, allele_string
```

Target fragment sizing:

```text
initial target: 8,192 rows per fragment-equivalent unit
benchmark range: 4,096 to 16,384 rows
```

Position boundary invariant:

```text
Do not split one position_key across cold fragment-equivalent lookup units.
```

Rationale:

- Current cold Parquet uses 8,192-row row groups and already shows tight point-lookup locality.
- chr1 cold has only 1.049 rows per unique position on average.
- Batched point lookup should narrow quickly from a `position_key` BTREE lookup to candidate rows.

## Indexing

### Required

Create a scalar `BTREE` index on:

```text
position_key
```

Use it for:

- `position_key = value`
- `position_key IN (...)`
- `position_key BETWEEN min AND max` for prefetch or range batching

### Benchmark-Gated

Evaluate a scalar or label-list index on:

```text
variant_keys
```

Use only if it outperforms the current sidecar variant Bloom plus `position_key` lookup flow. Until then, keep the sidecar Bloom as the cold negative-filter path.

Do not index dense, high-cardinality text fields such as:

```text
variation_name
dbsnp_ids
```

The `everything` path reads them as payload fields, not as filters.

## Lance 2.1 Baseline Encoding

Lance 2.1 is the stable baseline. It supports improved integer and string compression and random access with nested fields. Packed struct use is limited to fixed-width children.

Use these general policies:

- Fixed-width keys and coordinates: top-level, bitpacking where useful, optional LZ4.
- Sparse flags: RLE/bitpacking, fixed packed structs where fields are read together.
- Strings: separate top-level fields with FSST where applicable plus LZ4 for warm and zstd for cold.
- All-null logical fields: constant-null layout.
- Lists: Lance list structural encoding, bitpacking for integer list values.

Compression defaults:

```text
warm: LZ4 or zstd level 1 for payload strings, chosen by benchmark
cold: zstd level 3 for sparse/string payload, LZ4 or none for position_key
```

## Lance 2.2 Enhanced Encoding

Lance 2.2 is the enhanced benchmark target. It adds useful capabilities for this dataset:

- Variable packed structs.
- Broader general block compression.
- RLE block encoding.
- Generalized constant layout.

Use 2.2 variable packed structs for groups that the `everything` path normally reads together:

```text
match_payload
identity_text
clinical_payload
af_1kg
af_gnomade
af_gnomadg
```

Do not pack `position_key` into these structs. It must remain top-level for indexing and lookup.

## Physical Field Groups

### Top-Level Key Fields

```text
position_key
start
end
chrom
```

`position_key` remains top-level in all layouts.

`chrom` is constant within chromosome fragments. In a combined dataset it should still be present logically, but encoded as a constant scalar in each fragment whenever possible.

`start` and `end` remain top-level in the default design because they are useful for validation and matching. In Lance 2.1, `end` may also be included in a fixed-width packed match struct if the runtime can still access it efficiently.

### Match Payload

Logical fields:

```text
allele_string
end
failed
```

2.1:

- `allele_string` remains separate and uses FSST plus LZ4 or zstd.
- `end` and `failed` may be packed together as fixed-width `match_fixed`.

2.2:

- Pack as `match_payload`.
- Use variable packed struct encoding.
- Compress strings with FSST plus general compression.

### Identity Text

Logical fields:

```text
variation_name
dbsnp_ids
```

2.1:

- Keep separate top-level strings.
- Avoid dictionary encoding because both are very high cardinality.
- Use FSST plus LZ4 for warm and FSST plus zstd for cold.

2.2:

- Pack into `identity_text`.
- Use variable packed struct encoding.
- Keep independent compression per child.

### Clinical Payload

Logical fields:

```text
clin_sig
clin_sig_allele
clinical_impact
pubmed
clinvar_ids
cosmic_ids
```

2.1:

- Keep separate top-level sparse strings.
- Use sparse nullable structure, FSST, and zstd for cold.

2.2:

- Pack into `clinical_payload`.
- Use variable packed struct encoding.
- Preserve `clinical_impact` as sparse, not constant-null, because it is not globally null.

### Variant Flags

Logical fields:

```text
failed
somatic
phenotype_or_disease
strand
```

2.1:

- Pack fixed-width fields where they are read together.
- Use RLE/bitpacking.

2.2:

- Same fixed-width strategy.
- Use 2.2 RLE block encoding where beneficial.

`strand` is nearly always null and appears only in `other_cold.parquet`, so chromosome fragments should use constant-null encoding when possible and sparse fixed-width encoding where values exist.

### Frequency Payloads

Split frequency fields into three physical groups.

`af_1kg`:

```text
AF
AFR
AMR
EAS
EUR
SAS
```

`af_gnomade`:

```text
gnomADe
gnomADe_AFR
gnomADe_AMR
gnomADe_ASJ
gnomADe_EAS
gnomADe_FIN
gnomADe_MID
gnomADe_NFE
gnomADe_REMAINING
gnomADe_SAS
```

`af_gnomadg`:

```text
gnomADg
gnomADg_AFR
gnomADg_AMI
gnomADg_AMR
gnomADg_ASJ
gnomADg_EAS
gnomADg_FIN
gnomADg_MID
gnomADg_NFE
gnomADg_REMAINING
gnomADg_SAS
```

2.1:

- Keep fields separate.
- Use FSST for string values.
- Use LZ4 or zstd level 1 for warm.
- Use zstd level 3 for cold.

2.2:

- Use one variable packed struct per frequency group.
- Keep child-level string compression.
- Benchmark `af_gnomadg` with and without packing because it is dense and large.

Rationale:

- 1000G fields are about 34.5% present in chr1 warm and 6.1% present in chr1 cold.
- gnomAD exome fields are about 4.4-5.7% present in chr1 warm and about 4.8-5.4% present in chr1 cold.
- gnomAD genome fields are about 76.1-76.6% present in chr1 warm and about 53.7% present in chr1 cold.

### Always-Null Logical Fields

Logical fields:

```text
minor_allele
minor_allele_freq
assembly_ids
gencode_ids
genebuild_ids
gnomade_ids
gnomadg_ids
polyphen_ids
refseq_ids
regbuild_ids
sift_ids
src_1000genomes_ids
```

Policy:

- Preserve the logical fields.
- Encode as constant-null.
- Do not materialize value buffers.
- Validate all-null status during cache build.

### Sparse Preserve-Only Fields

Logical fields:

```text
var_synonyms
hgmd_public_ids
```

Policy:

- Preserve as top-level sparse strings.
- Exclude from the `everything` hot projection.
- Use FSST plus zstd for cold and FSST plus LZ4 or zstd level 1 for warm.

### Constant Metadata Fields

Logical fields:

```text
species
assembly
cache_version
serializer_type
source_cache_path
source_file
source_assembly
source_clinvar
source_cosmic
source_dbsnp
source_gencode
source_genebuild
source_gnomade
source_gnomadg
source_hgmd_public
source_polyphen
source_refseq
source_regbuild
source_sift
source_src_1000genomes
```

Policy:

- Preserve as logical fields.
- Encode as constant scalars per fragment when possible.
- Exclude from the `everything` hot projection.
- Do not pack into hot row-major structs.

## Field Encoding Matrix

| Logical field(s) | Warm 2.1 | Warm 2.2 | Cold 2.1 | Cold 2.2 |
|---|---|---|---|---|
| `position_key` | top-level, sorted, bitpack, LZ4 or none | same | top-level, sorted, BTREE, smaller miniblock benchmark | same |
| `start` | top-level, bitpack, LZ4 | same | top-level, bitpack, zstd or LZ4 | same |
| `end` | top-level or fixed `match_fixed`, bitpack | in `match_payload` candidate if benchmark wins | top-level or fixed `match_fixed`, bitpack | in `match_payload` candidate if benchmark wins |
| `chrom` | constant scalar per chromosome fragment | same | constant scalar per chromosome fragment | same |
| `allele_string` | top-level, FSST, LZ4/zstd-1 | `match_payload`, FSST, LZ4/zstd-1 | top-level, FSST, zstd-3 | `match_payload`, FSST, zstd-3 |
| `variant_keys` | list structural, bitpack values, LZ4 | same | optional if replacing sidecar Bloom, list structural, zstd | same plus optional label-list benchmark |
| `failed` | sparse fixed, RLE/bitpack, optional fixed pack | `match_payload` or `variant_flags` | sparse fixed, RLE/bitpack | `match_payload` or `variant_flags` |
| `somatic` | `variant_flags`, RLE/bitpack | same | `variant_flags`, RLE/bitpack | 2.2 RLE block candidate |
| `phenotype_or_disease` | `variant_flags`, RLE/bitpack | same | `variant_flags`, RLE/bitpack | 2.2 RLE block candidate |
| `strand` | constant-null for normal chromosome fragments | same | constant-null except non-null fragments | sparse fixed/RLE where present |
| `variation_name` | top-level, FSST, no dictionary | `identity_text` | top-level, FSST, zstd-3, no dictionary | `identity_text`, zstd-3 |
| `dbsnp_ids` | top-level, FSST, no dictionary | `identity_text` | top-level, FSST, zstd-3, no dictionary | `identity_text`, zstd-3 |
| `clin_sig` | sparse string, FSST | `clinical_payload` | sparse string, FSST, zstd-3 | `clinical_payload`, zstd-3 |
| `clin_sig_allele` | sparse string, FSST | `clinical_payload` | sparse string, FSST, zstd-3 | `clinical_payload`, zstd-3 |
| `clinical_impact` | sparse string, FSST, not constant-null | `clinical_payload` | sparse string, FSST, zstd-3 | `clinical_payload`, zstd-3 |
| `pubmed` | sparse string, FSST | `clinical_payload` | sparse string, FSST, zstd-3 | `clinical_payload`, zstd-3 |
| `clinvar_ids` | sparse string, FSST | `clinical_payload` | sparse string, FSST, zstd-3 | `clinical_payload`, zstd-3 |
| `cosmic_ids` | sparse string, FSST | `clinical_payload` | sparse string, FSST, zstd-3 | `clinical_payload`, zstd-3 |
| `AF` | top-level sparse string, FSST | `af_1kg` | top-level sparse string, FSST, zstd-3 | `af_1kg`, zstd-3 |
| `AFR` | top-level sparse string, FSST | `af_1kg` | top-level sparse string, FSST, zstd-3 | `af_1kg`, zstd-3 |
| `AMR` | top-level sparse string, FSST | `af_1kg` | top-level sparse string, FSST, zstd-3 | `af_1kg`, zstd-3 |
| `EAS` | top-level sparse string, FSST | `af_1kg` | top-level sparse string, FSST, zstd-3 | `af_1kg`, zstd-3 |
| `EUR` | top-level sparse string, FSST | `af_1kg` | top-level sparse string, FSST, zstd-3 | `af_1kg`, zstd-3 |
| `SAS` | top-level sparse string, FSST | `af_1kg` | top-level sparse string, FSST, zstd-3 | `af_1kg`, zstd-3 |
| `gnomADe*` | top-level sparse strings, FSST | `af_gnomade` | top-level sparse strings, FSST, zstd-3 | `af_gnomade`, zstd-3 |
| `gnomADg*` | top-level strings, FSST, LZ4/zstd-1 | `af_gnomadg`, benchmark packed vs unpacked | top-level strings, FSST, zstd-3 | `af_gnomadg`, benchmark packed vs unpacked |
| `minor_allele` | constant-null | constant-null | constant-null | constant-null |
| `minor_allele_freq` | constant-null | constant-null | constant-null | constant-null |
| all-null ID/helper columns | constant-null | constant-null | constant-null | constant-null |
| `var_synonyms`, `hgmd_public_ids` | preserve-only sparse strings | same | preserve-only sparse strings, zstd-3 | same |
| constant source metadata | constant scalar, excluded from hot projection | same | constant scalar, excluded from hot projection | same |

`gnomADe*` means all fields in `af_gnomade`.

`gnomADg*` means all fields in `af_gnomadg`.

## Runtime Reconstruction

The runtime should expose the same logical columns regardless of physical packing.

For Lance 2.1:

- Most logical fields are direct top-level Lance fields.
- Fixed packed structs may require a projection adapter for child fields.

For Lance 2.2:

- Variable packed structs require a projection adapter that maps logical column names to packed child paths.
- The adapter must allow existing annotation code to request logical names such as `gnomADg_NFE` without knowing whether they live inside `af_gnomadg`.

The adapter should be isolated behind a small variation-cache reader interface:

```text
open_variation_lance(chrom, tier, logical_columns) -> RecordBatch stream
lookup_variation_lance(chrom, position_keys, logical_columns) -> RecordBatch stream
```

## Validation

Cache build validation must check:

- Logical schema contains all required fields.
- `position_key` is non-null.
- Cold fragments are sorted by `position_key`.
- Cold fragments do not split one `position_key` across fragment-equivalent lookup units.
- All-null fields configured as constant-null are actually all-null in builder input.
- `clinical_impact` is not treated as globally all-null.
- `minor_allele` and `minor_allele_freq` remain present in the logical schema.
- 2.2-only variable packed structs are not required to read 2.1 artifacts.

## Benchmark Plan

Benchmark warm and cold independently.

Warm full-scan benchmark:

```text
projection: everything variation columns
metric: rows/s, bytes read, decode CPU, wall time
compare:
  current Parquet warm
  Lance 2.1 warm
  Lance 2.2 warm unpacked
  Lance 2.2 warm packed
```

Cold batched point-lookup benchmark:

```text
batch size: 2,000 position_key values
projection: everything variation columns
metric: p50/p95 latency, rows/s, bytes read, I/O operations, decode CPU
compare:
  current Parquet cold + side indexes
  Lance 2.1 cold + position_key BTREE
  Lance 2.2 cold + position_key BTREE
  optional Lance variant_keys index
```

Required acceptance gates:

- Warm Lance must not regress full-scan throughput by more than 5% unless it reduces storage or cold lookup enough to justify the trade-off.
- Cold Lance must improve or match 2,000-key batched lookup latency against current indexed Parquet.
- 2.2 layout must beat 2.1 on either warm full scan or cold lookup before becoming the default enhanced artifact.
- VCF output body hash must match the current `indexed_parquet` path for the same input and options.

## Risks

### Lance 2.2 Stability

Lance 2.2 is an enhanced target, not the only production path. The stable 2.1 layout must remain readable and benchmarked.

### Packed Struct Overuse

Packed structs reduce random-access I/O when all children are read together, but they can hurt selective reads. This design packs only groups that the `everything` path reads together and keeps source metadata outside hot packed groups.

### Dense gnomAD Genome Fields

`gnomADg*` fields are dense and large, especially in warm. Packing may help point lookup but can hurt scans or partial projections. Benchmark packed and unpacked 2.2 variants.

### Index Replacement Risk

The current sidecar Bloom and position index are known to match the Parquet layout. Lance `variant_keys` indexing should be optional until it proves better for negative filtering.

## Open Implementation Choices

These are implementation choices to settle during planning, not design ambiguities:

- Whether a single Lance dataset can cleanly hold tier-specific fragment encodings.
- Exact warm fragment size after benchmarking 250k, 500k, and 1M rows.
- Exact cold fragment-equivalent size after benchmarking 4,096, 8,192, and 16,384 rows.
- Whether `gnomADg*` should remain unpacked in 2.2 if full-scan throughput is better.
- Whether `variant_keys` should be materialized in cold Lance or remain sidecar-only.
