# Tiered Warm Parquet Cache Design

## Goal

Build a tiered variation lookup path where Parquet remains the master VEPyR cache representation, a small AF-selected warm Parquet tier accelerates common variant probes, and Fjall remains the correctness-preserving fallback for probes not covered by the warm tier.

## Motivation

Profiling showed the variation lookup path spends most of its time in Fjall point lookups. A warm tier selected by allele frequency can remove most of those point lookups while preserving the existing Fjall backend for cold data and correctness fallback.

Measured on the full benchmark VCF over chromosomes 1-22 with selector `max_global_af >= 1%` plus `+/-1` positions:

| Metric | Value |
|---|---:|
| Source variation rows | 1,112,354,246 |
| Warm rows | 46,996,195 |
| Warm row share | 4.22% |
| Source Parquet size | 36,501.7 MiB |
| Warm Parquet size | 3,536.6 MiB |
| Warm size share | 9.69% |
| Total lookup probes | 8,999,356 |
| Warm-covered probes | 6,227,473 |
| Fjall fallback probes | 2,771,883 |
| Fjall probes avoided | 69.20% |

## Non-Goals

- Do not replace Fjall in the first implementation.
- Do not use cold Parquet for runtime lookup in the first implementation.
- Do not load full chromosome warm Parquet files into memory.
- Do not load all 78 variation columns for each active warm chunk.
- Do not change transcript/regulatory context loading in this change.

## Cache Files

For each chromosome, write two master Parquet files:

```text
variation/
  chr1_warm.parquet
  chr1_cold.parquet
  chr2_warm.parquet
  chr2_cold.parquet
  ...
```

The original `chrN.parquet` may remain for compatibility, but the tiered runtime reads only `chrN_warm.parquet` for acceleration and uses `variation.fjall` for fallback.

Warm/cold split is based on complete positions, not individual rows:

```text
hot position = max_global_af >= 0.01
warm position = hot position OR hot_position - 1 OR hot_position + 1
warm parquet = all rows whose start is a warm position
cold parquet = all remaining rows
```

`max_global_af` is:

```text
max(minor_allele_freq, AF, gnomADg, gnomADe)
```

The string frequency columns contain allele-frequency pairs such as `A:0.1,C:2.6e-05`; the selector uses the maximum numeric frequency in each field.

## Key Columns

Add generated key columns to `_warm` and `_cold` Parquet files.

### `position_key`

Purpose: completeness and fallback routing.

```text
position_key = numeric representation of chrom_code + start
```

It has the same semantics as the current Fjall position key. Runtime logic:

```text
position_key absent from active warm chunk -> Fjall fallback
position_key present but variant_key misses -> definitive warm miss, no Fjall
```

### `variant_keys`

Purpose: fast exact candidate lookup.

```text
variant_keys = List<UInt64>
```

The list contains one deterministic key for each alternate allele represented by the cache row's `allele_string`. Multi-allelic cache rows therefore keep one physical Parquet row and expose several lookup keys.

For short A/C/G/T alleles, use an echtvar-style compact packed key. For long, symbolic, or unsupported alleles, use a stable 64-bit hash over the normalized tuple:

```text
chrom_code, start, ref_allele, alt_allele
```

All hits must still verify `allele_string`, so hash collisions cannot create false positives.

### Hashing and HashMap Choice

Persistent `variant_keys` stored in Parquet must not depend on Rust `HashMap` hasher output. Use a stable portable `rapidhash` function, preferably the current stable V3 API, to generate hashed variant keys for long, symbolic, or otherwise unpackable alleles.

Runtime warm-chunk lookup should avoid hashing those keys a second time. The preferred index type is:

```rust
type VariantIndex = hashbrown::HashMap<
    u64,
    SmallVec<[u32; 1]>,
    nohash_hasher::BuildNoHashHasher<u64>,
>;
```

This is valid only because `variant_key` is already a high-quality stable 64-bit hash or a deliberately packed key with acceptable distribution. If benchmarks show clustering for packed short-allele keys, switch only the runtime map hasher to `rapidhash::fast::RandomState`:

```rust
type VariantIndex = hashbrown::HashMap<
    u64,
    SmallVec<[u32; 1]>,
    rapidhash::fast::RandomState,
>;
```

Do not use the default standard-library `HashMap` for this hot index unless benchmarks show it is equivalent. The implementation should benchmark `nohash_hasher`, `rapidhash::fast::RandomState`, and the existing `ahash` dependency on a representative warm chunk before making the default permanent.

## Row Groups As Chunks

Warm Parquet row groups are the runtime chunks. Row groups should be sorted by:

```text
position_key, start, end, variant_keys
```

The first implementation should use `250,000` warm rows per row group and benchmark `100,000`, `250,000`, and `500,000`.

Runtime loads one active warm chunk per chromosome lane, with one adjacent chunk allowed for boundary-crossing VCF buffers.

## Runtime Data Structure

Each chromosome lane owns its warm chunk state. No mutable warm chunk state is shared across chromosome lanes.

```rust
struct WarmChunkContext {
    row_group_id: usize,
    min_position_key: u64,
    max_position_key: u64,
    position_keys: Vec<u64>,
    variant_index: HashMap<u64, SmallVec<[u32; 1]>>,
    batch: RecordBatch,
    columns: WarmColumnIndices,
}
```

`position_keys` is sorted and deduplicated. Use binary search for the completeness check.

`variant_index` maps each `variant_key` to one or more row numbers in `batch`. Multiple rows are allowed because of collisions, duplicate sources, or multi-allelic representations.

`batch` is a projected Arrow `RecordBatch` for the active row group only. It must not include source metadata columns that are not used by lookup/CSQ output.

Initial runtime projection:

```text
position_key
variant_keys
chrom
start
end
allele_string
variation_name
failed
somatic
strand
minor_allele
minor_allele_freq
clin_sig
phenotype_or_disease
clinical_impact
clin_sig_allele
AF
gnomADg
gnomADe
```

Additional output columns should be added only when a correctness test proves they are needed.

## Probe Flow

For each VEP-compatible probe generated from a VCF row:

```text
compute position_key
compute variant_key

if active warm chunk covers position_key:
    if variant_key exists in variant_index:
        read candidate Arrow row(s)
        verify allele_string
        if verified:
            emit warm cache values
            skip Fjall for this probe

    if position_key exists in position_keys:
        definitive warm miss
        skip Fjall for this probe

fallback:
    call existing Fjall position lookup
```

Worst-case cold probe cost is two cheap in-memory checks plus the existing Fjall lookup. Warm probes avoid Fjall entirely.

## Expected Improvement

The AF>=1% plus `+/-1` selector should reduce Fjall point lookups by about 69% on the 4M-variant benchmark VCF.

Expected runtime impact:

| Scope | Expected improvement |
|---|---:|
| Fjall point lookup count | 65-70% fewer |
| Variation lookup stage | 1.8x-3.0x faster |
| End-to-end VCF annotation | 1.2x-1.6x faster |

For a current 49-57 second full benchmark run, a realistic first target is 35-45 seconds if warm chunk indexing overhead stays low. The minimum acceptance target is at least 15% end-to-end improvement with 0 mismatches and max RSS under 6 GB.

## Correctness Requirements

- Warm split must be complete by position.
- Warm hits must verify `allele_string`; key match alone is not sufficient.
- Warm definitive misses may skip Fjall only when `position_key` exists in the active warm chunk.
- If a VCF buffer crosses a row-group boundary, the runtime must either load the adjacent warm chunk or fallback to Fjall for positions outside the active chunk.
- Existing Fjall output must remain the fallback source of truth.
- Merged, Ensembl, and RefSeq caches must all validate with 0 mismatches before enabling by default.

## Metrics

Add profiling counters:

```text
warm_chunks_loaded
warm_chunk_rows
warm_chunk_load_s
warm_chunk_index_s
warm_variant_key_hits
warm_variant_key_misses
warm_verified_hits
warm_definitive_misses
warm_position_misses
fjall_fallbacks
fjall_fallback_saved
```

Report these under existing `VEP_KV_PROFILE=1` output.

## Rollout

Phase 1 is opt-in only:

```text
VEP_WARM_VARIATION_CACHE=1
```

or an equivalent VEPyR option after the Rust implementation is stable.

Auto-detection should require:

```text
variation/chrN_warm.parquet
variation.fjall/
```

If warm files are missing or malformed, the system must silently fall back to current Fjall-only behavior unless profiling/debug mode is enabled.
