# Indexed Parquet Cache Format Design

## Goal

Make `indexed_parquet` the default cache format for VEP annotation in both this repository and VEPyR. The format must run annotation without touching Fjall for variation or translation SIFT/PolyPhen lookup.

The legacy Fjall path remains available only through an explicit `cache_format=legacy_fjall` option.

## Motivation

The current fast path was assembled from several opt-in pieces: warm/cold variation Parquet, cold Parquet row-group/page pruning, variant Bloom sidecars, and compact Parquet SIFT lookup. Benchmarks showed that this path can match the Fjall baseline on chr4 with `forks=1` while producing identical VCF body hashes:

```text
Fjall baseline                              19.09s, 16,096 rows/s
Indexed Parquet variation + compact SIFT    19.06s, 16,125 rows/s
VCF body mismatches                         0
```

The implementation should turn that path into one coherent cache format rather than a set of environment-variable switches.

## Non-Goals

- Do not keep `use_fjall` as a public API flag.
- Do not keep the old VEP Parquet interval-join lookup path as a public cache format.
- Do not remove the generic `IntervalJoinExec` implementation from `datafusion-bio-function-ranges`; it is a separate ranges feature.
- Do not require environment variables for normal annotation.
- Do not add or use `region_bin`.
- Do not add or use `vepyr.echtvar_bin_bits`.
- Do not generate Fjall stores by default.

## Public API

### Rust

Introduce a `CacheFormat` enum with these public values:

```text
indexed_parquet
legacy_fjall
```

Default:

```text
cache_format = indexed_parquet
```

The current `backend=parquet|fjall` and `use_fjall` naming is misleading because the fast path currently enters the `fjall` backend name to reach `KvLookupExec`. Replace that naming at public call sites with `cache_format`.

The old VEP annotation `backend=parquet` route, which uses the Parquet table lookup path rather than the indexed warm/cold cache, should be removed from public annotation APIs. It may remain temporarily only as an internal test/developer compatibility path during migration.

### VEPyR

Expose the same parameter in Python and pyo3:

```python
cache_format: Literal["indexed_parquet", "legacy_fjall"] = "indexed_parquet"
```

Apply it consistently to cache generation and annotation. `forks > 0` must work with `indexed_parquet`; it must not force `legacy_fjall`.

### Legacy Compatibility

`legacy_fjall` means:

- generate and use `variation.fjall`
- generate and use `translation_sift.fjall`
- preserve the existing Fjall lookup semantics

`indexed_parquet` means:

- do not generate `variation.fjall`
- do not generate `translation_sift.fjall`
- annotation fails clearly if required indexed Parquet artifacts are missing
- annotation does not silently fall back to Fjall

## Indexed Parquet Cache Layout

The cache root for `indexed_parquet` must contain:

```text
variation/
  chr1_warm.parquet
  chr1_cold.parquet
  ...
variation.position_index/
  chr1.posidx
  ...
variation.variant_bloom_index/
  chr1.varbf
  ...
translation_core/
  chr1.parquet
  ...
translation_sift/
  chr1.parquet
  ...
```

`translation_sift/` keeps the old public directory name, but its schema changes under `indexed_parquet` to the compact lookup schema described below.

## Variation Schema

Both `chrN_warm.parquet` and `chrN_cold.parquet` use the same schema family.

Required source/cache columns:

```text
chrom: Utf8-compatible, non-null
start: Int64, non-null
end: Int64, non-null
variation_name: Utf8-compatible, nullable
allele_string: Utf8-compatible, non-null
failed: Int64-compatible, nullable or non-null
somatic: Int64-compatible, nullable
phenotype_or_disease: Int64-compatible, nullable
clin_sig: Utf8-compatible, nullable
clin_sig_allele: Utf8-compatible, nullable
pubmed: Utf8-compatible, nullable
AF and all AF_COL_NAMES columns requested by colocated output: Utf8-compatible, nullable
```

Required generated key columns:

```text
position_key: UInt64, non-null
variant_keys: List<UInt64>, non-null
```

The files must preserve every normal variation cache output column needed by annotation and requested `cache_columns`. The writer may drop internal helper columns that are not annotation outputs.

Explicitly forbidden columns and metadata:

```text
region_bin
vepyr.echtvar_bin_bits
```

## Variation Row-Group Layout

Warm variation:

```text
target row group rows: 500,000
sort order: position_key, start, end, allele_string
```

Cold variation:

```text
target row group rows: 8,192
data page row count target: 1,024
sort order: position_key, start, end, allele_string
```

Hard invariant for both warm and cold files:

```text
A single position_key must never be split across row groups.
```

The writer may exceed the target row group size if all pending rows share the same `position_key`. It must flush only at a `position_key` boundary.

Validation must fail if adjacent row groups share the same `position_key`:

```text
row_group[i].max(position_key) == row_group[i + 1].min(position_key)
```

This invariant replaces the earlier bin-boundary experiment. The layout is not a `bin20` layout.

## Variation Page Index and Statistics

`chrN_cold.parquet` must be written with:

```text
position_key min/max statistics
position_key column index
position_key offset index
```

Runtime cold lookup must use:

- row-group pruning from `position_key` min/max
- page pruning from the `position_key` page index when available
- row selection for the candidate page/range

The cache validator must fail `indexed_parquet` cold variation files when the `position_key` statistics or page indexes are missing.

Warm files must have `position_key` min/max statistics. Page indexes on warm files are allowed but not required for the first implementation because warm lookup treats row groups as chunks.

## Variation Bloom Filter

The required Bloom filter for variation is the sidecar:

```text
variation.variant_bloom_index/chrN.varbf
```

It is a Bloom filter over `variant_keys`, not `position_key`.

Runtime uses it before reading cold Parquet:

```text
warm not covered
  -> position index says cold may contain position
  -> variant Bloom says none of the candidate variant keys are present
  -> skip cold Parquet read
```

If the variant Bloom file is missing or unreadable under `cache_format=indexed_parquet`, annotation must return a clear cache-format error rather than falling back to Fjall.

Parquet column-level Bloom filters are optional in this design. The runtime must not depend on them for correctness or performance claims until a separate benchmark proves value.

## Variation Position Index

The required position index is:

```text
variation.position_index/chrN.posidx
```

It stores cold `position_key` presence. Runtime uses it to skip cold Parquet for positions not present in the cold file.

If the position index is missing under `cache_format=indexed_parquet`, annotation must return a clear cache-format error. The old behavior of generating a fallback index by scanning Parquet may remain as an explicit developer/debug path, but not as the normal default.

## Compact Translation SIFT Schema

Under `cache_format=indexed_parquet`, `translation_sift/` is the compact lookup table, not the old 6-column split table.

Required schema:

```text
transcript_id: Utf8-compatible, non-null
predictions: Binary, non-null
```

`predictions` uses the same serialized byte representation as the current `translation_sift.fjall` SIFT keyspace. It contains the sorted compact SIFT and PolyPhen predictions for one transcript.

Required logical constraints:

- one row per resolved `transcript_id`
- duplicate transcript rows resolved with the same semantics as the current `translation_sift.fjall` builder
- canonical chromosomes preferred over patch/alt sources
- when source class ties, keep the richer prediction matrix

The old full SIFT columns are not present in `indexed_parquet` `translation_sift/`:

```text
chrom
start
end
sift_predictions
polyphen_predictions
```

Those columns may continue to exist only in intermediate builder state or in `legacy_fjall`.

## Compact Translation SIFT Layout

Required layout:

```text
path: translation_sift/
target row group rows: 16
compression: uncompressed by default
sort order: transcript_id
```

The row group size of 16 is part of the default format because compact SIFT lookup reads sparse transcript ids. Runtime reader batch size is separate and may remain 8,192.

Page index:

- write `transcript_id` statistics and page index when supported by the Parquet writer
- do not require runtime to use the page index in the first implementation, because `SiftParquetStore` builds an in-memory exact transcript-id to row-ref index at open time

Bloom filters:

- no separate SIFT Bloom sidecar is required in the first implementation
- an exact in-memory `HashMap<transcript_id, row refs>` is the lookup index
- Parquet column-level Bloom filters on `transcript_id` are optional and must not be used as a required correctness dependency

Validation must fail if `translation_sift/` is missing or if any file lacks `transcript_id` or `predictions`.

## Runtime Data Flow

For variation:

```text
VCF row
  -> build VEP-compatible probe starts
  -> warm Parquet probe
  -> cold position index check
  -> variant Bloom check
  -> cold Parquet row-group/page read
  -> allele verification
  -> emit cache columns or nulls
```

For translation SIFT/PolyPhen:

```text
annotation buffer
  -> collect transcript ids
  -> SiftParquetStore::get_many(transcript_ids)
  -> grouped row-group reads from translation_sift/
  -> deserialize predictions
  -> annotate transcript consequences
```

For `cache_format=indexed_parquet`, neither flow may open or query Fjall.

## Environment Variable Cleanup

Normal annotation must not require these environment variables:

```text
VEP_WARM_VARIATION_CACHE
VEP_WARM_COLD_VARIATION_BACKEND
VEP_WARM_COLD_VARIATION_INDEX_MODE
VEP_WARM_VARIATION_DIR
VEP_VARIATION_COLD_DIR
VEP_VARIATION_POSITION_INDEX_DIR
VEP_VARIATION_BLOOM_INDEX_DIR
VEP_SIFT_LOOKUP_PARQUET_DIR
VEP_COLD_PARQUET_LOAD_PAGE_INDEX
```

The implementation should replace normal-use env routing with `cache_format` and cache-root discovery.

Developer-only overrides may remain temporarily for benchmarks, but they must not be needed in VEPyR or documented as normal operation.

## Cache Builder Responsibilities

For `cache_format=indexed_parquet`, the cache builder must:

1. Build `translation_core/`.
2. Build compact `translation_sift/` with schema `transcript_id, predictions`.
3. Build warm/cold variation Parquet.
4. Drop `region_bin` if present in source input.
5. Add `position_key` and `variant_keys`.
6. Write cold variation with target 8,192 row groups and 1,024-row data pages.
7. Preserve the no-split `position_key` row-group invariant.
8. Write `variation.position_index/`.
9. Write `variation.variant_bloom_index/`.
10. Validate required indexed artifacts before reporting success.
11. Skip `variation.fjall` and `translation_sift.fjall`.

For `cache_format=legacy_fjall`, the cache builder may keep the existing legacy outputs and Fjall stores.

## Implementation References

The implementation should reuse existing library code rather than shelling out to examples.

Variation builder/library references:

- `datafusion/bio-function-vep/src/warm_cache/build.rs`
  - current warm/cold variation split
  - current `PositionAlignedWriter`
  - current `.posidx` generation from cold Parquet
- `datafusion/bio-function-vep/src/kv_cache/position_index.rs`
  - `PositionIndex::from_parquet`
  - `PositionIndex::from_parquet_files`
  - `PositionIndex::write_to_path`
- `datafusion/bio-function-vep/src/kv_cache/variant_bloom_index.rs`
  - `VariantBloomIndex::from_parquet`
  - `VariantBloomIndex::write_to_path`

Example prototypes to fold into library code:

- `datafusion/bio-function-vep/examples/rewrite_cold_variation_layout.rs`
  - reuse the `data_page_row_count` writer behavior
  - reuse the 8,192-row cold row-group layout behavior
  - reuse position-key-aligned row-group flushing
  - do not reuse `--bin-bits`, `genomic_bin`, or `vepyr.echtvar_bin_bits`
- `datafusion/bio-function-vep/examples/build_variation_position_index.rs`
  - keep only as a CLI wrapper over library code
- `datafusion/bio-function-vep/examples/build_variation_variant_bloom_index.rs`
  - keep only as a CLI wrapper over library code

Compact SIFT references:

- `datafusion/bio-function-vep/src/kv_cache/sift_parquet_store.rs`
  - runtime compact parquet lookup store
  - required compact schema reader: `transcript_id`, `predictions`
- `datafusion/bio-function-vep/examples/build_sift_lookup_parquet.rs`
  - fold the compact writer behavior into the cache builder
  - keep the 16-row group default
  - write to `translation_sift/`, not `translation_sift_lookup/`
- `datafusion/bio-function-vep/src/kv_cache/sift_store.rs`
  - reuse the existing Fjall-compatible prediction serialization/deserialization semantics

The examples are useful prototypes and benchmark tools. They must not be required steps in normal cache generation.

## Annotation Responsibilities

For `cache_format=indexed_parquet`, annotation must:

1. Discover all indexed Parquet artifacts from the cache root.
2. Open variation warm/cold Parquet and required side indexes.
3. Open compact `translation_sift/` through `SiftParquetStore`.
4. Reject missing indexed artifacts with clear errors.
5. Avoid opening `variation.fjall`.
6. Avoid opening `translation_sift.fjall`.
7. Keep VCF output identical to `legacy_fjall` for supported caches.

For `cache_format=legacy_fjall`, annotation may use the current Fjall stores.

## VEPyR Changes

Python API:

```python
build_cache(..., cache_format: str = "indexed_parquet", ...)
annotate(..., cache_format: str = "indexed_parquet", ...)
```

Rust pyo3 API:

```text
cache_format = "indexed_parquet"
```

Remove public `use_fjall` and `build_fjall`. Existing behavior should be expressible as:

```python
cache_format="legacy_fjall"
```

The VEPyR wrapper must stop forcing Fjall when `forks > 0` or when warm variation is enabled.

VEPyR wiring requirements:

- add `cache_format` to `src/vepyr/__init__.py`, `src/vepyr/_core.pyi`, and `src/lib.rs`
- thread `cache_format` through `src/annotate.rs` into the Rust annotation call
- thread `cache_format` through cache generation into this repository's `CacheBuilder`
- default `build_cache()` and `annotate()` to `indexed_parquet`
- remove the `forks > 0 requires use_fjall=True` guard
- remove the `warm_variation_cache=True requires use_fjall=True` guard
- remove public warm/cold artifact override parameters that are no longer normal API:
  - `warm_variation_cache`
  - `warm_variation_dir`
  - `variation_cold_dir`
  - `variation_position_index_dir`
  - `variation_bloom_index_dir`
- remove public `use_fjall`
- remove public `build_fjall`
- update Python docstrings and type stubs to describe `cache_format`
- ensure no removed parameter becomes a silent no-op

Normal VEPyR annotation should discover indexed artifacts from the cache root. Developer-only environment overrides may remain in the Rust crate during migration, but they must not be part of the VEPyR public API.

## Testing Requirements

Unit tests:

- parse `CacheFormat`
- reject unknown cache format values
- compact `translation_sift/` schema is exactly `transcript_id, predictions`
- compact SIFT row group default is 16
- variation output schema drops `region_bin`
- no `vepyr.echtvar_bin_bits` metadata is written
- cold variation row groups do not split one `position_key`
- validator detects adjacent row groups sharing a `position_key`
- validator detects missing `position_key` page index
- validator detects missing variant Bloom sidecar

Integration tests:

- build a small indexed cache with no Fjall outputs
- annotate with `cache_format=indexed_parquet`
- assert no Fjall directory is required
- assert compact SIFT predictions match legacy Fjall predictions
- assert VCF body hash matches `cache_format=legacy_fjall` on the same fixture
- run annotation with `forks=1` and `forks=0`

VEPyR smoke tests:

- `vepyr.build_cache(cache_format="indexed_parquet")` succeeds on the small fixture cache
- `vepyr.build_cache()` defaults to `indexed_parquet`
- generated cache contains:
  - `variation/chrN_warm.parquet`
  - `variation/chrN_cold.parquet`
  - `variation.position_index/chrN.posidx`
  - `variation.variant_bloom_index/chrN.varbf`
  - compact `translation_sift/chrN.parquet`
- generated cache does not contain:
  - `variation.fjall`
  - `translation_sift.fjall`
  - `translation_sift_lookup/`
- compact `translation_sift/` schema is exactly `transcript_id,predictions`
- generated variation Parquet has no `region_bin`
- generated variation Parquet has no `vepyr.echtvar_bin_bits`
- cold variation `position_key` page index is present
- adjacent cold row groups never share a `position_key`
- `vepyr.annotate(cache_format="indexed_parquet", forks=0)` succeeds
- `vepyr.annotate(cache_format="indexed_parquet", forks=1)` succeeds
- default `vepyr.annotate()` uses `indexed_parquet`
- indexed output VCF body hash matches the expected golden fixture
- if the fixture also builds `legacy_fjall`, indexed output body hash matches `legacy_fjall`
- removed parameters are rejected by Python instead of silently ignored:
  - `use_fjall`
  - `build_fjall`
  - `warm_variation_cache`
  - `warm_variation_dir`
  - `variation_cold_dir`
  - `variation_position_index_dir`
  - `variation_bloom_index_dir`

Benchmark acceptance:

- chr4 indexed Parquet runtime should remain close to the measured Fjall baseline
- generated VCF body hash must match Fjall baseline
- detailed metrics must show zero variation Fjall point gets and no translation SIFT Fjall opens

## Migration Notes

Existing caches with old full-schema `translation_sift/` are not valid `indexed_parquet` caches. They can be used through `legacy_fjall` or rebuilt.

Existing caches with `region_bin` in variation Parquet should still be readable by legacy/dev tools, but new `indexed_parquet` cache generation must omit it.

Existing code paths that discover `translation_sift_lookup/` should be migrated to `translation_sift/` for `indexed_parquet`. The indexed cache builder must not create a `translation_sift_lookup/` output.

## Open Risks

- Removing bin-boundary flushing changes row-group boundaries compared with the previous experimental `8k_page8_bin20` artifact. Because runtime never used `region_bin` or `vepyr.echtvar_bin_bits`, this should be safe, but the cleaned layout must be benchmarked before release.
- Replacing old `translation_sift/` with compact `translation_sift/` can affect tools that expected the old 6-column table. Those tools should either move to `legacy_fjall` compatibility or use an explicitly named intermediate artifact.
- Some Parquet writer APIs may not expose column-level Bloom filters uniformly. This design requires the sidecar variant Bloom and treats Parquet column Bloom filters as optional.

## Acceptance Criteria

The implementation is complete when:

```text
cache_format=indexed_parquet is the default in Rust and VEPyR
cache generation writes variation indexed Parquet and compact translation_sift/
cache generation does not write Fjall stores by default
annotation runs without Fjall for variation or SIFT
region_bin is absent from generated indexed variation Parquet
vepyr.echtvar_bin_bits is absent from generated indexed variation Parquet metadata
cold variation uses 8,192 target row groups and 1,024-row data pages
no position_key spans row-group boundaries
translation_sift/ uses transcript_id,predictions with 16-row groups
VCF body hash matches legacy_fjall on benchmark fixtures
chr4 performance is close to the established Fjall baseline
```
