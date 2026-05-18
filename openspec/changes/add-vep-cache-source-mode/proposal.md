# Add Explicit VEP Cache Source Mode

## Why

`annotate_vep()` currently models Ensembl-vs-merged behavior with an optional
`merged` boolean in `options_json`. That is not expressive enough for
RefSeq-only caches and allows ambiguous execution when a cache table does not
declare whether it represents Ensembl, merged, or RefSeq data.

The companion `datafusion-bio-formats` change will expose
`bio.vep.cache_source_type` schema metadata from an explicit
`cache_source_type` provider/export option. `bio-function-vep` must consume that
metadata and apply Ensembl VEP-compatible source-specific transcript semantics.

## What Changes

- Add an explicit `CacheSourceType`/source-mode contract with valid values
  `ensembl`, `merged`, and `refseq`.
- Require cache table schema metadata `bio.vep.cache_source_type`; missing or
  invalid values are planning/execution errors.
- Remove unsupported `merged = true` compatibility handling from `options_json`.
- Replace boolean merged transcript filtering with source-mode-aware filtering.
- Support RefSeq-only transcript IDs, including Ensembl VEP's mitochondrial
  RefSeq cache exceptions.
- Populate CSQ `SOURCE` only for merged cache mode, not RefSeq-only mode.
- Extend RefSeq hydration helpers so RefSeq-only caches get the same transcript
  sequence handling as merged RefSeq rows.

## Impact

### Affected Specs
- **NEW**: `vep-cache-source-mode` — source-mode contract for VEP annotation

### Affected Code
- `datafusion/bio-function-vep/src/config.rs` — source-mode enum and metadata key constants if kept in config
- `datafusion/bio-function-vep/src/annotate_table_function.rs` — validate cache source metadata and reject `merged`
- `datafusion/bio-function-vep/src/annotate_provider.rs` — store source mode, use it for transcript filtering and CSQ `SOURCE`
- `datafusion/bio-function-vep/src/transcript_consequence.rs` — replace `is_vep_transcript(id, merged)` with source-mode-aware filtering
- `datafusion/bio-function-vep/src/schema_contract.rs` or a new `cache_source.rs` — shared validation helpers if extraction should be isolated
- `datafusion/bio-function-vep/tests/` or existing unit modules — source-mode, RefSeq, and rejected-option tests

### Breaking Changes

`options_json` no longer accepts `{"merged": true}`. Callers must register or
export cache tables with schema metadata
`bio.vep.cache_source_type = "merged"`.
