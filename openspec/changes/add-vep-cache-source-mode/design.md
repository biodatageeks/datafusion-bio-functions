# VEP Cache Source Mode Design

## Context

`bio-function-vep` already supports Ensembl transcript annotation and a partial
merged-cache mode via `options_json` key `merged`. RefSeq-only caches need
different behavior: VEP accepts RefSeq stable IDs, preserves odd mitochondrial
RefSeq transcript names, applies RefSeq hydration logic, and does not populate
CSQ `SOURCE` merely because the cache is RefSeq-only.

The companion `datafusion-bio-formats` change writes cache source mode into
Arrow schema metadata from an explicit provider/export option:

```
bio.vep.cache_source_type = ensembl | merged | refseq
```

This repo consumes that metadata. It must not infer cache source mode from
paths, table names, or option booleans.

Research update:
- `datafusion-bio-formats` commit
  `ca91f694e4db145ff0f8c1118b92a7f4c6f9aee3` adds
  `CacheSourceType`, exports `VEP_CACHE_SOURCE_TYPE_METADATA_KEY`, requires
  `EnsemblCacheOptions.cache_source_type`, and writes
  `bio.vep.cache_source_type` on variation, transcript, exon, translation,
  translation split, regulatory, and motif schemas.
- Ensembl VEP release/115 derives cache source type from `--merged`/`--refseq`
  flags, rejects `homo_sapiens_refseq` and `homo_sapiens_merged` as species
  names, filters RefSeq transcripts only for RefSeq caches or merged-cache rows
  whose `_source_cache` is `RefSeq`, and exposes CSQ `SOURCE` through the
  merged flag field set only.

## Goals / Non-Goals

Goals:
- Require explicit source-mode metadata for VEP cache tables.
- Remove `merged` boolean compatibility from `options_json`.
- Match Ensembl VEP transcript inclusion behavior for Ensembl, merged, and
  RefSeq caches.
- Keep CSQ `SOURCE` output merged-only.
- Keep the change local to `bio-function-vep` semantics and tests.

Non-goals:
- Registering or exporting native cache tables. That belongs to
  `datafusion-bio-formats`.
- Inferring source mode from directory names such as `homo_sapiens_refseq`.
- Adding a new `annotate_vep()` positional argument for source mode.
- Reworking cache storage backends or parquet layout.

## Source Mode Contract

Add a small internal enum:

```rust
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum CacheSourceType {
    Ensembl,
    Merged,
    RefSeq,
}
```

Parsing is strict:

```rust
impl CacheSourceType {
    pub const METADATA_KEY: &'static str = "bio.vep.cache_source_type";

    pub fn parse(raw: &str) -> Result<Self> {
        match raw {
            "ensembl" => Ok(Self::Ensembl),
            "merged" => Ok(Self::Merged),
            "refseq" => Ok(Self::RefSeq),
            other => Err(DataFusionError::Plan(format!(
                "annotate_vep(): unsupported bio.vep.cache_source_type '{other}'; expected one of ensembl, merged, refseq"
            ))),
        }
    }

    pub fn from_schema(schema: &Schema) -> Result<Self> {
        let raw = schema.metadata().get(Self::METADATA_KEY).ok_or_else(|| {
            DataFusionError::Plan(
                "annotate_vep(): cache table schema is missing bio.vep.cache_source_type; expected one of ensembl, merged, refseq".to_string(),
            )
        })?;
        Self::parse(raw)
    }
}
```

The helper can live in a new
`datafusion/bio-function-vep/src/cache_source.rs` module to avoid adding more
state to already-large `annotate_provider.rs`.

## Validation Boundary

Validate source mode at table-function planning time where possible:

1. `AnnotateFunction::call()` rejects `options_json` containing top-level
   source selector keys `merged` or `refseq` before creating the provider.
2. For the current partitioned-cache path API, `AnnotateFunction::call()`
   synchronously reads the Arrow schema metadata from the first
   `variation/*.parquet` file under `cache_source`. This is required because
   `AnnotateProvider::schema()` must already know whether RefSeq and `SOURCE`
   output fields are present during planning.
3. It stores `CacheSourceType` on `AnnotateProvider`.

Context tables resolved later by `AnnotateProvider` should either:
- carry the same `bio.vep.cache_source_type`, or
- return a clear error if they are VEP cache tables and the metadata is missing
  or inconsistent.

The variation table source mode is authoritative for provider output schema.
Transcript and other context table metadata must match it. If any context table
metadata differs, return an error naming both values.

## Transcript Filtering

Replace boolean source-mode selection in `TranscriptSelectionFlags` and retire
the unused `transcript_consequence::is_vep_transcript(id, merged)` helper or
replace it with source-aware helper functions that operate on
`TranscriptFeature`.

Rules:

| Source mode | Accepted rows |
|---|---|
| `Ensembl` | `stable_id` starts with `ENST` |
| `RefSeq` | RefSeq stable IDs and RefSeq mitochondrial cache names |
| `Merged` | `ENST` rows plus rows whose normalized row source is `RefSeq` and match RefSeq rules |

RefSeq stable ID rules should match Ensembl VEP:
- standard RefSeq IDs: `^[A-Z]{2}_\d+`
- mitochondrial exceptions: `^\d{4}$|^(rna-)?[A-Z0-9]{3,}$`

For merged mode, RefSeq ID and mitochondrial exceptions should be accepted only
for rows whose source normalizes to `RefSeq`, so a non-RefSeq row with a
gene-like ID is not accidentally included. This matches Ensembl VEP
`AnnotationType::Transcript::filter_transcript()`, which applies the RefSeq
whitelist for `source_type eq 'refseq'` or for `source_type eq 'merged'` rows
whose `_source_cache` is `RefSeq`.

## CSQ SOURCE Output

Replace boolean `merged` output gating:

```rust
let source_val = if merged { source } else { "" };
```

with:

```rust
let source_val = if source_type == CacheSourceType::Merged {
    source
} else {
    ""
};
```

This preserves VEP's merged-cache-only `SOURCE` behavior and avoids incorrectly
emitting `RefSeq` for RefSeq-only caches.

## RefSeq Hydration

Update `is_refseq_transcript_for_hydration()` so it accepts:
- normalized row source `RefSeq`
- standard RefSeq stable IDs (`NM_`, `NR_`, `XM_`, `XR_`)
- RefSeq mitochondrial exceptions when source mode is `RefSeq` or row source is
  normalized to `RefSeq`

This keeps edited transcript and CDS reconstruction paths available for
RefSeq-only caches.

## Options JSON Migration

`options_json` remains for existing non-source settings such as:
- `failed`
- `reference_fasta_path`
- `transcripts_table`
- `exons_table`
- `distance`
- `extended_probes`
- `partitioned`

But source mode is not read from `options_json`. If `merged` appears in
`options_json`, return:

```
annotate_vep(): options_json key 'merged' is unsupported; register/export cache tables with schema metadata bio.vep.cache_source_type='merged'
```

If `refseq` appears in `options_json`, return:

```
annotate_vep(): options_json key 'refseq' is unsupported; register/export cache tables with schema metadata bio.vep.cache_source_type='refseq'
```

If `cache_source_type` appears in `options_json`, ignore it as an authority and
continue requiring schema metadata. Tests should cover this to prevent a quiet
fallback.

## Cache Builder / Dependency Update

The optional `cache-builder` feature currently depends on
`datafusion-bio-format-ensembl-cache` and uses split translation schema helpers.
After updating that dependency to
`ca91f694e4db145ff0f8c1118b92a7f4c6f9aee3`, callers must pass an explicit
`CacheSourceType` into `EnsemblCacheOptions` and split translation schema
constructors. `CacheBuilder` should therefore require or store an explicit
source type instead of guessing from cache paths.

## Testing Strategy

Unit tests:
- parse `ensembl`, `merged`, `refseq`; reject missing/invalid metadata
- reject `{"merged": true}` and `{"merged": false}`
- source-specific transcript inclusion for `ENST`, `NM_`, `NR_`, `XM_`, `XR_`,
  `4540`, `COX3`, and `rna-TRNK`
- SOURCE output is populated only in merged mode
- hydration eligibility includes RefSeq mitochondrial IDs

Integration tests:
- in-memory MemTable schemas with metadata for each source mode
- gated real-cache tests using an environment variable such as
  `DF_BIO_FUNCTION_REFSEQ_CACHE=/Users/mwiewior/workspace/data_vepyr/homo_sapiens_refseq/115_GRCh38`
  after the companion format repo can export/register that cache with explicit
  source metadata
