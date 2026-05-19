# Add VEP Cache Source Mode Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Make `annotate_vep()` consume explicit VEP cache source mode metadata (`bio.vep.cache_source_type = ensembl | merged | refseq`), remove source selection from `options_json`, and match Ensembl VEP behavior for Ensembl, merged, and RefSeq caches.

**Architecture:** Add a small internal cache-source module, resolve source mode from the variation parquet schema during table-function planning, carry `CacheSourceType` through `AnnotateProvider`, validate all cache table schemas against the same mode, and make transcript filtering, RefSeq hydration, and CSQ field selection source-mode aware.

**Tech Stack:** Rust, Apache Arrow/DataFusion, Parquet, `datafusion-bio-function-vep`, optional `datafusion-bio-format-ensembl-cache` dependency pinned to BioFormats commit `ca91f694e4db145ff0f8c1118b92a7f4c6f9aee3`.

---

## Evidence Checked

- `datafusion-bio-formats` commit `ca91f694e4db145ff0f8c1118b92a7f4c6f9aee3` adds `CacheSourceType`, exports `VEP_CACHE_SOURCE_TYPE_METADATA_KEY`, requires `EnsemblCacheOptions.cache_source_type`, and writes `bio.vep.cache_source_type` on all VEP cache schemas.
- Ensembl VEP `release/115` derives cache source mode from `--merged` and `--refseq`, rejects source mode embedded in species names, applies RefSeq filtering only for RefSeq caches or merged rows whose `_source_cache` is `RefSeq`, and exposes CSQ `SOURCE` only through the merged field set.
- Current `bio-function-vep` builds the output schema in `AnnotateProvider::new(...)`, so source mode must be known before provider construction. For the current partitioned-cache path API, this means reading the Arrow schema metadata from a `variation/*.parquet` file during `AnnotateFunction::call()`.
- Current `options_json` handling still reads source selectors (`merged`, `refseq`) in several places. Those reads must be removed or converted into explicit errors.

## Files To Change

- `datafusion/bio-function-vep/Cargo.toml`
- `datafusion/bio-function-vep/src/lib.rs`
- `datafusion/bio-function-vep/src/cache_source.rs` (new)
- `datafusion/bio-function-vep/src/annotate_table_function.rs`
- `datafusion/bio-function-vep/src/annotate_provider.rs`
- `datafusion/bio-function-vep/src/transcript_consequence.rs`
- `datafusion/bio-function-vep/src/cache_builder.rs`
- `datafusion/bio-function-vep/src/vcf_sink.rs`
- `datafusion/bio-function-vep/examples/annotate_vep_golden_bench.rs`
- Existing tests in `datafusion/bio-function-vep/src/annotate_table_function.rs` and any source-specific unit tests close to the changed helpers.

## Implementation Steps

### 1. Add the Cache Source Contract

- [ ] 1.1 Create `datafusion/bio-function-vep/src/cache_source.rs`.
- [ ] 1.2 Add a local enum instead of depending on the optional BioFormats crate from normal function code:

```rust
use std::fs::File;
use std::path::PathBuf;
use std::str::FromStr;

use arrow::datatypes::Schema;
use datafusion::common::{DataFusionError, Result};
use parquet::arrow::arrow_reader::{ArrowReaderMetadata, ArrowReaderOptions};

pub(crate) const CACHE_SOURCE_METADATA_KEY: &str = "bio.vep.cache_source_type";

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub(crate) enum CacheSourceType {
    Ensembl,
    Merged,
    RefSeq,
}

impl CacheSourceType {
    pub(crate) fn as_str(self) -> &'static str {
        match self {
            Self::Ensembl => "ensembl",
            Self::Merged => "merged",
            Self::RefSeq => "refseq",
        }
    }

    pub(crate) fn from_schema(schema: &Schema) -> Result<Self> {
        let raw = schema.metadata().get(CACHE_SOURCE_METADATA_KEY).ok_or_else(|| {
            DataFusionError::Plan(format!(
                "annotate_vep(): cache table schema is missing {CACHE_SOURCE_METADATA_KEY}; expected one of ensembl, merged, refseq"
            ))
        })?;
        raw.parse()
    }

    pub(crate) fn from_partitioned_cache_source(cache_source: &str) -> Result<Self> {
        let parquet = first_variation_parquet(cache_source)?;
        let file = File::open(&parquet).map_err(|err| {
            DataFusionError::Execution(format!(
                "annotate_vep(): failed to open cache variation parquet '{}': {err}",
                parquet.display()
            ))
        })?;
        let metadata = ArrowReaderMetadata::load(&file, ArrowReaderOptions::default()).map_err(|err| {
            DataFusionError::Execution(format!(
                "annotate_vep(): failed to read Arrow schema metadata from cache variation parquet '{}': {err}",
                parquet.display()
            ))
        })?;
        Self::from_schema(metadata.schema().as_ref())
    }
}

impl FromStr for CacheSourceType {
    type Err = DataFusionError;

    fn from_str(raw: &str) -> Result<Self> {
        match raw {
            "ensembl" => Ok(Self::Ensembl),
            "merged" => Ok(Self::Merged),
            "refseq" => Ok(Self::RefSeq),
            other => Err(DataFusionError::Plan(format!(
                "annotate_vep(): unsupported {CACHE_SOURCE_METADATA_KEY} '{other}'; expected one of ensembl, merged, refseq"
            ))),
        }
    }
}
```

- [ ] 1.3 Implement `first_variation_parquet(cache_source: &str) -> Result<PathBuf>` using the existing partitioned layout convention. It should search under `<cache_source>/variation` recursively enough for the existing layout and return the first `.parquet` file in deterministic sorted order.
- [ ] 1.4 Add `pub(crate) mod cache_source;` to `datafusion/bio-function-vep/src/lib.rs`.
- [ ] 1.5 Add unit tests for accepted values, missing metadata, invalid values, and case sensitivity.

Implementation note: keep this enum local because `datafusion-bio-format-ensembl-cache` is optional behind `cache-builder`; normal `annotate_vep()` builds should not gain a required dependency on that crate.

### 2. Reject Legacy Source Selectors at Table-Function Planning

- [ ] 2.1 In `annotate_table_function.rs`, add a helper that inspects only top-level JSON object keys:

```rust
fn options_json_has_key(options_json: Option<&str>, key: &str) -> Result<bool> {
    let Some(raw) = options_json else {
        return Ok(false);
    };
    let Some(value) = parse_options_json_value(raw)? else {
        return Ok(false);
    };
    Ok(value.as_object().is_some_and(|obj| obj.contains_key(key)))
}

fn reject_options_json_source_selectors(options_json: Option<&str>) -> Result<()> {
    for key in ["merged", "refseq"] {
        if options_json_has_key(options_json, key)? {
            return Err(DataFusionError::Plan(format!(
                "annotate_vep(): options_json key '{key}' is unsupported; register/export cache tables with schema metadata {CACHE_SOURCE_METADATA_KEY}='{key}'"
            )));
        }
    }
    Ok(())
}
```

- [ ] 2.2 If the existing JSON parsing helper is private to another module, either move the selector rejection helper there or parse with `serde_json::from_str::<serde_json::Value>()` locally and preserve existing invalid JSON error style.
- [ ] 2.3 In `AnnotateFunction::call()`, call `reject_options_json_source_selectors(options_json.as_deref())?` immediately after extracting `options_json`.
- [ ] 2.4 Still ignore `options_json.cache_source_type` as an authority. Missing schema metadata must remain an error even if that option exists.
- [ ] 2.5 Add tests for `{"merged":true}`, `{"merged":false}`, `{"refseq":true}`, and `{"refseq":false}`.

The exact error strings expected by tests:

```text
annotate_vep(): options_json key 'merged' is unsupported; register/export cache tables with schema metadata bio.vep.cache_source_type='merged'
annotate_vep(): options_json key 'refseq' is unsupported; register/export cache tables with schema metadata bio.vep.cache_source_type='refseq'
```

### 3. Resolve Source Mode Before Provider Construction

- [ ] 3.1 In `AnnotateFunction::call()`, after validating `options_json` and before `AnnotateProvider::new(...)`, read the mode:

```rust
let cache_source_type = CacheSourceType::from_partitioned_cache_source(&cache_source)?;
```

- [ ] 3.2 Pass the mode into provider construction:

```rust
let provider = AnnotateProvider::new(
    session,
    vcf_table,
    cache_source,
    backend,
    cache_source_type,
    options_json,
    vcf_schema,
)?;
```

- [ ] 3.3 Update every `AnnotateProvider::new(...)` call site, including temporary providers created inside context setup.
- [ ] 3.4 Add missing-metadata and invalid-metadata tests that fail before execution because `schema()` must already know the source-mode-specific CSQ fields.
- [ ] 3.5 Add a test where `options_json` contains `{"cache_source_type":"refseq"}` but parquet metadata is missing; it must still fail for missing metadata.

### 4. Store Source Mode on `AnnotateProvider`

- [ ] 4.1 In `annotate_provider.rs`, import `crate::cache_source::CacheSourceType`.
- [ ] 4.2 Add `cache_source_type: CacheSourceType` to `AnnotateProvider`.
- [ ] 4.3 Change `TranscriptSelectionFlags` to carry source mode:

```rust
struct TranscriptSelectionFlags {
    cache_source_type: CacheSourceType,
    gencode_basic: bool,
    gencode_primary: bool,
    all_refseq: bool,
    exclude_predicted: bool,
}
```

- [ ] 4.4 Change construction from `TranscriptSelectionFlags::from_options_json(options_json)` to:

```rust
let transcript_selection =
    TranscriptSelectionFlags::from_options_json(cache_source_type, options_json)?;
```

- [ ] 4.5 Remove `TranscriptSourceMode`. Use `CacheSourceType` everywhere instead.
- [ ] 4.6 Update `Debug` output to include `cache_source_type`.
- [ ] 4.7 Add `cache_source_type: CacheSourceType` to `ContigAnnotationConfig` and pass it through `hydrate_window(...)`, lookup construction, and any temporary provider setup.

### 5. Remove Source-Mode Reads from `options_json`

- [ ] 5.1 Remove all source-mode reads like:

```rust
parse_json_bool_option(opts, "merged")
parse_json_bool_option(opts, "refseq")
```

- [ ] 5.2 Keep non-source options unchanged: `failed`, `reference_fasta_path`, `transcripts_table`, `exons_table`, `distance`, `extended_probes`, `partitioned`, `all_refseq`, `exclude_predicted`, `gencode_basic`, `gencode_primary`, and backend options.
- [ ] 5.3 Update `TranscriptSelectionFlags::from_options_json(...)` so `all_refseq` and `exclude_predicted` are only meaningful for RefSeq-capable modes:

```rust
match cache_source_type {
    CacheSourceType::Ensembl if all_refseq || exclude_predicted => {
        return Err(DataFusionError::Plan(
            "annotate_vep(): all_refseq/exclude_predicted require cache schema metadata bio.vep.cache_source_type='refseq' or 'merged'".to_string(),
        ));
    }
    _ => {}
}
```

- [ ] 5.4 Run the source search at the end of the implementation and confirm no source mode is read from `options_json`:

```bash
rg 'parse_json_bool_option\(opts, "(merged|refseq)"|"\\"merged\\""|"\\"refseq\\""' \
  datafusion/bio-function-vep/src \
  datafusion/bio-function-vep/examples
```

Any remaining hits must be tests asserting rejection, fixture JSON, or user-facing error strings.

### 6. Validate All Cache Table Metadata

- [ ] 6.1 Add an async helper on `AnnotateProvider`:

```rust
async fn validate_registered_cache_source(&self, table: &str, role: &str) -> Result<()> {
    let df = self.session.table(table).await?;
    let schema = df.schema();
    let actual = CacheSourceType::from_schema(schema.as_arrow().as_ref())?;
    if actual != self.cache_source_type {
        return Err(DataFusionError::Plan(format!(
            "annotate_vep(): {role} table '{table}' has {CACHE_SOURCE_METADATA_KEY}='{}' but variation cache has {CACHE_SOURCE_METADATA_KEY}='{}'",
            actual.as_str(),
            self.cache_source_type.as_str(),
        )));
    }
    Ok(())
}
```

- [ ] 6.2 If `DFSchemaRef::as_arrow()` returns a borrowed schema rather than `SchemaRef` in the local DataFusion version, adjust the `.as_ref()` part to match the current API.
- [ ] 6.3 Call this helper after registering or resolving each cache table used for a scan:
  - variation table for the current contig
  - transcript table
  - exon table
  - translation table
  - split translation tables
  - regulatory table
  - motif table
- [ ] 6.4 Missing metadata should be an error for these VEP cache tables. Do not silently skip missing metadata for backward compatibility.
- [ ] 6.5 Add mismatch tests where variation metadata is `refseq` and transcript metadata is `ensembl`, and where transcript metadata is missing.

### 7. Make Transcript Filtering Match Ensembl VEP

- [ ] 7.1 Replace `TranscriptSelectionFlags::source_mode` logic with `cache_source_type`.
- [ ] 7.2 Keep the existing prefilters for failed, gencode, `all_refseq`, and `exclude_predicted`, but base source-mode inclusion on `CacheSourceType`.
- [ ] 7.3 Add source helpers near the existing transcript filtering helpers:

```rust
fn is_standard_refseq_accession(id: &str) -> bool {
    let bytes = id.as_bytes();
    bytes.len() >= 4 && bytes[0].is_ascii_uppercase() && bytes[1].is_ascii_uppercase() && bytes[2] == b'_' && bytes[3].is_ascii_digit()
}

fn is_mitochondrial_refseq_exception(tx: &TranscriptFeature) -> bool {
    is_mitochondrial_contig(&tx.seq_region_name)
        && (is_mitochondrial_refseq_name(&tx.transcript_id)
            || tx
                .display_xref
                .as_deref()
                .is_some_and(is_mitochondrial_refseq_name))
}

fn is_default_refseq_transcript_id(tx: &TranscriptFeature) -> bool {
    is_standard_refseq_accession(&tx.transcript_id)
        || is_mitochondrial_refseq_exception(tx)
        || tx
            .display_xref
            .as_deref()
            .is_some_and(is_standard_refseq_accession)
}
```

- [ ] 7.4 Implement merged-row source detection using the normalized row source already present in the provider. Prefer `_source_cache`/`source_cache` when available; fall back to the exported source label only if that is the existing row field:

```rust
fn row_source_is_refseq(tx: &TranscriptFeature) -> bool {
    normalize_source_label(tx.source_cache.as_deref().or(tx.source.as_deref()))
        .is_some_and(|source| source == "RefSeq")
}
```

- [ ] 7.5 Update `passes_transcript_selection(...)` to this shape:

```rust
match selection.cache_source_type {
    CacheSourceType::Ensembl => true,
    CacheSourceType::RefSeq => {
        selection.all_refseq || is_default_refseq_transcript_id(tx)
    }
    CacheSourceType::Merged => {
        if row_source_is_refseq(tx) {
            selection.all_refseq || is_default_refseq_transcript_id(tx)
        } else {
            true
        }
    }
}
```

- [ ] 7.6 Keep `exclude_predicted` as an earlier filter for `XM_`/`XR_` in RefSeq-capable modes.
- [ ] 7.7 Retire `transcript_consequence::is_vep_transcript(id, merged)` if it remains unused. If tests or future call sites need it, replace it with a `TranscriptFeature` plus `CacheSourceType` signature.
- [ ] 7.8 Add unit tests for Ensembl, RefSeq, merged, non-`ENST` Ensembl-side IDs, numeric MT IDs, `COX3`, `rna-TRNK`, and a non-RefSeq merged row whose ID would otherwise match the MT exception pattern.

### 8. Drive CSQ Fields from Source Mode

- [ ] 8.1 Update annotation-column selection so RefSeq-specific fields appear for `CacheSourceType::RefSeq` and `CacheSourceType::Merged`.
- [ ] 8.2 Ensure the CSQ `SOURCE` field appears only for `CacheSourceType::Merged`.
- [ ] 8.3 Replace boolean gates like:

```rust
let source_val = if merged { source } else { "" };
```

with:

```rust
let source_val = if self.cache_source_type == CacheSourceType::Merged {
    source
} else {
    ""
};
```

- [ ] 8.4 Update `CsqPlaceholderLayout::for_mode(...)` and benchmark field comparison setup to derive `refseq_fields` and `source_field` from `CacheSourceType`:

```rust
let refseq_fields = matches!(cache_source_type, CacheSourceType::RefSeq | CacheSourceType::Merged);
let source_field = cache_source_type == CacheSourceType::Merged;
```

- [ ] 8.5 Add tests proving merged mode emits `SOURCE`, while RefSeq-only and Ensembl-only modes leave it absent or empty according to the established schema contract.

### 9. Extend RefSeq Hydration Eligibility

- [ ] 9.1 Add `cache_source_type: CacheSourceType` to:
  - `hydrate_window(...)`
  - `hydrate_refseq_translation_cds_from_reference(...)`
  - `hydrate_transcript_cdna_from_reference(...)`
  - `is_refseq_transcript_for_hydration(...)`
- [ ] 9.2 Implement hydration eligibility as:

```rust
fn is_refseq_transcript_for_hydration(tx: &TranscriptFeature, cache_source_type: CacheSourceType) -> bool {
    row_source_is_refseq(tx)
        || is_standard_refseq_accession(&tx.transcript_id)
        || matches!(cache_source_type, CacheSourceType::RefSeq)
            && is_mitochondrial_refseq_exception(tx)
}
```

- [ ] 9.3 For merged mode, allow MT exceptions when `row_source_is_refseq(tx)` is true.
- [ ] 9.4 Add tests proving `4540`, `COX3`, and `rna-TRNK` are hydration-eligible in RefSeq mode and for merged RefSeq rows.

### 10. Update Test Fixtures to Carry Metadata

- [ ] 10.1 Add helpers in `annotate_table_function.rs` tests:

```rust
fn schema_with_cache_source(fields: Vec<Field>, source: CacheSourceType) -> Arc<Schema> {
    let mut metadata = HashMap::new();
    metadata.insert(CACHE_SOURCE_METADATA_KEY.to_string(), source.as_str().to_string());
    Arc::new(Schema::new(fields).with_metadata(metadata))
}

fn batch_with_cache_source(batch: &RecordBatch, source: CacheSourceType) -> RecordBatch {
    let schema = batch.schema();
    let mut metadata = schema.metadata().clone();
    metadata.insert(CACHE_SOURCE_METADATA_KEY.to_string(), source.as_str().to_string());
    let schema = Arc::new(schema.as_ref().clone().with_metadata(metadata));
    RecordBatch::try_new(schema, batch.columns().to_vec()).expect("metadata-only schema update")
}
```

- [ ] 10.2 Add `write_batch_to_cache_with_source(...)` and `write_batch_to_chrom_with_source(...)`.
- [ ] 10.3 Keep existing fixture helpers as Ensembl defaults by forwarding to the source-aware helpers with `CacheSourceType::Ensembl`.
- [ ] 10.4 Add a dedicated helper that writes parquet without source metadata for missing-metadata tests.
- [ ] 10.5 Update RefSeq and merged tests to use source metadata instead of `options_json` source booleans.

### 11. Update Optional Cache Builder Integration

- [ ] 11.1 In `datafusion/bio-function-vep/Cargo.toml`, update the optional dependency:

```toml
datafusion-bio-format-ensembl-cache = { git = "https://github.com/DataFusion-Bio/datafusion-bio-formats.git", rev = "ca91f694e4db145ff0f8c1118b92a7f4c6f9aee3", optional = true }
```

- [ ] 11.2 In `cache_builder.rs`, import the BioFormats enum under the `cache-builder` feature:

```rust
use datafusion_bio_format_ensembl_cache::CacheSourceType as BioFormatsCacheSourceType;
```

- [ ] 11.3 Add a source-type field to `CacheBuilder` and require it in construction, or add a builder setter that must be called before execution:

```rust
pub struct CacheBuilder {
    cache_root: PathBuf,
    output_dir: PathBuf,
    cache_source_type: BioFormatsCacheSourceType,
}
```

- [ ] 11.4 Pass it to BioFormats options:

```rust
let mut options = EnsemblCacheOptions::new(cache_root)
    .with_cache_source_type(self.cache_source_type);
options.target_partitions = Some(partitions);
```

- [ ] 11.5 Update split translation schema calls:

```rust
translation_core_schema(false, self.cache_source_type)
translation_sift_schema(false, self.cache_source_type)
```

- [ ] 11.6 Update cache-builder tests and examples to pass `BioFormatsCacheSourceType::Ensembl` unless they intentionally build RefSeq or merged caches.

### 12. Update VCF Sink and Golden Benchmark Options

- [ ] 12.1 In `vcf_sink.rs`, stop serializing `refseq` and `merged` into `options_json`.
- [ ] 12.2 If `AnnotateVcfConfig` still exposes `refseq` and `merged`, either:
  - return a clear error that source mode now comes from cache schema metadata, or
  - deprecate the fields and ignore them only after tests verify they are not serialized.
- [ ] 12.3 Prefer an explicit error for now, because silently ignoring these fields can hide a user migration bug.
- [ ] 12.4 In `examples/annotate_vep_golden_bench.rs`, remove only our `annotate_vep()` `options_json` source flags. Keep Docker VEP `--refseq` and `--merged` flags for expected-output generation.

### 13. Verification

- [ ] 13.1 Run formatting:

```bash
cargo fmt
```

- [ ] 13.2 Run focused tests:

```bash
cargo test -p datafusion-bio-function-vep cache_source
cargo test -p datafusion-bio-function-vep is_vep_transcript
cargo test -p datafusion-bio-function-vep test_annotate_vep_rejects_options_json_source_selectors
cargo test -p datafusion-bio-function-vep test_annotate_vep_refseq_and_merged_modes_emit_refseq_specific_csq_fields
```

- [ ] 13.3 Run the broader package tests:

```bash
cargo test -p datafusion-bio-function-vep annotate_vep
```

- [ ] 13.4 Run cache-builder checks after updating the BioFormats dependency:

```bash
cargo test -p datafusion-bio-function-vep --features cache-builder cache_builder
cargo clippy -p datafusion-bio-function-vep --all-targets --features cache-builder -- -D warnings
```

- [ ] 13.5 Confirm source mode is no longer read from `options_json`:

```bash
rg 'parse_json_bool_option\(opts, "(merged|refseq)"|"\\"merged\\""|"\\"refseq\\""' \
  datafusion/bio-function-vep/src \
  datafusion/bio-function-vep/examples
```

- [ ] 13.6 Validate the OpenSpec change:

```bash
openspec validate add-vep-cache-source-mode --strict
```

## Commit Boundary

Recommended implementation commits:

1. Add `cache_source.rs`, planning-time metadata read, provider source-mode field, and tests for missing/invalid metadata.
2. Remove `options_json` source selectors and update VCF sink/golden benchmark option generation.
3. Replace transcript filtering, CSQ field selection, and RefSeq hydration with source-mode-aware behavior.
4. Update fixtures and integration tests for Ensembl, RefSeq, and merged modes.
5. Update optional cache-builder dependency and API calls to BioFormats commit `ca91f694e4db145ff0f8c1118b92a7f4c6f9aee3`.
