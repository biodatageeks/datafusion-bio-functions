# Implementation Tasks

## 0. Research Inputs

- [x] 0.1 Inspect `datafusion-bio-formats` commit `ca91f694e4db145ff0f8c1118b92a7f4c6f9aee3`
- [x] 0.2 Inspect Ensembl VEP release/115 source-mode, RefSeq transcript-filter, MT exception, and CSQ field behavior
- [x] 0.3 Create detailed implementation handoff in `implementation-plan.md`

## 1. Cache Source Mode Contract

### 1.1 Add source-mode helper module
- [x] 1.1.1 Create `datafusion/bio-function-vep/src/cache_source.rs`
- [x] 1.1.2 Define `CacheSourceType::{Ensembl, Merged, RefSeq}`
- [x] 1.1.3 Add constant `CACHE_SOURCE_METADATA_KEY = "bio.vep.cache_source_type"`
- [x] 1.1.4 Implement strict parser accepting only `ensembl`, `merged`, and `refseq`
- [x] 1.1.5 Implement `CacheSourceType::from_schema(&Schema) -> Result<CacheSourceType>`
- [x] 1.1.6 Add unit tests for accepted values, missing metadata, invalid metadata, and case sensitivity

### 1.2 Export module from crate root
- [x] 1.2.1 Add `pub(crate) mod cache_source;` to `datafusion/bio-function-vep/src/lib.rs`
- [x] 1.2.2 Use the helper from `annotate_table_function.rs`, `annotate_provider.rs`, and `transcript_consequence.rs`
- [x] 1.2.3 Read source mode from partitioned parquet variation schema metadata at `annotate_vep()` planning time

## 2. Table Function Validation

### 2.1 Reject unsupported legacy source-selector options
- [x] 2.1.1 Add an `options_json_has_key(options, key)` helper near existing JSON option parsing
- [x] 2.1.2 In `AnnotateFunction::call()`, reject any `options_json` containing top-level `merged` or `refseq`
- [x] 2.1.3 Return error text: `annotate_vep(): options_json key 'merged' is unsupported; register/export cache tables with schema metadata bio.vep.cache_source_type='merged'`
- [x] 2.1.4 Return matching error text for `refseq` pointing to `bio.vep.cache_source_type='refseq'`
- [x] 2.1.5 Add tests for `{"merged":true}`, `{"merged":false}`, `{"refseq":true}`, and `{"refseq":false}`

### 2.2 Resolve cache source mode from schema metadata
- [x] 2.2.1 In `AnnotateFunction::call()`, resolve the main cache source table/provider schema
- [x] 2.2.2 Read `bio.vep.cache_source_type` from the main cache schema using `CacheSourceType::from_schema`
- [x] 2.2.3 Pass `CacheSourceType` into `AnnotateProvider::new(...)`
- [x] 2.2.4 Add tests proving missing metadata fails before provider execution
- [x] 2.2.5 Add tests proving `options_json.cache_source_type` does not replace schema metadata

## 3. AnnotateProvider Source Mode Wiring

### 3.1 Store source mode on provider
- [x] 3.1.1 Add `cache_source_type: CacheSourceType` field to `AnnotateProvider`
- [x] 3.1.2 Update `AnnotateProvider::new(...)` signature and call sites
- [x] 3.1.3 Update `Debug` output to include source mode

### 3.2 Remove boolean merged reads
- [x] 3.2.1 Remove `parse_json_bool_option(opts, "merged")` use from `scan_with_transcript_engine_partitioned`
- [x] 3.2.2 Remove `parse_json_bool_option(opts, "merged")` use from non-partitioned transcript context loading
- [x] 3.2.3 Remove `parse_json_bool_option(opts, "merged")` use from CSQ output construction
- [x] 3.2.4 Replace each removed use with `self.cache_source_type`

### 3.3 Validate context source metadata
- [x] 3.3.1 When variation/context tables are resolved, read schema metadata
- [x] 3.3.2 If context table metadata exists and differs from provider `cache_source_type`, return an error naming both values
- [x] 3.3.3 If context table metadata is missing, return a clear error requiring `bio.vep.cache_source_type`
- [x] 3.3.4 Add tests for mismatched main/transcript source metadata

## 4. Source-Specific Transcript Filtering

### 4.1 Replace `is_vep_transcript`
- [x] 4.1.1 Change `is_vep_transcript(id: &str, merged: bool)` to source-mode-aware filtering
- [x] 4.1.2 Prefer signature `is_vep_transcript(tx: &TranscriptFeature, source_type: CacheSourceType) -> bool`
- [x] 4.1.3 Implement Ensembl mode: accept only IDs starting `ENST`
- [x] 4.1.4 Implement RefSeq mode: accept standard RefSeq IDs matching `^[A-Z]{2}_\d+`
- [x] 4.1.5 Implement RefSeq mitochondrial exceptions: accept `^\d{4}$` and `^(rna-)?[A-Z0-9]{3,}$`
- [x] 4.1.6 Implement merged mode: accept `ENST` rows plus RefSeq rows matching RefSeq rules when row source normalizes to `RefSeq`

### 4.2 Update transcript filtering call sites
- [x] 4.2.1 Update partitioned context loading filter to call source-mode-aware `is_vep_transcript`
- [x] 4.2.2 Update non-partitioned context loading filter to call source-mode-aware `is_vep_transcript`
- [x] 4.2.3 Keep HGNC backfill behavior for merged caches; do not apply merged-only HGNC fill to RefSeq-only mode unless existing data requires it

### 4.3 Add transcript filter tests
- [x] 4.3.1 Test Ensembl mode includes `ENST00000367770` and excludes `NM_000546.6`
- [x] 4.3.2 Test RefSeq mode includes `NM_000546.6`, `NR_123456.1`, `XM_011520402.2`, and `XR_001734695.1`
- [x] 4.3.3 Test RefSeq mode includes mitochondrial names `4540`, `COX3`, and `rna-TRNK`
- [x] 4.3.4 Test RefSeq mode excludes `ENST00000367770`
- [x] 4.3.5 Test merged mode includes Ensembl rows and source-normalized RefSeq rows
- [x] 4.3.6 Test merged mode does not include gene-like non-RefSeq rows solely because their ID matches the mitochondrial exception regex

## 5. CSQ SOURCE Output

### 5.1 Gate SOURCE by source mode
- [x] 5.1.1 Replace `let source_val = if merged { source } else { "" };`
- [x] 5.1.2 Use `self.cache_source_type == CacheSourceType::Merged` as the only condition for emitting SOURCE
- [x] 5.1.3 Add unit or integration tests for SOURCE output in Ensembl, RefSeq, and merged modes

## 6. RefSeq Hydration

### 6.1 Extend hydration eligibility
- [x] 6.1.1 Update `is_refseq_transcript_for_hydration()` to accept `CacheSourceType`
- [x] 6.1.2 Keep current row-source normalization behavior for `RefSeq`
- [x] 6.1.3 Keep standard RefSeq ID detection for `NM_`, `NR_`, `XM_`, and `XR_`
- [x] 6.1.4 Add mitochondrial RefSeq ID detection when source mode is `RefSeq` or row source normalizes to `RefSeq`

### 6.2 Update hydration call sites
- [x] 6.2.1 Pass `self.cache_source_type` into `hydrate_refseq_translation_cds_from_reference`
- [x] 6.2.2 Pass source mode through any helper layers needed by `is_refseq_transcript_for_hydration`
- [x] 6.2.3 Add tests proving `4540`, `COX3`, and `rna-TRNK` are hydration-eligible in RefSeq mode

## 7. Integration and Regression Tests

### 7.1 Update existing MemTable fixtures
- [x] 7.1.1 Add `bio.vep.cache_source_type = "ensembl"` metadata to existing annotation test schemas that represent Ensembl-only caches
- [x] 7.1.2 Add helper for constructing test schemas with cache source metadata
- [x] 7.1.3 Verify existing annotation tests pass after the explicit metadata requirement

### 7.2 Add source-mode integration tests
- [ ] 7.2.1 Create an in-memory RefSeq transcript table with IDs `NM_000546.6`, `NR_123456.1`, `4540`, and `rna-TRNK`
- [x] 7.2.2 Verify `annotate_vep()` includes RefSeq transcripts when metadata is `refseq`
- [x] 7.2.3 Verify `annotate_vep()` excludes RefSeq transcripts when metadata is `ensembl`
- [x] 7.2.4 Verify merged mode emits SOURCE and RefSeq mode leaves SOURCE empty

### 7.3 Add gated real-cache validation
- [ ] 7.3.1 Add ignored/gated test controlled by `DF_BIO_FUNCTION_REFSEQ_CACHE`
- [ ] 7.3.2 Use the local RefSeq cache exported/registered with explicit source metadata
- [ ] 7.3.3 Compare representative chr22 and MT variant outputs against Ensembl VEP

## 8. Verification

### 8.1 Local checks
- [x] 8.1.1 Run `cargo fmt`
- [x] 8.1.2 Run `cargo test -p datafusion-bio-function-vep cache_source`
- [x] 8.1.3 Run `cargo test -p datafusion-bio-function-vep is_vep_transcript`
- [x] 8.1.4 Run `cargo test -p datafusion-bio-function-vep annotate_vep`
- [x] 8.1.5 Run `cargo clippy -p datafusion-bio-function-vep --all-targets -- -D warnings`

### 8.2 Cross-repo validation
- [x] 8.2.1 Update optional `datafusion-bio-format-ensembl-cache` dependency to `ca91f694e4db145ff0f8c1118b92a7f4c6f9aee3`
- [ ] 8.2.2 Run a RefSeq chr22 smoke annotation with `cache_source_type = refseq`
- [ ] 8.2.3 Run a merged-mode smoke annotation with `cache_source_type = merged`
- [x] 8.2.4 Confirm no code path reads `options_json.merged`
- [x] 8.2.5 Confirm no code path reads `options_json.refseq` as source mode
