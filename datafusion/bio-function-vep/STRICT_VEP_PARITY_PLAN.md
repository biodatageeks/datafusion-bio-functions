# Strict Ensembl VEP Parity Plan

Goal: get `datafusion/bio-function-vep` to zero mismatches on the README chr1 golden benchmark for all 74 CSQ columns, with zero extra CSQ entries, by implementing Ensembl VEP semantics exactly rather than tuning heuristics.

Date: March 12, 2026

## Hard Constraints

- No workaround logic that only fits the current chr1 benchmark.
- No heuristic fallback when Ensembl VEP has a defined algorithm or data flow.
- A review of the existing Rust implementation is required before rewriting any parity-sensitive area so the current behavior, deltas, and hidden coupling are documented first.
- Any behavior change must be traceable to a specific Ensembl VEP source path or to the underlying Ensembl Variation code used by VEP.
- If the required upstream behavior lives outside the current `ensembl-vep` checkout, the missing upstream module must be brought in or reimplemented from source-equivalent logic before the parity work is considered complete.
- Every parity-sensitive method that is added or modified must include a method-level comment pointing to the exact Ensembl VEP GitHub permalink that defines the implemented behavior.
- The traceability comment must include a specific file and line range, not just a repo root, directory, or issue link.
- If one Rust method implements logic derived from multiple upstream methods, the comment must list each upstream permalink and state which part of the behavior each link covers.
- Comments must use stable GitHub blob permalinks pinned to the Ensembl VEP commit/release being targeted so the reference does not drift.
- After Phase 1, the plan document must be updated with implementation findings:
  - remaining discrepancies
  - any workarounds still present
  - any newly discovered divergences from Ensembl VEP or its upstream dependencies
  - any source-traceability gaps that blocked exact parity

## Scope

Primary target:

- README chr1 golden benchmark command produces:
  - `0` mismatches for all 74 CSQ fields
  - `0` extra CSQ entries
  - `0` missing CSQ entries
  - `100%` `most_severe_consequence`
  - `100%` term-set parity

Secondary target:

- Close the hidden parity gap for HGVS behavior that the current 74-field benchmark does not exercise because the harness does not pass `--hgvs --fasta`.

## Current Baseline

Validated on this checkout with:

```bash
cargo run -p datafusion-bio-function-vep --example annotate_vep_golden_bench --release -- \
  vep-benchmark/data/HG002_chr1.vcf.gz \
  /Users/mwiewior/research/data/vep/115_GRCh38_variation_1.parquet \
  parquet \
  0 \
  /tmp/vep_golden \
  vep-benchmark/data/HG002_chr1.vcf.gz \
  /tmp/annotate_vep_golden_bench_curr \
  /Users/mwiewior/research/data/vep \
  --steps=ensembl,datafusion \
  --extended-probes
```

Observed baseline:

- `323,430` variants compared
- `2,997,504` CSQ entries compared
- `53/74` fields perfect
- `21/74` fields imperfect
- `15` extra CSQ entries on our side

Imperfect fields:

- `Consequence`
- `INTRON`
- `CDS_position`
- `Protein_position`
- `Amino_acids`
- `Existing_variation`
- `FLAGS`
- `AF`
- `AFR_AF`
- `AMR_AF`
- `EAS_AF`
- `EUR_AF`
- `SAS_AF`
- `gnomADe_AF`
- `gnomADg_AF`
- `MAX_AF`
- `MAX_AF_POPS`
- `CLIN_SIG`
- `SOMATIC`
- `PHENO`
- `PUBMED`

## Upstream Source of Truth

These files define the behavior that the Rust implementation must match.

### Existing/co-located variant semantics

- `modules/Bio/EnsEMBL/VEP/AnnotationType/Variation.pm`
  - `compare_existing`
  - shifted and unshifted allele matching
- `modules/Bio/EnsEMBL/VEP/InputBuffer.pm`
  - overlap against both shifted and unshifted coordinates
- `modules/Bio/EnsEMBL/VEP/OutputFactory.pm`
  - `add_colocated_variant_info`
  - `add_colocated_frequency_data`
- `modules/Bio/EnsEMBL/VEP/AnnotationSource/Cache/BaseCacheVariation.pm`
  - allele-to-frequency mapping used before output
- external dependency:
  - `Bio::EnsEMBL::Variation::Utils::Sequence::get_matched_variant_alleles`

### Transcript consequence and field serialization semantics

- `modules/Bio/EnsEMBL/VEP/OutputFactory.pm`
  - transcript flags
  - `Amino_acids`
  - `CDS_position`
  - `Protein_position`
  - HGVS output
- underlying Ensembl Variation transcript-variation logic used by VEP
  - splice predicates
  - exon/intron numbering
  - coding consequence classification
  - peptide allele string generation
  - boundary shift handling

## Main Diagnosis

### 1. Co-located variant handling is not modeled like VEP

Current Rust code:

- builds a position bucket of nearby cache rows
- derives `Existing_variation`, AF, and clinical metadata from that bucket
- uses local fallback rules to choose AF when the main joined row is empty

Relevant files:

- [annotate_provider.rs](/Users/mwiewior/research/git/datafusion-bio-functions/datafusion/bio-function-vep/src/annotate_provider.rs#L386)
- [annotate_provider.rs](/Users/mwiewior/research/git/datafusion-bio-functions/datafusion/bio-function-vep/src/annotate_provider.rs#L1790)
- [variant_lookup_exec.rs](/Users/mwiewior/research/git/datafusion-bio-functions/datafusion/bio-function-vep/src/variant_lookup_exec.rs#L87)
- [variant_lookup_exec.rs](/Users/mwiewior/research/git/datafusion-bio-functions/datafusion/bio-function-vep/src/variant_lookup_exec.rs#L772)

VEP behavior:

- computes `matched_alleles` per existing variant
- checks both shifted and unshifted input alleles
- filters IDs and metadata by matched allele
- computes AF and `MAX_AF` from matched allele frequency data, not by scanning nearby rows

Relevant upstream files:

- [Variation.pm](/Users/mwiewior/research/git/ensembl-vep/modules/Bio/EnsEMBL/VEP/AnnotationType/Variation.pm#L146)
- [InputBuffer.pm](/Users/mwiewior/research/git/ensembl-vep/modules/Bio/EnsEMBL/VEP/InputBuffer.pm#L312)
- [OutputFactory.pm](/Users/mwiewior/research/git/ensembl-vep/modules/Bio/EnsEMBL/VEP/OutputFactory.pm#L1005)
- [OutputFactory.pm](/Users/mwiewior/research/git/ensembl-vep/modules/Bio/EnsEMBL/VEP/OutputFactory.pm#L1139)
- [BaseCacheVariation.pm](/Users/mwiewior/research/git/ensembl-vep/modules/Bio/EnsEMBL/VEP/AnnotationSource/Cache/BaseCacheVariation.pm#L185)
- [BaseCacheVariation.pm](/Users/mwiewior/research/git/ensembl-vep/modules/Bio/EnsEMBL/VEP/AnnotationSource/Cache/BaseCacheVariation.pm#L398)

This is the root cause for the largest mismatch block:

- `Existing_variation`
- all mismatching AF fields
- `MAX_AF`
- `MAX_AF_POPS`
- most `CLIN_SIG`
- most `SOMATIC`
- most `PHENO`
- most `PUBMED`

### 2. Transcript consequence engine still contains approximation logic

Current Rust code still uses:

- heuristic start/stop classification
- post-hoc parent-term stripping as a safety net
- custom intron/splice overlap approximations
- custom peptide formatting rules when exact codon classification misses

Relevant files:

- [transcript_consequence.rs](/Users/mwiewior/research/git/datafusion-bio-functions/datafusion/bio-function-vep/src/transcript_consequence.rs#L600)
- [transcript_consequence.rs](/Users/mwiewior/research/git/datafusion-bio-functions/datafusion/bio-function-vep/src/transcript_consequence.rs#L1048)
- [transcript_consequence.rs](/Users/mwiewior/research/git/datafusion-bio-functions/datafusion/bio-function-vep/src/transcript_consequence.rs#L1943)

This is the root cause for:

- `Consequence`
- `INTRON`
- `CDS_position`
- `Protein_position`
- `Amino_acids`
- the `15` extra CSQ entries

### 3. FLAGS ordering is still synthesized in some cases

Current Rust code reconstructs `FLAGS` when ordered source data is unavailable.

Relevant file:

- [transcript_consequence.rs](/Users/mwiewior/research/git/datafusion-bio-functions/datafusion/bio-function-vep/src/transcript_consequence.rs#L3012)

VEP emits attribute order directly from transcript attributes.

Relevant upstream file:

- [OutputFactory.pm](/Users/mwiewior/research/git/ensembl-vep/modules/Bio/EnsEMBL/VEP/OutputFactory.pm#L1425)

### 4. The current benchmark does not prove HGVS parity

The benchmark requests `HGVSc` and `HGVSp` in the CSQ field list but does not enable `--hgvs`, and our serializer currently emits empty strings for those fields in transcript CSQ assembly.

Relevant files:

- [annotate_vep_golden_bench.rs](/Users/mwiewior/research/git/datafusion-bio-functions/datafusion/bio-function-vep/examples/annotate_vep_golden_bench.rs#L514)
- [annotate_provider.rs](/Users/mwiewior/research/git/datafusion-bio-functions/datafusion/bio-function-vep/src/annotate_provider.rs#L2028)
- [OutputFactory.pm](/Users/mwiewior/research/git/ensembl-vep/modules/Bio/EnsEMBL/VEP/OutputFactory.pm#L1696)

This does not block the README 74-field target, but it does block true end-to-end parity with VEP.

## Execution Plan

## Phase -1: Existing Code Review and Traceability Inventory

Purpose: review the current Rust implementation before changing semantics, and establish a per-method mapping to Ensembl VEP source.

Tasks:

1. Review each existing parity-sensitive Rust module before modification:
   - `variant_lookup_exec.rs`
   - `annotate_provider.rs`
   - `transcript_consequence.rs`
   - `hgvs.rs`
   - any helper modules involved in allele normalization or CSQ serialization
2. For each parity-sensitive function or method, record:
   - current Rust behavior
   - known mismatch classes it contributes to
   - intended Ensembl VEP source method(s)
   - whether the current implementation is exact, partial, heuristic, or unrelated
3. Add or update method comments so every parity-sensitive method has a traceability block with:
   - Ensembl VEP GitHub permalink
   - exact file and line range
   - short note on what semantic contract is being implemented
4. If no exact upstream source has been identified for a method yet:
   - mark the method as unresolved
   - do not rewrite it as “VEP parity complete”
5. Create a review checklist so no parity-sensitive method is merged without:
   - implementation review
   - upstream source traceability comment
   - validation evidence

Acceptance:

- every parity-sensitive method has an explicit upstream traceability comment
- every planned rewrite starts from a reviewed Rust baseline rather than direct replacement

## Phase 0: Lock the Validation Surface

Purpose: prevent false positives while the implementation changes.

Tasks:

1. Treat the chr1 benchmark command above as the release gate.
2. Add a machine-readable comparison summary for:
   - per-field mismatch counts
   - extra CSQ entries
   - missing CSQ entries
3. Add a second gate that fails if:
   - `ours_total_entries != golden_total_entries`
   - `ours_only != 0`
   - `golden_only != 0`
4. Preserve a named list of the currently known mismatch loci as fixed regression fixtures.
5. Add a second benchmark mode for future HGVS validation:
   - VEP with `--hgvs --fasta`
   - Rust with real HGVS serialization enabled

Acceptance:

- benchmark failure becomes deterministic and directly attributable to field-level or entry-level regressions

## Phase 1: Port VEP Existing-Variant Matching Exactly

Purpose: eliminate all co-located ID/frequency/clinical mismatches by replacing the current bucket/fallback model with VEP’s matched-allele model.

### 1.1 Introduce explicit matched-allele state

Replace the current co-located aggregation model with a structure equivalent to VEP’s per-existing-variant object:

- existing variant identity
- cache coordinates
- allele string
- matched allele pairs
- matched through shifted input
- matched through unshifted input
- per-row frequency data
- per-row clinical metadata

Required code changes:

- replace the coarse `ColocatedData` / `ColocatedCacheEntry` abstraction
- store enough information to reproduce `matched_alleles`, not just “this nearby row may apply”
- add traceability comments on the replacement methods referencing the exact Ensembl VEP GitHub lines for `compare_existing` and the matched-allele flow

### 1.2 Carry shifted and unshifted input state through lookup

Current build-side row data only keeps normalized coordinates plus original start.

Required additions:

- input allele string in VEP representation
- unshifted allele string
- shifted allele string
- unshifted start and end
- strand-aware allele matching state if needed by the upstream algorithm

Required code changes:

- extend [variant_lookup_exec.rs](/Users/mwiewior/research/git/datafusion-bio-functions/datafusion/bio-function-vep/src/variant_lookup_exec.rs#L314)
- change candidate collection so it mirrors [InputBuffer.pm](/Users/mwiewior/research/git/ensembl-vep/modules/Bio/EnsEMBL/VEP/InputBuffer.pm#L312)
- annotate the relevant lookup methods with the exact upstream GitHub permalinks for shifted and unshifted overlap behavior

### 1.3 Port `get_matched_variant_alleles` semantics

This is the core dependency for exact matching.

Tasks:

1. Obtain the exact upstream implementation from Ensembl Variation for release 115 semantics.
2. Reimplement it in Rust without simplification.
3. Validate it against a fixture suite covering:
   - SNVs
   - multi-allelic substitutions
   - insertions
   - deletions
   - left/right trimmed indels
   - repeated-sequence shifted indels
   - null/special alleles

Critical rule:

- do not keep the current normalized-key shortcut as the primary source of truth
- the Rust implementation comment must cite the exact Ensembl Variation or Ensembl VEP GitHub lines that define the matching algorithm

### 1.4 Rebuild `Existing_variation` from matched variants only

Current code exposes all nearby IDs through `existing_variation`.

Required behavior:

- include only variants whose `matched_alleles` contain the current allele or shifted allele, matching [OutputFactory.pm](/Users/mwiewior/research/git/ensembl-vep/modules/Bio/EnsEMBL/VEP/OutputFactory.pm#L1005)
- preserve VEP’s ordering:
  - non-somatic first
  - prefix rank: `rs`, then HGMD-like IDs, then COSMIC-like IDs, then others

Expected impact:

- fix `Existing_variation`
- fix some `SOMATIC` ordering cases

### 1.5 Port AF and MAX_AF exactly

Replace the current logic in [annotate_provider.rs](/Users/mwiewior/research/git/datafusion-bio-functions/datafusion/bio-function-vep/src/annotate_provider.rs#L1811) and [annotate_provider.rs](/Users/mwiewior/research/git/datafusion-bio-functions/datafusion/bio-function-vep/src/annotate_provider.rs#L1856) with VEP-equivalent behavior:

1. Build allele-to-frequency maps per existing variant exactly like [BaseCacheVariation.pm](/Users/mwiewior/research/git/ensembl-vep/modules/Bio/EnsEMBL/VEP/AnnotationSource/Cache/BaseCacheVariation.pm#L185)
2. Apply matched-allele selection exactly like [OutputFactory.pm](/Users/mwiewior/research/git/ensembl-vep/modules/Bio/EnsEMBL/VEP/OutputFactory.pm#L1139)
3. Interpolate only where VEP does:
   - global `AF`
   - exactly one remaining allele
   - biallelic context
4. Compute `MAX_AF` and `MAX_AF_POPS` only from the non-combined populations VEP uses
5. Preserve VEP population names exactly

Expected impact:

- fix `AF`
- fix `AFR_AF`
- fix `AMR_AF`
- fix `EAS_AF`
- fix `EUR_AF`
- fix `SAS_AF`
- fix `gnomADe_AF`
- fix `gnomADg_AF`
- fix `MAX_AF`
- fix `MAX_AF_POPS`

### 1.6 Port clinical metadata filtering exactly

Replace the current filtered-bucket logic with matched-variant logic mirroring [OutputFactory.pm](/Users/mwiewior/research/git/ensembl-vep/modules/Bio/EnsEMBL/VEP/OutputFactory.pm#L1005):

- `CLIN_SIG`
- `SOMATIC`
- `PHENO`
- `PUBMED`

Special cases:

- `clin_sig_allele`
- compound clinical labels
- suppression of all-zero arrays
- duplicate removal rules
- COSMIC/HGMD special handling only if upstream logic actually does it in this path

Expected impact:

- fix `CLIN_SIG`
- fix `SOMATIC`
- fix `PHENO`
- fix `PUBMED`

### 1.7 Post-Phase-1 discrepancy review and plan refresh

Purpose: capture what was actually learned during the existing-variant parity work before moving to transcript consequence parity.

Tasks:

1. Re-run the chr1 benchmark immediately after Phase 1 lands.
2. Compare the new report against the pre-Phase-1 baseline.
3. Update this plan document with a dedicated findings section covering:
   - which mismatch groups were fully resolved
   - which mismatch groups remain
   - any behavior that still depends on a workaround or approximation
   - any divergence discovered between:
     - Rust code and Ensembl VEP
     - Ensembl VEP and underlying Ensembl Variation code
     - cache data shape and what Ensembl VEP expects at runtime
4. For every remaining mismatch group, record:
   - exact affected fields
   - example loci
   - responsible Rust methods
   - corresponding Ensembl VEP source references
   - whether the issue is:
     - missing ported logic
     - incorrect port
     - cache/input data limitation
     - unresolved upstream dependency
5. Explicitly list any workarounds that still exist after Phase 1.
6. If Phase 1 reveals that any earlier assumptions in this plan were wrong, revise the next phases before Phase 2 implementation starts.

Deliverable:

- a committed update to this markdown plan summarizing post-Phase-1 findings

Acceptance:

- the plan is refreshed with Phase-1 findings before Phase 2 begins
- remaining discrepancies and any surviving workarounds are explicitly documented, not inferred from benchmark output alone

### Phase 1 Interim Findings: current implementation pass on March 12, 2026

Status:

- the first strict matched-alleles port landed and the chr1 benchmark was rerun
- this is a real semantic improvement, not a benchmark-specific workaround
- Phase 1 is not complete yet because the existing-variant block is reduced but not eliminated

Validation command:

```bash
cargo run -p datafusion-bio-function-vep --example annotate_vep_golden_bench --release -- \
  vep-benchmark/data/HG002_chr1.vcf.gz \
  /Users/mwiewior/research/data/vep/115_GRCh38_variation_1.parquet \
  parquet \
  0 \
  /tmp/vep_golden \
  vep-benchmark/data/HG002_chr1.vcf.gz \
  /tmp/annotate_vep_golden_bench_phase1a \
  /Users/mwiewior/research/data/vep \
  --steps=ensembl,datafusion \
  --extended-probes
```

Observed result:

- rows still align exactly: `323,430` golden, `323,430` ours
- CSQ entry diff is unchanged: `15` extra transcript entries, `0` missing entries
- perfect fields remain `53/74`
- transcript/FLAGS mismatch set is unchanged
- existing-variant mismatch block dropped sharply:
  - `Existing_variation`: `311,616 -> 25,519`
  - `AF`: `290,664 -> 5,147`
  - `AFR_AF`: `290,664 -> 5,147`
  - `AMR_AF`: `290,655 -> 5,147`
  - `EAS_AF`: `290,613 -> 5,147`
  - `EUR_AF`: `290,663 -> 5,147`
  - `SAS_AF`: `290,664 -> 5,147`
  - `gnomADe_AF`: `44,099 -> 2,300`
  - `gnomADg_AF`: `442,577 -> 25,469`
  - `MAX_AF`: `444,018 -> 26,860`
  - `MAX_AF_POPS`: `442,085 -> 29,636`
  - `CLIN_SIG`: `5,653 -> 118`
  - `SOMATIC`: `3,201 -> 339`
  - `PHENO`: `15,346 -> 1,071`
  - `PUBMED`: `7,931 -> 260`

What the new findings say:

- the matched-alleles data model is the correct direction and resolves the bulk of the co-located divergence
- the remaining Phase 1 errors are mostly over-inclusion, not under-inclusion
- representative failures now look like:
  - VEP emits empty co-located fields, Rust still emits one matched existing variant and its AF/clinical data
  - VEP emits `0&1`, Rust emits `0&0&1`, meaning one extra matched existing variant survives filtering

New discrepancy diagnosis discovered during this pass:

- Rust currently treats every indel as if both:
  - a shifted/minimized input allele string and
  - an unshifted/original input allele string
  should always participate in `compare_existing`
- that is too broad
- upstream VEP output uses `shift_hash->{alt_orig_allele_string}` only when shift metadata exists, and `compare_existing` only consults `unshifted_*` state when that state is actually defined
- current Rust still conflates:
  - ordinary VCF-to-VEP allele minimization
  - true VEP shift-state / unshifted-state semantics

Current conclusion for remaining Phase 1 work:

- Phase 1 must replace the current generic `unshifted_allele_string = raw VCF ref/alt` assumption with the exact upstream condition under which VEP defines shift/unshift metadata
- until that is fixed, the remaining `Existing_variation`, AF, `MAX_AF*`, `CLIN_SIG`, `SOMATIC`, `PHENO`, and `PUBMED` mismatches will not reach zero
- no Phase 2 transcript work should be treated as complete Phase 1 closure

Additional validated findings after later March 12, 2026 reruns:

- full workspace unit tests were rerun and passed after the cache-schema contract change:

```bash
cargo test --workspace --lib
```

- the failed-row filter port produced a much stronger interim state before the parser-space refactor:
  - report: `/tmp/annotate_vep_golden_bench_phase1d/HG002_chr1_0_comparison_report.txt`
  - perfect fields: `63/74`
  - remaining imperfect fields: `11`
  - key existing/co-located fields at that point:
    - `Existing_variation`: `120`
    - `gnomADg_AF`: `120`
    - `MAX_AF_POPS`: `120`
    - `MAX_AF`: `1,680`
    - `gnomADe_AF`: `23`
  - all 1000G AF fields were exact at that stage
  - all clinical fields were exact at that stage

- the parser-style `compare_existing()` refactor then regressed parity:
  - report: `/tmp/annotate_vep_golden_bench_phase1e/HG002_chr1_0_comparison_report.txt`
  - perfect fields: `55/74`
  - remaining imperfect fields: `19`
  - key existing/co-located fields after the refactor:
    - `Existing_variation`: `341`
    - `AF/AFR_AF/AMR_AF/EAS_AF/EUR_AF/SAS_AF`: `15` each
    - `gnomADe_AF`: `2`
    - `gnomADg_AF`: `324`
    - `MAX_AF`: `1,884`
    - `MAX_AF_POPS`: `324`
    - `SOMATIC`: `12`
    - `PHENO`: `12`

- a later parser-space sink/key cleanup restored the refactor to the same `phase1e` benchmark state but did not recover the `phase1d` gains:
  - report: `/tmp/annotate_vep_golden_bench_phase1g/HG002_chr1_0_comparison_report.txt`
  - perfect fields: `55/74`
  - remaining imperfect fields: `19`
  - conclusion:
    - parser-space candidate overlap and per-allele sink keys are necessary to avoid cross-allele contamination
    - those fixes are not sufficient to restore VEP parity on their own

- a later exact serializer/traceability pass removed two non-semantic divergences without changing the remaining colocated-match diagnosis:
  - report: `/tmp/annotate_vep_golden_bench_phase1j/HG002_chr1_0_comparison_report.txt`
  - perfect fields: `56/74`
  - remaining imperfect fields: `18`
  - exact fixes validated on chr1:
    - `FLAGS`: `428 -> 0` after parsing ordered `cds_*` attributes from `raw_object_json` instead of reconstructing canonical order from booleans
    - `MAX_AF`: `1,884 -> 324` after preserving the raw winning frequency string instead of reformatting small values away from VEP's scientific notation
  - remaining colocated failures after this pass are all true matching/allele-identity gaps, not serialization noise:
    - `Existing_variation`: `341`
    - `gnomADg_AF`: `324`
    - `MAX_AF`: `324`
    - `MAX_AF_POPS`: `324`
    - `AF/AFR_AF/AMR_AF/EAS_AF/EUR_AF/SAS_AF`: `15` each
    - `gnomADe_AF`: `2`
    - `SOMATIC`: `12`
    - `PHENO`: `12`

- the next strict compare-space port closed the previously identified output-allele gap and restored a much stronger colocated baseline:
  - report: `/tmp/annotate_vep_golden_bench_phase1j/HG002_chr1_0_comparison_report.txt`
  - annotate time: `92.7s`
  - perfect fields: `64/74`
  - remaining imperfect fields: `10`
  - exact improvements validated on chr1:
    - `Existing_variation`: `341 -> 120`
    - `AF/AFR_AF/AMR_AF/EAS_AF/EUR_AF/SAS_AF`: `15 -> 0`
    - `gnomADg_AF`: `324 -> 120`
    - `MAX_AF`: `324 -> 120`
    - `MAX_AF_POPS`: `324 -> 120`
    - `SOMATIC`: `12 -> 0`
    - `PHENO`: `12 -> 0`
    - `gnomADe_AF`: `2 -> 23`
  - unchanged independent transcript-side gaps:
    - `Consequence`: `15`
    - `INTRON`: `19`
    - `CDS_position`: `2`
    - `Protein_position`: `2`
    - `Amino_acids`: `43`
    - `15` extra transcript CSQ entries
  - current colocated gap is now tightly localized:
    - VEP still suppresses a subset of repeat-shifted existing variants that Rust emits
    - representative loci from the current report:
      - `chr1:105272550` / golden suppresses `rs60126136`
      - `chr1:98448445` / golden suppresses `rs1170278308`
      - `chr1:56092823` / golden suppresses `rs1180305122`
      - `chr1:27971599` / Rust emits `gnomADe_AF=0` where golden stays empty
  - concrete code-level finding from this pass:
    - `annotate_provider.rs` still calls `variant_fields(..., None, ...)` and `frequency_fields(..., None, ...)`
    - the helper methods already model a second allele identity for VEP shift-state filtering, but the live CSQ path never supplies it
    - this confirms that true upstream shift metadata is still missing in the assembled Rust variation-feature state, even though the compare-space match step is now much closer to VEP

New discrepancy diagnosis discovered during the later reruns:

- there are two remaining strict Phase 1 gaps, both confirmed against upstream source and benchmark loci:
- the earlier Gap B diagnosis is now substantially resolved:
  - passing minimized compare-space alleles into `compare_existing()` eliminated the broader downstream output-allele mismatch block
  - the remaining colocated failures are no longer a general parser-vs-output allele identity problem
  - they now point to the narrower missing-shift-metadata problem in the live CSQ assembly path
- there is one remaining strict Phase 1 gap plus the already-separate transcript consequence work:
  - Gap A: real upstream shift-state on the input variation feature is still missing in Rust
    - representative locus: `chr1:203883527` / golden `rs146407839`
    - cache row is a shifted repeat representation at `203883528..203883536` with `allele_string=TGTATTTAT/T`
    - upstream Perl `get_matched_variant_alleles()` returns no match for both:
      - parser-trimmed input `TGTATTTATTTA/-` at `203883528`
      - raw VCF input `TTGTATTTATTTA/T` at `203883527`
    - therefore golden parity at this locus requires a third, shifted input state that Rust does not yet materialize
  - Gap A now has a second validated manifestation in the output layer:
    - representative loci: `chr1:105272550`, `chr1:98448445`, `chr1:56092823`, `chr1:27971599`
    - Rust already stores secondary-allele-aware filtering helpers in `annotate_provider.rs`, but the live CSQ assembly still passes `None`
    - upstream `add_colocated_variant_info()` filters on `hash->{Allele}` plus `vf->{shifted_allele_string}`
    - upstream `add_colocated_frequency_data()` filters on `hash->{Allele}` plus `shift_hash->{alt_orig_allele_string}`
    - until Rust materializes and threads the exact shift-hash state into those calls, a small repeat-shifted colocated mismatch set will survive

Implication for the Phase 1 closure criteria:

- the earlier assumption that the `phase1e` regression was only a sink-key or coordinate-key mismatch was incomplete
- Phase 1 must now port the actual upstream shift-state model, not just parser-space `compare_existing()`
- after that, Phase 1 must thread the exact shift-state alleles into colocated field emission the same way upstream `OutputFactory.pm` does
- no transcript-consequence work should proceed as if Phase 1 were closed until both of those parity gaps are eliminated

Phase 1 acceptance:

- all remaining co-located field mismatches are eliminated
- no new regressions in the previously perfect 53 fields

## Phase 2: Replace Transcript Consequence Approximations with VEP Semantics

Purpose: eliminate the remaining consequence and extra-entry mismatches by removing heuristic logic and porting exact transcript consequence behavior.

### 2.1 Remove heuristic consequence branches

The following logic must stop being authoritative:

- start/stop heuristics
- “safety net” intron stripping at splice sites
- inferred coding consequences when exact coding classification fails
- custom frameshift amino-acid formatting rules not grounded in VEP’s peptide allele string

Relevant current code:

- [transcript_consequence.rs](/Users/mwiewior/research/git/datafusion-bio-functions/datafusion/bio-function-vep/src/transcript_consequence.rs#L1048)
- [transcript_consequence.rs](/Users/mwiewior/research/git/datafusion-bio-functions/datafusion/bio-function-vep/src/transcript_consequence.rs#L1943)

Traceability requirement:

- every replacement method in the transcript consequence path must carry a method comment with the exact Ensembl VEP GitHub permalink for the predicate or serializer behavior it matches

### 2.2 Port transcript overlap and predicate semantics

Tasks:

1. Match VEP transcript hit generation exactly.
2. Match exon/intron membership exactly for:
   - insertions at exon boundaries
   - deletions spanning exon/intron boundaries
   - retained-intron and non-coding biotypes
3. Match splice acceptor/donor windows exactly.
4. Match polypyrimidine tract behavior exactly.
5. Match when `intron_variant` co-emits with splice terms and when it does not.
6. Match UTR co-emission behavior for boundary indels.

Expected impact:

- fix `Consequence`
- fix `INTRON`
- remove the `15` extra CSQ entries

### 2.3 Port exact coding position semantics

Tasks:

1. Port CDS position and protein position behavior from the transcript-variation layer, not from local formatting rules.
2. Match boundary-shift handling equivalent to:
   - `tv->{_boundary_shift}`
   - `format_coords`
   - translation start/end logic
3. Match incomplete CDS formatting for `?-N`
4. Ensure coding position output remains correct when consequence classification falls on splice/coding boundaries

Expected impact:

- fix `CDS_position`
- fix `Protein_position`

### 2.4 Port exact amino-acid serialization

Tasks:

1. Reproduce VEP `pep_allele_string` output semantics rather than deriving a local simplified string.
2. Validate exact formatting for:
   - frameshift insertions
   - frameshift deletions
   - preserved first affected amino acid cases
   - codon-boundary indels
   - retained stop/start cases

Expected impact:

- fix `Amino_acids`

Phase 2 acceptance:

- `Consequence`, `INTRON`, `CDS_position`, `Protein_position`, and `Amino_acids` all reach zero mismatches
- `ours_only == 0`

## Phase 3: Preserve Exact FLAGS Ordering

Purpose: remove the final non-semantic serializer mismatch.

Tasks:

1. Stop relying on boolean reconstruction for `FLAGS` whenever ordered transcript attributes are not preserved.
2. Prefer one of:
   - an ordered attribute list supplied by the cache format
   - ordered extraction from `raw_object_json`
3. Match VEP’s emitted attribute order exactly, including the case where the same logical flags appear in different orders across transcripts.

Relevant current file:

- [transcript_consequence.rs](/Users/mwiewior/research/git/datafusion-bio-functions/datafusion/bio-function-vep/src/transcript_consequence.rs#L3012)

Expected impact:

- fix `FLAGS`

Phase 3 acceptance:

- `FLAGS` reaches zero mismatches on chr1

## Phase 4: Close the HGVS Parity Gap

Purpose: reach real VEP parity beyond the current README benchmark.

Tasks:

1. Enable a second golden benchmark mode that runs VEP with:
   - `--hgvs`
   - `--fasta`
2. Stop hardcoding blank `HGVSc` and `HGVSp` in transcript CSQ assembly.
3. Use the already-computed transcript consequence HGVS fields when enabled.
4. If the current HGVS implementation differs from VEP:
   - port strand handling
   - port transcript/protein accession formatting
   - port offset handling
   - port escaping behavior

Relevant current files:

- [annotate_vep_golden_bench.rs](/Users/mwiewior/research/git/datafusion-bio-functions/datafusion/bio-function-vep/examples/annotate_vep_golden_bench.rs#L514)
- [annotate_provider.rs](/Users/mwiewior/research/git/datafusion-bio-functions/datafusion/bio-function-vep/src/annotate_provider.rs#L2028)
- [OutputFactory.pm](/Users/mwiewior/research/git/ensembl-vep/modules/Bio/EnsEMBL/VEP/OutputFactory.pm#L1678)

Acceptance:

- HGVS benchmark reaches zero mismatches when HGVS is enabled

## Traceability Comment Format

Required minimum format for parity-sensitive methods:

```rust
/// VEP parity:
/// - https://github.com/Ensembl/ensembl-vep/blob/<commit>/modules/Bio/EnsEMBL/VEP/<file>.pm#Lx-Ly
/// - Contract: short statement of the exact VEP behavior implemented here.
```

If multiple upstream sources are required:

```rust
/// VEP parity:
/// - https://github.com/Ensembl/ensembl-vep/blob/<commit>/modules/Bio/EnsEMBL/VEP/<file>.pm#Lx-Ly
///   Role: candidate selection
/// - https://github.com/Ensembl/ensembl-vep/blob/<commit>/modules/Bio/EnsEMBL/VEP/<file>.pm#La-Lb
///   Role: output serialization
```

Rules:

- use GitHub blob permalinks
- pin to the exact commit or release tag under implementation
- include file and line ranges
- keep the contract statement short and specific
- update the traceability comment whenever the Rust method’s upstream source mapping changes

## Testing Plan

### Unit tests

Add exact fixture tests for the currently mismatching loci:

- `1:182297818`
- `1:225123725`
- `1:212909295`
- `1:235470575`
- `1:167426237`
- `1:15178152`
- `1:44206707`
- `1:9244919`
- `1:56249200`
- `1:175160811`
- `1:225524663`
- `1:101495223`

Each fixture should assert:

- exact term list
- exact `Existing_variation`
- exact AF fields
- exact clinical fields
- exact transcript count

### Property-style tests

For allele matching and frequency selection:

- shifted vs unshifted indels
- repeated-sequence indels
- allele-string permutations
- biallelic interpolation for global `AF`
- no interpolation for gnomAD sub-populations

### Full-run regression tests

Run after each phase:

1. chr22 1000-variant benchmark
2. full chr1 benchmark
3. update this plan with findings after Phase 1 before proceeding
4. optional HGVS-enabled benchmark once Phase 4 starts

## Acceptance Criteria

The strict parity work is complete only when all of the following hold:

1. README chr1 benchmark:
   - `0` field mismatches across all 74 columns
   - `0` extra entries
   - `0` missing entries
2. No heuristic consequence paths remain in the hot path for cases VEP defines exactly.
3. Existing/co-located output is driven by explicit `matched_alleles` semantics.
4. Every parity-sensitive method has a review record and a method-level traceability comment with exact Ensembl VEP GitHub file and line references.
5. `FLAGS` ordering is sourced from ordered upstream-equivalent data.
6. A second HGVS-enabled benchmark path exists for validating real full-field parity.

## Recommended Execution Order

1. Phase -1: existing code review and traceability inventory
2. Phase 0: validation hardening
3. Phase 1: exact existing/co-located matching and frequency logic
4. Phase 2: exact transcript consequence and transcript-count parity
5. Phase 3: exact `FLAGS`
6. Phase 4: HGVS parity beyond the current README procedure

## Expected Outcome by Phase

- After Phase 1:
  - 15 large co-located field mismatch groups should collapse
- After Phase 2:
  - the 5 consequence/position/amino-acid fields should collapse
  - the 15 extra entries should disappear
- After Phase 3:
  - `FLAGS` should collapse
- After Phase 4:
  - parity is no longer limited by the current benchmark’s blank HGVS fields

## Non-Goals

- optimizing runtime before semantic parity is reached
- preserving current workaround behavior for backward compatibility
- adding new consequence features unrelated to Ensembl VEP 115 parity
