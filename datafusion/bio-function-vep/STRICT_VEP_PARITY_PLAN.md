# Strict Ensembl VEP Parity Plan

Goal: get `datafusion/bio-function-vep` to zero mismatches on the README chr1 golden benchmark for all 74 CSQ columns, with zero extra CSQ entries, by implementing Ensembl VEP semantics exactly rather than tuning heuristics.

Date: March 13, 2026

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

## Current Verified Status

Validated on this checkout with:

```bash
cargo run -p datafusion-bio-function-vep --example annotate_vep_golden_bench --release -- \
  vep-benchmark/data/HG002_chr1.vcf.gz \
  /Users/mwiewior/research/data/vep/115_GRCh38_variation_1.parquet \
  parquet \
  0 \
  /tmp/vep_golden \
  vep-benchmark/data/HG002_chr1.vcf.gz \
  /tmp/annotate_vep_golden_bench_release_entry_diag_fix2 \
  /Users/mwiewior/research/data/vep \
  --steps=ensembl,datafusion \
  --extended-probes
```

Observed result:

- `323,430` variants compared
- `2,997,504` CSQ entries compared
- `74/74` fields perfect
- `0` extra CSQ entries on our side
- `0` missing CSQ entries
- `100%` most-severe parity
- `100%` term-set parity
- release annotate time `49.783s`

**Phase 4 HGVS parity: ACHIEVED** (March 14, 2026)

- Non-merged (Ensembl-only) chr1 benchmark: **74/74 fields at zero mismatches** with `--hgvs --fasta`
- Merged (Ensembl+RefSeq) chr1 benchmark: 72/74 fields, ~105 remaining mismatches (all RefSeq-specific)
- Performance: 51s without HGVS, 79s with HGVS (matches pre-HGVS baseline)

Remaining open work:

- merged (RefSeq) HGVS parity — edited transcript shifting, cdna_mapper_segments
- release-115 flag-surface parity beyond the fixed README field set
- VEP-style on/off gating for feature-expanding and output-gating flags
- dynamic CSQ field/header parity for flags not covered by the current 74-column benchmark

Important caveat:

- the current `74/74` result proves parity for the non-merged benchmark profile, not full release-115 CLI flag parity
- the current Rust path still behaves like a fixed-schema benchmark serializer in several places where VEP conditionally adds or removes columns based on flags

## Reproducible Benchmark Commands

Important harness behavior:

- `annotate_vep_golden_bench` only skips Docker VEP when the expected golden VCF already exists inside the chosen `work_dir`.
- For full `chr1` runs with `sample_limit=0`, the expected filename is:
  - `HG002_chr1_0_vep115_golden.vcf`
- So the reusable workflow is:
  - create a dedicated `work_dir`
  - copy the correct golden VCF into that directory with that exact filename
  - run the benchmark with `--steps=datafusion`

### 1. README chr1 benchmark against the saved default golden

```bash
mkdir -p /tmp/annotate_vep_default_compare_repo_golden
cp \
  /Users/mwiewior/research/git/datafusion-bio-functions/vep-benchmark/data/output/default/HG002_chr1_0_vep115_golden.vcf \
  /tmp/annotate_vep_default_compare_repo_golden/HG002_chr1_0_vep115_golden.vcf

cargo run --release -p datafusion-bio-function-vep --example annotate_vep_golden_bench -- \
  /Users/mwiewior/research/git/datafusion-bio-functions/vep-benchmark/data/HG002_chr1.vcf.gz \
  /Users/mwiewior/research/data/vep/115_GRCh38_variation_1.parquet \
  parquet \
  0 \
  /Users/mwiewior/research/data/vep/homo_sapiens/115_GRCh38 \
  /Users/mwiewior/research/git/datafusion-bio-functions/vep-benchmark/data/HG002_chr1.vcf.gz \
  /tmp/annotate_vep_default_compare_repo_golden \
  /Users/mwiewior/research/data/vep \
  --steps=datafusion \
  --extended-probes
```

Expected artifact paths:

- `/tmp/annotate_vep_default_compare_repo_golden/HG002_chr1_0_comparison_report.txt`
- `/tmp/annotate_vep_default_compare_repo_golden/HG002_chr1_0_discrepancies.txt`

### 2. HGVS non-merged chr1 benchmark against the saved golden

```bash
mkdir -p /tmp/annotate_vep_hgvs_nonmerge_compare
cp \
  /Users/mwiewior/research/git/datafusion-bio-functions/vep-benchmark/data/output/hgvsc_nonmerge/HG002_chr1_0_vep115_golden.vcf \
  /tmp/annotate_vep_hgvs_nonmerge_compare/HG002_chr1_0_vep115_golden.vcf

cargo run --release -p datafusion-bio-function-vep --example annotate_vep_golden_bench -- \
  --hgvs \
  --extended-probes \
  --reference-fasta-path=/tmp/ensembl_release115_official_fasta/Homo_sapiens.GRCh38.dna.toplevel.fa.gz \
  /Users/mwiewior/research/git/datafusion-bio-functions/vep-benchmark/data/HG002_chr1.vcf.gz \
  /Users/mwiewior/research/data/vep/115_GRCh38_variation_1.parquet \
  parquet \
  0 \
  /Users/mwiewior/research/data/vep/homo_sapiens/115_GRCh38 \
  /Users/mwiewior/research/git/datafusion-bio-functions/vep-benchmark/data/HG002_chr1.vcf.gz \
  /tmp/annotate_vep_hgvs_nonmerge_compare \
  /Users/mwiewior/research/data/vep \
  --steps=datafusion
```

### 3. HGVS merged chr1 benchmark against the saved golden

```bash
mkdir -p /tmp/annotate_vep_hgvs_merged_compare
cp \
  /Users/mwiewior/research/git/datafusion-bio-functions/vep-benchmark/data/output/hgvsc_merged/HG002_chr1_0_vep115_golden.vcf \
  /tmp/annotate_vep_hgvs_merged_compare/HG002_chr1_0_vep115_golden.vcf

cargo run --release -p datafusion-bio-function-vep --example annotate_vep_golden_bench -- \
  --merged \
  --hgvs \
  --extended-probes \
  --reference-fasta-path=/tmp/ensembl_release115_official_fasta/Homo_sapiens.GRCh38.dna.toplevel.fa.gz \
  /Users/mwiewior/research/git/datafusion-bio-functions/vep-benchmark/data/HG002_chr1.vcf.gz \
  /Users/mwiewior/research/data/vep/115_GRCh38_variation_1.parquet \
  parquet \
  0 \
  /Users/mwiewior/research/data/vep/homo_sapiens_merged/115_GRCh38 \
  /Users/mwiewior/research/git/datafusion-bio-functions/vep-benchmark/data/HG002_chr1.vcf.gz \
  /tmp/annotate_vep_hgvs_merged_compare \
  /Users/mwiewior/research/data/vep \
  --steps=datafusion
```

### 4. Fresh Docker VEP golden generation for a new work directory

Only use this when intentionally regenerating the golden VCF.

```bash
cargo run --release -p datafusion-bio-function-vep --example annotate_vep_golden_bench -- \
  --merged \
  --hgvs \
  --extended-probes \
  --reference-fasta-path=/tmp/ensembl_release115_official_fasta/Homo_sapiens.GRCh38.dna.toplevel.fa.gz \
  /Users/mwiewior/research/git/datafusion-bio-functions/vep-benchmark/data/HG002_chr1.vcf.gz \
  /Users/mwiewior/research/data/vep/115_GRCh38_variation_1.parquet \
  parquet \
  0 \
  /Users/mwiewior/research/data/vep/homo_sapiens_merged/115_GRCh38 \
  /Users/mwiewior/research/git/datafusion-bio-functions/vep-benchmark/data/HG002_chr1.vcf.gz \
  /tmp/annotate_vep_hgvs_merged_fresh \
  /Users/mwiewior/research/data/vep \
  --steps=ensembl,datafusion
```

### 5. Unit test commands used while iterating on HGVS parity

```bash
cargo test -p datafusion-bio-function-vep --lib format_hgvsc_
cargo test -p datafusion-bio-function-vep --lib test_format_hgvsp_
cargo test -p datafusion-bio-function-vep --lib
cargo test --workspace --lib
```

## Release 115 Flag-Parity Audit

This section supersedes the earlier rough split between "flags with side effects" and "pure output flags". For Ensembl VEP release 115, the exact target is:

- match which flags change processing or row generation
- match which flags imply other flags
- match which flags only gate output fields
- match which CSQ columns exist at all for a given flag set

### Corrections to the Initial Categorisation

- `--pubmed` is not a pure output-only flag in release 115. `Config.pm` groups `pubmed` with `check_existing`, `af`, `af_1kg`, `af_gnomad`, `af_gnomade`, `af_gnomadg`, and `max_af`, all of which imply `--check_existing`.
- `--protein` is primarily an output-gating flag in release 115. It controls `ENSP` emission in `OutputFactory`; it does not independently change transcript consequence classification.
- `--canonical` only becomes selection-sensitive when VEP pick-order modes are active. In the current Rust parity scope, it should be treated as output gating unless and until `pick`, `pick_allele`, `per_gene`, or related modes are implemented.
- `--mirna` is not just a generic annotation column switch. It enables miRNA structural output in `OutputFactory` and depends on transcript attributes and overlap state.
- `--af_esp` should not be treated as part of the current release-115 target surface unless explicitly added. It still appears in some compatibility/incompatibility code paths, but it is not part of the active release-115 option set we are currently matching.

### Flag Families That Must Match VEP Semantics

#### Processing-side and row-set-changing flags

- `--hgvs`, `--hgvsc`, `--hgvsp`, `--shift_hgvs`
  - VEP effect:
    - enables transcript/protein HGVS generation
    - requires FASTA in offline mode
    - changes shifted allele/coordinate behavior through HGVS-specific 3-prime shifting
  - Current Rust status:
    - partially supported
    - wrapper-level flag parsing and `shift_hgvs` control exist
    - exact Ensembl Variation HGVS algorithm is still incomplete
- `--regulatory`
  - VEP effect:
    - expands the overlap set beyond transcripts
    - adds regulatory-feature and motif-feature consequence rows
    - can also be enabled indirectly when a plugin requests regulatory or motif features
  - Current Rust status:
    - regulatory, motif, and miRNA context engines exist
    - benchmark path already exercises regulatory output
    - explicit top-level VEP-style `regulatory` option gating is not yet fully modeled
- `--check_existing`, `--af`, `--af_1kg`, `--af_gnomad`, `--af_gnomade`, `--af_gnomadg`, `--max_af`, `--pubmed`
  - VEP effect:
    - turns on co-located/existing-variant lookup
    - enables AF/clinical/publication output from matched existing variants
    - `af_*`, `max_af`, and `pubmed` all imply `check_existing`
  - Current Rust status:
    - exact co-located matching semantics are in place for the README benchmark profile
    - implication handling exists for `af`, `af_1kg`, `af_gnomade`, `af_gnomadg`, `max_af`, and `pubmed`
    - release-115 alias parity for `af_gnomad` still needs to be added

#### Output-gating flags already represented in the current 74-field schema

- `--symbol`, `--biotype`, `--numbers`, `--ccds`, `--canonical`, `--mane`, `--mane_select`, `--tsl`, `--protein`, `--uniprot`, `--variant_class`, `--gene_phenotype`
  - VEP effect:
    - add or remove specific CSQ fields from the header and output
    - do not change consequence logic on their own
  - Current Rust status:
    - most corresponding columns already exist in the fixed 74-field schema
    - current serializer still emits many of them unconditionally in benchmark mode instead of gating them by option
    - this is benchmark parity, not yet CLI flag parity

#### Output-gating flags whose columns are still missing from the current Rust CSQ schema

- `--appris` -> `APPRIS`
- `--domains` -> `DOMAINS`
- `--sift b` -> `SIFT`
- `--polyphen b` -> `PolyPhen`
- `--mirna` -> `miRNA`

Current Rust status:

- these are not part of the current 74-column schema in `golden_benchmark.rs`
- some of the underlying data or logic already exists
  - miRNA overlap logic exists in the transcript engine
  - protein-feature overlap support exists upstream but is not yet serialized as `DOMAINS`
- output/header parity for these fields is still open

### Source-of-Truth References for Flag Behavior

- `Config.pm`
  - `check_existing`, AF flags, `max_af`, `pubmed`: `modules/Bio/EnsEMBL/VEP/Config.pm#L124-L136`
  - transcript/output flags: `modules/Bio/EnsEMBL/VEP/Config.pm#L177-L220`
  - default "everything" profile enabling `hgvs`, `symbol`, `numbers`, `domains`, `regulatory`, `canonical`, `protein`, `biotype`, `uniprot`, `tsl`, `variant_class`, `sift`, `polyphen`, `ccds`, `gene_phenotype`, `mane`, `mirna`: `modules/Bio/EnsEMBL/VEP/Config.pm#L350-L381`
  - flag implication rules: `modules/Bio/EnsEMBL/VEP/Config.pm#L439-L489`
  - pick-order interaction for `canonical`, `appris`, `tsl`, `biotype`, `ccds`, `mane`: `modules/Bio/EnsEMBL/VEP/Config.pm#L305-L306`
- `Runner.pm`
  - offline HGVS gating and `shift_hgvs=0` semantics: `modules/Bio/EnsEMBL/VEP/Runner.pm#L726-L738`, `modules/Bio/EnsEMBL/VEP/Runner.pm#L771-L773`
  - plugin-driven regulatory enablement: `modules/Bio/EnsEMBL/VEP/Runner.pm#L959-L962`
- `BaseRunner.pm`
  - process-wide `shift_hgvs` default wiring: `modules/Bio/EnsEMBL/VEP/BaseRunner.pm#L491-L496`
- `OutputFactory.pm`
  - numbers, domains, symbol, CCDS, ENSP, uniprot, canonical, biotype, gene phenotype, MANE, TSL, APPRIS, miRNA: `modules/Bio/EnsEMBL/VEP/OutputFactory.pm#L1434-L1610`
  - HGVS, SIFT, PolyPhen output gating: `modules/Bio/EnsEMBL/VEP/OutputFactory.pm#L1693-L1794`

### Current Rust Gaps Relative to That Surface

- fixed 74-field CSQ schema is still hard-coded in `golden_benchmark.rs`
- many benchmark-profile transcript fields are emitted even when the corresponding VEP flag would normally be off
- there is no complete release-115 option parser for all output-gating flags listed above
- there is no exact field-list/header parity layer yet for non-benchmark flag combinations
- missing output columns still need serializer, tests, and benchmark coverage:
  - `APPRIS`
  - `DOMAINS`
  - `SIFT`
  - `PolyPhen`
  - `miRNA`

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

### HGVS semantics used by VEP

- `modules/Bio/EnsEMBL/VEP/Config.pm`
  - `hgvs`
  - `hgvsc`
  - `hgvsp`
  - `hgvsp_use_prediction`
  - `remove_hgvsp_version`
  - `shift_hgvs`
- `modules/Bio/EnsEMBL/VEP/Runner.pm`
  - offline HGVS + FASTA gating
  - `shift_hgvs=0` forcing `shift_3prime=0`
- `modules/Bio/EnsEMBL/VEP/BaseRunner.pm`
  - offline/cache-backed HGVS shift defaults
- `modules/Bio/EnsEMBL/VEP/OutputFactory.pm`
  - HGVS emission
  - URI escaping
  - output-level version stripping / predicted-format toggles
- external dependency:
  - `Bio::EnsEMBL::Variation::TranscriptVariationAllele`
    - `hgvs_transcript`
    - `hgvs_protein`
    - `_return_3prime`
    - `_clip_alleles`
    - `_get_cDNA_position`
    - `_get_hgvs_protein_type`
    - `_get_hgvs_peptides`
    - `_get_fs_peptides`
    - `_get_alternate_cds`
    - `_stop_loss_extra_AA`
  - `Bio::EnsEMBL::Variation::Utils::Sequence`
    - `hgvs_variant_notation`
    - `format_hgvs_string`

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

### 4. ~~The current benchmark does not prove HGVS parity~~ RESOLVED

**HGVS parity achieved on non-merged benchmark** (March 14, 2026).

The non-merged chr1 benchmark with `--hgvs --fasta` now shows **74/74 fields at zero mismatches** across 2,997,504 CSQ entries. The HGVS generation path ports exact Ensembl Variation algorithms including `hgvs_transcript()`, `hgvs_protein()`, `perform_shift()`, `_stop_loss_extra_AA()`, `_check_for_peptide_duplication()`, `_get_alternate_cds()`, and `genomic2pep()` with full traceability to VEP release/115 source.

Relevant files:

- [annotate_vep_golden_bench.rs](/Users/mwiewior/research/git/datafusion-bio-functions/datafusion/bio-function-vep/examples/annotate_vep_golden_bench.rs#L514)
- [annotate_provider.rs](/Users/mwiewior/research/git/datafusion-bio-functions/datafusion/bio-function-vep/src/annotate_provider.rs#L1176)
- [annotate_provider.rs](/Users/mwiewior/research/git/datafusion-bio-functions/datafusion/bio-function-vep/src/annotate_provider.rs#L2338)
- [annotate_table_function.rs](/Users/mwiewior/research/git/datafusion-bio-functions/datafusion/bio-function-vep/src/annotate_table_function.rs#L1850)
- [hgvs.rs](/Users/mwiewior/research/git/datafusion-bio-functions/datafusion/bio-function-vep/src/hgvs.rs#L48)
- [hgvs.rs](/Users/mwiewior/research/git/datafusion-bio-functions/datafusion/bio-function-vep/src/hgvs.rs#L290)
- [variant_lookup_exec.rs](/Users/mwiewior/research/git/datafusion-bio-functions/datafusion/bio-function-vep/src/variant_lookup_exec.rs#L907)
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

1. Finish the HGVS harness and setup validation first:
   - enable the second golden benchmark mode that runs VEP with `--hgvs --fasta`
   - require an indexed reference FASTA at setup time rather than failing later during execution
   - make the HGVS integration test use a valid self-contained indexed FASTA fixture
2. Keep wrapper semantics aligned with VEP:
   - `hgvs` enables both `hgvsc` and `hgvsp`
   - offline/cache-backed HGVS requires FASTA
   - `shift_hgvs=0` disables HGVS-only 3-prime shifting by forcing `shift_3prime=0`
   - `remove_hgvsp_version`, `hgvsp_use_prediction`, and output-layer escaping follow VEP behavior exactly
3. Replace the simplified transcript HGVS formatter with a source-equivalent Ensembl Variation port:
   - port `hgvs_variant_notation`
   - port `format_hgvs_string`
   - port `_clip_alleles`
   - port `_get_cDNA_position`
   - port transcript numbering (`c.` vs `n.`), intronic offsets, 5' UTR negatives, and 3' UTR `*` coordinates
   - port shifted insertion/deletion handling via `_return_3prime` / `shift_hash`
4. Replace the simplified protein HGVS formatter with a source-equivalent Ensembl Variation port:
   - port `_get_hgvs_protein_type`
   - port `_get_hgvs_peptides`
   - port `_get_fs_peptides`
   - port `_get_alternate_cds`
   - port `_stop_loss_extra_AA`
   - port peptide duplication / delins / frameshift / stop-loss / start-loss formatting
5. Add upstream-style regression coverage before relying on full chr1:
   - synonymous `=`
   - intronic HGVS on both strands
   - 5' UTR and 3' UTR transcript coordinates
   - insertion vs duplication
   - deletion vs delins
   - frameshift with exact first changed AA
   - stop-loss extension
   - version stripping / predicted-format toggles / output escaping

Relevant current files:

- [annotate_vep_golden_bench.rs](/Users/mwiewior/research/git/datafusion-bio-functions/datafusion/bio-function-vep/examples/annotate_vep_golden_bench.rs#L514)
- [annotate_provider.rs](/Users/mwiewior/research/git/datafusion-bio-functions/datafusion/bio-function-vep/src/annotate_provider.rs#L1176)
- [annotate_provider.rs](/Users/mwiewior/research/git/datafusion-bio-functions/datafusion/bio-function-vep/src/annotate_provider.rs#L2338)
- [annotate_table_function.rs](/Users/mwiewior/research/git/datafusion-bio-functions/datafusion/bio-function-vep/src/annotate_table_function.rs#L1850)
- [hgvs.rs](/Users/mwiewior/research/git/datafusion-bio-functions/datafusion/bio-function-vep/src/hgvs.rs#L48)
- [hgvs.rs](/Users/mwiewior/research/git/datafusion-bio-functions/datafusion/bio-function-vep/src/hgvs.rs#L290)
- [variant_lookup_exec.rs](/Users/mwiewior/research/git/datafusion-bio-functions/datafusion/bio-function-vep/src/variant_lookup_exec.rs#L907)
- [OutputFactory.pm](/Users/mwiewior/research/git/ensembl-vep/modules/Bio/EnsEMBL/VEP/OutputFactory.pm#L1678)

Acceptance:

- HGVS integration tests pass with a real indexed FASTA fixture
- HGVS benchmark reaches zero mismatches when HGVS is enabled
- no simplified HGVS-only heuristic formatter remains in the hot path

## Phase 5: Release-115 Flag Surface and Output-Schema Parity

Purpose: move from "README benchmark parity" to "VEP option-surface parity" for the supported release-115 flags.

Tasks:

1. Add a real release-115 option matrix to `options_json` parsing and benchmark harness wiring:
   - `regulatory`
   - `symbol`
   - `biotype`
   - `numbers`
   - `ccds`
   - `canonical`
   - `protein`
   - `uniprot`
   - `mane`
   - `mane_select`
   - `tsl`
   - `appris`
   - `variant_class`
   - `gene_phenotype`
   - `mirna`
   - `domains`
   - `sift`
   - `polyphen`
   - release-115 alias `af_gnomad` -> `af_gnomade`
2. Split processing-side flags from output-gating flags in the Rust implementation:
   - processing-side behavior must change lookup/context loading exactly when VEP does
   - output-only flags must not silently stay "always on" in serializer code
3. Make CSQ field/header generation flag-sensitive like VEP rather than permanently fixed to the README 74-column contract:
   - disabled fields must disappear from the emitted CSQ field order when VEP would omit them
   - enabled fields must appear in VEP order
   - comparison tooling must support multiple field profiles, not only the fixed 74-field profile
4. Implement missing release-115 output columns:
   - `APPRIS`
   - `DOMAINS`
   - `SIFT`
   - `PolyPhen`
   - `miRNA`
5. Add exact integration benchmarks against Docker VEP for at least these profiles:
   - default non-HGVS profile
   - `--regulatory`
   - `--check_existing --af --af_1kg --af_gnomade --af_gnomadg --max_af --pubmed`
   - transcript-output profile:
     - `--symbol --biotype --numbers --ccds --canonical --protein --uniprot --mane --tsl --variant_class --gene_phenotype`
   - extra output profile:
     - `--appris --domains --sift b --polyphen b --mirna`
   - HGVS profile:
     - `--hgvs --fasta`
     - `--hgvs --shift_hgvs 0 --fasta`

Acceptance:

- for every supported profile, our CSQ header field list matches VEP exactly
- enabling/disabling a flag changes row count, field set, and coordinates exactly when VEP does
- no benchmark-only always-on serializer path remains for fields controlled by VEP flags
- the missing release-115 columns above are emitted and tested

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
7. The release-115 flag implication matrix is matched for all supported flags in scope.
8. Output-gating flags add or remove CSQ columns exactly as VEP does, rather than relying on a permanently fixed schema.
9. Missing release-115 output fields in the audited scope (`APPRIS`, `DOMAINS`, `SIFT`, `PolyPhen`, `miRNA`) are implemented and benchmarked.

Status on March 13, 2026:

- Criteria `1` through `5` are satisfied on chr1.
- Criterion `6` remains open.
  - wrapper-level HGVS emission and indexed-FASTA validation are in place
  - exact Ensembl Variation HGVS generation is still incomplete
  - the HGVS-enabled chr1 benchmark still has real mismatches
- Criteria `7` through `9` are not yet satisfied.
  - flag implication coverage is still partial
  - fixed-schema benchmark output still masks missing output-gating parity
  - several release-115 output columns are still absent from the Rust CSQ schema

## Recommended Execution Order

1. Phase -1: existing code review and traceability inventory
2. Phase 0: validation hardening
3. Phase 1: exact existing/co-located matching and frequency logic
4. Phase 2: exact transcript consequence and transcript-count parity
5. Phase 3: exact `FLAGS`
6. Phase 4: HGVS parity beyond the current README procedure
7. Phase 5: release-115 flag surface and dynamic output-schema parity

## Expected Outcome by Phase

- After Phase 1:
  - completed
- After Phase 2:
  - completed
- After Phase 3:
  - completed
- After Phase 4:
  - parity is no longer limited by the current benchmark’s unexercised or simplified HGVS path
- After Phase 5:
  - the implementation can be compared against Ensembl VEP under multiple flag profiles rather than only the fixed README profile

## Phase 4 HGVS Parity Status (updated March 14, 2026)

### Non-merged (Ensembl-only) benchmark — near-zero

Benchmark: chr1 non-merged, 2,997,504 CSQ entries, `--hgvs --fasta`, Docker VEP 115.2.
**72/74 fields at zero mismatches. HGVSc: 2, HGVSp: 4. Total: 6 mismatches (99.9998%).**

Remaining 6 mismatches:
- 2 HGVSc: spurious HGVSc where VEP emits empty (ENST00000948304, ENST00000673836)
- 3 HGVSp: `extTerN` where VEP says `extTer?` (VEP `_trim_incomplete_codon` trims the alt CDS before UTR append)
- 1 HGVSp: insertion flanking off-by-2 at chr1:248638949 (protein position mapper difference)

### Merged (Ensembl+RefSeq) benchmark

Benchmark: chr1 merged, 4,737,090 CSQ entries, `--hgvs --shift_hgvs 1 --fasta`.
72/74 fields at zero mismatches. Only HGVSc (106) and HGVSp (5) remain.

### Progression

| Run | HGVSc mismatches | HGVSp mismatches | Total |
|-----|-----------------|-----------------|-------|
| mapper_fix (earliest) | 16,215 | 114 | 16,329 |
| resample | 5,636 | 114 | 5,750 |
| resume1 | 584 | 114 | 698 |
| resume2 | 860 | 114 | 974 |
| **Current** | **106** | **125** | **231** |

### Fixes applied

1. **Genomic shift `seq_strand=1`** — VEP `_genomic_shift()` always passes `seq_strand=1` to `perform_shift()`. We were passing the transcript strand, causing `hgvs_reverse=true` for reverse-strand transcripts and wrong allele rotation.
   - [TranscriptVariationAllele.pm#L431](https://github.com/Ensembl/ensembl-variation/blob/release/115/modules/Bio/EnsEMBL/Variation/TranscriptVariationAllele.pm#L431)
2. **3' UTR extension for stop-loss/frameshift** — VEP `_get_alternate_cds()` appends `_three_prime_utr()` before translating so `_stop_loss_extra_AA()` can find the new stop codon in the UTR.
   - [TranscriptVariationAllele.pm#L2335-L2372](https://github.com/Ensembl/ensembl-variation/blob/release/115/modules/Bio/EnsEMBL/Variation/TranscriptVariationAllele.pm#L2335-L2372)
   - [TranscriptVariationAllele.pm#L2406-L2461](https://github.com/Ensembl/ensembl-variation/blob/release/115/modules/Bio/EnsEMBL/Variation/TranscriptVariationAllele.pm#L2406-L2461)
3. **Protein dup before shift** — VEP `hgvs_protein()` calls `_check_for_peptide_duplication()` BEFORE `_shift_3prime()`. We had the order reversed.
   - [TranscriptVariationAllele.pm#L1700-L1758](https://github.com/Ensembl/ensembl-variation/blob/release/115/modules/Bio/EnsEMBL/Variation/TranscriptVariationAllele.pm#L1700-L1758)
4. **Transcript-level shift via `perform_shift`** — VEP `_return_3prime()` calls `perform_shift()` with `hgvs_reverse` flag for edited RefSeq transcripts. We were using simple `rotate_left`/`rotate_right`.
   - [TranscriptVariationAllele.pm#L237-L243](https://github.com/Ensembl/ensembl-variation/blob/release/115/modules/Bio/EnsEMBL/Variation/TranscriptVariationAllele.pm#L237-L243)
   - [TranscriptVariationAllele.pm#L291-L351](https://github.com/Ensembl/ensembl-variation/blob/release/115/modules/Bio/EnsEMBL/Variation/TranscriptVariationAllele.pm#L291-L351)

### ACHIEVED: 74/74 zero mismatches on non-merged (Ensembl-only) benchmark

**Date: March 14, 2026**

Non-merged (Ensembl-only) chr1 benchmark with Docker VEP 115.2:
- **74/74 CSQ fields at zero mismatches**
- **2,997,504 CSQ entries compared**
- **100.000% accuracy across all fields**
- **0 missing entries, 0 extra entries**

Performance (chr1, 323K variants, VEP-only cache with promoted columns):
- Without HGVS: **51s** (matches pre-HGVS baseline)
- With HGVS: **79s** (+55% HGVS overhead, inherent cost of coordinate mapping for 3M CSQ entries)

### Fixes applied (with Ensembl VEP traceability)

1. **Genomic shift `seq_strand=1`** — VEP `_genomic_shift()` always passes `seq_strand=1` to `perform_shift()`
   - [TranscriptVariationAllele.pm#L431](https://github.com/Ensembl/ensembl-variation/blob/release/115/modules/Bio/EnsEMBL/Variation/TranscriptVariationAllele.pm#L431)
2. **3' UTR extension for stop-loss/frameshift** — VEP `_get_alternate_cds()` appends `_three_prime_utr()` before translating
   - [TranscriptVariationAllele.pm#L2335-L2372](https://github.com/Ensembl/ensembl-variation/blob/release/115/modules/Bio/EnsEMBL/Variation/TranscriptVariationAllele.pm#L2335-L2372)
3. **Protein dup before shift** — VEP `hgvs_protein()` calls `_check_for_peptide_duplication()` BEFORE `_shift_3prime()`
   - [TranscriptVariationAllele.pm#L1700-L1758](https://github.com/Ensembl/ensembl-variation/blob/release/115/modules/Bio/EnsEMBL/Variation/TranscriptVariationAllele.pm#L1700-L1758)
4. **Transcript-level shift via `perform_shift`** — VEP `_return_3prime()` calls `perform_shift()` with `hgvs_reverse` flag
   - [TranscriptVariationAllele.pm#L237-L243](https://github.com/Ensembl/ensembl-variation/blob/release/115/modules/Bio/EnsEMBL/Variation/TranscriptVariationAllele.pm#L237-L243)
   - [TranscriptVariationAllele.pm#L291-L351](https://github.com/Ensembl/ensembl-variation/blob/release/115/modules/Bio/EnsEMBL/Variation/TranscriptVariationAllele.pm#L291-L351)
5. **`stop_loss_extra_aa` ref_len** — VEP `_peptide()` cached peptide excludes terminal `*`; use `trim_end_matches('*').len()`
   - [TranscriptVariationAllele.pm#L2430-L2432](https://github.com/Ensembl/ensembl-variation/blob/release/115/modules/Bio/EnsEMBL/Variation/TranscriptVariationAllele.pm#L2430-L2432)
   - [BaseTranscriptVariation.pm#L1282-L1291](https://github.com/Ensembl/ensembl-variation/blob/release/115/modules/Bio/EnsEMBL/Variation/BaseTranscriptVariation.pm#L1282-L1291)
6. **Port VEP `genomic2pep` for insertion protein position** — map BOTH flanking bases independently via `int((cds_1based + 2) / 3)`; HGVS start = higher position (from seq_region_start mapping), end = lower
   - [TranscriptMapper.pm#L451-L487](https://github.com/Ensembl/ensembl/blob/release/115/modules/Bio/EnsEMBL/TranscriptMapper.pm#L451-L487)
   - [TranscriptVariationAllele.pm#L1680-L1682](https://github.com/Ensembl/ensembl-variation/blob/release/115/modules/Bio/EnsEMBL/Variation/TranscriptVariationAllele.pm#L1680-L1682)
   - [BaseTranscriptVariation.pm#L467-L499](https://github.com/Ensembl/ensembl-variation/blob/release/115/modules/Bio/EnsEMBL/Variation/BaseTranscriptVariation.pm#L467-L499)
7. **LoF biotype UTR suppression** — VEP `_three_prime_utr()` returns undef for `protein_coding_LoF` transcripts
   - [BaseTranscriptVariation.pm#L1106-L1116](https://github.com/Ensembl/ensembl-variation/blob/release/115/modules/Bio/EnsEMBL/Variation/BaseTranscriptVariation.pm#L1106-L1116)
8. **Boundary-spanning deletion HGVSc suppression** — VEP returns undef when variant extends beyond transcript boundaries
   - [TranscriptVariationAllele.pm#L1416](https://github.com/Ensembl/ensembl-variation/blob/release/115/modules/Bio/EnsEMBL/Variation/TranscriptVariationAllele.pm#L1416)
9. **Full cDNA hydration from FASTA** — reconstruct transcript cDNA from exon coordinates + reference FASTA for UTR extension
   - [TranscriptVariationAllele.pm#L2412-L2418](https://github.com/Ensembl/ensembl-variation/blob/release/115/modules/Bio/EnsEMBL/Variation/TranscriptVariationAllele.pm#L2412-L2418)

### Performance optimizations applied

1. **Promoted parquet columns** — `translateable_seq`, `cdna_mapper_segments`, `bam_edit_status`, `has_non_polya_rna_edit`, `spliced_seq`, `flags_str` moved from `raw_object_json` blob to top-level columns, eliminating 7× JSON parsing of ~25KB per transcript (biodatageeks/datafusion-bio-formats#125, #126)
2. **Batched FASTA reads** — read entire transcript span in one FASTA query instead of per-exon reads (~8.7x reduction)
3. **Targeted cDNA hydration** — only hydrate transcripts whose CDS overlaps indel variants or whose stop codon overlaps any variant
4. **SNV shift skip** — skip `build_hgvs_genomic_shift` for SNVs/MNVs (84% of variants)

### Merged (Ensembl+RefSeq) benchmark status

Merged chr1 benchmark with `--hgvs --merged`:
- 72/74 fields at zero mismatches
- HGVSc: ~100 mismatches (RefSeq-specific: edited transcript shifting, cdna_mapper_segments differences)
- HGVSp: ~5 mismatches (RefSeq transcripts missing spliced_seq)
- All non-HGVS fields: zero mismatches

Remaining merged mismatches are ALL on RefSeq transcripts (NM_, NR_, XM_, XR_) and require:
- Porting VEP's transcript-level `_return_3prime()` shifting with correct pre/post sequence extraction
- Resolving `cdna_mapper_segments` differences for 3 specific RefSeq transcripts

### Test coverage

365 unit tests covering:
- `stop_loss_extra_aa`: internal stops (LoF), no new stop, frameshift, `trim_end_matches('*')`
- `perform_shift_ensembl`: forward/reverse, `hgvs_reverse` flag, genomic `seq_strand=1`
- `clip_alleles`: negative strand prefix/suffix coordinate adjustments
- `check_for_peptide_duplication`: current position, fallback offsets, multi-residue
- `three_prime_utr_seq`: LoF suppression, spliced_seq/cdna_seq fallback
- `genomic_to_cds_index`: positive/negative strand, outside CDS
- `coding_segments`: strand-aware ordering
- `translate_protein_from_cds`: stop inclusion, incomplete codons, N bases
- Protein HGVS: deletion, missense, delins, frameshift, synonymous, start_lost

## Phase 5: `--everything` Flag Parity

**Status: blocked on cache columns (biodatageeks/datafusion-bio-formats#127)**

**Goal**: support `--everything` flag producing 80/80 zero mismatches on chr1 non-merged benchmark.

### VEP `--everything` CSQ field order (80 fields, confirmed via Docker VEP 115.2)

```
 1 Allele                 21 FLAGS                  41 miRNA
 2 Consequence            22 VARIANT_CLASS          42 HGVS_OFFSET
 3 IMPACT                 23 SYMBOL_SOURCE          43 AF
 4 SYMBOL                 24 HGNC_ID                44 AFR_AF
 5 Gene                   25 CANONICAL              45 AMR_AF
 6 Feature_type           26 MANE                   46 EAS_AF
 7 Feature                27 MANE_SELECT            47 EUR_AF
 8 BIOTYPE                28 MANE_PLUS_CLINICAL     48 SAS_AF
 9 EXON                   29 TSL                    49-58 gnomADe_*_AF
10 INTRON                 30 APPRIS                 59-69 gnomADg_*_AF
11 HGVSc                  31 CCDS                   70 MAX_AF
12 HGVSp                  32 ENSP                   71 MAX_AF_POPS
13 cDNA_position          33 SWISSPROT              72 CLIN_SIG
14 CDS_position           34 TREMBL                 73 SOMATIC
15 Protein_position       35 UNIPARC                74 PHENO
16 Amino_acids            36 UNIPROT_ISOFORM        75 PUBMED
17 Codons                 37 GENE_PHENO             76 MOTIF_NAME
18 Existing_variation     38 SIFT                   77 MOTIF_POS
19 DISTANCE               39 PolyPhen               78 HIGH_INF_POS
20 STRAND                 40 DOMAINS                79 MOTIF_SCORE_CHANGE
                                                    80 TRANSCRIPTION_FACTORS
```

### Differences from current 74-field schema

| Change | Details |
|--------|---------|
| **New fields (6)** | `APPRIS` (30), `SIFT` (38), `PolyPhen` (39), `DOMAINS` (40), `miRNA` (41), `HGVS_OFFSET` (42) |
| **New field** | `MANE` (26) — generic MANE flag, separate from MANE_SELECT/MANE_PLUS_CLINICAL |
| **Removed** | `SOURCE` — not present in --everything output |
| **Reordered** | `VARIANT_CLASS` moved from after batch1 to position 22 (after FLAGS) |
| **Reordered** | `MOTIF_*` fields moved from positions 24-28 to end (76-80) |
| **Renamed** | gnomAD sub-population fields gained `_AF` suffix (e.g., `gnomADe_AFR` → `gnomADe_AFR_AF`) |

### VEP `--everything` enables (Config.pm L346-374)

`sift=b`, `polyphen=b`, `ccds`, `hgvs`, `symbol`, `numbers`, `domains`, `regulatory`, `canonical`, `protein`, `biotype`, `af`, `af_1kg`, `af_gnomade`, `af_gnomadg`, `max_af`, `pubmed`, `uniprot`, `mane`, `tsl`, `appris`, `variant_class`, `gene_phenotype`, `mirna`

### New cache columns needed (biodatageeks/datafusion-bio-formats#127)

| Column | Table | Type | Source |
|--------|-------|------|--------|
| `appris` | transcript | `Utf8` | Transcript attribute code `appris` |
| `protein_features` | translation | `List<Struct<analysis,hseqname,start,end>>` | `_variation_effect_feature_cache.protein_features` |
| `sift_predictions` | translation | `Map<Utf8, Utf8>` or equivalent | Decoded from binary `ProteinFunctionPredictionMatrix` |
| `polyphen_predictions` | translation | `Map<Utf8, Utf8>` or equivalent | Decoded from binary `ProteinFunctionPredictionMatrix` |

### Implementation steps

**Step 1: Flag parsing** ✅ — Added `everything: bool` to `VepFlags` and `HgvsFlags`. When true, all sub-flags are enabled per Config.pm.

**Step 2: CSQ schema** ✅ — Added 80-field `CSQ_FIELD_NAMES_EVERYTHING` constant in `golden_benchmark.rs`. Updated CSQ assembly in `annotate_provider.rs` with three format paths (cache-hit, transcript engine, fallback) producing 80-field CSQ when `--everything` is active. Field reordering, `SOURCE` removal, `MOTIF_*` at end, and `MANE` generic field placeholder implemented.

**Step 3: Wire new CSQ columns**:
- **HGVS_OFFSET**: ✅ Wired from `HgvsGenomicShift.shift_length` per-transcript (strand-aware). Gated on `hgvsc && tc.hgvsc.is_some()` so offset is only emitted when HGVSc was actually computed.
- **MANE**: ✅ Emits `MANE_Select` when transcript has `mane_select`, `MANE_Plus_Clinical` when it has `mane_plus_clinical`, empty otherwise.
- **APPRIS**: ✅ Loaded from `appris` transcript parquet column. Format: `principal1`→`P1`, `alternative2`→`A2` (matches VEP OutputFactory.pm#L1563-L1570).
- **SIFT**: ✅ Loaded from `sift_predictions` translation parquet column (`List<Struct<position,amino_acid,prediction,score>>`). Lookup by (protein_position, alt_amino_acid) for single AA substitutions only. Format: `prediction(score)` with spaces→underscores (matches VEP `--sift b` mode).
- **PolyPhen**: ✅ Loaded from `polyphen_predictions` translation parquet column (same struct). Format: `prediction(score)` (matches VEP `--polyphen b` mode).
- **DOMAINS**: ✅ Loaded from `protein_features` translation parquet column (`List<Struct<analysis,hseqname,start,end>>`). Overlap check: variant protein_position vs feature [start,end]. Format: `analysis:hseqname` with `[\s;=]`→`_`, joined with `&`. Fixed after biodatageeks/datafusion-bio-formats#128 populated `analysis` column.
- **miRNA**: ✅ Shows 100% in benchmark (no miRNA transcripts in chr1 1000-variant sample). `ncrna_structure` column now available in cache for full structure-aware overlap when miRNA transcripts are encountered.
- **gnomAD sub-pops**: ✅ Fixed: when `flags.everything` is true, all AF columns emit in CSQ (overrides `emit_in_csq: false` for gnomAD sub-populations).

**Step 4: Benchmark** ✅ — Added `--everything` CLI flag to benchmark harness. When active, passes `--everything` to Docker VEP and uses `CSQ_FIELD_NAMES_EVERYTHING` (80 fields) for per-field comparison via `compare_csq_fields_with_names()`.

**Step 5: Parameterized comparison** ✅ — Added `compare_csq_fields_with_names()` that accepts custom field name lists, keeping backward compatibility with `compare_csq_fields()` (defaults to 74-field `CSQ_FIELD_NAMES`).

### Benchmark result: 80/80 zero mismatches ✅

chr1, 1000 variants, non-merged, `--everything` mode: **34,741 CSQ entries compared, all 80 fields at 100.0% match.**

### Known issue: SIFT/PolyPhen memory explosion (#38)

`load_translations()` eagerly loads ALL SIFT/PolyPhen prediction matrices into memory. On chr1 this consumes ~20GB+:
- 22,832 translations × ~11,228 SIFT entries each = ~256M entries (~10GB)
- PolyPhen similar = another ~10GB
- Total: ~26GB observed during benchmark

Fix: load predictions lazily or only for transcripts that overlap input variants. See biodatageeks/datafusion-bio-functions#38.

### VEP traceability

- `--everything` definition: [Config.pm#L346-L374](https://github.com/Ensembl/ensembl-vep/blob/release/115/modules/Bio/EnsEMBL/VEP/Config.pm#L346-L374)
- Flag implications: [Config.pm#L436-L441](https://github.com/Ensembl/ensembl-vep/blob/release/115/modules/Bio/EnsEMBL/VEP/Config.pm#L436-L441)
- APPRIS output: [OutputFactory.pm#L1563-L1570](https://github.com/Ensembl/ensembl-vep/blob/release/115/modules/Bio/EnsEMBL/VEP/OutputFactory.pm#L1563-L1570)
- DOMAINS output: [OutputFactory.pm#L1448-L1466](https://github.com/Ensembl/ensembl-vep/blob/release/115/modules/Bio/EnsEMBL/VEP/OutputFactory.pm#L1448-L1466)
- SIFT/PolyPhen output: [OutputFactory.pm#L1746-L1799](https://github.com/Ensembl/ensembl-vep/blob/release/115/modules/Bio/EnsEMBL/VEP/OutputFactory.pm#L1746-L1799)
- miRNA output: [OutputFactory.pm#L1572-L1612](https://github.com/Ensembl/ensembl-vep/blob/release/115/modules/Bio/EnsEMBL/VEP/OutputFactory.pm#L1572-L1612)
- CSQ field order: [Constants.pm#L66-L138](https://github.com/Ensembl/ensembl-vep/blob/release/115/modules/Bio/EnsEMBL/VEP/Constants.pm#L66-L138)

### Critical files

- `annotate_provider.rs` — flag parsing, CSQ assembly, column loading
- `transcript_consequence.rs` — TranscriptConsequence/Feature structs, miRNA logic
- `golden_benchmark.rs` — CSQ_FIELD_NAMES
- `annotate_vep_golden_bench.rs` — benchmark harness

## Non-Goals

- optimizing runtime before semantic parity is reached
- preserving current workaround behavior for backward compatibility
- adding new consequence features unrelated to Ensembl VEP 115 parity
