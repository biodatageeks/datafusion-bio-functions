## Why

`annotate_vep` currently supports only the `flag_pick_allele_gene` subset of Ensembl VEP's pick-mode surface. To compare reliably against Ensembl VEP runs that use `--pick`, `--pick_allele`, `--per_gene`, `--pick_allele_gene`, or their flagging variants, the Rust VEP implementation needs the same grouping, filtering, and `PICK` emission behavior.

The downstream `vepyr` e2e flow normalizes benchmark input with `bcftools norm -m -both`, so the first validation target can use split biallelic rows while the upstream implementation still models VEP's allele grouping semantics.

## What Changes

- Generalize the existing `PickFlags` / `mark_flag_pick_allele_gene` path into a shared VEP pick-mode engine.
- Parse and support these `options_json` booleans:
  - `pick`
  - `pick_allele`
  - `per_gene`
  - `pick_allele_gene`
  - `flag_pick`
  - `flag_pick_allele`
  - `flag_pick_allele_gene`
- Keep the existing `pick_order` parser and ranking criteria, and use it for every pick/flag-pick mode.
- Apply VEP release 115 mode precedence when multiple pick modes are enabled.
- For filtering modes, emit only selected CSQ/typed annotation entries.
- For flagging modes, retain all entries and add the standalone `PICK` field after `FLAGS`, marking selected entries with `PICK=1`.
- Bypass the cache-hit synthetic CSQ fast path for all pick/flag-pick modes because ranking requires transcript-level consequence candidates.
- Expose the full pick-mode option surface in the VCF sink and downstream `vepyr` wrapper.
- Add e2e comparison profiles for at least the normalized merged-cache pick modes used by the current benchmark workflow.

## Capabilities

### New Capabilities

- `vep-pick-modes`: Ensembl VEP-compatible pick and flag-pick modes for transcript consequence selection, CSQ shaping, typed annotation columns, and downstream VCF output.

### Modified Capabilities

- None. No archived base specs exist yet in `openspec/specs/`; this change introduces a new capability delta.

## Impact

- Affected upstream code:
  - `datafusion/bio-function-vep/src/annotate_provider.rs`
  - `datafusion/bio-function-vep/src/transcript_consequence.rs` if assignment-level allele metadata is needed beyond the current normalized-row path
  - `datafusion/bio-function-vep/src/golden_benchmark.rs`
  - `datafusion/bio-function-vep/src/vcf_sink.rs`
  - `datafusion/bio-function-vep/src/annotate_table_function.rs` tests
- Affected downstream code:
  - `/Users/mwiewior/research/git/vepyr/src/annotate.rs`
  - `/Users/mwiewior/research/git/vepyr/src/lib.rs`
  - `/Users/mwiewior/research/git/vepyr/e2e-testing/scripts/run_annotation_fast.py`
  - `/Users/mwiewior/research/git/vepyr/e2e-testing/scripts/run_annotation.py`
- No new external dependencies.
- No breaking change for existing default `annotate_vep` output when no pick mode is enabled.
