## 1. Upstream Picker Configuration

- [x] 1.1 Replace `PickFlags.flag_pick_allele_gene` with `PickMode` plus `pick_order`
- [x] 1.2 Parse `pick`, `pick_allele`, `per_gene`, `pick_allele_gene`, `flag_pick`, `flag_pick_allele`, and `flag_pick_allele_gene` from `options_json`
- [x] 1.3 Preserve existing `pick_order` parsing and error messages for invalid or empty criteria
- [x] 1.4 Implement VEP release 115 precedence when multiple pick modes are enabled
- [x] 1.5 Add unit tests for mode parsing, default mode, and multi-flag precedence

## 2. Upstream Selection Engine

- [x] 2.1 Replace `mark_flag_pick_allele_gene` with a generic `apply_pick_mode` / selection helper
- [x] 2.2 Implement all-assignment grouping for `pick` and `flag_pick`
- [x] 2.3 Implement allele grouping for `pick_allele` and `flag_pick_allele`
- [x] 2.4 Implement per-gene grouping for `per_gene`
- [x] 2.5 Implement allele-plus-gene grouping for `pick_allele_gene` and `flag_pick_allele_gene`
- [x] 2.6 Preserve VEP behavior for non-transcript assignments in per-gene modes
- [x] 2.7 Apply filtering before CSQ sorting and typed-column population
- [x] 2.8 Add unit tests for every mode using small synthetic assignments

## 3. CSQ Schema and Annotation Output

- [x] 3.1 Compute `include_pick_output` from flag modes only
- [x] 3.2 Ensure `PICK` is inserted after `FLAGS` in CSQ headers for `flag_pick`, `flag_pick_allele`, and `flag_pick_allele_gene`
- [x] 3.3 Ensure filtering modes do not emit a `PICK` CSQ field by default
- [x] 3.4 Generalize cache-hit fast-path bypass from `flag_pick_allele_gene` to all non-`None` pick modes
- [x] 3.5 Add `annotate_table_function` tests for CSQ entry counts and selected features per mode
- [x] 3.6 Add tests that typed list columns stay aligned with CSQ entries after filtering and flagging
- [x] 3.7 Decide and document `most_severe_consequence` behavior for filter modes

## 4. VCF Sink and Benchmark Helpers

- [x] 4.1 Add all pick/flag-pick booleans to `AnnotateVcfConfig`
- [x] 4.2 Serialize all enabled pick/flag-pick booleans in `AnnotateVcfConfig::to_options_json`
- [x] 4.3 Update `vcf_sink` tests for CSQ header layout in flag modes
- [x] 4.4 Update golden comparison helpers if they need to select CSQ field names by pick mode instead of a single boolean

## 5. Downstream vepyr Integration

- [x] 5.1 Update `/Users/mwiewior/research/git/vepyr` dependency pin after upstream commit lands
- [x] 5.2 Add `pick`, `pick_allele`, `per_gene`, `pick_allele_gene`, `flag_pick`, and `flag_pick_allele` parameters to `vepyr.annotate`
- [x] 5.3 Keep existing `flag_pick_allele_gene` and `pick_order` downstream API behavior
- [x] 5.4 Pass all new parameters through to `AnnotateVcfConfig`
- [x] 5.5 Add or update downstream tests/checks for Python API parameter forwarding

## 6. E2E Validation

- [x] 6.1 Add `run_annotation_fast.py` cache profiles for at least `merged_pick`, `merged_pick_allele`, `merged_per_gene`, and `merged_pick_allele_gene`
- [x] 6.2 Document that the fast e2e runner normalizes input with `bcftools norm -m -both`
- [ ] 6.3 Run `chr22` normalized comparisons against Ensembl VEP for flag and filtering modes
- [ ] 6.4 Verify `PICK` field parity for flag modes
- [ ] 6.5 Verify selected CSQ entry parity for filtering modes
- [ ] 6.6 Capture report JSON paths and summary metrics in the implementation PR

## 7. Final Verification

- [x] 7.1 Run `cargo fmt --all -- --check`
- [x] 7.2 Run `cargo test -p datafusion-bio-function-vep --lib`
- [x] 7.3 Run `cargo check -p datafusion-bio-function-vep --all-targets`
- [x] 7.4 Run downstream `cargo check` in `/Users/mwiewior/research/git/vepyr`
- [ ] 7.5 Push upstream and downstream commits
