## Context

Ensembl VEP release 115 implements all pick and flag-pick modes in `Bio::EnsEMBL::VEP::OutputFactory`. The modes share one ranking function (`pick_worst_VariationFeatureOverlapAllele`) and differ mainly by grouping policy and whether selected entries replace or merely mark the full consequence set.

The current Rust implementation already has the difficult ranking pieces:
- `PickCriterion` and `pick_order` parsing
- VEP-style default order: `mane_select,mane_plus_clinical,canonical,appris,tsl,biotype,ccds,rank,length,ensembl,refseq`
- progressive tie filtering in `pick_worst_assignment`
- transcript length, APPRIS, TSL, source, rank, MANE, canonical, CCDS, and biotype scoring
- `flag_pick_allele_gene` marking with a standalone CSQ `PICK` field

The missing piece is the full option matrix: all grouping modes, filtering modes, flagging modes, schema/header behavior, and downstream API exposure.

The current `vepyr` fast e2e runner normalizes input with `bcftools norm -m -both` before annotation. That means benchmark rows are split into one ALT allele per row, so `pick` and `pick_allele` behave the same in the current e2e path, as do `per_gene` and `pick_allele_gene`. The upstream implementation should still use VEP-compatible allele grouping so it does not encode the normalized-input assumption into the core picker.

## Goals / Non-Goals

**Goals:**
- Match Ensembl VEP release 115 pick/flag-pick behavior for small variants in `annotate_vep`.
- Reuse the existing ranking criteria implementation instead of introducing a second picker.
- Preserve default output when no pick mode is enabled.
- Keep CSQ and typed list columns aligned after filtering or flagging.
- Expose the full option surface through `AnnotateVcfConfig` and `vepyr.annotate`.
- Validate against normalized `bcftools norm -m -both` e2e inputs first.

**Non-Goals:**
- Implement `--most_severe` or `--summary`; those are separate output-shaping modes.
- Change consequence classification or HGVS generation.
- Guarantee byte-for-byte VCF parity with Ensembl VEP; comparison remains CSQ field/entry parity.
- Add unsplit multi-allelic e2e parity as a blocking requirement for the first implementation. The core grouping design must not prevent it later.
- Add new cache formats or dependencies.

## Decisions

### Decision 1: Represent pick behavior as a single `PickMode`

Introduce:

```rust
enum PickMode {
    None,
    Pick,
    PickAllele,
    PerGene,
    PickAlleleGene,
    FlagPick,
    FlagPickAllele,
    FlagPickAlleleGene,
}
```

`PickFlags` becomes a configuration container with `mode: PickMode` and `pick_order: Vec<PickCriterion>`.

When multiple mode booleans are present, use Ensembl VEP release 115 precedence from `OutputFactory.pm`:

1. `pick`
2. `pick_allele`
3. `per_gene`
4. `pick_allele_gene`
5. `flag_pick`
6. `flag_pick_allele`
7. `flag_pick_allele_gene`

Rationale: VEP's implementation uses an `if`/`elsif` chain rather than rejecting combinations. Mirroring that precedence avoids surprising comparison differences for users who pass multiple flags.

### Decision 2: Use one selector for all modes

Replace `mark_flag_pick_allele_gene` with a generic function that computes selected assignment indices:

```rust
fn apply_pick_mode(
    assignments: Vec<TranscriptConsequence>,
    ctx: &PreparedContext<'_>,
    pick_flags: &PickFlags,
    row_allele: &str,
) -> Vec<TranscriptConsequence>
```

The function should:
- build groups according to `PickMode`
- use the existing `pick_worst_assignment` for each group
- for filter modes, return only selected entries plus VEP-retained non-transcript entries for per-gene modes
- for flag modes, set `picked = true` on selected entries and return all entries

The grouping policies are:
- `Pick` / `FlagPick`: all assignments for the row
- `PickAllele` / `FlagPickAllele`: assignments grouped by CSQ allele
- `PerGene`: transcript assignments grouped by gene, with non-transcript assignments retained
- `PickAlleleGene` / `FlagPickAlleleGene`: first group by CSQ allele, then apply per-gene selection within each allele group

For current normalized e2e inputs, all assignments in a row share the same `row_allele`. If future unsplit multi-allelic parsing produces multiple assignment alleles in one row, the grouping key must come from the assignment's CSQ allele, not from the row-level fallback.

### Decision 3: Filter before CSQ serialization and typed-column population

Filtering modes must modify `row_assignments` before `sorted_indices` is built. That keeps:
- CSQ entries
- typed list columns (`Feature`, `Consequence`, `PICK`, etc.)
- feature ordering

aligned by construction.

Flag modes keep all assignments but set `picked = true` before serialization and typed-column population.

### Decision 4: `PICK` field is emitted only for flag modes

Ensembl VEP adds `PICK` as an output field for `flag_pick`, `flag_pick_allele`, and `flag_pick_allele_gene`. Plain filtering modes (`pick`, `pick_allele`, `per_gene`, `pick_allele_gene`) reduce the consequence set and do not require a `PICK` field.

Therefore:
- `include_pick_output` should be true only for `PickMode::Flag*`.
- `golden_benchmark::csq_field_names_for_mode(..., include_pick)` remains the source of CSQ field layout.
- The VCF CSQ header must include `|FLAGS|PICK|...` only for flag modes.

### Decision 5: All pick modes bypass the synthetic cache-hit path

The cache-hit fast path can synthesize a CSQ row from variation-cache consequence strings, but pick modes require ranking transcript/regulatory/motif assignment candidates. All modes except `None` must force transcript-engine evaluation when CSQ or typed annotation columns are requested.

This generalizes the current `flag_pick_allele_gene` bypass.

### Decision 6: Downstream exposure is a thin pass-through

Add booleans to `AnnotateVcfConfig` and `vepyr.annotate` rather than adding a single stringly typed mode. This matches Ensembl VEP CLI names and keeps e2e profiles readable:

```python
vepyr.annotate(..., pick=True)
vepyr.annotate(..., flag_pick_allele_gene=True, pick_order="...")
```

The Rust side remains responsible for VEP precedence if multiple booleans are true.

### Decision 7: E2E validation starts from normalized input

The e2e scripts already call `bcftools norm -m -both` unless `--no-normalize` is passed. Initial validation should explicitly record that:
- `pick` and `pick_allele` are equivalent on normalized input
- `per_gene` and `pick_allele_gene` are equivalent on normalized input
- full unsplit multi-ALT parity can be added later without blocking this capability

## Risks / Trade-offs

- **Risk: multiple pick flags are passed together** -> Mitigation: mirror VEP release 115 precedence and add tests documenting the precedence.
- **Risk: non-transcript entries in per-gene modes differ from VEP** -> Mitigation: port Ensembl's `pick_VariationFeatureOverlapAllele_per_gene` behavior: non-transcript assignments are retained, and flag per-gene modes mark retained non-transcript entries.
- **Risk: filtered typed columns drift from CSQ ordering** -> Mitigation: apply filtering before `sorted_indices` and reuse the same `sorted_indices` for CSQ and typed columns.
- **Risk: current rows are biallelic but future parser behavior changes** -> Mitigation: design grouping around CSQ allele value, using row-level allele only as the current fallback.
- **Risk: e2e comparison hides unsplit multi-allelic bugs** -> Mitigation: document normalized e2e scope and add a small unit/integration fixture for allele grouping even before full unsplit VCF e2e.

## Migration Plan

1. Implement the upstream generalized picker behind existing option parsing.
2. Preserve `flag_pick_allele_gene` behavior and tests as a compatibility baseline.
3. Add new unit tests for each pick mode and precedence.
4. Add `AnnotateVcfConfig` fields and update `vepyr` dependency pin.
5. Add downstream wrapper parameters and e2e cache profiles.
6. Run local upstream tests, downstream `cargo check`, and at least one normalized `chr22` e2e profile against Ensembl VEP.

Rollback is straightforward: disable the new options in downstream profiles or revert the upstream change; default `annotate_vep` behavior is unchanged when no pick mode is set.

## Open Questions

- Should `most_severe_consequence` in `annotate_vep` remain based on all computed assignments for filter modes, or should it be recomputed from emitted assignments? Ensembl VEP VCF output has no equivalent scalar column, so this is a local API decision.
- Should we add full unsplit multi-allelic VCF e2e parity now, or defer until `datafusion-bio-format-vcf` row semantics for multi-ALT records are explicitly specified?
