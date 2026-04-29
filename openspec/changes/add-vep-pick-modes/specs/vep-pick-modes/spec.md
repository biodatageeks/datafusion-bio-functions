## ADDED Requirements

### Requirement: VEP Pick Mode Option Parsing
The system SHALL parse Ensembl VEP-compatible pick and flag-pick options from `annotate_vep` `options_json`.

The supported boolean options MUST include:
- `pick`
- `pick_allele`
- `per_gene`
- `pick_allele_gene`
- `flag_pick`
- `flag_pick_allele`
- `flag_pick_allele_gene`

The system MUST continue to parse `pick_order` as a comma-separated list of VEP ranking criteria.

#### Scenario: No pick mode enabled
- **WHEN** no pick or flag-pick option is set
- **THEN** the system uses normal unfiltered consequence output
- **AND** the CSQ schema does not include `PICK`

#### Scenario: Single pick mode enabled
- **WHEN** `options_json` contains `{"pick": true}`
- **THEN** the system selects `PickMode::Pick`
- **AND** it uses the configured or default `pick_order`

#### Scenario: Multiple pick modes enabled
- **WHEN** more than one pick or flag-pick option is set
- **THEN** the system applies Ensembl VEP release 115 precedence: `pick`, `pick_allele`, `per_gene`, `pick_allele_gene`, `flag_pick`, `flag_pick_allele`, `flag_pick_allele_gene`
- **AND** only the highest-precedence enabled mode controls output shaping

#### Scenario: Invalid pick order criterion
- **WHEN** `pick_order` contains an unsupported criterion
- **THEN** `annotate_vep` fails with an execution error naming the invalid criterion

#### Scenario: Empty pick order
- **WHEN** `pick_order` parses to no criteria
- **THEN** `annotate_vep` fails with an execution error stating that `pick_order` must contain at least one criterion

### Requirement: Shared VEP Ranking Engine
The system SHALL use one shared ranking engine for every pick and flag-pick mode.

The ranking engine MUST:
- use the default order `mane_select,mane_plus_clinical,canonical,appris,tsl,biotype,ccds,rank,length,ensembl,refseq`
- allow callers to override the order via `pick_order`
- apply criteria progressively, retaining only tied best candidates after each criterion
- select the first remaining candidate after all criteria if ties remain
- score lower values as better, matching Ensembl VEP release 115
- use VEP-compatible transcript length scoring where longer transcripts rank better
- rank `rank` by the most severe Sequence Ontology consequence on the assignment

#### Scenario: Default ranking selects one candidate
- **WHEN** multiple transcript consequences are candidates for one pick group
- **AND** default `pick_order` yields a unique best candidate
- **THEN** the system selects that candidate

#### Scenario: Custom pick order changes winner
- **WHEN** `pick_order` changes the leading criterion from `mane_select` to `rank`
- **THEN** the system ranks candidates by consequence severity before later criteria
- **AND** the selected candidate can differ from the default-order winner

#### Scenario: Ties fall through to later criteria
- **WHEN** two candidates tie on the first configured criterion
- **THEN** the system keeps both candidates
- **AND** evaluates the next configured criterion only on that tied subset

#### Scenario: Ties remain after all criteria
- **WHEN** candidates tie on every configured criterion
- **THEN** the system selects the first candidate in deterministic feature-order input order

### Requirement: Filtering Pick Modes
The system SHALL support VEP-compatible filtering pick modes that reduce emitted consequence entries.

Filtering modes MUST NOT add a `PICK` field to CSQ by default.

#### Scenario: Pick one consequence per variant
- **WHEN** `pick` is enabled
- **THEN** the system emits exactly one selected consequence entry for the input variant
- **AND** all non-selected consequence entries are omitted from CSQ and typed list columns

#### Scenario: Pick one consequence per allele
- **WHEN** `pick_allele` is enabled
- **THEN** the system groups candidate consequences by CSQ allele value
- **AND** emits exactly one selected consequence entry per allele group

#### Scenario: Pick one consequence per gene
- **WHEN** `per_gene` is enabled
- **THEN** the system groups transcript consequence candidates by gene stable ID
- **AND** emits one selected transcript consequence per gene group
- **AND** retains non-transcript consequence entries that Ensembl VEP retains outside transcript gene groups

#### Scenario: Pick one consequence per allele and gene
- **WHEN** `pick_allele_gene` is enabled
- **THEN** the system groups candidates by CSQ allele value
- **AND** within each allele group, groups transcript candidates by gene stable ID
- **AND** emits one selected transcript consequence per allele-plus-gene group
- **AND** retains non-transcript consequence entries for each allele group according to VEP per-gene behavior

#### Scenario: Normalized biallelic e2e input
- **WHEN** input has been normalized with `bcftools norm -m -both`
- **THEN** each VCF row has one ALT allele
- **AND** `pick` and `pick_allele` select the same number of entries for that row
- **AND** `per_gene` and `pick_allele_gene` select the same transcript groups for that row

### Requirement: Flag Pick Modes
The system SHALL support VEP-compatible flag-pick modes that retain all consequence entries and mark selected entries with `PICK=1`.

Flag-pick modes MUST include a standalone `PICK` CSQ field immediately after `FLAGS`, and typed annotation output MUST expose a corresponding `PICK` list column aligned with other transcript-level list columns.

#### Scenario: Flag one consequence per variant
- **WHEN** `flag_pick` is enabled
- **THEN** the system retains every consequence entry
- **AND** marks exactly one selected entry with `PICK=1`
- **AND** leaves unselected entries with an empty CSQ `PICK` field and null typed `PICK` value

#### Scenario: Flag one consequence per allele
- **WHEN** `flag_pick_allele` is enabled
- **THEN** the system retains every consequence entry
- **AND** marks exactly one selected entry per CSQ allele group with `PICK=1`

#### Scenario: Flag one consequence per allele and gene
- **WHEN** `flag_pick_allele_gene` is enabled
- **THEN** the system retains every consequence entry
- **AND** marks one selected transcript consequence per allele-plus-gene group with `PICK=1`
- **AND** marks retained non-transcript entries according to VEP per-gene flagging behavior

#### Scenario: CSQ header includes PICK for flag modes
- **WHEN** any flag-pick mode is enabled
- **THEN** the CSQ header field list includes `PICK`
- **AND** `PICK` appears immediately after `FLAGS`

### Requirement: Pick Modes Use Transcript Consequence Evaluation
The system SHALL evaluate transcript consequences for all pick and flag-pick modes instead of using synthetic variation-cache CSQ rows.

#### Scenario: Cache hit with pick mode
- **WHEN** an input variant has a variation-cache consequence hit
- **AND** any pick or flag-pick mode is enabled
- **THEN** the system bypasses the synthetic cache-hit CSQ fast path
- **AND** computes transcript/regulatory/motif assignments through the transcript consequence engine

#### Scenario: Cache hit without pick mode
- **WHEN** an input variant has a variation-cache consequence hit
- **AND** no pick or flag-pick mode is enabled
- **THEN** the system MAY use the existing synthetic cache-hit CSQ fast path

### Requirement: CSQ and Typed Column Alignment
The system SHALL keep CSQ entries and typed list annotation columns aligned after applying pick or flag-pick modes.

#### Scenario: Filtering mode alignment
- **WHEN** a filtering pick mode removes non-selected assignments
- **THEN** the Nth emitted CSQ entry corresponds to the Nth element in typed list columns such as `Feature`, `Consequence`, `BIOTYPE`, and `PICK` when present

#### Scenario: Flag mode alignment
- **WHEN** a flag-pick mode retains all assignments
- **THEN** typed list column order matches CSQ entry order
- **AND** the typed `PICK` list marks the same entries that have CSQ `PICK=1`

#### Scenario: Deterministic output order
- **WHEN** selected assignments are emitted
- **THEN** the system orders entries using the existing VEP-compatible CSQ sort order: transcript entries first, regulatory entries next, motif entries next, and deterministic feature ID ordering within feature type

### Requirement: Downstream VCF and Python API Exposure
The system SHALL expose the full pick and flag-pick option surface through the VCF sink and downstream `vepyr` wrapper.

#### Scenario: Rust VCF sink emits pick options
- **WHEN** an `AnnotateVcfConfig` enables any pick or flag-pick option
- **THEN** `AnnotateVcfConfig::to_options_json` serializes the corresponding boolean option
- **AND** serializes `pick_order` when provided

#### Scenario: Python annotate passes pick options
- **WHEN** `vepyr.annotate` is called with a pick or flag-pick keyword argument
- **THEN** the downstream wrapper passes the option through to the Rust VCF annotation config

#### Scenario: Existing downstream behavior remains compatible
- **WHEN** existing code calls `vepyr.annotate(..., flag_pick_allele_gene=True, pick_order=...)`
- **THEN** the call continues to work
- **AND** produces the same PICK semantics as before unless corrected for VEP parity by this change

### Requirement: Normalized E2E Pick-Mode Comparison
The system SHALL provide e2e comparison coverage for pick and flag-pick modes using the existing normalized VCF workflow.

#### Scenario: Fast runner normalizes input
- **WHEN** `e2e-testing/scripts/run_annotation_fast.py` is run without `--no-normalize`
- **THEN** it normalizes the input VCF with `bcftools norm -m -both`
- **AND** annotates the chromosome-extracted normalized VCF

#### Scenario: Flag-pick comparison reports PICK parity
- **WHEN** a fast e2e profile enables a flag-pick mode
- **THEN** the generated comparison report includes shared CSQ field comparison for `PICK`
- **AND** `PICK` mismatches are counted in `field_mismatch_counts` when present

#### Scenario: Filtering pick comparison reports selected-entry parity
- **WHEN** a fast e2e profile enables a filtering pick mode
- **THEN** the generated comparison report compares selected CSQ entries against the matching Ensembl VEP output
- **AND** reports CSQ entry-count mismatches when selected entry counts differ
