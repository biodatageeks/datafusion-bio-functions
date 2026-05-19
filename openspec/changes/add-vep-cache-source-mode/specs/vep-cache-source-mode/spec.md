## ADDED Requirements

### Requirement: Explicit Cache Source Metadata
The system SHALL require each VEP cache table consumed by `annotate_vep()` to
declare schema metadata key `bio.vep.cache_source_type` with value `ensembl`,
`merged`, or `refseq`.

#### Scenario: Ensembl cache metadata is accepted
- **WHEN** `annotate_vep()` reads a transcript, variation, or context table whose schema metadata contains `bio.vep.cache_source_type = "ensembl"`
- **THEN** the table is accepted as an Ensembl cache table
- **AND** source-specific transcript filtering uses Ensembl semantics

#### Scenario: RefSeq cache metadata is accepted
- **WHEN** `annotate_vep()` reads a transcript table whose schema metadata contains `bio.vep.cache_source_type = "refseq"`
- **THEN** the table is accepted as a RefSeq cache table
- **AND** source-specific transcript filtering uses RefSeq semantics

#### Scenario: Merged cache metadata is accepted
- **WHEN** `annotate_vep()` reads a transcript table whose schema metadata contains `bio.vep.cache_source_type = "merged"`
- **THEN** the table is accepted as a merged cache table
- **AND** source-specific transcript filtering uses merged-cache semantics

#### Scenario: Missing source metadata is rejected
- **WHEN** `annotate_vep()` is called with a cache table that does not expose `bio.vep.cache_source_type`
- **THEN** the system returns a clear error requiring `ensembl`, `merged`, or `refseq`
- **AND** the system does not infer source type from table names, file paths, or cache directory names

#### Scenario: Invalid source metadata is rejected
- **WHEN** `annotate_vep()` is called with `bio.vep.cache_source_type = "foo"`
- **THEN** the system returns a clear error listing the allowed values `ensembl`, `merged`, and `refseq`

### Requirement: No Options JSON Source Selector
The system SHALL NOT accept `merged = true`, `merged = false`, `refseq = true`,
or `refseq = false` in `annotate_vep()` `options_json` as a source-mode
selector.

#### Scenario: Legacy merged boolean is rejected
- **WHEN** a user calls `annotate_vep('vcf', 'cache', 'parquet', '{"merged":true}')`
- **THEN** the system returns a clear error explaining that `merged` is unsupported
- **AND** the error instructs the user to provide cache schema metadata `bio.vep.cache_source_type = "merged"` instead

#### Scenario: Legacy refseq boolean is rejected
- **WHEN** a user calls `annotate_vep('vcf', 'cache', 'parquet', '{"refseq":true}')`
- **THEN** the system returns a clear error explaining that `refseq` is unsupported as an `options_json` source selector
- **AND** the error instructs the user to provide cache schema metadata `bio.vep.cache_source_type = "refseq"` instead

#### Scenario: Source mode is not read from options_json
- **WHEN** a user calls `annotate_vep()` with `options_json` containing `cache_source_type = "refseq"`
- **THEN** the system ignores that value as a source-mode authority
- **AND** the system uses only schema metadata `bio.vep.cache_source_type`
- **AND** missing schema metadata remains an error

### Requirement: Source-Specific Transcript Filtering
The system SHALL filter transcripts according to the explicit cache source mode,
matching Ensembl VEP cache behavior for Ensembl, merged, and RefSeq caches.

#### Scenario: Ensembl mode accepts stable transcripts
- **WHEN** cache source mode is `ensembl`
- **THEN** transcripts with non-empty stable IDs are included after common gencode filters
- **AND** transcript IDs are not required to begin `ENST`

#### Scenario: RefSeq mode accepts canonical RefSeq transcripts
- **WHEN** cache source mode is `refseq`
- **THEN** transcripts with stable IDs beginning `NM_`, `NR_`, `XM_`, and `XR_` are included
- **AND** Ensembl-only transcript IDs beginning `ENST` are excluded

#### Scenario: RefSeq mode keeps mitochondrial RefSeq names
- **WHEN** cache source mode is `refseq`
- **AND** a mitochondrial RefSeq transcript has an Ensembl VEP-accepted ID such as `4540`, `COX3`, or `rna-TRNK`
- **THEN** the transcript is included

#### Scenario: Merged mode accepts Ensembl and RefSeq rows
- **WHEN** cache source mode is `merged`
- **THEN** Ensembl transcript rows beginning `ENST` are included
- **AND** RefSeq rows beginning `NM_`, `NR_`, `XM_`, or `XR_` are included
- **AND** mitochondrial RefSeq rows accepted by Ensembl VEP are included when the row source is RefSeq
- **AND** non-RefSeq-source rows with non-empty stable IDs are included even when their IDs do not begin `ENST`

### Requirement: Source-Specific CSQ SOURCE Output
The system SHALL populate the CSQ `SOURCE` field only for merged cache mode.

#### Scenario: Merged mode populates SOURCE
- **WHEN** cache source mode is `merged`
- **AND** an annotated transcript row has source label `RefSeq` or `Ensembl`
- **THEN** the CSQ `SOURCE` field contains that row source value

#### Scenario: RefSeq mode leaves SOURCE empty
- **WHEN** cache source mode is `refseq`
- **AND** an annotated transcript row has source label `RefSeq`
- **THEN** the CSQ `SOURCE` field is empty

#### Scenario: Ensembl mode leaves SOURCE empty
- **WHEN** cache source mode is `ensembl`
- **AND** an annotated transcript row has source label `Ensembl`
- **THEN** the CSQ `SOURCE` field is empty

### Requirement: RefSeq Hydration Support
The system SHALL apply RefSeq transcript hydration logic to RefSeq cache rows
identified by either normalized row source or RefSeq-compatible stable ID.

#### Scenario: RefSeq stable IDs are hydrated
- **WHEN** a RefSeq transcript row has stable ID beginning `NM_`, `NR_`, `XM_`, or `XR_`
- **AND** the annotation path is configured with `reference_fasta_path`
- **THEN** RefSeq CDS/spliced sequence hydration considers that transcript eligible

#### Scenario: RefSeq mitochondrial IDs are hydrated
- **WHEN** a RefSeq mitochondrial transcript row has stable ID such as `4540`, `COX3`, or `rna-TRNK`
- **AND** the row source normalizes to `RefSeq` or cache source mode is `refseq`
- **THEN** RefSeq hydration considers that transcript eligible

### Requirement: Cache Source Mode Validation Tests
The system SHALL include tests that cover accepted source modes, rejected legacy
options, RefSeq transcript filters, and CSQ `SOURCE` behavior.

#### Scenario: Unit tests cover source mode parsing
- **WHEN** the test suite runs
- **THEN** it verifies accepted metadata values `ensembl`, `merged`, and `refseq`
- **AND** it verifies missing/invalid metadata produces clear errors

#### Scenario: Unit tests cover RefSeq transcript IDs
- **WHEN** the test suite runs
- **THEN** it verifies `NM_`, `NR_`, `XM_`, `XR_`, numeric mitochondrial IDs, and `rna-*` IDs are handled according to source mode

#### Scenario: Integration tests compare RefSeq behavior
- **WHEN** a gated local RefSeq cache fixture is available
- **THEN** the test suite compares representative chr22 and MT annotations against Ensembl VEP output
