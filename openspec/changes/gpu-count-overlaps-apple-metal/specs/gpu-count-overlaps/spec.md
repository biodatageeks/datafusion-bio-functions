## ADDED Requirements

### Requirement: Optional Apple GPU Count Overlaps Backend
The system SHALL provide an optional Apple GPU backend for `count_overlaps` that preserves the existing table-function contract and returns exactly one count value per right-side input row.

The backend MUST:
- Be gated behind an `apple-gpu` cargo feature
- Compile out cleanly when the feature is disabled
- Be available only on supported macOS Apple GPU targets
- Keep the existing CPU COITree implementation as the default correctness baseline and fallback
- Preserve the current `count_overlaps` output schema
- Preserve right input row order
- Return `Int64` counts matching the CPU path

#### Scenario: CPU remains available without Apple GPU feature
- **WHEN** the crate is built without `apple-gpu`
- **THEN** `count_overlaps` uses the existing CPU implementation
- **AND** no Metal dependency is required

#### Scenario: Unsupported platform fallback
- **WHEN** `bio.count_overlaps_backend = 'auto'`
- **AND** the process is not running on a supported macOS Apple GPU target
- **THEN** `count_overlaps` uses the CPU backend
- **AND** execution succeeds with existing semantics

#### Scenario: Forced GPU unavailable
- **WHEN** `bio.count_overlaps_backend = 'apple_gpu'`
- **AND** Metal device initialization fails before any output batch is emitted
- **THEN** the implementation MAY fall back to CPU with a structured fallback metric
- **AND** the output MUST match CPU results

### Requirement: Endpoint Rank Count Algorithm
The GPU backend SHALL compute count-overlaps using sorted endpoint ranks rather than interval-tree traversal.

For each contig, the backend MUST maintain:
- a sorted array of left interval starts
- a sorted array of left interval ends
- an offset and length for the contig's slice in each array

For each right query interval, the backend MUST compute:
- `started = count(left.start <= query.end)`
- `ended_prior = count(left.end < query.start)`
- `count = started - ended_prior`

#### Scenario: Contained intervals are counted correctly
- **WHEN** left intervals include intervals nested inside other intervals
- **AND** a right interval overlaps both the outer and inner intervals
- **THEN** the GPU backend returns the same count as the CPU COITree backend

#### Scenario: Duplicate intervals are counted independently
- **WHEN** the left table contains duplicate intervals on the same contig
- **AND** a right interval overlaps those duplicates
- **THEN** each duplicate contributes one to the count

#### Scenario: Missing contig returns zero
- **WHEN** a right row references a contig absent from the left index
- **THEN** the GPU backend returns count `0` for that row

### Requirement: Strict and Weak Semantics Match CPU
The GPU backend SHALL implement strict and weak interval semantics exactly as the existing CPU path does.

The backend MUST:
- Apply weak mode without shrinking query intervals
- Apply strict mode by using `query.start + 1` and `query.end - 1`
- Return zero for strict-mode queries that become empty after adjustment
- Match CPU results for point intervals, endpoint-touching intervals, and normal intervals

#### Scenario: Weak endpoint touch counts as overlap
- **WHEN** weak mode is used
- **AND** a left interval endpoint touches a right interval endpoint according to current CPU semantics
- **THEN** the GPU backend returns the same nonzero count as the CPU backend

#### Scenario: Strict endpoint touch does not count
- **WHEN** strict mode is used
- **AND** two intervals only touch at an endpoint
- **THEN** the GPU backend returns the same count as the CPU backend

#### Scenario: Strict point query becomes empty
- **WHEN** strict mode is used
- **AND** the adjusted right interval has `start > end`
- **THEN** the GPU backend returns count `0`

### Requirement: Backend Selection Configuration
The system SHALL expose a session-level backend selection option for `count_overlaps`.

The option MUST support:
- `auto`
- `cpu`
- `apple_gpu`

The default MUST be `auto`.

`auto` mode MUST choose CPU when:
- the Apple GPU feature is not compiled
- the platform is unsupported
- the workload is below benchmark-derived GPU thresholds
- input types or null semantics are unsupported by the GPU backend
- GPU initialization or allocation fails before output starts

#### Scenario: Explicit CPU mode
- **WHEN** `SET bio.count_overlaps_backend = 'cpu'`
- **THEN** `count_overlaps` uses the CPU backend even on Apple GPU hardware

#### Scenario: Auto mode below threshold
- **WHEN** `bio.count_overlaps_backend = 'auto'`
- **AND** the input is below the measured GPU crossover threshold
- **THEN** `count_overlaps` uses the CPU backend

#### Scenario: Auto mode eligible large workload
- **WHEN** `bio.count_overlaps_backend = 'auto'`
- **AND** the platform, feature, input types, and workload size are GPU-eligible
- **THEN** `count_overlaps` uses the Apple GPU backend

### Requirement: Merge Path Rank-Sweep Optimization
The system SHALL support a future Merge Path rank-sweep optimization for large sorted or sort-beneficial `count_overlaps` workloads.

The optimization MUST:
- Compute the same endpoint ranks as the binary-rank GPU kernel
- Use Merge Path partitioning to split sorted endpoint merges into independent threadgroup work
- Preserve original right-row output order through row ID scatter
- Be selected only when benchmarks show it is faster than the binary-rank GPU path

The optimization MUST NOT:
- Be required for initial Apple GPU backend correctness
- Sort every small or unsorted right batch unconditionally
- Materialize overlap pairs
- Maintain an active interval list

#### Scenario: Sorted right partition uses rank sweep
- **WHEN** a right partition is already sorted by the endpoint needed for a rank pass
- **AND** the workload exceeds the Merge Path threshold
- **THEN** the backend MAY use Merge Path rank-sweep for that rank pass
- **AND** the final counts MUST match the binary-rank GPU path

#### Scenario: Random small right batch avoids rank sweep
- **WHEN** a right batch is small or randomly ordered
- **THEN** the backend MUST NOT sort it solely to force Merge Path execution unless benchmarks prove that path is faster for the configured threshold

### Requirement: Correctness and Performance Validation
The implementation SHALL validate the Apple GPU backend against the CPU backend before enabling automatic GPU dispatch.

Validation MUST include:
- deterministic fixture tests
- nested/contained interval tests
- strict and weak semantics tests
- multi-partition right-side tests
- output row-order tests
- ignored Apple GPU integration tests
- benchmark-derived auto-dispatch thresholds

#### Scenario: CPU and GPU agree on deterministic fixtures
- **WHEN** GPU tests are run on supported Apple hardware
- **THEN** GPU `count_overlaps` output matches CPU output exactly for all deterministic fixtures

#### Scenario: Auto dispatch threshold is measured
- **WHEN** the GPU backend is enabled in `auto` mode
- **THEN** the threshold for selecting GPU is based on benchmark results documented in the implementation
- **AND** small workloads continue to use CPU when CPU is faster
