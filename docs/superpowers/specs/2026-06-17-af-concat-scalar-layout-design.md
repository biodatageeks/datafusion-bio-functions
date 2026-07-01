# AF concat-scalar layout — design (2026-06-17)

Replaces the fullzip `List<Utf8>` AF bundle with **3 concatenated scalar `Utf8` columns**
encoded miniblock+zstd, to cut read amplification and storage on the variation point-lookup path.

Related: [[af-bundle-stock-lance-fullzip]] (current state), [[af-bundle-fullzip-benchmark]],
[[vep-annotation-bottlenecks]]. Supersedes the fullzip-List approach in
`docs/superpowers/plans/2026-06-17-af-bundle-resume.md`.

## Motivation

The 27 AF columns are read together on every variant point-lookup; storing them separately causes
read amplification (many column regions per row). The fullzip `List` bundle solved column-count
but is large (FSST, no block zstd) and amplifies under the real multi-column take.

Measured (chr1, 735,017 real row_ids = 0.83% of 88.2M rows, 79-batch replay, full 21-col
projection, single-thread):

| | separate-27 (miniblock+zstd) | fullzip-3 List (FSST) | **concat-3 scalar (miniblock+zstd 4 KB)** |
|---|---|---|---|
| AF columns | 27 | 3 | 3 |
| AF on disk | 3.20 GB | 5.74 GB | **1.50 GB** |
| disk bytes read (full take) | 1.99 GB | 1.90 GB | **0.52 GB** |
| IOPS (full take) | 191.9k | 107.8k | **80.3k** |
| wall (full take, warm) | 11.65 s | 4.84 s | 5.35 s |

Ordering for the point-lookup workload: **concat-3 ≫ fullzip-3 ≫ separate-27**. Separate columns
(the pre-bundle base) are worst — 27 column regions ⇒ 192k IOPS, 1.99 GB read. Concat-3 beats
fullzip-3 on size (3.8×) and read (3.6×) at ~equal speed, and beats separate-27 on every axis
(3.8× read, 2.4× IOPS, 2.2× wall, 2.1× disk).

Concat wins on size, bytes-read, IOPS, and cold wall; loses only ~0.5 s warm (decode). The
remaining ~15× amplification vs the 0.83% row fraction is **structural to columnar point lookups**
(chunk granularity, flat across 1–8 KB minichunk) — breaking it needs a row-major point store, out
of scope here (documented as the next-step lever).

## Design

Per AF group (`af_global` 6, `af_gnomade` 10, `af_gnomadg` 11), store **one scalar `Utf8` value per
row**: the members joined by **`|`** (ASCII Unit Separator), empty field for an absent member,
**null** when all members absent. Positional and fixed-width (always N-1 separators), so decode
splits into exactly N members.

- **Encoding**: the concat fields use the shared `lance_field_metadata()` preset (miniblock + zstd
  level 3 + **minichunk 4 KB** + dict settings) — byte-for-byte consistent with the other string
  columns (`variation_name`, `dbsnp_ids`, …). No fullzip, no dict-divisor, no `List` levels → no
  large-chunk bug, stock lance, no fork.
- **Separator**: `|` is verified collision-free — 0 of 613M non-empty AF values contain it
  (`:` is the `allele:freq` delimiter = 100%, `,` = 8%; `|`/`&`/`;`/`/`/tab = 0%). Readable and
  text-inspectable. **Assert at build time** that no AF value contains `|` (done inline during
  `concat_group`, no extra pass); fail loudly if source data ever changes.
- **Large values**: AF members can be up to ~25 KB (STR/repeat-expansion alleles); `split('|')` is
  length-agnostic and Lance stores an over-minichunk value as its own chunk — handled, no special
  casing.
- **Replace** the fullzip path entirely (remove `bundle_group`/`unbundle_group`/`list_item_field`/
  `bundled_list_field`/dict-divisor).

## Components (all in `lance_cache/af_bundle.rs` unless noted)

Write:
- `concat_field(name)` → scalar `Utf8` field with `lance_field_metadata()`.
- `concat_group(members) -> StringArray` (prototyped) — join with `|`, null when all-absent.
- `bundle_schema` / `bundle_af_columns` → emit the 3 concat columns (drop the env gate; default).
- build-time `|` collision assertion.

Read:
- `split_group(col: &StringArray, width) -> Vec<StringArray>` — split each value on `|` into N
  positional members (empty → `""`/absent), null parent → all-absent.
- `unbundle_af_columns` → split the 3 concat columns into the 27 logical `Utf8` AF columns.
- `bundle_projection` (`variation_runtime.rs`) → map the 27 logical AF names to the 3 concat column
  names (same names, scalar type).

## Phases

1. **Code** — write-side default concat; read-side `split_group` + `unbundle_af_columns`;
   projection; build assertion; update tests (`lance_lookup_round_trips_bundled_af_columns`,
   `bundle_unbundle_round_trips_values`) + new `concat_group`↔`split_group` round-trip.
2. **Data** — rebuild `variation.lance/chr1` with concat (context datasets unchanged).
3. **Validate e2e** — rebuild vepyr (path dep), run `run_annotation_fast.py chr1 --backend lance`
   **without `--skip-compare`** → gate: **0 CSQ/AF mismatches** vs VEP reference; confirm size +
   lookup timing.
4. **Cleanup** — remove fullzip-List code + diagnostic scaffolding (rowid dump in
   `variation_runtime.rs`, `examples/frag_rows.rs`, `examples/af_row_size.rs`,
   `examples/af_sep_check.rs`); revert sandbox
   `Cargo.toml` patch removal is kept (stock). Commit.

## Gates
- 0 CSQ/AF mismatches (parity, the correctness gate for split-on-read).
- AF on disk ≈ 1.5 GB (chr1); whole dataset ≈ 4.3 GB.
- e2e completes on stock lance, single-thread, no panic.

## Out of scope (follow-ups)
- **Non-AF minichunk tuning** — non-AF columns now dominate reads (304 MB of 524 MB); bumping their
  minichunk is the likely next win (one-line `lance_field_metadata` change + the same bench).
- **Row-major point store** (cold-tier point-store) — the only thing that breaks the ~15× columnar
  point-lookup floor.
- Genome-wide variation rebuild (only chr1 in this plan).
