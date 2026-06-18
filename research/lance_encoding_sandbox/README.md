# Lance Encoding Sandbox

Research-only harness for building and benchmarking chr-local Lance variation-cache layouts.

```bash
RUSTFLAGS="-C target-cpu=native" cargo run --release --manifest-path research/lance_encoding_sandbox/crates/lance_sandbox/Cargo.toml -- \
  run --config research/lance_encoding_sandbox/configs/current_v21_zstd3.toml
```

The sandbox writes generated datasets and reports under the config's `dataset.output_root`.
Benchmark examples use `RUSTFLAGS="-C target-cpu=native"` so the Rust CLI is compiled for the local CPU before measuring Lance paths.

## Commands

- `build --config <file>`: build Lance dataset and indexes.
- `inspect --config <file>`: inspect the generated dataset's actual schema, encodings, and sizes.
- `inspect-path --dataset-path <path> [--lance-version 2.1]`: inspect any Lance dataset path and emit JSON to stdout. The PyLance notebook uses this command so field metadata, compressed bytes, page counts, and observed encodings come from Rust Lance file metadata instead of PyLance schema guesses.
- `bench --config <file> [--positions-file <file>]`: run warm full scan, cold lookup-batch sweep, Lance row-id scan + `take_rows`, external position-to-row-id sidecar + `take_rows`, and any configured heed payload sidecar probes; if the configured positions file is missing, derive it first. An explicit `--positions-file` must already exist. Progress is printed while each phase runs, `reports/benchmark.partial.json` is refreshed after each completed phase, and the final result is printed with `termimad` table formatting before `reports/benchmark.json` is written.
- `run --config <file> [--positions-file <file>]`: build, inspect, benchmark, and write a summary. The optional positions override has the same behavior as `bench`.
- `extract-keys --config <file>`: derive chr-local cold `position UInt32` lookup keys from the configured VCF and external indexes using the same extended VEPyr probe-start expansion as the cold-path side benchmark.
- `sample-cold --config <file> [--limits 10000,100000] [--output-dir <dir>]`: scan the built cold tier, deduplicate positions, and write evenly sampled lookup files such as `chr1_cold_sample_10k_positions_u32.txt` and `chr1_cold_sample_100k_positions_u32.txt`.

Example benchmark inputs:

```bash
RUSTFLAGS="-C target-cpu=native" cargo run --release --manifest-path research/lance_encoding_sandbox/crates/lance_sandbox/Cargo.toml -- \
  extract-keys --config research/lance_encoding_sandbox/configs/current_v21_zstd3.toml

RUSTFLAGS="-C target-cpu=native" cargo run --release --manifest-path research/lance_encoding_sandbox/crates/lance_sandbox/Cargo.toml -- \
  sample-cold --config research/lance_encoding_sandbox/configs/current_v21_zstd3.toml --limits 10000,100000

RUSTFLAGS="-C target-cpu=native" cargo run --release --manifest-path research/lance_encoding_sandbox/crates/lance_sandbox/Cargo.toml -- \
  bench --config research/lance_encoding_sandbox/configs/current_v21_zstd3.toml \
  --positions-file research/lance_encoding_sandbox/inputs/chr1_cold_sample_10k_positions_u32.txt

RUSTFLAGS="-C target-cpu=native" cargo run --release --manifest-path research/lance_encoding_sandbox/crates/lance_sandbox/Cargo.toml -- \
  bench --config research/lance_encoding_sandbox/configs/current_v21_zstd3.toml \
  --positions-file research/lance_encoding_sandbox/inputs/chr1_cold_sample_100k_positions_u32.txt
```

Targeted Lance 2.2 struct-packing experiment:

```bash
RUSTFLAGS="-C target-cpu=native" cargo run --release --manifest-path research/lance_encoding_sandbox/crates/lance_sandbox/Cargo.toml -- \
  run --config research/lance_encoding_sandbox/configs/packed_targeted_v22.toml \
  --positions-file research/lance_encoding_sandbox/inputs/chr1_cold_sample_10k_positions_u32.txt
```

The benchmark resolves the logical VEPyr `--everything` workload against both flat schemas and packed schemas. If a required logical field is not top-level but exists as a child of a top-level struct, the physical projection reads that struct parent once.

## PyLance Notebook

The PyLance notebook builds Lance datasets from the variation Parquet files, applies per-field Arrow metadata, creates scalar indexes, and renders inspect output with `rich` tables. Field metadata and observed encodings are collected by invoking the Rust sandbox CLI's `inspect-path` command against the generated dataset.

```bash
uv venv research/lance_encoding_sandbox/.venv-pylance
uv pip install --python research/lance_encoding_sandbox/.venv-pylance/bin/python \
  -r research/lance_encoding_sandbox/requirements-pylance.txt
research/lance_encoding_sandbox/.venv-pylance/bin/python -m ipykernel install \
  --user --name lance-pylance-sandbox --display-name "Lance PyLance Sandbox"
```

Open `research/lance_encoding_sandbox/notebooks/pylance_variation_builder.ipynb` with the `Lance PyLance Sandbox` kernel. The notebook defaults to `RUN_BUILD = False` and `RUN_INDEXES = False` so opening or executing it does not accidentally build the full chr1 dataset. Set those flags to `True` after editing `FIELD_METADATA`; use `ROW_LIMIT_PER_TIER = 10000` for a smoke build before a full run. Set `STRUCT_PACKING_ENABLED = True` to build the targeted Lance 2.2 packed-struct schema used by `packed_targeted_v22.toml`.

The notebook uses `PARQUET_BATCH_SIZE` for the input PyArrow batches. PyLance passes these batches into the Lance writer, so this stays modest for nested/list fields while avoiding very small Python batches. It also sets `LANCE_MINIBLOCK_MAX_VALUES = "4096"` in `os.environ` immediately before writing the dataset; set that notebook variable to another string value for sweeps, or `None` to remove the env var and use Lance's built-in default.

By default, the notebook writes warm and cold tiers in separate Lance transactions (`WRITE_TIERS_SEPARATELY = True`). This keeps tier fragments physically separated and avoids a Lance 2.1 writer failure seen when one combined stream crosses from warm `variant_keys` list pages into cold all-null `variant_keys` pages. Reducing `max_rows_per_group` or the PyArrow input batch size to `1024` does not address that boundary condition.

## Config Parameters

For Lance encoding metadata, "Lance default" means the sandbox does not attach the corresponding `lance-encoding:*` field metadata and Lance chooses its built-in behavior for the field type and file version.
Unknown TOML keys are rejected, so removed parameters fail fast instead of being silently ignored.

### Dataset Parameters

`name`

Names the experiment directory under `dataset.output_root`. Lance default: none; this is a sandbox output naming field.

`cache_root`

Points at the VEPyr merged cache root containing `variation` Parquet input files. Lance default: none; this is a sandbox source path.

`chrom`

Selects the chromosome to build and benchmark. Lance default: none; this is a sandbox source filter. The current benchmark inputs are chr1-oriented.

`output_root`

Root directory for generated Lance datasets and reports. Lance default: none; this is a sandbox output path.

`lance_version`

Selects Lance file version `2.1` or `2.2`. Lance `WriteParams.data_storage_version` defaults to `None`, which resolves to Lance `2.1`; the sandbox requires this field to keep every run explicit.

`batch_size = 65536`

Controls how many rows are read from the source Parquet files per Arrow record batch while building the sandbox dataset. This is a sandbox/Parquet reader setting, not a Lance `WriteParams` setting, so there is no Lance write default for it. It affects build memory/throughput, not the benchmark scanner output batch size.

`overwrite = true`

Allows `build` or `run` to delete and rebuild the existing output dataset at `<output_root>/<dataset.name>/<chrom>.lance`. Lance `WriteParams` defaults to `WriteMode::Create`, which fails if a dataset already exists. The sandbox config default is `false`; the provided configs set it to `true` for repeatable experiments.

### Lookup Key Parameters

`mode = "position"`

Uses a chromosome-local integer lookup key instead of the global encoded `position_key`. This is sandbox logic, not a Lance parameter, so Lance has no default. The preferred sandbox mode is `position` because each dataset is already chromosome-specific.

`column = "position"`

Names the generated lookup column. Lance has no default column name. In `position` mode, the build step adds this column as `UInt32` and drops `position_key` from the optimized layout.

`type = "uint32"`

Declares the lookup key type. Lance has no default for sandbox lookup keys. `mode = "position"` requires `uint32`; compatibility mode `position_key` requires `uint64`.

Position source

The source is no longer configurable. In `position` mode the builder always derives `position` from the Arrow `start` column and fails if `start` is missing, null, negative, or outside `UInt32`.

### Index Parameters

`position_btree = true`

Creates a Lance BTree scalar index on the configured lookup key column. Lance default: no secondary index is created unless requested.

`tier_bitmap = true`

Creates a Lance bitmap scalar index on `tier`. Lance default: no secondary index is created unless requested.

### Benchmark Parameters

`positions_file`

Text file containing chr-local `UInt32` lookup positions, one per line. Lance default: none; this is a sandbox benchmark input. The provided chr1 e2e input is produced by `extract-keys`; it expands each VCF row into VEPyr probe starts, removes warm-covered positions, checks the external position index, and then applies the external variant bloom admission check.

`vcf_path`

VCF used by `extract-keys` to derive lookup positions. Lance default: none; this is only used by the key extraction helper.

`lookup_batch_sizes = [238, 512, 1024, 2000, 4096, 8192, 16384]`

Defines the cold lookup sweep. This is sandbox benchmark logic, not a Lance scanner default. For each value, the benchmark sorts/deduplicates the positions when enabled, splits them into chunks of this many keys, and runs one Lance scan per chunk using `tier = 1 AND position IN (...)`. Smaller values mimic more fragmented request streams; larger values test global batching across the chromosome.

`sort_position_keys = true`

Sorts and deduplicates cold positions before building lookup chunks. This is sandbox benchmark logic, not a Lance scanner parameter, so Lance has no default. Sorting keeps each benchmark deterministic and tends to improve physical locality compared with probing keys in VCF/input order.

### Sidecar Parameters

`payload_backends = ["heed_payload_zstd"]`

Enables optional LMDB/heed payload caches that map chr-local `position UInt32` keys directly to decoded `--everything` payload structs. Values are `EverythingPayload { rows: Vec<EverythingRow> }`, and `EverythingRow` stores the 47 logical fields required by the sandbox `--everything` workload. Duplicate positions are preserved by storing multiple rows under the same position key. Supported values are `heed_payload_raw` and `heed_payload_zstd`; the main v2.1 config enables only `heed_payload_zstd` because raw payloads can be much larger for full chr1.

The heed payload benchmark bypasses Lance reads entirely after the payload sidecar is built. Its benchmark row reports LMDB payload decode time in `Resolve s`; Lance IO counters remain zero because no Lance scanner or `take_rows` call is used for that row.

`heed_map_size_gib = 64`

Sets the LMDB map size for heed payload sidecars. This reserves virtual address space and caps the maximum LMDB environment size; it does not preallocate that amount on disk. The sandbox default is `64` GiB. Increase it before building raw payload sidecars if LMDB reports `MDB_MAP_FULL`.

### Inspect Parameters

`include_physical_columns = true`

Reports physical Lance column metadata, including observed encodings, compressed bytes, and page counts when Lance exposes enough metadata to map them. Lance default: none; this controls sandbox reporting only. Sandbox default: `true`.

`include_index_sizes = true`

Reports Lance index sizes and metadata where available. Lance default: none; this controls sandbox reporting only. Sandbox default: `true`.

### Encoding Parameters

`structural = "miniblock"`

Sets `lance-encoding:structural-encoding`. Lance default: metadata unset.

`compression = "zstd"`

Sets `lance-encoding:compression`. Lance default: metadata unset.

`compression_level = 3`

Sets `lance-encoding:compression-level`. Lance default: metadata unset.

`dict_values_compression = "zstd"`

Sets `lance-encoding:dict-values-compression`. Lance default: metadata unset.

`dict_values_compression_level = 3`

Sets `lance-encoding:dict-values-compression-level`. Lance default: metadata unset.

`rle_threshold = 0.95`

Sets `lance-encoding:rle-threshold`. Lance default: metadata unset.

`dict_size_ratio = 0.99`

Sets `lance-encoding:dict-size-ratio`. Lance default: metadata unset.

`dict_divisor = 1`

Sets `lance-encoding:dict-divisor`. Lance default: metadata unset.

`minichunk_size = 4096`

Sets the Lance `lance-encoding:minichunk-size` field metadata applied by the sandbox builder. Lance's metadata default is unset; for binary miniblock encoding, Lance then falls back to the `LANCE_BINARY_MINIBLOCK_CHUNK_SIZE` environment variable or `4096` bytes if that is unset/invalid. The sandbox configs set `4096` explicitly to match that default fallback because the config schema currently requires a numeric value. In Lance 2.1, values `>= 32768` are rejected and the default is used; Lance 2.2 supports larger values. Larger values generally improve compression and scan efficiency by amortizing metadata over more rows; smaller values can improve point-lookup locality because fewer unrelated rows are tied to a selected binary/FSST minichunk.

Field-specific `[fields.<name>.encoding]` blocks may override any of the encoding fields above for one logical field. Lance default: metadata unset for fields without overrides.

### Struct Packing Parameters

`[structs.<name>]`

Defines an optional Lance 2.2 struct-packing group. Lance default: no struct packing.

`enabled = true`

Enables the group. Sandbox default: `false`.

`packed_metadata = true`

Applies the group's selected encoding metadata to the packed struct field in addition to child fields. Lance default: no sandbox metadata propagation. Sandbox default: `true`.

`fields = ["AF", "..."]`

Lists top-level logical fields to move into the struct. Lance default: no field packing.

`packed_targeted_v22.toml`

Adds a benchmark-only Lance 2.2 profile that keeps `position`, `tier`, `chrom`, `start`, `variant_keys`, and constant/null helper fields top-level while packing always-read `--everything` payload groups such as `identity_text`, `clinical_payload`, `freq_gnomade`, and split `freq_gnomadg_*` structs. Use `inspect` after building and confirm the observed encodings include `VariablePackedStruct`; otherwise the config only changed the logical Arrow schema and did not produce the intended Lance 2.2 physical layout.

Known upstream limitation: Lance 7.0.0, also reproduced with 6.0.1, can panic when a Lance 2.2 packed struct contains nullable children that request per-value compression. The observed panic is `Per-value compression not yet supported for block type: Nullable` from `rust/lance-encoding/src/compression.rs`. GitHub issue searches for the exact panic and for packed nullable compression currently return no matching Lance issue:

- <https://github.com/lance-format/lance/issues?q=is%3Aissue+%22Per-value+compression+not+yet+supported+for+block+type%3A+Nullable%22>
- <https://github.com/lance-format/lance/issues?q=is%3Aissue+%22lance-encoding%3Apacked%22+Nullable+compression>
- <https://github.com/lance-format/lance/issues?q=is%3Aissue+packed+struct+nullable+compression>

Related upstream work exists for Lance format 2.1/2.2 compression and packed structs, but it is not the same bug: the Lance 2.1 umbrella tracks compression and structural encodings, issue #2862 tracks variable-width packed structs, and issue #3111 tracks nullable fixed-size-list zipped compression. Until Lance handles nullable blocks in the packed/per-value compression path, sandbox packed-struct configs should either avoid nullable compressed children or strip compression metadata from nullable children inside packed structs before building.

### Lance Environment Knob

`LANCE_MINIBLOCK_MAX_VALUES`

This is not a sandbox config field; it is a Lance process environment variable used by the primitive miniblock encoder. Lance's default is `4096` values, and values are clamped to `1..=4096`. It caps the number of primitive values in a miniblock chunk, while Lance also caps the compressed chunk payload to slightly under `8 KiB`. No sandbox parameter maps directly to this value today. It is related to `minichunk_size` only at the broad encoding-granularity level: `LANCE_MINIBLOCK_MAX_VALUES` limits primitive miniblock value count, while `minichunk_size` is field metadata consumed by binary/FSST miniblock chunking.

## Fixed Sandbox Behavior

- Warm and cold writes use Lance's default `WriteParams.max_rows_per_file`, currently `1_048_576` rows in Lance 7.0.0. There are no `warm_fragment_rows` or `cold_fragment_rows` parameters because this sandbox should follow Lance defaults unless a specific file-size sweep is added later.
- Warm, cold, and key-extraction scans use Lance's default scanner batch size. Lance uses `LANCE_DEFAULT_BATCH_SIZE` when that environment variable is valid; otherwise it falls back to `max(object_store.block_size / 4, 8192)`, commonly `16384` with a 64KiB object-store block size. There are no `scan_batch_size` or `warm_scan_batch_size` parameters.
- Benchmarks run a key-only cold index prewarm before the measured cold lookup sweep. The prewarm uses the first configured `lookup_batch_sizes` value and projects only the lookup key column, so it warms Lance scalar-index/tier-filter paths without preloading the full `--everything` payload. The prewarm result is reported separately and is not included in `cold_results`.
- Benchmarks also compare two row-id materialization paths. The BTree path reads Lance scalar-index `page_data.lance` once, builds a sorted `position UInt32 -> row_id` map, and then calls `Dataset::take_rows` for the `--everything` projection. The external flat sidecar at `reports/position_row_ids.bin` stores the same map for reuse. Existing sidecars are validated against the BTree page-data pairs on load; fresh sidecars are built atomically and validated against a dataset scan before being written.
- Optional heed payload sidecars live under `reports/payload_heed_*.lmdb`. They store the full 47-field `--everything` row payload per position, not Lance row ids, so their measured lookup rows do not call `Dataset::take_rows`.
- Position lookup in `mode = "position"` is always derived from `start`; there is no configurable `key.source`.
- Benchmarks always read the VEPyr `--everything` logical projection and fail if any required field cannot be resolved to a top-level column or packed struct. There is no configurable benchmark projection.
- The removed inspect toggles were reporting-only switches with no effect on dataset layout; the active inspect controls are `include_physical_columns` and `include_index_sizes`.

## Notes

- Only Lance `2.1` and `2.2` are supported.
- Optimized sandbox configs use `position UInt32` derived from `start`.
