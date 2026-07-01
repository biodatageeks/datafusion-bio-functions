# Regenerate all variation-cache datasets into the chr1 bundled layout — Design (2026-06-18)

## Goal

Make every per-contig Lance dataset under
`/Users/mwiewior/workspace/data_vepyr/115_GRCh38_merged/variation.lance/` share a
single, uniform schema — the **chr1 bundled layout** — so the whole variation cache
can be read as one polars `LazyFrame` (a trivial vertical `concat`), and so the
annotation pipeline reads one consistent format across all contigs.

The original trigger was: *"can we read all the per-chromosome datasets into a single
polars LazyFrame?"* The blocker is that the datasets currently have **two different
schemas**. This design eliminates that split by regenerating the stale datasets.

## Background / current state

- `variation.lance/` holds **463 separate Lance datasets**: 25 main contigs
  (`chr1`–`chr22`, `chrX`, `chrY`, `chrMT`) plus ~438 patch/scaffold contigs
  (`HG109_PATCH.lance`, `GL000008.2.lance`, …), each its own dataset with its own
  `_indices`. Total size **70 GB**; the volume has **~52 GB free**.
- The datasets have **exactly 2 distinct schemas**:
  - **NEW (chr1 only)** — 21 columns: `chrom, start, end, allele_string, failed,
    variation_name, clin_sig, clin_sig_allele, clinical_impact, phenotype_or_disease,
    pubmed, somatic, minor_allele, minor_allele_freq, clinvar_ids, cosmic_ids,
    dbsnp_ids, tier (int8), af_global (Utf8), af_gnomade (Utf8), af_gnomadg (Utf8)`.
    The 27 individual population AF columns are **bundled** into the 3 pipe-joined
    `Utf8` columns; a derived `tier (int8)` column drives warm/cold ordering.
  - **OLD (the other 462)** — 44 columns: **no `tier`**, ~30 exploded population AF
    columns (`AF, AFR, AMR, EAS, EUR, SAS, gnomADe, gnomADe_AFR … gnomADg_REMAINING`).

### Compatibility — confirmed both sides (no code changes needed)

- **Producer.** `build_lance_variation_chrom`
  (`datafusion/bio-function-vep/examples/build_lance_variation_chrom.rs` →
  `src/lance_cache/build.rs:152`) emits the chr1 bundled layout **unconditionally**.
  `streaming_lance_output_schema` (build.rs:1179) bundles via
  `af_bundle::bundle_schema`, and `write_variation_tier_pass` (build.rs:359) bundles
  each batch + appends `tier`. The 3 AF columns are scalar **`Utf8`** (miniblock+zstd),
  not List/struct (af_bundle.rs module doc + `concat_field`). No flag produces the old
  exploded layout — that came from the *separate legacy* path
  `examples/build_lance_variation_cache.rs` →
  `warm_cache::lance_variation::write_merged_lance_variation_dataset`. So the 462 stale
  datasets are simply pre-bundling builds.
- **Consumer.** `vepyr` (`/Users/mwiewior/research/git/vepyr`) has **no
  variation-reader of its own**. It depends on this repo via a **local PATH
  dependency** (`vepyr/Cargo.toml:32` → `../datafusion-bio-functions/datafusion/
  bio-function-vep`, features `cache-builder,lance-cache`) and forwards
  `cache_format="lance"` + an options JSON to the upstream runtime. All
  bundled→CSQ unbundling lives upstream (`lance_cache/af_bundle.rs`,
  `variation_runtime.rs::resolve_and_take`) — the same path that already reads
  `chr1.lance` green. vepyr's `_core.abi3.so` was rebuilt against current upstream src.
  Regenerating all contigs to the bundled layout will **not** break vepyr wiring; the
  only references to exploded population names in vepyr are VEP **CSQ output** field
  expectations in `tests/_golden_suite.py`, validated by the e2e parity gate below.

### Naming round-trip (verified)

`canonical_chrom_label` (`manifest.rs:97`): `1..22 → chr1..chr22`, `X/Y → chrX/chrY`,
`M/MT → chrMT`, all other labels (patches/scaffolds) returned unchanged.
`dataset_dir_name` = canonical label + `.lance` (non-`[A-Za-z0-9_.-]` bytes
hex-encoded). The builder maps `--chrom` → `source_chrom` by stripping a leading `chr`
(build.rs:157). Therefore, for each existing `<name>.lance`, invoking
`--chrom <name>` rebuilds the same dataset. All 463 current names are clean
(alnum / `_` / `.` / `-`), so the round-trip holds.

## Approach (chosen: A — shell driver, no Rust changes)

The per-contig builder function rebuilds its DataFusion context internally on each
call, so a single-process internal loop saves little wall-clock over invoking the
existing binary per contig. The shell driver is therefore the simplest robust option:
zero new Rust code, naturally resumable, manifest stays consistent incrementally,
exact 463 scope, disk-safe. (Fallback if per-contig context registration proves to
dominate wall-clock: promote to a checkpointed single-process binary — "Approach C".)

### Component 1 — driver script `scripts/regenerate_variation_cache.sh`

1. Build once:
   `cargo build --release --example build_lance_variation_chrom --features lance-cache,cache-builder`.
2. Enumerate contigs: list immediate `*.lance` subdirs of `variation.lance/`, strip the
   `.lance` suffix → the `--chrom` values (463 of them).
3. For each contig, before rebuilding, inspect the existing dataset schema; if it
   already has `tier` **and** `af_global` (new layout) and `--force` is not set,
   **skip** (so `chr1` is skipped and a crashed run resumes cheaply).
4. Otherwise run:
   ```
   ./target/release/examples/build_lance_variation_chrom \
     --cache-root  /Users/mwiewior/workspace/data_vepyr/homo_sapiens_merged/115_GRCh38 \
     --output-dir  /Users/mwiewior/workspace/data_vepyr/115_GRCh38_merged \
     --chrom <name> --cache-source-type merged --partitions 8 --overwrite
   ```
   The binary removes-then-rewrites only that dataset and **upserts its manifest
   entry**.
5. Append `<contig> OK rows=<n> elapsed=<s>` or `<contig> FAIL <err>` to
   `regen_progress.log`. Failures do not abort the run; they are summarized at the end
   and the script is re-runnable.

### Data flow

merged Ensembl source (`homo_sapiens_merged/115_GRCh38/<source_chrom>/`) →
DataFusion scan → tier-classify (`WARM_AF_THRESHOLD=0.01`, `WARM_POSITION_RADIUS=1`) →
bundle 27 AF cols → 3 `Utf8` → Lance **V2_2** write + `start` BTree & `tier` bitmap
indices → `variation.lance/<contig>.lance` → manifest upsert.
Consumer: vepyr → upstream `resolve_and_take` (`bundle_projection` take + `unbundle`) →
CSQ output. Unchanged from the chr1 path.

## Disk & safety

- **In-place, per-contig overwrite** is the only disk-viable option (a full
  side-by-side rebuild needs ~70 GB; only ~52 GB free). The builder removes each
  dataset immediately before rewriting it, so peak extra disk ≈ one contig (a few GB
  for the largest), well under the 52 GB headroom. The bundled layout is *smaller*
  than the exploded one for the 462, so net cache size **drops**.
- This overwrites **production in place** (approved). Resumability + per-contig
  manifest upsert keep the cache consistent: at any point each dataset is either fully
  old or fully new, and its manifest entry matches.

## Error handling

- **Per-contig failure** (0 rows, missing source partition): the single-contig binary
  errors on 0 rows; the driver logs it, continues, and reports at the end. The 463
  existing contigs all had rows, so 0-row results are flagged as anomalies, not
  silently skipped.
- **Crash mid-contig:** the next run's `--overwrite` removes the partial dataset first.
- **Resume:** schema check skips already-migrated contigs unless `--force`.

## Validation (gating — must pass before declaring done)

1. **Staged dry-run first:** run only 3 contigs — one patch + `chrMT` + `chr21` — to
   confirm naming round-trip, output schema, disk behavior, and e2e parity **before**
   the full 463 run.
2. **Row-count parity:** capture each contig's old row count before regenerating;
   assert new == old (no data loss).
3. **Schema uniformity:** assert all 463 dataset schemas equal chr1's 21-column schema
   (presence of `tier int8`, `af_global/af_gnomade/af_gnomadg` `Utf8`, absence of the
   exploded population columns).
4. **Single LazyFrame (the original ask):**
   `pl.concat([scan(d) for d in all_463], how="vertical")` builds without schema
   conflict and its total row count equals the sum of per-dataset counts.
5. **e2e annotation parity:** `cd vepyr && env -u VIRTUAL_ENV -u CONDA_PREFIX uv run
   python e2e-testing/scripts/run_annotation_fast.py <contig> --cache merged --backend
   lance` on `chr21` + one patch; compare CSQ/AF vs VEP golden (0 mismatches). Confirms
   bundled→CSQ unbundling holds beyond chr1.

## Non-goals

- **No builder code changes** — the builder already emits the target layout.
- **The 2.4 GB AF-size gate** from `2026-06-17-af-bundle-resume.md` — that gate belongs
  to the abandoned fullzip `List<Utf8>` experiment, not the shipped scalar-`Utf8`
  layout. Out of scope here.
- **The ~1453 empty source contigs** not currently cached — scope is exactly the 463
  existing datasets.
- **The legacy merged-Parquet builder** (`build_lance_variation_cache.rs`) — untouched.

## Deliverables

- `scripts/regenerate_variation_cache.sh` — the driver (enumerate, skip-migrated,
  per-contig build with `--overwrite`, checkpoint log, end-of-run summary).
- A validation script (or section of the driver) implementing checks 2–4 above.
- A short run log / note recording the staged dry-run result and the final full-run
  summary (counts, failures, before/after total cache size).
