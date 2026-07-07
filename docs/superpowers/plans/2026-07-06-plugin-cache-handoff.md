# Custom VEP Plugin Caches — Handoff

**Date:** 2026-07-06
**Branch:** `feat/plugin-cache-alphamissense` (pushed) — **PR #190**
**Companion repo:** `biodatageeks/vepyr-plugins` (private, pushed) — AlphaMissense manifest
**Status:** Complete & verified. AlphaMissense achieves **exact chr22 golden parity** (1,912/1,912 populated CSQ lines match, 0 mismatch, 0 over-emit).

---

## 1. What this delivers

A reusable subsystem for adding **custom VEP plugins** (CADD/ClinVar/dbNSFP/AlphaMissense-style per-variant and per-transcript scores) as **frequency-tiered per-chromosome Parquet caches** whose values are emitted as VEP CSQ output fields — proven end-to-end with **AlphaMissense**.

Design spec: `docs/superpowers/specs/2026-07-05-custom-vep-plugin-caches-design.md`
Plans: `docs/superpowers/plans/2026-07-05-alphamissense-plugin-prototype.md` and `…-plugin-buffer-batched-lookup.md`

## 2. Architecture (module `datafusion/bio-function-vep/src/plugin_cache/`)

Build pipeline (declarative, no per-plugin Rust for the common case):

```
TOML source manifest (vepyr-plugins)
  → provider.rs        register raw source (CSV/TSV/Parquet); gzip→temp decompress
  → ingest_sql view    map raw cols → (chrom,start,end,allele_string,[match cols],values)
  → normalize.rs       canonical_contig UDF + coordinate shift (bare Ensembl, 1-based)
  → join.rs            LEFT JOIN variation shard on (chrom,start,allele_string) → tier (AF discarded)
  → write.rs           tiered shard plugin/<name>/<chrom>.parquet (reuses point_lookup_writer_properties)
  → cache_manifest.rs  emit manifest.json (value→CSQ map, match cols, per-chrom counts)
```

Runtime (buffer-batched, page-scoped — mirrors the variation lookup, `variation_lookup.rs` untouched):

```
registry.rs   PluginRegistry::open(cache_root, chrom)  (discover manifests, open shards)
lookup.rs     PluginLookup::take_buffer(&starts)        3-phase PageDir take, reads only the buffer's pages
              PluginBufferSlice::probe(start, allele, &[discriminator])  sync per-transcript filter
registry.rs   BufferSlices::probe_all(start, allele, &EngineAttrs) → Vec<PluginScalar>
csq.rs        amino_acid_change / format_scalar / field_suffix / empty_suffix
```

**Key files:**
- `source_manifest.rs` — manifest TOML types (`SourceManifest`, `MatchColumn`, `ValueColumn`, …).
- `provider.rs` — provider factory; gzip decompress-to-temp (DataFusion built w/o `compression` feature).
- `normalize.rs` — `canonical_contig` UDF + `wrap_normalization`.
- `join.rs` — `tier_sql` / `tiered_stream`.
- `write.rs` — `plugin_output_schema` + `PluginShardWriter`.
- `cache_manifest.rs` — `CacheManifest` + `discover_plugins`.
- `build.rs` — `build_plugin_chrom`.
- `lookup.rs` — `PluginLookup` (`take_buffer`), `PluginBufferSlice`, `PluginScalar`.
- `registry.rs` — `PluginRegistry`, `BufferSlices`, `EngineAttrs`.
- `csq.rs` — CSQ formatting helpers.

## 3. Cache format (matches the variation cache verbatim)

- Key columns: `chrom` (Utf8), `start`/`end` (UInt32, 1-based), `allele_string` (Utf8, `ref/alt`).
- Optional **match discriminator** column(s) after `allele_string` (§3.4) — e.g. AlphaMissense `protein_variant`.
- Then value columns, then derived `tier` (Int8; 0=warm if variation AF≥0.01, else 1=cold).
- Physical: `point_lookup_writer_properties(schema, &["tier","start"])` — ZSTD(3), no dictionary, 4 KiB/512-row pages, page-index stats, 1M-row groups, sorted `(tier,start)`.

## 4. Per-transcript matching (the crux for AlphaMissense/dbNSFP)

VEP AlphaMissense matches each transcript consequence by **amino-acid change** (`_aminoacid_changes_match`, `transcript_match=0`) and gates on **missense only**. So the plugin is per-(variant×transcript), not purely per-variant.

- Manifest declares `[[match_column]] column=protein_variant, engine_attr=amino_acid_change`.
- Runtime forms `amino_acid_change = {ref_aa}{Protein_position}{alt_aa}` (e.g. `W320R`) from the engine's `Amino_acids` (`"W/R"`) + `Protein_position` (`"320"`); probes `(start, allele, protein_variant)`.
- A non-missense line has no aa-change → `None` discriminator → probe miss → **empty output = the gate**.
- Per-variant plugins (CADD/ClinVar) declare no `match_column`: `probe(start, allele, &[])` → same value on every transcript line (VEP's variant-level behavior). Same code path, generic.

Data fact (why gating matters): on HG002 chr22, the 188 am-variants span 5,817 CSQ lines but only 1,912 are missense-populated; 3,905 must stay empty. A per-variant cache would over-populate all of them.

## 5. AlphaMissense prototype — how it was built & validated

- Source: `AlphaMissense_hg38.tsv.gz` (hg38, 1-based, chr-prefixed). Manifest in `vepyr-plugins/plugins/alphamissense/alphamissense.source.toml`.
- Built chr22 shard: `/tmp/plugin_cache/plugin/alphamissense/chr22.parquet` (1,481,489 rows).
- Golden: Ensembl VEP release 115, merged cache, `--everything --hgvs --plugin AlphaMissense` → `/Users/mwiewior/research/git/vepyr/sandbox/HG002_chr22_everything_hgvs_merged_am.vcf`. Runbook: `vepyr/e2e-testing/vep-docker.md`.
- **Parity: 1,912/1,912 populated CSQ lines match, 0 mismatch, 0 over-emit**, correct field order (`am_class`, `am_pathogenicity`) and gating. Independently re-verified.

### Reproduce (commands)

```bash
# 1. build chr22 slice (or full genome via full AlphaMissense_hg38.tsv.gz)
cargo run -q -p datafusion-bio-function-vep --example build_plugin -- \
  --manifest <vepyr-plugins>/plugins/alphamissense/alphamissense.source.toml \
  --source-path /Users/mwiewior/research/git/vepyr/sandbox/alphamissense/am_chr22.tsv.gz \
  --variation-shard /Users/mwiewior/workspace/data_vepyr/115_GRCh38_merged/variation/chr22.parquet \
  --out /tmp/plugin_cache --chrom 22

# 2. annotate with the plugin cache enabled
cargo run --release -p datafusion-bio-function-vep --example annotate_vcf -- \
  --input  /Users/mwiewior/research/git/vepyr/sandbox/HG002_chr22.vcf \
  --cache  /Users/mwiewior/workspace/data_vepyr/115_GRCh38_merged \
  --out    /tmp/vepyr_chr22_am.vcf \
  --fasta  /Users/mwiewior/workspace/data_vepyr/Homo_sapiens.GRCh38.dna.primary_assembly.fa \
  --plugin-cache /tmp/plugin_cache --everything --hgvs
```

Enable knob: `AnnotateVcfConfig.plugin_cache_root: Option<PathBuf>` (`None` = disabled, byte-identical). For `workers=1` it flows via `to_options_json` → the `annotate_vep` UDTF (`annotate_table_function.rs`); for `workers>1` via `drive_sharded_vcf_annotation`.

## 6. Key decisions & non-obvious deviations (all justified in-code)

1. **gzip → temp decompress** — DataFusion is built `default-features=false` **without `compression`** (xz2/liblzma link collision with noodles-cram, root `Cargo.toml`), so `register_csv` can't read `.gz`; `provider.rs::materialize_plain` decompresses to a temp file.
2. **Explicit value-column projection** in `wrap_normalization` (no `SELECT * EXCLUDE`) — DF-version robust.
3. **Buffer-batched lookup, not whole-file** — the first prototype loaded the whole shard into memory; replaced with `take_buffer` (page-scoped, reuses `parquet_cache::page_dir`) per the variation model.
4. **`builder.schema()` returns the FULL file schema in parquet 58** (not the projected one) — `take_buffer` derives the projected schema from the payload columns (`lookup.rs::projected_schema`); watch for this if you touch the take path.
5. **TOML ordering** — top-level scalar keys (`plugin_name`, `coordinate_system`, `ingest_sql`) MUST precede any `[[table]]` header or TOML absorbs them.
6. **Backend string** — `annotate_vep` now accepts `"parquet"` (canonical) with `"lance"` as a deprecated alias; the enum variant was renamed `Lance→Parquet`. Was rejecting `"parquet"` despite Parquet being the only backend.
7. **F32 formatting** — `format!("{v}")` (shortest round-trip) reproduces `0.2199`/`0.361` exactly; no Utf8 fallback needed. If a future plugin's source strings don't round-trip through f32, declare the value column `Utf8`.

## 7. The two bugs V1 (the golden gate) caught

- **Header/body width divergence** — header listed plugin fields (from manifests) but the body emitted none, because the executed `workers=1` path (SQL/`annotate_vep` UDTF, built from `options_json`) never carried `plugin_cache_root`. Fix: emit/parse it through `options_json`.
- **Field order** swapped vs VEP — fixed by ordering manifest `value_columns` as `am_class, am_pathogenicity`.

Both are why an end-to-end parity gate is non-negotiable; unit tests alone missed them.

## 8. Adding the next plugin

Follow the `vep-add-plugin` skill (`~/.claude/skills/vep-add-plugin/`):
1. Write `plugins/<name>/<name>.source.toml` in `vepyr-plugins` (provider, input schema, `ingest_sql`, `coordinate_system`, `value_columns`, optional `match_column`, tier).
2. New *format* only → add a `ProviderKind` arm in `provider.rs`.
3. Build per-chrom (`build_plugin` example); annotate + parity-gate (`annotate_vcf` example) before production.

Per-variant plugins (CADD/ClinVar) omit `match_column`. Per-transcript protein-scored plugins (dbNSFP) reuse the AlphaMissense pattern (`protein_variant ← amino_acid_change`).

## 9. Open follow-ups (non-blocking)

- **Full-genome AlphaMissense build** — pipeline is ready; point `--source-path` at the full `AlphaMissense_hg38.tsv.gz` and build all chroms (bgzip+tabix not needed for our build; VEP golden needs it).
- **True full-pipeline width assertion** — current guard is a `to_options_json` round-trip unit test; a fixture-cache e2e test would be stronger.
- **`indel-bearing plugin`** — the runtime probe keys on `start_val` (VCF POS); equal to `input_start` for SNVs (AlphaMissense is SNV-only). An indel plugin should use `input_start`.
- **`vepyr` Python plumbing** — expose `plugin_cache_root` through the PyO3 API so vepyr's CLI/notebook can enable plugins (only the Rust `annotate_to_vcf` knob exists today).
- **Non-goals** (spec §9): per-gene/transcript-ID plugins (LOEUF/pLI/G2P) and allele-agnostic per-position tracks are out of scope for this cache shape.

## 10. State of the world

- Main repo: PR #190 (`feat/plugin-cache-alphamissense`, ~28 commits). Supersedes spec-only PR #189.
- vepyr-plugins: pushed `master` (`23f5553`) with the AlphaMissense manifest.
- Tests: 803 crate + integration/roundtrip green; clippy `-D warnings` clean; `variation_lookup.rs` untouched; no-plugin runs byte-identical.
- Local build artifacts (not committed): `/tmp/plugin_cache/…` (chr22 shard), `/tmp/vepyr_chr22_am.vcf`, `/Users/mwiewior/research/git/vepyr/sandbox/alphamissense/` (persisted AM tabix + slice + .pm).
