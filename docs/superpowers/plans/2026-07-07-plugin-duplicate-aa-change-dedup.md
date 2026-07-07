# Handoff: plugin cache duplicate aa-change dedup (first-in-file wins)

> **SHIPPED (2026-07-07, commit 3f5087b):** the **Option B** alternative was
> implemented, not the "Option A" (`_src_ord` + `QUALIFY`) plan headlined below.
> See `datafusion/bio-function-vep/src/plugin_cache/dedup.rs` (`dedup_keep_first`)
> and its wiring in `build.rs` (`build_plugin_chrom`): an in-stream `HashSet`
> keep-first over a single-partition, file-ordered normalized stream, run before
> the tier join. No `_src_ord` column. Full chr1–22 rebuilt; merged_am gate 17→0.

**Date:** 2026-07-07
**Repos/branches:** `datafusion-bio-functions` @ `feat/plugin-cache-alphamissense` (PR #190, pushed through `402ec64`); `vepyr` @ `feat/plugin-cache` (pushed `79b3970`); `vepyr-plugins` @ `master` (pushed `e4bd38d`).
**Status:** Root cause fully diagnosed + VEP behavior confirmed in code. Fix (Option A) **NOT yet implemented**.

## The finding

Full chr1–22 AlphaMissense parity gate report
`vepyr/e2e-testing/reports/fast_chr1_22_merged_am_summary_20260707_0818.md`:
**87/88 CSQ fields at 100%; 17 `am_pathogenicity` mismatches** (all other fields incl. `am_class` = 0). Exactly **two positions**:

- `chr3:184257577 C>T` (aa-change `H101Y`) — vepyr `0.0898` vs VEP `0.0431` (×4 lines)
- `chr8:12191811 C>A` (aa-change `D43Y`) — vepyr `0.1336` vs VEP `0.1206` (×13 lines)

## Root cause — overlapping-gene ambiguity, NOT a coord/gating bug

At each position two genes/UniProts map to the **same genomic variant with the same
amino-acid change** but **different scores**. AlphaMissense source rows
(`tabix AlphaMissense_hg38.bgz.tsv.gz chr3:184257577-184257577`):

```
chr3 184257577 C T  P0DPD7 ENST00000324557.8 H101Y 0.0431 likely_benign   <- VEP uses this (first)
chr3 184257577 C T  P0DPD8 ENST00000402825.7 H101Y 0.0898 likely_benign   <- vepyr used this
chr8 12191811  C A  P0C5J1 ENST00000533852.6 D43Y  0.1206 likely_benign   <- VEP
chr8 12191811  C A  Q8N7N1 ENST00000448228.6 D43Y  0.1336 likely_benign   <- vepyr
```

vepyr's match discriminator is `template = "{ref_aa}{Protein_position}{alt_aa}"` =
`H101Y` / `D43Y` — **identical across both UniProts**. So both source rows land in the
cache with the same key `(start, allele_string, protein_variant)` but different scores.
On probe, vepyr returns whichever the shard yields first; the shard's `(tier, start)`
sort does not preserve source order, so it returns the *second* row. VEP returns the
**first row in file order** for *all* transcript lines at the variant (note every
mismatch shows the same VEP value across all 4/13 lines → VEP is not resolving
per-transcript, it just takes the first aa-change match).

## VEP behavior — CONFIRMED in code (first-in-file wins)

`vepyr/sandbox/alphamissense/AlphaMissense.pm`, `run` method (~lines 236–282):

```perl
my @data = @{$self->get_data($vf->{chr}, $vf->{start} - 2, $vf->{end})};  # tabix/file order
foreach (@data) {
    next unless $self->_aminoacid_changes_match($tva, $_->{result}->{am_protein_variant});
    if ($self->{transcript_match}) { ... next if transcript mismatch ... }   # default 0 → skipped
    my $matches = get_matched_variant_alleles(...);   # ref/alt/pos/strand equivalence only
    ...
    if (@$matches) { return \%res; }   # <-- RETURNS on the FIRST matching record
}
return {};
```

- Iterates records in tabix/file order, `return`s on the **first** aa-change match.
- `transcript_match` default **0** → no transcript disambiguation (our manifest also
  doesn't set it — consistent).
- `get_matched_variant_alleles` is only allele-equivalence (indel/strand), not a
  scoring tiebreak.

So the target rule to replicate: **for a duplicate `(start, allele, protein_variant)`
key, keep the first source-order row.**

## The fix — Option A (generic, faithful) — TO IMPLEMENT

DataFusion's pipeline (CSV scan → `ingest_sql` → `normalize` → LEFT-JOIN tier →
sort-by-`start`) does **not** preserve source file order (joins/multi-partition scans
reorder; final shard sorts by `(tier, start)`). So source order must be materialized
explicitly:

1. **Single-partition read** of the raw source (`target_partitions = 1`) so the CSV
   scan yields rows in file order.
2. **Inject a monotonic `_src_ord` Int64 column at read time** — do it in Rust (read
   the batches in `provider.rs::register_sources`, append an incrementing column,
   register that as `plugin_<name>_src`). DataFusion has no guaranteed implicit row id
   (`ROW_NUMBER() OVER ()` order is undefined across partitions), so inject explicitly.
3. **Dedup keeping min `_src_ord` per key**, e.g. in the build SQL:
   ```sql
   QUALIFY ROW_NUMBER() OVER (
     PARTITION BY chrom, start, end, allele_string, <match columns>
     ORDER BY _src_ord) = 1
   ```
   (or a Rust arg-min pass), then **drop `_src_ord`** before writing the shard.

This is plugin-generic (any plugin, any discriminator set) and matches VEP's
first-in-file rule regardless of downstream reordering. `<match columns>` =
`SourceManifest::match_column_names()` (empty for per-variant plugins → dedup on
`(chrom,start,end,allele_string)`).

### Where it goes
- `plugin_cache/provider.rs::register_sources` / `materialize_plain` — single-partition
  read + append `_src_ord`. (Note: this path already returns `Vec<TempPath>`; keep that.)
- `plugin_cache/build.rs::build_plugin_chrom` and/or `join.rs` — add the dedup
  (`QUALIFY`/window) in the normalized→tiered SQL, and ensure `_src_ord` is projected
  through `wrap_normalization` and dropped before `plugin_output_schema` write.
- `plugin_cache/normalize.rs::wrap_normalization` — must carry `_src_ord` through.

### Alternative (Option B) and cross-check
- Option B: Rust `HashSet` keep-first pass over the single-partition **normalized**
  stream *before* the tier join (drop later dups). Simpler but must run before any
  join/sort.
- Cross-check: in both cases VEP's pick = alphabetically-smaller `uniprot_id`
  (`P0DPD7<P0DPD8`, `P0C5J1<Q8N7N1`), suggesting AM is globally sorted by
  pos-then-transcript/UniProt, i.e. first-in-file ≡ smallest UniProt. Do **not** rely on
  this (plugin-specific, needs `uniprot_id` in ingest); `_src_ord` is the assumption-free
  equivalent.

## Verification plan

**Order matters — the cache MUST be rebuilt before re-checking parity.** The current
`data_vepyr/plugin_cache` chr3/chr8 shards were built by the *buggy* build and still
contain both duplicate rows, so re-running the gate against them will still show the
17 mismatches. The fix is a **build-time** dedup, so it only takes effect after the
affected shards are regenerated. So: implement the fix → **rebuild the chr3 and chr8
shards** (step 2) → then verify parity on **chr3 and chr8** (step 3). Confirming the
fix works == those two chromosomes going 17 → 0 `am_pathogenicity` mismatches.

1. Unit test: build a synthetic source with two rows sharing
   `(chrom,start,allele_string, protein_variant)` but different value + different
   `_src_ord`; assert the shard keeps the first and the probe returns its value.
2. Rebuild the affected shards with the new build:
   ```
   cargo run -q -p datafusion-bio-function-vep --example build_plugin -- \
     --manifest /Users/mwiewior/workspace/vepyr-plugins/plugins/alphamissense/alphamissense.source.toml \
     --source-path <full AlphaMissense_hg38.tsv.gz> \
     --variation-cache-dir /Users/mwiewior/workspace/data_vepyr/115_GRCh38_merged \
     --out /Users/mwiewior/workspace/data_vepyr/plugin_cache --chrom 3   # and 8
   ```
   (Full source is `vepyr/sandbox/alphamissense/AlphaMissense_hg38.tsv.gz`.)
3. Re-run the gate for chr3/chr8 (or full):
   ```
   cd vepyr/e2e-testing/scripts && uv run python run_annotation_fast.py chr3 --profile merged_am --force
   uv run python run_annotation_fast.py chr8 --profile merged_am --force
   ```
   Expect **17 → 0** `am_pathogenicity` mismatches; all 88 fields 100%.
4. Full crate tests + clippy `-D warnings`; commit on `feat/plugin-cache-alphamissense`.

## Context / already done on this PR
- PR #190 review rounds C1–C3, N1–N3, R1–R2 + doc fixes all resolved (see git log).
- The multi-allelic Codex finding (`review 4641949793`) was assessed a **false positive**
  (engine is biallelic-input by contract) and replied to on the PR.
- Disk: `data_vepyr` volume runs tight (~12 GiB free after the 23 GB golden); a full
  rebuild of all chroms is unnecessary — only chr3 + chr8 shards need rebuilding to
  verify the fix.
