# Plan: Populate VEP CSQ Fields Beyond the Default 29

## Status Summary

**Completed** (commit `99a12a9`): All 29 default CSQ fields match VEP at 100% on full chr22 (715,133 entries). Fields that VEP leaves empty by default (HGVSc, HGVSp, Existing_variation, SOURCE, motif fields) are correctly left empty. HGVS computation code exists in `hgvs.rs` but is not emitted until `--hgvs` mode is enabled.

**Upstream update** (`datafusion-bio-formats` issue #121 comment, 2026-03-10; PR #75, commit `5888ba9`): `datafusion-bio-format-ensembl-cache` now exposes the following promoted top-level columns in both text and storable decoding paths:
- transcript: `gene_phenotype`, `ccds`, `swissprot`, `trembl`, `uniparc`, `uniprot_isoform`, `cds_start_nf`, `cds_end_nf`, `mature_mirna_regions`
- motif: `transcription_factors`

Downstream consequence: `datafusion-bio-functions` should prefer direct Arrow projection of those columns over `raw_object_json` parsing. `raw_object_json` should remain only for fields that are still intentionally unpromoted (for example SIFT / PolyPhen / DOMAINS) or as a short-term fallback during rollout. One caveat remains: upstream does not expose a `flags_str` / raw `FLAGS` order-preserving column, so exact VEP encounter order still needs either canonical reconstruction or a JSON fallback.

## VEP Field Discovery

VEP's default output (with `--regulatory`, no extra flags) produces **28 fields** (our 29 minus SOURCE). The `--everything` flag expands to **80 fields**. Each extra field requires a specific VEP flag to enable.

### Cache Parquet Schema (available top-level columns)

| Parquet File | Key Columns Available |
|-------------|----------------------|
| `variation` | `variation_name`, `AF`, `AFR`..`SAS`, `gnomADe`+pops, `gnomADg`+pops, `clin_sig`, `somatic`, `phenotype_or_disease`, `pubmed`, `minor_allele`, `minor_allele_freq` |
| `transcript` | `stable_id`, `version`, `is_canonical`, `tsl`, `mane_select`, `mane_plus_clinical`, `translation_stable_id`, `source`, `gene_phenotype`, `ccds`, `swissprot`, `trembl`, `uniparc`, `uniprot_isoform`, `cds_start_nf`, `cds_end_nf`, `mature_mirna_regions`, `raw_object_json` (still needed only for exact `FLAGS` encounter order fallback plus `_variation_effect_feature_cache.protein_function_predictions` / `_variation_effect_feature_cache.protein_features`) |
| `translation` | `stable_id`, `version`, `transcript_id`, `translation_seq`, `cds_sequence` |
| `motif` | `motif_id`, `strand`, `score`, `binding_matrix`, `transcription_factors`, `raw_object_json` (chr22 has 0 rows; JSON only as fallback) |

---

## Batch 1: Trivial — Top-Level Cache Columns, No JSON Parsing

**Fields**: VARIANT_CLASS, CANONICAL, TSL, MANE_SELECT, MANE_PLUS_CLINICAL, ENSP, GENE_PHENO, CCDS, SWISSPROT, TREMBL, UNIPARC, UNIPROT_ISOFORM, TRANSCRIPTION_FACTORS

| CSQ Field | VEP Flag | Cache Source | Format | Notes |
|-----------|----------|-------------|--------|-------|
| VARIANT_CLASS | `--variant_class` | **Computed from REF/ALT** | `"SNV"`, `"deletion"`, `"insertion"`, `"indel"`, `"substitution"` | No cache needed; pure string logic on allele lengths |
| CANONICAL | `--canonical` | `transcript.is_canonical` (bool) | `"YES"` or empty | Already loaded in `load_transcripts` schema |
| TSL | `--tsl` | `transcript.tsl` (int32) | `"1"`..`"5"` or empty | Already in parquet as top-level column |
| MANE_SELECT | `--mane` | `transcript.mane_select` (string) | e.g. `"NM_000185.4"` | Already in parquet |
| MANE_PLUS_CLINICAL | `--mane` | `transcript.mane_plus_clinical` (string) | e.g. `"NM_000185.4"` | Already in parquet |
| ENSP | `--protein` | `transcript.translation_stable_id` (string) | e.g. `"ENSP00000215727"` | Already in parquet; we already extract `translation.stable_id` |
| GENE_PHENO | `--gene_phenotype` | `transcript.gene_phenotype` (bool) | `"1"` or empty | Promoted upstream; absent attribute materializes as `false` |
| CCDS | `--ccds` | `transcript.ccds` (string) | e.g. `"CCDS13783.1"` | Promoted upstream |
| SWISSPROT | `--uniprot` | `transcript.swissprot` (string) | e.g. `"P05546.236"` | Promoted upstream |
| TREMBL | `--uniprot` | `transcript.trembl` (string) | e.g. `"-"` (dash = none) | Promoted upstream |
| UNIPARC | `--uniprot` | `transcript.uniparc` (string) | e.g. `"UPI000012C603"` | Promoted upstream |
| UNIPROT_ISOFORM | `--uniprot` | `transcript.uniprot_isoform` (string) | e.g. `"P05546-1"` | Promoted upstream |
| TRANSCRIPTION_FACTORS | default | `motif.transcription_factors` (string) | e.g. `"CUX1"` | Promoted upstream; alias differences are normalized upstream (`_transcription_factors` / `transcription_factors`) |

### Implementation

1. Add fields to `TranscriptFeature`: `is_canonical: bool`, `tsl: Option<i32>`, `mane_select: Option<String>`, `mane_plus_clinical: Option<String>`, `translation_stable_id: Option<String>`, `gene_phenotype: bool`, `ccds: Option<String>`, `swissprot: Option<String>`, `trembl: Option<String>`, `uniparc: Option<String>`, `uniprot_isoform: Option<String>`
2. Add `transcription_factors: Option<String>` to `MotifFeature`
3. Extract directly in `load_transcripts` / `load_motif_features` from top-level cache columns
4. Add to `TranscriptConsequence` or emit directly in the CSQ format string from `ctx.transcripts[idx]` / motif context
5. VARIANT_CLASS: compute in `annotate_batch_with_transcript_engine` from REF/ALT before the per-transcript loop
6. Gate each field on an options_json flag where VEP makes it opt-in (for example `"canonical": true`, `"gene_phenotype": true`, `"ccds": true`, `"uniprot": true`)
7. Keep JSON reads only as a short-term fallback during rollout, not as the primary implementation path

---

## Batch 2: Runtime Refactor — Replace Existing JSON Helpers with Promoted Columns

These are not new CSQ output fields by themselves, but they directly affect the transcript consequence engine and remove two current `raw_object_json` parsing helpers from the hot path.

| Runtime Need | Cache Source | Current Downstream Helper | Semantics / Caveat |
|--------------|--------------|---------------------------|--------------------|
| Mature miRNA genomic regions | `transcript.mature_mirna_regions` (`List<Struct<start: Int64, end: Int64>>`) | `parse_mirna_regions_from_json()` | `NULL` for non-miRNA transcripts; populated only for miRNA transcripts |
| CDS incomplete flags for runtime logic | `transcript.cds_start_nf`, `transcript.cds_end_nf` (bool) | `parse_transcript_flags()` | Materialized booleans; when absent upstream the value is `false` |
| Exact `FLAGS` string ordering | not promoted upstream | `parse_transcript_flags()` | If exact VEP attribute encounter order still matters, keep a JSON fallback or reconstruct canonically from booleans |

### Implementation

1. Extend Arrow projection in `load_transcripts` to request `mature_mirna_regions`, `cds_start_nf`, and `cds_end_nf`
2. Populate `TranscriptFeature` from those columns instead of `parse_mirna_regions_from_json()` / `parse_transcript_flags()` for runtime logic
3. Decide the `flags_str` strategy explicitly: reconstruct canonically from booleans, or keep `raw_object_json` only for that one string fallback
4. Keep `raw_object_json` parsing for genuinely unpromoted payloads only (SIFT / PolyPhen / DOMAINS), not transcript flags or miRNA regions

---

## Batch 3: Medium — Allele-Aware Frequency Extraction from Variation Cache

**Fields**: Existing_variation, AF, AFR_AF..SAS_AF, gnomADe_AF+pops, gnomADg_AF+pops, MAX_AF, MAX_AF_POPS, CLIN_SIG, SOMATIC, PHENO, PUBMED

### Sub-group 3A: Existing_variation + Clinical (per-variant, not per-transcript)

| CSQ Field | VEP Flag | Cache Column | Notes |
|-----------|----------|-------------|-------|
| Existing_variation | `--check_existing` (auto with cache) | `variation.variation_name` | Requires co-located variant lookup with allele matching; VEP emits per-transcript |
| CLIN_SIG | `--check_existing` | `variation.clin_sig` | Already exposed as `cache_clin_sig` output column |
| SOMATIC | `--check_existing` | `variation.somatic` | Already in cache |
| PHENO | `--check_existing` | `variation.phenotype_or_disease` | Already in cache |
| PUBMED | `--check_existing` | `variation.pubmed` | Already in cache |

### Sub-group 3B: Allele Frequencies (31 fields)

| CSQ Field Group | VEP Flag | Cache Columns | Format in Cache |
|----------------|----------|---------------|-----------------|
| AF | `--af` / `--af_1kg` | `variation.AF` | `"allele:freq"` e.g. `"T:0.0034"` |
| AFR_AF..SAS_AF (5) | `--af_1kg` | `variation.AFR`..`variation.SAS` | Same format |
| gnomADe_AF + 9 pops | `--af_gnomade` | `variation.gnomADe`, `variation.gnomADe_AFR`.. | Same format |
| gnomADg_AF + 9 pops | `--af_gnomadg` | `variation.gnomADg`, `variation.gnomADg_AFR`.. | Same format; 99.2% populated |
| MAX_AF | `--max_af` | **Derived** from all above | Max frequency across all populations |
| MAX_AF_POPS | `--max_af` | **Derived** | Population name(s) with max AF |

### Implementation

1. AF extraction requires: match variant allele against `"allele:freq"` string → extract freq for matching allele
2. Extend `lookup_variants` join to carry AF columns from variation cache into the annotated batch
3. Parse `"allele:freq"` format per-allele in `annotate_batch_with_transcript_engine`
4. MAX_AF/MAX_AF_POPS: compute across all extracted AF values
5. Gate on options_json flags: `"af": true`, `"af_1kg": true`, `"af_gnomade": true`, `"af_gnomadg": true`, `"max_af": true`

---

## Batch 4: Hard — Requires New Processing Logic

### 4A: HGVSc + HGVSp (code exists, needs fixes + opt-in wiring)

| CSQ Field | VEP Flag | Status |
|-----------|----------|--------|
| HGVSc | `--hgvs` + `--fasta` | Code in `hgvs.rs` computes exonic + intronic notation. **Needs fixes**: use `n.` for non-coding transcripts (vs `c.`); reverse-complement alleles for negative strand |
| HGVSp | `--hgvs` + `--fasta` | Code in `hgvs.rs` handles missense/synonymous/stop/frameshift/inframe. Currently computed but not emitted |
| HGVS_OFFSET | `--hgvs` | Not yet implemented. Reports position shift when VEP normalizes the variant |

**Remaining work:**
- [ ] Fix `n.` prefix for non-coding transcripts (check `is_non_coding_biotype()`)
- [ ] Reverse-complement alleles in HGVSc for negative-strand transcripts
- [ ] Add options_json flag `"hgvs": true` to gate emission
- [ ] Wire `tc.hgvsc` / `tc.hgvsp` into CSQ format when flag is set

### 4B: SIFT + PolyPhen (compressed binary prediction matrices)

| CSQ Field | VEP Flag | Cache Location |
|-----------|----------|---------------|
| SIFT | `--sift b` | `transcript.raw_object_json._variation_effect_feature_cache.protein_function_predictions.sift` |
| PolyPhen | `--polyphen b` | Same path, `polyphen_humdiv` key |

**Complexity**: Prediction data is a gzip-compressed `Bio::EnsEMBL::Variation::ProteinFunctionPredictionMatrix` — a binary matrix indexed by (amino acid position, substitution). Need to:
1. Decompress gzip blob from JSON
2. Decode VEP's custom binary matrix format
3. Look up (protein_position, alt_amino_acid) → (score, prediction_label)
4. Format: `"tolerated(0.23)"` / `"deleterious(0.01)"` for SIFT; `"benign(0.001)"` / `"possibly_damaging(0.5)"` / `"probably_damaging(0.95)"` for PolyPhen

**Alternatively**: Could extract these during cache ETL into top-level columns per-position, but that would massively expand the parquet size.

### 4C: DOMAINS (protein feature overlap)

| CSQ Field | VEP Flag | Cache Location |
|-----------|----------|---------------|
| DOMAINS | `--domains` | `transcript.raw_object_json._variation_effect_feature_cache.protein_features` |

**Implementation**: JSON array of `Bio::EnsEMBL::ProteinFeature` objects with `start`, `end`, `hseqname`, `analysis._display_label`. Check if variant protein position overlaps any feature range. Format: `"source:id"` e.g. `"PDB-ENSP_mappings:1jmj.A&Pfam:PF00079"`.

### 4D: Motif Fields (no test data on chr22)

| CSQ Field | VEP Flag | Difficulty |
|-----------|----------|-----------|
| MOTIF_NAME | default | Easy — `motif.motif_id` |
| MOTIF_POS | default | Easy — positional computation |
| HIGH_INF_POS | default | Hard — requires PWM information content scoring |
| MOTIF_SCORE_CHANGE | default | Hard — requires PWM `score(alt) - score(ref)` |

`TRANSCRIPTION_FACTORS` moved out of this batch: upstream now exposes it as top-level `motif.transcription_factors`, so downstream should project it directly in Batch 1 instead of reading motif JSON.

**Blocker**: Only chr22 motif parquet exists and has 0 rows. Need to generate motif parquets for other chromosomes to test. `binding_matrix` column contains the PWM ID but the actual matrix weights may need to come from the VEP cache's binding_matrix table.

---

## Implementation Order

```
DONE ─────────────────────────────────────────────────────────────────
  Phase 1A  Existing_variation wiring (cache-hit path)           [x]
  Phase 1B  SOURCE (gated on --merged)                           [x]
  Phase 2A  TranscriptFeature.version, TranslationFeature.stable_id/version  [x]
  Phase 2B  hgvs.rs: format_hgvsc (exonic + intronic)           [x]
  Phase 2C  hgvs.rs: format_hgvsp (all variant types)           [x]
  Phase 2D  TranscriptConsequence.hgvsc/hgvsp computed           [x]
  Benchmark All 29 fields 100% match on chr22 (715k entries)     [x]

Batch 2 — Runtime refactor from promoted columns ──────────────────
  [x] mature_mirna_regions  read from List<Struct<start,end>>
  [x] cds_start_nf/cds_end_nf  read from promoted booleans
  [x] flags_str             canonical reconstruction from booleans
  [x] Removed raw_object_json parsing + serde_json dep           4cedfc8

Batch 1 — Top-level cache columns (CSQ positions 29-40) ──────────
  [x] VARIANT_CLASS         classify_variant() from VEP-minimized REF/ALT
  [x] CANONICAL             transcript.is_canonical → "YES" or ""
  [x] TSL                   transcript.tsl (upstream fix: bio-formats#122)
  [x] MANE_SELECT           transcript.mane_select
  [x] MANE_PLUS_CLINICAL    transcript.mane_plus_clinical
  [x] ENSP                  transcript.translation_stable_id
  [x] GENE_PHENO            transcript.gene_phenotype → "1" or ""
  [x] CCDS                  transcript.ccds
  [x] SWISSPROT             transcript.swissprot
  [x] TREMBL                transcript.trembl
  [x] UNIPARC               transcript.uniparc
  [x] UNIPROT_ISOFORM       transcript.uniprot_isoform
  [x] TRANSCRIPTION_FACTORS motif.transcription_factors (pos 27, 0 motif rows on chr22)
  Benchmark All 41 fields 100% match on chr22 (2421 entries)     [x]

TODO ─────────────────────────────────────────────────────────────────
Batch 3 — Medium (allele-aware freq extraction) ─────────────────────
  [ ] Existing_variation    co-located variant lookup (per-transcript CSQ)
  [ ] CLIN_SIG/SOMATIC/PHENO/PUBMED  variation cache (already exposed)
  [ ] AF + 1000G pops (6)  variation.AF/AFR/AMR/EAS/EUR/SAS
  [ ] gnomADe + pops (10)  variation.gnomADe + 9 sub-populations
  [ ] gnomADg + pops (10)  variation.gnomADg + 9 sub-populations
  [ ] MAX_AF + MAX_AF_POPS  derived from all AF fields

Batch 4 — Hard (new processing logic) ───────────────────────────────
  [ ] HGVSc/HGVSp          fix n. prefix + strand RC; gate on --hgvs
  [ ] HGVS_OFFSET           position shift tracking
  [ ] SIFT                  binary prediction matrix decompression
  [ ] PolyPhen              binary prediction matrix decompression
  [ ] DOMAINS               protein feature overlap check
  [ ] Motif fields (4)      needs non-chr22 cache + PWM scoring engine
```

## Critical Files

| File | Batch 1 | Batch 2 | Batch 3 | Batch 4 |
|------|---------|---------|---------|---------|
| `src/transcript_consequence.rs` | Add transcript/motif CSQ metadata fields | Consume promoted miRNA/NF context fields | — | HGVSc fixes |
| `src/annotate_provider.rs` | Extract new top-level transcript/motif columns; emit in CSQ | Retire JSON helpers except optional `flags_str` fallback | Extend lookup join for AF cols; parse allele:freq | Wire HGVSc/HGVSp emission |
| `src/hgvs.rs` | — | — | — | Fix n./c. prefix, strand RC |
| `src/golden_benchmark.rs` | Extend CSQ_FIELD_NAMES if fields are added | — | — | — |
| `examples/annotate_vep_golden_bench.rs` | Add VEP flags for comparison | — | — | Add --hgvs + --fasta |
