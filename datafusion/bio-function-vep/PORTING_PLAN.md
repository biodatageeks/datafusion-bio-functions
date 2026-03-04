# Ensembl-VEP Porting Plan for `datafusion-bio-function-vep`

This document defines the migration plan for Ensembl-VEP-style consequence annotation using unified `parquet` and `fjall` cache backends in this repository.

Reference inputs used for parity planning:
- `ensembl-vep` release 115 behavior
- `gdl-annotations-infra` commit `de03fa05122cf6af5faeceeda6ac81eb76282a92` (`modules/rust/parquet_annotator_rust`)

## Feature Comparison

| Feature | Ensembl-VEP (release 115) | `bio-functions-vep` (current) | Port Target |
|---|---|---|---|
| Runtime model | Perl CLI pipeline | DataFusion SQL UDF/UDTF | DataFusion SQL-native annotation engine |
| Known variant lookup | Yes (`--check_existing`, colocated matching) | Yes (`lookup_variants`) | Keep and integrate with consequence output |
| Match behavior | VEP allele/coordinate normalization | `exact`, `exact_or_colocated_ids`, `exact_or_vep_existing`, optional `extended_probes` | Preserve existing modes, use as `Existing_variation` source |
| Consequence engine | Full transcript/regulatory/motif consequence evaluation | Partial transcript/exon phase-2 baseline + fallback path | Full SO term evaluation with rank + impact |
| Consequence term coverage | 41/41 | Partial runtime subset (transcript/exon-focused; regulatory/SV missing) | 41/41 (or explicit unsupported list guarded by tests) |
| `most_severe_consequence` | Yes | Ranked SO term when transcript/exon context is available; fallback placeholder (`sequence_variant`) for cache-hit-only mode | Fully computed rank-min consequence |
| CSQ field generation | Full CSQ payload | Not yet | VEP-compatible CSQ fields (phased) |
| Backends | VEP caches/filesystem model | Parquet + Fjall (lookup path) | Unified backend contract for all consequence artifacts |
| SQL composability | No | Yes | Yes (primary interface) |
| Distributed query integration | No native DataFusion integration | Native DataFusion plans | Keep DataFusion-first execution path |

## Unified Cache Contract (Parquet + Fjall)

Target: both backends must expose the same logical artifacts and key semantics. Physical layout differs by backend.

### Required Logical Artifacts

| Artifact | Purpose | Minimum Required Columns | Parquet Layout (target) | Fjall Namespace (target) |
|---|---|---|---|---|
| `variation` | Existing IDs, clinical flags, allele strings | `chrom,start,end,variation_name,allele_string,clin_sig,somatic,phenotype_or_disease,pubmed` | `<root>/<chrom>/variation/*.parquet` | `variation:{chrom}:{start}:{end}` |
| `transcripts` | Transcript interval + coding context | `chrom,start,end,transcript_id,stable_id,gene_id,strand,biotype,cds_start,cds_end,is_canonical` | `<root>/<chrom>/transcripts/*.parquet` | `tx:{chrom}:{bin}` |
| `exons` | Exon/intron boundaries and exon order | `transcript_id,exon_number,start,end,strand` | `<root>/<chrom>/exons/*.parquet` | `exon:{transcript_id}` |
| `translations` | CDS/protein mapping | `transcript_id,protein_id,cdna_len,cds_len,protein_len,translation_seq` | `<root>/<chrom>/translations/*.parquet` | `translation:{transcript_id}` |
| `transcript_attribs` | Transcript ranking and flags | `transcript_id,mane,appris,tsl,ccds,canonical_flag,sift,polyphen,domains` | `<root>/<chrom>/transcript_attribs/*.parquet` | `tx_attr:{transcript_id}` |
| `gene_xrefs` | Symbol and HGNC mapping | `gene_id,symbol,symbol_source,hgnc_id` | `<root>/<chrom>/gene_xrefs/*.parquet` | `gene_xref:{gene_id}` |
| `gene_phenotypes` | Gene phenotype flag | `gene_id,gene_pheno` | `<root>/<chrom>/gene_phenotypes/*.parquet` | `gene_pheno:{gene_id}` |
| `population_af` | AF / MAX_AF / MAX_AF_POPS | `chrom,start,end,allele_key,AF,...,MAX_AF,MAX_AF_POPS` | `<root>/<chrom>/population_allele_frequencies/*.parquet` | `pop_af:{chrom}:{start}:{end}:{allele_key}` |
| `regulatory_features` | Regulatory overlap and feature metadata | `chrom,start,end,stable_id,biotype,feature_type` | `<root>/<chrom>/regulatory_features/*.parquet` | `reg:{chrom}:{bin}` |
| `motif_features` | TF motif overlap + score context | `chrom,start,end,motif_id,motif_name,tf_name,pwm_or_score_blob` | `<root>/<chrom>/motif_features/*.parquet` | `motif:{chrom}:{bin}` |
| `mirna_features` | Mature miRNA regions | `chrom,start,end,mirna_id,arm,transcript_id` | `<root>/<chrom>/mirna_features/*.parquet` | `mirna:{chrom}:{bin}` |
| `sv_features` | Structural-event overlap deltas for ablation/amplification terms | `chrom,start,end,event_type,feature_id,feature_kind` | `<root>/<chrom>/sv_features/*.parquet` | `sv:{chrom}:{bin}` |

### Backend Normalization Rules

| Rule | Requirement |
|---|---|
| Coordinates | Normalize to 1-based closed intervals at the logical contract boundary |
| Chromosome naming | Normalize aliases (`chr1` vs `1`) before key lookup |
| Allele matching | Preserve current `match_allele` + `match_allele_relaxed` behavior |
| Consequence ranking | Use VEP rank ordering for `most_severe_consequence` |
| Nullability | Missing artifact lookups must degrade to partial annotation, not hard fail |

## Consequence Coverage Matrix (All 41 SO Terms)

Legend:
- `TX`: `transcripts`
- `EXON`: `exons`
- `SEQ`: `translations` (+ CDS/codon context)
- `VAR`: `variation`
- `REG`: `regulatory_features`
- `MOTIF`: `motif_features`
- `MIRNA`: `mirna_features`
- `SV`: `sv_features`

| # | SO Term | Required Artifacts | Notes |
|---|---|---|---|
| 1 | `transcript_ablation` | `TX`, `SV` | Requires structural-event span vs full transcript overlap |
| 2 | `splice_acceptor_variant` | `TX`, `EXON` | Intron acceptor positions (-1/-2) |
| 3 | `splice_donor_variant` | `TX`, `EXON` | Intron donor positions (+1/+2) |
| 4 | `stop_gained` | `TX`, `EXON`, `SEQ` | Codon translation required |
| 5 | `frameshift_variant` | `TX`, `EXON`, `SEQ` | Coding indel frame delta |
| 6 | `stop_lost` | `TX`, `EXON`, `SEQ` | Stop codon destruction |
| 7 | `start_lost` | `TX`, `EXON`, `SEQ` | Start codon destruction |
| 8 | `transcript_amplification` | `TX`, `SV` | Structural amplification event |
| 9 | `feature_elongation` | `TX`, `SV` | Structural elongation event |
| 10 | `feature_truncation` | `TX`, `SV` | Structural truncation event |
| 11 | `inframe_insertion` | `TX`, `EXON`, `SEQ` | Coding indel, frame-preserving |
| 12 | `inframe_deletion` | `TX`, `EXON`, `SEQ` | Coding indel, frame-preserving |
| 13 | `missense_variant` | `TX`, `EXON`, `SEQ` | Amino-acid substitution |
| 14 | `protein_altering_variant` | `TX`, `EXON`, `SEQ` | Generic protein change fallback |
| 15 | `splice_donor_5th_base_variant` | `TX`, `EXON` | Donor +5 rule |
| 16 | `splice_region_variant` | `TX`, `EXON` | Donor/acceptor surrounding windows |
| 17 | `splice_donor_region_variant` | `TX`, `EXON` | Donor +3..+6 rule |
| 18 | `splice_polypyrimidine_tract_variant` | `TX`, `EXON` | Acceptor polypyrimidine window |
| 19 | `incomplete_terminal_codon_variant` | `TX`, `EXON`, `SEQ` | Incomplete terminal codon logic |
| 20 | `start_retained_variant` | `TX`, `EXON`, `SEQ` | Start codon preserved |
| 21 | `stop_retained_variant` | `TX`, `EXON`, `SEQ` | Stop codon preserved |
| 22 | `synonymous_variant` | `TX`, `EXON`, `SEQ` | AA unchanged |
| 23 | `coding_sequence_variant` | `TX`, `EXON`, `SEQ` | Generic coding fallback |
| 24 | `mature_miRNA_variant` | `MIRNA` | Mature miRNA coordinates/arm mapping |
| 25 | `5_prime_UTR_variant` | `TX`, `EXON` | UTR localization, strand-aware |
| 26 | `3_prime_UTR_variant` | `TX`, `EXON` | UTR localization, strand-aware |
| 27 | `non_coding_transcript_exon_variant` | `TX`, `EXON` | Exonic non-coding transcript |
| 28 | `intron_variant` | `TX`, `EXON` | Intronic non-splice fallback |
| 29 | `NMD_transcript_variant` | `TX` | Biotype/annotation flag-driven |
| 30 | `non_coding_transcript_variant` | `TX` | Biotype-driven |
| 31 | `coding_transcript_variant` | `TX` | Coding transcript generic fallback |
| 32 | `upstream_gene_variant` | `TX` | Distance window (strand-aware) |
| 33 | `downstream_gene_variant` | `TX` | Distance window (strand-aware) |
| 34 | `TFBS_ablation` | `REG`, `MOTIF`, `SV` | Structural TFBS removal |
| 35 | `TFBS_amplification` | `REG`, `MOTIF`, `SV` | Structural TFBS amplification |
| 36 | `TF_binding_site_variant` | `REG`, `MOTIF` | TF motif overlap |
| 37 | `regulatory_region_ablation` | `REG`, `SV` | Structural removal of regulatory region |
| 38 | `regulatory_region_amplification` | `REG`, `SV` | Structural amplification of regulatory region |
| 39 | `regulatory_region_variant` | `REG` | Region overlap |
| 40 | `intergenic_variant` | `TX` | No transcript/regulatory hit fallback |
| 41 | `sequence_variant` | (fallback) | Base fallback when no specific term applies |

## Porting Sequence

1. `annotate_vep` API and backend abstraction (phase 1, scaffold)
2. Shared cache readers for `TX/EXON/SEQ/VAR` (SNV + small indel consequences)
3. Most-severe ranking + CSQ assembly + SQL-visible output parity tests
4. Regulatory + motif + miRNA consequences
5. Structural-event consequences (`SV`-dependent terms)
6. Performance tuning (projection pushdown, batching, fjall key-locality)

## Progress Status (Updated March 4, 2026)

| Step | Status | Notes |
|---|---|---|
| 1. `annotate_vep` API/backend abstraction | Complete | Unified `annotate_vep(vcf_table, cache_source, backend[, options_json])`; parquet/fjall cache source resolution implemented |
| 2. Shared readers + transcript/exon baseline | In progress | `lookup_variants` integration complete; transcript/exon-driven consequence engine wired when context tables are provided (`<cache>_transcripts`, `<cache>_exons`, or `options_json` overrides) |
| 3. Most-severe ranking + CSQ assembly + SQL tests | In progress | SO rank model (41 terms) implemented; `most_severe_consequence` computed from ranked terms in transcript/exon mode; fallback cache-hit CSQ remains for no-context mode |
| 4. Regulatory/motif/miRNA consequences | Not started | Artifacts and matrix defined; runtime evaluation not yet merged |
| 5. Structural-event consequences | Not started | `SV`-dependent terms (`ablation`, `amplification`, truncation/elongation) not yet merged |
| 6. Performance tuning | Not started | Need profiling and optimizations after consequence coverage is broader |

### Remaining Gaps (Next Priorities)

1. Add translation/codon-level effects for `missense`, `synonymous`, `stop_gained`, `stop_retained`, `start_retained`, and incomplete terminal codon logic.
2. Add regulatory/motif/miRNA table readers and consequence assignment path.
3. Add structural-event interpretation using `sv_features`.
4. Expand golden parity tests to include transcript/regulatory/SV representative cases and assert deterministic CSQ ordering.

## Acceptance Criteria

| Area | Required Check |
|---|---|
| SQL API | `annotate_vep(vcf_table, cache_source, backend[, options_json])` works for both backends |
| Consequence parity | Golden tests against VEP reference outputs for representative SNV/indel/regulatory variants |
| SO completeness | Matrix above covered (or explicit, tested unsupported terms) |
| Determinism | Stable `most_severe_consequence` ranking across runs/backends |
| Regression safety | Existing `lookup_variants` behavior unchanged |
