# VEP Plugins → vepyr: Systematic Port Feasibility & Difficulty Analysis

**Date:** 2026-07-10
**Status:** Analysis / brainstorming output (pre-decision)
**Repos in scope:** `datafusion-bio-functions` (plugin framework, primary), `datafusion-bio-formats` (source providers), `vepyr` (Python surface), `biodatageeks/vepyr-plugins` (manifest catalog)
**Excluded from this analysis (owned by someone else):** CADD, AlphaMissense, SpliceAI, dbNSFP (+ ClinVar, which is not a `.pm` plugin — handled as a `--custom` VCF preset)

> Method note: this is a brainstorming/feasibility pass, not a plan. The
> "obra superpowers brainstorming skill" is not installed in this environment;
> the analysis instead applies the same decompose → classify → adversarially-
> pressure-test methodology and lands the output as a `superpowers/specs` doc,
> matching the repo's existing workflow.

---

## 0. TL;DR verdict

**Porting the VEP plugin catalogue to vepyr is feasible but sharply bimodal in difficulty.** It is *not* one project — it is **five projects wearing one name**, and their difficulty differs by ~2 orders of magnitude.

- **~30 of the 78 in-scope plugins (the per-variant / per-transcript *scoring* class) are already essentially "done by construction."** The `plugin_cache` subsystem was built exactly for them. Porting one = **authoring a TOML source manifest** in `vepyr-plugins` + acquiring the source data + running the parity gate. **No Rust.** This is a genuine assembly line and is the systematic, high-ROI core of the effort.
- **The remaining ~48 plugins are the hard tail** and are *not* systematizable with today's framework. They split into three groups that each need a **new subsystem** before any of their members can be ported:
  1. **Gene / transcript-ID-keyed** plugins (LOEUF, pLI, LoFtool, GO, G2P, …) — need a gene-keyed lookup mechanism (explicit non-goal today).
  2. **Algorithmic / compute-from-sequence** plugins (Downstream, NMD, UTRAnnotator, MaxEntScan, …) — need per-plugin Rust reimplementation *inside the annotation engine*, plus exact parity against a Perl algorithm.
  3. **Online API / DB-service** plugins (LD, ProtVar, GeneBe, LOVD, …) — structurally incompatible with vepyr's offline-by-design contract; out of scope unless reframed as snapshot caches.
- **Testing is the cheap part for the easy class and the expensive part for the hard class.** The §8 parity gate (golden diff of the plugin's CSQ fields vs real Ensembl VEP, chr-scoped) is already defined and is fully templatable for data-lookup plugins. Algorithmic plugins need bespoke edge-case corpora and are where correctness risk concentrates.

**Bottom line:** commit to the ~30 data-lookup plugins as a repeatable factory (weeks, not months, of *engineering*; the real cost is data hosting + parity data generation). Treat the gene-keyed subsystem as one well-scoped follow-up that unlocks ~12 more cheaply. Treat the algorithmic and online classes as individually-scoped features, not a "port everything" sweep — several (UTRAnnotator, MaxEntScan, NMD) are multi-week engine projects on their own.

---

## 1. What already exists (the thing that makes this feasible)

The engine already ships a **complete, declarative plugin subsystem** — `datafusion/bio-function-vep/src/plugin_cache/` (~3.3k LOC across 14 files), plus:

- A **source-manifest format** (`source_manifest.rs`): per-plugin TOML declaring provider(s) (`vcf`/`csv`/`tsv`/`parquet`/`bed`), input schema, an `ingest_sql` SELECT that maps raw columns → the shared key `(chrom, start, end, allele_string)` + value columns, coordinate system, and value→CSQ-field mapping.
- A **per-transcript match discriminator** (`[[match_column]]` + `template`, e.g. `"{ref_aa}{Protein_position}{alt_aa}"`) that reproduces VEP's missense-gated, amino-acid-matched scoring (the AlphaMissense/dbNSFP/REVEL pattern) — a *no-match emits VEP's empty value*, so missense gating "falls out for free."
- A **build pipeline** (`build.rs`, `builder.rs`, `join.rs`): normalize contig/coordinate → LEFT-join the variation cache for allele frequency → derive warm/cold tier → write per-chrom, PageDir-indexed Parquet shards that reuse the variation cache's tuned point-lookup writer properties verbatim.
- A **runtime lookup** (`lookup.rs`, `registry.rs`, `provider.rs`): buffer-batched, page-scoped point lookup per plugin, CSQ injection per transcript consequence — the variation-lookup model cloned, variation path untouched.
- A **Python surface**: `vepyr.build_plugin_cache(plugin, version, …)` already resolves `plugins/<plugin>/<plugin>.source.toml` from `biodatageeks/vepyr-plugins` at a pinned tag; `annotate(..., plugin_cache_root=…)` wires the built caches into output.
- A **defined parity strategy** (design §8): unit tests for the freq-join tiering + contig/coord normalization, a synthetic round-trip test, and the **decisive gate** — chr-scoped golden diff of the plugin's CSQ fields against real Ensembl VEP.

**Consequence:** for the class the framework targets, "porting a plugin" is a *data + config + test* task, not a *code* task. That is the single most important finding in this analysis.

### The framework's own declared scope (design §9)

- **Target class:** per-variant scoring — tabular VCF/TSV/BED keyed by `chrom/pos/ref/alt`. **Directly supported.**
- **Non-goal:** per-gene / per-transcript-ID plugins (key by ID, not position) — "would need a separate gene-keyed mechanism."
- **Secondary:** per-interval / overlap plugins (conservation, tracks) — structurally supported by `start`/`end`, but frequency-tiering degrades to the interval start.
- **Escape hatch:** a sibling-crate Rust normalizer referenced by name, only when `ingest_sql` cannot express the transform.

This analysis inherits that scope boundary as the primary axis of difficulty.

---

## 2. Classification of all 78 in-scope plugins

Every plugin is assigned to a **bucket** by its *data/compute mechanism*, which is what determines porting cost. Difficulty is **Trivial / Low / Medium / High / Very High**. "New infra" flags a plugin that cannot be ported until a not-yet-built subsystem exists.

> Confidence: bucket assignment is high-confidence for the well-known scoring/
> constraint/algorithmic plugins and mechanism-inferred for a handful of niche
> ones (marked ⚠). Confirming each ⚠ by reading its `.pm` header is itself the
> first checklist item of the per-plugin playbook (§4), so no classification
> here is load-bearing without that cheap confirmation.

### Bucket A — Per-variant scoring, data-file lookup → **framework READY, port = manifest only. LOW.**

Keyed by `chrom/pos/ref/alt`; one score row per variant. Pure `ingest_sql`, no `match_column`.

| Plugin | Source format | Difficulty | Notes |
|---|---|---|---|
| REVEL | tabix TSV | Low | Highest-value missense score; often used per-transcript → may want `match_column` (Bucket B) |
| BayesDel | tabix TSV | Low | |
| CAPICE | tabix TSV | Low | |
| ClinPred | tabix TSV | Low | |
| EVE | tabix VCF | Low | |
| Enformer | tabix | Low | |
| FATHMM_MKL | tabix TSV | Low | |
| Geno2MP | tabix VCF | Low | phenotype counts |
| Mastermind | tabix VCF | Low | literature citations |
| MaveDB | tabix | Low | multiplex assay scores |
| gnomADMt | tabix VCF | Low | mitochondrial AF |
| satMutMPRA | tabix | Low | saturation-mutagenesis MPRA |
| SpliceVault | tabix | Low | |
| mutfunc ⚠ | tabix | Low | |
| AVADA ⚠ | tabix | Low | literature-derived variants |
| SubsetVCF | VCF | Low | annotate from an arbitrary VCF (≈ `--custom`) |
| FATHMM ⚠ | tabix/DB | Low–Med | some versions call a service (→ Bucket F) |

### Bucket B — Per-transcript / missense-gated scoring → **framework READY (`match_column`). LOW–MEDIUM.**

Score depends on the amino-acid change / transcript; uses the `[[match_column]]` template.

| Plugin | Source | Difficulty | Notes |
|---|---|---|---|
| PrimateAI | tabix TSV | Low–Med | aa-matched |
| MPC | tabix | Low–Med | missense |
| VARITY | tabix | Low–Med | missense |
| TranscriptAnnotator ⚠ | tabix | Medium | generic per-transcript annotator |
| Carol / Condel | derived | Low–Med | combine SIFT+PolyPhen — *engine-attribute*, may not need a data file (borderline Bucket E, but trivial math) |

### Bucket C — Position / interval tracks (no allele, or region-scored) → **structurally supported; tiering degrades; may need a provider. MEDIUM.**

| Plugin | Source | Difficulty | Notes |
|---|---|---|---|
| Conservation | bigwig/GERP | Medium | position track; needs bigwig provider **or** pre-converted TSV |
| gnomADc | tabix/bigwig | Low–Med | constraint/coverage track |
| MTR | tabix | Low–Med | missense tolerance ratio (regional) |
| dbscSNV | tabix TSV | Low | splice score, position-keyed |
| PromoterAI ⚠ | tabix | Low–Med | |
| RiboseqORFs ⚠ | tabix/GFF | Medium | ORF overlap |
| FunMotifs ⚠ | tabix | Low–Med | |
| ReferenceQuality ⚠ | bigwig/region | Medium | needs region/track provider |
| gnomADMt (dup) | — | — | listed in A |
| AncestralAllele | FASTA | Medium | **new provider**: per-position base from an ancestral-genome FASTA |
| FlagLRG | region map | Low–Med | overlap flag |

### Bucket D — Gene / transcript-ID keyed → **NEW gene-keyed subsystem required. HIGH (subsystem) → then LOW–MED each.**

Not position-keyed; the variation freq-join and PageDir point-lookup **do not apply**. Explicit non-goal today. All of these unlock together once a gene-ID-keyed lookup exists.

| Plugin | Keyed by | Difficulty | Notes |
|---|---|---|---|
| LOEUF | gene | High → Med | gnomAD constraint |
| pLI | gene | High → Med | constraint |
| LoFtool | gene | High → Med | per-gene percentile |
| DosageSensitivity | gene | High → Med | haploinsufficiency |
| GO | gene/transcript | High → Med | Gene Ontology terms |
| G2P | gene panel | High | + panel/inheritance logic beyond a lookup |
| MechPredict ⚠ | gene | High → Med | |
| PhenotypeOrthologous ⚠ | gene/ortholog | High → Med | |
| Phenotypes | gene/region | Medium | GFF of phenotype features (partly Bucket C) |
| GWAS ⚠ | position/gene | Medium | catalog associations |
| IntAct ⚠ | variant/gene | Medium | interaction data (partly A) |
| OpenTargets ⚠ | target/disease | Medium | tabix or API (partly A/F) |
| neXtProt ⚠ | protein/gene | High | protein annotations, may be DB-backed (partly F) |

### Bucket E — Algorithmic compute (from transcript/sequence, no external score file) → **per-plugin Rust in the engine. TRIVIAL → VERY HIGH.**

Each reimplements a Perl algorithm. No manifest can express these; they hook the transcript-consequence engine. Wide internal spread:

| Plugin | Difficulty | Notes |
|---|---|---|
| SingleLetterAA | Trivial | 3-letter → 1-letter AA in HGVSp (string map) |
| Blosum62 | Trivial–Low | static substitution matrix keyed on aa change |
| TSSDistance | Low | distance to TSS (geometry; needs transcript coords) |
| NearestGene | Low–Med | overlaps VEP `--nearest`; already roadmap P2 |
| NearestExonJB | Medium | nearest exon junction (geometry) |
| SameCodon | Medium | flag variants sharing a codon |
| SpliceRegion | Medium | extended splice-region SO terms (geometry) |
| HGVSReferenceBase | Low–Med | reference base for HGVS |
| HGVSIntronOffset | Medium | intron-offset HGVS math |
| RefSeqHGVS | Med–High | RefSeq-based HGVS (needs RefSeq transcripts) |
| CSN | Med–High | Clinical Sequencing Nomenclature notation |
| ProteinSeqs | Medium | emit wild-type + mutant protein sequences |
| Downstream | High | predict downstream protein after frameshift (sequence synthesis) |
| DeNovo | Medium | de-novo flag from trio genotypes (needs sample/GT plumbing) |
| NMD | High | nonsense-mediated-decay prediction (transcript-structure rules) |
| GeneSplicer | High | wraps the GeneSplicer binary; reimplement or shell out |
| MaxEntScan | High | MaxEnt splice model + score matrices (algorithm + data) |
| UTRAnnotator | Very High | 5′UTR uORF/uAUG effect engine — a mini-annotator itself |
| PolyPhen_SIFT | Low / skip | already native in vepyr core (redundant) |

### Bucket F — Online API / DB service → **incompatible with offline design. OUT OF SCOPE (or snapshot-cache reframing). N/A–HIGH.**

| Plugin | Mechanism | Disposition |
|---|---|---|
| LD | Ensembl DB / REST | Out of scope; needs LD snapshot |
| ProtVar | REST API | Out of scope |
| GeneBe | REST API (ACMG) | Out of scope |
| LOVD | REST API | Out of scope |
| PON_P2 | web tool | Out of scope |
| DAS | Distributed Annotation System | Deprecated; skip |
| GXA | Expression Atlas (live/gene) | Out of scope (partly D) |
| Paralogues | Compara DB + data | High; needs paralogue snapshot (partly D) |

### Bucket G — Output filter / modifier / utility → **needs an output-hook (not built). LOW–MEDIUM.**

These don't add scores; they drop/flag/format output rows. vepyr's `pick`/`filter` machinery is the closest existing seam.

| Plugin | Difficulty | Notes |
|---|---|---|
| NonSynonymousFilter | Low | drop non-missense lines |
| CCDSFilter | Low | restrict to CCDS transcripts |
| RankFilter | Low | rank-threshold filter |
| LocalID | Low | inject a local ID column |
| Draw | Medium | transcript diagram (visualization; niche) |

### Count by bucket (78 total)

| Bucket | ~Count | Framework status | Dominant cost |
|---|---|---|---|
| A — per-variant lookup | ~17 | **Ready** | data hosting + parity |
| B — per-transcript lookup | ~5 | **Ready** | data + aa-match template |
| C — position/interval track | ~10 | Mostly ready (± provider) | tiering caveat, a provider or two |
| D — gene/ID-keyed | ~13 | **New subsystem** | one-time gene-keyed lookup, then cheap |
| E — algorithmic | ~19 | **Per-plugin engine code** | Rust + exact Perl parity |
| F — online/DB | ~8 | Out of scope | reframe as snapshot or skip |
| G — output filter | ~5 | Needs output hook | modest |

(Buckets overlap at the margins — REVEL A↔B, OpenTargets A↔D↔F — so counts are ±2.)

---

## 3. Difficulty drivers (what actually costs time)

For the *systematic* class (A/B, and C once a provider exists), engineering is **not** the bottleneck. The real costs, in order:

1. **Source-data acquisition & hosting.** Each plugin's data is a separate multi-GB download with its own license, versioning, and coordinate quirks (`chr1` vs `1`, 0- vs 1-based, INFO-array packing, multi-allelic split). The framework already standardizes contig/coordinate and gives an `ingest_sql` + `unnest` escape hatch, but *sourcing and license-clearing each file* is real, per-plugin, non-code work.
2. **Parity-data generation.** The decisive gate needs **real Ensembl VEP output for the same plugin** on a test region — i.e. you must run VEP-with-the-Perl-plugin once per plugin to get golden output. That means installing each plugin's Perl deps + data at least once. This is the hidden long pole of "systematic testing."
3. **The aa-change/missense edge cases** for Bucket B (indels, MNVs, stop-gain, selenocysteine, multi-isoform genes) — the `match_column` template handles the happy path; parity failures cluster here.

For the hard class (D/E/F), engineering *is* the bottleneck and dominates everything above.

---

## 4. The systematic porting playbook (Bucket A/B/C)

This is the part that *is* an assembly line. Per plugin:

1. **Confirm mechanism** — read the `.pm` header; verify it's a positional data-lookup (resolves any ⚠). ~15 min. Kick to D/E/F if not.
2. **Acquire + stage source data** for the test chrom(s).
3. **Author `plugins/<name>/<name>.source.toml`** in `vepyr-plugins`: provider + params + schema, `ingest_sql` (col mapping only — the builder owns contig/coord), value→CSQ mapping, `[[match_column]]` if per-transcript, tier policy.
4. **Build** the per-chrom cache via `vepyr.build_plugin_cache(...)`.
5. **Parity-gate:** run real Ensembl VEP with the plugin on the test region → run vepyr with the cache → diff the plugin's CSQ fields → require 100% (plus w1-vs-w4 byte-identity).
6. **Contribute** the manifest (+ golden fixture) as a PR to `vepyr-plugins`. No change to the engine repo unless a *new format/provider* is needed.

**Templatable artifacts to build once** (the real force-multiplier, ~a few days):
- A **manifest scaffolder** (`vep-add-plugin` skill, already designed in §7 of the plugin-cache design) — stub a `.source.toml` from a couple of prompts.
- A **generic parity harness** parameterized by `(plugin, csq_fields, test_region)` — one script runs VEP+plugin, runs vepyr, diffs. Reused verbatim for every A/B/C plugin. This is what turns "test each plugin" from N bespoke efforts into N config rows.
- A **provider gap list** (bigwig, ancestral-FASTA) — a handful of Bucket-C plugins share these; add each once in `datafusion-bio-formats` and several plugins unlock.

**Rough effort (framework unchanged):** ~0.5–2 engineer-days per A/B plugin, *dominated by data + parity, not code*. ~30 plugins → a **factory measured in a few engineer-weeks of coding + a larger, parallelizable tail of data/license/parity work.**

---

## 5. The non-systematic tail (Bucket D/E/F/G)

- **Bucket D (gene-keyed), ~13 plugins.** One-time subsystem: a gene/transcript-ID-keyed lookup (a small keyed table joined on the engine's gene/transcript ID, emitted per consequence). Estimate **~1–2 engineer-weeks** for the mechanism; then each plugin is a **Low–Medium** manifest-like add (~1–3 days). Best ROI in the tail: one build unlocks LOEUF, pLI, LoFtool, DosageSensitivity, GO, and more.
- **Bucket E (algorithmic), ~19 plugins.** No shared unlock — each is its own feature. Internal spread is enormous: SingleLetterAA/Blosum62/TSSDistance are **hours–days**; Downstream/NMD/GeneSplicer/MaxEntScan are **1–4 weeks each**; **UTRAnnotator is a multi-week engine sub-project.** Parity is the risk: matching a Perl algorithm bit-for-bit across indels, strands, and incomplete-CDS edge cases. Recommend cherry-picking the high-value, low-cost ones (SpliceRegion, NearestGene, ProteinSeqs, the HGVS helpers) and treating the heavyweights as individually-justified roadmap items.
- **Bucket F (online), ~8 plugins.** Structurally against vepyr's offline-by-design contract. **Do not port as-is.** The only in-spirit path is reframing a service as a *snapshot cache* (e.g. an LD snapshot, a ProtVar dump) — which turns it back into a Bucket A/D data problem, at the cost of freshness. Otherwise: out of scope.
- **Bucket G (filters), ~5 plugins.** Need an output-row filter/modifier hook. vepyr's existing `pick`/`filter`/`most_severe` machinery is the natural home; adding filter predicates there is **Low–Medium** and shared across all five.

---

## 6. Systematic testing strategy

The framework already names the gate (design §8); systematizing it means building the harness *once* and running it as data:

1. **Tier 0 — unit (Rust, per mechanism, not per plugin):** freq-join tiering table-tests, `canonical_contig` + coordinate-shift tests, synthetic shard round-trip (hit/miss/tier). Already the pattern in-repo.
2. **Tier 1 — golden parity (per plugin, the decisive gate):** the generic `(plugin, csq_fields, region)` harness diffing vepyr vs real Ensembl VEP. **Reused verbatim** for all A/B/C plugins → testing is O(config), not O(code). This is the single most important test investment.
3. **Tier 2 — byte-identity across workers:** w1 vs w4 body-identical (existing parallel-annotation invariant) — free, catches concurrency bugs in the lookup.
4. **Tier 3 — algorithmic corpora (Bucket E only):** hand-built edge-case VCFs per algorithm (frameshift near stop for Downstream, uORF variants for UTRAnnotator, splice-boundary indels for MaxEntScan). This is bespoke and is where test cost is real.

**Hidden dependency (call it out loud):** Tier 1 requires a working Ensembl VEP + each plugin's Perl deps + data installed *once* to generate golden output. Standing up that reference environment is a prerequisite line item, not an afterthought — and it is the true gating cost of "systematic testing," more than writing any diff code.

---

## 7. Recommended sequencing

1. **Wave 0 — factory (a few days).** Build the `vep-add-plugin` scaffolder + the generic parity harness + stand up the reference VEP environment. Nothing ports faster than this pays back.
2. **Wave 1 — Bucket A/B high-value scoring (~20 plugins).** REVEL, PrimateAI, MPC, VARITY, BayesDel, CAPICE, ClinPred, EVE, FATHMM_MKL, Mastermind, Geno2MP, MaveDB, SpliceVault, dbscSNV, gnomADMt, satMutMPRA, Enformer, … Pure manifest work on the factory. Parallelizable across people; each is independent.
3. **Wave 2 — Bucket C providers + tracks (~8 plugins).** Add bigwig + ancestral-FASTA providers once; unlock Conservation, gnomADc, MTR, AncestralAllele, ReferenceQuality, etc.
4. **Wave 3 — Bucket D subsystem (~13 plugins).** Build gene-ID-keyed lookup once; then LOEUF, pLI, LoFtool, DosageSensitivity, GO, … as cheap adds.
5. **Wave 4 — Bucket E cherry-picks + Bucket G filters.** SingleLetterAA, Blosum62, TSSDistance, SpliceRegion, NearestGene, ProteinSeqs, HGVS helpers; the five filter plugins via the output hook.
6. **Wave 5 (opt-in, per-item justified) — Bucket E heavyweights.** UTRAnnotator, MaxEntScan, NMD, GeneSplicer, Downstream — each a standalone project with its own parity corpus.
7. **Never (as-is) — Bucket F.** Only via snapshot reframing, item by item.

---

## 8. Risks & pressure-tests (adversarial)

- **"~30 easy plugins" assumes data availability & licensing.** Some sources are gated/registration-only or non-redistributable; the manifest catalog is public but the *data* isn't shippable. Mitigation: manifests reference user-provided paths; `vepyr-plugins` ships only TOML + fixtures (already the design). Cost is per-user download UX, not our code — but it *does* cap how "one command, zero setup" the story can be.
- **Parity may not hit 100% for allele-edge-cases** even in the easy class (multi-allelic, MNV, spanning-deletion, indel left-alignment). The `match_column` template covers the missense happy path; expect a real debugging tail here per plugin. This is why Tier-1 parity is mandatory, not optional.
- **Bucket boundaries are fuzzy.** REVEL/OpenTargets/FATHMM/neXtProt straddle buckets; a ⚠ plugin may reveal an API dependency (→ Bucket F) only when you read its header. The §4 step-1 header read is the cheap guard.
- **Gene-keyed and algorithmic estimates are order-of-magnitude, not commitments.** UTRAnnotator especially could balloon; scope it as its own spec before committing.
- **Reference-VEP environment is a real prerequisite.** No golden output ⇒ no parity gate ⇒ no trustworthy port. Budget it explicitly in Wave 0.

---

## 9. Answer to the question asked

- **Is a systematic port feasible?** *Yes for the class the framework was built for (~30 plugins), where it is a genuine, low-risk assembly line.* Partially for another ~13 behind one well-scoped subsystem. Selectively for the algorithmic tail. No, as-is, for the ~8 online plugins.
- **How hard?** Bimodal. The systematic core is **LOW** difficulty and **high** throughput once the Wave-0 factory exists — the cost is data + parity ops, not engineering. The tail is where difficulty lives: **one MEDIUM subsystem** (gene-keyed) unlocks a dozen cheaply, but **a handful of HIGH/VERY-HIGH algorithmic plugins** (UTRAnnotator, MaxEntScan, NMD, GeneSplicer, Downstream) are multi-week projects each and should never be lumped into a "port them all" estimate.
- **The one-line strategy:** *Build the factory, run the ~30-plugin assembly line, invest once in the gene-keyed subsystem, and treat every algorithmic/online plugin as an individually-justified feature — never as a line item in a bulk port.*
