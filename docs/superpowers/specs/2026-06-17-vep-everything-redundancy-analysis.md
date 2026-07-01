# VEP `--everything` flow — redundancy analysis vs Ensembl VEP (2026-06-17)

Companion to `docs/superpowers/plans/2026-06-17-vep-engine-perf-handoff.md`.
Goal: gain **confidence** that our transcript-consequence engine does not do redundant work
that VEP avoids — i.e. distinguish *genuine reclaimable redundancy* from *work VEP also does*
(which is therefore correctness-load, not waste). Every claim is backed by file:line on **both**
sides. Conservative throughout: per-call costs in our engine are already 1–2 µs, so impact is
bounded by **volume**, not per-call savings.

## Reference checkout (evidence base)
The actual consequence + HGVS engine is **not** in `ensembl-vep` (the driver) — it lives in two
sibling Perl APIs that VEP depends on. All three are checked out at `release/115(.2)`:

| Repo | Path | Role |
|---|---|---|
| `ensembl-vep` | `/Users/mwiewior/research/git/ensembl-vep` (rel 115.2) | driver: option sets, OutputFactory, pick, colocated/AF assembly |
| `ensembl-variation` | `/Users/mwiewior/research/git/ensembl-variation` (rel 115) | `TranscriptVariation(Allele)`, HGVS, consequence predicates |
| `ensembl` core | `/Users/mwiewior/research/git/ensembl` (rel 115) | `Transcript`, `TranscriptMapper` (coord table) |

## What `--everything` actually turns on
`ensembl-vep/modules/Bio/EnsEMBL/VEP/Config.pm:348-374` — the `everything` OPTION_SET sets:
`sift=b, polyphen=b, ccds, hgvs(→hgvsc+hgvsp at :378-383), symbol, numbers, domains, regulatory,
canonical, protein, biotype, af, af_1kg, af_gnomade, af_gnomadg, max_af, pubmed, uniprot, mane,
tsl, appris, variant_class, gene_phenotype, mirna`.
**Critically: it sets NO `pick` / `most_severe` / `flag_pick*`.** So `--everything` =
*full HGVSc/HGVSp for every overlapping transcript, no filtering* — exactly the multi-isoform,
all-transcript mode our e2e parity gate validates (4.74M rows, 0 mismatches). This frames the whole
analysis: **under the profiled mode VEP itself does NOT skip transcripts.**

---

## The headline structural difference (root cause of C, D, and #2)

VEP and our engine have **inverted loop nesting**:

```
VEP (ensembl-vep AnnotationType/Transcript.pm:101-147)
  for each TRANSCRIPT in region buffer:        # outer  (:101)
      build/reuse geometry ONCE on the Transcript object
      for each VARIANT overlapping it:         # inner  (:121)
          TranscriptVariation->new(-transcript=>$tr, -variation_feature=>$vf)   (:135)

OURS (annotate_provider.rs:5875 → transcript_consequence.rs:1149)
  for each ROW = (position, single ALT):       # outer  (annotate_provider.rs:5875)
      COITree.query(...)                        # transcript_consequence.rs:1149
      for each TRANSCRIPT hit:                  # inner
          rebuild exon sort / coord table / cds bounds  PER (variant×transcript)
```

Because VEP's transcript is the **outer** loop, every piece of per-transcript geometry is built by
the *first* overlapping variant and reused by *all the rest* (and all their alleles). Our engine,
with the variant outer, re-derives that geometry on every (variant × transcript) candidate — 4.66M
times for chr1. This single inversion is the source of redundancies **C**, **D**, and the per-allele
issue **#2**. It is *not* a correctness difference (parity gate is green) — it is pure recompute.

---

## Redundancy-by-redundancy verdict

### A. Materialize-then-pick — **0 under `--everything`; ~3–5 s only if `--pick` is used**
- **VEP picks BEFORE materializing.** `OutputFactory.pm`: `get_all_VariationFeatureOverlapAlleles`
  calls `filter_VariationFeatureOverlapAlleles` (`:513`); under `--pick` that returns a single
  survivor `[pick_worst_VariationFeatureOverlapAllele($vfoas)]` (`:584-586`) **before** the
  per-allele `..._to_output_hash` loop (`:388-403`) that reads `hgvs_transcript`/`hgvs_protein`
  (`:1699`, `:1708`). HGVS is additionally **lazy + memoized** (`TranscriptVariationAllele.pm:1320`
  hgvs_transcript guard, `:1614` hgvs_protein guard) — losers never trigger it.
- **Under `--everything` (no pick)** `filter_…` falls through every pick branch and returns the full
  list (`OutputFactory.pm:584-625` → unfiltered `return`), so VEP *also* materializes HGVS for all
  ~12 transcripts. **We match VEP here.**
- **Our engine** computes HGVSc (`transcript_consequence.rs:1368`) + HGVSp (`:1409`) and `push`es
  (`:1435`) unconditionally; `picked` (`:587`) is set later in `collapse_pick_sort`
  (`annotate_provider.rs:2569` `pick_worst_assignment`) — i.e. **materialize-then-pick**.
- **Verdict:** redundancy A is **real but conditional**. In the profiled all-transcript mode it
  reclaims **0** (VEP does the same work). It only bites if the *product* runs `--pick`/most-severe,
  where it would reclaim most of hgvsc 5.49 s + hgvsp 0.18 s. **This is a semantic/product decision
  (open question #1), not a bug.** If pick is ever the target mode, restructure to classify+rank →
  pick → materialize HGVS only for the winner, matching VEP's order.

### B. HGVSc on non-cDNA outputs — **genuine, but small; bounded by the up/down/intergenic share of `coord_map`**
- **VEP hard-gates** the whole HGVSc/cDNA/protein block behind `if($pre->{within_feature})`
  (`OutputFactory.pm:1669`, closes `:1718`); `hgvs_transcript` (`:1699`) is *never reached* for
  upstream/downstream (consequences flagged `within_feature => 0` at `Config.pm:914,931`) or
  intergenic (not a `TranscriptVariationAllele` at all). The output write is itself guarded
  `$hash->{HGVSc} = $hgvs_t if $hgvs_t` (`:1702`).
- **Our engine** gates the *overlap* materialization at `transcript_consequence.rs:1194`
  (`if variant_overlaps_tx`), but the handoff's measured `hgvsc_calls == tx_outputs == 3.85M`
  indicates `format_hgvsc` (`hgvs.rs`) still runs on outputs with no cDNA position and pays
  `coord_map` before returning empty. VEP avoids that cost entirely via the pre-gate.
- **Verdict:** real redundancy = the `coord_map` time spent on up/down/intergenic outputs that can
  never produce HGVSc. **Conservative impact: a slice of the 4.14 s `coord_map`** (size = fraction of
  the 3.85M outputs that are non-cDNA). Fix = a VEP-style `within_feature`/`has-cDNA-position`
  short-circuit *before* `coord_map`. **First confirm the counter** — re-profile with
  `hgvsc_detail_profile` and read how many `format_hgvsc` calls return `None`.

### C. Per-transcript geometry recomputed per (variant × transcript) — **the biggest genuine, behavior-preserving win**
- **VEP builds each of these exactly ONCE per transcript**, hung on the Transcript object's
  `{_variation_effect_feature_cache}` (reused by every overlapping variant *and* allele):
  - cdna↔genomic↔pep **mapper / coord table**: `BaseTranscriptVariation.pm:1079`
    `{mapper} ||= get_TranscriptMapper`; mapper itself memoized `Transcript.pm:2270`
    `{transcript_mapper} ||=`, table built once in ctor `TranscriptMapper.pm:96-159`.
  - **sorted exons / introns / interval trees**: `BaseTranscriptVariation.pm:796` `{sorted_exons}||=`,
    `:775` introns, `:1038-1067` exon/intron interval trees; exon array `Transcript.pm:1457-1463`.
  - **CDS bounds**: `cdna_coding_start/end` memoized `Transcript.pm:946-995`/`1016-1062`,
    captured into the mapper at `TranscriptMapper.pm:120-121`.
- **Our engine recomputes per call:** `add_intron_splice_terms` does `tx_exons.to_vec(); sort_by_key`
  every call (`transcript_consequence.rs:3131-3132`); `exon_segments` rebuilds + `genomic_to_cdna_index`
  **linear-scans** the segment list (`:7270-7299`); `which_exon_str`/`which_intron_str`/
  `compute_cdna_position` (`:1218-1220`); `coding_cdna_bounds` recomputed inside
  `shift_to_hgvs_coding_coordinates` (`hgvs.rs:1505-1546`). Static *scalar* fields (`cds_start/end`,
  `cdna_coding_start/end`, `cdna_mapper_segments`) **are** already cached on `TranscriptFeature`
  (`transcript_consequence.rs:194-286`) — but the *derivation paths above don't use them*; they
  re-derive from raw exons.
- **Verdict:** **genuine redundancy, fully behavior-preserving** (parity gate proves equality).
  Fix = build a sorted `Vec<TranscriptCdnaCoord>` + cds bounds + sortedness flag **once per
  transcript** (lazily, on first candidate), cache on the per-batch `PreparedContext`, and
  binary-search instead of linear scan — mirroring VEP's `{_variation_effect_feature_cache}`.
  **Conservative impact: ~2–4 s** (bounded because per-call is 1–2 µs over 4.66M candidates; the win
  is removing the rebuild+sort+linear-scan, helping both the 76k fast path and 11k fallback).

  **STATUS (2026-06-17): implemented + unit-verified; e2e parity & profile PENDING.**
  Implemented on branch `vep-transcript-geometry-cache` (commits `8a7d4ca`, `25cd70d`) — but with a
  *different, lower-risk mechanism* than the original plan threaded a param through ~59 call sites.
  Instead, a lazily-populated `OnceLock` cache (`GeometryCache` newtype) lives on
  `TranscriptFeature` so every `&TranscriptFeature` leaf function — and thus the engine hot path —
  hits it with **no parameter threading** (~10 edits). `transcript_cdna_coords` memoizes via the
  field; `genomic_to_cdna_index_for_transcript` reads the cached table
  (`cdna_index_from_sorted_coords`) instead of rebuilding `exon_segments`. This mirrors VEP's
  `{_variation_effect_feature_cache}` first-touch semantics even more closely than the plan.
  - **Verified (local):** clippy clean; new unit tests assert `cached == uncached` and that the
    cached lookup matches the linear `genomic_to_cdna_index` for both strands; full lib suite
    **880 passed / 31 failed**, where the 31 are pre-existing missing-fixture/env failures
    **identical on baseline** (0 new failures). 197 hgvs/cdna/parity unit tests exercise the
    changed paths.
  - **NOT yet verified (blocked):** the authoritative e2e parity gate (4.74M-row compare) and the
    profile delta vs the ~2–4 s estimate. The `vepyr` e2e/profiling toolchain runs a **prebuilt
    `polars-bio` wheel git-pinned to a fixed rev of this crate**, so it does not contain the local
    change. Running it requires rebuilding the `polars-bio` extension against this working tree
    (`maturin develop` with a path-override) and reinstalling into vepyr's env — a heavy multi-repo
    build. Until that runs: **do not claim the ~2–4 s win as measured**, and do not merge past the
    parity gate (this touches HGVSc/cDNA-position correctness).

### D. Fallback string round-trip — **~0.3–0.8 s; plus a bigger why-do-we-fall-back question**
- **VEP keeps coords as integers** in the prebuilt mapper and projects via `genomic2cdna` (no
  string round-trip): `BaseTranscriptVariation.pm:483` `{_cdna_coords} ||= [...]`.
- **Our fallback** (2.0 s, 3 µs/call, 18% of hgvsc calls) round-trips through `String`:
  `notation_to_hgvsc_coords` → `(String,String)`, then `shift_to_hgvs_coding_coordinates` *re-parses*
  via `split_hgvs_coord` and *recomputes* `coding_cdna_bounds` (`hgvs.rs:1505-1546`).
- **Verdict:** real but modest. Fix = keep coords as integers + intron-offset until final format.
  **Conservative impact: ~0.3–0.8 s.** Larger, cheaper lever: **why do 18% fall back?** Routing more
  of those to the ~1 µs fast path beats optimizing the fallback. (Subsumed by C — the shared coord
  table removes the recompute regardless of path.)

### E. Colocated / AF field assembly — **mostly matches VEP; the genuine Rust-only waste is the string clones**
- **VEP also re-sorts and re-parses per output row** — this is *not* extra redundancy on our part:
  - sorts `@{$vf->{existing}}` every row (`OutputFactory.pm:1032-1041`, called per-vfoa at `:1383`);
  - sorts the ~27 AF group keys every row (`:1161`);
  - re-parses each freq string (`"A:0.001,G:0.999"`) per row, rebuilding `%freq_data`/`%remaining`
    and re-deriving the minor allele (`:1171-1196`).
  - allele→freq matching is a **linear** `grep {$_->{a_allele} eq $this_allele}` (`:1045`, `:1154`).
  - VEP *does* memoize one thing once (`get_matched_variant_alleles` →
    `AnnotationType/Variation.pm:200-203`), but **not** the freq parsing or the two sorts.
- **Our engine** (`annotate_provider.rs`): `sorted_entries()` allocs+sorts a `Vec` twice per row
  (`:1923-1931`); `frequency_fields` builds a `HashMap`+`HashSet` per (entry × ~27 AF cols)
  (`:2071-2075`) and clones each freq 3–4× (`:2082-2124`).
- **Verdict:** the per-row sort + freq parse **match VEP's (already non-optimal) behavior** — so
  changing them is an *improvement over VEP*, not closing a gap. The **3–4× string clone has no VEP
  analog** (Perl uses refcounted COW strings) → that part is **pure Rust-side overhead** and the
  safest win. **Conservative impact: ~1–2 s**, of which the clone elimination is the
  highest-confidence slice. Hoistable beyond VEP: pre-sort entries + parse AF strings into a map
  **once** at `ColocatedData` build (cache-read), linear `&str` scan + format only the chosen allele.

### #2. Per-ALT-allele recompute — **real for multi-allelic sites; a SUBSET of C**
- **VEP does NOT recompute the genomic→cDNA mapping per ALT.** One `TranscriptVariation` per
  (transcript × VariationFeature) holds lightweight per-allele objects carrying only the allele
  sequence (`VariationFeatureOverlap.pm:412-473`); coords are memoized once on the TV
  (`BaseTranscriptVariation.pm:483` `{_cdna_coords} ||=`, `:142` `{cdna_start}`), shared by all
  alleles; per-allele work is only sequence-dependent (consequence prediction, codon/aa change,
  HGVS allele string).
- **Our engine** runs once per **row = (position, single ALT)** (`annotate_provider.rs:5875` row loop;
  multi-allelic sites are pre-split to one row per ALT — no inner allele loop). So an N-allele site
  re-runs the COITree query (`transcript_consequence.rs:1149`) + all per-transcript geometry N times.
- **Verdict:** genuine redundancy, but **bounded by the multi-allelic fraction** (small in human WGS).
  **Fixing C subsumes it**: a transcript-keyed geometry cache on the `PreparedContext` is reused
  across rows regardless of whether they are different positions or different ALTs at one position.
  No separate fix needed — note it as a free rider on C.

---

## Conservative impact summary

| # | Redundancy | VEP avoids it? | Genuine vs matches VEP | Conservative reclaim | Risk |
|---|---|---|---|---|---|
| C | per-(variant×transcript) geometry rebuild | **yes** (once/transcript) | **genuine** | **~2–4 s** | low (behavior-preserving) |
| E | AF/colocated per-row clones | partly (clones have no VEP analog) | clones **genuine**; sort/parse match VEP | **~1–2 s** (clones = safest slice) | low |
| B | HGVSc coord_map on non-cDNA outputs | **yes** (`within_feature` pre-gate) | **genuine** (confirm counter first) | slice of 4.14 s coord_map | low |
| D | fallback string round-trip | yes (integer coords) | genuine (subset of C) | ~0.3–0.8 s | low |
| #2 | per-ALT recompute | **yes** (shared TV coords) | genuine, **subset of C** | folded into C | low |
| A | materialize-then-pick | yes (pick before materialize) | **0 under `--everything`** | 0 now / ~3–5 s iff `--pick` | n/a (product decision) |

**Low-risk behavior-preserving set (C + E-clones + B + D):** realistically **~3–5 s** off the 25 s
engine — consistent with the handoff's independent estimate. **A is 0 in the profiled all-transcript
mode**; the order-of-magnitude lever remains *reducing transcripts-per-variant via pick mode*, which
is a semantic/product choice, not a redundancy fix.

## Confidence statement (the "conservative look" the task asked for)
- **High confidence, evidence-backed:** C, D, #2 are genuine recompute that VEP structurally avoids
  by caching geometry once per transcript and sharing the TranscriptVariation across alleles
  (`_variation_effect_feature_cache`, `{_cdna_coords}`). E's *string clones* are pure Rust overhead.
- **High confidence it is NOT redundancy:** the 12–14× transcript fan-out and per-transcript HGVS
  under `--everything` — VEP does the identical work (parity gate: 4.74M rows, 0 mismatches). Do not
  "optimize" by dropping transcripts unless the product switches to `--pick`.
- **Confirm before trusting the number:** B's exact size (how many `format_hgvsc` calls return
  `None`) and E's sort/parse hoist (matches VEP, so it's an over-VEP improvement, lower certainty on
  its share). Trust no figure until re-profiled with the `vep-perf-profiling` skill, single-thread,
  all levels, AND the e2e parity gate WITHOUT `--skip-compare` (C, B, D touch HGVSc correctness).

## Recommended order of attack (highest confidence first)
1. **C — transcript-keyed geometry cache on `PreparedContext`** (sorted exon/coord table + cds bounds,
   binary search). Subsumes #2 and the recompute half of D. Re-profile + parity gate.
2. **E — eliminate the 3–4× freq string clones** (Rust-only waste); then hoist the entry sort + AF
   parse to `ColocatedData` build (over-VEP improvement).
3. **B — `within_feature`/has-cDNA pre-gate before `coord_map`** once the counter confirms the share.
4. **D — integer coords through to final format** (largely free after C).
5. **A — only if the product mandates `--pick`** (then restructure to pick-before-materialize).
