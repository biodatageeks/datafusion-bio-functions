# VEP Plugin Port Factory (Wave 0) — Design

**Date:** 2026-07-13
**Status:** Design, approved in brainstorming; not yet implemented
**Supersedes the runbook in:** `HANDOFF veppluginsport.md` (2026-07-12) — see §1, several of its
premises are false against the code.
**Builds on:**
- `docs/superpowers/specs/2026-07-05-custom-vep-plugin-caches-design.md` (engine design)
- `docs/superpowers/specs/2026-07-06-vepyr-plugin-cache-integration-design.md`
- `docs/superpowers/specs/2026-07-10-vep-plugins-systematic-port-feasibility.md` (bucket analysis;
  on branch `dev-test-plugins-feasibility`)

---

## 0. Summary

The feasibility analysis concluded that ~30 of the 78 in-scope VEP plugins are a genuine assembly
line: porting one means authoring a TOML manifest, not writing Rust. That conclusion holds. What
does *not* exist yet is the assembly line itself — the thing that makes plugin #2 cost a day
instead of a week.

This spec designs **Wave 0: the factory**. Four components:

1. **A region-sliced mini-cache fixture** — because the real cache is 34 GB and CI cannot hold it.
2. **A parity harness** with a blame-attribution rule — the decisive gate, reusable across every
   plugin as configuration rather than code.
3. **Engine work** — wiring the already-existing VCF/BED readers into `plugin_cache` (they are
   `NotImplemented` today) plus four small hardening fixes.
4. **A manifest scaffolder skill**.

**Definition of done:** `parity --check` passes for three clients, each proving a different axis —
AlphaMissense (the harness is correct), REVEL (a newcomer can port a TSV plugin from a manifest
alone), and one VCF-sourced plugin (the new provider wiring works end-to-end).

---

## 1. Corrections to the prior handoff

The 2026-07-12 handoff proposed a copy-paste runbook. Four of its premises are false against the
code, and following it would waste a session:

| Handoff claim | Reality |
|---|---|
| `vepyr-plugins` is inaccessible; the blocker is repo access (§5) | It is cloned locally and public to us. It has a **working AlphaMissense manifest on `master`** — a real, parity-passed reference. |
| "Python surface (already shipped): `vepyr.build_plugin_cache(...)`" (§2) | **Does not exist.** `grep -r plugin_cache` in the `vepyr` repo returns nothing. The builder is `PluginCacheBuilder` (Rust) driven by `cargo run --example build_plugin --features parquet-cache`. |
| Seed manifests for REVEL/PrimateAI (§6.2–6.3), to be pasted in | Written blind. They declare `has_header = true` and **no `schema` block** — but `provider.rs::csv_schema` builds the Arrow schema *solely* from `[source.csv].schema`, so an absent block yields a zero-field schema. The only working manifest (AlphaMissense) uses `has_header = false` + an explicit `schema`. Discard the seeds; regenerate from that pattern. |
| Bucket A is "framework READY, port = manifest only" | Only for `csv`/`tsv`/`parquet`. `provider.rs:111-116` returns `NotImplemented` for **`vcf`** and **`bed`**, and roughly a third of Bucket A's sources are tabix VCF (EVE, Mastermind, Geno2MP, gnomADMt, SubsetVCF). |

Two facts the handoff missed, both favourable:

- **A reference Ensembl VEP is already installed** (`ensembl-vep/`, cache v115, GRCh38 FASTA). The
  feasibility spec called standing this up "the true gating cost of systematic testing". It is
  largely paid.
- **The format readers already exist.** `datafusion-bio-formats` ships `bio-format-vcf` (TableProvider
  with INFO projection, indexed read, projection pushdown), `bio-format-bed`, and `bio-format-bbi`
  (bigWig/bigBed). So the missing providers are a **wiring** job, not a format-implementation job.
  This also cheapens Bucket C, which the feasibility spec priced as MEDIUM on the assumption that a
  bigwig provider had to be built.

---

## 2. Base branch — `master-sitekwb` (and why it is not `dev-test`)

`plugin_cache` lives only on `master` (14 files). `dev-test` has **zero** of them, so nothing in this
spec compiles there.

The obvious move — "sync `dev-test` with `master`" — is **not** the mechanical merge it looks like.
The branches diverged at `v0.10.0` and both moved:

- `master` took `a6e19ad` (#181, *Lance-only cache, grid-aligned parallel annotation*), which rewrote
  `annotate_provider.rs` (**+7130/−2576**) and **deleted** `variant_lookup_exec.rs`.
- `dev-test` has **45 commits** of its own, including an engine feature that exists nowhere else:
  multi-ALT CSQ per-allele expansion (`PerAltCtx`, `vcf_to_vep_allele_multi`). `master` still treats
  `alt` as a single allele.

So merging them means **porting a feature into a rewritten hot path**: 149 files merge cleanly, and
the 5 conflicting hunks are all CSQ per-allele assembly. Get it wrong and nothing crashes — every
multi-allelic variant is just silently mis-annotated. That deserves its own PR and its own parity
test, and it has **nothing to do with the plugin factory**.

**Decision (2026-07-13): decouple.** All factory work is based on **`master-sitekwb`** (cut from
`master`, so it already has `plugin_cache`). It is treated as our `main`: every PR targets it, and
**`master`/`main` are never committed to**. The multi-ALT port from `dev-test` is tracked separately
and blocks nothing here.

---

## 3. Goals / non-goals

**Goals.** Make porting a Bucket A/B plugin cost a manifest plus a config row. Make the parity gate
executable and reproducible by someone who is not the engine's author. Make a failing gate say
*which* thing is broken.

**Non-goals for this spec** (each is a separate project; recorded in §9 with the plugins that will
demand them):
- Gene/transcript-ID-keyed lookup (Bucket D).
- Interval/overlap lookup — the cache is point-lookup only; `end` is stored but never read.
- `match_column` OR/fallback semantics.
- Widening `ValueType` beyond `Utf8`/`Float32`/`Int32`.
- The bigwig provider (no client in Wave 0; wiring is now known to be cheap — Wave 2).

---

## 4. Architecture: three repos, three roles

Ownership today is smeared across repos: the CSQ comparator lives in `vepyr/e2e-testing/scripts/`,
the builder is a cargo example in `datafusion-bio-functions`, and the manifests live in
`vepyr-plugins`. A PR adding a plugin is therefore not testable by itself. The fix:

- **`datafusion-bio-functions` — the engine.** Mechanisms only: providers, builder, lookup, template.
  It knows nothing about any specific plugin. Tests are per-mechanism unit tests.
- **`vepyr` — the toolkit.** Exposes the two things the handoff *assumed* existed:
  - `vepyr.build_plugin_cache(...)` — a PyO3 binding over `PluginCacheBuilder`.
  - `vepyr.parity` — the CSQ comparator, lifted out of
    `e2e-testing/scripts/run_annotation_fast.py:384` (`compare_vcfs`) into an importable module.
    One comparator, two consumers: the WGS e2e suite and plugin parity. It already handles the hard
    parts (merge-join on the variant key, CSQ entry-count and entry-order semantics, per-field
    mismatch counts with examples).
- **`vepyr-plugins` — the catalogue.** Data and configuration; no engine code. Per plugin: manifest,
  source slice, golden, parity declaration. Its CI depends on `pip install vepyr` plus the
  mini-cache. **It never compiles Rust.**

---

## 5. Component 1 — the mini-cache fixture

The full cache is 34 GB; `variation/chr22.parquet` alone is 541 MB and the FASTA is 2.9 GB. Hosted
CI cannot rebuild the vepyr side, and none of it is committable. So slice once, narrowly.

- **Region: `chr22:22,000,000–23,500,000`** (1.5 Mb). Not arbitrary — it contains the locus on which
  AlphaMissense parity was manually confirmed (`chr22:22,893,742 C>G`), so the regression in §8 is
  meaningful.
- **Slicer:** `scripts/slice_cache.py` in `vepyr`. Filters **every** cache table (`variation`,
  `transcript`, `exon`, `translation_core`, `translation_sift`, `motif`, `regulatory`) to the region
  and cuts the matching FASTA window. Expected output: **tens of MB** (order 50–100 MB).
- **Distribution: a GitHub release asset**, restored in CI via `actions/cache`. Not git-lfs, not a
  private bucket.
- **Input VCF:** the HG002 benchmark sliced to the region — a few hundred variants, tens of KB,
  **committed** to `vepyr-plugins`.

**Stated honestly:** the "1,912/1,912" figure from the manual AlphaMissense parity **cannot** be
reproduced in CI — that was full chr22. What CI reproduces is *parity on the mini-cache region*. The
full-chr22 run remains a nightly/manual check on a machine with the 34 GB cache.

---

## 6. Component 2 — the parity harness

One tool, parameterised per plugin by `plugins/<name>/parity.toml`. Two modes:

**`--refresh-golden` (local; needs Perl + the plugin's data).** Runs the **real Ensembl VEP** with
`--plugin <X>` over the region and writes `plugins/<name>/golden/<name>.vcf`. This is the *only*
place Perl is required, and it runs when the plugin or its data version changes — not per PR.

**`--check` (CI; hermetic).** Takes the committed golden, builds the plugin cache with
`vepyr.build_plugin_cache` against the mini-cache, annotates the committed input VCF, and diffs with
`vepyr.parity`.

### 6.1 The gate, and the blame-attribution rule

A naive diff conflates three failure sources, because plugin fields are *derived from engine
attributes*: the AlphaMissense discriminator is `{ref_aa}{Protein_position}{alt_aa}`, so if vepyr's
core annotation disagrees with VEP about the transcript or the amino-acid change, the plugin field
comes out empty or different — and a naive harness reports "REVEL fails parity" when REVEL's manifest
is perfect and the core is at fault. The predictable consequence is someone "fixing" parity by
loosening the gate.

This is not hypothetical: `e2e-testing/reports/` contains `mismatches_parquet_vs_vep_classified.tsv`
and `issue88_remaining_unresolved.md` — the core *has* known divergences from VEP.

So from the same pair of files (golden VEP-with-plugin ↔ vepyr-with-plugin-cache), and with no extra
runs, compute two independent verdicts:

1. **Core agreement.** Diff the CSQ fields the discriminator depends on: `Feature`, `Consequence`,
   `Amino_acids`, `Protein_position`, `ref`/`alt`. This is *not* the gate — it is a property of the
   core, guarded by the existing e2e suite. Here it only determines the set of lines on which vepyr
   and VEP already agree.
2. **The port gate.** Diff the plugin's CSQ fields **restricted to that agreed subset**. If core and
   VEP agree about the transcript and the amino-acid change, then a plugin-field mismatch is
   unambiguously the manifest's or `plugin_cache`'s fault. Nothing else.

**Pass criteria:** 100% agreement on the plugin's fields over the agreed subset; **zero
over-emission** (vepyr must not populate a field where VEP left it empty — one of the two bugs the
manual AlphaMissense e2e caught that unit tests did not); and **w1-vs-w4 body-identity** (free, and
it catches races in the lookup).

Lines excluded for core drift are reported **loudly and separately** (`excluded: core drift, N
lines`) — never silently folded into zero. A rising exclusion rate is a signal about the core.

A `--strict` flag makes core drift fail the build too; that is the mode for PRs against the engine.

### 6.2 One harness, two subjects under test

The same tool tests different things depending on what changed:

- **PR in `vepyr-plugins` (a new manifest)** → the subject is the **manifest**. Engine and vepyr are
  pinned; run that one plugin.
- **PR in `datafusion-bio-functions` (an engine change)** → the subject is the **engine**. Manifests
  are pinned; run the **whole plugin corpus** as a regression net.

The corpus, built once for the plugins' sake, becomes a free safety net for the engine.

### 6.3 Redistribution

Source slices are often not redistributable (REVEL is non-commercial; AlphaMissense is CC-BY-NC).
`parity.toml` therefore carries `redistributable = true|false`. When `false`, the plugin's test is
**skipped in public CI with an explicit, visible message** — never silently — and runs only in the
nightly job on a machine that holds the data. Without this field the factory would either violate a
licence or quietly go green on empty tests.

---

## 7. Component 3 — engine changes

### 7.1 Provider wiring (not implementation)

`ProviderKind` (`plugin_cache/source_manifest.rs:21-29`) already *declares* `vcf`, `csv`, `tsv`,
`parquet`, `bed` — so a manifest with `provider = "vcf"` passes validation today and only explodes at
build time, in `provider.rs:111-116`, with `NotImplemented`.

Wire `bio-format-vcf` under `ProviderKind::Vcf` and `bio-format-bed` under `Bed`. The manifest gains
a `[source.vcf]` block declaring which INFO fields to project — exactly analogous to the existing
`[source.csv].schema`. INFO arrays (per-allele values) are handled in `ingest_sql` via the existing
`unnest` escape hatch.

### 7.2 Manifest hardening

Four small fixes, each of which saves hours of debugging:

- **`#[serde(deny_unknown_fields)]` on `SourceManifest`.** Today a typo — or the `[tier]` block that
  the handoff still documents and the code does not have — is **silently ignored**. The manifest
  "works" and quietly does something else.
- **`--overwrite` stops being a no-op** (`builder.rs:106` literally reads `let _ = self.overwrite;`).
- **`--source-path` applies to all sources**, not just `sources.first_mut()` — required for
  multi-part sources.
- **`--chrom` accepts repetition** (today only the first occurrence is read).

---

## 8. Component 4 — the scaffolder, and the repo layout

**Scaffolder:** a `vep-add-plugin` **skill**, not a CLI — because half the work is *reading the
`.pm` header and deciding which bucket the plugin belongs to*, which is model work, not parser work.
Input: plugin name + data URL. Output: a `<name>.source.toml` stub on the AlphaMissense pattern, a
`parity.toml`, an empty `golden/`, and the §4 playbook checklist. It mechanically enforces exactly
what the handoff's blind seeds got wrong: an explicit `schema` instead of `has_header`, the TOML
scalar-ordering trap (top-level keys must precede any `[[table]]` header), and a `path` that the
build driver overrides.

**Layout in `vepyr-plugins`:**

```
plugins/<name>/
  <name>.source.toml     # manifest (exists today for alphamissense)
  parity.toml            # csq_fields, region, Perl plugin version, redistributable
  fixtures/<name>.slice  # source slice for the region (if redistributable)
  golden/<name>.vcf      # real VEP --plugin <X> output (written by --refresh-golden)
```

---

## 9. Known engine limits (deliberately out of scope, with their future clients)

Recording these so they are not rediscovered painfully:

- **`ValueType` has only `Utf8`/`Float32`/`Int32`** — no Int64/Float64. Will bite any plugin needing
  double precision.
- **`match_column` supports multiple columns but only conjunctive exact match** — no OR/fallback, so
  VEP's "match on transcript_id, else on aa-change" cannot be expressed. Also: the multi-column path
  has **no test**, and transcript matching is bare string equality with **no version normalisation**
  (`ENST00000123.4` ≠ `ENST00000123`).
- **Point lookup only.** `end` is in the shard schema but never read; `PluginLookup` pages on `start`
  alone. Bucket C interval plugins need real overlap support.
- **The probe keys on `input_start` (VCF POS)**, which may be wrong for indels. AlphaMissense is
  SNV-only, so this has never been exercised.

---

## 10. Acceptance criteria

Wave 0 is done when `parity --check` passes for three clients, each proving a different axis:

1. **AlphaMissense** — the harness reproduces parity on the mini-cache region. This proves *the
   harness is correct* (we have an independently confirmed result to reproduce). Its manifest is
   **not touched** — it is owned elsewhere; we use it strictly as a fixture.
2. **REVEL** — the first new client, TSV source. This proves the thesis "a port is just a manifest"
   holds for someone who did not write the engine.
3. **Mastermind or gnomADMt** — VCF source. This proves the newly wired provider works end-to-end
   rather than merely compiling.

---

## 11. Sequencing

0. Cut `master-sitekwb` from `master` (§2). Every PR below targets it.
1. Engine: provider wiring + hardening (§7). PR → `master-sitekwb`.
2. `vepyr`: `build_plugin_cache` binding + `vepyr.parity` module + `slice_cache.py`; publish the
   mini-cache release asset (§5).
3. `vepyr-plugins`: harness `--check`/`--refresh-golden`, `parity.toml`, CI (§6).
4. The three clients (§10). AlphaMissense first — it is the harness's own test.
5. Only then Wave 1 (the ~20 Bucket A/B scoring plugins) as pure manifest work.

---

## 12. Risks

- **The mini-cache slicer must be consistent across tables.** A transcript sliced away while its
  variants remain would produce phantom core drift and poison the attribution rule. The slicer needs
  its own test: annotate the region with the mini-cache and with the full cache, and require
  identical output.
- **Licence-gated sources cap how hermetic CI can be.** Mitigated by `redistributable = false` +
  loud skips, but it means the public CI corpus will be a subset of the real one.
- **Parity may not reach 100% on allele edge cases** (multi-allelic, MNV, spanning deletion, indel
  left-alignment) even in the easy class. The attribution rule keeps these honest instead of letting
  them be papered over.
- **The scaffolder can encode today's manifest shape and then rot.** It must be regenerated from the
  `SourceManifest` struct, not from a copied example — otherwise it will keep emitting `[tier]`-style
  ghosts.
