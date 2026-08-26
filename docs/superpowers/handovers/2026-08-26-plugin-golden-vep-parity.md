# Golden-VEP plugin parity — findings handover (2026-08-26)

The acceptance gate from `2026-08-25-plugin-cache-testing-handover.md` §6 has been
run. This records what was analysed, what was run, what mismatched, and the
traps — several of which fail silently and will cost the next person a day.

Engine: `86fd72d` (this repo, `plugins-fixes`).
Manifests: `vepyr-plugins` `c9fe785`, `f9cfa30`, `2545530` (local, unpushed).
Docs + reference script: `vepyr` `plugins-fixes` (local, unpushed).
Narrative report: <https://claude.ai/code/artifact/28178c87-b954-40b0-adee-dfd2933681de>

---

## 1. Result

All 38 plugin CSQ fields report **100.00%** on both contigs.

| contig | variants | CSQ entries/field | mismatch rows | format-only | unpaired | CSQ count/order |
|---|---:|---:|---:|---:|---:|---:|
| chr21 | 55,812 | 895,482 | **2** | 0 | 0 | 0 |
| chr1 | 323,430 | 5,688,867 | **78** | 0 | 0 | 0 |

Byte parity (chr21, `md5_concordance.py`): **55,811 of 55,812 records byte-identical**.
Header 272 vs 266 lines — see §7.

Every remaining row is the upstream VEP HGVS-shift defect (§5.4), not ours.

**The caches were never wrong.** In every failing case the correct value was
physically present in the shard and the lookup asked the wrong question. No
shard was rebuilt for correctness — only to change column *types* (§5.5).

---

## 2. Files analysed

### Ensembl (ground truth — read, not inferred)

Extract from the image with
`docker run --rm ensemblorg/ensembl-vep:release_116.0 bash -lc 'sed -n "/^sub X/,/^}/p" <file>'`.

| file | why it mattered |
|---|---|
| `Bio/EnsEMBL/Variation/Utils/Sequence.pm` → `get_matched_variant_alleles` | the actual matcher: minimises **both** sides, compares `(ref, alt, pos)` |
| … → `trim_sequences` | prefix-then-suffix vs suffix-then-prefix; `empty_to_dash` |
| … → `_get_trim_directions` | returns `[0,1]` when either allele is longer than 1 |
| `Plugins/CADD.pm` `run`/`parse_data` | delegates to `get_matched_variant_alleles`; fetches `[start-2, end]`; anchor-shifts data rows |
| `Plugins/SpliceAI.pm` `run` | **imports** the matcher but never calls it — exact `eq` comparison |
| `Plugins/dbNSFP.pm` `run` | exact `pos`+`alt`; `return {} unless $vf->{start} == $vf->{end}` (SNV-only) |
| `Plugins/AlphaMissense.pm` `run` | calls `get_matched_variant_alleles` |

Plugin `.pm` files are also mirrored in the Drive sources folder, which is where
they were taken from for the container's `--dir_plugins`.

### Engine

`src/allele.rs` (already contained a full `get_matched_variant_alleles` port —
`trim_sequences_ensembl`, `trim_directions`), `plugin_cache/{registry,lookup,
cache_manifest,source_manifest,build,normalize}.rs`, `annotate_provider.rs`
(probe call sites ~6425/6473, buffer take-set ~5470), `vcf_sink.rs`
(`csq_header_description`, raw header lines).

### Comparator (`vepyr` PR #48)

`comparison/{compare,profiles,report,cli}.py`, `md5_concordance.py`,
`verify_parity_gate.py`.

---

## 3. Commands run

### Reference (the part that did not exist)

`e2e-testing/scripts/build_vep_plugin_reference.sh <chrom>` in `vepyr` captures
this end to end. Core call:

```bash
docker run --rm --user "$(id -u):$(id -g)" -v "$DATA":/data \
  ensemblorg/ensembl-vep:release_116.0 \
  vep --cache --cache_version 116 --dir_cache /data --offline --merged \
      --everything --no_stats --force_overwrite --vcf \
      --fasta /data/input/Homo_sapiens.GRCh38.dna.primary_assembly.fa \
      --input_file /data/vep116_chr21/input/HG002_norm_chr21.vcf.gz \
      --output_file /data/vep116_chr21/output/HG002_chr21_5plugins_vep116.vcf \
      --dir_plugins /data/vep116_chr21/plugins \
      --custom …/clinvar_chr21.vcf.gz,ClinVar,vcf,exact,0,CLNSIG,CLNREVSTAT,CLNDN,CLNVC,CLNVI \
      --plugin SpliceAI,snv=…/spliceai_chr21.vcf.gz,indel=…/spliceai_chr21.vcf.gz \
      --plugin CADD,snv=…/cadd_all_chr21.tsv,indels=…/cadd_all_chr21.tsv \
      --plugin AlphaMissense,file=…/alphamissense_chr21.tsv.gz \
      --plugin dbNSFP,…/dbNSFP5.3.1a_grch38_chr21.gz,<19 cols>
```

Runtime: chr21 ~14 min, chr1 ~75 min (native arm64, no emulation).

Source slices via `rclone serve http` + `tabix` (handover 2026-08-25 §3).
Contig naming: `chr21` for AlphaMissense, bare `21` for the rest.

### Shard rebuilds (only for the `Utf8` retype)

```python
vepyr.build_plugin_cache("cadd", "f9cfa30", source_path=".../cadd_all_chr21.tsv",
    cache_dir="$DATA/cache/116_GRCh38_merged", plugin_cache_root="$DATA/plugin_cache",
    chroms=["21"], plugins_repo="$HOME/workspace/vepyr-plugins", overwrite=True)
```

`version` is a **git ref** materialised via `git worktree`, so manifest edits
must be **committed** first — working-tree edits are invisible.

Times: spliceai chr21 17 s / chr1 209 s; cadd chr21 45 s / chr1 296 s.
Row counts reproduced the originals exactly (cadd chr21 121,809,062 / warm
399,098; chr1 699,768,898 / warm 1,929,017).

### Comparison

```bash
uv run python run_comparison.py --release 116 --profile merged_plugins \
  --chroms 21 --plugin-cache $DATA/plugin_cache \
  --vep $DATA/vep116_chr21/output/HG002_chr21_5plugins_vep116.vcf.gz \
  --workers 4 --force

python md5_concordance.py --pair <vep>.vcf <vepyr>.vcf --mode strict --explain
```

Isolating an upstream behaviour — run VEP on a handful of variants with
`--fields Location,Allele,Feature,Consequence,CADD_RAW` and bisect flags. This
is how §5.4 was pinned; do this before theorising.

---

## 4. Mismatches found

### chr21 (per CADD field, before fixes: 232)

| # | variant | shape | entries | cause |
|---|---|---|---:|---|
| 1 | chr21:13973877 `TTGTGTGTGTGTG>GTGTGTGTGTGTG` | equal-length MNV | 85 | A |
| 2 | chr21:26032805 `AAT>TAT` | equal-length MNV | 68 | A |
| 3 | chr21:26062230 `AAC>ACAC` | insertion | 66 | A (indel) |
| 4 | chr21:28149937 `CGT>TGT` | MNV | 5 | A |
| 5 | chr21:39327921 `AAATG>GAATG` | MNV | 3 | A |
| 6 | chr21:27136796 `CA>AA` | MNV | 2 | A |
| 7 | chr21:18829623 `AATAT>TATAT` | MNV | 1 | A |
| 8 | chr21:26989261 `CTC…>TTC…` (33 bp) | MNV | 1 | A |
| 9 | chr21:33549183 `A>ATTGT` | insertion | 1 | **upstream** |

Plus 1,533,782 format-only entries over six `Float32` fields, and 32 missing
header lines. **Now: 2 rows (variant 9 only).**

### chr1 (before the ClinVar fix: 138 rows)

| field group | rows | variants | cause |
|---|---:|---:|---|
| CADD_RAW/PHRED | 39 each | 5 insertions (13779899, 77558660, 152709213, 111290866, 217771848) | **upstream** |
| ClinVar ×6 | 10 each | chr1:65364614 `GT>TT` | misclassified `allele_match` |

**Now: 78 rows (CADD only).**

---

## 5. Causes, and what fixed them

### 5.1 Probe key not reduced — and trim *order* decides the row

`vcf_to_vep_input_allele` strips only a single anchor base. `bcftools norm
-m -both` splits multi-allelics **without re-trimming**, so non-minimal records
are routine. Added `allele::plugin_probe_allele` (left-first `trim_sequences`),
consulted **only on a miss** so no existing hit can change.

Order is load-bearing, not cosmetic:

| VCF | left-first (Ensembl) | right-first |
|---|---|---|
| `CGTGTGT/CGTGT` | `GT/-` @ 13836153 — no row, empty (correct) | `GT/-` @ 13836149 — another variant's score |
| `AAC/ACAC` | `-/C` @ 26062231 — the row VEP reports | same |

### 5.2 Allele semantics are per-plugin

`AlleleMatch::{Exact,Minimised}` in the manifest. CADD, AlphaMissense and
`--custom` (ClinVar) minimise; SpliceAI and dbNSFP compare verbatim.

> The `exact` in `--custom <file>,ClinVar,vcf,exact,0,…` is the **overlap mode,
> not the allele rule**. VEP matches `chr1:65364614 GT>TT` to ClinVar's `G>T`
> row 1258041. Misreading this cost 10 variants × 6 fields on chr1.

### 5.3 CSQ field order

Ensembl emits plugin blocks in **command-line order**, each plugin's fields
**sorted by name**, `--custom` last in declared order. We emitted plugins
alphabetically in declaration order → 99.3% of records differed, invisible to
the name-keyed comparator. Pinned via `csq_rank` + `field_order`.

### 5.4 Upstream: HGVS shifting breaks CADD's fetch window

`CADD.pm` fetches `[VF.start-2, VF.end]` from the tabix file. HGVS 3′ shifting
mutates those VF coordinates, so the window moves off the row and the score
vanishes — for the shifted (frameshift) transcript only.

```
--hgvs                 → chr1 80/139 hits   (39 lost = exactly our mismatches)
--hgvs --shift_hgvs 0  → chr1 119/139 hits
plain (no --hgvs)      → chr1 119/139 hits
```

A cosmetic naming feature corrupting a data lookup. vepyr matches
`--shift_hgvs 0` and matches CADD being a per-variant score. **Not reproduced.**
Worth reporting upstream.

### 5.5 Float re-formatting — retyping alone is not enough

Six `Float32` columns (CADD RAW/PHRED, SpliceAI DS_×4) were re-rendered with
Rust's shortest round-trip. Setting them to `Utf8` changed the stored type but
**not the value**, because `ingest_sql` still did `CAST(raw_score AS FLOAT)` —
the padding was destroyed before storage. Selecting the raw columns verbatim
fixed it. Now 0 on both contigs.

Not "avoid Float32": `am_pathogenicity` stays `Float32` and is clean, because
AlphaMissense doesn't zero-pad. **Don't re-format a source that pads.**

### 5.6 Missing plugin header lines

`description` per value column, emitted as `##<FIELD>=…` from `vcf_sink.rs`.
All 32 byte-identical.

---

## 6. Observations — silent failure modes

Ranked by how long each cost.

1. **`ProjectionMask::leaves` ignores leaf order.** Permuting the projection to
   match reordered field names moves the names and leaves the values put; every
   field pairs with a sibling's value. CADD fell to 0.73%, ledger 4.6M rows.
   **Permute the returned values, not the projection.**
2. **A rebuild silently reverts anything only patched into `manifest.json`.**
   `allele_match` was hand-patched, lost on the next build, and CADD jumped back
   to 232. Settings belong in the source TOML.
3. **dbNSFP parses its version from the *filename*.** A slice named
   `dbnsfp_chr21.tsv.gz` is rejected, all 19 fields disappear from the CSQ
   header, and VEP still **exits 0**. The name must contain e.g. `5.3.1a`.
   The reference script now asserts 38 plugin fields for this reason.
4. **`--database 0` is not valid input.** It appears in VEP's own
   `##VEP-command-line` header but parses as a stray positional. `--offline`
   covers it.
5. **CADD's slice must be plain text and comment-free.** Its manifest declares
   no `compression` (a `.gz` is read as binary), and `tabix -h` leaves a
   `## CADD` line that breaks the CSV reader. Concatenate SNV+indel with
   comments stripped. chr1 combined = 20 GB — watch disk.
6. **The summary markdown hides format-only and order-only counts.**
   `report.py:512` / `cli.py:571` compute "fields at 100%" from real mismatches
   alone. Read plugin results from the **per-contig JSON**
   (`field_format_mismatch_counts`). This is very likely how the earlier
   "37/38 at 100%" chr22 figure arose — that run also predates identity pairing.
7. **`build_plugin_cache(version=…)` is a git ref**, materialised via
   `git worktree`. Uncommitted manifest edits are invisible.
8. **A chrom-scoped rebuild rewrites the whole manifest**, dropping other
   chroms' entries. Re-add them, or the sibling contig silently disappears.

### Method note

Five of my explanations were overturned by the next measurement: the
`vcf_to_vep_allele` one-liner (a no-op — it doesn't suffix-trim MNVs), the
"equal-length only" restriction (the real rule was trim order), "cause B is a
start-convention bug" (it was the same defect as A), "VEP is just inconsistent"
(it's HGVS shifting), and ClinVar's `exact`. Each died to reading Ensembl's
source or running VEP directly. **Read the `.pm`, or run the flag bisection.
Do not reason from mismatch shapes.**

Also: chr1 caught a defect chr21 could not. Two contigs found four defects.

---

## 7. Open items

1. **Header parity** — 272 vs 266. The gap is six `##INFO=<ID=ClinVar…>`
   definitions VEP writes for `--custom`, each embedding the **absolute source
   path**, so no implementation can match them portably. They also declare INFO
   keys that never appear as INFO fields (values ride in CSQ). Suggest
   `md5_concordance.py` treat them as provenance rather than emitting untrue
   INFO definitions to win a digest.
2. **chr2–20, 22, X** — untouched. Given chr1's yield, expect more.
3. **Push the companion commits** — `vepyr-plugins` (3) and `vepyr` (docs +
   `build_vep_plugin_reference.sh`) are local only.
4. **`verify_parity_gate.py` still refuses plugin profiles** by design; the
   plugin field contract is unpinned there.
5. **Report the HGVS-shift/CADD-window bug upstream.**
6. **Other contigs' shards still carry the old `Float32` columns** — any contig
   beyond chr1/chr21 must be rebuilt from the current manifests, or its
   manifest will disagree with its columns.

## 8. Artifacts on disk

| what | where |
|---|---|
| chr21 / chr1 references (+ `.gz`, `.tbi`) | `$DATA/vep116_chr{21,1}/output/` |
| normalized inputs | `$DATA/vep116_chr21/input/HG002_norm.vcf.gz` (whole-genome), per-contig slices alongside |
| plugin caches (all 5, chr1+chr21) | `$DATA/plugin_cache/plugin/<name>/` |
| reports + mismatch ledgers | `vepyr/e2e-testing/reports/fast_chr{1,21}_merged_plugins_116_*` |
| prior-session plain slices (reusable) | `$DATA/plugin_input/_slices/` |
