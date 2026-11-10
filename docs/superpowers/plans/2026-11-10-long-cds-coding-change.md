# Long-CDS Coding Change Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Coding-change classification stops doing work proportional to the whole coding sequence for every variant-transcript pair, so dense variants on very long transcripts (TTN) cost about what they cost on short ones.

**Architecture:** `classify_coding_change_with_spliced_span_cached` builds `mutated`, a full copy of the CDS, for every pair; `substitute_peptide_window` then copies the whole reference peptide; every indel translates the whole mutated CDS. This plan first measures which of these dominates on the TTN block, then replaces them one at a time: a lazy `MutatedCds` view, an `AltPeptide` view that stores only the changed codons, and indel translation that reuses the unchanged reference prefix. Each fast path sits behind `CodingFastPaths`, and an equivalence test compares every path with the current code on a variant matrix.

**Tech Stack:** Rust 2024; changes in `datafusion/bio-function-vep/src/transcript_consequence.rs` and `src/transcript_consequence/reference_memo.rs`; `/usr/bin/sample` and `atos` for the profile.

**Spec:** `docs/superpowers/specs/2026-11-10-parallel-vcf-output-bottlenecks-design.md` (sections 1.4, 3.4, 4)

## Global Constraints

- Branch `perf/long-cds-coding-change` off `origin/master`, after plan 1 (it removes chr2 as the scheduling straggler, so this plan is judged on CPU).
- Output bytes identical for every input. `CodingClassification` is identical field by field (compared through `Debug`).
- Fast paths apply only where the substitution fast path applies today: a `ReferenceEntry::Ready` reference with `local_substitution_safe`, not patched for an RNA edit, an ACGT-only allele, no RefSeq edit gap and no spliced intron gap. Everything else keeps the current code.
- `VEP_CODING_FAST_PATHS=0` disables all three in the same binary; the value is read once per process.
- Gate targets (spec section 4): TTN block (41,279 records, `2:178525989-178830802`), 1 worker, CPU time at least 15% lower (median of three runs per build); full gate with `--target-pct HG002=-2,ClinVar=1`.
- Expected size: coding-change classification with translation is 33% of busy samples on ClinVar chr2 and chr2 is 12% of merged ClinVar CPU, so the whole-dataset gain is a few percent at most. Stop after Task 1 if its table shows no item at 5%.
- Line numbers refer to `68d0f7b`.

## Review Focus

- A variant in the first codon and one in the last (stop) codon, on both strands: start_lost, start_retained, stop_lost and stop_retained fields identical (matrix rows in Task 2).
- An insertion or deletion that is not a multiple of three, including one whose new stop lies in the 3' UTR: `protein_hgvs` and `alt_translation_extension` identical (matrix rows in Task 2).
- A transcript whose reference is not `local_substitution_safe` (CDS length not a multiple of three, an `N` in the CDS): fast paths not taken, result identical (matrix rows in Task 2; a counter test in Task 5).
- The vertebrate mitochondrial codon table (`CodonTable::VertebrateMitochondrial`): identical results (matrix rows in Task 2; ClinVar MT in the gate).
- A reference declined by the 128 MiB reference budget (`ReferenceEntry::Declined`): the function still returns the same classification through the owned fallback (test in Task 3).

---

### Task 1: TTN block benchmark and profile

**Files:**
- Modify: `datafusion/bio-function-vep/PERFORMANCE.md` (new section "Long-CDS coding change", baseline table)

**Interfaces:**
- Consumes: `build_venv.sh` (plan 0).
- Produces: `$G/p4/ttn.vcf.gz`; a table with one row per item below, giving its share of busy samples on the annotation thread. Tasks 3, 4 and 5 are done only for items at 5% or more.

- [ ] **Step 1: Build the input and a profiling venv**

```bash
G="$HOME/research/benchmarks/local-gate"; L=datafusion/bio-function-vep/scripts/local_gate; D="$HOME/workspace/data_vepyr"
mkdir -p "$G/p4"; git worktree add --detach "$G/p4/functions-baseline" origin/master
tabix -h "$D/input/clinvar/clinvar_20261011.norm.vcf.gz" 2:178525989-178830802 | bgzip > "$G/p4/ttn.vcf.gz"
tabix -p vcf "$G/p4/ttn.vcf.gz"; zgrep -vc '^#' "$G/p4/ttn.vcf.gz"
$L/build_venv.sh "$G/p4/functions-baseline" "$G/p4/venv-baseline"
CARGO_PROFILE_RELEASE_DEBUG=line-tables-only $L/build_venv.sh "$G/p4/functions-baseline" "$G/p4/venv-profile"
```
Expected: 41279 records.

- [ ] **Step 2: Baseline timing**

Three runs of
```bash
/usr/bin/time -l "$G/p4/venv-baseline/bin/python" -P -c 'import sys,vepyr
vepyr.annotate(sys.argv[1], sys.argv[2], reference_fasta=sys.argv[3], output_vcf=sys.argv[4], workers=1, compression="plain", preserve_record_layout=True, show_progress=False)' \
  "$G/p4/ttn.vcf.gz" "$D/cache/116_GRCh38_merged" "$D/input/Homo_sapiens.GRCh38.dna.primary_assembly.fa" "$G/p4/ttn.out.vcf"
```
Record user + sys seconds and wall for each; keep the median and `grep -v '^#' "$G/p4/ttn.out.vcf" | md5`.

- [ ] **Step 3: Profile**

Start the same command with `venv-profile` in the background and run `/usr/bin/sample <pid> 600 1 -mayDie -file "$G/p4/ttn.sample.txt"`. For the busiest thread, sum samples whose stack contains each symbol, and resolve the self-time offsets of `classify_coding_change_with_spliced_span_cached` to source lines with `atos -o <path to _core*.so> -l <load address from the report> <address>`.

| Item | Symbol or lines | Task |
|---|---|---|
| CDS copy | `_platform_memmove` directly under `classify_coding_change_with_spliced_span_cached`, lines 10032-10043 | 3 |
| Peptide copy | `reference_memo::substitute_peptide_window` | 4 |
| Indel translation | `translate_protein_for_coding_classification`, `translate_protein_from_cds_in` | 5 |
| HGVS alternate translation | `alternate_translation_for_vep_hgvs`, `translated_alt_protein_for_hgvs`, `translate_protein_string` | 5 |
| Classifier self time | the three source lines with the most self samples | note only |

- [ ] **Step 4: Record and decide**

Write the timing and the table into `PERFORMANCE.md`. List which of Tasks 3-5 qualify (share at least 5%). If none does, commit the section with that conclusion and stop.

```bash
git commit -am "docs(vep): TTN block baseline profile"
```

### Task 2: Switch and equivalence matrix

**Files:**
- Modify: `datafusion/bio-function-vep/src/transcript_consequence.rs` (`classify_coding_change_with_spliced_span_cached` at 9824)
- Test: `datafusion/bio-function-vep/src/transcript_consequence/tests/coding_fast_paths.rs` (new; add the `mod` line next to `allocation_free`)

**Interfaces:**
- Produces:
  - `#[derive(Clone, Copy, Debug, PartialEq, Eq)] pub(crate) struct CodingFastPaths { pub lazy_cds: bool, pub peptide_view: bool, pub indel_prefix: bool }` with `const ALL` and `const NONE`, and `fn from_env_value(value: Option<&str>) -> Self` (`NONE` for a value that `pipeline_trace::enabled_from_env_value` rejects, `ALL` otherwise and when unset).
  - `fn classify_coding_change_with_spliced_span_cached_using(…same parameters…, fast: CodingFastPaths) -> Option<CodingClassification>`; the existing function calls it with a process-wide `OnceLock<CodingFastPaths>` read from `VEP_CODING_FAST_PATHS`.
  - Test helpers in the new file: `fn long_transcript(strand: i8, codons: usize, table: CodonTable) -> (TranscriptFeature, Vec<ExonFeature>, TranslationFeature)` (three exons, ACGT CDS from a fixed-seed generator, start `ATG`, one terminal stop, 300 nt 3' UTR) and `fn matrix(tx: &TranscriptFeature) -> Vec<VariantInput>`.

- [ ] **Step 1: Write the failing test**

```rust
#[test]
fn fast_paths_match_the_reference_implementation() {
    for table in [CodonTable::Standard, CodonTable::VertebrateMitochondrial] {
        for strand in [1, -1] {
            let (tx, exons, tl) = long_transcript(strand, 3_000, table);
            let exon_refs: Vec<&ExonFeature> = exons.iter().collect();
            let memo = ready_memo(&tx, &tl);      // a ReferenceMemo with an ample budget
            for variant in matrix(&tx) {
                let slow = classify_coding_change_with_spliced_span_cached_using(
                    &tx, &exon_refs, Some(&tl), &variant, VepSemantics::V116, None, Some(&memo), CodingFastPaths::NONE);
                for fast in [CodingFastPaths::ALL,
                             CodingFastPaths { lazy_cds: true, ..CodingFastPaths::NONE },
                             CodingFastPaths { peptide_view: true, ..CodingFastPaths::NONE },
                             CodingFastPaths { indel_prefix: true, ..CodingFastPaths::NONE }] {
                    let got = classify_coding_change_with_spliced_span_cached_using(
                        &tx, &exon_refs, Some(&tl), &variant, VepSemantics::V116, None, Some(&memo), fast);
                    assert_eq!(format!("{got:?}"), format!("{slow:?}"), "{variant:?} {fast:?} strand {strand}");
                }
            }
        }
    }
}

#[test]
fn coding_fast_paths_env() {
    assert_eq!(CodingFastPaths::from_env_value(None), CodingFastPaths::ALL);
    assert_eq!(CodingFastPaths::from_env_value(Some("0")), CodingFastPaths::NONE);
}
```
`matrix` returns, at CDS codons 1, 2, 1500, 2999 and 3000 (the stop): an SNV at each codon position; a 2-base MNV; a 3-base MNV across two codons; insertions of 1, 2, 3 and 4 bases; deletions of 1, 2, 3 and 6 bases; a deletion-insertion (`delins`) of 2 for 5 bases; an SNV with ALT `N`; and one deletion spanning an exon boundary. Repeat the whole test body for `VepSemantics::V115`. Add a second transcript per strand whose CDS has one `N` (not `local_substitution_safe`) and run the same matrix on it.

- [ ] **Step 2: Run it to verify it fails**

Run: `cargo test -p datafusion-bio-function-vep --lib coding_fast_paths`
Expected: compile error, `CodingFastPaths` not found.

- [ ] **Step 3: Implement the switch only**

Add the struct and the `_using` function; the flags are accepted and not yet read, so every arm runs the current code.

- [ ] **Step 4: Run the test**

Run: `cargo test -p datafusion-bio-function-vep --lib coding_fast_paths`
Expected: PASS (all arms identical by construction). The matrix must classify at least one variant as each of: synonymous, missense, stop_gained, frameshift, in-frame insertion, in-frame deletion, start_lost, stop_lost; assert that in the test so the oracle cannot go vacuous.

- [ ] **Step 5: Commit**

```bash
git commit -am "test(vep): equivalence matrix and switch for coding-change fast paths"
```

### Task 3: Lazy mutated CDS (if Task 1 lists it)

**Files:**
- Modify: `datafusion/bio-function-vep/src/transcript_consequence.rs` (lines 10032-10043 and every later use of `mutated` in the function)
- Modify: `datafusion/bio-function-vep/src/transcript_consequence/reference_memo.rs` (`substitute_peptide_window` at 219)

**Interfaces:**
- Consumes: `CodingFastPaths::lazy_cds` (Task 2).
- Produces:
  - `pub(super) struct MutatedCds<'a> { reference: &'a [u8], start: usize, end: usize, alt: &'a [u8], full: std::cell::OnceCell<Vec<u8>> }` (`end` inclusive)
  - `impl<'a> MutatedCds<'a> { fn new(reference: &'a [u8], start: usize, end: usize, alt: &'a [u8]) -> Self; fn len(&self) -> usize; fn codon(&self, codon_index: usize) -> Option<[u8; 3]>; fn window(&self, range: std::ops::Range<usize>) -> Vec<u8>; fn starts_with_atg(&self) -> bool; fn tail(&self, offset: usize) -> Vec<u8>; fn full(&self) -> &[u8]; fn materialised(&self) -> bool }`
  - `substitute_peptide_window(table: CodonTable, reference: &[char], mutated: &MutatedCds<'_>, start: usize, end: usize, memo: &ReferenceMemo) -> Option<Vec<char>>`

- [ ] **Step 1: Write the failing tests**

```rust
#[test]
fn mutated_cds_view_matches_the_materialised_sequence() {
    let reference = b"ATGAAACCCGGGTTTTAG";
    for (start, end, alt) in [(3usize, 3usize, &b"C"[..]), (3, 5, b""), (6, 5, b"GG"), (0, 2, b"ATA"), (15, 17, b"TGGA")] {
        let view = MutatedCds::new(reference, start, end, alt);
        let mut owned = reference[..start].to_vec();
        owned.extend_from_slice(alt);
        owned.extend_from_slice(&reference[end + 1..]);
        assert_eq!(view.len(), owned.len());
        assert_eq!(view.full(), owned.as_slice());
        assert_eq!(view.starts_with_atg(), owned.starts_with(b"ATG"));
        for ci in 0..owned.len() / 3 {
            assert_eq!(view.codon(ci).unwrap().as_slice(), &owned[ci * 3..ci * 3 + 3]);
        }
        assert_eq!(view.tail(4), owned[4..].to_vec());
        assert_eq!(view.window(2..7), owned[2..7].to_vec());
    }
}

#[test]
fn a_local_substitution_does_not_materialise_the_cds() {
    // SNV at codon 1500 of long_transcript(1, 3_000, Standard), CodingFastPaths { lazy_cds: true, .. NONE }.
    // A cfg(test) counter incremented in MutatedCds::full() stays at 0 for this call.
}

#[test]
fn a_declined_reference_still_classifies() {
    // Same SNV with a ReferenceMemo whose budget is 0 (ReferenceEntry::Declined):
    // the result equals the CodingFastPaths::NONE result.
}
```
(`(6, 5, b"GG")` is an insertion: `end + 1 == start`.)

- [ ] **Step 2: Run to verify they fail**

Run: `cargo test -p datafusion-bio-function-vep --lib mutated_cds a_local_substitution a_declined_reference`
Expected: compile error, `MutatedCds` not found.

- [ ] **Step 3: Implement**

With `fast.lazy_cds`, construct `MutatedCds` in place of the `Vec`. `substitute_peptide_window` reads codons through `codon()`; `mutated_cds_starts_with_atg` becomes `starts_with_atg()`; the local alternate codon window uses `window()`. Every other consumer takes `mutated.full()`. Without the flag, call `full()` once right after construction, which reproduces today's cost.

- [ ] **Step 4: Run the tests**

Run: `cargo test -p datafusion-bio-function-vep --lib mutated_cds coding_fast_paths a_local_substitution a_declined_reference`
Expected: PASS.

- [ ] **Step 5: Commit**

```bash
git commit -am "perf(vep): build the mutated CDS only when a consumer needs all of it"
```

### Task 4: Peptide view for local substitutions (if Task 1 lists it)

**Files:**
- Modify: `datafusion/bio-function-vep/src/transcript_consequence/reference_memo.rs` (`substitute_peptide_window`)
- Modify: `datafusion/bio-function-vep/src/transcript_consequence.rs` (every use of `new_aas` after line 10063; `build_protein_hgvs_data`)

**Interfaces:**
- Consumes: `MutatedCds` (Task 3), `CodingFastPaths::peptide_view` (Task 2).
- Produces:
  - `pub(super) enum AltPeptide<'a> { Window { reference: &'a [char], first: usize, replacement: smallvec::SmallVec<[char; 4]> }, Owned(Vec<char>) }`
  - `impl AltPeptide<'_> { fn len(&self) -> usize; fn get(&self, index: usize) -> Option<char>; fn first(&self) -> Option<char>; fn iter(&self) -> impl Iterator<Item = char> + '_; fn range_string(&self, range: std::ops::Range<usize>) -> String; fn position_from(&self, from: usize, pred: impl Fn(char) -> bool) -> Option<usize>; fn to_vec(&self) -> Vec<char> }`
  - `substitute_peptide_window(..) -> Option<AltPeptide<'a>>` returning `Window`.
  - `build_protein_hgvs_data` takes `new_aas: &AltPeptide<'_>`.

- [ ] **Step 1: Write the failing tests**

```rust
#[test]
fn alt_peptide_window_behaves_like_the_owned_peptide() {
    let reference: Vec<char> = "MKPGF*".chars().collect();
    let window = AltPeptide::Window { reference: &reference, first: 2, replacement: smallvec::smallvec!['L', 'W'] };
    let owned = AltPeptide::Owned("MKLWF*".chars().collect());
    assert_eq!(window.len(), owned.len());
    assert_eq!(window.to_vec(), owned.to_vec());
    assert_eq!(window.first(), Some('M'));
    for i in 0..7 { assert_eq!(window.get(i), owned.get(i)); }
    assert_eq!(window.range_string(1..5), "KLWF");
    assert_eq!(window.position_from(0, |aa| aa == '*'), Some(5));
    assert_eq!(window.iter().collect::<String>(), "MKLWF*");
}

#[test]
fn a_local_substitution_does_not_copy_the_peptide() {
    // SNV at codon 1500, CodingFastPaths { lazy_cds: true, peptide_view: true, indel_prefix: false }:
    // substitute_peptide_window returns AltPeptide::Window whose replacement has length 1.
}
```

- [ ] **Step 2: Run to verify they fail**

Run: `cargo test -p datafusion-bio-function-vep --lib alt_peptide a_local_substitution_does_not_copy`
Expected: compile error, `AltPeptide` not found.

- [ ] **Step 3: Implement**

`new_aas` becomes an `AltPeptide`; the fallback branch wraps its `Vec<char>` in `Owned`. Replace each slice or index on `new_aas` by the matching method. Where a consumer needs a contiguous `&[char]` for more than a window, call `to_vec()` there; do not change that consumer's logic. Without `fast.peptide_view`, convert the `Window` to `Owned` immediately.

- [ ] **Step 4: Run the tests**

Run: `cargo test -p datafusion-bio-function-vep --lib alt_peptide coding_fast_paths && cargo test -p datafusion-bio-function-vep`
Expected: PASS.

- [ ] **Step 5: Commit**

```bash
git commit -am "perf(vep): represent a substituted peptide as a window over the reference"
```

### Task 5: Indels reuse the reference prefix (if Task 1 lists it)

**Files:**
- Modify: `datafusion/bio-function-vep/src/transcript_consequence.rs` (the fallback branch at 10072-10077; `translated_alt_protein_for_hgvs` at 7420)

**Interfaces:**
- Consumes: `MutatedCds` (Task 3), `CodingFastPaths::indel_prefix` (Task 2), `CodingReference::aas`, `local_substitution_safe`.
- Produces:
  - `fn translate_alt_from_codon(table: CodonTable, reference_aas: &[char], mutated: &MutatedCds<'_>, first_codon: usize) -> Option<Vec<char>>` — `reference_aas[..first_codon]` followed by the translation of `mutated.tail(first_codon * 3)`; `None` when `first_codon > reference_aas.len()`.
  - `fn translated_alt_protein_for_hgvs_from(tx: &TranscriptFeature, reference_protein: &str, mutated: &MutatedCds<'_>, first_codon: usize, semantics: VepSemantics) -> Option<String>` — the same string `translated_alt_protein_for_hgvs(tx, mutated.full(), semantics)` returns, built from `reference_protein[..first_codon]` plus the translation of the remaining CDS (trimmed as that function trims it) and the 3' UTR.

- [ ] **Step 1: Write the failing tests**

```rust
#[test]
fn translate_alt_from_codon_equals_the_full_translation() {
    let (tx, _exons, tl) = long_transcript(1, 600, CodonTable::Standard);
    let cds = reference_translateable_seq_for_vep(&tx, Some(&tl)).unwrap();
    let reference: Vec<char> = translate_protein_from_cds(cds.as_bytes()).unwrap();
    for (start, end, alt) in [(900usize, 900usize, &b""[..]), (901, 900, b"A"), (300, 305, b"GG"), (1797, 1799, b"")] {
        let view = MutatedCds::new(cds.as_bytes(), start, end, alt);
        let full = translate_protein_from_cds(view.full()).unwrap();
        assert_eq!(translate_alt_from_codon(CodonTable::Standard, &reference, &view, start / 3).unwrap(), full);
    }
}

#[test]
fn hgvs_alternate_protein_from_prefix_equals_the_full_translation() {
    // Same four edits, both VepSemantics values: translated_alt_protein_for_hgvs_from(..)
    // equals translated_alt_protein_for_hgvs(&tx, view.full(), semantics).
}

#[test]
fn an_unsafe_reference_takes_the_full_translation() {
    // The transcript with an N in its CDS, a 1-base deletion, CodingFastPaths::ALL:
    // the cfg(test) counter in translate_alt_from_codon stays at 0 and the result
    // equals the CodingFastPaths::NONE result.
}
```

- [ ] **Step 2: Run to verify they fail**

Run: `cargo test -p datafusion-bio-function-vep --lib translate_alt_from_codon hgvs_alternate_protein_from_prefix an_unsafe_reference`
Expected: compile error, `translate_alt_from_codon` not found.

- [ ] **Step 3: Implement**

In the fallback branch, when `fast.indel_prefix` and the constraints of Global Constraints hold, call `translate_alt_from_codon(tx.codon_table, &old_aas, &mutated, start_idx / 3)`; otherwise the current call. For the frameshift and stop-lost closure, use `translated_alt_protein_for_hgvs_from` under the same conditions, with the reference protein text taken from the memo.

- [ ] **Step 4: Run the tests**

Run: `cargo test -p datafusion-bio-function-vep`
Expected: PASS, including `coding_fast_paths`.

- [ ] **Step 5: Commit**

```bash
git commit -am "perf(vep): translate indels from the first affected codon"
```

### Task 6: TTN block check and the gate

**Files:**
- Modify: `datafusion/bio-function-vep/PERFORMANCE.md` (results in "Long-CDS coding change")

**Interfaces:**
- Consumes: `build_venv.sh`, `run_arm.sh`, `gate.sh` (plan 0); `$G/p4/ttn.vcf.gz`, the baseline MD5 and timing (Task 1).

- [ ] **Step 1: Lint and unit tests**

Run: `cargo fmt --all -- --check && cargo clippy --all-targets --all-features -- -D warnings && cargo test --all`
Expected: all pass.

- [ ] **Step 2: TTN block**

Build `$G/p4/venv-candidate` from this checkout. Run Task 1's command three times with it, and once with `VEP_CODING_FAST_PATHS=0`.
Expected: all four body MD5s equal the baseline MD5; median user + sys at most 85% of the baseline median; the `VEP_CODING_FAST_PATHS=0` run within 3% of the baseline.

- [ ] **Step 3: Smoke**

```bash
$L/run_arm.sh "$G/p4/venv-candidate" "$G/p4/smoke" clinvar 8 parity --chroms "2 17 MT"
$L/run_arm.sh "$G/p4/venv-candidate" "$G/p4/smoke" hg002 8 parity --chroms chr22
```
Expected: `PASS` on every row.

- [ ] **Step 4: Full gate**

```bash
$L/gate.sh "$G/p4/functions-baseline" "$PWD" "$G/p4/gate" --target-pct HG002=-2,ClinVar=1
cat "$G/p4/gate/verdict.txt"
```
Expected: every rule `PASS`.

- [ ] **Step 5: Document, commit and open the PR**

Add to `PERFORMANCE.md`: the before/after TTN table (CPU, wall, the Task 1 shares re-measured), the gate totals and verdict, build identities, and which of Tasks 3-5 were implemented.

```bash
git commit -am "docs(vep): long-CDS coding change measurements"
```
