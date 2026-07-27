# Porting an Ensembl-VEP test — the contract

This directory is the audit trail for the Ensembl-VEP → vepyr test port. Everything in it
exists to answer one question mechanically: **did any Perl assertion get dropped on the way
across?**

A ported test that passes proves nothing about faithfulness. It proves the Rust code you
*did* write works. It says nothing about the Perl subtests you did not notice. That gap is
what the gate in [`../port_ledger.rs`](../port_ledger.rs) closes, and this README is its
contract.

Design rationale:
`porting-tests/docs/superpowers/specs/2026-07-27-test-port-forward-port-to-master-design.md`.

---

## What a port is

**One Perl `.t` file.** That is the audit object and the unit of review:

```
tests/port/perl/<Test>.t              vendored upstream source, byte-for-byte
tests/port/perl/<Test>.t.provenance   where it came from, which ref, which sha256
tests/port/<Test>.ledger.toml         one [[row]] per Perl assertion, classified
tests/port_<something>.rs             the Rust analogues the ledger's rows name
```

**One PR per `.t`.** A reviewer opening that PR sees the complete Perl source, the complete
classification, and every Rust assertion derived from it. Nothing about faithfulness
requires opening a second PR.

Worked example: [`Utils.ledger.toml`](Utils.ledger.toml) (44 rows) with
[`../port_utils.rs`](../port_utils.rs).

## How to add one

1. **Vendor the `.t` verbatim.** Byte-for-byte — every ledger row's `line` refers to the
   upstream file, so a single added header comment would shift all of them.

   ```bash
   cp …/ensembl-vep/t/<Test>.t datafusion/bio-function-vep/tests/port/perl/<Test>.t
   shasum -a 256 datafusion/bio-function-vep/tests/port/perl/<Test>.t
   ```

2. **Write the `.provenance` sidecar** next to it: upstream URL, path, git ref + commit,
   sha256, vendoring date, licence. Provenance goes *beside* the file, never inside it.
   The root `NOTICE` already carries the Apache-2.0 attribution for this whole directory;
   extend it only if you vendor from a new project.

3. **Get the machine numbering.** Add a `[[row]]` for `n = 1` only, set
   `expected_assertions` to a deliberate lie (e.g. `1`), and run the gate. It prints one
   line per unaccounted assertion, with its ordinal, its function and its line:

   ```bash
   cargo test -p datafusion-bio-function-vep --test port_ledger --offline
   ```

   Let those failures drive the ledger. Do not hand-count the `.t`; the scanner's numbering
   *is* the numbering, and disagreeing with it is a bug report against the scanner.

4. **Classify every row**, then write the Rust analogues the `unit-port` rows name.

5. **Green means row-complete, not done.** Read §"What the gate does not check" before you
   call it finished.

## The three classes

Every row carries `n`, `line`, `desc` and `class`. What else it needs depends on the class:

| `class` | also requires | means |
|---|---|---|
| `unit-port` | `rust = "<path from crate root>::<fn>"` | a Rust test covers this assertion |
| `architectural-no-analogue` | `reason = "…"` | vepyr cannot have this by design |
| `blocked-future-work` | `issue = "biodatageeks/vepyr#NN"` | it should exist; it does not yet |

`rust` is resolved against a textual index of every `#[test]` / `#[tokio::test]` fn in
`src/` and `tests/`. A comment claiming coverage does not create an index entry, so the
claim cannot be faked; nor can a `#[ignore]`d test stand in for a real one unless the
attribute carries both an issue reference and a human reason
(`#[ignore = "vepyr#42: needs a 2 GB fixture"]`).

**`reason` is prose a reviewer has to agree with.** "Not applicable", "N/A" and "no
analogue" are not reasons — they restate the class. A reason says *what the Perl subtest
pins down* and *why no vepyr code path can exhibit it*, and points at the analysis in
`porting-tests/detailed_plans/<Test>.md`. The gate only checks that the field is non-empty;
the field's job is to be argued with in review. An unjustified claim of impossibility is a
coverage gap wearing a label.

**Reclassifying down is correct; deleting is not.** If a Rust analogue cannot be
forward-ported because the API it targeted is gone, or if it turns out the existing test
does not actually cover the case, move the row to `blocked-future-work` and file the issue.
Rows 3 and 6 of `Utils.ledger.toml` are the live example: the in-`src` test they pointed at
asserts four of `format_coords`' six cases, so two rows were demoted rather than left
claiming coverage that is not there. A silently dropped row is the exact failure this gate
exists to prevent.

## Why `expected_assertions` exists

It is a **checksum on the scanner**, not a target. It is checked twice — against what the
enumerator found in the vendored Perl, and against how many `[[row]]` entries the ledger
carries.

If it ever mismatches, exactly one of two things happened:

- the **vendored file changed** — re-audit it, and update the `.provenance` sha256; or
- the **scanner mis-read it** — fix the scanner.

**Never just update the number.** Doing so converts a detected coverage hole into an
undetected one, which is strictly worse than not having the gate at all.

## What the gate does **not** check

It verifies that an assertion is **accounted for**. It never verifies that an assertion is
**strong**.

A `unit-port` row is satisfied by a `#[test] fn` that exists, is uniquely named, and is not
mutely `#[ignore]`d. That fn may assert nothing at all. It may assert something true of the
Rust standard library rather than of vepyr — `Utils.ledger.toml` rows 8 and 9 are exactly
that case, and say so in the test file: vepyr has no central `convert_arrayref` analogue,
so the tests pin the separator contract its call sites rely on rather than vepyr code. They
are honest, and they are also nearly vacuous. The gate cannot tell the difference.

Two things close that hole, and neither is CI:

- **Review.** The `desc` and `reason` fields exist so a human can disagree. A ledger is an
  argument, not a receipt.
- **The periodic tier on ii-hpc** — `cargo-mutants` scoped to the components a port
  declares it covers, plus an advisory LLM-as-judge over the
  `(Perl subtest, ledger row, Rust test)` triple, plus the real-VEP oracle cross-check for
  CSQ-observable ports. See §5 of the design spec. A port that passes the gate but leaves
  mutants alive in its owning component is a hollow test.

So: green is necessary, never sufficient.

## Two things that will bite you

**A row counts a *statement*, not a runtime TAP result.** The scanner's rule is "an
assertion call at a statement boundary", which is what makes multi-line argument lists safe.
It also means one statement inside a loop emits several TAP results but earns exactly one
row. Three sites in the corpus do this:

| site | statements | TAP results | why |
|---|---:|---:|---|
| `AnnotationSource.t:97` | 1 | 7 | `is_deeply(…) for @$inputs` |
| `Config.t:146` | 1 | 3 | assertion inside a `foreach` |
| `AnnotationSource_File.t:110` | 1 | 5 | `is()` inside a 5-format `foreach` |

That is not a defect, but the row's `desc` must say so — write "…for each of the 7 inputs",
not just the first case, or the next reader will think six results vanished.

**The mapping is many-to-many.** One Rust test file may serve rows from several ledgers, and
one ledger's rows may point into several files (`Utils.ledger.toml` names both
`src/transcript_consequence.rs` and `tests/port_utils.rs`). Both are fine. The audit unit is
the `.t`, not the `.rs`. The gate does enforce the reverse direction: every
`tests/port_*.rs` must be named by at least one ledger row, so a port cannot land by simply
not writing a ledger.

## The corpus sweep

Beyond the vendored files, the enumerator is pinned against **all 49 upstream `.t` files**,
with each per-file count cross-checked against the hand-written tables in
`porting-tests/detailed_plans/` and every disagreement adjudicated in writing. It reads an
`ensembl-vep` checkout, so it is `#[ignore]`d and CI stays offline:

```bash
cargo test -p datafusion-bio-function-vep --test port_ledger --offline \
  -- --ignored --nocapture corpus::
```

Those counts are **frozen**. If a change to the lexer moves one, this test fails — and that
is the point: the scanner's one intolerable failure mode is silently under-counting, because
an assertion it never sees is never required to have a ledger row, and the gate then reports
success while coverage is incomplete.

So if a count moves: **prove which side is right before editing the table.** Read the Perl
at the disputed line, decide whether the old or the new count is correct, and record the
adjudication in the entry's `verdict` string. Re-baselining the table to make a build pass
disarms the only defence against a scanner that has started lying.

## Checklist before opening the PR

```bash
cargo test  -p datafusion-bio-function-vep --test port_ledger --offline   # gate green
cargo test  -p datafusion-bio-function-vep --offline --no-fail-fast       # suite green
cargo clippy -p datafusion-bio-function-vep --all-targets --offline
cargo fmt --check
```

- [ ] `.t` vendored byte-for-byte, `.provenance` sidecar written, sha256 recorded in both
      the sidecar and the ledger.
- [ ] `expected_assertions` equals what the scanner found — and you did not edit it to make
      that true.
- [ ] Every `architectural-no-analogue` `reason` is an argument, not a restatement.
- [ ] Every `blocked-future-work` `issue` is a real, filed issue number.
- [ ] Every `desc` on a loop-backed row says how many cases it stands for.
- [ ] You deliberately broke the ledger once and watched the gate fail. A gate nobody
      watched bite is not known to work.
