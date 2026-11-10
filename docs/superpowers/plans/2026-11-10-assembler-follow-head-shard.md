# Assembler Follows the Head Shard Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** For plain and gzip VCF output, the assembler copies the bytes already written to the lowest unpublished shard while its run is still annotating, so no large copy is left when the last run finishes.

**Architecture:** Shard writers append through a `BufWriter<File>`, so every byte present in a shard file is final. `run_assembler_thread` already owns the output writer and receives one-shard `ContigShardRange` messages in id order (plan 1). Between messages it now reads whatever the next shard file holds and appends it; when that shard is published it drains the rest and deletes the file. BGZF assembly is untouched.

**Tech Stack:** Rust 2024, `std::sync::mpsc`, `std::fs`; all changes in `datafusion/bio-function-vep/src/vcf_sink.rs`.

**Spec:** `docs/superpowers/specs/2026-11-10-parallel-vcf-output-bottlenecks-design.md` (sections 1.2, 3.2, 4)

## Entry criterion

Read the residual-tail table that plan 1 (Task 6, step 7) wrote into `PERFORMANCE.md`. Implement this plan only if `tail / wall >= 0.03` on at least two of ClinVar contigs 1, 2 and 17. Otherwise add one sentence to that section saying the plan was not needed, with the three ratios, and stop.

## Global Constraints

- Branch `perf/assembler-follow-head` off `origin/master` after plan 1 is merged.
- Output bytes identical to plan 1's output, for plain, gzip and BGZF.
- `<output>.part` is discarded unless the driver set `finalize`, exactly as today.
- `VEP_ASSEMBLER_FOLLOW=0` restores copy-on-publish in the same binary.
- Gate targets (spec section 4), against `origin/master` with plan 1: 8-worker sum of `time_s` at least 3% lower on ClinVar and 2% lower on HG002.
- Line numbers refer to `68d0f7b` plus plan 1.

## Review Focus

- A shard published before the follower ever opened it (a fast run): copied whole on publish (test 2 in Task 1).
- A shard that is empty at publish: no bytes, file removed, next id followed (test 2).
- A run that fails after its shard was partly followed: `finalize` stays false and `<output>.part` is removed (test 3).
- A range covering several shards (the static path under `VEP_SHARD_RUN_POOL=0` still sends whole-contig ranges): all of them copied in id order (test 4).
- The follower must not spin: it waits on the channel with a timeout, and a busy loop would cost a core for the whole run (reviewer checks the wait; `FOLLOW_POLL` is 5 ms).

---

### Task 1: Follower

**Files:**
- Modify: `datafusion/bio-function-vep/src/vcf_sink.rs` (`run_assembler_thread` at 1478; tests from 2690)

**Interfaces:**
- Consumes: `ContigShardRange`, `write_vcf_body_chunk`, `VcfLocalWriter` (existing).
- Produces:
  - `const FOLLOW_POLL: std::time::Duration = std::time::Duration::from_millis(5);`
  - `fn assembler_follow_enabled(env: Option<&str>) -> bool` — `true` when unset. Variable: `VEP_ASSEMBLER_FOLLOW`.
  - `struct ShardFollower { tempdir: PathBuf, next_id: usize, open: Option<std::fs::File>, bytes: u64 }`
  - `impl ShardFollower { fn new(tempdir: PathBuf) -> Self; fn follow(&mut self, writer: &mut VcfLocalWriter, buffer: &mut [u8]) -> Result<u64>; fn publish(&mut self, range: &ContigShardRange, writer: &mut VcfLocalWriter, buffer: &mut [u8]) -> Result<u64> }`
    `follow` appends what `partition_{next_id}` currently holds past the bytes already taken and returns the count (0 when the file does not exist yet). `publish` finishes every id of the range in order (draining an open file, or copying an unopened one whole), removes each file and leaves `next_id == range.end_shard_id`.

- [ ] **Step 1: Write the failing tests**

All four use a temp dir, a `Plain` writer primed with a one-line header and the existing test helper that writes `partition_{id:04}.vcf.body`.

```rust
#[test]
fn follower_appends_bytes_before_the_shard_is_published() {
    // shard 0 holds "a\n" x 1000; follow() returns 2000 and the writer's file grew by 2000.
    // Append "b\n" x 1000 to the shard; publish(0..1) returns 2000; shard file is gone;
    // body == "a\n" x 1000 + "b\n" x 1000; follower.next_id == 1.
}

#[test]
fn follower_publish_copies_unopened_and_empty_shards() {
    // shards 0 ("x\n"), 1 (empty), 2 ("y\n"); no follow() call.
    // publish(0..1), publish(1..2), publish(2..3) return 2, 0, 2; body == "x\ny\n"; no shard file left.
}

#[test]
fn assembler_discards_a_followed_output_when_not_finalized() {
    // run_assembler_thread with shard 0 present; wait (<= 2 s) until <output>.part
    // is longer than the header; drop the sender with finalize == false.
    // The thread returns Ok(0); neither <output> nor <output>.part exists.
}

#[test]
fn follower_publish_handles_a_multi_shard_range() {
    // shards 0 ("a\n"), 1 ("b\n"), 2 ("c\n"); follow() once (takes shard 0's bytes);
    // publish(0..3) returns 4; body == "a\nb\nc\n"; next_id == 3.
}

#[test]
fn assembler_follow_defaults_on_and_zero_disables() {
    assert!(assembler_follow_enabled(None));
    assert!(!assembler_follow_enabled(Some("0")));
}
```

- [ ] **Step 2: Run the tests to verify they fail**

Run: `cargo test -p datafusion-bio-function-vep --lib follower assembler_follow assembler_discards_a_followed`
Expected: compile error, `ShardFollower` not found.

- [ ] **Step 3: Implement**

`run_assembler_thread` loop: `match rx.recv_timeout(FOLLOW_POLL)` — `Ok(range)`: `follower.publish(..)`, add `range.output_lines`, emit the `vcf_assemble range` trace; `Err(Timeout)`: `follower.follow(..)` when following is enabled; `Err(Disconnected)`: leave the loop and finalize exactly as today. With following disabled, use the blocking `rx.recv()` and `publish` only. `follow` and `publish` read with the thread's existing 8 MiB buffer.

- [ ] **Step 4: Run the tests**

Run: `cargo test -p datafusion-bio-function-vep --lib assembler follower`
Expected: PASS, including the unchanged `assembler_plain_appends_ranges_in_order_and_counts_lines`, `assembler_bgzf_raw_concat_decodes_in_order` and `assembler_discards_partial_output_when_not_finalized`.

- [ ] **Step 5: Commit**

```bash
git commit -am "perf(vep): assembler follows the head shard for plain output"
```

### Task 2: Smoke and gate

**Files:**
- Modify: `datafusion/bio-function-vep/PERFORMANCE.md` (extend "Shard run pool" with "Assembler follow")

**Interfaces:**
- Consumes: `build_venv.sh`, `run_arm.sh`, `gate.sh` (plan 0).

- [ ] **Step 1: Lint and unit tests**

Run: `cargo fmt --all -- --check && cargo clippy --all-targets --all-features -- -D warnings && cargo test --all`
Expected: all pass.

- [ ] **Step 2: Smoke, including gzip**

```bash
G="$HOME/research/benchmarks/local-gate"; L=datafusion/bio-function-vep/scripts/local_gate
git worktree add --detach "$G/p2/functions-baseline" origin/master
$L/build_venv.sh "$PWD" "$G/p2/venv-candidate"
$L/run_arm.sh "$G/p2/venv-candidate" "$G/p2/smoke/follow" hg002 8 parity --chroms chr22
$L/run_arm.sh "$G/p2/venv-candidate" "$G/p2/smoke/nofollow" hg002 8 time --chroms chr22 VEP_ASSEMBLER_FOLLOW=0
$L/run_arm.sh "$G/p2/venv-candidate" "$G/p2/smoke/follow" clinvar 8 parity --chroms "21 Y MT"
```
Then one gzip run: annotate `tabix -h … HG002_normalized.vcf.gz chr22` with `compression="gzip"` and compare `gzip -dc out.vcf.gz | grep -v '^#' | md5` with the `follow` arm's `body_md5`.
Expected: `PASS` verdicts; one `body_md5` per contig across `follow`, `nofollow` and gzip.

- [ ] **Step 3: Full gate**

```bash
$L/gate.sh "$G/p2/functions-baseline" "$PWD" "$G/p2/gate" --target-pct HG002=2,ClinVar=3 --bgzf-identity
cat "$G/p2/gate/verdict.txt"
```
Expected: every rule `PASS`.

- [ ] **Step 4: Residual tail after the change**

Repeat plan 1's Task 6, step 7 on ClinVar contigs 1, 2 and 17 and add the before/after `tail / wall` table to `PERFORMANCE.md`, with both build identities and the gate totals.

- [ ] **Step 5: Commit and open the PR**

```bash
git commit -am "docs(vep): assembler follow measurements"
```
