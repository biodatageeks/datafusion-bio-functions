# Local Parity and Performance Gate Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Three scripts that build vepyr against any functions checkout and run vepyr's e2e runner per chromosome on HG002 and ClinVar with the merged cache, producing one parity, identity and performance verdict.

**Architecture:** `build_venv.sh` builds a wheel into a throwaway venv from a dedicated vepyr worktree whose `Cargo.toml` points at a local functions checkout. `run_arm.sh` runs `e2e-testing/scripts/run_comparison.py` once per chromosome under `/usr/bin/time -l`, records one TSV row per chromosome and deletes the VCFs. `report.py` reads arm directories and applies the rules of the spec, section 4.

**Tech Stack:** bash, Python 3.12 standard library, `uv`, Rust 1.91.1, vepyr `e2e-testing/scripts/run_comparison.py`, `tabix`, `bgzip`.

**Spec:** `docs/superpowers/specs/2026-11-10-parallel-vcf-output-bottlenecks-design.md` (section 4)

## Global Constraints

- Branch `chore/local-gate` off `origin/master` of `datafusion-bio-functions`; scripts live in `datafusion/bio-function-vep/scripts/local_gate/`.
- The scripts never write inside the user's vepyr checkout. All state goes under `GATE_ROOT` (default `$HOME/research/benchmarks/local-gate`).
- Defaults: `VEPYR_REPO` = `<functions repo root>/../vepyr`, `VEPYR_REF` = `origin/master`, `DATA_VEPYR_DIR` = `$HOME/workspace/data_vepyr`.
- Build command, verbatim: `RUSTUP_TOOLCHAIN=1.91.1 RUSTFLAGS='-C target-cpu=native' CARGO_INCREMENTAL=0 uv sync --frozen --no-cache --no-default-groups --no-editable --reinstall-package vepyr`, with `CARGO_TARGET_DIR` and `UV_PROJECT_ENVIRONMENT` set per venv.
- Guards before every chromosome: at least 20 GiB free on the arm directory's volume; background CPU (`ps -A -o %cpu=` summed) under 200%, waited for up to 180 s.
- Every `VEP_*` and `MIMALLOC_*` variable inherited from the shell is unset before a run; only the `ENV=VALUE` pairs given to `run_arm.sh` are applied.
- Deletions use `"${VAR:?}/..."` targets only.

## Review Focus

- A run that crashes mid-chromosome must leave a row with a non-zero exit code, not a missing row: `report.py` fails on missing or failed rows (test in Task 3).
- A chromosome whose report JSON is missing `annotation.time_s` must fail the report, not count as 0 s (test in Task 3).
- BGZF arms must hash the decompressed body, or identity against plain arms is meaningless (Task 2, step 3).
- A stale vepyr build worktree with a leftover path dependency must be reset before each build (Task 1, step 2).
- Two arms with different chromosome sets must fail the comparison, not compare the intersection (test in Task 3).

---

### Task 1: `build_venv.sh`

**Files:**
- Create: `datafusion/bio-function-vep/scripts/local_gate/build_venv.sh`

**Interfaces:**
- Produces: `build_venv.sh <functions-checkout> <venv-dir>`; leaves `<venv-dir>/bin/python` with vepyr installed and `<venv-dir>/GATE_IDENTITY.json` with keys `functions_sha`, `functions_dirty` (bool), `vepyr_sha`, `extension_sha256`, `toolchain`, `rustflags`. The vepyr build worktree is `$GATE_ROOT/vepyr-build`; later tasks run the e2e scripts from it.

- [ ] **Step 1: Write the script header and argument checks**

`set -euo pipefail`; exit 2 with a usage line unless two arguments are given and `<functions-checkout>/datafusion/bio-function-vep/Cargo.toml` exists.

- [ ] **Step 2: Create or reset the vepyr build worktree**

If `$GATE_ROOT/vepyr-build` is missing: `GIT_LFS_SKIP_SMUDGE=1 git -C "$VEPYR_REPO" worktree add --detach "$GATE_ROOT/vepyr-build" "$VEPYR_REF"`. Otherwise: `git -C "$GATE_ROOT/vepyr-build" checkout -- Cargo.toml Cargo.lock` then `GIT_LFS_SKIP_SMUDGE=1 git -C "$GATE_ROOT/vepyr-build" checkout --detach "$VEPYR_REF"`.

- [ ] **Step 3: Point the dependency at the functions checkout**

Replace the single line starting with `datafusion-bio-function-vep = ` in the worktree's `Cargo.toml` by
`datafusion-bio-function-vep = { path = "<absolute functions-checkout>/datafusion/bio-function-vep", features = ["cache-builder"] }`.
Fail if the line count matching `^datafusion-bio-function-vep = ` is not exactly 1.

- [ ] **Step 4: Build and record identity**

Run the build command from Global Constraints in the worktree with `CARGO_TARGET_DIR="$GATE_ROOT/target/$(basename <venv-dir>)"` and `UV_PROJECT_ENVIRONMENT=<venv-dir>`, logging to `<venv-dir>.build.log`. Write `GATE_IDENTITY.json`; `extension_sha256` is the SHA-256 of the one `_core*.so` under `<venv-dir>/lib/python*/site-packages/vepyr/`.

- [ ] **Step 5: Verify**

Run:
```bash
git -C . worktree add --detach "$HOME/research/benchmarks/local-gate/functions-baseline" origin/master
datafusion/bio-function-vep/scripts/local_gate/build_venv.sh \
  "$HOME/research/benchmarks/local-gate/functions-baseline" \
  "$HOME/research/benchmarks/local-gate/venv-baseline"
"$HOME/research/benchmarks/local-gate/venv-baseline/bin/python" -P -c 'import vepyr; print(vepyr.__version__)'
```
Expected: a version string, and `GATE_IDENTITY.json` with `"functions_dirty": false`.

- [ ] **Step 6: Commit**

```bash
git add datafusion/bio-function-vep/scripts/local_gate/build_venv.sh
git commit -m "chore(vep): local gate build script"
```

### Task 2: `run_arm.sh`

**Files:**
- Create: `datafusion/bio-function-vep/scripts/local_gate/run_arm.sh`

**Interfaces:**
- Consumes: `<venv-dir>/bin/python`, `$GATE_ROOT/vepyr-build/e2e-testing/scripts/run_comparison.py` from Task 1.
- Produces: `run_arm.sh <venv-dir> <arm-dir> <hg002|clinvar> <workers> <parity|time> [--bgzf] [--chroms "<c1> <c2> ..."] [ENV=VALUE ...]`. Appends to `<arm-dir>/contigs.tsv` with header
  `dataset	contig	workers	mode	compression	exit	wall_s	rss_gib	time_s	records	body_md5	verdict	env`
  and writes `<arm-dir>/arm.json` (`venv` identity copied from `GATE_IDENTITY.json`, arguments, start time). It is resumable: a `(dataset, contig)` already present with `exit` 0 is skipped.

- [ ] **Step 1: Dataset table**

`hg002`: profile `merged`, suffix `merged`, `--vcf "$DATA_VEPYR_DIR/input/HG002_normalized.vcf.gz" --no-normalize`, md5 mode `both`, contigs `chr1` … `chr22`. `clinvar`: profile `clinvar_merged`, suffix `clinvar_merged`, `--no-normalize`, md5 mode `strict-clin-sig`, contigs `1` … `22 X Y MT`. `--chroms` overrides the contig list.

- [ ] **Step 2: Per-chromosome run**

After the guards, run from `$GATE_ROOT/vepyr-build/e2e-testing/scripts`:
```bash
/usr/bin/time -l env "${UNSET[@]}" DATA_VEPYR_DIR="$DATA_VEPYR_DIR" "${EXTRA_ENV[@]}" \
  "$VENV/bin/python" -P run_comparison.py --release 116 --profile "$PROFILE" \
  --chroms "$CONTIG" --workers "$WORKERS" --output-dir "$ARM" "${INPUT_ARGS[@]}" \
  "${MODE_ARGS[@]}" ${BGZF:+--bgzf}
```
`MODE_ARGS` is `--comparison-mode md5 --md5-mode <mode>` for `parity` and `--skip-compare` for `time`. Stdout and stderr go to `<arm-dir>/logs/<dataset>.<contig>.log`.

- [ ] **Step 3: Record the row, then delete the VCFs**

`wall_s` and `rss_gib` come from the `real` and `maximum resident set size` lines of `/usr/bin/time -l`. `time_s` is `annotation.time_s` and `records` is `annotation.output_variants` from `<arm-dir>/reports/fast_chr<N>_<suffix>_116_report.json` (`<N>` is the contig without a `chr` prefix); write `NA` when the file or key is missing. `body_md5` is the MD5 of the lines not starting with `#` of `<arm-dir>/results/116/fast_chr<N>/vepyr_parquet_chr<N>_<suffix>.vcf`, or of `bgzip -dc` of the `.vcf.gz` under `--bgzf`. `verdict` is `PASS` when every `md5 …: body` line in the log says `PASS`, `FAIL` when any says otherwise, and `-` in `time` mode. Then remove `"${ARM:?}/results/116/fast_chr$N"/vepyr_* "${ARM:?}/results/116/fast_chr$N"/vep_*`.

- [ ] **Step 4: Smoke run**

Run:
```bash
G="$HOME/research/benchmarks/local-gate"
datafusion/bio-function-vep/scripts/local_gate/run_arm.sh "$G/venv-baseline" "$G/smoke/parity" hg002 8 parity --chroms chr22
datafusion/bio-function-vep/scripts/local_gate/run_arm.sh "$G/venv-baseline" "$G/smoke/time1" hg002 1 time --chroms chr22
datafusion/bio-function-vep/scripts/local_gate/run_arm.sh "$G/venv-baseline" "$G/smoke/bgzf" hg002 8 time --bgzf --chroms chr22
cat "$G"/smoke/*/contigs.tsv
```
Expected: three rows with `exit` 0, `records` 50861, the same `body_md5` in all three, `verdict` `PASS` on the first, and no `vepyr_*` or `vep_*` file left under `$G/smoke`.

- [ ] **Step 5: Commit**

```bash
git add datafusion/bio-function-vep/scripts/local_gate/run_arm.sh
git commit -m "chore(vep): local gate per-chromosome e2e runner"
```

### Task 3: `report.py`

**Files:**
- Create: `datafusion/bio-function-vep/scripts/local_gate/report.py`
- Test: `datafusion/bio-function-vep/scripts/local_gate/test_report.py`

**Interfaces:**
- Consumes: `contigs.tsv` rows from Task 2.
- Produces: `report.py --baseline DIR [DIR ...] --candidate DIR [DIR ...] [--parity DIR ...] [--identity DIR ...] [--baseline-1w DIR --candidate-1w DIR] --target-pct HG002=<p>,ClinVar=<p>`; exit 0 when every rule holds, 1 otherwise, 2 on unusable input. Importable functions:
  `load_arm(path) -> list[dict]`, `totals(rows) -> dict[str, float]` (dataset to summed `time_s`),
  `check(baseline, candidate, parity, identity, baseline_1w, candidate_1w, targets) -> list[tuple[str, bool, str]]` (rule name, passed, detail).

- [ ] **Step 1: Write the failing tests**

In `test_report.py` (unittest, arm directories built in a temp dir by a helper `arm(rows)`):

```python
def test_meets_target_and_passes(self):
    base = [arm(rows("hg002", {"chr1": 10.0, "chr2": 8.0})) for _ in range(2)]
    cand = [arm(rows("hg002", {"chr1": 8.0, "chr2": 7.0})) for _ in range(2)]
    results = dict((n, ok) for n, ok, _ in check(base, cand, [], [], None, None, {"hg002": 10.0}))
    self.assertTrue(all(results.values()))

def test_target_missed(self):          # 18.0 -> 17.5 is 2.8%, target 10%
    ...; self.assertFalse(results["target hg002"])

def test_single_contig_regression(self):   # chr2 8.0 -> 9.0: +12.5% and +1.0 s
    ...; self.assertFalse(results["no contig regression hg002"])

def test_small_absolute_regression_is_allowed(self):  # chr2 1.0 -> 1.3: +30% but 0.3 s
    ...; self.assertTrue(results["no contig regression hg002"])

def test_rss_over_ten_percent_fails(self):  # max rss 4.0 -> 4.5
    ...; self.assertFalse(results["peak rss hg002"])

def test_identity_mismatch_fails(self):     # one arm has a different body_md5 for chr2
    ...; self.assertFalse(results["identity"])

def test_parity_fail_verdict_fails(self):
    ...; self.assertFalse(results["parity"])

def test_failed_or_missing_row_fails(self):  # exit 1 row; and a contig absent from one arm
    ...; self.assertFalse(results["complete"])

def test_missing_time_fails(self):           # time_s == "NA"
    ...; self.assertFalse(results["complete"])

def test_one_worker_regression_over_three_percent_fails(self):  # 100.0 -> 104.0
    ...; self.assertFalse(results["1 worker hg002"])
```

- [ ] **Step 2: Run the tests to verify they fail**

Run: `python3 -m unittest datafusion/bio-function-vep/scripts/local_gate/test_report.py -v`
Expected: FAIL with `ImportError` for `report`.

- [ ] **Step 3: Implement `report.py`**

Rules, by name: `complete` (every arm has the same contig set per dataset, every row `exit` 0 and numeric `time_s`); `parity` (every row of every `--parity` arm has `verdict` `PASS`); `identity` (per `(dataset, contig)`, one distinct `body_md5` across `--identity`, `--parity` and `--candidate` arms and the first `--baseline` arm); `target <dataset>` (mean candidate total at most `(1 - p/100)` of mean baseline total); `no contig regression <dataset>` (per contig, mean candidate `time_s` not above both `1.10 ×` and `+0.5 s` of mean baseline); `peak rss <dataset>` (max candidate `rss_gib` at most `1.10 ×` max baseline); `1 worker <dataset>` (candidate total at most `1.03 ×` baseline). Print a per-contig table (baseline mean, candidate mean, ratio) and one `PASS`/`FAIL` line per rule.

- [ ] **Step 4: Run the tests to verify they pass**

Run: `python3 -m unittest datafusion/bio-function-vep/scripts/local_gate/test_report.py -v`
Expected: 10 tests, OK.

- [ ] **Step 5: Commit**

```bash
git add datafusion/bio-function-vep/scripts/local_gate/report.py datafusion/bio-function-vep/scripts/local_gate/test_report.py
git commit -m "chore(vep): local gate report with parity, identity and performance rules"
```

### Task 4: `gate.sh` and the README

**Files:**
- Create: `datafusion/bio-function-vep/scripts/local_gate/gate.sh`
- Create: `datafusion/bio-function-vep/scripts/local_gate/README.md`

**Interfaces:**
- Consumes: Tasks 1-3.
- Produces: `gate.sh <baseline-functions-checkout> <candidate-functions-checkout> <gate-dir> --target-pct HG002=<p>,ClinVar=<p> [--bgzf-identity] [--candidate-env ENV=VALUE ...]`. It builds `venv-baseline` and `venv-candidate` under `<gate-dir>` and runs, in this order, for `hg002` then `clinvar`:
  1. `candidate` warm-up, 8 workers, `time` (discarded: arm `warmup`);
  2. `baseline` r1, `candidate` r1, `candidate` r2, `baseline` r2, 8 workers, `time`;
  3. `candidate` parity, 8 workers;
  4. `candidate` 4 workers and 1 worker, `baseline` 1 worker, `time`;
  5. with `--bgzf-identity`: `candidate` 8 workers, `time`, `--bgzf`;
  then `report.py` with the 1-worker arms, the 4-worker and BGZF arms as `--identity`, writing `<gate-dir>/verdict.txt`.

- [ ] **Step 1: Write `gate.sh`** as the sequence above; every arm directory is `<gate-dir>/<dataset>/<arm>`; the script is resumable because `run_arm.sh` is.

- [ ] **Step 2: Write the README**: purpose, the three defaults, one full `gate.sh` example, expected duration (about two hours), the 20 GiB disk guard, and the statement that the CI `/parity-wgs` run on the vepyr pin PR remains the release gate.

- [ ] **Step 3: Self-gate (baseline against itself)**

Run:
```bash
G="$HOME/research/benchmarks/local-gate"
datafusion/bio-function-vep/scripts/local_gate/gate.sh "$G/functions-baseline" "$G/functions-baseline" \
  "$G/self" --target-pct HG002=-3,ClinVar=-3 --bgzf-identity
```
Expected: `verdict.txt` with every rule `PASS` (a target of -3 allows 3% noise), 47 `PASS` parity rows, and identical `body_md5` across all arms. Record the two 8-worker totals per dataset; they are the run-to-run noise the later plans compare against.

- [ ] **Step 4: Commit and open the PR**

```bash
git add datafusion/bio-function-vep/scripts/local_gate/gate.sh datafusion/bio-function-vep/scripts/local_gate/README.md
git commit -m "chore(vep): local gate driver and README"
```
