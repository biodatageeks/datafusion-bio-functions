# Regenerate Variation Cache into chr1 Bundled Layout — Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Regenerate all 463 per-contig Lance datasets under `variation.lance/` into the uniform chr1 bundled layout (21 cols: `tier` + 3 scalar-`Utf8` AF columns) so the whole cache reads as one polars `LazyFrame`.

**Architecture:** A bash driver loops the existing `build_lance_variation_chrom` release binary over every contig, in-place with `--overwrite`. No Rust changes — the builder already emits the target layout. Python helpers handle schema-kind detection, a pre-flight row-count snapshot, and post-run validation (schema uniformity + single-LazyFrame + row-count parity).

**Tech Stack:** Rust (existing example binary, features `lance-cache,cache-builder`), Bash driver, Python 3 (`lance`, `polars`, `pyarrow`) for snapshot/validation.

## Global Constraints

- Cache root (source): `/Users/mwiewior/workspace/data_vepyr/homo_sapiens_merged/115_GRCh38`
- Output dir: `/Users/mwiewior/workspace/data_vepyr/115_GRCh38_merged`
- Entity dir: `<output>/variation.lance` (holds `<contig>.lance` subdirs + `chrom_manifest` file)
- `--cache-source-type merged`, `--partitions 8`
- **In-place regeneration only** — `variation.lance` is 70 GB; volume has ~52 GB free. No side-by-side copy. Peak extra disk ≈ one contig.
- Target schema (chr1, 21 cols): `chrom, start, end, allele_string, failed, variation_name, clin_sig, clin_sig_allele, clinical_impact, phenotype_or_disease, pubmed, somatic, minor_allele, minor_allele_freq, clinvar_ids, cosmic_ids, dbsnp_ids, tier (int8), af_global (Utf8), af_gnomade (Utf8), af_gnomadg (Utf8)`.
- New-layout detector: a dataset is "new" iff its schema contains `tier` AND `af_global`.
- Naming: pass `--chrom <name>` where `<name>` = the existing dataset dir name minus `.lance` (round-trips for all 463).
- No builder code changes. Legacy `build_lance_variation_cache.rs` untouched.

---

### Task 1: Schema-kind detector helper

A tiny Python helper the driver uses to decide skip-vs-rebuild and that validation reuses.

**Files:**
- Create: `scripts/variation_cache_schema_kind.py`

**Interfaces:**
- Produces CLI: `python scripts/variation_cache_schema_kind.py <dataset.lance>` → prints exactly one of `new` / `old` / `missing` to stdout; exit 0 for `new`/`old`, exit 2 for `missing`.

- [ ] **Step 1: Write the helper**

```python
#!/usr/bin/env python3
"""Classify a variation Lance dataset's schema as new (bundled chr1 layout),
old (exploded population AF columns), or missing."""
import sys
import lance

def kind(path: str) -> str:
    try:
        ds = lance.dataset(path)
    except Exception:
        return "missing"
    names = {f.name for f in ds.schema}
    if "tier" in names and "af_global" in names:
        return "new"
    return "old"

if __name__ == "__main__":
    if len(sys.argv) != 2:
        print("usage: variation_cache_schema_kind.py <dataset.lance>", file=sys.stderr)
        sys.exit(64)
    k = kind(sys.argv[1])
    print(k)
    sys.exit(0 if k in ("new", "old") else 2)
```

- [ ] **Step 2: Verify against the two known datasets**

Run:
```bash
python scripts/variation_cache_schema_kind.py /Users/mwiewior/workspace/data_vepyr/115_GRCh38_merged/variation.lance/chr1.lance
python scripts/variation_cache_schema_kind.py /Users/mwiewior/workspace/data_vepyr/115_GRCh38_merged/variation.lance/chrMT.lance
```
Expected: first prints `new`, second prints `old`.

- [ ] **Step 3: Commit**

```bash
git add scripts/variation_cache_schema_kind.py
git commit -m "feat: variation cache schema-kind detector"
```

---

### Task 2: Pre-flight row-count snapshot

Capture every contig's current row count BEFORE any overwrite, so we can prove no data loss afterward.

**Files:**
- Create: `scripts/variation_cache_snapshot_counts.py`

**Interfaces:**
- Produces CLI: `python scripts/variation_cache_snapshot_counts.py <entity_dir> <out.json>` → writes `{ "<contig>": <rows>, ... }` for every `*.lance` subdir, prints `N contigs, total rows = M`.

- [ ] **Step 1: Write the snapshot script**

```python
#!/usr/bin/env python3
"""Snapshot per-contig row counts of every dataset in a variation entity dir."""
import json
import sys
from pathlib import Path
import lance

def main(entity_dir: str, out_path: str) -> None:
    root = Path(entity_dir)
    counts = {}
    for ds_dir in sorted(root.glob("*.lance")):
        name = ds_dir.name[: -len(".lance")]
        counts[name] = lance.dataset(str(ds_dir)).count_rows()
    Path(out_path).write_text(json.dumps(counts, indent=2, sort_keys=True))
    total = sum(counts.values())
    print(f"{len(counts)} contigs, total rows = {total}")

if __name__ == "__main__":
    if len(sys.argv) != 3:
        print("usage: variation_cache_snapshot_counts.py <entity_dir> <out.json>", file=sys.stderr)
        sys.exit(64)
    main(sys.argv[1], sys.argv[2])
```

- [ ] **Step 2: Run the snapshot against production (this is the real pre-flight)**

Run:
```bash
python scripts/variation_cache_snapshot_counts.py \
  /Users/mwiewior/workspace/data_vepyr/115_GRCh38_merged/variation.lance \
  /Users/mwiewior/workspace/data_vepyr/115_GRCh38_merged/variation_counts_before.json
```
Expected: prints `463 contigs, total rows = <M>` and writes the JSON. Confirm the file exists and `chr1` is present with ~88153966.

- [ ] **Step 3: Record the baseline total disk size**

Run:
```bash
du -sh /Users/mwiewior/workspace/data_vepyr/115_GRCh38_merged/variation.lance
```
Expected: ~70G. Note the number in `regen_progress.log` later.

- [ ] **Step 4: Commit the script (NOT the snapshot JSON — it lives in the data dir, outside the repo)**

```bash
git add scripts/variation_cache_snapshot_counts.py
git commit -m "feat: variation cache pre-flight row-count snapshot"
```

---

### Task 3: Driver script

The core deliverable: enumerate contigs, skip already-migrated, rebuild each in-place with `--overwrite`, checkpoint to a log.

**Files:**
- Create: `scripts/regenerate_variation_cache.sh`

**Interfaces:**
- Consumes: `scripts/variation_cache_schema_kind.py` (Task 1).
- Produces CLI: `scripts/regenerate_variation_cache.sh [contig ...]`. With no args, processes all `*.lance`. With args, processes only the named contigs (used by the dry-run in Task 4). Env: `FORCE=1` rebuilds even already-`new` datasets. Appends one line per contig to `<output>/regen_progress.log`.

- [ ] **Step 1: Write the driver**

```bash
#!/usr/bin/env bash
set -euo pipefail

CACHE_ROOT="/Users/mwiewior/workspace/data_vepyr/homo_sapiens_merged/115_GRCh38"
OUTPUT_DIR="/Users/mwiewior/workspace/data_vepyr/115_GRCh38_merged"
ENTITY_DIR="$OUTPUT_DIR/variation.lance"
BIN="./target/release/examples/build_lance_variation_chrom"
KIND="python scripts/variation_cache_schema_kind.py"
LOG="$OUTPUT_DIR/regen_progress.log"
PARTITIONS=8
FORCE="${FORCE:-0}"

[ -x "$BIN" ] || { echo "binary not built: $BIN (run cargo build first)" >&2; exit 1; }

if [ "$#" -gt 0 ]; then
  CONTIGS=("$@")
else
  mapfile -t CONTIGS < <(find "$ENTITY_DIR" -maxdepth 1 -name '*.lance' -type d \
    -exec basename {} .lance \; | sort)
fi

echo "=== regen start $(date) contigs=${#CONTIGS[@]} force=$FORCE ===" | tee -a "$LOG"
fail=0
for contig in "${CONTIGS[@]}"; do
  ds="$ENTITY_DIR/$contig.lance"
  if [ "$FORCE" != "1" ] && [ "$($KIND "$ds" 2>/dev/null || true)" = "new" ]; then
    echo "$contig SKIP already-new" | tee -a "$LOG"
    continue
  fi
  start=$(date +%s)
  if "$BIN" --cache-root "$CACHE_ROOT" --output-dir "$OUTPUT_DIR" \
       --chrom "$contig" --cache-source-type merged \
       --partitions "$PARTITIONS" --overwrite >>"$LOG" 2>&1; then
    elapsed=$(( $(date +%s) - start ))
    rows=$($KIND "$ds" >/dev/null 2>&1 && python -c "import lance;print(lance.dataset('$ds').count_rows())")
    echo "$contig OK rows=$rows elapsed=${elapsed}s" | tee -a "$LOG"
  else
    echo "$contig FAIL (see $LOG)" | tee -a "$LOG"
    fail=$((fail+1))
  fi
done
echo "=== regen done $(date) failures=$fail ===" | tee -a "$LOG"
[ "$fail" -eq 0 ]
```

- [ ] **Step 2: Make executable + shellcheck (if available)**

Run:
```bash
chmod +x scripts/regenerate_variation_cache.sh
command -v shellcheck >/dev/null && shellcheck scripts/regenerate_variation_cache.sh || echo "shellcheck not installed, skipping"
```
Expected: executable bit set; shellcheck (if present) reports no errors.

- [ ] **Step 3: Build the binary (prerequisite for any run)**

Run:
```bash
cargo build --release --example build_lance_variation_chrom --features lance-cache,cache-builder
```
Expected: compiles clean; `./target/release/examples/build_lance_variation_chrom` exists.

- [ ] **Step 4: Smoke-test the SKIP path on chr1 (no rebuild — chr1 is already new)**

Run:
```bash
scripts/regenerate_variation_cache.sh chr1
```
Expected: log line `chr1 SKIP already-new`; chr1.lance untouched (still `new`, ~88153966 rows). This proves the resume/skip guard works without overwriting anything.

- [ ] **Step 5: Commit**

```bash
git add scripts/regenerate_variation_cache.sh
git commit -m "feat: in-place variation cache regeneration driver"
```

---

### Task 4: Staged dry-run on 3 representative contigs + e2e parity gate

Validate end-to-end on a small, fast, representative sample BEFORE the full sweep. Pick one patch, `chrMT` (tiny main, non-numeric mapping), and `chr21` (real-size main).

**Files:**
- Modify (data, not repo): overwrites `chrMT.lance`, `chr21.lance`, and one patch dataset in place.

**Interfaces:**
- Consumes: Task 3 driver, Task 1 detector.

- [ ] **Step 1: Pick the patch contig and capture the 3 baseline counts**

Run:
```bash
PATCH=$(find /Users/mwiewior/workspace/data_vepyr/115_GRCh38_merged/variation.lance \
  -maxdepth 1 -name '*PATCH.lance' -type d -exec basename {} .lance \; | sort | head -1)
echo "patch=$PATCH"
for c in "$PATCH" chrMT chr21; do
  python -c "import lance;print('$c', lance.dataset('/Users/mwiewior/workspace/data_vepyr/115_GRCh38_merged/variation.lance/$c.lance').count_rows())"
done
```
Expected: prints the chosen patch name and the three current (old-layout) row counts. Record them.

- [ ] **Step 2: Regenerate the 3 contigs**

Run:
```bash
scripts/regenerate_variation_cache.sh "$PATCH" chrMT chr21
```
Expected: three `OK rows=...` log lines; non-zero exit only on failure.

- [ ] **Step 3: Verify each is now the new layout with unchanged row count**

Run:
```bash
for c in "$PATCH" chrMT chr21; do
  k=$(python scripts/variation_cache_schema_kind.py /Users/mwiewior/workspace/data_vepyr/115_GRCh38_merged/variation.lance/$c.lance)
  n=$(python -c "import lance;print(lance.dataset('/Users/mwiewior/workspace/data_vepyr/115_GRCh38_merged/variation.lance/$c.lance').count_rows())")
  echo "$c kind=$k rows=$n"
done
```
Expected: each `kind=new`; each `rows` equals its Step 1 baseline.

- [ ] **Step 4: e2e annotation parity on chr21 against VEP golden**

Run:
```bash
cd /Users/mwiewior/research/git/vepyr && env -u VIRTUAL_ENV -u CONDA_PREFIX \
  uv run python e2e-testing/scripts/run_annotation_fast.py chr21 --cache merged --backend lance
```
Expected: 0 CSQ/AF mismatches (same gate chr1 already passes). This confirms bundled→CSQ unbundling is correct for a freshly-regenerated contig. If it fails, STOP — do not run the full sweep; investigate unbundle/build before proceeding.

- [ ] **Step 5: Record the dry-run outcome (no code commit — data + log only)**

Append a note to `regen_progress.log`:
```bash
echo "=== dry-run OK: $PATCH, chrMT, chr21 regenerated, parity green $(date) ===" \
  >> /Users/mwiewior/workspace/data_vepyr/115_GRCh38_merged/regen_progress.log
```
Expected: line appended. No git commit (the scripts are already committed; only data changed).

---

### Task 5: Full sweep over all remaining contigs

With the dry-run green, regenerate everything. Already-new contigs (chr1 + the 3 dry-run ones) are skipped automatically.

**Files:**
- Modify (data, not repo): overwrites the remaining ~459 datasets in place.

- [ ] **Step 1: Launch the full sweep (long-running — run in background)**

Run:
```bash
scripts/regenerate_variation_cache.sh
```
Expected: `=== regen start ... contigs=463 ... ===`, a stream of `OK`/`SKIP` lines, ending `=== regen done ... failures=0 ===`. (Several hours; chr1-sized contigs dominate. Patches are seconds each.)

- [ ] **Step 2: Confirm zero failures**

Run:
```bash
grep -c ' FAIL' /Users/mwiewior/workspace/data_vepyr/115_GRCh38_merged/regen_progress.log
tail -3 /Users/mwiewior/workspace/data_vepyr/115_GRCh38_merged/regen_progress.log
```
Expected: `FAIL` count is `0`; last line shows `failures=0`. If any failed, re-run the driver (it resumes by skipping already-new datasets and retrying the rest).

- [ ] **Step 3: Record post-run disk size**

Run:
```bash
du -sh /Users/mwiewior/workspace/data_vepyr/115_GRCh38_merged/variation.lance
```
Expected: ≤ 70G (bundled layout is smaller than exploded; net should drop). Note it in the log.

---

### Task 6: Validation — schema uniformity, single LazyFrame, row-count parity

Final gate proving the original goal is met. Also delivers the reusable single-LazyFrame reader.

**Files:**
- Create: `scripts/variation_lazyframe.py` (reusable reader — the original ask)
- Create: `scripts/variation_cache_validate.py` (the gate)

**Interfaces:**
- `variation_lazyframe.py` produces: `scan_all_variation(entity_dir: str) -> polars.LazyFrame` — vertical concat of all `*.lance` datasets.
- `variation_cache_validate.py` consumes `scan_all_variation` and `variation_counts_before.json` (Task 2).

- [ ] **Step 1: Write the reusable single-LazyFrame reader**

```python
#!/usr/bin/env python3
"""Read every per-contig variation Lance dataset as one polars LazyFrame."""
from pathlib import Path
import lance
import polars as pl

def scan_all_variation(entity_dir: str) -> pl.LazyFrame:
    frames = []
    for ds_dir in sorted(Path(entity_dir).glob("*.lance")):
        ds = lance.dataset(str(ds_dir))
        frames.append(pl.scan_pyarrow_dataset(ds.scanner()))
    if not frames:
        raise SystemExit(f"no .lance datasets under {entity_dir}")
    return pl.concat(frames, how="vertical")

if __name__ == "__main__":
    import sys
    lf = scan_all_variation(sys.argv[1])
    print(lf.select(pl.len()).collect().item())
```

- [ ] **Step 2: Write the validation gate**

```python
#!/usr/bin/env python3
"""Validate the regenerated variation cache: uniform schema, single LazyFrame,
row-count parity vs the pre-flight snapshot."""
import json
import sys
from pathlib import Path
import lance
import polars as pl
from variation_lazyframe import scan_all_variation

TARGET = {
    "chrom", "start", "end", "allele_string", "failed", "variation_name",
    "clin_sig", "clin_sig_allele", "clinical_impact", "phenotype_or_disease",
    "pubmed", "somatic", "minor_allele", "minor_allele_freq", "clinvar_ids",
    "cosmic_ids", "dbsnp_ids", "tier", "af_global", "af_gnomade", "af_gnomadg",
}

def main(entity_dir: str, before_json: str) -> None:
    root = Path(entity_dir)
    before = json.loads(Path(before_json).read_text())
    errors = []

    # 1. schema uniformity
    per_contig_rows = {}
    for ds_dir in sorted(root.glob("*.lance")):
        name = ds_dir.name[: -len(".lance")]
        ds = lance.dataset(str(ds_dir))
        names = {f.name for f in ds.schema}
        if names != TARGET:
            errors.append(f"{name}: schema {sorted(names ^ TARGET)} differs from target")
        per_contig_rows[name] = ds.count_rows()

    # 2. row-count parity vs snapshot
    for name, n_before in before.items():
        n_after = per_contig_rows.get(name)
        if n_after != n_before:
            errors.append(f"{name}: rows {n_before} -> {n_after} (mismatch)")

    # 3. single LazyFrame builds and total matches sum of per-dataset counts
    total_lf = scan_all_variation(entity_dir).select(pl.len()).collect().item()
    total_sum = sum(per_contig_rows.values())
    if total_lf != total_sum:
        errors.append(f"LazyFrame total {total_lf} != sum-of-counts {total_sum}")

    if errors:
        print("FAIL:")
        for e in errors:
            print("  -", e)
        sys.exit(1)
    print(f"OK: {len(per_contig_rows)} contigs, uniform schema, total rows = {total_sum}")

if __name__ == "__main__":
    main(sys.argv[1], sys.argv[2])
```

- [ ] **Step 3: Run the gate**

Run:
```bash
cd scripts && python variation_cache_validate.py \
  /Users/mwiewior/workspace/data_vepyr/115_GRCh38_merged/variation.lance \
  /Users/mwiewior/workspace/data_vepyr/115_GRCh38_merged/variation_counts_before.json
```
Expected: `OK: 463 contigs, uniform schema, total rows = <M>` (M equals the Task 2 baseline total). Any `FAIL:` line lists the offending contig — fix by re-running the driver on it (`scripts/regenerate_variation_cache.sh <contig>` with `FORCE=1`).

- [ ] **Step 4: Commit the reader + validator**

```bash
git add scripts/variation_lazyframe.py scripts/variation_cache_validate.py
git commit -m "feat: variation cache validation + single-LazyFrame reader"
```

---

## Self-Review

**Spec coverage:**
- Driver (Approach A, in-place, per-contig overwrite, manifest upsert via binary, checkpoint log) → Task 3. ✓
- Skip-migrated / resume → Task 1 detector + Task 3 guard. ✓
- Disk-safety (in-place only) → enforced by design; baseline+post size recorded (Tasks 2, 5). ✓
- Staged dry-run (patch + chrMT + chr21) before full run → Task 4. ✓
- e2e annotation parity → Task 4 Step 4. ✓
- Row-count parity → Task 2 snapshot + Task 6 check. ✓
- Schema uniformity → Task 6. ✓
- Single LazyFrame (original ask) → Task 6 `variation_lazyframe.py`. ✓
- Non-goals (no builder changes, AF-size gate excluded, empty contigs excluded, legacy builder untouched) → respected; no task touches them. ✓

**Placeholder scan:** No TBD/TODO; every script step shows full content; commands have expected output. ✓

**Type consistency:** `scan_all_variation(entity_dir)` defined in Task 6 Step 1, consumed in Step 2 with the same name/signature. `variation_cache_schema_kind.py` CLI contract (`new`/`old`/`missing`) used identically in Tasks 1, 3, 4. `variation_counts_before.json` produced in Task 2, consumed in Task 6. ✓

**Note on naming:** the validator imports `variation_lazyframe` as a sibling module, so Task 6 Step 3 runs from inside `scripts/` (`cd scripts && python …`). Consistent across the plan.
