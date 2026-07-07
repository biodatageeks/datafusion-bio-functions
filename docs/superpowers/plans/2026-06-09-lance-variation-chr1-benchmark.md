# Lance Variation chr1 Benchmark Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Build a reproducible chr1 benchmark/report harness for the approved Lance 2.1/2.2 variation-cache layout.

**Architecture:** Add a focused Python benchmark script under `research/` that materializes warm and cold Lance datasets from existing chr1 Parquet, supports Lance 2.1 top-level and Lance 2.2 packed variants, benchmarks warm full scans and cold 2,000-key point lookups, and writes JSON plus Markdown summary reports. Add unit tests for the helper logic that chooses logical columns, maps packed projections, and builds report rows.

**Tech Stack:** Python 3, PyArrow, pylance/lance Python package through `uv run --with pylance --with pyarrow`, pytest for helper tests, existing chr1 Parquet cache at `/Users/mwiewior/workspace/data_vepyr/115_GRCh38_merged`.

---

### Task 1: Test Benchmark Helper Contracts

**Files:**
- Create: `research/test_bench_lance_variation_chr1.py`
- Create later: `research/bench_lance_variation_chr1.py`

- [ ] **Step 1: Write failing tests**

Create `research/test_bench_lance_variation_chr1.py` with tests for the public helper functions the benchmark script will expose:

```python
from pathlib import Path

import pyarrow as pa

from research.bench_lance_variation_chr1 import (
    AF_GNOMADG_COLUMNS,
    ALWAYS_NULL_COLUMNS,
    EVERYTHING_VARIATION_COLUMNS,
    physical_projection_for_variant,
    result_markdown_table,
    schema_for_variant,
)


def test_everything_columns_preserve_known_null_logical_fields():
    assert "minor_allele" in EVERYTHING_VARIATION_COLUMNS
    assert "minor_allele_freq" in EVERYTHING_VARIATION_COLUMNS
    assert "minor_allele" in ALWAYS_NULL_COLUMNS
    assert "minor_allele_freq" in ALWAYS_NULL_COLUMNS


def test_21_projection_keeps_logical_top_level_names():
    logical = ["position_key", "allele_string", "gnomADg", "gnomADg_NFE"]
    assert physical_projection_for_variant("2.1-unpacked", logical) == logical


def test_22_packed_projection_maps_group_children_to_struct_once():
    logical = ["position_key", "gnomADg", "gnomADg_NFE", "AF", "AFR"]
    assert physical_projection_for_variant("2.2-packed", logical) == [
        "position_key",
        "af_gnomadg",
        "af_1kg",
    ]


def test_schema_for_22_packed_contains_packed_af_struct():
    source = pa.schema(
        [
            pa.field("position_key", pa.uint64(), nullable=False),
            pa.field("gnomADg", pa.string()),
            pa.field("gnomADg_NFE", pa.string()),
            pa.field("AF", pa.string()),
        ]
    )
    schema = schema_for_variant("2.2-packed", source, ["position_key", "gnomADg", "gnomADg_NFE", "AF"])
    assert "position_key" in schema.names
    assert "af_gnomadg" in schema.names
    assert schema.field("af_gnomadg").type == pa.struct(
        [pa.field(name, pa.string()) for name in ["gnomADg", "gnomADg_NFE"]]
    )
    assert schema.field("af_gnomadg").metadata[b"lance-encoding:packed"] == b"true"
    assert set(AF_GNOMADG_COLUMNS) >= {"gnomADg", "gnomADg_NFE"}


def test_result_markdown_table_formats_core_metrics():
    rows = [
        {
            "variant": "2.1-unpacked",
            "tier": "warm",
            "operation": "full_scan",
            "rows": 10,
            "seconds": 2.0,
            "rows_per_s": 5.0,
            "artifact_gib": 1.25,
        }
    ]
    markdown = result_markdown_table(rows)
    assert "| variant | tier | operation | rows | seconds | rows/s | artifact GiB |" in markdown
    assert "| 2.1-unpacked | warm | full_scan | 10 | 2.000 | 5 | 1.250 |" in markdown
```

- [ ] **Step 2: Run tests to verify they fail**

Run:

```bash
uv run --with pytest --with pyarrow pytest research/test_bench_lance_variation_chr1.py -q
```

Expected: FAIL with `ModuleNotFoundError: No module named 'research.bench_lance_variation_chr1'`.

### Task 2: Implement Benchmark Script

**Files:**
- Create: `research/bench_lance_variation_chr1.py`
- Test: `research/test_bench_lance_variation_chr1.py`

- [ ] **Step 1: Add minimal implementation for tested helpers**

Implement:

- `EVERYTHING_VARIATION_COLUMNS`
- `ALWAYS_NULL_COLUMNS`
- AF group constants
- `schema_for_variant()`
- `physical_projection_for_variant()`
- `result_markdown_table()`

The helper implementation must pass Task 1 before adding the heavier benchmark code.

- [ ] **Step 2: Run helper tests**

Run:

```bash
uv run --with pytest --with pyarrow pytest research/test_bench_lance_variation_chr1.py -q
```

Expected: PASS.

- [ ] **Step 3: Add Lance materialization and benchmark commands**

Extend `research/bench_lance_variation_chr1.py` with:

- `build_lance_dataset()` for warm/cold per variant.
- `benchmark_warm_full_scan()` using Parquet and Lance.
- `sample_cold_keys()` for deterministic 2,000-key cold samples.
- `benchmark_cold_point_lookup()` using `position_key IN (...)`.
- `write_reports()` for JSON and Markdown output.
- CLI flags:
  - `--cache-dir`
  - `--output-dir`
  - `--report-dir`
  - `--chrom`
  - `--variants`
  - `--cold-fragment-rows`
  - `--cold-sample-size`
  - `--force-build`
  - `--skip-build`

- [ ] **Step 4: Run smoke benchmark on limited rows**

Run:

```bash
uv run --with pylance --with pyarrow research/bench_lance_variation_chr1.py \
  --cache-dir /Users/mwiewior/workspace/data_vepyr/115_GRCh38_merged \
  --chrom chr1 \
  --variants 2.1-unpacked,2.2-packed \
  --output-dir /tmp/vepyr-lance-variation-smoke \
  --report-dir /tmp/vepyr-lance-variation-smoke-report \
  --cold-sample-size 128 \
  --row-limit 10000 \
  --force-build
```

Expected: exits 0 and writes:

```text
/tmp/vepyr-lance-variation-smoke-report/chr1_lance_variation_benchmark.json
/tmp/vepyr-lance-variation-smoke-report/chr1_lance_variation_benchmark.md
```

### Task 3: Run chr1 Benchmarks and Commit

**Files:**
- Modify: `research/bench_lance_variation_chr1.py`
- Create: `research/reports/chr1_lance_variation_benchmark.md`
- Create: `research/reports/chr1_lance_variation_benchmark.json`

- [ ] **Step 1: Run full chr1 benchmark**

Run:

```bash
uv run --with pylance --with pyarrow research/bench_lance_variation_chr1.py \
  --cache-dir /Users/mwiewior/workspace/data_vepyr/115_GRCh38_merged \
  --chrom chr1 \
  --variants 2.1-unpacked,2.2-packed \
  --output-dir /Users/mwiewior/workspace/data_vepyr/lance_variation_chr1_bench \
  --report-dir research/reports \
  --cold-fragment-rows 8192 \
  --cold-sample-size 2000 \
  --force-build
```

Expected: exits 0 and writes the chr1 report files under `research/reports/`.

- [ ] **Step 2: Run final tests**

Run:

```bash
uv run --with pytest --with pyarrow pytest research/test_bench_lance_variation_chr1.py -q
```

Expected: PASS.

- [ ] **Step 3: Commit implementation and report**

Run:

```bash
git add \
  docs/superpowers/plans/2026-06-09-lance-variation-chr1-benchmark.md \
  research/bench_lance_variation_chr1.py \
  research/test_bench_lance_variation_chr1.py \
  research/reports/chr1_lance_variation_benchmark.md \
  research/reports/chr1_lance_variation_benchmark.json
git commit -m "bench: add lance chr1 variation benchmark"
```

Expected: commit succeeds after hooks.

## Self-Review

Spec coverage:

- Lance 2.1 and 2.2 are both represented by benchmark variants.
- Warm full-scan and cold 2,000-key lookup are both measured.
- The logical `everything` projection preserves all-null fields.
- 2.2 packed structs are benchmarked as an enhanced variant.
- Results are written to a summary report.

Known scope boundary:

- This plan does not replace the Rust annotation runtime with Lance. It produces measured chr1 evidence first so the runtime migration can be planned against actual Lance performance.
