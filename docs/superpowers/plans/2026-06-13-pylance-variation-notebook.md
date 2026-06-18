# PyLance Variation Notebook Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Build a Python-only notebook that creates a Lance variation dataset from chr Parquet files with configurable per-field Lance encoding metadata, scalar indexes, and Rich table inspection.

**Architecture:** The notebook owns the PyArrow/PyLance workflow and writes datasets under the existing Lance encoding sandbox output root. It streams warm and cold Parquet batches, adds `tier` and chr-local `position UInt32`, applies Arrow field metadata, creates Lance scalar indexes, and inspects exposed logical/index/physical metadata with Python tables. Rust is not used by the notebook.

**Tech Stack:** `uv`, Python, Jupyter, PyArrow, pandas, pylance/import `lance`, Rich.

---

### Task 1: Create Python Environment

**Files:**
- Create: `research/lance_encoding_sandbox/requirements-pylance.txt`
- Create: `research/lance_encoding_sandbox/.venv-pylance/`

- [ ] **Step 1: Add dependency file**

Create `research/lance_encoding_sandbox/requirements-pylance.txt` with:

```text
pylance
pyarrow
pandas
rich
jupyter
ipykernel
nbformat
tabulate
```

- [ ] **Step 2: Create and populate uv venv**

Run:

```bash
uv venv research/lance_encoding_sandbox/.venv-pylance
uv pip install --python research/lance_encoding_sandbox/.venv-pylance/bin/python \
  -r research/lance_encoding_sandbox/requirements-pylance.txt
```

Expected: `python -c "import lance, pyarrow, pandas, rich"` succeeds in the venv.

### Task 2: Create Notebook

**Files:**
- Create: `research/lance_encoding_sandbox/notebooks/pylance_variation_builder.ipynb`

- [ ] **Step 1: Add config cells**

Notebook config defines source Parquet paths, output Lance path, file version, overwrite mode, batch size, default metadata, and `FIELD_METADATA`.

- [ ] **Step 2: Add schema and batch helpers**

Helpers apply Arrow field metadata, add `tier`, derive `position UInt32` from `start`, and stream warm/cold Parquet batches into `lance.write_dataset`.

- [ ] **Step 3: Add index cells**

Notebook opens the dataset and creates a BTree index on `position` and bitmap index on `tier` using the PyLance APIs exposed by the installed package.

- [ ] **Step 4: Add inspection cells**

Notebook renders Rich tables for schema metadata, indexes, file sizes, and physical fragments/column metadata where PyLance exposes it. It writes `reports/inspect_python.json` and Markdown tables. Page-level fields that PyLance does not expose are marked `not exposed by pylance`.

### Task 3: Verify

**Files:**
- Verify: `research/lance_encoding_sandbox/notebooks/pylance_variation_builder.ipynb`
- Verify: `research/lance_encoding_sandbox/requirements-pylance.txt`

- [ ] **Step 1: Validate notebook JSON**

Run:

```bash
python -m json.tool research/lance_encoding_sandbox/notebooks/pylance_variation_builder.ipynb >/dev/null
```

Expected: exit code 0.

- [ ] **Step 2: Validate dependency imports**

Run:

```bash
research/lance_encoding_sandbox/.venv-pylance/bin/python - <<'PY'
import lance, pyarrow, pandas, rich
print("ok")
PY
```

Expected: prints `ok`.

- [ ] **Step 3: Validate core notebook helper syntax**

Run the notebook helper cells with the venv Python without building the full chr1 dataset.

Expected: helper functions import and compile, and API discovery prints supported Lance index methods.
