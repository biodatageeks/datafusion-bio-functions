# Lance 2.2 Struct Packing Benchmark Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Add a targeted Lance 2.2 packed-struct benchmark profile and make the Rust benchmark plus PyLance notebook handle both flat and packed schemas for the same VEPyr `--everything` logical workload.

**Architecture:** Keep lookup keys, tier, and other index/filter fields top-level. Pack only logical fields that are always read together in the `--everything` path into Lance 2.2 struct groups. The benchmark resolves logical fields to top-level fields or packed parent fields and reports both requested logical fields and physical projections.

**Tech Stack:** Rust sandbox CLI, Lance 7.0.0, Arrow, TOML configs, Jupyter/PyLance notebook, Python `rich`.

---

### Task 1: Projection Compatibility Tests

**Files:**
- Modify: `research/lance_encoding_sandbox/crates/lance_sandbox/src/build.rs`

- [ ] **Step 1: Add failing tests**

Add tests in the existing `#[cfg(test)] mod tests` in `build.rs` that construct small Lance schemas and verify:

```rust
#[test]
fn everything_projection_keeps_flat_fields_top_level() {
    let config = test_config("2.1", "");
    let schema = lance_schema(&[
        ("position", DataType::UInt32),
        ("tier", DataType::UInt8),
        ("variation_name", DataType::Utf8),
        ("dbsnp_ids", DataType::Utf8),
    ]);

    let projection = physical_projection_for_everything(&config, &schema).unwrap();

    assert!(projection.contains(&"position".to_string()));
    assert!(projection.contains(&"variation_name".to_string()));
    assert!(projection.contains(&"dbsnp_ids".to_string()));
}

#[test]
fn everything_projection_deduplicates_packed_parent_fields() {
    let config = test_config(
        "2.2",
        r#"
[structs.identity_text]
enabled = true
packed_metadata = true
fields = ["variation_name", "dbsnp_ids"]
"#,
    );
    let schema = lance_schema(&[
        ("position", DataType::UInt32),
        ("tier", DataType::UInt8),
        ("identity_text", DataType::Struct(Fields::from(vec![
            Field::new("variation_name", DataType::Utf8, true),
            Field::new("dbsnp_ids", DataType::Utf8, true),
        ]))),
    ]);

    let projection = physical_projection_for_everything(&config, &schema).unwrap();

    assert!(projection.contains(&"position".to_string()));
    assert_eq!(
        projection.iter().filter(|name| name.as_str() == "identity_text").count(),
        1
    );
    assert!(!projection.contains(&"variation_name".to_string()));
    assert!(!projection.contains(&"dbsnp_ids".to_string()));
}
```

- [ ] **Step 2: Run red test**

Run:

```bash
cargo test --manifest-path research/lance_encoding_sandbox/crates/lance_sandbox/Cargo.toml build::tests::everything_projection_deduplicates_packed_parent_fields
```

Expected before implementation: the packed-schema test fails if logical-to-physical schema handling is incomplete.

- [ ] **Step 3: Implement resolver cleanup if needed**

Update `physical_projection_for_everything` only if the tests reveal missing behavior. The resolver must preserve top-level lookup fields and emit each packed parent at most once.

- [ ] **Step 4: Run green tests**

Run:

```bash
cargo test --manifest-path research/lance_encoding_sandbox/crates/lance_sandbox/Cargo.toml build::tests::everything_projection
```

Expected: both projection tests pass.

### Task 2: Targeted Lance 2.2 Benchmark Config

**Files:**
- Create: `research/lance_encoding_sandbox/configs/packed_targeted_v22.toml`
- Modify: `research/lance_encoding_sandbox/README.md`

- [ ] **Step 1: Add config**

Create a Lance 2.2 config with `dataset.name = "packed_targeted_v22_position_u32"`, `lance_version = "2.2"`, the same index and benchmark settings as the current 2.1 profile, and these struct groups:

```toml
[structs.match_payload]
enabled = true
packed_metadata = true
fields = ["allele_string", "end", "failed"]

[structs.identity_text]
enabled = true
packed_metadata = true
fields = ["variation_name", "dbsnp_ids"]

[structs.clinical_payload]
enabled = true
packed_metadata = true
fields = ["clin_sig", "clin_sig_allele", "clinical_impact", "pubmed", "clinvar_ids", "cosmic_ids"]

[structs.variant_flags]
enabled = true
packed_metadata = true
fields = ["somatic", "phenotype_or_disease", "strand"]

[structs.freq_1kg]
enabled = true
packed_metadata = true
fields = ["AF", "AFR", "AMR", "EAS", "EUR", "SAS"]

[structs.freq_gnomade]
enabled = true
packed_metadata = true
fields = ["gnomADe", "gnomADe_AFR", "gnomADe_AMR", "gnomADe_ASJ", "gnomADe_EAS", "gnomADe_FIN", "gnomADe_NFE", "gnomADe_REMAINING", "gnomADe_SAS", "gnomADe_MID"]

[structs.freq_gnomadg_core]
enabled = true
packed_metadata = true
fields = ["gnomADg", "gnomADg_AFR", "gnomADg_NFE", "gnomADg_AMR"]

[structs.freq_gnomadg_tail]
enabled = true
packed_metadata = true
fields = ["gnomADg_AMI", "gnomADg_ASJ", "gnomADg_EAS", "gnomADg_FIN", "gnomADg_MID", "gnomADg_SAS", "gnomADg_REMAINING"]
```

- [ ] **Step 2: Document command**

Add a README example showing:

```bash
RUSTFLAGS="-C target-cpu=native" cargo run --release --manifest-path research/lance_encoding_sandbox/crates/lance_sandbox/Cargo.toml -- \
  run --config research/lance_encoding_sandbox/configs/packed_targeted_v22.toml \
  --positions-file research/lance_encoding_sandbox/inputs/chr1_cold_sample_10k_positions_u32.txt
```

- [ ] **Step 3: Parse config**

Run:

```bash
cargo test --manifest-path research/lance_encoding_sandbox/crates/lance_sandbox/Cargo.toml config::tests::checked_in_configs_parse
```

Expected: the new config parses and validates as Lance 2.2.

### Task 3: Notebook Struct Packing Support

**Files:**
- Modify: `research/lance_encoding_sandbox/notebooks/pylance_variation_builder.ipynb`

- [ ] **Step 1: Add notebook config constants**

Add:

```python
STRUCT_PACKING_ENABLED = False
STRUCT_PACKING_GROUPS = {
    'match_payload': ['allele_string', 'end', 'failed'],
    'identity_text': ['variation_name', 'dbsnp_ids'],
    'clinical_payload': ['clin_sig', 'clin_sig_allele', 'clinical_impact', 'pubmed', 'clinvar_ids', 'cosmic_ids'],
    'variant_flags': ['somatic', 'phenotype_or_disease', 'strand'],
    'freq_1kg': ['AF', 'AFR', 'AMR', 'EAS', 'EUR', 'SAS'],
    'freq_gnomade': ['gnomADe', 'gnomADe_AFR', 'gnomADe_AMR', 'gnomADe_ASJ', 'gnomADe_EAS', 'gnomADe_FIN', 'gnomADe_NFE', 'gnomADe_REMAINING', 'gnomADe_SAS', 'gnomADe_MID'],
    'freq_gnomadg_core': ['gnomADg', 'gnomADg_AFR', 'gnomADg_NFE', 'gnomADg_AMR'],
    'freq_gnomadg_tail': ['gnomADg_AMI', 'gnomADg_ASJ', 'gnomADg_EAS', 'gnomADg_FIN', 'gnomADg_MID', 'gnomADg_SAS', 'gnomADg_REMAINING'],
}
```

- [ ] **Step 2: Pack schema and batches**

Update `build_target_schema` to replace configured child fields with `pa.struct` parents when `STRUCT_PACKING_ENABLED` is true. Update `transform_batch` to emit `pa.StructArray.from_arrays(...)` for packed parent fields and keep ungrouped fields unchanged.

- [ ] **Step 3: Inspect/benchmark generated dataset**

The existing Rust `inspect-path` and `bench --dataset-path` cells should continue using `OUTPUT_DATASET`, so they must work for flat and packed notebook datasets.

- [ ] **Step 4: Validate notebook syntax**

Run a Python JSON/AST check over every code cell:

```bash
python - <<'PY'
import ast, json
from pathlib import Path
nb = json.loads(Path('research/lance_encoding_sandbox/notebooks/pylance_variation_builder.ipynb').read_text())
for idx, cell in enumerate(nb['cells']):
    if cell.get('cell_type') == 'code':
        ast.parse(''.join(cell.get('source', [])), filename=f'cell-{idx}')
print('ok')
PY
```

Expected: `ok`.

### Task 4: Full Verification

**Files:**
- All modified sandbox files.

- [ ] **Step 1: Format Rust**

Run:

```bash
cargo fmt --manifest-path research/lance_encoding_sandbox/crates/lance_sandbox/Cargo.toml --check
```

- [ ] **Step 2: Run Rust tests**

Run:

```bash
cargo test --manifest-path research/lance_encoding_sandbox/crates/lance_sandbox/Cargo.toml
```

- [ ] **Step 3: Release build**

Run:

```bash
RUSTFLAGS="-C target-cpu=native" cargo build --release --manifest-path research/lance_encoding_sandbox/crates/lance_sandbox/Cargo.toml
```

- [ ] **Step 4: Report outcome**

Summarize:

- Added config path.
- Notebook toggle name.
- Whether the Rust benchmark now resolves physical projections for flat and packed schemas.
- Verification commands and results.
