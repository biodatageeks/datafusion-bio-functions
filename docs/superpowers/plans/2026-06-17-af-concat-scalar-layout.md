# AF concat-scalar layout — Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Replace the fullzip `List<Utf8>` AF bundle with 3 concatenated scalar `Utf8` columns
(miniblock+zstd, `|`-joined) so variant point-lookups read ~3.6× fewer bytes and AF storage shrinks 3.8×.

**Architecture:** Each AF group (`af_global`/`af_gnomade`/`af_gnomadg`) becomes one scalar `Utf8`
column: members joined by `|`, empty field = absent member, null = all-absent. Encoded with the
shared `lance_field_metadata()` string preset (miniblock + zstd, 4 KB minichunk) — same as every
other string column, so no fullzip / no `List` levels / no large-chunk bug, on stock lance. Write
side concatenates; read side splits back into the 27 logical AF columns.

**Tech Stack:** Rust, Apache Arrow 58, DataFusion 53, Lance 7.0.0 (stock, crates.io).

**Spec:** `docs/superpowers/specs/2026-06-17-af-concat-scalar-layout-design.md`

---

## File structure

- `datafusion/bio-function-vep/src/lance_cache/af_bundle.rs` — all bundle/unbundle logic (write
  concat, read split, schema, projection, separator constant + assertion, tests). The only file
  with substantive changes.
- `datafusion/bio-function-vep/src/lance_cache/variation_runtime.rs` — remove the diagnostic
  row-id dump; the lance round-trip test lives here and must keep passing.
- Delete diagnostic examples: `examples/frag_rows.rs`, `examples/af_row_size.rs`,
  `examples/af_sep_check.rs`.

Current `af_bundle.rs` already contains a prototype `concat_field`/`concat_group` (env-gated) plus
the fullzip-List functions and a `VEP_AF_NOBUNDLE` baseline gate. This plan promotes concat to the
only path and deletes the rest.

---

### Task 1: Separator `|`, concat assertion, and `split_group` (pure logic + unit test)

**Files:**
- Modify: `datafusion/bio-function-vep/src/lance_cache/af_bundle.rs` (`AF_CONCAT_SEP` line 201,
  `concat_group` lines 222-249, `bundle_af_columns` concat branch)
- Test: same file, `mod tests`

- [ ] **Step 1: Write the failing round-trip test**

Add to `mod tests` in `af_bundle.rs`:

```rust
    #[test]
    fn concat_split_round_trips_values() {
        // 6-member group, 4 rows: full, partial, all-empty, all-null
        let members_owned = vec![
            sa(&[Some("A:0.1"), Some("x"), Some(""), None]),
            sa(&[Some("A:0.2"), None, Some(""), None]),
            sa(&[Some("A:0.3"), Some(""), Some(""), None]),
            sa(&[Some("A:0.4"), Some("y"), None, None]),
            sa(&[Some("A:0.5"), Some(""), Some(""), None]),
            sa(&[Some("A:0.6"), Some("z"), Some(""), None]),
        ];
        let members: Vec<&StringArray> = members_owned.iter().collect();
        let concat = concat_group(&members).unwrap();
        assert!(!concat.is_null(0));
        assert!(!concat.is_null(1));
        assert!(concat.is_null(2)); // all-empty -> null
        assert!(concat.is_null(3)); // all-null  -> null

        let back = split_group(&concat, 6);
        assert_eq!(back.len(), 6);
        let expect: Vec<Vec<&str>> = vec![
            vec!["A:0.1", "A:0.2", "A:0.3", "A:0.4", "A:0.5", "A:0.6"],
            vec!["x", "", "", "y", "", "z"],
            vec!["", "", "", "", "", ""],
            vec!["", "", "", "", "", ""],
        ];
        for row in 0..4 {
            for col in 0..6 {
                assert_eq!(back[col].value(row), expect[row][col], "row {row} col {col}");
            }
        }
    }

    #[test]
    fn concat_group_rejects_separator_in_value() {
        let members_owned = vec![sa(&[Some("a|b")]), sa(&[Some("c")])];
        let members: Vec<&StringArray> = members_owned.iter().collect();
        assert!(concat_group(&members).is_err());
    }
```

- [ ] **Step 2: Run to verify it fails**

Run: `cargo test -p datafusion-bio-function-vep --features lance-cache concat_split_round_trips_values concat_group_rejects_separator_in_value`
Expected: FAIL — `split_group` not found / `concat_group` returns `StringArray` not `Result`.

- [ ] **Step 3: Change the separator, make `concat_group` assert + return `Result`, add `split_group`**

Change line 201 from:

```rust
const AF_CONCAT_SEP: char = '\u{1f}';
```
to:
```rust
// `|` is verified collision-free in AF values (0 of 613M; `:`=allele/freq delimiter, `,`=8%).
// The build-time check in `concat_group` guards against future source changes.
const AF_CONCAT_SEP: char = '|';
```

Replace `concat_group` (lines 222-249) with:

```rust
/// Concatenate one group's members into a single Utf8 column (`|`-joined, positional,
/// null when all members absent). Errors if any member value contains the separator.
pub fn concat_group(members: &[&StringArray]) -> Result<StringArray> {
    let n = members[0].len();
    let width = members.len();
    let mut b = StringBuilder::with_capacity(n, n * 16);
    let mut s = String::new();
    for r in 0..n {
        s.clear();
        let mut all_absent = true;
        for (c, m) in members.iter().enumerate().take(width) {
            if c > 0 {
                s.push(AF_CONCAT_SEP);
            }
            if !m.is_null(r) {
                let v = m.value(r);
                if !v.is_empty() {
                    if v.contains(AF_CONCAT_SEP) {
                        return Err(DataFusionError::Execution(format!(
                            "AF value contains the '{AF_CONCAT_SEP}' bundle separator: {v:?}"
                        )));
                    }
                    all_absent = false;
                    s.push_str(v);
                }
            }
        }
        if all_absent {
            b.append_null();
        } else {
            b.append_value(&s);
        }
    }
    Ok(b.finish())
}

/// Split one bundled concat column back into its `width` member `Utf8` columns (positional;
/// absent member -> "", null parent -> "" for every member, matching downstream expectations).
pub fn split_group(col: &StringArray, width: usize) -> Vec<StringArray> {
    let n = col.len();
    let mut builders: Vec<StringBuilder> =
        (0..width).map(|_| StringBuilder::with_capacity(n, n * 4)).collect();
    for r in 0..n {
        if col.is_null(r) {
            for b in builders.iter_mut() {
                b.append_value("");
            }
        } else {
            let mut parts = col.value(r).split(AF_CONCAT_SEP);
            for b in builders.iter_mut() {
                b.append_value(parts.next().unwrap_or(""));
            }
        }
    }
    builders.into_iter().map(|mut b| b.finish()).collect()
}
```

Update the one caller in `bundle_af_columns` (the concat branch) from
`cols.push(Arc::new(concat_group(&arrays)) as ArrayRef);` to:

```rust
            cols.push(Arc::new(concat_group(&arrays)?) as ArrayRef);
```

- [ ] **Step 4: Run to verify it passes**

Run: `cargo test -p datafusion-bio-function-vep --features lance-cache concat_split_round_trips_values concat_group_rejects_separator_in_value`
Expected: PASS (2 tests).

- [ ] **Step 5: Commit**

```bash
git add datafusion/bio-function-vep/src/lance_cache/af_bundle.rs
git commit -m "feat(vep): concat AF separator '|' + split_group + collision assertion"
```

---

### Task 2: Make concat the only write path

**Files:**
- Modify: `af_bundle.rs` — `concat_field` (lines 211-220), `bundle_schema` (253-274),
  `bundle_af_columns` (289-318); add `lance_field_metadata` import.

- [ ] **Step 1: Use the shared string preset for concat fields, drop the minichunk env helper**

Add to the `use` of schema items at the top of `af_bundle.rs` (it already imports
`LANCE_STRUCTURAL_ENCODING_*`): add `lance_field_metadata`.

Replace `concat_field` (lines 211-220) with:

```rust
/// Field for a concatenated AF group: scalar `Utf8` encoded exactly like the other string
/// columns (miniblock + zstd + 4 KB minichunk via the shared preset) for full consistency.
fn concat_field(name: &str) -> Field {
    Field::new(name, DataType::Utf8, true).with_metadata(lance_field_metadata())
}
```

Delete `af_minichunk_meta` (lines 207-209) and `concat_enabled` (lines 203-205).

- [ ] **Step 2: Make `bundle_schema` always emit concat fields**

Replace `bundle_schema` (lines 253-274) with:

```rust
/// Schema with the 27 AF `Utf8` fields replaced by 3 concatenated scalar `Utf8` group fields
/// (appended after the non-AF fields, preserving non-AF order + schema metadata).
pub fn bundle_schema(schema: &Schema) -> Schema {
    let af: HashSet<&str> = af_column_order().into_iter().collect();
    let mut fields: Vec<Arc<Field>> = schema
        .fields()
        .iter()
        .filter(|f| !af.contains(f.name().as_str()))
        .cloned()
        .collect();
    for (name, _) in AF_GROUPS {
        fields.push(Arc::new(concat_field(name)));
    }
    Schema::new_with_metadata(fields, schema.metadata().clone())
}
```

- [ ] **Step 3: Make `bundle_af_columns` always concatenate**

Replace `bundle_af_columns` (lines 289-318) with:

```rust
pub fn bundle_af_columns(batch: &RecordBatch) -> Result<RecordBatch> {
    let af: HashSet<&str> = af_column_order().into_iter().collect();
    let schema = batch.schema();
    let mut cols: Vec<ArrayRef> = Vec::new();
    for (i, f) in schema.fields().iter().enumerate() {
        if !af.contains(f.name().as_str()) {
            cols.push(batch.column(i).clone());
        }
    }
    for (_, members) in AF_GROUPS {
        let arrays = members
            .iter()
            .map(|m| string_col(batch, m))
            .collect::<Result<Vec<_>>>()?;
        cols.push(Arc::new(concat_group(&arrays)?) as ArrayRef);
    }
    let out_schema = Arc::new(bundle_schema(schema.as_ref()));
    RecordBatch::try_new(out_schema, cols)
        .map_err(|e| DataFusionError::ArrowError(Box::new(e), None))
}
```

- [ ] **Step 4: Verify it compiles (fullzip helpers now unused — that's expected, removed in Task 4)**

Run: `cargo build -p datafusion-bio-function-vep --features lance-cache,cache-builder`
Expected: builds (warnings about unused `bundle_group`/`bundled_list_field`/etc. are OK).

- [ ] **Step 5: Commit**

```bash
git add datafusion/bio-function-vep/src/lance_cache/af_bundle.rs
git commit -m "feat(vep): make concat the only AF write path (4 KB shared preset)"
```

---

### Task 3: Rewrite the read path to split

**Files:**
- Modify: `af_bundle.rs` — `unbundle_af_columns` (lines 320-347)

- [ ] **Step 1: Write the failing read-path test**

Add to `mod tests`:

```rust
    #[test]
    fn unbundle_af_columns_splits_concat() {
        use datafusion::arrow::datatypes::{Field, Schema};
        // build a batch with one non-AF col + the 3 concat group cols
        let chrom = sa(&[Some("chr1"), Some("chr1")]);
        let af_global = StringArray::from(vec![Some("A:0.1|||||"), None]);
        let af_gnomade = StringArray::from(vec![Some("g|||||||||"), None]);
        let af_gnomadg = StringArray::from(vec![Some("h||||||||||"), None]);
        let schema = Arc::new(Schema::new(vec![
            Field::new("chrom", DataType::Utf8, true),
            Field::new("af_global", DataType::Utf8, true),
            Field::new("af_gnomade", DataType::Utf8, true),
            Field::new("af_gnomadg", DataType::Utf8, true),
        ]));
        let batch = RecordBatch::try_new(
            schema,
            vec![
                Arc::new(chrom),
                Arc::new(af_global),
                Arc::new(af_gnomade),
                Arc::new(af_gnomadg),
            ],
        )
        .unwrap();
        let out = unbundle_af_columns(&batch).unwrap();
        // 1 non-AF + 27 AF = 28 columns
        assert_eq!(out.num_columns(), 28);
        let af_idx = out.schema().index_of("AF").unwrap();
        let af = out
            .column(af_idx)
            .as_any()
            .downcast_ref::<StringArray>()
            .unwrap();
        assert_eq!(af.value(0), "A:0.1");
        assert_eq!(af.value(1), ""); // null parent -> ""
    }
```

- [ ] **Step 2: Run to verify it fails**

Run: `cargo test -p datafusion-bio-function-vep --features lance-cache unbundle_af_columns_splits_concat`
Expected: FAIL — `unbundle_af_columns` still expects a `ListArray` (downcast error).

- [ ] **Step 3: Rewrite `unbundle_af_columns` to split scalar Utf8**

Replace `unbundle_af_columns` (lines 320-347) with:

```rust
/// Read side: expand each bundled concat `Utf8` group column back into its member `Utf8`
/// columns. Pass-through for non-group columns; robust to 0..3 groups present.
pub fn unbundle_af_columns(batch: &RecordBatch) -> Result<RecordBatch> {
    let schema = batch.schema();
    let mut fields: Vec<Arc<Field>> = Vec::new();
    let mut cols: Vec<ArrayRef> = Vec::new();
    for (i, f) in schema.fields().iter().enumerate() {
        if let Some((_, members)) = AF_GROUPS.iter().find(|(n, _)| *n == f.name()) {
            let concat = batch
                .column(i)
                .as_any()
                .downcast_ref::<StringArray>()
                .ok_or_else(|| {
                    DataFusionError::Execution(format!(
                        "bundled column '{}' must be Utf8",
                        f.name()
                    ))
                })?;
            for (m, arr) in members.iter().zip(split_group(concat, members.len())) {
                fields.push(Arc::new(Field::new(*m, DataType::Utf8, true)));
                cols.push(Arc::new(arr) as ArrayRef);
            }
        } else {
            fields.push(f.clone());
            cols.push(batch.column(i).clone());
        }
    }
    RecordBatch::try_new(
        Arc::new(Schema::new_with_metadata(fields, schema.metadata().clone())),
        cols,
    )
    .map_err(|e| DataFusionError::ArrowError(Box::new(e), None))
}
```

- [ ] **Step 4: Run to verify it passes**

Run: `cargo test -p datafusion-bio-function-vep --features lance-cache unbundle_af_columns_splits_concat`
Expected: PASS.

- [ ] **Step 5: Commit**

```bash
git add datafusion/bio-function-vep/src/lance_cache/af_bundle.rs
git commit -m "feat(vep): split concat AF columns back to 27 on read"
```

---

### Task 4: Remove the dead fullzip-List machinery and fix the lance round-trip test

**Files:**
- Modify: `af_bundle.rs` — delete `bundle_group` (84-131), `unbundle_group` (133-161),
  `list_item_field` (163-179), `bundled_list_field` (181-194); delete the old
  `bundle_unbundle_round_trips_values` test (it uses `bundle_group`/`unbundle_group`).
- Modify: `variation_runtime.rs` — `lance_lookup_round_trips_bundled_af_columns` test (must pass
  unchanged with concat, since it goes through `bundle_af_columns` + `take_rows` + `unbundle`).

- [ ] **Step 1: Delete the four fullzip functions and the old List round-trip test**

Remove `bundle_group`, `unbundle_group`, `list_item_field`, `bundled_list_field` from
`af_bundle.rs`, and delete the `bundle_unbundle_round_trips_values` test. Remove now-unused
imports (`ListArray`, `OffsetBuffer`, `LANCE_STRUCTURAL_ENCODING_FULLZIP`,
`LANCE_STRUCTURAL_ENCODING_KEY` if no longer referenced).

- [ ] **Step 2: Run the full lance_cache test suite**

Run: `cargo test -p datafusion-bio-function-vep --features lance-cache,cache-builder --lib lance_cache::`
Expected: PASS, including `lance_cache::variation_runtime::tests::lance_lookup_round_trips_bundled_af_columns` (it now round-trips through concat).

- [ ] **Step 3: If `lance_lookup_round_trips_bundled_af_columns` fails**, it asserts the 27 AF
  columns come back with exact values after `take_rows`+`unbundle`. With concat, absent members
  return `""` (not null) — if the test asserted nulls, update its expectations to `""` for absent
  members (matching `split_group`). Show the diff in the commit.

- [ ] **Step 4: Clippy + fmt**

Run: `cargo clippy -p datafusion-bio-function-vep --features lance-cache,cache-builder -- -D warnings && cargo fmt`
Expected: clean.

- [ ] **Step 5: Commit**

```bash
git add datafusion/bio-function-vep/src/lance_cache/af_bundle.rs datafusion/bio-function-vep/src/lance_cache/variation_runtime.rs
git commit -m "refactor(vep): drop fullzip-List AF bundle, concat is the only layout"
```

---

### Task 5: Remove diagnostic scaffolding

**Files:**
- Modify: `variation_runtime.rs` — delete the `VEP_DUMP_ROWIDS` block (the `if let Some(path) =
  std::env::var_os("VEP_DUMP_ROWIDS")` dump added for benchmarking).
- Delete: `datafusion/bio-function-vep/examples/frag_rows.rs`,
  `datafusion/bio-function-vep/examples/af_row_size.rs`,
  `datafusion/bio-function-vep/examples/af_sep_check.rs`.

- [ ] **Step 1: Remove the row-id dump from `resolve_and_take`**

In `variation_runtime.rs`, delete the whole `if let Some(path) = std::env::var_os("VEP_DUMP_ROWIDS") { ... }` block (added in this session for capturing benchmark ids).

- [ ] **Step 2: Delete the diagnostic example files**

```bash
git rm datafusion/bio-function-vep/examples/frag_rows.rs \
       datafusion/bio-function-vep/examples/af_row_size.rs \
       datafusion/bio-function-vep/examples/af_sep_check.rs
```

- [ ] **Step 3: Verify build + tests still green**

Run: `cargo test -p datafusion-bio-function-vep --features lance-cache,cache-builder --lib lance_cache::`
Expected: PASS.

- [ ] **Step 4: Commit**

```bash
git add datafusion/bio-function-vep/src/lance_cache/variation_runtime.rs
git commit -m "chore(vep): remove AF-bundle benchmark scaffolding"
```

---

### Task 6: Rebuild cache, rebuild vepyr, run e2e parity (acceptance gate)

**Files:** none (data + validation). This is the gate that proves split-on-read is correct.

- [ ] **Step 1: Build the chr1 variation cache with concat**

```bash
cd /Users/mwiewior/research/git/datafusion-bio-functions
cargo build --release --example build_lance_variation_chrom --features lance-cache,cache-builder
./target/release/examples/build_lance_variation_chrom \
  --cache-root /Users/mwiewior/workspace/data_vepyr/homo_sapiens_merged/115_GRCh38 \
  --output-dir /Users/mwiewior/workspace/data_vepyr/115_GRCh38_merged \
  --chrom chr1 --cache-source-type merged --partitions 8 --overwrite
```
Expected: `wrote chrom=chr1 ... rows=88153966`. (Context datasets untouched.)

- [ ] **Step 2: Confirm AF columns encoded as expected**

```bash
SB=research/lance_encoding_sandbox/crates/lance_sandbox/target/release/lance_sandbox
$SB inspect-path --dataset-path /Users/mwiewior/workspace/data_vepyr/115_GRCh38_merged/variation.lance/chr1.lance --position-field start --lance-version 2.2 2>/dev/null | grep -A2 'af_'
```
Expected: `af_global`/`af_gnomade`/`af_gnomadg` are scalar `Utf8`, `MiniBlock` + `zstd(level=3)`,
AF total ≈ 1.5 GB.

- [ ] **Step 3: Rebuild the vepyr extension (picks up split-on-read via the path dep)**

```bash
cd /Users/mwiewior/research/git/vepyr
env -u VIRTUAL_ENV -u CONDA_PREFIX uv run maturin develop --release
```
Expected: `🛠 Installed vepyr-0.1.0`.

- [ ] **Step 4: Run the e2e WITH comparison (parity gate)**

```bash
cd /Users/mwiewior/research/git/vepyr/e2e-testing/scripts
LANCE_CPU_THREADS=1 LANCE_IO_THREADS=1 RAYON_NUM_THREADS=1 \
  env -u VIRTUAL_ENV -u CONDA_PREFIX uv run python run_annotation_fast.py \
  chr1 --cache merged --forks 0 --force --backend lance
```
Expected: completes (323,430 variants), no panic, and the comparison step reports
**0 CSQ/AF mismatches** vs the VEP reference. The AF concordance is the proof that split-on-read
reproduces the 27 columns exactly.

- [ ] **Step 5: If mismatches are AF-only**, inspect a mismatched variant's AF field: the likely
  cause is the absent-member representation (`""` vs null) differing from what the annotation
  formats. Compare `split_group` output for that variant against the source AF columns; adjust the
  absent-member normalization in `split_group` (e.g. `append_null()` instead of `append_value("")`)
  and re-run from Step 1. Do **not** proceed until mismatches are 0.

- [ ] **Step 6: Commit a short result note**

```bash
git add docs/superpowers/specs/2026-06-17-af-concat-scalar-layout-design.md
git commit -m "docs(vep): AF concat layout validated e2e (0 mismatches, AF ~1.5 GB)"
```

---

## Out of scope (follow-ups, not this plan)
- Non-AF minichunk tuning (the non-AF columns now dominate reads — likely the next win).
- Row-major / cold-tier point store (the only thing that breaks the columnar ~15× point-lookup floor).
- Genome-wide variation rebuild (this plan is chr1 only).
