# Engine Forks/Threads Cleanup → Single `workers` Knob (coordinated with vepyr)

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Remove the dead cross-contig `forks`/`chromosome_lanes` parallelism machinery and the serial `inline_lookup` path from the VEP engine, rename the surviving within-contig knob `threads` → `workers`, and collapse vepyr to a single `--workers` knob — all in one coordinated change verified byte-identical on chr1.

**Architecture:** vepyr is the **sole consumer** of `datafusion-bio-function-vep` and depends on it via a **local `path` dependency** (`vepyr/Cargo.toml:34`), so there is no git-rev re-pin — `maturin develop` rebuilds vepyr against the local engine. The engine's `threads` path and `forks` path are **mutually exclusive by validation** and share no execution logic, so the forks path is a clean excision. The serial case has no `if partitions==1` gate; the only thing forcing the inline path is the `inline_lookup` flag, so removing inline makes `workers==1` route through one spawned lookup worker (`target_partitions==1`).

**Tech Stack:** Rust (datafusion, tokio, pyo3), Python (vepyr PyO3 extension), maturin, cargo, pytest.

## Global Constraints

- **Single source of truth for the knob name:** after this change the only annotation-concurrency knob is `workers` — engine `AnnotateVcfConfig.workers` (renamed from `threads`), options-JSON key `"workers"` (renamed from `"threads"`), internal `annotation_workers` (renamed from `annotation_threads`), Python/CLI `workers`. No `forks`, no `threads`, no per-contig `workers`, no `chromosome_lanes`, no `inline_lookup`, no `contig_parallelism` anywhere.
- **`workers > 1` requires a tabix-indexed input** (`.tbi`/`.csi`) — enforced by `annotate_to_vcf`. `workers == 1` is the serial path (now via one spawned lookup worker).
- **`target_partitions` field stays** on `AnnotateVcfConfig` (cold-Parquet/Lance shard count); the provider still bumps it to `workers.max(target_partitions)`.
- **Behavior-preserving refactor.** The regression gate is a **byte-identical chr1 comparison** against a baseline captured BEFORE any edits (Task 0). Serial output must be unchanged; `workers>1` output must remain byte-identical to serial (the engine's existing property).
- **Accepted tradeoff (explicitly chosen):** removing the inline path means the serial (`workers==1`) case now spawns one tokio lookup task + uses an mpsc channel instead of inline polling — a small per-run overhead.
- After all tasks: in `datafusion-bio-functions` — `cargo test`, `cargo clippy -- -D warnings`, `cargo fmt -- --check` clean. In `vepyr` — `cargo test`, `cargo clippy -- -D warnings`, `cargo fmt -- --check`, `pytest tests/` green, and the chr1 e2e gate byte-identical.
- **`partitions`** (cache-build knob in `build_cache`) is out of scope — do not touch.
- All engine line numbers below are relative to current `HEAD` (`2fe78f6`). Verify the surrounding verbatim code matches before editing; if a span has drifted, search by the quoted code, not the line number.

## File Structure

**Engine — `/Users/mwiewior/research/git/datafusion-bio-functions/datafusion/bio-function-vep/src/`:**
- `vcf_sink.rs` — `AnnotateVcfConfig` (struct/Default/Debug/`to_options_json`), `VepConcurrencyPlan`, `annotate_to_vcf` validation + output dispatch, the 4 forks fan-out functions + `PartitionedVcfBody`, and the forks/scheduler unit tests.
- `annotate_provider.rs` — options parse block, `ContigAnnotationConfig`, `ContigAnnotationExec::new`, `partition_contigs_for_execution`, the inline-lookup machinery (`InlineLookupPartitionHandle`, `LookupPartitionHandle::Inline`, `inline_lookup_partition`, the two selection sites, the poll arm, the sharded `inline_seen` guard), and the test constructor.

**vepyr — `/Users/mwiewior/research/git/vepyr/`:**
- `src/annotate.rs`, `src/lib.rs`, `src/vepyr/_core.pyi`, `src/vepyr/__init__.py` — single `workers` knob; build `AnnotateVcfConfig { workers, target_partitions, .. }` with no `forks`/old-`workers` fields; emit options-JSON key `"workers"`.
- `tests/test_annotate.py`, `tests/test_build_cache.py` — update to the single-knob contract.
- `e2e-testing/scripts/run_annotation_fast.py`, `run_annotation_fast_all.py` — single `--workers` flag.
- `README.md`, `docs/quickstart.md`, `docs/performance.md` — docs.

> **Note:** this plan supersedes and absorbs `vepyr/docs/superpowers/plans/2026-06-21-single-workers-knob.md` (the vepyr-only version). The vepyr edits here differ because the engine field is renamed to `workers` (not left as `threads`).

---

## Task 0: Baseline capture + branch

**Files:** none (setup).

- [ ] **Step 1: Create working branches in both repos**

```bash
cd /Users/mwiewior/research/git/datafusion-bio-functions && git checkout -b engine-single-workers-knob
cd /Users/mwiewior/research/git/vepyr && git checkout -b single-workers-knob
```

- [ ] **Step 2: Build the current engine into vepyr and capture a chr1 baseline**

Run (uses the existing e2e driver against the current/old API; produces the reference output to diff against at the end):

```bash
cd /Users/mwiewior/research/git/vepyr
maturin develop --release 2>&1 | tail -3
python e2e-testing/scripts/run_annotation_fast.py chr1 --backend lance --skip-compare --force 2>&1 | tail -15
```

- [ ] **Step 3: Save the baseline output path**

Locate the annotated chr1 VCF the driver just wrote (printed in its output, typically under the cache profile's results dir). Copy it to a stable location and record its checksum:

```bash
# Replace <CHR1_OUTPUT> with the path printed by the driver above.
cp <CHR1_OUTPUT> /tmp/chr1_baseline.vcf
grep -v '^##' /tmp/chr1_baseline.vcf | md5 | tee /tmp/chr1_baseline.md5
```

Expected: a non-empty md5 recorded. This is the byte-identical gate target for Task 8. (We strip `##` header lines, which may carry version/timestamp noise, and compare data + column header.)

---

## Task 1: `AnnotateVcfConfig` — drop forks/old-workers, rename threads→workers

**Files:**
- Modify: `datafusion/bio-function-vep/src/vcf_sink.rs` (struct ~855-875, Default ~877-916, Debug ~918-931, `to_options_json` ~1010-1037)

**Interfaces:**
- Produces: `AnnotateVcfConfig { ..., target_partitions: usize, workers: usize, ... }` — `workers` is the within-contig knob (was `threads`); `forks` and the old per-contig `workers` are gone. `to_options_json` emits key `"workers"` (was `"threads"`) and no longer emits `"forks"`/`"contig_parallelism"`/`"inline_lookup"`.

- [ ] **Step 1: Edit the struct fields (vcf_sink.rs ~856-866)**

Replace:

```rust
    /// User-facing chromosome lane count. Missing or `Some(0)` is a strict
    /// single-lane path with no helper lookup task; `Some(n)` runs up to `n`
    /// contigs concurrently.
    pub forks: Option<usize>,
    /// Number of annotation lookup workers per active contig.
    pub workers: usize,
    /// Maximum number of independent cold-Parquet lookup shards per contig.
    pub target_partitions: usize,
    /// Number of parallel fused window-annotation workers within a contig
    /// (`1` = serial). The single within-contig parallelism knob.
    pub threads: usize,
```

with:

```rust
    /// Maximum number of independent cold-Parquet lookup shards per contig.
    pub target_partitions: usize,
    /// Number of parallel fused window-annotation pipelines within a contig
    /// (`1` = serial). The single annotation-concurrency knob.
    pub workers: usize,
```

- [ ] **Step 2: Edit the `Default` impl (vcf_sink.rs ~907-910)**

Replace:

```rust
            forks: Some(0),
            workers: 1,
            target_partitions: 1,
            threads: 1,
```

with:

```rust
            target_partitions: 1,
            workers: 1,
```

- [ ] **Step 3: Edit the `Debug` impl (vcf_sink.rs ~923-925)**

Replace:

```rust
            .field("forks", &self.forks)
            .field("workers", &self.workers)
            .field("target_partitions", &self.target_partitions)
```

with:

```rust
            .field("workers", &self.workers)
            .field("target_partitions", &self.target_partitions)
```

- [ ] **Step 4: Edit `to_options_json` (vcf_sink.rs ~1014-1036)**

Replace:

```rust
        if let Some(forks) = self.forks {
            opts.insert(
                "forks".into(),
                serde_json::Value::Number(serde_json::Number::from(if forks == 0 {
                    0
                } else {
                    self.workers.max(1)
                })),
            );
            opts.insert(
                "contig_parallelism".into(),
                serde_json::Value::Number(serde_json::Number::from(forks.max(1))),
            );
            opts.insert("inline_lookup".into(), serde_json::Value::Bool(forks == 0));
        }
        opts.insert(
            "target_partitions".into(),
            serde_json::Value::Number(serde_json::Number::from(self.target_partitions.max(1))),
        );
        opts.insert(
            "threads".into(),
            serde_json::Value::Number(serde_json::Number::from(self.threads.max(1))),
        );
```

with:

```rust
        opts.insert(
            "target_partitions".into(),
            serde_json::Value::Number(serde_json::Number::from(self.target_partitions.max(1))),
        );
        opts.insert(
            "workers".into(),
            serde_json::Value::Number(serde_json::Number::from(self.workers.max(1))),
        );
```

- [ ] **Step 5: Verify the crate still references the old fields elsewhere (expected: compile errors guide the next tasks)**

Run: `cd /Users/mwiewior/research/git/datafusion-bio-functions && cargo build -p datafusion-bio-function-vep 2>&1 | grep -E "error\[|no field|cannot find" | head -30`
Expected: errors at the `VepConcurrencyPlan`, validation, and provider sites — these are fixed in Tasks 2-6. (Do not commit yet; the crate is intentionally mid-refactor.)

---

## Task 2: `VepConcurrencyPlan` + `annotate_to_vcf` validation/dispatch/profile

**Files:**
- Modify: `datafusion/bio-function-vep/src/vcf_sink.rs` (`VepConcurrencyPlan` ~132-175, validation ~1301-1350, output dispatch ~1558-1604)

**Interfaces:**
- Consumes: `AnnotateVcfConfig.workers` from Task 1.
- Produces: `VepConcurrencyPlan { lookup_partitions, spawn_vcf_provider_open }` (no `chromosome_lanes`/`inline_lookup`); a single `if config.workers > 1 { sharded } else { serial }` output dispatch.

- [ ] **Step 1: Replace `VepConcurrencyPlan` struct + impl (vcf_sink.rs 132-175)**

Replace the entire block (struct 132-138 and `impl` 140-175) with:

```rust
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
struct VepConcurrencyPlan {
    lookup_partitions: usize,
    spawn_vcf_provider_open: bool,
}

impl VepConcurrencyPlan {
    fn from_config(config: &AnnotateVcfConfig) -> Self {
        // Single `workers` knob: N position-ordered lookup partitions feeding
        // N independent annotation range pipelines, one contig at a time.
        // Correctness for workers>1 relies on the indexed-input guard in
        // annotate_to_vcf (contiguous, position-ordered partitions).
        let lookup_partitions = config.workers.max(1);
        Self {
            lookup_partitions,
            spawn_vcf_provider_open: lookup_partitions > 1,
        }
    }
}
```

- [ ] **Step 2: Update validation in `annotate_to_vcf` (vcf_sink.rs 1301-1337)**

Replace:

```rust
    if config.workers == 0 {
        return Err(DataFusionError::Plan(
            "annotate_to_vcf(): workers must be a positive integer".to_string(),
        ));
    }
    if config.target_partitions == 0 {
        return Err(DataFusionError::Plan(
            "annotate_to_vcf(): target_partitions must be a positive integer".to_string(),
        ));
    }
    if config.refseq || config.merged {
        return Err(DataFusionError::Plan(
            "annotate_to_vcf(): refseq and merged config fields are unsupported; cache source mode must come from cache schema metadata bio.vep.cache_source_type".to_string(),
        ));
    }
    // Parallel annotation (threads>1) splits the contig into N position-range
    // lookup partitions. That requires a bgzipped + tabix-indexed input so the
    // VCF scan yields contiguous, position-ordered partitions; an unindexed
    // input would force a round-robin repartition that scrambles output order.
    if config.threads > 1
        && !std::path::Path::new(&format!("{input_vcf}.tbi")).exists()
        && !std::path::Path::new(&format!("{input_vcf}.csi")).exists()
    {
        return Err(DataFusionError::Plan(format!(
            "annotate_to_vcf(): threads>1 requires a tabix-indexed input (`{input_vcf}.tbi` or `.csi`); \
             bgzip + tabix the VCF, or run with threads=1"
        )));
    }
    // threads>1 (within-contig sharded output) and forks>0 (contig-lane
    // parallelism) are distinct, mutually exclusive knobs. Combining them would
    // give the engine multiple output partitions, each restarting shard ids
    // from 0 and colliding on shard filenames. Reject early with a clear message.
    if config.threads > 1 && config.forks.is_some_and(|forks| forks > 0) {
        return Err(DataFusionError::Plan(
            "annotate_to_vcf(): threads>1 and forks>0 are mutually exclusive".to_string(),
        ));
    }
```

with:

```rust
    if config.workers == 0 {
        return Err(DataFusionError::Plan(
            "annotate_to_vcf(): workers must be a positive integer".to_string(),
        ));
    }
    if config.target_partitions == 0 {
        return Err(DataFusionError::Plan(
            "annotate_to_vcf(): target_partitions must be a positive integer".to_string(),
        ));
    }
    if config.refseq || config.merged {
        return Err(DataFusionError::Plan(
            "annotate_to_vcf(): refseq and merged config fields are unsupported; cache source mode must come from cache schema metadata bio.vep.cache_source_type".to_string(),
        ));
    }
    // Parallel annotation (workers>1) splits the contig into N position-range
    // lookup partitions. That requires a bgzipped + tabix-indexed input so the
    // VCF scan yields contiguous, position-ordered partitions; an unindexed
    // input would force a round-robin repartition that scrambles output order.
    if config.workers > 1
        && !std::path::Path::new(&format!("{input_vcf}.tbi")).exists()
        && !std::path::Path::new(&format!("{input_vcf}.csi")).exists()
    {
        return Err(DataFusionError::Plan(format!(
            "annotate_to_vcf(): workers>1 requires a tabix-indexed input (`{input_vcf}.tbi` or `.csi`); \
             bgzip + tabix the VCF, or run with workers=1"
        )));
    }
```

- [ ] **Step 3: Update the profile log (vcf_sink.rs 1340-1350)**

Replace:

```rust
        eprintln!(
            "[VEP_PROFILE] concurrency_plan lookup_partitions={} chromosome_lanes={} workers={} cold_parquet_target_partitions={} inline_lookup={} spawn_vcf_provider_open={}",
            concurrency_plan.lookup_partitions,
            concurrency_plan.chromosome_lanes,
            config.workers,
            config.target_partitions,
            concurrency_plan.inline_lookup,
            concurrency_plan.spawn_vcf_provider_open
        );
```

with:

```rust
        eprintln!(
            "[VEP_PROFILE] concurrency_plan lookup_partitions={} workers={} cold_parquet_target_partitions={} spawn_vcf_provider_open={}",
            concurrency_plan.lookup_partitions,
            config.workers,
            config.target_partitions,
            concurrency_plan.spawn_vcf_provider_open
        );
```

- [ ] **Step 4: Collapse the output dispatch (vcf_sink.rs 1558-1604)**

Replace the condition `if config.threads > 1 {` (line 1558) with `if config.workers > 1 {`.

Then delete the entire `else if concurrency_plan.chromosome_lanes > 1 { ... }` arm (lines 1590-1604, from `} else if concurrency_plan.chromosome_lanes > 1 {` up to but not including the final `} else {`), so the structure becomes `if config.workers > 1 { <sharded> } else { <serial stream> }`. The verbatim arm to delete:

```rust
    } else if concurrency_plan.chromosome_lanes > 1 {
        total_rows = stream_partitioned_vcf_body(
            df,
            &mut writer,
            &pb,
            config,
            total_input,
            Arc::new(vcf_info_fields),
            Arc::new(unique_format_tags),
            Arc::new(sample_names),
            coordinate_zero_based,
            concurrency_plan.chromosome_lanes,
            &mut sink_profile,
        )
        .await?;
```

(Leave the preceding `if config.workers > 1 { ... drive_sharded_vcf_annotation ... }` body and the trailing `} else { ... serial stream ... }` intact.)

- [ ] **Step 5: Commit checkpoint is deferred** — the crate does not compile until Task 3-6 land (the forks fan-out functions and `chromosome_lanes` field still exist). Proceed directly to Task 3.

---

## Task 3: Delete the forks fan-out functions + `PartitionedVcfBody`

**Files:**
- Modify: `datafusion/bio-function-vep/src/vcf_sink.rs` (delete `PartitionedVcfBody` ~85-95 and the forks-only functions in the ~431-797 range)

**Interfaces:**
- Consumes: nothing new. After Task 2 deleted the only caller (`stream_partitioned_vcf_body` at the removed dispatch arm), these become unreferenced.

- [ ] **Step 1: Delete the forks-only functions**

Delete these functions in full (verify each is unreferenced via the grep in Step 2 first if unsure). Per the dead-code map they are contiguous-ish in `vcf_sink.rs`:
- `write_vcf_partition_body` (~431-498)
- `schedule_partition_body_jobs` (~500-552) — contains the `JoinSet` fan-out
- `spawn_partition_body_scheduler` (~554-580)
- `stream_vcf_partition_to_writer` (~581-682)
- `stream_partitioned_vcf_body` (~683-797)
- any `partition_body_job_limit` helper (referenced only by the deleted scheduler + its test)

Also delete the `PartitionedVcfBody` struct (vcf_sink.rs 84-95):

```rust
#[derive(Debug)]
struct PartitionedVcfBody {
    partition_id: usize,
    path: PathBuf,
    batches: usize,
    input_rows: usize,
    lines: usize,
    bytes: usize,
    stream_next: Duration,
    format_duration: Duration,
    temp_write_duration: Duration,
}
```

**Keep** `VcfBodyTempDir` (97-130) — it is used by the surviving `workers>1` sharded path.

- [ ] **Step 2: Confirm nothing references the deleted symbols**

Run: `cd /Users/mwiewior/research/git/datafusion-bio-functions && grep -nE "stream_partitioned_vcf_body|schedule_partition_body_jobs|spawn_partition_body_scheduler|write_vcf_partition_body|stream_vcf_partition_to_writer|PartitionedVcfBody|partition_body_job_limit" datafusion/bio-function-vep/src/vcf_sink.rs`
Expected: matches only inside the now-removed unit tests (handled in Task 6) — no production references. If any production reference remains, it was missed in Task 2; fix before continuing.

---

## Task 4: Provider parse + `ContigAnnotationConfig` + exec partitioning

**Files:**
- Modify: `datafusion/bio-function-vep/src/annotate_provider.rs` (parse 5229-5276, config build 5320-5354, struct 8388-8400, exec `new` 8459-8475 + second call site ~8550-8561, `partition_contigs_for_execution` 8565-8576, test ctor 12506-12512)

**Interfaces:**
- Consumes: options-JSON key `"workers"` (from Task 1's `to_options_json` and vepyr Task 7).
- Produces: `ContigAnnotationConfig` with `target_partitions`, `annotation_workers` (renamed), `vcf_shard_ctx` — no `chromosome_lanes`, no `inline_lookup`. `ContigAnnotationExec::new` no longer takes a `chromosome_lanes` arg (always single output partition).

- [ ] **Step 1: Replace the parse block (annotate_provider.rs 5229-5276)**

Replace:

```rust
        let worker_forks = self
            .options_json
            .as_deref()
            .and_then(|opts| Self::parse_json_i64_option(opts, "forks"))
            .and_then(|value| usize::try_from(value).ok())
            .unwrap_or(0);
        let target_partitions = self
            .options_json
            .as_deref()
            .and_then(|opts| Self::parse_json_i64_option(opts, "target_partitions"))
            .and_then(|value| usize::try_from(value).ok())
            .filter(|value| *value > 0)
            .unwrap_or(1);
        let chromosome_lanes = self
            .options_json
            .as_deref()
            .and_then(|opts| Self::parse_json_i64_option(opts, "contig_parallelism"))
            .and_then(|value| usize::try_from(value).ok())
            .filter(|value| *value > 0)
            .unwrap_or_else(|| if worker_forks > 0 { worker_forks } else { 1 });
        let chromosome_lanes = if fetch_limit.is_some() {
            1
        } else {
            chromosome_lanes
        };
        let inline_lookup = self
            .options_json
            .as_deref()
            .and_then(|opts| Self::parse_json_bool_option(opts, "inline_lookup"))
            .unwrap_or(worker_forks == 0);
        // Number of parallel window workers for the CPU-bound annotation step
        // within a contig. Additive to the partition-lookup knobs above:
        // `threads <= 1` keeps the serial inline annotation path.
        let annotation_threads = self
            .options_json
            .as_deref()
            .and_then(|opts| Self::parse_json_i64_option(opts, "threads"))
            .and_then(|value| usize::try_from(value).ok())
            .filter(|value| *value > 0)
            .unwrap_or(1);
        // The single `threads` knob drives BOTH stages: with threads>1 the
        // lookup runs as N spawned, position-ordered partitions (so it isn't a
        // serial bottleneck) and annotation runs as N parallel window workers.
        let (inline_lookup, target_partitions) = if annotation_threads > 1 {
            (false, annotation_threads.max(target_partitions))
        } else {
            (inline_lookup, target_partitions)
        };
```

with:

```rust
        let target_partitions = self
            .options_json
            .as_deref()
            .and_then(|opts| Self::parse_json_i64_option(opts, "target_partitions"))
            .and_then(|value| usize::try_from(value).ok())
            .filter(|value| *value > 0)
            .unwrap_or(1);
        // The single `workers` knob drives BOTH stages: with workers>1 the
        // lookup runs as N spawned, position-ordered partitions (so it isn't a
        // serial bottleneck) and annotation runs as N parallel window workers.
        let annotation_workers = self
            .options_json
            .as_deref()
            .and_then(|opts| Self::parse_json_i64_option(opts, "workers"))
            .and_then(|value| usize::try_from(value).ok())
            .filter(|value| *value > 0)
            .unwrap_or(1);
        let target_partitions = if annotation_workers > 1 {
            annotation_workers.max(target_partitions)
        } else {
            target_partitions
        };
```

- [ ] **Step 2: Update the `ContigAnnotationConfig` construction (annotate_provider.rs 5340-5342)**

Replace:

```rust
            chromosome_lanes,
            inline_lookup,
            annotation_threads,
```

with:

```rust
            annotation_workers,
```

- [ ] **Step 3: Update the `ContigAnnotationExec::new` call in the parse path (annotate_provider.rs ~5350-5354)**

The call currently passes `chromosome_lanes` as the 4th argument:

```rust
        let exec = ContigAnnotationExec::new(
            projected_schema,
            self.schema.clone(),
            contigs,
            chromosome_lanes,
```

Remove the `chromosome_lanes,` argument line so the call is `new(projected_schema, self.schema.clone(), contigs, <session>, ...)` matching the new signature in Step 5.

- [ ] **Step 4: Edit the `ContigAnnotationConfig` struct fields (annotate_provider.rs 8393-8400)**

Replace:

```rust
    /// Maximum number of active chromosome lanes in the enclosing execution plan.
    chromosome_lanes: usize,
    /// Poll lookup streams inline instead of spawning lookup tasks.
    inline_lookup: bool,
    /// Number of parallel window workers for the annotation step within a
    /// contig. `<= 1` keeps the serial inline path. Additive to the
    /// partition-lookup knobs above (lookup partitioning is unchanged).
    annotation_threads: usize,
```

with:

```rust
    /// Number of parallel window workers for the annotation step within a
    /// contig. `<= 1` is serial (one spawned lookup partition).
    annotation_workers: usize,
```

- [ ] **Step 5: Edit `ContigAnnotationExec::new` signature + body (annotate_provider.rs 8460-8469)**

Replace:

```rust
    fn new(
        projected_schema: SchemaRef,
        full_schema: SchemaRef,
        contigs: Vec<String>,
        chromosome_lanes: usize,
        session: Arc<SessionContext>,
        cache: Arc<PartitionedAnnotationCache>,
        config: ContigAnnotationConfig,
    ) -> Self {
        let contig_partitions = partition_contigs_for_execution(contigs, chromosome_lanes);
```

with:

```rust
    fn new(
        projected_schema: SchemaRef,
        full_schema: SchemaRef,
        contigs: Vec<String>,
        session: Arc<SessionContext>,
        cache: Arc<PartitionedAnnotationCache>,
        config: ContigAnnotationConfig,
    ) -> Self {
        let contig_partitions = partition_contigs_for_execution(contigs);
```

- [ ] **Step 6: Confirm the second `ContigAnnotationExec::new` call site is already correct (annotate_provider.rs ~8554-8561)**

Read lines 8548-8562. The re-exec call site:

```rust
        ... ContigAnnotationExec::new(
            self.projected_schema.clone(),
            self.full_schema.clone(),
            contigs,
            Arc::clone(&self.session),
            Arc::clone(&self.cache),
            self.config.clone(),
        )))
```

If this call still passes a `chromosome_lanes`/lane-count argument before `Arc::clone(&self.session)`, remove that argument. (After Step 5 the signature has 6 params.)

- [ ] **Step 7: Simplify `partition_contigs_for_execution` (annotate_provider.rs 8565-8576)**

Replace:

```rust
fn partition_contigs_for_execution(
    contigs: Vec<String>,
    requested_parallelism: usize,
) -> Vec<Vec<String>> {
    if contigs.is_empty() {
        return vec![Vec::new()];
    }
    if requested_parallelism <= 1 {
        return vec![contigs];
    }
    contigs.into_iter().map(|contig| vec![contig]).collect()
}
```

with:

```rust
fn partition_contigs_for_execution(contigs: Vec<String>) -> Vec<Vec<String>> {
    if contigs.is_empty() {
        return vec![Vec::new()];
    }
    vec![contigs]
}
```

- [ ] **Step 8: Update the test constructor (annotate_provider.rs 12510-12512)**

Replace:

```rust
            chromosome_lanes: 1,
            inline_lookup: false,
            annotation_threads: 1,
```

with:

```rust
            annotation_workers: 1,
```

- [ ] **Step 9: Rename remaining `annotation_threads` references**

Run: `cd /Users/mwiewior/research/git/datafusion-bio-functions && grep -rn "annotation_threads" datafusion/bio-function-vep/src/`
For each remaining hit (e.g. the sharded gate ~11387 `if config.annotation_threads > 1`, the inflight pool ~11626 `config.annotation_threads.max(1)`), replace `annotation_threads` with `annotation_workers`. These are mechanical identifier renames; do not change logic.

---

## Task 5: Remove the serial inline-lookup path

**Files:**
- Modify: `datafusion/bio-function-vep/src/annotate_provider.rs` (handle struct ~9092-9098, enum ~9100-9111, `inline_lookup_partition` ~9879-9894, poll arm ~11101-11136, selection sites ~11962-11980 and ~11990-12007)

**Interfaces:**
- Consumes: nothing new. Forces both lookup-selection sites onto the spawned worker; `workers==1` → one partition → one `spawn_lookup_partition_worker`.

- [ ] **Step 1: Force spawned in the partitioned selection site (annotate_provider.rs 11962-11980)**

Replace:

```rust
            if config.inline_lookup {
                handles.push_back(inline_lookup_partition(
                    Arc::clone(&plan),
                    Arc::clone(&task_ctx),
                    partition_id,
                    chrom.to_string(),
                    sink,
                )?);
            } else {
                handles.push_back(spawn_lookup_partition_worker(
                    Arc::clone(&plan),
                    Arc::clone(&task_ctx),
                    partition_id,
                    chrom.to_string(),
                    sink,
                    LOOKUP_PARTITION_QUEUE_BATCHES,
                    pipeline_profile.clone(),
                ));
            }
```

with:

```rust
            handles.push_back(spawn_lookup_partition_worker(
                Arc::clone(&plan),
                Arc::clone(&task_ctx),
                partition_id,
                chrom.to_string(),
                sink,
                LOOKUP_PARTITION_QUEUE_BATCHES,
                pipeline_profile.clone(),
            ));
```

- [ ] **Step 2: Force spawned in the fallback selection site (annotate_provider.rs 11990-12007)**

Replace:

```rust
        if config.inline_lookup {
            VecDeque::from([LookupPartitionHandle::Inline(InlineLookupPartitionHandle {
                partition_id: 0,
                chrom: chrom.to_string(),
                stream: lookup_stream,
                sink: fallback_coloc_sink,
                next_batch_id: 0,
            })])
        } else {
            VecDeque::from([spawn_lookup_stream_worker(
                lookup_stream,
                plan.schema(),
                chrom.to_string(),
                fallback_coloc_sink,
                LOOKUP_PARTITION_QUEUE_BATCHES,
                pipeline_profile.clone(),
            )])
        }
```

with:

```rust
        VecDeque::from([spawn_lookup_stream_worker(
            lookup_stream,
            plan.schema(),
            chrom.to_string(),
            fallback_coloc_sink,
            LOOKUP_PARTITION_QUEUE_BATCHES,
            pipeline_profile.clone(),
        )])
```

- [ ] **Step 3: Delete `inline_lookup_partition` (annotate_provider.rs 9879-9894)**

Delete the whole function:

```rust
fn inline_lookup_partition(
    plan: Arc<dyn ExecutionPlan>,
    task_ctx: Arc<TaskContext>,
    partition_id: usize,
    chrom: String,
    sink: ColocatedSink,
) -> Result<LookupPartitionHandle> {
    let stream = plan.execute(partition_id, task_ctx)?;
    Ok(LookupPartitionHandle::Inline(InlineLookupPartitionHandle {
        partition_id,
        chrom,
        stream,
        sink,
        next_batch_id: 0,
    }))
}
```

- [ ] **Step 4: Remove the `Inline` enum variant and its `partition_id()` arm (annotate_provider.rs 9100-9110)**

Delete the variant line `Inline(InlineLookupPartitionHandle),` (9102) and the `LookupPartitionHandle::Inline(handle) => handle.partition_id,` arm (9109) in the `partition_id()` match.

- [ ] **Step 5: Delete the `InlineLookupPartitionHandle` struct (annotate_provider.rs 9092-9098)**

```rust
struct InlineLookupPartitionHandle {
    partition_id: usize,
    chrom: String,
    stream: SendableRecordBatchStream,
    sink: ColocatedSink,
    next_batch_id: usize,
}
```

- [ ] **Step 6: Delete the `Inline` poll arm (annotate_provider.rs 11101-11136)**

Delete the entire `LookupPartitionHandle::Inline(active) => match active.stream.as_mut().poll_next(cx) { ... },` arm. The `Spawned` arm remains and covers all handles.

- [ ] **Step 7: Build and resolve residual inline references**

Run: `cd /Users/mwiewior/research/git/datafusion-bio-functions && cargo build -p datafusion-bio-function-vep 2>&1 | grep -E "error|warning: unused|never" | head -40`
Resolve any remaining errors. Two likely items:
- The sharded `inline_seen` guard (~11411-11453): with no `Inline` handles, `take_spawned_parts()` never returns `None` for a real handle. If the compiler flags `inline_seen` as unused, remove the `let mut inline_seen = false;`, the `None => { handle.abort(); inline_seen = true; }` arm (change the `match` to bind the `Some((lookup_rx, lookup_jh))` parts directly via `if let`/expect), and the `if inline_seen { ... }` error block. If it does not warn, leaving it is acceptable (dead but harmless).
- `parse_json_bool_option`: run `grep -rn "parse_json_bool_option" datafusion/bio-function-vep/src/` — if it now has no callers, delete its definition; if it has other callers, keep it.

---

## Task 6: Remove engine tests for forks/scheduler; engine green

**Files:**
- Modify: `datafusion/bio-function-vep/src/vcf_sink.rs` (tests ~1832-1907 and ~1982-2018)

- [ ] **Step 1: Delete the forks/plan/serialization unit tests**

Delete these `#[test]` functions in full (they assert the removed `forks`/`chromosome_lanes`/`inline_lookup`/`contig_parallelism` behavior):
- `test_forks_zero_derives_strict_inline_concurrency_plan` (~1832-1845)
- `test_forks_one_preserves_legacy_single_partition_plan` (~1847-1860)
- `test_forks_nonzero_uses_fork_count_for_chromosome_lanes` (~1862-1876)
- `test_to_options_json_emits_forks_plan_when_set` (~1878-1890)
- `test_to_options_json_keeps_workers_and_chromosome_lanes_separate` (~1892-1907)
- `test_partitioned_body_writer_schedules_partition_zero_as_a_job` (~1982-1987)
- `test_partition_body_scheduler_uses_full_job_limit` (~1989-2018)

- [ ] **Step 2: Add a minimal plan test reflecting the new contract**

In the same `#[cfg(test)] mod` (where the deleted tests lived), add:

```rust
    #[test]
    fn concurrency_plan_serial_uses_single_lookup_partition() {
        let config = AnnotateVcfConfig::default();
        let plan = VepConcurrencyPlan::from_config(&config);
        assert_eq!(plan.lookup_partitions, 1);
        assert!(!plan.spawn_vcf_provider_open);
    }

    #[test]
    fn concurrency_plan_workers_drive_lookup_partitions() {
        let config = AnnotateVcfConfig {
            workers: 8,
            ..AnnotateVcfConfig::default()
        };
        let plan = VepConcurrencyPlan::from_config(&config);
        assert_eq!(plan.lookup_partitions, 8);
        assert!(plan.spawn_vcf_provider_open);
    }

    #[test]
    fn to_options_json_emits_workers_not_threads_or_forks() {
        let config = AnnotateVcfConfig {
            workers: 4,
            ..AnnotateVcfConfig::default()
        };
        let json: serde_json::Value =
            serde_json::from_str(&config.to_options_json()).unwrap();
        assert_eq!(json["workers"], 4);
        assert!(json.get("threads").is_none());
        assert!(json.get("forks").is_none());
        assert!(json.get("contig_parallelism").is_none());
        assert!(json.get("inline_lookup").is_none());
    }
```

- [ ] **Step 3: Engine full gate**

Run: `cd /Users/mwiewior/research/git/datafusion-bio-functions && cargo test -p datafusion-bio-function-vep 2>&1 | tail -20 && cargo clippy -p datafusion-bio-function-vep --all-features -- -D warnings 2>&1 | tail -10 && cargo fmt -- --check`
Expected: tests pass, no clippy warnings, no fmt diff. (If `fmt --check` reports a diff, run `cargo fmt`.)

- [ ] **Step 4: Commit the engine cleanup**

```bash
cd /Users/mwiewior/research/git/datafusion-bio-functions
git add datafusion/bio-function-vep/src/vcf_sink.rs datafusion/bio-function-vep/src/annotate_provider.rs
git commit -m "refactor(vep): single workers knob; remove forks/chromosome_lanes + inline lookup

Delete the cross-contig forks/chromosome-lane fan-out path and the serial
inline-lookup path; rename the within-contig knob threads -> workers
(config field, options-JSON key, annotation_workers). Serial now routes
through one spawned lookup partition. vepyr is the sole consumer (local
path dep)."
```

---

## Task 7: vepyr Rust glue → single `workers` knob against the renamed engine

**Files:**
- Modify: `vepyr/src/annotate.rs`, `vepyr/src/lib.rs`, `vepyr/src/vepyr/_core.pyi`

**Interfaces:**
- Consumes: engine `AnnotateVcfConfig { workers, target_partitions, .. }` (no `forks`/old-`workers`/`threads`) and options-JSON key `"workers"`.
- Produces: `_core.annotate_vcf(..., on_batch_written=None)` and `_core.create_annotator(..., limit=None)` (no `forks`/`workers` positional args); both read the knob from the `"workers"` options-JSON key.

- [ ] **Step 1: Rewrite `annotate.rs` helpers (vepyr/src/annotate.rs 15-105)**

Replace the four helpers (`effective_session_partitions`, `effective_runtime_threads`, `runtime_for_parallelism`, `options_json_with_parallelism`) with:

```rust
fn worker_thread_count(workers: usize) -> usize {
    workers.max(1)
}

fn runtime_for_workers(workers: usize) -> PyResult<Arc<Runtime>> {
    let runtime = Builder::new_multi_thread()
        .worker_threads(worker_thread_count(workers))
        .build()
        .map_err(|e| pyo3::exceptions::PyRuntimeError::new_err(format!("{e}")))?;
    Ok(Arc::new(runtime))
}

fn normalize_options(options_json: &str) -> PyResult<(String, String)> {
    let mut opts: Value = serde_json::from_str(options_json).map_err(|e| {
        pyo3::exceptions::PyValueError::new_err(format!("Invalid options JSON: {e}"))
    })?;
    let object = opts
        .as_object_mut()
        .ok_or_else(|| pyo3::exceptions::PyValueError::new_err("options JSON must be an object"))?;
    let cache_format = object
        .get("cache_format")
        .and_then(|v| v.as_str())
        .unwrap_or("indexed_parquet")
        .to_string();
    if !matches!(
        cache_format.as_str(),
        "indexed_parquet" | "legacy_fjall" | "lance"
    ) {
        return Err(pyo3::exceptions::PyValueError::new_err(
            "cache_format must be 'indexed_parquet', 'legacy_fjall', or 'lance'",
        ));
    }
    object.insert("cache_format".to_string(), Value::from(cache_format.clone()));
    let options_json = serde_json::to_string(&opts).map_err(|e| {
        pyo3::exceptions::PyValueError::new_err(format!("Invalid options JSON: {e}"))
    })?;
    Ok((options_json, cache_format))
}

fn workers_from_options(opts: &Value) -> usize {
    opts.get("workers")
        .and_then(|v| v.as_u64())
        .and_then(|n| usize::try_from(n).ok())
        .filter(|n| *n > 0)
        .unwrap_or(1)
}
```

- [ ] **Step 2: Rewrite `annotate_to_vcf_file` (vepyr/src/annotate.rs 166-350)**

Change the signature to drop `forks`/`workers` (end it `on_batch_written: Option<Py<PyAny>>,) -> PyResult<usize> {`). Replace the options/runtime preamble (187-197) with:

```rust
    let (options_json, cache_format) = normalize_options(options_json)?;
    let opts: Value = serde_json::from_str(&options_json).map_err(|e| {
        pyo3::exceptions::PyValueError::new_err(format!("Invalid options JSON: {e}"))
    })?;
    let workers = workers_from_options(&opts);

    let backend = cache_format.as_str();
    let rt = runtime_for_workers(workers)?;
```

In the `AnnotateVcfConfig { ... }` literal, replace the parallelism fields (314-325: `forks`, `workers`, `threads`, `target_partitions`) with:

```rust
        // Single annotation-concurrency knob.
        workers,
        target_partitions: 1,
```

- [ ] **Step 3: Rewrite `create_streaming_annotator` (vepyr/src/annotate.rs 353-433)**

Change the signature to drop `forks`/`workers` (end it `limit: Option<usize>,) -> PyResult<StreamingAnnotator> {`). Replace the preamble (364-376) with:

```rust
    let (options_json, cache_format) = normalize_options(options_json)?;
    let opts: Value = serde_json::from_str(&options_json).map_err(|e| {
        pyo3::exceptions::PyValueError::new_err(format!("Invalid options JSON: {e}"))
    })?;
    let workers = workers_from_options(&opts);
    let rt = runtime_for_workers(workers)?;

    let (stream, schema) = rt.block_on(async {
        let backend = cache_format.as_str();
        let session_partitions = worker_thread_count(workers);
```

(Delete the now-removed `let _opts: Value = ...` line that followed the old options call.)

- [ ] **Step 4: Rewrite the `annotate.rs` unit tests (vepyr/src/annotate.rs 435-492)**

Replace the entire `#[cfg(test)] mod tests` block with:

```rust
#[cfg(test)]
mod tests {
    use super::{normalize_options, worker_thread_count, workers_from_options};

    #[test]
    fn normalize_preserves_workers_and_cache_format() {
        let (json, fmt) =
            normalize_options(r#"{"cache_format":"lance","workers":4}"#).unwrap();
        let opts: serde_json::Value = serde_json::from_str(&json).unwrap();
        assert_eq!(fmt, "lance");
        assert_eq!(opts["workers"], 4);
        assert_eq!(workers_from_options(&opts), 4);
    }

    #[test]
    fn default_cache_format_and_workers_when_absent() {
        let (_json, fmt) = normalize_options(r#"{}"#).unwrap();
        assert_eq!(fmt, "indexed_parquet");
        let opts: serde_json::Value = serde_json::from_str(r#"{}"#).unwrap();
        assert_eq!(workers_from_options(&opts), 1);
    }

    #[test]
    fn invalid_cache_format_is_rejected() {
        pyo3::Python::initialize();
        let err = normalize_options(r#"{"cache_format":"fjall"}"#).unwrap_err();
        assert!(err.to_string().contains("cache_format"));
    }

    #[test]
    fn worker_thread_count_is_at_least_one() {
        assert_eq!(worker_thread_count(0), 1);
        assert_eq!(worker_thread_count(8), 8);
    }
}
```

- [ ] **Step 5: Update PyO3 signatures (vepyr/src/lib.rs 89-144)**

Replace `annotate_vcf` and `create_annotator` with the no-`forks`/`workers` versions (identical to Task 1 Steps 6 of the superseded vepyr-only plan):

```rust
#[pyfunction]
#[pyo3(signature = (vcf_path, cache_dir, output_path, options_json, show_progress=true, compression="", on_batch_written=None))]
#[allow(clippy::too_many_arguments)]
fn annotate_vcf(
    py: Python<'_>,
    vcf_path: &str,
    cache_dir: &str,
    output_path: &str,
    options_json: &str,
    show_progress: bool,
    compression: &str,
    on_batch_written: Option<Py<PyAny>>,
) -> PyResult<usize> {
    annotate::annotate_to_vcf_file(
        py, vcf_path, cache_dir, output_path, options_json, show_progress, compression,
        on_batch_written,
    )
}

#[pyfunction]
#[pyo3(signature = (vcf_path, cache_dir, options_json, skip_csq=true, limit=None))]
fn create_annotator(
    py: Python<'_>,
    vcf_path: &str,
    cache_dir: &str,
    options_json: &str,
    skip_csq: bool,
    limit: Option<usize>,
) -> PyResult<annotate::StreamingAnnotator> {
    annotate::create_streaming_annotator(py, vcf_path, cache_dir, options_json, skip_csq, limit)
}
```

- [ ] **Step 6: Update `_core.pyi`** — drop `forks`/`workers` from the `annotate_vcf` and `create_annotator` stubs (as in the superseded plan's Task 1 Step 7).

- [ ] **Step 7: Build the vepyr Rust side against the local engine**

Run: `cd /Users/mwiewior/research/git/vepyr && cargo test --lib annotate 2>&1 | tail -20 && cargo clippy --lib -- -D warnings 2>&1 | tail -10 && cargo fmt -- --check`
Expected: 4 annotate tests pass; clippy/fmt clean.

---

## Task 8: vepyr Python + CLI + docs + e2e gate

**Files:**
- Modify: `vepyr/src/vepyr/__init__.py`, `tests/test_annotate.py`, `tests/test_build_cache.py`, `e2e-testing/scripts/run_annotation_fast.py`, `run_annotation_fast_all.py`, `README.md`, `docs/quickstart.md`, `docs/performance.md`

**Interfaces:**
- Consumes: the rebuilt `_core` from Task 7 and the renamed `"workers"` options key.

- [ ] **Step 1: Rebuild the extension**

Run: `cd /Users/mwiewior/research/git/vepyr && maturin develop 2>&1 | tail -5`

- [ ] **Step 2: Apply the Python `annotate()` + tests + CLI + docs edits**

Apply Task 2 Steps 2-9, Task 3 (all steps), and Task 4 (all steps) **from the superseded plan** `vepyr/docs/superpowers/plans/2026-06-21-single-workers-knob.md`, with this single difference: wherever that plan sets the options key `opts["threads"] = workers`, use `opts["workers"] = workers` instead (the engine now parses `"workers"`). Concretely in `src/vepyr/__init__.py`:

```python
    if workers > 1:
        # Single annotation-concurrency knob: N within-contig fused pipelines.
        # Requires a tabix-indexed (bgzip+.tbi) input VCF.
        opts["workers"] = workers
```

And in `tests/test_annotate.py`, the assertions that the superseded plan writes as `seen[...]["threads"] == 4` become `seen[...]["workers"] == 4`.

- [ ] **Step 3: Python test gate**

Run: `cd /Users/mwiewior/research/git/vepyr && python -m pytest tests/test_annotate.py tests/test_run_annotation_fast.py 2>&1 | tail -20`
Expected: PASS (golden-cache integration tests may skip).

- [ ] **Step 4: Stray-reference sweep across both repos**

Run:
```bash
cd /Users/mwiewior/research/git/vepyr && grep -rnE "forks|--threads|\"threads\"|opts\[.threads.\]" --include="*.py" --include="*.rs" --include="*.md" . | grep -viE "worker_threads|rt-multi-thread|num_threads" | grep -v "/target/" | grep -v "docs/superpowers/plans/"
cd /Users/mwiewior/research/git/datafusion-bio-functions && grep -rnE "\.threads|\bforks\b|chromosome_lanes|inline_lookup|contig_parallelism|annotation_threads" datafusion/bio-function-vep/src/ | grep -viE "worker_threads|num_threads|os_threads"
```
Expected: no functional matches (only unrelated `worker_threads` / doc-plan references). Fix any straggler.

- [ ] **Step 5: chr1 e2e byte-identical gate**

Run:
```bash
cd /Users/mwiewior/research/git/vepyr
maturin develop --release 2>&1 | tail -3
python e2e-testing/scripts/run_annotation_fast.py chr1 --backend lance --skip-compare --force 2>&1 | tail -15
# Replace <CHR1_OUTPUT> with the path printed above (same as Task 0 Step 3).
grep -v '^##' <CHR1_OUTPUT> | md5 | tee /tmp/chr1_after.md5
diff /tmp/chr1_baseline.md5 /tmp/chr1_after.md5 && echo "BYTE-IDENTICAL (serial parity OK)"
```
Expected: `BYTE-IDENTICAL (serial parity OK)` — the refactor preserved serial output exactly.

- [ ] **Step 6: workers>1 parity gate**

Run (the chr1 input is bgzipped+tabixed, so `workers>1` is allowed):
```bash
cd /Users/mwiewior/research/git/vepyr
python e2e-testing/scripts/run_annotation_fast.py chr1 --backend lance --workers 4 --skip-compare --force 2>&1 | tail -10
# Replace <CHR1_OUTPUT> with the path printed above.
grep -v '^##' <CHR1_OUTPUT> | md5 | tee /tmp/chr1_w4.md5
diff /tmp/chr1_baseline.md5 /tmp/chr1_w4.md5 && echo "BYTE-IDENTICAL (workers=4 parity OK)"
```
Expected: `BYTE-IDENTICAL (workers=4 parity OK)` — parallel output matches serial.

- [ ] **Step 7: Full lint + test gate (both repos)**

Run:
```bash
cd /Users/mwiewior/research/git/datafusion-bio-functions && cargo test 2>&1 | tail -8 && cargo clippy --all-features -- -D warnings 2>&1 | tail -5 && cargo fmt -- --check
cd /Users/mwiewior/research/git/vepyr && cargo clippy -- -D warnings 2>&1 | tail -5 && cargo fmt -- --check && python -m pytest tests/ 2>&1 | tail -12
```
Expected: all green (pre-existing skips allowed).

- [ ] **Step 8: Commit vepyr**

```bash
cd /Users/mwiewior/research/git/vepyr
git add src/ tests/ e2e-testing/scripts/run_annotation_fast.py e2e-testing/scripts/run_annotation_fast_all.py README.md docs/quickstart.md docs/performance.md
git commit -m "refactor(vepyr): single --workers knob against renamed engine

Collapse forks/workers/threads to one --workers knob; emit options-JSON
'workers' key consumed by the cleaned engine. chr1 byte-identical
(serial and workers=4)."
```

---

## Self-Review

- **Spec coverage:** remove forks/chromosome_lanes (Tasks 2-4), remove inline path (Task 5), rename threads→workers across engine field + JSON key + annotation_workers (Tasks 1,4), remove old per-contig workers (Task 1), vepyr single knob + emit `"workers"` (Tasks 7-8), docs (Task 8), no re-pin needed (local path dep — Global Constraints), byte-identical gate serial + workers>1 (Tasks 0,8). `partitions` untouched. Covered.
- **Type/name consistency:** engine field `workers` ← `threads`; JSON key `"workers"` ← `"threads"` (emitted in `to_options_json` Task 1 Step 4, parsed in Task 4 Step 1, produced by vepyr Task 8 Step 2); internal `annotation_workers` ← `annotation_threads` (Task 4 Steps 1,2,4,8,9). `VepConcurrencyPlan` fields reduced to `lookup_partitions`/`spawn_vcf_provider_open` (Task 2) and the profile log + dispatch reference only those (Task 2 Steps 3,4). `partition_contigs_for_execution` arity reduced from 2→1 and both call sites updated (Task 4 Steps 3,5,6,7). vepyr `workers_from_options`/`normalize_options`/`runtime_for_workers`/`worker_thread_count` defined Task 7 Step 1, used Steps 2-3, tested Step 4.
- **Placeholder scan:** edits give verbatim old→new; the two read-then-edit spans (second `ContigAnnotationExec::new` call site Task 4 Step 6; the `inline_seen` guard Task 5 Step 7) are bounded with exact symbols and a build-error driver. No TBD/"handle errors".
- **Risk note:** the `inline_seen`/`take_spawned_parts` simplification (Task 5 Step 7) is the only judgement call; the fallback ("leave dead but harmless") is stated, and the byte-identical gate (Task 8 Steps 5-6) catches any behavioral regression.
