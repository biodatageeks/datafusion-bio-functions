# CADD tier-join reservation leak — RESOLVED (2026-08-25)

> **Status: fixed in `92f55f9`.** All 5 plugins build chr21 and validate.
> CADD: 121,809,062 rows in 2.1 min, `HashJoinInput` peak **0.70 GB** (was
> 8.02 GB at an 8 GiB pool and 20.00 GB at 20 GiB). Kept as a record of the
> mechanism and of what was ruled out — §5 in particular is worth reading
> before touching this join again.

## 0. What it was, and the piece that was missing

The tier join reserved **exactly whatever memory pool it was given** to collect
a build side of **~180 MB**, while RSS stayed flat — memory reserved but never
allocated, so no pool size helped.

The cause was `get_record_batch_memory_size` charging a *shared* Arrow buffer
once per batch that references it (it de-duplicates only *within* one batch).
Two layers of sharing had to be removed, and only finding both fixes it:

1. **`GROUP BY` output batches share buffers** (2d6c01d). Measured 31.7x
   over-count — see §4. Fixed by materializing the probe with a `concat` per
   batch. This dropped the increment from 592.7 MB to 18.2 MB but did **not**
   fix the build.
2. **The MemTable scan re-splits oversized batches** (92f55f9). DataFusion 53
   wraps a MemTable in `BatchSplitStream`, which slices any parent larger than
   the execution batch size into **zero-copy children**. The materialized probe
   was chunked at 1M rows, so its 11 owned parents became ~1300 slices that all
   shared 11 buffers, and the over-count returned — that is precisely what the
   18.2 MB increment was: one parent's buffer, charged per 8192-row slice.
   Fixed by chunking to `ctx.copied_config().batch_size()` so the scan has
   nothing left to split, plus forcing a real copy in
   `concat_probe_parts_owned` (Arrow's single-input `concat_batches` is
   zero-copy and would otherwise keep the sharing).

The lesson worth keeping: **an over-count that shrinks but does not vanish means
a second sharing layer, not a partial fix.** The 18.2 MB figure was measured
and recorded a full round before anyone recognised it as "one parent per
execution-sized slice".

## 0b. The test that killed hypotheses

Any explanation had to account for a large reservation/RSS gap. Three plausible
theories died on it before the real one survived — see §5. Measure before
theorising; this bug was unusually good at producing convincing wrong answers.

---

## 1. Locations

| Thing | Path / ref |
|---|---|
| Engine worktree (the bug lives here) | `datafusion-bio-functions`, branch `plugins-fixes`, PR #217, HEAD `2d6c01d` |
| Tier join | `datafusion/bio-function-vep/src/plugin_cache/join.rs` |
| Build driver + retry | `datafusion/bio-function-vep/src/plugin_cache/build.rs` |
| Memory tracer | `datafusion/bio-function-vep/src/plugin_cache/mem_trace.rs` |
| vepyr (PyO3 consumer, pins the engine by rev) | `vepyr`, branch `plugin-cache-build-validation`, PR #49 |
| Plugin manifests | `vepyr-plugins`, branch `fix/adding-a-plugin-sort-guidance`, PR #4 |
| Repro scripts | `scripts/debug/` (this repo) |
| Variation cache | `~/workspace/data_vepyr/cache/116_GRCh38_merged` |
| Raw plugin sources | Google Drive folder `1ZT3g31I0LXepORF_dusy47AuXOnbRlZ_`, rclone remote `gdrive-mw` |

`vepyr` pins the engine by git rev, so **every engine change needs a rev bump +
`uv sync --reinstall-package vepyr`** (~3 min). Bare `cargo build` on vepyr fails
to link (PyO3, no Python symbols) — that is expected, not a regression.

---

## 2. Reproduce

### 2.1 Sources without downloading 168 GB

All six source files are bgzip+tabix indexed. `rclone serve http` exposes Drive
over local HTTP with range support, so tabix pulls one chromosome out of an
87 GB file in ~50 s:

```bash
rclone serve http gdrive-mw: \
  --drive-root-folder-id 1ZT3g31I0LXepORF_dusy47AuXOnbRlZ_ \
  --addr localhost:18080 --read-only &
curl -s -o /dev/null -w '%{http_code}\n' http://localhost:18080/   # expect 200
```

Contig naming differs per source and a wrong prefix yields a silently empty
slice: **`chr`-prefixed for alphamissense, bare for the other four.** Verify with
`tabix -l <url> | head -3`.

### 2.2 Slice chr21

```bash
S=~/workspace/data_vepyr/plugin_input/_slices; mkdir -p $S
tabix http://localhost:18080/cadd/whole_genome_SNVs.tsv.gz 21 > $S/cadd_chr21_snv.tsv
tabix http://localhost:18080/cadd/gnomad.genomes.r4.0.indel.tsv.gz 21 > $S/cadd_chr21_indel.tsv
```

~3.5 GB + 58 MB, ~60 s total.

### 2.3 Build (fails)

```bash
cd <vepyr worktree>
VIRTUAL_ENV= RUST_LOG=info VEPYR_TRACE_MEMORY=1 VEPYR_EXPLAIN_TIER_JOIN=1 \
  ./.venv/bin/python <engine>/scripts/debug/build_plugin_chrom.py cadd 21
```

Expected failure:

```
plugin_cache[cadd/21]: tier join did not preserve row order (121809062 plugin
                       rows this chrom) -- retrying with an explicit sort
Resources exhausted: Failed to allocate additional 18.2 MB for HashJoinInput
with 8.0 GB already allocated - 1420.1 KB remain available for the total pool
```

`RUST_LOG=info` is required or the instrumentation prints nothing.

### 2.4 Turn the dials

```bash
VEPYR_PLUGIN_RETRY_MEMORY_MIB=20480   # pool size; it will take exactly this much
VEPYR_TRACE_MEMORY=1                  # per-reservation peaks (mem_trace.rs)
VEPYR_EXPLAIN_TIER_JOIN=1             # physical plan (join.rs)
```

---

## 3. Established, with evidence

**The retry path is where it fails, and the retry is not the cause.**
`tiered_stream` and `tiered_stream_sorted` are the same `tiered_stream_impl`
differing only by `ORDER BY` (`join.rs:41-58`). The fast path
(`build.rs`, `SessionContext::new_with_config`) uses the default **unbounded**
pool; the retry (`new_with_config_rt`) is the only path with a bound. So the
retry does not use more memory — it is the only path that counts.

CADD reaches the retry legitimately: its two `[[source]]` parts are `UNION ALL`ed,
so the join sees one sorted block after another and `assert_start_monotonic`
trips. That is the guard working.

**The collected side is the probe, and it is small.** Measured directly
(`scripts/debug/measure_probe_size.py`):

```
variation rows:      14,574,753
plugin rows (SNV):  120,265,857
PROBE VIEW rows:     10,806,229   <- 4 narrow columns, ~180 MB
```

Plan after the 2d6c01d fix:

```
HashJoinExec: mode=CollectLeft, join_type=Right
  <projection> -> DataSourceExec: partitions=1, partition_sizes=[11]   <- LEFT, COLLECTED (probe)
  DataSourceExec: partitions=1, partition_sizes=[14870]                <- RIGHT (plugin)
```

`CollectLeft` collects the left input. 11 batches, ~180 MB.

**The reservation tracks the ceiling, not the data.**

| pool | `HashJoinInput#3` peak | peak RSS |
|---:|---:|---:|
| 8 GiB | 8.02 GB | 16.9 GB |
| 20 GiB | **20.00 GB** (exactly the pool) | 16.9 GB |

Reserved memory exceeds RSS, so it is not backed by allocation.

**Only CADD is affected, because only CADD's plan swaps.** `join_type=Right`
appears because 120M plugin rows dwarf the 10.8M-row probe. dbNSFP plans
`join_type=Left` and collects its own rows — a genuine 1557 MB working set that
the 8 GiB default fixes properly.

---

## 4. Fixed along the way (keep these)

**`GROUP BY` cross-batch buffer over-count — 2d6c01d.**
`GroupedHashAggregateStream` emits batches sharing underlying buffers.
`HashJoinExec` reserves per batch via `get_record_batch_memory_size`, which
de-duplicates buffers only *within* one batch, so a shared buffer is charged once
per referencing batch. Measured (`scripts/debug/measure_batch_buffer_sharing.py`):

| source | batches | accounted | distinct buffers | over-count |
|---|---:|---:|---:|---:|
| raw parquet scan | 1793 | 277.9 MB | 277.9 MB | 1.0x |
| same scan + `GROUP BY` | 1791 | 5756.8 MB | 181.4 MB | **31.7x** |

Materializing the probe (a `concat` per batch → unshared buffers) dropped CADD's
increment from a repeated **592.7 MB** to a smooth **18.2 MB** and cut the probe
from 1791 batches to 11. Related upstream: arrow-rs#6439, which DataFusion's own
comment cites as why `get_record_batch_memory_size` exists — it does not cover
the cross-batch case.

**Retry pool budget — 39fc35c.** 2 GiB → 8 GiB default, plus
`VEPYR_PLUGIN_RETRY_MEMORY_MIB`. Fixes dbNSFP. Does not and cannot fix CADD.

**`schema_force_view_types=false` — aa877f3.** Halves the shared buffer
(592.7 → 304.2 MB) and drops `CAST(... AS Utf8View)` from the join keys. A real
improvement, but it was **not** the fix — it changed the buffer's width, not the
sharing.

---

## 5. Ruled out — do not re-tread

| Theory | Killed by |
|---|---|
| The `ORDER BY` sort is expensive | Retry plan == fast plan + `ORDER BY`; the error names `HashJoinInput`, not a sort reservation |
| The plugin is always the collected side | Plan dumps: `Left` for dbNSFP (plugin collected), `Right` for CADD (probe collected) |
| It is a scale problem | SpliceAI retries at 31.7M rows and succeeds; dbNSFP fails at 836k |
| A bigger pool fixes it | 2 / 8 / 20 GiB all fail; the reservation grows to match the pool exactly |
| `Utf8View` width is the cause | Disabling view types halved the increment, pattern unchanged |
| `GROUP BY` buffer sharing is the whole cause | Fixed in 2d6c01d; increment fell to 18.2 MB but the leak persisted — it was one of **two** sharing layers (§0) |
| The semi-join (b561629) is the consumer | `HashJoinInput#0` is 0.95 GB and **stable** across pool sizes |

---

## 6. Where to look next

`collect_left_input` (`datafusion-physical-plan-53.0.0/src/joins/hash_join/exec.rs`)
grows the reservation once per build batch (`:1874`) and once for the hash map
(`:~1935`), and `left_fut.try_once` (`:1299`) should make the collect happen once.
That accounts for ~180 MB and cannot account for 20 GB — so the extra growth
comes from somewhere else. Unexplored angles:

1. **Is `collect_left_input` running more than once?** Log entry/exit, or count
   `MemoryConsumer::new("HashJoinInput")` registrations. `mem_trace` keys on
   consumer identity, so a second collect would appear as a *new* `#n` — it does
   not, which argues against this but does not exclude a reused consumer.
2. **Is the probe side reserving under the build's consumer?** `stream.rs` is
   the probe loop; check whether anything there touches the build reservation.
3. **`join_type=Right` + `CollectLeft` specifically.** CADD is the only plugin
   that plans this combination and the only one that leaks. Force the swap off
   (make the probe look larger, or disable the optimizer's join-selection rule)
   and see whether the leak follows the swap. **This is the highest-value test.**
4. **Reproduce outside vepyr.** Drive the same plan from a small Rust test or
   `datafusion-cli` against `variation/chr21.parquet` + the CADD slice, under a
   `FairSpillPool`. If it reproduces, this is an upstream DataFusion bug and
   should be filed with that reproducer.

If it is upstream, the near-term workaround is to keep CADD off the retry path —
feed the join pre-sorted input so `assert_start_monotonic` never trips. Note the
correct sort is **not** `sort -k2,2n`: CADD's `ingest_sql` anchor-shifts indels to
`start = pos+1` while the SNVs at that position stay at `pos`, so the tie must
break on shifted-ness (unshifted first). See `vepyr-plugins` PR #4.

---

## 7. Working state

**chr21, built and validated** (`(tier, start)` ascending per shard; SpliceAI's
`assume_unique` checked exhaustively over all 31.7M keys):

| plugin | rows | shard | notes |
|---|---:|---:|---|
| clinvar | 49,937 | 1 MB | |
| alphamissense | 698,535 | 6 MB | |
| spliceai | 31,683,675 | 238 MB | `assume_unique` exhaustively verified |
| dbnsfp | 836,391 | 41 MB | |
| cadd | 121,809,062 | 1040 MB | `assume_unique` exhaustively verified over all 121.8M keys |

CADD chr21 post-fix: 124.8 s total (16.1 s read+normalize+dedup, 105.1 s
tier-join+write including the sorted retry), peak RSS 17.0 GB. Reservation
peaks: `ExternalSorter` 7.29 GB (real, spillable work), `HashJoinInput` 0.70 GB.

Full-genome runs have not been attempted. `notebooks/build_plugin_caches.ipynb`
in the vepyr branch drives the whole matrix; note it still downloads sources
rather than using the `rclone serve http` trick in §2.1, which is worth porting.

---

## 8. 2026-08-25 follow-up: remove the false HashJoin input and final merge

Two independent costs remained in the chr1 partition matrix:

1. The probe filter used an inner join against
   `SELECT DISTINCT start, allele_string FROM plugin`. The distinct aggregate
   emitted shared-buffer batches immediately below `HashJoinExec`, so its exact
   0.25 GiB AlphaMissense estimate became a 7.9 GiB runtime reservation and an
   inevitable Hash-to-SortMerge retry.
2. The final `(tier, start)` shard was encoded as warm and cold temporary
   Parquet files, then both files were decoded and encoded again into the final
   shard. That serial tail cost about 129 seconds for SpliceAI and 142 seconds
   for CADD at `target_partitions=4`.

The production path now registers a narrow `(start, allele_string)` plugin-key
MemTable and filters variation with a true `LEFT SEMI JOIN`. A semi join cannot
fan out variation rows when plugin keys repeat, so the distinct aggregate is
unnecessary. DataFusion still selects Hash or SortMerge from actual plan
statistics; no plugin, chromosome, partition-count, or join-strategy threshold
was added.

The final query now uses `ORDER BY tier, start` and streams directly into one
atomic `.build.tmp` shard. The storage-boundary check validates every emitted
key lexicographically before the batch is written. This removes the two tier
filters, two intermediate encodes, two intermediate decodes, and final
re-encode.

Release + `target-cpu=native`, chr1, `target_partitions=4`:

| plugin | before | after | wall delta | RSS before → after |
|---|---:|---:|---:|---:|
| ClinVar | 4.26 s | 4.08 s | -4.2% | 1.98 → 1.95 GiB |
| AlphaMissense | 15.35 s | 8.08 s | -47.4% | 5.71 → 3.01 GiB |
| dbNSFP | 75.29 s | 57.84 s | -23.2% | 11.83 → 11.81 GiB |
| SpliceAI | 357.38 s | 230.27 s | -35.6% | 31.75 → 27.09 GiB |
| CADD | 382.34 s | 221.89 s | -42.0% | 29.12 → 28.90 GiB |

AlphaMissense's probe HashJoin now reserves 0.23 GiB and completes on the first
plan. dbNSFP likewise completes its first HashJoin at 0.39 GiB. SpliceAI and
CADD report missing projected byte statistics and select SortMerge before
execution, so neither wastes time on a doomed HashJoin.

All five new shards match the previous p4 shards on row count plus both XOR and
sum aggregates of DuckDB whole-row hashes. Warm/cold manifest counts also match
exactly. Logs and outputs are under
`/Users/mwiewior/workspace/data_vepyr/plugin_partition_bench_direct_p4`.
