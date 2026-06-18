# Cold-Tier Point Store — Design Note (investigation seed)

> Not a plan yet. Captures the structural problem + the gate measurement that decides
> whether to build a non-columnar cold store. Companion to the size-optimization work
> (encoding tuning) and the (rejected/modest) late-materialization plan.

## Problem (measured, chr1 HG002, single-path lance, true 1-CPU)

The cold tier is read with massive **columnar miniblock amplification**:
- cold rows taken: **66,461** (9% of the 735,017 taken); they sit at ~46,785 distinct positions.
- distinct cold miniblocks decompressed: **11,147 of ~20,632 (~54%)** → ~6 useful rows per
  4096-row block (~0.15% useful).
- cold decode ≈ **8.35s** of the ~12s pure decode (disk-cold e2e `cold_tier_load` ~10.1s
  post-CacheOnly); **266 µs/row vs 8.2 µs/row hot** (32× worse per row).
- cold content is ~49% gnomADg-family frequencies; cold rows are sparse/null-heavy.

This amplification is **inherent to columnar + sparse, scattered point access**, not a tuning bug.

## Why a re-sort / re-tier does NOT fix it (ruled out)

- The cold query set is **unpredictable and genome-wide**, so any position-based sort still
  scatters a query across most miniblocks.
- "Cluster cold rows near hot variants": most cold reads *are* near-hot (86% within ±32 bp),
  but near-hot cold rows are ~24.7M (29% of cold) — still position-scattered → still thousands
  of miniblocks. No win.
- "Promote near-hot cold into hot": measured **backfire** (hot blocks are 5× costlier each —
  full common-variant frequency payloads — and the workload already reads ~all 871 of them;
  promotion only adds expensive hot blocks). Sweet spot was ±1–2 bp ≈ −4%, reversing past that.
- "Smaller cold miniblocks": measured **strictly worse** (per-block overhead + lost zstd context;
  256–2048 monotonically slower + 2.4× file bloat; 4096 is the floor, ≥8192 flat).

## The structural fix: a point-lookup cold store

Replace the cold tier's columnar miniblock layout with a **row-major / KV point store** keyed by
the lance row identity (or position), so a lookup fetches **exactly the requested rows** with no
4096-row block decompress. Hot tier stays dense Lance columnar (it's efficient — 8.2 µs/row,
reads ~all its blocks anyway). This is the *correct* shape for sparse point access and is what
the original two-tier design (warm-fjall + cold) was reaching for; the build config still lists
`heed_payload_zstd`, and the sandbox has `payload_sidecar.rs` / `row_sidecar.rs` prototypes.

## GO/NO-GO gate (do this before building anything)

Measure, on the **real 66,461 cold row_ids** (already dumped in `/tmp/chr1_rowids.txt`; filter to
tier==1 via the global-index → tier map used in the float/late-mat analyses):
1. **Columnar baseline:** time `ds.take(cold_row_ids, [read columns])` at true 1-CPU (≈ the cold
   slice of the current decode).
2. **Point-store candidate:** load the same rows' payloads from a heed/fjall (or simple row-major
   mmap) store keyed by row_id; time the point-gets at 1-CPU.
3. **GO** if the point store is materially faster (target: collapse the 32× cold/hot per-row gap —
   i.e. cold approaching hot's per-row cost). **NO-GO** if KV per-get overhead ≈ the columnar
   amplification it removes.

The sandbox already builds heed/fjall payload sidecars (`payload_sidecar`, `PayloadBench`) — reuse
that harness; no production code needed for the gate.

## Trade-offs / risks (if GO)

- **Two storage formats** (Lance hot + KV cold) → a new cold read path + a build step that writes
  the cold payload store, and **row-id stability** between the lance dataset and the sidecar
  (rebuild in lockstep, same as `.posidx`/`.varbf`).
- **Parity re-verification** (chr1+chr4 100% CSQ) on the new path.
- **Compaction/format maintenance** of the KV store.
- Composes with: **encoding/size tuning** (orthogonal — shrinks whatever stays columnar) and
  makes **late-materialization** moot for cold (point store already fetches exact rows).

## Recommended sequencing
1. Ship **encoding/size tuning** first (safest, low-risk IO win) — active work.
2. Run this **gate measurement** (cheap, sandbox-only) to decide if the cold point store is worth
   the architectural change.
3. Only if GO: brainstorm → spec → plan the two-format cold store.
