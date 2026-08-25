# Plugin-cache build validation — testing handover (2026-08-25)

How the five plugin caches get built from their raw sources and checked, what
has actually been run, and what is left. The engine-side memory investigation
that this testing uncovered lives in
`2026-08-24-cadd-tier-join-reservation-leak-RESOLVED.md`; this doc is the
*testing* half and does not repeat it.

---

## 1. Why this exists

`datafusion-bio-functions` PR #217 rewrote the plugin-cache build path
(streaming shard write, `assume_unique`, VCF/BED providers). Unit tests cover
the invariants in isolation; they cannot tell you whether a real 120M-row
chromosome builds, or whether the shard it writes is correct. Two reviews on
#217 asked for a real-chromosome check before merge. This is that check.

It is also the only exercise that touches several of the PR's changes on real
data at all — the VCF provider, `assume_unique`, the ordering guard and its
retry, and the streaming write at chromosome scale.

## 2. State

| scope | status |
|---|---|
| chr21, all 5 plugins | **done** — build + shard invariants pass |
| chr1, all 5 plugins | **done** — at `target_partitions` 4 and 8, shard equivalence verified (see the leak doc §8) |
| chr2–22, X | **not started** |
| Golden-VEP CSQ parity | **not started** — see §6, this is the real acceptance gate |

chr21 shard invariants (this session):

| plugin | rows | shard | `(tier,start)` | `assume_unique` |
|---|---:|---:|---|---|
| clinvar | 49,937 | 1 MB | OK | — |
| alphamissense | 698,535 | 6 MB | OK | — |
| spliceai | 31,683,675 | 238 MB | OK | exhaustive, 31.7M keys |
| dbnsfp | 836,391 | 41 MB | OK | — |
| cadd | 121,809,062 | 1040 MB | OK | exhaustive, 121.8M keys |

## 3. Getting the sources without downloading 168 GB

All six source files are bgzip+tabix indexed, so a chromosome can be pulled
without fetching the file. `rclone serve http` makes Drive range-readable:

```bash
rclone serve http gdrive-mw: \
  --drive-root-folder-id 1ZT3g31I0LXepORF_dusy47AuXOnbRlZ_ \
  --addr localhost:18080 --read-only &
curl -s -o /dev/null -w '%{http_code}\n' http://localhost:18080/    # expect 200
tabix http://localhost:18080/cadd/whole_genome_SNVs.tsv.gz 21 > cadd_chr21_snv.tsv
```

CADD's chr21 slice comes out of an 87 GB remote file in ~50 s. The full source
set is ~168 GB against ~230 GB free, so **without this you cannot hold two large
sources at once** and the whole matrix becomes a disk-juggling exercise.

**Contig naming differs per source and a wrong prefix yields a silently empty
slice, not an error.** `chr`-prefixed for alphamissense; bare for clinvar, cadd,
dbnsfp, spliceai. Verify with `tabix -l <url> | head -3`.

## 4. Building and checking one chromosome

```bash
python scripts/debug/build_plugin_chrom.py <plugin> <chrom>
```

Env: `RUST_LOG=info` (**required** — all instrumentation logs at info),
`VEPYR_TRACE_MEMORY=1` (per-reservation peaks + batch buffer-sharing stats),
`VEPYR_EXPLAIN_TIER_JOIN=1` (physical plan),
`VEPYR_PLUGIN_RETRY_MEMORY_MIB=<n>`.

Every shard must satisfy, and these are what the checks assert:

1. **`start` ascending within each tier.** `PageDir::resolve_ranges` binary
   searches inside a tier, so an unsorted run silently misses rows that are
   present — a wrong-answer bug, not a crash.
2. **Tiers grouped** warm-then-cold.
3. **No duplicate probe keys** for an `assume_unique` plugin. A duplicate means
   the manifest flag is wrong and the runtime `HashMap` keeps the *last* row,
   inverting VEP's first-in-file rule. Check this **exhaustively** — the
   builder's own check samples only the first 2M distinct keys, so a build
   passing is not evidence the flag is right.

`scripts/debug/measure_probe_size.py` and
`measure_batch_buffer_sharing.py` are the two measurement tools that the memory
investigation needed; keep them for the next scale problem.

## 5. Running the matrix

`notebooks/build_plugin_caches.ipynb` (vepyr branch
`plugin-cache-build-validation`) drives all five, smallest first
(clinvar 192 MB → cadd 88 GB), downloading each source, building, then deleting
it. **It predates the `rclone serve http` approach in §3 and still downloads —
porting that in is the single biggest improvement available to it.**

Order matters and is not arbitrary: clinvar is 192 MB and exercises the VCF
provider, so a broken pipeline surfaces in minutes rather than after an 88 GB
download.

Two mistakes from this session worth not repeating:

- **Do not pipe per-plugin output through `tail -N`.** I did, and it truncated
  the retry/slice diagnostics, which cost a wrong conclusion ("the retry never
  fired") that took a rerun to correct.
- **The peak-reservation table does not print when a build fails inside
  `tiered_stream_sorted`** rather than the write — the report call sits after
  it. Worth moving if you are debugging a failure there.

## 6. The gate that has not been run

Everything above checks that a shard is *internally* well-formed. None of it
checks that the annotations are *right*. The real acceptance gate is the
golden-VEP CSQ comparison: build the caches, annotate HG002 with all five
plugins attached, and compare per-field against an Ensembl VEP run with the
same plugins.

That harness already exists — `run_comparison.py --profile merged_plugins` on
vepyr `plugins-fixes` (PR #48), which adds the profile and the representation-
only field equality the plugin fields need. **`merged_plugins` has never been
run.** Until it is, "the plugin cache builds" is the only claim supported.

Note `verify_parity_gate.py` deliberately refuses plugin profiles: it pins the
Ensembl core CSQ contract, and a plugin run emits fields outside it. Extending
the gate to plugin fields is unstarted work.

## 7. Open items

1. **Golden-VEP parity for `merged_plugins`** (§6) — the actual gate.
2. **chr2–22, X.** chr1 and chr21 are the extremes of size; the middle is
   presumably uneventful but unverified.
3. **Port `rclone serve http` into the notebook** (§3).
4. **Plugin-field parity gate** — `verify_parity_gate.py` cannot cover plugin
   profiles today.
5. **BED provider has no real-data coverage.** No BED manifest exists in
   vepyr-plugins, so `ProviderKind::Bed` is exercised only by the schema-pin
   unit test added in `84e5cfe`.

## 8. Related

| what | where |
|---|---|
| Engine + memory investigation | `datafusion-bio-functions` PR #217, branch `plugins-fixes` |
| Multi-part `source_path`, notebook | `vepyr` PR #49, branch `plugin-cache-build-validation` |
| Comparator + `merged_plugins` profile | `vepyr` PR #48, branch `plugins-fixes` |
| CADD as two `[[source]]` parts | `vepyr-plugins` PR #4 |
| Memory/leak investigation | `2026-08-24-cadd-tier-join-reservation-leak-RESOLVED.md` |
