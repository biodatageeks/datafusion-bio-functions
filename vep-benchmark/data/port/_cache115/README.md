# _cache115 — Tier 1 v115 fixture

Tier 1 fast-cache for the Perl→Rust VEP test port. Built from upstream
Ensembl VEP v115 (`homo_sapiens`, `GRCh38`). Used by every `port_*.rs`
that declares `Cache tier: fast` in its `porting-tests/detailed_plans/`.

**Spec**: `porting-tests/testing.md` §1 (Tier 1).
**Region design**: see [`REGION_DESIGN.md`](REGION_DESIGN.md).

## Contents

| Path | Content |
|---|---|
| `parquet/115_GRCh38_vep/{transcript,variation,regulatory,motif_feature}/{21,MT}.parquet` | Per-entity parquet, chr21 + MT slice. Used by the vepyr engine. |
| `parquet/115_GRCh38_vep/variation.fjall/` + `translation_sift.fjall/` | fjall KV stores (optional; vepyr can run parquet-only) |
| `native_cache/homo_sapiens/115_GRCh38/{21,MT,info.txt}` | Native Ensembl cache subset, chr21 + MT. Used by the Perl-VEP oracle inside the Docker container in `gen_golden.sh`. |
| `reference.fa` + `.fai` | FASTA slice (chr21 + MT) |
| `ORACLE_IMAGE` | Pinned `ensemblorg/ensembl-vep@sha256:<digest>` for `gen_golden.sh` |
| `REGION_DESIGN.md` | Why these regions, what's covered, B1-B3 outcomes |

## How to rebuild

Pre-reqs (per `scripts/port/build_v115_parquet.sh`):
1. Whole-genome v115 vepyr cache at `$VEPYR_WHOLE_GENOME_PARQUET`
   (default: `/Users/wojtek/Documents/vepyr/_cache_v115/parquet/115_GRCh38_vep`).
   Build with one call: `vepyr.build_cache(release=115, cache_dir=...)`.
2. Native v115 Ensembl cache at
   `/Users/wojtek/Documents/vepyr/_cache_v115/homo_sapiens/115_GRCh38/`
   (auto-produced as a by-product of `vepyr.build_cache`).
3. Reference FASTA (primary assembly) at
   `/Users/wojtek/Documents/vepyr/_cache_v115/Homo_sapiens.GRCh38.dna.primary_assembly.fa[.fai]`.
4. `samtools` on `PATH`.

Then:

```bash
cd datafusion-bio-functions
./scripts/port/build_v115_parquet.sh
```

Idempotent: re-runs are no-ops if `parquet/` already exists.

## How to use from a port

```rust
mod port_common;

const HARD_FIELDS: &[&str] = &[
    "Consequence", "Feature", "Gene", "BIOTYPE", "IMPACT",
    "STRAND", "SIFT", "PolyPhen", "NEAREST", "SYMBOL",
];

#[tokio::test(flavor = "multi_thread")]
async fn port_cache_transcript_csq_matches_golden() {
    let ran = port_common::run_and_compare_csq("cache_transcript", HARD_FIELDS)
        .await
        .unwrap();
    assert!(ran, "fixtures must be present (git lfs pull + build_v115_parquet.sh)");
}
```

Where `cache_transcript` matches the directory under
`datafusion/bio-function-vep/tests/data/port/cache_transcript/`.
The harness automatically wires the v115 parquet + reference paths to this
fixture (`port_common::run_and_compare_csq`).
