
---

## FINAL RESULTS (2026-08-28) — engine `f81f442`

### chr1-22 regression gate, re-run after the defect D fix

| metric | required | actual |
|---|---:|---:|
| chromosomes `body=ok` | 22/22 | **22/22** |
| records | 4,096,123 | **4,096,123** |
| CSQ entries | 69,299,753 | **69,299,753** |
| mismatches, every category | 0 | **0** |

Zero value mismatches, zero order mismatches, zero CSQ count/order mismatches, zero
variants-on-one-side-only, on every chromosome.

### chrX / chrY

| | before | after |
|---|---|---|
| chrX | 148 value + 4 order, strict DIFF | **0 value**, 4 order (defect E), strict FAIL on 4 records |
| chrY | 216 value, strict DIFF | **0 value, 0 order, strict PASS — byte-identical** |

**All 368 field-value mismatches across chrX and chrY are closed.** The only residue is 4 records
of `DOMAINS` ordering on a single transcript (defect E), deliberately not fixed — see above.

### Engine commits, all pushed to `fix/chrxy-core-annotation-parity`

| commit | defect | rows |
|---|---|---:|
| `76a6008` | A — EXON range | 324 combined |
| `25313d8` | A — INTRON range | |
| `40e662c` | `SoTerm::tier()` | — |
| `cc92486` | B — transcript_ablation via the generic tier gate | 37 |
| `bd04392` | C — intron_variant at a 2 bp intron | 2 |
| `f81f442` | D — HGVSp across a dropped insertion flank | 1 |

978 unit tests, clippy clean, `cargo fmt` clean.
