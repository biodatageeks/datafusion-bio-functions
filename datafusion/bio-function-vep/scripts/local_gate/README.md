# Local parity and performance gate

Checks a change to `datafusion-bio-function-vep` on one machine before a vepyr pin PR:
output parity with Ensembl VEP, output identity across worker counts, and run time, on
HG002 and ClinVar with the release 116 merged cache.

The scripts build vepyr against a functions checkout through a **path dependency** and
call vepyr's own `e2e-testing/scripts/run_comparison.py` once per chromosome. A path
build never qualifies a release: the `/parity-wgs` run on the vepyr pin PR remains the
release gate.

## Requirements

- A vepyr checkout next to this repository (`../vepyr`), or `VEPYR_REPO`.
- Data under `$HOME/workspace/data_vepyr`, or `DATA_VEPYR_DIR`: the merged cache
  (`cache/116_GRCh38_merged`), `input/HG002_normalized.vcf.gz`,
  `input/clinvar/clinvar_20261011.norm.vcf.gz`, the reference FASTA and the two VEP
  goldens in `output/116/`. Docker is not needed.
- `uv`, Rust 1.91.1, `tabix`, `bgzip`.
- 20 GiB of free disk during every run. Plain ClinVar output is up to 7 GB per
  chromosome, so each chromosome's VCFs are deleted before the next one starts.
- An otherwise idle machine: a run waits up to three minutes for background CPU to drop
  below 200%.

State goes under `GATE_ROOT` (default `$HOME/research/benchmarks/local-gate`); nothing
is written inside the vepyr checkout except one registered worktree, `vepyr-build`.

## Full gate

```bash
git worktree add --detach "$HOME/research/benchmarks/local-gate/functions-baseline" origin/master
datafusion/bio-function-vep/scripts/local_gate/gate.sh \
  "$HOME/research/benchmarks/local-gate/functions-baseline" "$PWD" \
  "$HOME/research/benchmarks/local-gate/my-change" \
  --target-pct HG002=5,ClinVar=15 --bgzf-identity
```

It builds both sides, then for HG002 and for ClinVar runs a discarded warm-up, the
8-worker timing arms in the order baseline, candidate, candidate, baseline, the
candidate's parity run, the candidate at 4 workers and at 1 worker, the baseline at
1 worker and, with `--bgzf-identity`, the candidate with BGZF output. It ends with the
report in `<gate-dir>/verdict.txt` and exits 0 only if every rule passes. Allow about
three hours: besides annotation, the runner verifies its input on every call, about
12 s per ClinVar chromosome. The command is resumable: repeat it after an interruption.

`--target-pct` is the required improvement of the 8-worker total in percent. A negative
value allows that much slowdown: `HG002=-3` accepts up to 3% slower.
`--candidate-env VEP_X=1 ...` sets switches for the candidate arms only.

## Rules

| Rule | Meaning |
|---|---|
| `complete` | every chromosome finished and has a time; all arms cover the same chromosomes |
| `parity` | the candidate's body md5 matches the VEP golden on every chromosome (`strict` and `canonical` for HG002, `strict-clin-sig` for ClinVar) |
| `identity` | one raw body MD5 per chromosome across the candidate at 1, 4 and 8 workers, its BGZF output and the baseline |
| `target <dataset>` | mean candidate total at most `(1 - p/100)` of the mean baseline total |
| `no contig regression <dataset>` | no chromosome slower by more than 10% and more than 0.5 s |
| `peak rss <dataset>` | largest candidate RSS at most 10% above the largest baseline RSS |
| `1 worker <dataset>` | candidate total at most 3% above the baseline total |

Times are the runner's `annotation.time_s` (the `vepyr.annotate` call), which it rounds
to 0.1 s. RSS and wall time are those of the whole runner process.

## Single pieces

```bash
L=datafusion/bio-function-vep/scripts/local_gate; G="$HOME/research/benchmarks/local-gate"
$L/build_venv.sh "$PWD" "$G/venv-candidate"                      # build one side
$L/run_arm.sh "$G/venv-candidate" "$G/try/parity" hg002 8 parity --chroms chr22
$L/run_arm.sh "$G/venv-candidate" "$G/try/off" hg002 8 time --chroms chr22 VEP_SOME_SWITCH=0
python3 $L/report.py --baseline "$G/try/off" --candidate "$G/try/parity" --target-pct HG002=-50
python3 $L/test_report.py                                        # unit tests of the rules
```

Each arm directory holds `contigs.tsv` (one row per chromosome: exit code, wall time,
RSS, annotation time, records, raw body MD5, verdict, switches), `arm.json` (build
identity and arguments) and `logs/`.

`GATE_HG002_CHROMS` and `GATE_CLINVAR_CHROMS` restrict `gate.sh` to a few chromosomes to
rehearse the driver; a verdict from a restricted run is not a gate result.
