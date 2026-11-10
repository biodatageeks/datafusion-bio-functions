#!/usr/bin/env python3
"""Verdict of a local gate: parity, identity and performance rules over arm directories.

An arm directory holds `contigs.tsv`, one row per chromosome, written by `run_arm.sh`.

    report.py --baseline DIR [DIR ...] --candidate DIR [DIR ...]
              [--parity DIR ...] [--identity DIR ...]
              [--baseline-1w DIR ...] [--candidate-1w DIR ...]
              --target-pct HG002=<p>,ClinVar=<p>

Exit status: 0 when every rule holds, 1 when one fails, 2 on unusable input.

Rules (per dataset where a dataset is named):
  complete              no arm is empty, every row finished with exit 0 and a numeric
                        time, all arms of a dataset cover the same chromosomes and every
                        group of arms covers every dataset
  quiet                 every timing run started with background CPU at or below 200%
  parity                every row of every --parity arm has verdict PASS
  identity              one raw body MD5 per chromosome across the candidate, parity and
                        identity arms and the first baseline arm
  target <d>            mean candidate total <= (1 - p/100) x mean baseline total; a
                        dataset without a target fails
  no contig regression  no chromosome slower by more than 10% and by more than 0.5 s
  peak rss <d>          largest candidate RSS <= 1.10 x largest baseline RSS
  1 worker <d>          candidate total <= 1.03 x baseline total
"""

import argparse
import csv
import os
import sys

COLUMNS = [
    "dataset",
    "contig",
    "workers",
    "mode",
    "compression",
    "exit",
    "wall_s",
    "rss_gib",
    "time_s",
    "records",
    "body_md5",
    "verdict",
    "env",
    "bg_cpu",
]
MAX_BG_CPU = 200.0
CONTIG_REGRESSION_RATIO = 1.10
CONTIG_REGRESSION_SECONDS = 0.5
RSS_RATIO = 1.10
ONE_WORKER_RATIO = 1.03


def load_arm(path):
    """Rows of `<path>/contigs.tsv` as dicts of strings."""
    table = os.path.join(path, "contigs.tsv")
    with open(table, newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        missing = [c for c in COLUMNS if c not in (reader.fieldnames or [])]
        if missing:
            raise ValueError(f"{table}: missing columns {missing}")
        # A chromosome that failed and was run again has two rows: the last one counts.
        latest = {}
        for row in reader:
            latest[(row["dataset"], row["contig"])] = row
        return list(latest.values())


def _number(value):
    try:
        return float(value)
    except (TypeError, ValueError):
        return None


def totals(rows):
    """Summed `time_s` per dataset; rows without a numeric time are left out."""
    out = {}
    for row in rows:
        seconds = _number(row["time_s"])
        if seconds is not None:
            out[row["dataset"]] = out.get(row["dataset"], 0.0) + seconds
    return out


def _datasets(arms):
    seen = []
    for arm in arms:
        for row in arm:
            if row["dataset"] not in seen:
                seen.append(row["dataset"])
    return seen


def _of_dataset(arms, dataset):
    """The arms that hold `dataset`, each reduced to that dataset's rows."""
    out = []
    for arm in arms:
        rows = [row for row in arm if row["dataset"] == dataset]
        if rows:
            out.append(rows)
    return out


def _mean(values):
    return sum(values) / len(values)


def _mean_total(arms, dataset):
    sums = [totals(arm).get(dataset, 0.0) for arm in _of_dataset(arms, dataset)]
    return _mean(sums) if sums else None


def _contig_means(arms, dataset):
    by_contig = {}
    for arm in _of_dataset(arms, dataset):
        for row in arm:
            seconds = _number(row["time_s"])
            if seconds is not None:
                by_contig.setdefault(row["contig"], []).append(seconds)
    return {contig: _mean(values) for contig, values in by_contig.items()}


def _peak_rss(arms, dataset):
    values = [_number(row["rss_gib"]) for arm in _of_dataset(arms, dataset) for row in arm]
    values = [v for v in values if v is not None]
    return max(values) if values else None


def _complete(groups):
    """Problems that make a comparison meaningless. `groups` maps a group name to its arms."""
    problems = []
    arms = [arm for group in groups.values() for arm in group]
    for name, group in groups.items():
        for arm in group:
            if not arm:
                problems.append(f"{name}: an arm has no rows")
    for arm in arms:
        for row in arm:
            if row["exit"] != "0":
                problems.append(f"{row['dataset']} {row['contig']}: exit {row['exit']}")
            elif _number(row["time_s"]) is None:
                problems.append(f"{row['dataset']} {row['contig']}: no annotation time")
    datasets = _datasets(arms)
    for dataset in datasets:
        sets = [frozenset(row["contig"] for row in arm) for arm in _of_dataset(arms, dataset)]
        if len(set(sets)) > 1:
            union = frozenset().union(*sets)
            short = min(sets, key=len)
            problems.append(f"{dataset}: arms cover different chromosomes (e.g. missing {sorted(union - short)})")
    for name, group in groups.items():
        missing = [d for d in datasets if group and d not in _datasets(group)]
        if missing:
            problems.append(f"{name}: no arm for {', '.join(missing)}")
    return problems


def _busy(arms):
    problems = []
    for arm in arms:
        for row in arm:
            load = _number(row.get("bg_cpu"))
            if load is None:
                problems.append(f"{row['dataset']} {row['contig']}: background CPU not recorded")
            elif load > MAX_BG_CPU:
                problems.append(f"{row['dataset']} {row['contig']}: background CPU {load:.0f}%")
    return problems


def _identity(arms):
    digests = {}
    for arm in arms:
        for row in arm:
            digests.setdefault((row["dataset"], row["contig"]), set()).add(row["body_md5"])
    problems = []
    for (dataset, contig), values in sorted(digests.items()):
        if len(values) != 1 or not next(iter(values)) or next(iter(values)) == "NA":
            problems.append(f"{dataset} {contig}: {len(values)} distinct body MD5")
    return problems, len(digests)


def check(baseline, candidate, parity, identity, baseline_1w, candidate_1w, targets):
    """Apply every rule. Returns `[(rule name, passed, detail), ...]`.

    `baseline`, `candidate`, `parity` and `identity` are lists of arms (an arm is a list
    of rows). `baseline_1w` and `candidate_1w` are lists of arms or None. `targets` maps a
    dataset name to the required improvement in percent.
    """
    baseline_1w = baseline_1w or []
    candidate_1w = candidate_1w or []
    results = []

    problems = _complete(
        {
            "baseline": baseline,
            "candidate": candidate,
            "parity": parity,
            "identity": identity,
            "baseline 1 worker": baseline_1w,
            "candidate 1 worker": candidate_1w,
        }
    )
    results.append(("complete", not problems, "; ".join(problems[:5]) or "all rows finished"))

    problems = _busy(baseline + candidate + baseline_1w + candidate_1w)
    results.append(("quiet", not problems, "; ".join(problems[:5]) or "every timing run started on a quiet machine"))

    if parity:
        failed = [
            f"{row['dataset']} {row['contig']}: {row['verdict']}"
            for arm in parity
            for row in arm
            if row["verdict"] != "PASS"
        ]
        count = sum(len(arm) for arm in parity)
        results.append(("parity", not failed, "; ".join(failed[:5]) or f"{count} chromosomes PASS"))

    reference = []
    for dataset in _datasets(baseline):
        reference.extend(_of_dataset(baseline, dataset)[:1])
    problems, count = _identity(candidate + parity + identity + candidate_1w + reference)
    results.append(("identity", not problems, "; ".join(problems[:5]) or f"{count} chromosomes, one body MD5 each"))

    for dataset in _datasets(baseline + candidate):
        base_total = _mean_total(baseline, dataset)
        cand_total = _mean_total(candidate, dataset)
        if base_total is None or cand_total is None or base_total <= 0:
            results.append((f"target {dataset}", False, "dataset missing from the baseline or the candidate"))
            continue
        change = 100.0 * (1.0 - cand_total / base_total)
        direction = "faster" if change >= 0 else "slower"
        detail = f"baseline {base_total:.1f} s, candidate {cand_total:.1f} s, {abs(change):.1f}% {direction}"
        if dataset in targets:
            limit = (1.0 - targets[dataset] / 100.0) * base_total
            results.append((f"target {dataset}", cand_total <= limit, f"{detail}; target {targets[dataset]:g}%"))
        else:
            results.append((f"target {dataset}", False, f"{detail}; no target given"))

        base_contigs = _contig_means(baseline, dataset)
        cand_contigs = _contig_means(candidate, dataset)
        slower = [
            f"{contig} {base_contigs[contig]:.2f} -> {seconds:.2f} s"
            for contig, seconds in cand_contigs.items()
            if contig in base_contigs
            and seconds > CONTIG_REGRESSION_RATIO * base_contigs[contig]
            and seconds > base_contigs[contig] + CONTIG_REGRESSION_SECONDS
        ]
        results.append((f"no contig regression {dataset}", not slower, "; ".join(slower[:5]) or "none"))

        base_rss = _peak_rss(baseline, dataset)
        cand_rss = _peak_rss(candidate, dataset)
        if base_rss is None or cand_rss is None:
            results.append((f"peak rss {dataset}", False, "no RSS recorded"))
        else:
            results.append(
                (
                    f"peak rss {dataset}",
                    cand_rss <= RSS_RATIO * base_rss,
                    f"baseline {base_rss:.2f} GiB, candidate {cand_rss:.2f} GiB",
                )
            )

    for dataset in _datasets(baseline_1w + candidate_1w):
        base_total = _mean_total(baseline_1w, dataset)
        cand_total = _mean_total(candidate_1w, dataset)
        if base_total is None or cand_total is None or base_total <= 0:
            results.append((f"1 worker {dataset}", False, "dataset missing from one side"))
            continue
        results.append(
            (
                f"1 worker {dataset}",
                cand_total <= ONE_WORKER_RATIO * base_total,
                f"baseline {base_total:.1f} s, candidate {cand_total:.1f} s",
            )
        )
    return results


def parse_targets(text):
    targets = {}
    for item in filter(None, (part.strip() for part in (text or "").split(","))):
        name, sep, value = item.partition("=")
        if not sep or not name:
            raise ValueError(f"--target-pct item '{item}' is not NAME=PERCENT")
        targets[name.strip().lower()] = float(value)
    return targets


def _print_tables(baseline, candidate):
    for dataset in _datasets(baseline + candidate):
        base = _contig_means(baseline, dataset)
        cand = _contig_means(candidate, dataset)
        print(f"\n{dataset}: annotation time per chromosome, seconds (mean over arms)")
        print(f"{'contig':>8s} {'baseline':>9s} {'candidate':>9s} {'ratio':>6s}")
        for contig in base:
            if contig in cand:
                ratio = cand[contig] / base[contig] if base[contig] else float("nan")
                print(f"{contig:>8s} {base[contig]:9.2f} {cand[contig]:9.2f} {ratio:6.2f}")
        base_total, cand_total = _mean_total(baseline, dataset), _mean_total(candidate, dataset)
        if base_total and cand_total:
            print(f"{'total':>8s} {base_total:9.1f} {cand_total:9.1f} {cand_total / base_total:6.2f}")


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--baseline", nargs="+", required=True, metavar="DIR")
    parser.add_argument("--candidate", nargs="+", required=True, metavar="DIR")
    parser.add_argument("--parity", nargs="*", default=[], metavar="DIR")
    parser.add_argument("--identity", nargs="*", default=[], metavar="DIR")
    parser.add_argument("--baseline-1w", nargs="*", default=[], metavar="DIR", dest="baseline_1w")
    parser.add_argument("--candidate-1w", nargs="*", default=[], metavar="DIR", dest="candidate_1w")
    parser.add_argument("--target-pct", default="", metavar="NAME=PCT[,NAME=PCT]", dest="target_pct")
    try:
        args = parser.parse_args(argv)
    except SystemExit as stop:
        return 2 if stop.code else 0
    try:
        targets = parse_targets(args.target_pct)
        groups = {
            name: [load_arm(path) for path in getattr(args, name)]
            for name in ("baseline", "candidate", "parity", "identity", "baseline_1w", "candidate_1w")
        }
    except (OSError, ValueError) as error:
        print(f"error: {error}", file=sys.stderr)
        return 2
    unknown = sorted(set(targets) - set(_datasets(groups["baseline"] + groups["candidate"])))
    if unknown:
        print(f"error: --target-pct names no dataset of this comparison: {', '.join(unknown)}", file=sys.stderr)
        return 2

    _print_tables(groups["baseline"], groups["candidate"])
    results = check(
        groups["baseline"],
        groups["candidate"],
        groups["parity"],
        groups["identity"],
        groups["baseline_1w"],
        groups["candidate_1w"],
        targets,
    )
    print()
    for name, passed, detail in results:
        print(f"{'PASS' if passed else 'FAIL'}  {name}: {detail}")
    passed = all(ok for _, ok, _ in results)
    print(f"VERDICT: {'PASS' if passed else 'FAIL'}")
    return 0 if passed else 1


if __name__ == "__main__":
    sys.exit(main())
