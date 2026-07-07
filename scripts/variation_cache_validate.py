#!/usr/bin/env python3
"""Validate the regenerated variation cache: uniform schema, single LazyFrame,
row-count parity vs the pre-flight snapshot."""
import json
import sys
from pathlib import Path
import lance
import polars as pl
from variation_lazyframe import scan_all_variation

TARGET = {
    "chrom", "start", "end", "allele_string", "failed", "variation_name",
    "clin_sig", "clin_sig_allele", "clinical_impact", "phenotype_or_disease",
    "pubmed", "somatic", "minor_allele", "minor_allele_freq", "clinvar_ids",
    "cosmic_ids", "dbsnp_ids", "tier", "af_global", "af_gnomade", "af_gnomadg",
}


def main(entity_dir: str, before_json: str) -> None:
    root = Path(entity_dir)
    before = json.loads(Path(before_json).read_text())
    errors = []

    # 1. schema uniformity
    per_contig_rows = {}
    for ds_dir in sorted(root.glob("*.lance")):
        name = ds_dir.name[: -len(".lance")]
        ds = lance.dataset(str(ds_dir))
        names = {f.name for f in ds.schema}
        if names != TARGET:
            errors.append(f"{name}: schema {sorted(names ^ TARGET)} differs from target")
        per_contig_rows[name] = ds.count_rows()

    # 2. row-count parity vs snapshot
    for name, n_before in before.items():
        n_after = per_contig_rows.get(name)
        if n_after != n_before:
            errors.append(f"{name}: rows {n_before} -> {n_after} (mismatch)")

    # 3. single LazyFrame builds and total matches sum of per-dataset counts
    total_lf = scan_all_variation(entity_dir).select(pl.len()).collect().item()
    total_sum = sum(per_contig_rows.values())
    if total_lf != total_sum:
        errors.append(f"LazyFrame total {total_lf} != sum-of-counts {total_sum}")

    if errors:
        print("FAIL:")
        for e in errors:
            print("  -", e)
        sys.exit(1)
    print(f"OK: {len(per_contig_rows)} contigs, uniform schema, total rows = {total_sum}")


if __name__ == "__main__":
    main(sys.argv[1], sys.argv[2])
