#!/usr/bin/env python3
"""Classify a variation Lance dataset's schema as new (bundled chr1 layout),
old (exploded population AF columns), or missing."""
import sys
import lance


def kind(path: str) -> str:
    try:
        ds = lance.dataset(path)
    except Exception:
        return "missing"
    names = {f.name for f in ds.schema}
    if "tier" in names and "af_global" in names:
        return "new"
    return "old"


if __name__ == "__main__":
    if len(sys.argv) != 2:
        print("usage: variation_cache_schema_kind.py <dataset.lance>", file=sys.stderr)
        sys.exit(64)
    k = kind(sys.argv[1])
    print(k)
    sys.exit(0 if k in ("new", "old") else 2)
