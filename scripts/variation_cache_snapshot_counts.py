#!/usr/bin/env python3
"""Snapshot per-contig row counts of every dataset in a variation entity dir."""
import json
import sys
from pathlib import Path
import lance


def main(entity_dir: str, out_path: str) -> None:
    root = Path(entity_dir)
    counts = {}
    for ds_dir in sorted(root.glob("*.lance")):
        name = ds_dir.name[: -len(".lance")]
        counts[name] = lance.dataset(str(ds_dir)).count_rows()
    Path(out_path).write_text(json.dumps(counts, indent=2, sort_keys=True))
    total = sum(counts.values())
    print(f"{len(counts)} contigs, total rows = {total}")


if __name__ == "__main__":
    if len(sys.argv) != 3:
        print(
            "usage: variation_cache_snapshot_counts.py <entity_dir> <out.json>",
            file=sys.stderr,
        )
        sys.exit(64)
    main(sys.argv[1], sys.argv[2])
