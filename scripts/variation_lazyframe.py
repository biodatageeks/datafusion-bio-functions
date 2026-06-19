#!/usr/bin/env python3
"""Read every per-contig variation Lance dataset as one polars LazyFrame."""
from pathlib import Path
import lance
import polars as pl


def scan_all_variation(entity_dir: str) -> pl.LazyFrame:
    frames = []
    for ds_dir in sorted(Path(entity_dir).glob("*.lance")):
        ds = lance.dataset(str(ds_dir))
        frames.append(pl.scan_pyarrow_dataset(ds))
    if not frames:
        raise SystemExit(f"no .lance datasets under {entity_dir}")
    return pl.concat(frames, how="vertical")


if __name__ == "__main__":
    import sys

    lf = scan_all_variation(sys.argv[1])
    print(lf.select(pl.len()).collect().item())
