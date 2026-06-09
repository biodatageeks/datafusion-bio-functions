#!/usr/bin/env python3
"""Build and benchmark Lance chr1 variation-cache layouts.

This is a focused research harness for the VEP `everything` variation-cache
projection. It compares the current Parquet artifacts with Lance 2.1 and 2.2
physical layouts for:

* warm-cache full scans over the hot `everything` projection
* cold-cache batched point lookups on `position_key`

Run with pylance 7.0.0:

    uv run --with pylance==7.0.0 --with pyarrow \
      research/bench_lance_variation_chr1.py \
      --cache-dir /Users/mwiewior/workspace/data_vepyr/115_GRCh38_merged \
      --chrom chr1 \
      --variants 2.1-unpacked,2.2-packed \
      --cold-fragment-rows 512,1024,2048,4096,8192,16384 \
      --cold-sample-size 2000 \
      --force-build
"""

from __future__ import annotations

import argparse
import json
import os
import platform
import shutil
import time
from collections import OrderedDict
from collections.abc import Iterable, Iterator, Sequence
from dataclasses import asdict, dataclass
from pathlib import Path
from typing import Any

import lance
import pyarrow as pa
import pyarrow.dataset as ds
import pyarrow.parquet as pq


GIB = 1024**3
DEFAULT_CACHE_DIR = Path("/Users/mwiewior/workspace/data_vepyr/115_GRCh38_merged")
DEFAULT_OUTPUT_DIR = Path("/Users/mwiewior/workspace/data_vepyr/lance_variation_chr1_bench")
DEFAULT_REPORT_DIR = Path("research/reports")

AF_1KG_COLUMNS = ["AF", "AFR", "AMR", "EAS", "EUR", "SAS"]
AF_GNOMADE_COLUMNS = [
    "gnomADe",
    "gnomADe_AFR",
    "gnomADe_AMR",
    "gnomADe_ASJ",
    "gnomADe_EAS",
    "gnomADe_FIN",
    "gnomADe_MID",
    "gnomADe_NFE",
    "gnomADe_REMAINING",
    "gnomADe_SAS",
]
AF_GNOMADG_COLUMNS = [
    "gnomADg",
    "gnomADg_AFR",
    "gnomADg_AMI",
    "gnomADg_AMR",
    "gnomADg_ASJ",
    "gnomADg_EAS",
    "gnomADg_FIN",
    "gnomADg_MID",
    "gnomADg_NFE",
    "gnomADg_REMAINING",
    "gnomADg_SAS",
]

ALWAYS_NULL_COLUMNS = [
    "minor_allele",
    "minor_allele_freq",
    "assembly_ids",
    "gencode_ids",
    "genebuild_ids",
    "gnomade_ids",
    "gnomadg_ids",
    "polyphen_ids",
    "refseq_ids",
    "regbuild_ids",
    "sift_ids",
    "src_1000genomes_ids",
]

# Current hot projection: Rust WARM_RUNTIME_COLUMNS plus the cold match key
# columns. This is intentionally not all fields in the source Parquet; source
# metadata and preserve-only xrefs are not read by the current everything path.
EVERYTHING_VARIATION_COLUMNS = [
    "position_key",
    "variant_keys",
    "chrom",
    "start",
    "end",
    "allele_string",
    "variation_name",
    "failed",
    "somatic",
    "strand",
    "minor_allele",
    "minor_allele_freq",
    "clin_sig",
    "phenotype_or_disease",
    "clinical_impact",
    "clin_sig_allele",
    *AF_1KG_COLUMNS,
    *AF_GNOMADG_COLUMNS,
    *AF_GNOMADE_COLUMNS,
    "clinvar_ids",
    "cosmic_ids",
    "dbsnp_ids",
    "pubmed",
]

PACKED_STRUCT_GROUPS: OrderedDict[str, list[str]] = OrderedDict(
    [
        ("match_payload", ["allele_string", "end", "failed"]),
        ("identity_text", ["variation_name", "dbsnp_ids"]),
        (
            "clinical_payload",
            ["clin_sig", "clin_sig_allele", "clinical_impact", "pubmed", "clinvar_ids", "cosmic_ids"],
        ),
        ("variant_flags", ["somatic", "phenotype_or_disease", "strand"]),
        ("af_1kg", AF_1KG_COLUMNS),
        ("af_gnomade", AF_GNOMADE_COLUMNS),
        ("af_gnomadg", AF_GNOMADG_COLUMNS),
    ]
)
CHILD_TO_STRUCT = {
    child: struct_name
    for struct_name, children in PACKED_STRUCT_GROUPS.items()
    for child in children
}

LANCE_MINIBLOCK_ZSTD3_METADATA = {
    b"lance-encoding:structural-encoding": b"miniblock",
    b"lance-encoding:compression": b"zstd",
    b"lance-encoding:compression-level": b"3",
    b"lance-encoding:dict-values-compression": b"zstd",
    b"lance-encoding:dict-values-compression-level": b"3",
    b"lance-encoding:rle-threshold": b"0.95",
    b"lance-encoding:dict-size-ratio": b"0.99",
    b"lance-encoding:dict-divisor": b"1",
    b"lance-encoding:minichunk-size": b"16384",
}
CONSTANT_NULL_METADATA = {
    b"lance-encoding:structural-encoding": b"constant",
    b"lance-encoding:constant": b"null",
}
PACKED_STRUCT_METADATA = {
    b"lance-encoding:packed": b"true",
    b"lance-encoding:structural-encoding": b"miniblock",
    b"lance-encoding:compression": b"zstd",
    b"lance-encoding:compression-level": b"3",
}


@dataclass
class BuildStats:
    variant: str
    tier: str
    path: str
    storage_version: str
    fragment_rows: int
    rows: int
    build_seconds: float
    index_seconds: float
    artifact_bytes: int
    built: bool


@dataclass
class BenchStats:
    variant: str
    tier: str
    operation: str
    fragment_rows: int | None
    rows: int
    seconds: float
    rows_per_s: float
    artifact_gib: float
    columns: int
    path: str
    notes: str = ""


def storage_version_for_variant(variant: str) -> str:
    if variant == "2.1-unpacked":
        return "2.1"
    if variant == "2.2-packed":
        return "2.2"
    raise ValueError(f"unknown variant {variant!r}")


def validate_variant(variant: str) -> None:
    storage_version_for_variant(variant)


def with_metadata(field: pa.Field, metadata: dict[bytes, bytes]) -> pa.Field:
    merged = dict(field.metadata or {})
    merged.update(metadata)
    return field.with_metadata(merged)


def source_field_or_default(source_schema: pa.Schema, name: str) -> pa.Field:
    if name in source_schema.names:
        return source_schema.field(name)
    if name == "position_key":
        return pa.field(name, pa.uint64(), nullable=False)
    if name in {"start", "end"}:
        return pa.field(name, pa.int64())
    if name == "failed":
        return pa.field(name, pa.bool_())
    if name == "variant_keys":
        return pa.field(name, pa.list_(pa.string()))
    return pa.field(name, pa.string())


def encoded_field(source_schema: pa.Schema, name: str) -> pa.Field:
    field = source_field_or_default(source_schema, name)
    if name in ALWAYS_NULL_COLUMNS:
        return with_metadata(field.with_nullable(True), CONSTANT_NULL_METADATA)
    return with_metadata(field, LANCE_MINIBLOCK_ZSTD3_METADATA)


def schema_for_variant(
    variant: str,
    source_schema: pa.Schema,
    logical_columns: Sequence[str],
) -> pa.Schema:
    """Return the physical Lance schema for a logical projection."""

    validate_variant(variant)
    logical_set = set(logical_columns)

    if variant == "2.1-unpacked":
        return pa.schema([encoded_field(source_schema, name) for name in logical_columns])

    fields: list[pa.Field] = []
    emitted_structs: set[str] = set()
    for name in logical_columns:
        struct_name = CHILD_TO_STRUCT.get(name)
        if struct_name is None:
            fields.append(encoded_field(source_schema, name))
            continue
        if struct_name in emitted_structs:
            continue
        children = [child for child in PACKED_STRUCT_GROUPS[struct_name] if child in logical_set]
        child_fields = [source_field_or_default(source_schema, child).with_nullable(True) for child in children]
        fields.append(pa.field(struct_name, pa.struct(child_fields)).with_metadata(PACKED_STRUCT_METADATA))
        emitted_structs.add(struct_name)

    return pa.schema(fields)


def writer_schema_for_variant(
    variant: str,
    source_schema: pa.Schema,
    logical_columns: Sequence[str],
    unsafe_packed_metadata: bool = False,
) -> pa.Schema:
    """Return a pylance-7.0.0 writer-safe schema.

    The intended 2.2 schema marks struct fields with
    `lance-encoding:packed=true`. In pylance/lance 7.0.0 that metadata panics
    when any packed struct child is all-null in the written batch. The benchmark
    still writes the 2.2 physical struct layout, but strips only that unsafe
    metadata key by default so the run can complete reproducibly.
    """

    schema = schema_for_variant(variant, source_schema, logical_columns)
    if variant != "2.2-packed" or unsafe_packed_metadata:
        return schema

    fields: list[pa.Field] = []
    for field in schema:
        if field.name not in PACKED_STRUCT_GROUPS:
            fields.append(field)
            continue
        metadata = dict(field.metadata or {})
        metadata.pop(b"lance-encoding:packed", None)
        fields.append(field.with_metadata(metadata))
    return pa.schema(fields, metadata=schema.metadata)


def physical_projection_for_variant(variant: str, logical_columns: Sequence[str]) -> list[str]:
    """Map logical `everything` names to Lance physical projection names."""

    validate_variant(variant)
    if variant == "2.1-unpacked":
        return list(logical_columns)

    projection: list[str] = []
    emitted_structs: set[str] = set()
    for name in logical_columns:
        struct_name = CHILD_TO_STRUCT.get(name)
        if struct_name is None:
            projection.append(name)
        elif struct_name not in emitted_structs:
            projection.append(struct_name)
            emitted_structs.add(struct_name)
    return projection


def logical_columns_available(source_schema: pa.Schema, include_missing_variant_keys: bool) -> list[str]:
    available = set(source_schema.names)
    columns = []
    for name in EVERYTHING_VARIATION_COLUMNS:
        if name in available or (include_missing_variant_keys and name == "variant_keys"):
            columns.append(name)
    return columns


def table_batches(
    path: Path,
    columns: Sequence[str],
    batch_size: int,
    row_limit: int | None,
) -> Iterator[pa.RecordBatch]:
    dataset = ds.dataset(path, format="parquet")
    remaining = row_limit
    for batch in dataset.to_batches(columns=list(columns), batch_size=batch_size):
        if remaining is not None:
            if remaining <= 0:
                break
            if batch.num_rows > remaining:
                batch = batch.slice(0, remaining)
            remaining -= batch.num_rows
        yield batch


def array_for_field(
    batch_columns: dict[str, pa.Array],
    field: pa.Field,
    rows: int,
) -> pa.Array:
    if field.name in batch_columns:
        array = batch_columns[field.name]
        if array.type.equals(field.type):
            return array
        return array.cast(field.type)
    return pa.nulls(rows, type=field.type)


def struct_array_for_group(
    batch_columns: dict[str, pa.Array],
    source_schema: pa.Schema,
    struct_name: str,
    struct_field: pa.Field,
    rows: int,
) -> pa.StructArray:
    children: list[pa.Array] = []
    child_names = [child.name for child in struct_field.type]
    for child_name in child_names:
        child_field = source_field_or_default(source_schema, child_name).with_nullable(True)
        children.append(array_for_field(batch_columns, child_field, rows))
    return pa.StructArray.from_arrays(children, names=child_names)


def align_batch_to_schema(
    batch: pa.RecordBatch,
    source_schema: pa.Schema,
    lance_schema: pa.Schema,
    variant: str,
) -> pa.RecordBatch:
    rows = batch.num_rows
    batch_columns = {name: batch.column(i) for i, name in enumerate(batch.schema.names)}
    arrays: list[pa.Array] = []
    for field in lance_schema:
        if variant == "2.2-packed" and field.name in PACKED_STRUCT_GROUPS:
            arrays.append(struct_array_for_group(batch_columns, source_schema, field.name, field, rows))
        else:
            arrays.append(array_for_field(batch_columns, field, rows))
    return pa.RecordBatch.from_arrays(arrays, schema=lance_schema)


def aligned_batches(
    source_path: Path,
    source_schema: pa.Schema,
    lance_schema: pa.Schema,
    variant: str,
    logical_columns: Sequence[str],
    batch_size: int,
    row_limit: int | None,
) -> Iterator[pa.RecordBatch]:
    source_columns = [name for name in logical_columns if name in source_schema.names]
    for batch in table_batches(source_path, source_columns, batch_size, row_limit):
        yield align_batch_to_schema(batch, source_schema, lance_schema, variant)


def dir_bytes(path: Path) -> int:
    total = 0
    for root, _, files in os.walk(path):
        for name in files:
            total += (Path(root) / name).stat().st_size
    return total


def parquet_rows(path: Path) -> int:
    return pq.ParquetFile(path).metadata.num_rows


def row_count_for_limit(path: Path, row_limit: int | None) -> int:
    rows = parquet_rows(path)
    return min(rows, row_limit) if row_limit is not None else rows


def create_reader_from_batches(schema: pa.Schema, batches: Iterable[pa.RecordBatch]) -> pa.RecordBatchReader:
    return pa.RecordBatchReader.from_batches(schema, batches)


def build_lance_dataset(
    source_path: Path,
    lance_path: Path,
    variant: str,
    tier: str,
    fragment_rows: int,
    row_group_rows: int,
    batch_size: int,
    row_limit: int | None,
    force: bool,
    skip_build: bool,
    unsafe_packed_metadata: bool = False,
) -> BuildStats:
    if force and lance_path.exists():
        shutil.rmtree(lance_path)

    storage_version = storage_version_for_variant(variant)
    source_schema = pq.ParquetFile(source_path).schema_arrow
    logical_columns = logical_columns_available(source_schema, include_missing_variant_keys=True)
    lance_schema = writer_schema_for_variant(variant, source_schema, logical_columns, unsafe_packed_metadata)

    if lance_path.exists() or skip_build:
        if not lance_path.exists():
            raise FileNotFoundError(f"--skip-build requested but {lance_path} does not exist")
        dataset = lance.dataset(lance_path)
        index_started = time.perf_counter()
        if not dataset.has_index:
            dataset.create_scalar_index("position_key", "BTREE")
        index_seconds = time.perf_counter() - index_started
        return BuildStats(
            variant=variant,
            tier=tier,
            path=str(lance_path),
            storage_version=storage_version,
            fragment_rows=fragment_rows,
            rows=dataset.count_rows(),
            build_seconds=0.0,
            index_seconds=index_seconds,
            artifact_bytes=dir_bytes(lance_path),
            built=False,
        )

    lance_path.parent.mkdir(parents=True, exist_ok=True)
    batches = aligned_batches(
        source_path,
        source_schema,
        lance_schema,
        variant,
        logical_columns,
        batch_size,
        row_limit,
    )
    reader = create_reader_from_batches(lance_schema, batches)
    started = time.perf_counter()
    rows_per_file = fragment_rows if tier == "warm" else max(fragment_rows, 1_000_000)
    lance.write_dataset(
        reader,
        lance_path,
        schema=lance_schema,
        mode="create",
        data_storage_version=storage_version,
        max_rows_per_file=rows_per_file,
        max_rows_per_group=row_group_rows,
    )
    build_seconds = time.perf_counter() - started

    dataset = lance.dataset(lance_path)
    index_started = time.perf_counter()
    dataset.create_scalar_index("position_key", "BTREE")
    index_seconds = time.perf_counter() - index_started

    return BuildStats(
        variant=variant,
        tier=tier,
        path=str(lance_path),
        storage_version=storage_version,
        fragment_rows=fragment_rows,
        rows=dataset.count_rows(),
        build_seconds=build_seconds,
        index_seconds=index_seconds,
        artifact_bytes=dir_bytes(lance_path),
        built=True,
    )


def timed_table_read(callable_read: Any) -> tuple[float, int]:
    started = time.perf_counter()
    table = callable_read()
    seconds = time.perf_counter() - started
    return seconds, table.num_rows


def bench_parquet_full_scan(
    source_path: Path,
    logical_columns: Sequence[str],
    row_limit: int | None,
) -> BenchStats:
    columns = [name for name in logical_columns if name in pq.ParquetFile(source_path).schema_arrow.names]

    def read() -> pa.Table:
        dataset = ds.dataset(source_path, format="parquet")
        return dataset.head(row_limit, columns=columns) if row_limit is not None else dataset.to_table(columns=columns)

    seconds, rows = timed_table_read(read)
    return BenchStats(
        variant="parquet-current",
        tier="warm",
        operation="full_scan",
        fragment_rows=None,
        rows=rows,
        seconds=seconds,
        rows_per_s=rows / seconds if seconds > 0 else 0.0,
        artifact_gib=source_path.stat().st_size / GIB,
        columns=len(columns),
        path=str(source_path),
    )


def bench_lance_full_scan(
    lance_path: Path,
    variant: str,
    logical_columns: Sequence[str],
    row_limit: int | None,
) -> BenchStats:
    projection = physical_projection_for_variant(variant, logical_columns)

    def read() -> pa.Table:
        dataset = lance.dataset(lance_path)
        return dataset.to_table(columns=projection, limit=row_limit)

    seconds, rows = timed_table_read(read)
    return BenchStats(
        variant=variant,
        tier="warm",
        operation="full_scan",
        fragment_rows=None,
        rows=rows,
        seconds=seconds,
        rows_per_s=rows / seconds if seconds > 0 else 0.0,
        artifact_gib=dir_bytes(lance_path) / GIB,
        columns=len(projection),
        path=str(lance_path),
    )


def sample_cold_keys(source_path: Path, sample_size: int, row_limit: int | None) -> list[int]:
    row_count = row_count_for_limit(source_path, row_limit)
    stride = max(1, row_count // sample_size)
    keys: list[int] = []
    seen: set[int] = set()
    row_index = 0
    for batch in table_batches(source_path, ["position_key"], batch_size=262_144, row_limit=row_limit):
        column = batch.column(0)
        for value in column.to_pylist():
            if row_index % stride == 0 and value not in seen:
                seen.add(value)
                keys.append(value)
                if len(keys) >= sample_size:
                    return keys
            row_index += 1
    return keys


def lance_filter(keys: Sequence[int]) -> str:
    return "position_key IN (" + ",".join(str(key) for key in keys) + ")"


def bench_parquet_point_lookup(
    source_path: Path,
    logical_columns: Sequence[str],
    keys: Sequence[int],
    row_limit: int | None,
) -> BenchStats:
    source_schema = pq.ParquetFile(source_path).schema_arrow
    columns = [name for name in logical_columns if name in source_schema.names]
    key_array = pa.array(keys, type=pa.uint64())
    filter_expr = ds.field("position_key").isin(key_array)
    if row_limit is not None:
        max_key = max(keys) if keys else 0
        filter_expr = filter_expr & (ds.field("position_key") <= max_key)

    def read() -> pa.Table:
        dataset = ds.dataset(source_path, format="parquet")
        return dataset.to_table(
            columns=columns,
            filter=filter_expr,
            batch_readahead=16,
            fragment_readahead=8,
        )

    seconds, rows = timed_table_read(read)
    return BenchStats(
        variant="parquet-current",
        tier="cold",
        operation=f"point_lookup_{len(keys)}",
        fragment_rows=None,
        rows=rows,
        seconds=seconds,
        rows_per_s=rows / seconds if seconds > 0 else 0.0,
        artifact_gib=source_path.stat().st_size / GIB,
        columns=len(columns),
        path=str(source_path),
    )


def bench_lance_point_lookup(
    lance_path: Path,
    variant: str,
    fragment_rows: int,
    logical_columns: Sequence[str],
    keys: Sequence[int],
) -> BenchStats:
    projection = physical_projection_for_variant(variant, logical_columns)

    def read() -> pa.Table:
        dataset = lance.dataset(lance_path)
        return dataset.to_table(
            columns=projection,
            filter=lance_filter(keys),
            use_scalar_index=True,
            prefilter=True,
            late_materialization=True,
        )

    seconds, rows = timed_table_read(read)
    return BenchStats(
        variant=variant,
        tier="cold",
        operation=f"point_lookup_{len(keys)}",
        fragment_rows=fragment_rows,
        rows=rows,
        seconds=seconds,
        rows_per_s=rows / seconds if seconds > 0 else 0.0,
        artifact_gib=dir_bytes(lance_path) / GIB,
        columns=len(projection),
        path=str(lance_path),
    )


def result_markdown_table(rows: Sequence[dict[str, Any]]) -> str:
    lines = [
        "| variant | tier | operation | rows | seconds | rows/s | artifact GiB |",
        "|---|---|---:|---:|---:|---:|---:|",
    ]
    for row in rows:
        lines.append(
            f"| {row['variant']} | {row['tier']} | {row['operation']} | "
            f"{int(row['rows'])} | {float(row['seconds']):.3f} | "
            f"{float(row['rows_per_s']):.0f} | {float(row['artifact_gib']):.3f} |"
        )
    return "\n".join(lines)


def build_markdown_report(report: dict[str, Any]) -> str:
    config = report["config"]
    cold_rows = [
        row
        for row in report["results"]
        if row["tier"] == "cold" and row["variant"] != "parquet-current"
    ]
    lines = [
        f"# Lance Variation chr1 Benchmark",
        "",
        "## Summary",
        "",
        f"- Lance Python package: `{report['environment']['lance_version']}`",
        f"- PyArrow: `{report['environment']['pyarrow_version']}`",
        f"- Cache dir: `{config['cache_dir']}`",
        f"- Chromosome: `{config['chrom']}`",
        f"- Variants: `{', '.join(config['variants'])}`",
        f"- Cold sample size: `{config['cold_sample_size']}` position keys",
        f"- Cold fragment row range: `{', '.join(str(v) for v in config['cold_fragment_rows'])}`",
        f"- Row limit: `{config['row_limit'] if config['row_limit'] is not None else 'full chr1'}`",
        "",
        "## Results",
        "",
        result_markdown_table(report["results"]),
        "",
        "## Cold Fragment Detail",
        "",
        "| variant | row-group rows | lookup rows | seconds | rows/s | artifact GiB |",
        "|---|---:|---:|---:|---:|---:|",
    ]
    for row in cold_rows:
        lines.append(
            f"| {row['variant']} | {row['fragment_rows']} | {row['rows']} | "
            f"{row['seconds']:.3f} | {row['rows_per_s']:.0f} | {row['artifact_gib']:.3f} |"
        )

    parquet_cold = [row for row in report["results"] if row["tier"] == "cold" and row["variant"] == "parquet-current"]
    if parquet_cold:
        baseline = parquet_cold[0]
        lines.extend(
            [
                "",
                f"Parquet cold baseline: `{baseline['seconds']:.3f}s` for `{baseline['rows']}` rows "
                f"from `{config['cold_sample_size']}` sampled keys.",
                "",
            ]
        )
    else:
        lines.append("")

    lines.extend(
        [
        "## Build Artifacts",
        "",
        "| variant | tier | fragment rows | rows | build s | index s | artifact GiB | built |",
        "|---|---|---:|---:|---:|---:|---:|---|",
        ]
    )
    for row in report["builds"]:
        lines.append(
            f"| {row['variant']} | {row['tier']} | {row['fragment_rows']} | {row['rows']} | "
            f"{row['build_seconds']:.3f} | {row['index_seconds']:.3f} | "
            f"{row['artifact_bytes'] / GIB:.3f} | {row['built']} |"
        )
    lines.extend(
        [
            "",
            "## Layout Notes",
            "",
            "- `2.1-unpacked` keeps the hot projection as top-level Lance columns using storage version 2.1.",
            "- `2.2-packed` keeps `position_key`, `variant_keys`, `chrom`, and `start` top-level and packs match, identity, clinical, flag, and AF groups into structs.",
            "- With pylance/lance 7.0.0, the writer-safe default strips only the `lance-encoding:packed=true` metadata key because it panics on packed structs with all-null children; pass `--unsafe-packed-metadata` to reproduce the intended metadata path.",
            "- `minor_allele` and `minor_allele_freq` remain logical fields and are marked as constant-null candidates.",
            "- Cold Lance lookups use a `BTREE` scalar index on `position_key` with `use_scalar_index=True`.",
            "",
            "## Sources",
            "",
            "- Lance encoding strategy: <https://lance.org/format/file/encoding/>",
            "- Lance versioning: <https://lance.org/format/file/versioning/>",
            "- Lance table format: <https://lance.org/format/table/>",
            "- LanceDB scalar indexes: <https://docs.lancedb.com/indexing/scalar-index>",
        ]
    )
    if report.get("errors"):
        lines.extend(["", "## Errors", ""])
        for error in report["errors"]:
            lines.append(f"- `{error['variant']}` `{error['tier']}` `{error.get('fragment_rows')}`: {error['error']}")
    lines.append("")
    return "\n".join(lines)


def parse_csv(value: str) -> list[str]:
    return [item.strip() for item in value.split(",") if item.strip()]


def parse_fragment_rows(value: str) -> list[int]:
    if ":" in value:
        start_s, stop_s = value.split(":", 1)
        start = int(start_s)
        stop = int(stop_s)
        rows = []
        current = start
        while current <= stop:
            rows.append(current)
            current *= 2
        return rows
    return [int(item) for item in parse_csv(value)]


def run(args: argparse.Namespace) -> dict[str, Any]:
    variation_dir = args.cache_dir / "variation"
    warm_path = variation_dir / f"{args.chrom}_warm.parquet"
    cold_path = variation_dir / f"{args.chrom}_cold.parquet"
    if not warm_path.exists():
        raise FileNotFoundError(warm_path)
    if not cold_path.exists():
        raise FileNotFoundError(cold_path)

    variants = parse_csv(args.variants)
    for variant in variants:
        validate_variant(variant)
    cold_fragment_rows = parse_fragment_rows(args.cold_fragment_rows)

    args.output_dir.mkdir(parents=True, exist_ok=True)
    args.report_dir.mkdir(parents=True, exist_ok=True)

    warm_schema = pq.ParquetFile(warm_path).schema_arrow
    cold_schema = pq.ParquetFile(cold_path).schema_arrow
    warm_logical_columns = logical_columns_available(warm_schema, include_missing_variant_keys=True)
    cold_logical_columns = logical_columns_available(cold_schema, include_missing_variant_keys=True)

    print(f"Sampling {args.cold_sample_size} cold position keys from {cold_path}")
    cold_keys = sample_cold_keys(cold_path, args.cold_sample_size, args.row_limit)
    print(f"Sampled {len(cold_keys)} cold position keys")

    builds: list[dict[str, Any]] = []
    results: list[dict[str, Any]] = []
    errors: list[dict[str, Any]] = []

    print("Benchmarking current Parquet warm full scan")
    results.append(asdict(bench_parquet_full_scan(warm_path, warm_logical_columns, args.row_limit)))
    print("Benchmarking current Parquet cold point lookup")
    results.append(asdict(bench_parquet_point_lookup(cold_path, cold_logical_columns, cold_keys, args.row_limit)))

    for variant in variants:
        warm_lance_path = args.output_dir / f"{args.chrom}_warm_{variant}.lance"
        try:
            print(f"Building {variant} warm Lance dataset at {warm_lance_path}")
            warm_build = build_lance_dataset(
                warm_path,
                warm_lance_path,
                variant,
                "warm",
                args.warm_fragment_rows,
                args.warm_row_group_rows,
                args.batch_size,
                args.row_limit,
                args.force_build,
                args.skip_build,
                args.unsafe_packed_metadata,
            )
            builds.append(asdict(warm_build))
            print(f"Benchmarking {variant} warm full scan")
            results.append(asdict(bench_lance_full_scan(warm_lance_path, variant, warm_logical_columns, args.row_limit)))
        except Exception as exc:  # pragma: no cover - report operational failures.
            errors.append({"variant": variant, "tier": "warm", "fragment_rows": args.warm_fragment_rows, "error": repr(exc)})
            print(f"ERROR warm {variant}: {exc!r}")

        for fragment_rows in cold_fragment_rows:
            cold_lance_path = args.output_dir / f"{args.chrom}_cold_{variant}_fr{fragment_rows}.lance"
            try:
                print(f"Building {variant} cold Lance dataset fr={fragment_rows} at {cold_lance_path}")
                cold_build = build_lance_dataset(
                    cold_path,
                    cold_lance_path,
                    variant,
                    "cold",
                    fragment_rows,
                    fragment_rows,
                    args.batch_size,
                    args.row_limit,
                    args.force_build,
                    args.skip_build,
                    args.unsafe_packed_metadata,
                )
                builds.append(asdict(cold_build))
                print(f"Benchmarking {variant} cold point lookup fr={fragment_rows}")
                results.append(
                    asdict(
                        bench_lance_point_lookup(
                            cold_lance_path,
                            variant,
                            fragment_rows,
                            cold_logical_columns,
                            cold_keys,
                        )
                    )
                )
            except Exception as exc:  # pragma: no cover - report operational failures.
                errors.append({"variant": variant, "tier": "cold", "fragment_rows": fragment_rows, "error": repr(exc)})
                print(f"ERROR cold {variant} fr={fragment_rows}: {exc!r}")

    report = {
        "config": {
            "cache_dir": str(args.cache_dir),
            "output_dir": str(args.output_dir),
            "report_dir": str(args.report_dir),
            "chrom": args.chrom,
            "variants": variants,
            "cold_fragment_rows": cold_fragment_rows,
            "warm_fragment_rows": args.warm_fragment_rows,
            "cold_sample_size": args.cold_sample_size,
            "row_limit": args.row_limit,
            "batch_size": args.batch_size,
            "unsafe_packed_metadata": args.unsafe_packed_metadata,
        },
        "environment": {
            "lance_version": lance.__version__,
            "pyarrow_version": pa.__version__,
            "python": platform.python_version(),
            "platform": platform.platform(),
        },
        "source": {
            "warm_path": str(warm_path),
            "warm_rows": row_count_for_limit(warm_path, args.row_limit),
            "warm_total_rows": parquet_rows(warm_path),
            "cold_path": str(cold_path),
            "cold_rows": row_count_for_limit(cold_path, args.row_limit),
            "cold_total_rows": parquet_rows(cold_path),
            "warm_projection": warm_logical_columns,
            "cold_projection": cold_logical_columns,
            "cold_sampled_keys": len(cold_keys),
        },
        "builds": builds,
        "results": results,
        "errors": errors,
    }
    json_path = args.report_dir / f"{args.chrom}_lance_variation_benchmark.json"
    md_path = args.report_dir / f"{args.chrom}_lance_variation_benchmark.md"
    json_path.write_text(json.dumps(report, indent=2), encoding="utf-8")
    md_path.write_text(build_markdown_report(report), encoding="utf-8")
    print(f"Wrote {json_path}")
    print(f"Wrote {md_path}")
    return report


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--cache-dir", type=Path, default=DEFAULT_CACHE_DIR)
    parser.add_argument("--output-dir", type=Path, default=DEFAULT_OUTPUT_DIR)
    parser.add_argument("--report-dir", type=Path, default=DEFAULT_REPORT_DIR)
    parser.add_argument("--chrom", default="chr1")
    parser.add_argument("--variants", default="2.1-unpacked,2.2-packed")
    parser.add_argument("--cold-fragment-rows", default="512:16384")
    parser.add_argument("--warm-fragment-rows", type=int, default=500_000)
    parser.add_argument("--warm-row-group-rows", type=int, default=262_144)
    parser.add_argument("--cold-row-group-rows", type=int, default=8192)
    parser.add_argument("--cold-sample-size", type=int, default=2000)
    parser.add_argument("--batch-size", type=int, default=65_536)
    parser.add_argument("--row-limit", type=int)
    parser.add_argument("--force-build", action="store_true")
    parser.add_argument("--skip-build", action="store_true")
    parser.add_argument("--unsafe-packed-metadata", action="store_true")
    args = parser.parse_args()
    run(args)


if __name__ == "__main__":
    main()
