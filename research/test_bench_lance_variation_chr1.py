from __future__ import annotations

import pyarrow as pa

from research.bench_lance_variation_chr1 import (
    AF_GNOMADG_COLUMNS,
    ALWAYS_NULL_COLUMNS,
    EVERYTHING_VARIATION_COLUMNS,
    build_probe_starts,
    physical_projection_for_variant,
    position_key,
    result_markdown_table,
    schema_for_variant,
)


def test_everything_columns_preserve_known_null_logical_fields() -> None:
    assert "minor_allele" in EVERYTHING_VARIATION_COLUMNS
    assert "minor_allele_freq" in EVERYTHING_VARIATION_COLUMNS
    assert "minor_allele" in ALWAYS_NULL_COLUMNS
    assert "minor_allele_freq" in ALWAYS_NULL_COLUMNS


def test_21_projection_keeps_logical_top_level_names() -> None:
    logical = ["position_key", "allele_string", "gnomADg", "gnomADg_NFE"]
    assert physical_projection_for_variant("2.1-unpacked", logical) == logical


def test_22_packed_projection_maps_group_children_to_struct_once() -> None:
    logical = ["position_key", "gnomADg", "gnomADg_NFE", "AF", "AFR"]
    assert physical_projection_for_variant("2.2-packed", logical) == [
        "position_key",
        "af_gnomadg",
        "af_1kg",
    ]


def test_schema_for_22_packed_contains_packed_af_struct() -> None:
    source = pa.schema(
        [
            pa.field("position_key", pa.uint64(), nullable=False),
            pa.field("gnomADg", pa.string()),
            pa.field("gnomADg_NFE", pa.string()),
            pa.field("AF", pa.string()),
        ]
    )
    schema = schema_for_variant(
        "2.2-packed",
        source,
        ["position_key", "gnomADg", "gnomADg_NFE", "AF"],
    )
    assert "position_key" in schema.names
    assert "af_gnomadg" in schema.names
    assert schema.field("af_gnomadg").type == pa.struct(
        [pa.field(name, pa.string()) for name in ["gnomADg", "gnomADg_NFE"]]
    )
    assert schema.field("af_gnomadg").metadata[b"lance-encoding:packed"] == b"true"
    assert set(AF_GNOMADG_COLUMNS) >= {"gnomADg", "gnomADg_NFE"}


def test_result_markdown_table_formats_core_metrics() -> None:
    rows = [
        {
            "variant": "2.1-unpacked",
            "tier": "warm",
            "operation": "full_scan",
            "rows": 10,
            "seconds": 2.0,
            "rows_per_s": 5.0,
            "artifact_gib": 1.25,
        }
    ]
    markdown = result_markdown_table(rows)
    assert "| variant | tier | operation | rows | seconds | rows/s | artifact GiB |" in markdown
    assert "| 2.1-unpacked | warm | full_scan | 10 | 2.000 | 5 | 1.250 |" in markdown


def test_vcf_probe_helpers_match_extended_probe_shape() -> None:
    assert position_key("chr1", 1) == (1 << 48) | 1
    assert build_probe_starts(10, 10, "A", "G", extended=False) == [10]
    assert build_probe_starts(10, 10, "A", "G", extended=True) == [10, 11]
