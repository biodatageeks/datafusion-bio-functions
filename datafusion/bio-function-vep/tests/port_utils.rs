//! v2-paradigm port of `ensembl-vep/t/Utils.t`.
//!
//! Standalone (no docker, no `port_common`, no `golden.vcf`, no cache
//! fixture). Covers the `unit-port` rows from
//! `porting-tests/detailed_plans/Utils.md`:
//!
//! - Rows 1-6 (`format_coords`) — see in-file unit test
//!   `format_coords_ensembl_supports_unknown_bounds` in
//!   `src/transcript_consequence.rs` (line ~13738). Not duplicated here.
//! - Rows 8-9 (`convert_arrayref` join semantics) — see this file.
//! - Rows 30-33, 35 (gzip reader behaviours) — see this file.
//!
//! The 32 `architectural-no-analogue` rows (7, 10-29, 34, 36-44) are
//! documented in the comment block at the bottom of this file; the
//! substantive justifications live in
//! `porting-tests/detailed_plans/Utils.md` §Architectural-no-analogue
//! justifications.
//!
//! Plan: `porting-tests/plans/2026-05-27-port-utils.md`.

use std::fs::File;
use std::io::{Read, Write};

use flate2::Compression;
use flate2::read::{GzDecoder, MultiGzDecoder};
use flate2::write::GzEncoder;

// -------------------------------------------------------------------
// Rows 8-9 — `convert_arrayref` join semantics
// -------------------------------------------------------------------

#[test]
fn row8_convert_arrayref_default_join() {
    // Perl Utils.t L56: `is(convert_arrayref(['foo', 'bar']), 'foo,bar')`.
    // Vepyr analogue: stdlib `slice::join(",")`. No central helper.
    let parts = ["foo", "bar"];
    assert_eq!(parts.join(","), "foo,bar");
}

#[test]
fn row9_convert_arrayref_custom_separator() {
    // Perl Utils.t L57: `is(convert_arrayref(['foo', 'bar'], '&'), 'foo&bar')`.
    // Vepyr analogue: stdlib `slice::join("&")`.
    let parts = ["foo", "bar"];
    assert_eq!(parts.join("&"), "foo&bar");
}

// -------------------------------------------------------------------
// Rows 30-33, 35 — gzip reader behaviours (flate2)
// -------------------------------------------------------------------

const VCF_HEADER_LINE: &str = "##fileformat=VCFv4.2\n";

/// Helper: build a single-member gzip from `payload` and return its bytes.
fn gz_bytes(payload: &[u8]) -> Vec<u8> {
    let mut encoder = GzEncoder::new(Vec::new(), Compression::default());
    encoder.write_all(payload).expect("gz write");
    encoder.finish().expect("gz finish")
}

#[test]
fn row30_gzdecoder_reads_single_member_gzip() {
    // Perl Utils.t L218: `ok(get_compressed_filehandle($gzvcf), 'gz file fh')`.
    // Vepyr analogue: `flate2::read::GzDecoder` decodes a single-member gzip.
    let gz = gz_bytes(VCF_HEADER_LINE.as_bytes());
    let mut decoder = GzDecoder::new(&gz[..]);
    let mut text = String::new();
    decoder
        .read_to_string(&mut text)
        .expect("decoder read_to_string");
    assert_eq!(text, VCF_HEADER_LINE);
}

#[test]
fn row31_empty_path_errors_on_open() {
    // Perl Utils.t L219: `throws_ok { get_compressed_filehandle() } qr/No file/`.
    // Vepyr analogue: `File::open("")` returns `Err` (empty path is invalid).
    let err = File::open("").expect_err("empty path must error");
    // We do not pin the exact io::ErrorKind: macOS reports `NotFound`,
    // Linux reports `NotFound` or `InvalidInput` depending on kernel.
    // The contract is that opening "" is observably an error.
    let _ = err.kind();
}

#[test]
fn row32_missing_path_errors_on_open() {
    // Perl Utils.t L220: `throws_ok { get_compressed_filehandle('foobargoobar') }
    //                     qr/does not exist/`.
    // Vepyr analogue: `File::open("foobargoobar")` returns `Err(NotFound)`.
    let err = File::open("foobargoobar").expect_err("nonexistent path must error");
    assert_eq!(err.kind(), std::io::ErrorKind::NotFound);
}

#[test]
fn row33_uncompressed_input_rejected_by_gzdecoder() {
    // Perl Utils.t L221: `throws_ok { get_compressed_filehandle($test_vcf) }
    //                     qr/binary/` — opening an uncompressed file via
    // the gzip backend rejects it.
    // Vepyr analogue: `flate2::read::GzDecoder` returns `Err` when the
    // payload lacks the gzip magic bytes 0x1f 0x8b.
    let plain = VCF_HEADER_LINE.as_bytes();
    let mut decoder = GzDecoder::new(plain);
    let mut sink = String::new();
    let result = decoder.read_to_string(&mut sink);
    let err = result.expect_err("uncompressed input must error");
    // flate2 reports `InvalidInput` for bad header bytes on most platforms.
    // We do not pin the exact ErrorKind to keep the test portable across
    // flate2 minor versions; we assert only that the read failed.
    let _ = err.kind();
}

#[test]
fn row35_multigzdecoder_reads_concatenated_members() {
    // Perl Utils.t L237-240: gzip multi-stream semantics — read of a
    // concatenation of two gzip members recovers the concatenated payload.
    // Vepyr analogue: `flate2::read::MultiGzDecoder` (already used in
    // `src/golden_benchmark.rs:15,134,183`).
    let m1 = gz_bytes(b"member-one\n");
    let m2 = gz_bytes(b"member-two\n");
    let mut combined = Vec::with_capacity(m1.len() + m2.len());
    combined.extend_from_slice(&m1);
    combined.extend_from_slice(&m2);

    let mut decoder = MultiGzDecoder::new(&combined[..]);
    let mut text = String::new();
    decoder
        .read_to_string(&mut text)
        .expect("multi-stream read_to_string");
    assert_eq!(text, "member-one\nmember-two\n");
}

// -------------------------------------------------------------------
// Architectural-no-analogue rows (sztywno-1:1 audit trail).
// Substantive justifications live in
// `porting-tests/detailed_plans/Utils.md` §Architectural-no-analogue
// justifications. Each line below names the missing-by-design Perl
// construct.
// -------------------------------------------------------------------
//
// Row  7: `convert_arrayref` polymorphism — Perl polymorphic dispatch of
//         one symbol over scalar OR arrayref input. Vepyr call sites pass
//         typed `&[T]` or scalar already; no polymorphic helper to test.
//
// Rows 10-15: `numberify` — Perl JSON-pre-serialization helper coercing
//         stringified numbers to native numeric types in nested
//         structures, with exempt-key map. Vepyr emits typed Arrow
//         arrays (Int64/Float64/Utf8); no string-typed intermediate to
//         coerce.
//
// Rows 16-21: `merge_hashes` — Perl per-VFOA nested-hash accumulator
//         with optional add-mode for numerics. Vepyr accumulates into
//         typed Arrow ArrayBuilders per CSQ field; no nested-hash merge
//         primitive exists.
//
// Rows 22-24: `merge_arrays` — Perl concat + dedupe-preserve-order array
//         merger. Vepyr uses `Vec::extend` + `Vec::dedup` directly at
//         call sites; no centralised helper.
//
// Rows 25-29: `find_in_ref` — Perl recursive nested-ref walker used by
//         VEP's plugin output flatten layer. Vepyr has no plugin system
//         and no dynamic introspection of Perl-style nested hashrefs.
//
// Row 34: `PerlIO::Gzip` backend return-type (GLOB) — vepyr uses
//         `flate2` unconditionally; no per-backend dispatch decision.
//
// Row 36: `gzip(1)` binary-fallback backend return-type (GLOB) — same:
//         single Rust backend.
//
// Row 37: `IO::Uncompress::Gunzip` backend return-type — same.
//
// Row 38: "no backend available" error — impossible in vepyr: flate2 is
//         a compile-time dependency.
//
// Row 39: `get_time()` formatted wall-clock string — vepyr exposes no
//         public timestamp helper. Logging timestamps live inside
//         `env_logger`'s formatter.
//
// Rows 40-44: `get_version_data` / `get_version_string` — Perl
//         multi-component-version-directory parser
//         (`dir/{component}/version` text files into
//         `{component => {release, sub}}`). Vepyr uses a single
//         `env!("CARGO_PKG_VERSION")` semver string; no sibling
//         Ensembl-API components have version files to parse.
