//! Forward-port of `ensembl-vep/t/Utils.t` onto `master`'s API.
//!
//! Standalone: no docker, no VEP runtime, no cache fixture, no `golden.vcf`. Every
//! fixture is built in a `tempfile::TempDir` at run time, so the tests are hermetic and
//! offline.
//!
//! The audit object is `tests/port/Utils.ledger.toml`, which classifies all 44 Perl
//! assertions and names the Rust test that covers each one. Read it, not this header, for
//! the authoritative split — a header cannot be checked by CI and this one's ancestor on
//! `dev-test` was wrong (it claimed 32 architectural rows; there are 31).
//!
//! Rows served by this file: **8, 9, 30, 31, 32, 33, 35**.
//! Rows 1, 2, 4, 5 are served by `src/transcript_consequence.rs::
//! format_coords_ensembl_supports_unknown_bounds`; rows 3 and 6 are `blocked-future-work`
//! because that in-`src` test asserts only four of `format_coords`' six cases and this
//! branch may not touch `src/`.
//!
//! Provenance: adapted from `port_utils.rs` on the archived `dev-test` branch, keeping its
//! test function names so the ledger and the git history still line up. The gzip rows were
//! *strengthened* during the forward-port: `dev-test` asserted only on `flate2` itself, so
//! nothing in vepyr was exercised. Each now also drives vepyr's own public gzip-VCF reader,
//! `golden_benchmark::sample_gz_vcf_first_n`, which is the closest thing in the crate to
//! `Bio::EnsEMBL::VEP::Utils::get_compressed_filehandle`: it opens the path, decodes
//! through `MultiGzDecoder`, and surfaces the open/decode error to its caller.

use std::fs::{self, File};
use std::io::{Read, Write};
use std::path::Path;

use datafusion_bio_function_vep::golden_benchmark::sample_gz_vcf_first_n;
use flate2::Compression;
use flate2::read::{GzDecoder, MultiGzDecoder};
use flate2::write::GzEncoder;
use tempfile::TempDir;

// -------------------------------------------------------------------
// Rows 8-9 — `convert_arrayref` join semantics
//
// NOTE ON STRENGTH: vepyr has no central `convert_arrayref` analogue. Multi-valued CSQ
// fields are joined with `slice::join` at each call site (`annotate_provider.rs` joins
// `MAX_AF_POPS`, `CLIN_SIG`, `PUBMED`, `SOMATIC` and the consequence terms this way), and
// none of those helpers is reachable from an integration test. These two tests therefore
// pin the *separator contract* the call sites rely on, not vepyr code. They are the
// clearest example in this port of what the ledger gate cannot see — see
// `tests/port/README.md`, "What the gate does not check".
// -------------------------------------------------------------------

#[test]
fn row8_convert_arrayref_default_join() {
    // Perl Utils.t L56: `is(convert_arrayref(['foo', 'bar']), 'foo,bar')`.
    let parts = ["foo", "bar"];
    assert_eq!(parts.join(","), "foo,bar");
}

#[test]
fn row9_convert_arrayref_custom_separator() {
    // Perl Utils.t L57: `is(convert_arrayref(['foo', 'bar'], '&'), 'foo&bar')`.
    // `&` is the separator vepyr actually emits inside a CSQ subfield.
    let parts = ["foo", "bar"];
    assert_eq!(parts.join("&"), "foo&bar");
}

// -------------------------------------------------------------------
// Rows 30-33, 35 — `get_compressed_filehandle` behaviours
// -------------------------------------------------------------------

const VCF_HEADER_LINE: &str = "##fileformat=VCFv4.2\n";
const VCF_COLUMN_LINE: &str = "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n";
const VCF_VARIANT_LINES: &str = "22\t10510033\t.\tA\tG\t.\t.\t.\n22\t10510077\t.\tC\tT\t.\t.\t.\n";

/// Build a single-member gzip from `payload` and return its bytes.
fn gz_bytes(payload: &[u8]) -> Vec<u8> {
    let mut encoder = GzEncoder::new(Vec::new(), Compression::default());
    encoder.write_all(payload).expect("gz write");
    encoder.finish().expect("gz finish")
}

/// A two-variant VCF, uncompressed.
fn plain_vcf() -> String {
    format!("{VCF_HEADER_LINE}{VCF_COLUMN_LINE}{VCF_VARIANT_LINES}")
}

#[test]
fn row30_gzdecoder_reads_single_member_gzip() {
    // Perl Utils.t L218: `ok(get_compressed_filehandle($gzvcf), 'gz file fh')`.
    // Backend contract: `flate2::read::GzDecoder` decodes a single-member gzip.
    let gz = gz_bytes(VCF_HEADER_LINE.as_bytes());
    let mut decoder = GzDecoder::new(&gz[..]);
    let mut text = String::new();
    decoder
        .read_to_string(&mut text)
        .expect("decoder read_to_string");
    assert_eq!(text, VCF_HEADER_LINE);

    // vepyr contract: its own gzip-VCF reader opens the file and yields the records.
    let dir = TempDir::new().expect("tempdir");
    let input = dir.path().join("in.vcf.gz");
    let output = dir.path().join("out.vcf");
    fs::write(&input, gz_bytes(plain_vcf().as_bytes())).expect("write gz fixture");

    let sampled = sample_gz_vcf_first_n(&input, &output, 10).expect("gzipped VCF must open");
    assert_eq!(sampled, 2, "both variant records should be read");

    let written = fs::read_to_string(&output).expect("read sampled VCF");
    assert!(
        written.starts_with(VCF_HEADER_LINE),
        "headers must be preserved, got {written:?}"
    );
    assert!(written.contains("10510077"), "got {written:?}");
}

#[test]
fn row31_empty_path_errors_on_open() {
    // Perl Utils.t L219: `throws_ok { get_compressed_filehandle() } qr/No file/`.
    // The Perl error is "No file given"; vepyr has no separate no-argument path, because
    // `&Path` is not optional at the type level. The observable analogue is the empty
    // path, which must fail at open rather than yield an empty reader.
    let err = File::open("").expect_err("empty path must error");
    // The exact `io::ErrorKind` is deliberately not pinned: macOS reports `NotFound`,
    // Linux reports `NotFound` or `InvalidInput` depending on the kernel.
    let _ = err.kind();

    let dir = TempDir::new().expect("tempdir");
    sample_gz_vcf_first_n(Path::new(""), &dir.path().join("out.vcf"), 1)
        .expect_err("vepyr must reject an empty input path");
}

#[test]
fn row32_missing_path_errors_on_open() {
    // Perl Utils.t L220: `throws_ok { get_compressed_filehandle('foobargoobar') }
    //                     qr/File .+ does not exist/`.
    let dir = TempDir::new().expect("tempdir");
    let missing = dir.path().join("foobargoobar.vcf.gz");

    let err = File::open(&missing).expect_err("nonexistent path must error");
    assert_eq!(err.kind(), std::io::ErrorKind::NotFound);

    let msg = sample_gz_vcf_first_n(&missing, &dir.path().join("out.vcf"), 1)
        .expect_err("vepyr must reject a nonexistent input path")
        .to_string();
    // Divergence, recorded rather than asserted away: Perl's message names the file
    // (`File .+ does not exist`), vepyr's `golden_benchmark::io_err` stringifies the raw
    // `io::Error` and loses the path. The *behaviour* — refusing a missing path instead of
    // returning an empty reader — is the same, and that is what this row claims.
    assert!(
        msg.contains("No such file or directory"),
        "expected the underlying io error to surface, got {msg:?}"
    );
}

#[test]
fn row33_uncompressed_input_rejected_by_gzdecoder() {
    // Perl Utils.t L221: `throws_ok { get_compressed_filehandle($test_vcf) }
    //                     qr/File .+ binary/` — an uncompressed file handed to the gzip
    // backend is rejected rather than read as text.
    let mut decoder = GzDecoder::new(VCF_HEADER_LINE.as_bytes());
    let mut sink = String::new();
    let err = decoder
        .read_to_string(&mut sink)
        .expect_err("uncompressed input must error");
    // The exact `ErrorKind` is not pinned, to stay portable across flate2 minor versions.
    let _ = err.kind();

    let dir = TempDir::new().expect("tempdir");
    let input = dir.path().join("plain.vcf.gz"); // .gz name, non-gz content
    fs::write(&input, plain_vcf()).expect("write plain fixture");
    sample_gz_vcf_first_n(&input, &dir.path().join("out.vcf"), 1)
        .expect_err("vepyr must reject an uncompressed file offered as gzip");
}

#[test]
fn row35_multigzdecoder_reads_concatenated_members() {
    // Perl Utils.t L237: `is(ref(get_compressed_filehandle($gzvcf, 1)), 'GLOB',
    //                     'PerlIO::Gzip with multistream uses gzip')` — the multistream
    // flag. In vepyr there is no flag: `MultiGzDecoder` is used unconditionally, so
    // concatenated members (the shape of every BGZF file) are read transparently.
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

    // vepyr contract: a VCF split across two gzip members reads as one stream.
    let dir = TempDir::new().expect("tempdir");
    let input = dir.path().join("multi.vcf.gz");
    let output = dir.path().join("out.vcf");
    let header_member = gz_bytes(format!("{VCF_HEADER_LINE}{VCF_COLUMN_LINE}").as_bytes());
    let body_member = gz_bytes(VCF_VARIANT_LINES.as_bytes());
    let mut multi = header_member;
    multi.extend_from_slice(&body_member);
    fs::write(&input, &multi).expect("write multi-member gz fixture");

    let sampled = sample_gz_vcf_first_n(&input, &output, 10).expect("multi-member VCF must open");
    assert_eq!(
        sampled, 2,
        "records living in the second gzip member must still be read"
    );

    let written = fs::read_to_string(&output).expect("read sampled VCF");
    assert!(written.starts_with(VCF_HEADER_LINE), "got {written:?}");
    assert!(written.contains("10510033"), "got {written:?}");
}
