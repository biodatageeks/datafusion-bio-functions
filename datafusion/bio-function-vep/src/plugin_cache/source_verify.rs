//! Build-input verification. Each resolved source file is hashed once, before
//! the first chromosome is ingested, and compared to the digest its manifest
//! declares (`path_md5`, else `md5`). A byte-parity-validated cache must not
//! silently drift from the input it was validated on — a truncated download,
//! a newer weekly release, a re-sorted copy — and publish a manifest that
//! claims the manifest's version regardless.

use std::io::Read;
use std::path::Path;
use std::str::FromStr;
use std::time::UNIX_EPOCH;

use datafusion::common::{DataFusionError, Result};
use log::{info, warn};
use md5::{Digest, Md5};

use crate::plugin_cache::cache_manifest::SourceRecord;
use crate::plugin_cache::source_manifest::{SourceManifest, SourceSpec};

/// What to do when a source's digest is declared.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Default)]
pub enum SourceVerification {
    /// Hash and fail the build on a mismatch (default).
    #[default]
    Strict,
    /// Hash, log a mismatch and record the actual digest, but keep building —
    /// for deliberate builds against a derived or re-compressed artifact.
    Warn,
    /// Do not hash. Provenance keys are still copied into the built manifest,
    /// without a `verified_md5`. For chromosome slices of a large source,
    /// whose digest can never match the whole file's.
    Skip,
}

impl SourceVerification {
    pub fn as_str(self) -> &'static str {
        match self {
            SourceVerification::Strict => "strict",
            SourceVerification::Warn => "warn",
            SourceVerification::Skip => "skip",
        }
    }
}

impl FromStr for SourceVerification {
    type Err = DataFusionError;

    fn from_str(s: &str) -> Result<Self> {
        match s {
            "strict" => Ok(SourceVerification::Strict),
            "warn" => Ok(SourceVerification::Warn),
            "skip" => Ok(SourceVerification::Skip),
            other => Err(DataFusionError::Execution(format!(
                "unknown source verification mode {other:?}; expected \"strict\", \"warn\" or \
                 \"skip\""
            ))),
        }
    }
}

/// Streaming MD5 of a file in bounded memory. Returns the lowercase hex digest
/// and the number of bytes hashed.
pub fn md5_file(path: &Path) -> Result<(String, u64)> {
    let mut file = std::fs::File::open(path).map_err(|e| {
        DataFusionError::Execution(format!("open source '{}' for hashing: {e}", path.display()))
    })?;
    let mut hasher = Md5::new();
    let mut buf = vec![0u8; 1 << 20];
    let mut size = 0u64;
    loop {
        let n = file.read(&mut buf).map_err(|e| {
            DataFusionError::Execution(format!("read source '{}' for hashing: {e}", path.display()))
        })?;
        if n == 0 {
            break;
        }
        hasher.update(&buf[..n]);
        size += n as u64;
    }
    Ok((format!("{:x}", hasher.finalize()), size))
}

/// `(size, mtime_ns)` of a file: with the file name, the fingerprint an
/// earlier verification is keyed by. Full mtime precision, so two writes
/// within one second are not confused.
fn fingerprint(path: &Path) -> Result<(u64, Option<i64>)> {
    let meta = std::fs::metadata(path).map_err(|e| {
        DataFusionError::Execution(format!("stat source '{}': {e}", path.display()))
    })?;
    let mtime_ns = meta
        .modified()
        .ok()
        .and_then(|t| t.duration_since(UNIX_EPOCH).ok())
        .and_then(|d| i64::try_from(d.as_nanos()).ok());
    Ok((meta.len(), mtime_ns))
}

/// Check that no verified source changed since it was fingerprinted. The
/// providers reopen each source path per chromosome, so a file replaced after
/// the hash but before or during ingestion would be published under a digest
/// it does not have; the build calls this after every chromosome and refuses
/// to write the manifest if a size or mtime moved. A replacement with the
/// same size and mtime is the same blind spot as the reuse rule's.
pub fn check_sources_unchanged(manifest: &SourceManifest, records: &[SourceRecord]) -> Result<()> {
    for record in records.iter().filter(|r| r.verified_md5.is_some()) {
        let Some(spec) = manifest.sources.iter().find(|s| s.part == record.part) else {
            continue;
        };
        let path = Path::new(&spec.path);
        let (size, mtime_ns) = fingerprint(path)?;
        if record.size != Some(size) || record.mtime_ns != mtime_ns {
            return Err(DataFusionError::Execution(format!(
                "plugin '{}' {}: {} changed while it was being ingested (size {:?} -> {size}, \
                 mtime_ns {:?} -> {mtime_ns:?}); the digest verified before the build no longer \
                 describes the data read, so no manifest was written. Rebuild from the stable \
                 file.",
                manifest.plugin_name,
                spec.label(),
                path.display(),
                record.size,
                record.mtime_ns,
            )));
        }
    }
    Ok(())
}

/// Verify every `[[source]]` of `manifest` under `mode` and return the
/// provenance records to write into the built cache's manifest.
///
/// `prior` are the records of an earlier build into the same plugin directory.
/// A source whose earlier record carries a digest for the same file name, size
/// and mtime is not re-hashed: that digest is reused and compared as if it had
/// just been computed, so a per-chromosome workflow that calls the builder once
/// per contig hashes an 87 GB input once, not once per call — in `Warn` mode
/// too, where the recorded digest deliberately differs from the declared one,
/// and in `Strict` mode a known mismatch fails without re-reading the file.
/// Like any mtime-keyed cache this cannot tell apart a replacement copied with
/// its timestamps preserved and the same byte length; a different name, size
/// or mtime always re-hashes.
pub fn verify_sources(
    manifest: &SourceManifest,
    mode: SourceVerification,
    prior: &[SourceRecord],
) -> Result<Vec<SourceRecord>> {
    manifest
        .sources
        .iter()
        .map(|spec| verify_source(&manifest.plugin_name, spec, mode, prior))
        .collect()
}

fn verify_source(
    plugin_name: &str,
    spec: &SourceSpec,
    mode: SourceVerification,
    prior: &[SourceRecord],
) -> Result<SourceRecord> {
    let mut record = SourceRecord::from_spec(spec);
    let label = spec.label();
    let Some(expected) = spec.expected_md5() else {
        info!("plugin '{plugin_name}' {label}: no md5 declared, source not verified");
        return Ok(record);
    };
    if mode == SourceVerification::Skip {
        info!("plugin '{plugin_name}' {label}: source verification skipped");
        return Ok(record);
    }

    let path = Path::new(&spec.path);
    let (size, mtime_ns) = fingerprint(path)?;
    let earlier = prior.iter().find_map(|p| {
        (p.part == spec.part
            && p.file == record.file
            && p.size == Some(size)
            && mtime_ns.is_some()
            && p.mtime_ns == mtime_ns)
            .then_some(p.verified_md5.as_deref())
            .flatten()
    });
    let (actual, hashed) = match earlier {
        Some(digest) => {
            info!(
                "plugin '{plugin_name}' {label}: md5 {digest} computed by an earlier build \
                 (same file name, size and mtime), not re-hashing"
            );
            (digest.to_string(), size)
        }
        None => {
            info!(
                "plugin '{plugin_name}' {label}: hashing {} ({size} bytes)",
                path.display()
            );
            md5_file(path)?
        }
    };
    record.verified_md5 = Some(actual.clone());
    record.size = Some(hashed);
    record.mtime_ns = mtime_ns;
    if actual == expected {
        info!("plugin '{plugin_name}' {label}: md5 {actual} matches the manifest");
        return Ok(record);
    }

    let declared_as = if spec.path_md5.is_some() {
        "path_md5"
    } else {
        "md5"
    };
    let upstream = spec
        .url
        .as_deref()
        .map(|u| format!(", upstream {u}"))
        .unwrap_or_default();
    let message = format!(
        "plugin '{plugin_name}' {label}: MD5 mismatch for {}: expected {expected} (manifest \
         {declared_as}{upstream}), got {actual} over {hashed} bytes. If this file is \
         deliberately different from the declared input (a chromosome slice, a \
         re-compression), build with source verification \"warn\" or \"skip\", or declare its \
         digest as path_md5 in the manifest.",
        path.display()
    );
    match mode {
        SourceVerification::Strict => Err(DataFusionError::Execution(message)),
        SourceVerification::Warn => {
            warn!("{message}");
            Ok(record)
        }
        SourceVerification::Skip => unreachable!("skip returns before hashing"),
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn md5_file_matches_known_vectors() {
        let dir = tempfile::tempdir().unwrap();
        let empty = dir.path().join("empty");
        std::fs::write(&empty, b"").unwrap();
        assert_eq!(
            md5_file(&empty).unwrap(),
            ("d41d8cd98f00b204e9800998ecf8427e".to_string(), 0)
        );
        let abc = dir.path().join("abc");
        std::fs::write(&abc, b"abc").unwrap();
        assert_eq!(
            md5_file(&abc).unwrap(),
            ("900150983cd24fb0d6963f7d28e17f72".to_string(), 3)
        );
    }

    #[test]
    fn md5_file_streams_past_one_buffer() {
        let dir = tempfile::tempdir().unwrap();
        let big = dir.path().join("big");
        let body = vec![b'a'; (1 << 20) + 12345];
        std::fs::write(&big, &body).unwrap();
        let (digest, size) = md5_file(&big).unwrap();
        assert_eq!(size, body.len() as u64);
        assert_eq!(digest, format!("{:x}", Md5::digest(&body)));
    }

    /// A one-source manifest over `path` declaring `md5`, plus the record an
    /// earlier build would have written for that file with digest `earlier`.
    fn reuse_fixture(path: &Path, md5: &str, earlier: &str) -> (SourceManifest, Vec<SourceRecord>) {
        let toml = format!(
            r##"
plugin_name = "demo"
coordinate_system = "1-based"
ingest_sql = "SELECT 1"

[[source]]
provider = "csv"
path = "{}"
url = "https://example.org/demo.tsv"
md5 = "{md5}"

[[value_columns]]
column = "x"
csq_field = "X"
type = "Utf8"
"##,
            path.display()
        );
        let manifest: SourceManifest = toml::from_str(&toml).unwrap();
        let (size, mtime_ns) = fingerprint(path).unwrap();
        let prior = vec![SourceRecord {
            verified_md5: Some(earlier.into()),
            size: Some(size),
            mtime_ns,
            ..SourceRecord::from_spec(&manifest.sources[0])
        }];
        (manifest, prior)
    }

    // The earlier build's digest is a property of the file, so it is reused
    // whatever it was — warn mode's deliberate mismatch included — and a
    // strict build with a known mismatch fails without re-reading the file.
    #[test]
    fn an_earlier_digest_is_reused_regardless_of_whether_it_matched() {
        let dir = tempfile::tempdir().unwrap();
        let src = dir.path().join("demo.tsv");
        std::fs::write(&src, b"abc").unwrap();
        let real = "900150983cd24fb0d6963f7d28e17f72";
        let recorded = "11111111111111111111111111111111";
        let declared = "22222222222222222222222222222222";

        // Warn: the recorded (mismatching) digest is carried forward, not `real`.
        let (manifest, prior) = reuse_fixture(&src, declared, recorded);
        let records = verify_sources(&manifest, SourceVerification::Warn, &prior).unwrap();
        assert_eq!(records[0].verified_md5.as_deref(), Some(recorded));

        // Strict: the known mismatch is reported from the record, not re-hashed.
        let error = verify_sources(&manifest, SourceVerification::Strict, &prior)
            .unwrap_err()
            .to_string();
        assert!(error.contains(recorded), "{error}");
        assert!(!error.contains(real), "{error}");

        // Strict with a recorded match is accepted without hashing either.
        let (manifest, prior) = reuse_fixture(&src, declared, declared);
        let records = verify_sources(&manifest, SourceVerification::Strict, &prior).unwrap();
        assert_eq!(records[0].verified_md5.as_deref(), Some(declared));

        // Without a usable record the file is hashed for real.
        let records = verify_sources(&manifest, SourceVerification::Warn, &[]).unwrap();
        assert_eq!(records[0].verified_md5.as_deref(), Some(real));
    }

    #[test]
    fn a_source_that_moves_after_verification_is_detected() {
        let dir = tempfile::tempdir().unwrap();
        let src = dir.path().join("demo.tsv");
        std::fs::write(&src, b"abc").unwrap();
        let (manifest, _) = reuse_fixture(&src, "900150983cd24fb0d6963f7d28e17f72", "x");
        let records = verify_sources(&manifest, SourceVerification::Strict, &[]).unwrap();
        check_sources_unchanged(&manifest, &records).unwrap();

        // Same bytes rewritten later: the mtime moves, and that is enough.
        let file = std::fs::File::options().write(true).open(&src).unwrap();
        file.set_modified(
            std::time::SystemTime::UNIX_EPOCH + std::time::Duration::from_secs(1_700_000_000),
        )
        .unwrap();
        drop(file);
        let error = check_sources_unchanged(&manifest, &records)
            .unwrap_err()
            .to_string();
        assert!(
            error.contains("changed while it was being ingested"),
            "{error}"
        );

        // A skipped source has no digest to protect and is not checked.
        let skipped = verify_sources(&manifest, SourceVerification::Skip, &[]).unwrap();
        std::fs::write(&src, b"abcd").unwrap();
        check_sources_unchanged(&manifest, &skipped).unwrap();
    }

    #[test]
    fn mode_parses_its_three_spellings_only() {
        assert_eq!(
            "strict".parse::<SourceVerification>().unwrap(),
            SourceVerification::Strict
        );
        assert_eq!(
            "warn".parse::<SourceVerification>().unwrap(),
            SourceVerification::Warn
        );
        assert_eq!(
            "skip".parse::<SourceVerification>().unwrap(),
            SourceVerification::Skip
        );
        let error = "Strict"
            .parse::<SourceVerification>()
            .unwrap_err()
            .to_string();
        assert!(
            error.contains("unknown source verification mode"),
            "{error}"
        );
        assert_eq!(SourceVerification::default(), SourceVerification::Strict);
    }
}
