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

/// `(size, mtime)` of a file: the fingerprint an earlier verification is keyed by.
fn fingerprint(path: &Path) -> Result<(u64, Option<i64>)> {
    let meta = std::fs::metadata(path).map_err(|e| {
        DataFusionError::Execution(format!("stat source '{}': {e}", path.display()))
    })?;
    let mtime = meta
        .modified()
        .ok()
        .and_then(|t| t.duration_since(UNIX_EPOCH).ok())
        .map(|d| d.as_secs() as i64);
    Ok((meta.len(), mtime))
}

/// Verify every `[[source]]` of `manifest` under `mode` and return the
/// provenance records to write into the built cache's manifest.
///
/// `prior` are the records of an earlier build into the same plugin directory.
/// A source whose earlier record carries the expected digest and the file's
/// current `(size, mtime)` is trusted without re-hashing, so a per-chromosome
/// workflow that calls the builder once per contig hashes an 87 GB input once,
/// not once per call.
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
    let (size, mtime) = fingerprint(path)?;
    let earlier = prior.iter().find(|p| {
        p.part == spec.part
            && p.verified_md5.as_deref() == Some(expected)
            && p.size == Some(size)
            && mtime.is_some()
            && p.mtime == mtime
    });
    if earlier.is_some() {
        info!(
            "plugin '{plugin_name}' {label}: md5 {expected} verified by an earlier build \
             (size and mtime unchanged), not re-hashing"
        );
        record.verified_md5 = Some(expected.to_string());
        record.size = Some(size);
        record.mtime = mtime;
        return Ok(record);
    }

    info!(
        "plugin '{plugin_name}' {label}: hashing {} ({size} bytes)",
        path.display()
    );
    let (actual, hashed) = md5_file(path)?;
    record.verified_md5 = Some(actual.clone());
    record.size = Some(hashed);
    record.mtime = mtime;
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
