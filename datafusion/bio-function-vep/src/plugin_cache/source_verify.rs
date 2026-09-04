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

/// What identifies a file between two looks at it without reading it: its
/// size and mtime, plus the replacement-sensitive inode and change time —
/// a copy can preserve size and mtime (`cp -p`), but an atomic replacement
/// is a new inode and any rewrite or metadata restore bumps the ctime.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
struct Fingerprint {
    size: u64,
    mtime_ns: Option<i64>,
    ino: Option<u64>,
    ctime_ns: Option<i64>,
}

impl Fingerprint {
    fn of(path: &Path) -> Result<Self> {
        let meta = std::fs::metadata(path).map_err(|e| {
            DataFusionError::Execution(format!("stat source '{}': {e}", path.display()))
        })?;
        let mtime_ns = meta
            .modified()
            .ok()
            .and_then(|t| t.duration_since(UNIX_EPOCH).ok())
            .and_then(|d| i64::try_from(d.as_nanos()).ok());
        #[cfg(unix)]
        let (ino, ctime_ns) = {
            use std::os::unix::fs::MetadataExt;
            (
                Some(meta.ino()),
                meta.ctime()
                    .checked_mul(1_000_000_000)
                    .and_then(|s| s.checked_add(meta.ctime_nsec())),
            )
        };
        #[cfg(windows)]
        let (ino, ctime_ns) = {
            use std::os::windows::fs::MetadataExt;
            // 100 ns intervals since 1601-01-01, as ns since the Unix epoch.
            const EPOCH_DIFF_100NS: i64 = 116_444_736_000_000_000;
            (
                None,
                i64::try_from(meta.creation_time())
                    .ok()
                    .and_then(|t| t.checked_sub(EPOCH_DIFF_100NS))
                    .and_then(|t| t.checked_mul(100)),
            )
        };
        #[cfg(not(any(unix, windows)))]
        let (ino, ctime_ns) = (None, None);
        Ok(Fingerprint {
            size: meta.len(),
            mtime_ns,
            ino,
            ctime_ns,
        })
    }

    /// True when `record` was taken from a file with this exact identity. A
    /// fingerprint without a change time can never match: without it, a
    /// replacement could not be told apart, so the file is re-hashed.
    fn matches(&self, record: &SourceRecord) -> bool {
        self.ctime_ns.is_some()
            && record.size == Some(self.size)
            && record.mtime_ns == self.mtime_ns
            && record.ino == self.ino
            && record.ctime_ns == self.ctime_ns
    }

    fn stamp(&self, record: &mut SourceRecord) {
        record.size = Some(self.size);
        record.mtime_ns = self.mtime_ns;
        record.ino = self.ino;
        record.ctime_ns = self.ctime_ns;
    }
}

/// Record `path`'s current identity on `record`, as a build would after
/// hashing it. For tests and tools that construct provenance records by hand.
pub fn stamp_fingerprint(path: &Path, record: &mut SourceRecord) -> Result<()> {
    Fingerprint::of(path)?.stamp(record);
    Ok(())
}

/// Check that no verified source changed since it was fingerprinted. The
/// providers reopen each source path per chromosome, so a file replaced after
/// the hash but before or during ingestion would be published under a digest
/// it does not have; the build calls this after each chromosome is ingested
/// and before its shard is committed, and refuses to go on if the identity
/// moved (size, mtime, inode or change time).
pub fn check_sources_unchanged(manifest: &SourceManifest, records: &[SourceRecord]) -> Result<()> {
    for record in records.iter().filter(|r| r.verified_md5.is_some()) {
        let Some(spec) = manifest.sources.iter().find(|s| s.part == record.part) else {
            continue;
        };
        let path = Path::new(&spec.path);
        let now = Fingerprint::of(path)?;
        if !now.matches(record) {
            return Err(DataFusionError::Execution(format!(
                "plugin '{}' {}: {} changed while it was being ingested (size {:?} -> {}, \
                 mtime_ns {:?} -> {:?}, ino {:?} -> {:?}, ctime_ns {:?} -> {:?}); the digest \
                 verified before the build no longer describes the data read, so the shard \
                 was not committed and no manifest was written. Rebuild from the stable file.",
                manifest.plugin_name,
                spec.label(),
                path.display(),
                record.size,
                now.size,
                record.mtime_ns,
                now.mtime_ns,
                record.ino,
                now.ino,
                record.ctime_ns,
                now.ctime_ns,
            )));
        }
    }
    Ok(())
}

/// Verify every `[[source]]` of `manifest` under `mode` and return the
/// provenance records to write into the built cache's manifest.
///
/// `prior` are the records of an earlier build into the same plugin directory.
/// A source whose earlier record carries a digest for a file of the same name
/// and identity (size, mtime, inode, change time) is not re-hashed: that
/// digest is reused and compared as if it had just been computed, so a
/// per-chromosome workflow that calls the builder once per contig hashes an
/// 87 GB input once, not once per call — in `Warn` mode too, where the
/// recorded digest deliberately differs from the declared one, and in `Strict`
/// mode a known mismatch fails without re-reading the file. A replacement or
/// rewrite is a new inode or a fresh change time and always re-hashes; only an
/// in-place overwrite that also forges the change time could pass, which
/// requires more than file ownership.
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
    let now = Fingerprint::of(path)?;
    let earlier = prior.iter().find_map(|p| {
        (p.part == spec.part && p.file == record.file && now.matches(p))
            .then_some(p.verified_md5.as_deref())
            .flatten()
    });
    let actual = match earlier {
        Some(digest) => {
            info!(
                "plugin '{plugin_name}' {label}: md5 {digest} computed by an earlier build over \
                 this same file (name, size, mtime, inode and ctime unchanged), not re-hashing"
            );
            digest.to_string()
        }
        None => {
            info!(
                "plugin '{plugin_name}' {label}: hashing {} ({} bytes)",
                path.display(),
                now.size
            );
            let (digest, hashed) = md5_file(path)?;
            if hashed != now.size {
                return Err(DataFusionError::Execution(format!(
                    "plugin '{plugin_name}' {label}: {} changed while it was being hashed \
                     ({} bytes at start, {hashed} read)",
                    path.display(),
                    now.size
                )));
            }
            digest
        }
    };
    let hashed = now.size;
    record.verified_md5 = Some(actual.clone());
    now.stamp(&mut record);
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
        let mut record = SourceRecord::from_spec(&manifest.sources[0]);
        record.verified_md5 = Some(earlier.into());
        Fingerprint::of(path).unwrap().stamp(&mut record);
        (manifest, vec![record])
    }

    /// Replace `path` atomically with `bytes`, preserving its mtime — the
    /// replacement a size/mtime fingerprint cannot see.
    fn replace_preserving_mtime(path: &Path, bytes: &[u8]) {
        let mtime = std::fs::metadata(path).unwrap().modified().unwrap();
        let tmp = path.with_extension("new");
        std::fs::write(&tmp, bytes).unwrap();
        std::fs::File::options()
            .write(true)
            .open(&tmp)
            .unwrap()
            .set_modified(mtime)
            .unwrap();
        std::fs::rename(&tmp, path).unwrap();
        assert_eq!(std::fs::metadata(path).unwrap().modified().unwrap(), mtime);
    }

    #[test]
    fn a_replacement_with_preserved_size_and_mtime_is_rehashed_and_detected() {
        let dir = tempfile::tempdir().unwrap();
        let src = dir.path().join("demo.tsv");
        std::fs::write(&src, b"abc").unwrap();
        let abc = "900150983cd24fb0d6963f7d28e17f72";
        let (manifest, prior) = reuse_fixture(&src, abc, abc);
        // Sanity: the record is trusted for the untouched file.
        let records = verify_sources(&manifest, SourceVerification::Strict, &prior).unwrap();
        assert_eq!(records[0].verified_md5.as_deref(), Some(abc));

        replace_preserving_mtime(&src, b"abd");
        // Same size, same mtime: only the inode / ctime give it away.
        let error = verify_sources(&manifest, SourceVerification::Strict, &prior)
            .unwrap_err()
            .to_string();
        assert!(error.contains("MD5 mismatch"), "{error}");
        assert!(
            error.contains("4911e516e5aa21d327512e0c8b197616"),
            "{error}"
        );
        // And a source replaced mid-build is caught before the shard commits.
        let error = check_sources_unchanged(&manifest, &records)
            .unwrap_err()
            .to_string();
        assert!(
            error.contains("changed while it was being ingested"),
            "{error}"
        );
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
