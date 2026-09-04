//! Build-input verification. Each resolved source file is hashed once, before
//! the first chromosome is ingested, and compared to the digest its manifest
//! declares (`md5`). A byte-parity-validated cache must not
//! silently drift from the input it was validated on — a truncated download,
//! a newer weekly release, a re-sorted copy — and publish a manifest that
//! claims the manifest's version regardless.

use std::io::{BufRead, Read};
use std::path::Path;
use std::str::FromStr;
use std::time::UNIX_EPOCH;

use datafusion::common::{DataFusionError, Result};
use log::{info, warn};
use md5::{Digest, Md5};
use noodles_csi::BinningIndex;

use crate::plugin_cache::cache_manifest::{IndexRecord, SourceRecord};
use crate::plugin_cache::source_manifest::{SourceIndex, SourceManifest, SourceSpec};

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

    /// True when a record's identity fields were taken from a file with this
    /// exact identity. A fingerprint without a change time can never match:
    /// without it, a replacement could not be told apart, so the file is
    /// re-hashed.
    fn matches_fields(
        &self,
        size: Option<u64>,
        mtime_ns: Option<i64>,
        ino: Option<u64>,
        ctime_ns: Option<i64>,
    ) -> bool {
        self.ctime_ns.is_some()
            && size == Some(self.size)
            && mtime_ns == self.mtime_ns
            && ino == self.ino
            && ctime_ns == self.ctime_ns
    }

    fn matches(&self, record: &SourceRecord) -> bool {
        self.matches_fields(record.size, record.mtime_ns, record.ino, record.ctime_ns)
    }

    fn matches_index(&self, record: &IndexRecord) -> bool {
        self.matches_fields(record.size, record.mtime_ns, record.ino, record.ctime_ns)
    }

    fn stamp(&self, record: &mut SourceRecord) {
        record.size = Some(self.size);
        record.mtime_ns = self.mtime_ns;
        record.ino = self.ino;
        record.ctime_ns = self.ctime_ns;
    }

    fn stamp_index(&self, record: &mut IndexRecord) {
        record.size = Some(self.size);
        record.mtime_ns = self.mtime_ns;
        record.ino = self.ino;
        record.ctime_ns = self.ctime_ns;
    }
}

/// The sibling index a tabix source is read through.
fn index_path(spec: &SourceSpec) -> Option<std::path::PathBuf> {
    (spec.index == Some(SourceIndex::Tabix)).then(|| format!("{}.tbi", spec.path).into())
}

/// The message for a digest that is not what the manifest declares.
fn mismatch_message(
    plugin_name: &str,
    label: &str,
    path: &Path,
    expected: &str,
    upstream: &str,
    actual: &str,
    hashed: u64,
) -> String {
    format!(
        "plugin '{plugin_name}' {label}: MD5 mismatch for {}: expected {expected} (manifest \
         md5{upstream}), got {actual} over {hashed} bytes. If this file is deliberately \
         different from the declared input (a chromosome slice, a re-compression prepared \
         as the plugin's README describes), build with source verification \"warn\" (records \
         the digest found) or \"skip\".",
        path.display()
    )
}

/// Check that a tabix index describes this data file. A stale index — one
/// built from another version of the data — still parses and still names
/// contigs, and the provider would follow it into the wrong bytes, so for
/// every contig the index names, the record its first chunk points at is read
/// from the data and must carry that contig in the index's own contig column.
/// One block per contig is decompressed; nothing else is read.
pub fn validate_tabix_index(data: &Path, index: &Path) -> Result<()> {
    let open = |p: &Path| {
        std::fs::File::open(p)
            .map_err(|e| DataFusionError::Execution(format!("open '{}': {e}", p.display())))
    };
    let idx = noodles_tabix::io::Reader::new(open(index)?)
        .read_index()
        .map_err(|e| {
            DataFusionError::Execution(format!("read tabix index '{}': {e}", index.display()))
        })?;
    let header = idx.header().ok_or_else(|| {
        DataFusionError::Execution(format!("tabix index '{}' has no header", index.display()))
    })?;
    let contig_column = header.reference_sequence_name_index();
    let comment = header.line_comment_prefix();
    let mut reader = noodles_bgzf::io::Reader::new(open(data)?);
    let mismatch = |name: &[u8], detail: String| {
        DataFusionError::Execution(format!(
            "tabix index '{}' does not describe '{}': it names contig '{}' but {detail}. The \
             index was built from a different version of the data; rebuild it with `tabix` \
             from the verified file.",
            index.display(),
            data.display(),
            String::from_utf8_lossy(name),
        ))
    };
    for (name, reference) in header
        .reference_sequence_names()
        .iter()
        .zip(idx.reference_sequences())
    {
        let Some(first) = reference
            .bins()
            .values()
            .flat_map(|bin| bin.chunks())
            .map(|chunk| chunk.start())
            .min()
        else {
            continue;
        };
        let offset = u64::from(first);
        // Seek to the block only, then bound the in-block offset by the block's
        // decompressed length ourselves: the reader would otherwise index past
        // the block (and panic) for an offset that belongs to other data.
        let (block, within): (u64, u16) = first.into();
        let block_start = noodles_bgzf::VirtualPosition::try_from((block, 0)).map_err(|e| {
            mismatch(
                name,
                format!("its first chunk ({offset:#x}) is not a position: {e}"),
            )
        })?;
        reader.seek(block_start).map_err(|e| {
            mismatch(
                name,
                format!("its first chunk ({offset:#x}) is not at a BGZF block of the data: {e}"),
            )
        })?;
        let block_len = reader
            .fill_buf()
            .map_err(|e| {
                mismatch(
                    name,
                    format!("the block at its first chunk ({offset:#x}) is unreadable: {e}"),
                )
            })?
            .len();
        if usize::from(within) > block_len {
            return Err(mismatch(
                name,
                format!(
                    "its first chunk ({offset:#x}) points {within} bytes into a block that holds \
                     only {block_len}"
                ),
            ));
        }
        reader.consume(usize::from(within));
        let mut line = Vec::new();
        // A chunk starts at a record, but tolerate a comment line or two in
        // case an indexer recorded the preceding header block.
        for _ in 0..4 {
            line.clear();
            let n = reader.read_until(b'\n', &mut line).map_err(|e| {
                mismatch(
                    name,
                    format!("the data at its first chunk ({offset:#x}) is unreadable: {e}"),
                )
            })?;
            if n == 0 {
                return Err(mismatch(
                    name,
                    format!("its first chunk ({offset:#x}) is past the end"),
                ));
            }
            if line.first() != Some(&comment) {
                break;
            }
        }
        let found = line
            .strip_suffix(b"\n")
            .unwrap_or(&line)
            .split(|&b| b == b'\t')
            .nth(contig_column)
            .unwrap_or(&[]);
        if found != &name[..] {
            return Err(mismatch(
                name,
                format!(
                    "the record at its first chunk ({offset:#x}) reads contig '{}'",
                    String::from_utf8_lossy(found)
                ),
            ));
        }
    }
    Ok(())
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
    for record in records {
        let Some(spec) = manifest.sources.iter().find(|s| s.part == record.part) else {
            continue;
        };
        let changed = |what: &str, path: &Path, before: [String; 4], now: &Fingerprint| {
            DataFusionError::Execution(format!(
                "plugin '{}' {}: {what} {} changed while it was being ingested (size {} -> {}, \
                 mtime_ns {} -> {:?}, ino {} -> {:?}, ctime_ns {} -> {:?}); the digest \
                 verified before the build no longer describes the data read, so no shard was \
                 committed and no manifest was written. Rebuild from the stable file.",
                manifest.plugin_name,
                spec.label(),
                path.display(),
                before[0],
                now.size,
                before[1],
                now.mtime_ns,
                before[2],
                now.ino,
                before[3],
                now.ctime_ns,
            ))
        };
        let path = Path::new(&spec.path);
        // The data and the index are verified independently (a manifest may
        // declare no digest for the data and still have its index hashed), so
        // each is checked on its own evidence.
        if record.verified_md5.is_some()
            && let now = Fingerprint::of(path)?
            && !now.matches(record)
        {
            return Err(changed(
                "source",
                path,
                [
                    format!("{:?}", record.size),
                    format!("{:?}", record.mtime_ns),
                    format!("{:?}", record.ino),
                    format!("{:?}", record.ctime_ns),
                ],
                &now,
            ));
        }
        if let (Some(index), Some(index_path)) = (&record.index, index_path(spec))
            && index.verified_md5.is_some()
        {
            let now = Fingerprint::of(&index_path)?;
            if !now.matches_index(index) {
                return Err(changed(
                    "index",
                    &index_path,
                    [
                        format!("{:?}", index.size),
                        format!("{:?}", index.mtime_ns),
                        format!("{:?}", index.ino),
                        format!("{:?}", index.ctime_ns),
                    ],
                    &now,
                ));
            }
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
    if mode == SourceVerification::Skip {
        info!("plugin '{plugin_name}' {label}: source verification skipped");
        return Ok(record);
    }
    let upstream = spec
        .url
        .as_deref()
        .map(|u| format!(", upstream {u}"))
        .unwrap_or_default();
    let earlier = prior
        .iter()
        .find(|p| p.part == spec.part && p.file == record.file);

    // The index decides which records each chromosome reads, so it is hashed
    // (it is small) and recorded whenever the source has one.
    if let (Some(index), Some(index_path)) = (record.index.as_mut(), index_path(spec)) {
        let now = Fingerprint::of(&index_path)?;
        let actual = match earlier
            .and_then(|p| p.index.as_ref())
            .filter(|i| i.file == index.file && now.matches_index(i))
            .and_then(|i| i.verified_md5.clone())
        {
            Some(digest) => digest,
            None => md5_file(&index_path)?.0,
        };
        index.verified_md5 = Some(actual);
        now.stamp_index(index);
        // No digest is declared for an index, so its correspondence with the
        // data is checked structurally instead; a mismatch is never a
        // deliberate derived input, so it fails in warn mode too.
        validate_tabix_index(Path::new(&spec.path), &index_path)?;
    }

    let Some(expected) = spec.expected_md5() else {
        info!("plugin '{plugin_name}' {label}: no md5 declared, source not verified");
        return Ok(record);
    };

    let path = Path::new(&spec.path);
    let now = Fingerprint::of(path)?;
    let reused = earlier
        .filter(|p| now.matches(p))
        .and_then(|p| p.verified_md5.as_deref());
    let actual = match reused {
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
    record.verified_md5 = Some(actual.clone());
    now.stamp(&mut record);
    if actual == expected {
        info!("plugin '{plugin_name}' {label}: md5 {actual} matches the manifest");
        return Ok(record);
    }

    let message = mismatch_message(
        plugin_name,
        &label,
        path,
        expected,
        &upstream,
        &actual,
        now.size,
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

    /// Write `rows` (chrom, 1-based pos, rest) as BGZF at `path` with a sibling
    /// tabix index over columns 1 and 2.
    fn write_bgzf_tabix(path: &Path, rows: &[(&str, usize, &str)]) {
        use noodles_core_tabix::Position;
        use noodles_csi::binning_index::index::reference_sequence::bin::Chunk;
        use std::io::Write;
        let mut writer = noodles_bgzf::io::Writer::new(std::fs::File::create(path).unwrap());
        let mut indexer = noodles_tabix::index::Indexer::default();
        indexer.set_header(
            noodles_csi::binning_index::index::header::Builder::vcf()
                .set_line_comment_prefix(b'#')
                .build(),
        );
        let mut start = writer.virtual_position();
        for &(chrom, pos, rest) in rows {
            writeln!(writer, "{chrom}\t{pos}\t{rest}").unwrap();
            let end = writer.virtual_position();
            let p = Position::try_from(pos).unwrap();
            indexer
                .add_record(chrom, p, p, Chunk::new(start, end))
                .unwrap();
            start = end;
        }
        writer.finish().unwrap();
        let index = std::fs::File::create(format!("{}.tbi", path.display())).unwrap();
        noodles_tabix::io::Writer::new(index)
            .write_index(&indexer.build())
            .unwrap();
    }

    #[test]
    fn a_tabix_index_from_other_data_is_rejected() {
        let dir = tempfile::tempdir().unwrap();
        let good = dir.path().join("good.vcf.gz");
        write_bgzf_tabix(&good, &[("1", 100, "a"), ("2", 200, "b")]);
        validate_tabix_index(&good, &good.with_extension("gz.tbi")).unwrap();

        // The same contigs, but a different (longer) data file: its index
        // points where `good` has other bytes or no bytes at all.
        let other = dir.path().join("other.vcf.gz");
        let rows: Vec<(&str, usize, &str)> = (0..2000)
            .map(|i| ("1", i + 1, "padding padding padding padding padding"))
            .chain([("2", 5000, "x")])
            .collect();
        write_bgzf_tabix(&other, &rows);
        let error = validate_tabix_index(&good, &other.with_extension("gz.tbi"))
            .unwrap_err()
            .to_string();
        assert!(error.contains("does not describe"), "{error}");
        assert!(error.contains("contig '2'"), "{error}");

        // And a stale index that only names contig 1 over data that has both
        // still checks what it names — the omission is the data digest's job.
        let partial = dir.path().join("partial.vcf.gz");
        write_bgzf_tabix(&partial, &[("1", 100, "a")]);
        validate_tabix_index(&partial, &partial.with_extension("gz.tbi")).unwrap();
    }

    /// A tabix VCF source with its `.tbi`, declaring `md5` for the data.
    fn indexed_fixture(dir: &Path, md5: &str) -> (SourceManifest, Path2) {
        let src = dir.join("demo.vcf.gz");
        write_bgzf_tabix(&src, &[("1", 100, "a"), ("2", 200, "b")]);
        let tbi = dir.join("demo.vcf.gz.tbi");
        let toml = format!(
            r##"
plugin_name = "demo"
coordinate_system = "1-based"
ingest_sql = "SELECT 1"

[[source]]
provider = "vcf"
path = "{}"
url = "https://example.org/demo.vcf.gz"
md5 = "{md5}"
index = "tabix"

[[value_columns]]
column = "x"
csq_field = "X"
type = "Utf8"
"##,
            src.display()
        );
        (toml::from_str(&toml).unwrap(), Path2 { src, tbi })
    }

    struct Path2 {
        src: std::path::PathBuf,
        tbi: std::path::PathBuf,
    }

    #[test]
    fn a_tabix_index_is_hashed_recorded_and_checked_like_the_data() {
        let dir = tempfile::tempdir().unwrap();
        let (probe, paths) = indexed_fixture(dir.path(), "00000000000000000000000000000000");
        let (abc, _) = md5_file(&paths.src).unwrap();
        let (index, index_size) = md5_file(&paths.tbi).unwrap();
        let mut manifest = probe;
        manifest.sources[0].md5 = Some(abc.clone());
        let records = verify_sources(&manifest, SourceVerification::Strict, &[]).unwrap();
        let rec = &records[0];
        assert_eq!(rec.verified_md5.as_deref(), Some(abc.as_str()));
        let idx = rec.index.as_ref().unwrap();
        assert_eq!(idx.file, "demo.vcf.gz.tbi");
        assert_eq!(idx.verified_md5.as_deref(), Some(index.as_str()));
        assert_eq!(idx.size, Some(index_size));
        assert!(idx.ctime_ns.is_some());
        check_sources_unchanged(&manifest, &records).unwrap();

        // The index is hashed even when the data declares no digest…
        let mut undeclared = manifest.clone();
        undeclared.sources[0].md5 = None;
        let records2 = verify_sources(&undeclared, SourceVerification::Strict, &[]).unwrap();
        assert_eq!(records2[0].verified_md5, None);
        assert_eq!(
            records2[0].index.as_ref().unwrap().verified_md5.as_deref(),
            Some(index.as_str())
        );
        // …and never in skip mode.
        let skipped = verify_sources(&manifest, SourceVerification::Skip, &[]).unwrap();
        let idx = skipped[0].index.as_ref().unwrap();
        assert_eq!(idx.file, "demo.vcf.gz.tbi");
        assert_eq!(idx.verified_md5, None);

        // An index replaced mid-build is a different file: the pre-commit
        // check refuses the build — also when the data itself declared no
        // digest and only the index was verified. The replacement is a
        // different length: an in-place overwrite of the same length can land
        // in the same coarse kernel clock tick as the fingerprint on Linux,
        // leaving mtime and ctime identical, and this test is about the check,
        // not about that window.
        let other = dir.path().join("other.vcf.gz");
        write_bgzf_tabix(&other, &[("1", 100, "a"), ("2", 200, "b"), ("3", 300, "c")]);
        std::fs::copy(other.with_extension("gz.tbi"), &paths.tbi).unwrap();
        for recs in [&records, &records2] {
            let error = check_sources_unchanged(&manifest, recs)
                .unwrap_err()
                .to_string();
            assert!(error.contains("index"), "{error}");
            assert!(
                error.contains("changed while it was being ingested"),
                "{error}"
            );
        }
        // And a fresh verification refuses an index that describes other data,
        // in warn mode as well: a stale index is never a deliberate input.
        for mode in [SourceVerification::Strict, SourceVerification::Warn] {
            let error = verify_sources(&manifest, mode, &[])
                .unwrap_err()
                .to_string();
            assert!(error.contains("does not describe"), "{error}");
        }
        let _ = paths.src;
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
