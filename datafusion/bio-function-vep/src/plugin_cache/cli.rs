//! Argument handling for the `build_plugin` driver.
//!
//! Lives in the crate (not in `examples/build_plugin.rs`) so it is unit-testable:
//! the example is a thin `parse → load_manifest → build → print` shell. An earlier
//! revision kept the argv handling in the example, where nothing could exercise it —
//! and a documented `--overwrite` flag ended up never being parsed, silently.
//!
//! Every malformed invocation is an error. In particular a valued flag with no
//! value (`… --chrom` as the last argument) is rejected rather than dropped: a
//! dangling `--chrom` would otherwise silently widen the build from one chromosome
//! to every chromosome.

use std::path::PathBuf;

use datafusion::common::{DataFusionError, Result};

use crate::plugin_cache::builder::PluginCacheBuilder;
use crate::plugin_cache::cache_manifest::CacheManifest;
use crate::plugin_cache::source_manifest::SourceManifest;

/// Flags that consume the following argv token as their value.
const VALUED_FLAGS: [&str; 5] = [
    "--manifest",
    "--variation-cache-dir",
    "--out",
    "--source-path",
    "--chrom",
];

/// Flags that stand alone (no value).
const BARE_FLAGS: [&str; 1] = ["--overwrite"];

/// Parsed `build_plugin` command line.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct PluginBuildArgs {
    /// `--manifest <path>`: the source manifest TOML.
    pub manifest: PathBuf,
    /// `--variation-cache-dir <dir>`: cache root containing `variation/`.
    pub variation_cache_dir: PathBuf,
    /// `--out <dir>`: output cache root; shards land in `<out>/plugin/<name>/`.
    pub out: PathBuf,
    /// `--source-path` values, verbatim and in order: either a bare `<path>`
    /// (single-source manifests only) or `<part>=<path>`.
    pub source_paths: Vec<String>,
    /// `--chrom <c>`, repeatable. Empty = build every chrom in the variation cache.
    pub chroms: Vec<String>,
    /// `--overwrite`: start from an empty chrom list (a clean rebuild) instead of
    /// UPSERTing into the previous `manifest.json`.
    pub overwrite: bool,
}

/// Assign a single-valued flag, rejecting a second occurrence (last-wins or
/// first-wins would silently ignore the other value).
fn set_once(slot: &mut Option<String>, value: &str, flag: &str) -> Result<()> {
    match slot {
        Some(prev) => Err(DataFusionError::Execution(format!(
            "{flag} was given twice ('{prev}' then '{value}'); it accepts a single value"
        ))),
        None => {
            *slot = Some(value.to_string());
            Ok(())
        }
    }
}

/// Unwrap a required flag into a path.
fn required(slot: Option<String>, flag: &str) -> Result<PathBuf> {
    slot.map(PathBuf::from)
        .ok_or_else(|| DataFusionError::Execution(format!("{flag} is required")))
}

impl PluginBuildArgs {
    /// Parse `build_plugin`'s arguments (argv **without** the program name).
    ///
    /// Errors on: an unrecognized flag, a valued flag with no value (it is the last
    /// argument, or is immediately followed by another flag), a repeated
    /// single-valued flag, and a missing required flag.
    pub fn parse(args: &[String]) -> Result<Self> {
        let mut manifest: Option<String> = None;
        let mut variation_cache_dir: Option<String> = None;
        let mut out: Option<String> = None;
        let mut source_paths: Vec<String> = Vec::new();
        let mut chroms: Vec<String> = Vec::new();
        let mut overwrite = false;

        let mut i = 0usize;
        while i < args.len() {
            let flag = args[i].as_str();

            if flag == "--overwrite" {
                overwrite = true;
                i += 1;
                continue;
            }
            if !VALUED_FLAGS.contains(&flag) {
                return Err(DataFusionError::Execution(format!(
                    "unrecognized argument '{flag}'; expected one of {}, {}",
                    VALUED_FLAGS.join(", "),
                    BARE_FLAGS.join(", ")
                )));
            }

            // A valued flag must be followed by a value. Silently dropping a dangling
            // flag is how a trailing `--chrom` turns a one-chromosome build into an
            // all-chromosome build without a word of warning.
            let value = match args.get(i + 1) {
                None => {
                    return Err(DataFusionError::Execution(format!(
                        "{flag} requires a value but is the last argument"
                    )));
                }
                Some(next) if next.starts_with("--") => {
                    return Err(DataFusionError::Execution(format!(
                        "{flag} requires a value but is followed by the flag '{next}'"
                    )));
                }
                Some(v) => v.as_str(),
            };

            match flag {
                "--manifest" => set_once(&mut manifest, value, flag)?,
                "--variation-cache-dir" => set_once(&mut variation_cache_dir, value, flag)?,
                "--out" => set_once(&mut out, value, flag)?,
                "--source-path" => source_paths.push(value.to_string()),
                "--chrom" => chroms.push(value.to_string()),
                // Unreachable: the `VALUED_FLAGS` guard above admits nothing else.
                other => {
                    return Err(DataFusionError::Execution(format!(
                        "internal: unhandled valued flag '{other}'"
                    )));
                }
            }
            i += 2;
        }

        Ok(Self {
            manifest: required(manifest, "--manifest")?,
            variation_cache_dir: required(variation_cache_dir, "--variation-cache-dir")?,
            out: required(out, "--out")?,
            source_paths,
            chroms,
            overwrite,
        })
    }

    /// Redirect each `[[source]]`'s `path` from a `--source-path` override.
    ///
    /// A bare `--source-path <path>` is accepted only when the manifest declares a
    /// single `[[source]]`; a multi-source manifest requires `<part>=<path>`, once
    /// per part. An unknown part is an error (a typo'd part would otherwise leave
    /// that source silently pointing at its manifest default).
    pub fn apply_source_overrides(&self, manifest: &mut SourceManifest) -> Result<()> {
        for spec in &self.source_paths {
            match spec.split_once('=') {
                // --source-path <part>=<path>
                Some((part, path)) => {
                    let target = manifest
                        .sources
                        .iter_mut()
                        .find(|s| s.part.as_deref() == Some(part))
                        .ok_or_else(|| {
                            DataFusionError::Execution(format!(
                                "--source-path '{part}=...': no [[source]] with part = \"{part}\""
                            ))
                        })?;
                    target.path = path.to_string();
                }
                // --source-path <path> — unambiguous only for a single-source manifest
                None => {
                    if manifest.sources.len() != 1 {
                        return Err(DataFusionError::Execution(format!(
                            "--source-path <path> is ambiguous: the manifest has {} sources; \
                             use --source-path <part>=<path> (parts: {})",
                            manifest.sources.len(),
                            manifest
                                .sources
                                .iter()
                                .map(|s| s.part.as_deref().unwrap_or("<none>"))
                                .collect::<Vec<_>>()
                                .join(", ")
                        )));
                    }
                    manifest.sources[0].path = spec.clone();
                }
            }
        }
        Ok(())
    }

    /// Load `--manifest` (which validates it) and apply the `--source-path` overrides.
    pub fn load_manifest(&self) -> Result<SourceManifest> {
        let mut manifest = SourceManifest::load(&self.manifest)?;
        self.apply_source_overrides(&mut manifest)?;
        Ok(manifest)
    }

    /// The manifest's file name, as recorded in the built cache manifest.
    fn manifest_file(&self) -> String {
        self.manifest
            .file_name()
            .unwrap_or(self.manifest.as_os_str())
            .to_string_lossy()
            .into_owned()
    }

    /// Build every requested chromosome and write `plugin/<name>/manifest.json`.
    pub async fn build(&self, manifest: &SourceManifest) -> Result<CacheManifest> {
        let mut builder = PluginCacheBuilder::new(
            manifest,
            self.manifest_file(),
            &self.variation_cache_dir,
            &self.out,
        );
        if !self.chroms.is_empty() {
            builder = builder.with_chrom_filter(self.chroms.clone());
        }
        builder.build_all().await
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::path::Path;

    fn argv(args: &[&str]) -> Vec<String> {
        args.iter().map(|s| s.to_string()).collect()
    }

    /// The three required flags, so each test can add only what it exercises.
    fn base(extra: &[&str]) -> Vec<String> {
        let mut v = argv(&[
            "--manifest",
            "/m.toml",
            "--variation-cache-dir",
            "/cache",
            "--out",
            "/out",
        ]);
        v.extend(argv(extra));
        v
    }

    /// `SourceManifest::load` requires a file; tests materialize one from a string.
    fn write_manifest(dir: &Path, body: &str) -> PathBuf {
        let path = dir.join("demo.source.toml");
        std::fs::write(&path, body).unwrap();
        path
    }

    #[test]
    fn parses_the_required_flags() {
        let a = PluginBuildArgs::parse(&base(&[])).unwrap();
        assert_eq!(a.manifest, PathBuf::from("/m.toml"));
        assert_eq!(a.variation_cache_dir, PathBuf::from("/cache"));
        assert_eq!(a.out, PathBuf::from("/out"));
        assert!(a.chroms.is_empty());
        assert!(a.source_paths.is_empty());
        assert!(!a.overwrite);
    }

    // The regression this module exists for: `--overwrite` was documented in the
    // example's header but never parsed, so `PluginCacheBuilder::with_overwrite`
    // had no caller outside `#[cfg(test)]` and the flag did nothing.
    #[test]
    fn parses_overwrite_as_a_bare_flag() {
        let a = PluginBuildArgs::parse(&base(&["--overwrite"])).unwrap();
        assert!(a.overwrite, "--overwrite must be parsed");
    }

    #[test]
    fn overwrite_defaults_to_false() {
        assert!(!PluginBuildArgs::parse(&base(&[])).unwrap().overwrite);
    }

    // `--overwrite` takes no value, so it must not swallow the token after it.
    #[test]
    fn overwrite_does_not_consume_the_next_argument() {
        let a = PluginBuildArgs::parse(&base(&["--overwrite", "--chrom", "21"])).unwrap();
        assert!(a.overwrite);
        assert_eq!(a.chroms, ["21"]);
    }

    #[test]
    fn chrom_repeats_and_keeps_order() {
        let a = PluginBuildArgs::parse(&base(&["--chrom", "21", "--chrom", "22"])).unwrap();
        assert_eq!(
            a.chroms,
            ["21", "22"],
            "every --chrom must be collected, in order"
        );
    }

    #[test]
    fn source_path_repeats() {
        let a = PluginBuildArgs::parse(&base(&[
            "--source-path",
            "snv=/a.tsv",
            "--source-path",
            "indel=/b.tsv",
        ]))
        .unwrap();
        assert_eq!(a.source_paths, ["snv=/a.tsv", "indel=/b.tsv"]);
    }

    // A trailing valued flag used to be silently dropped: `… --chrom` then meant
    // "no chrom filter", i.e. build EVERY chromosome. That is a silent widening.
    #[test]
    fn trailing_flag_without_a_value_is_an_error() {
        let err = PluginBuildArgs::parse(&base(&["--chrom"]))
            .expect_err("a dangling --chrom must not be silently dropped");
        let msg = err.to_string();
        assert!(msg.contains("--chrom"), "must name the flag: {msg}");
        assert!(msg.contains("last argument"), "must say why: {msg}");
    }

    #[test]
    fn valued_flag_followed_by_another_flag_is_an_error() {
        let err = PluginBuildArgs::parse(&base(&["--chrom", "--overwrite"]))
            .expect_err("--chrom --overwrite must not read '--overwrite' as a chromosome");
        assert!(err.to_string().contains("--overwrite"), "{err}");
    }

    #[test]
    fn unrecognized_flag_is_an_error() {
        let err = PluginBuildArgs::parse(&base(&["--chrome", "21"]))
            .expect_err("a typo'd flag must be rejected, not ignored");
        assert!(err.to_string().contains("--chrome"), "{err}");
    }

    #[test]
    fn repeated_single_valued_flag_is_an_error() {
        let err = PluginBuildArgs::parse(&base(&["--out", "/other"]))
            .expect_err("--out twice must not silently keep just one of them");
        assert!(err.to_string().contains("--out"), "{err}");
    }

    #[test]
    fn missing_required_flag_is_an_error() {
        for missing in ["--manifest", "--variation-cache-dir", "--out"] {
            let kept: Vec<String> = base(&[])
                .chunks(2)
                .filter(|pair| pair[0] != missing)
                .flat_map(<[String]>::to_vec)
                .collect();
            let err = PluginBuildArgs::parse(&kept).unwrap_err();
            assert!(
                err.to_string().contains(missing),
                "dropping {missing} must be an error naming it, got: {err}"
            );
        }
    }

    const TWO_SOURCE: &str = r##"
plugin_name = "demo"
coordinate_system = "1-based"
ingest_sql = "SELECT 1"

[[source]]
part = "snv"
provider = "csv"
path = "snv-default.tsv"
  [source.csv]
  schema = [{ name = "chrom", type = "Utf8" }]

[[source]]
part = "indel"
provider = "csv"
path = "indel-default.tsv"
  [source.csv]
  schema = [{ name = "chrom", type = "Utf8" }]

[[value_columns]]
column = "s"
csq_field = "S"
type = "Float32"
"##;

    const ONE_SOURCE: &str = r##"
plugin_name = "demo"
coordinate_system = "1-based"
ingest_sql = "SELECT 1"

[[source]]
provider = "csv"
path = "default.tsv"
  [source.csv]
  schema = [{ name = "chrom", type = "Utf8" }]

[[value_columns]]
column = "s"
csq_field = "S"
type = "Float32"
"##;

    /// Args for `manifest_path` plus the given extras.
    fn args_for(manifest_path: &Path, extra: &[&str]) -> PluginBuildArgs {
        let mut v = argv(&[
            "--manifest",
            manifest_path.to_str().unwrap(),
            "--variation-cache-dir",
            "/cache",
            "--out",
            "/out",
        ]);
        v.extend(argv(extra));
        PluginBuildArgs::parse(&v).unwrap()
    }

    #[test]
    fn bare_source_path_targets_the_only_source() {
        let dir = tempfile::tempdir().unwrap();
        let path = write_manifest(dir.path(), ONE_SOURCE);
        let args = args_for(&path, &["--source-path", "/real/input.tsv"]);
        let m = args.load_manifest().unwrap();
        assert_eq!(m.sources[0].path, "/real/input.tsv");
    }

    // A bare path cannot say WHICH source it means, so on a multi-source manifest it
    // must fail loudly and name the parts rather than silently patching the first.
    #[test]
    fn bare_source_path_on_a_multi_source_manifest_is_an_error() {
        let dir = tempfile::tempdir().unwrap();
        let path = write_manifest(dir.path(), TWO_SOURCE);
        let args = args_for(&path, &["--source-path", "/real/input.tsv"]);
        let err = args
            .load_manifest()
            .expect_err("a bare --source-path is ambiguous with 2 sources");
        let msg = err.to_string();
        assert!(msg.contains("ambiguous"), "{msg}");
        assert!(msg.contains("snv"), "must list the parts: {msg}");
        assert!(msg.contains("indel"), "must list the parts: {msg}");
    }

    #[test]
    fn part_scoped_source_path_targets_the_right_source() {
        let dir = tempfile::tempdir().unwrap();
        let path = write_manifest(dir.path(), TWO_SOURCE);
        let args = args_for(&path, &["--source-path", "indel=/real/indel.tsv"]);
        let m = args.load_manifest().unwrap();
        assert_eq!(
            m.sources[0].path, "snv-default.tsv",
            "the un-targeted source must keep its manifest path"
        );
        assert_eq!(m.sources[1].path, "/real/indel.tsv");
    }

    // Every --source-path must land, not just the first (the pre-fix `arg()` helper
    // read only the first occurrence).
    #[test]
    fn every_part_scoped_source_path_is_applied() {
        let dir = tempfile::tempdir().unwrap();
        let path = write_manifest(dir.path(), TWO_SOURCE);
        let args = args_for(
            &path,
            &[
                "--source-path",
                "snv=/real/snv.tsv",
                "--source-path",
                "indel=/real/indel.tsv",
            ],
        );
        let m = args.load_manifest().unwrap();
        assert_eq!(m.sources[0].path, "/real/snv.tsv");
        assert_eq!(m.sources[1].path, "/real/indel.tsv");
    }

    #[test]
    fn unknown_part_in_source_path_is_an_error() {
        let dir = tempfile::tempdir().unwrap();
        let path = write_manifest(dir.path(), TWO_SOURCE);
        let args = args_for(&path, &["--source-path", "sv=/real/sv.tsv"]);
        let err = args
            .load_manifest()
            .expect_err("an unknown part must not be silently ignored");
        let msg = err.to_string();
        assert!(msg.contains("sv"), "must name the offending part: {msg}");
    }
}
