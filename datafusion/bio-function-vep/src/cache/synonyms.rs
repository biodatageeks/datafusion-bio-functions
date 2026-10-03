//! Native VEP chromosome synonyms, preserved beside converted entity shards.

use std::collections::{BTreeSet, HashMap};
use std::path::Path;

use datafusion::common::{DataFusionError, Result};

pub(crate) const CHROM_SYNONYMS_FILE: &str = "chr_synonyms.txt";

#[derive(Debug, Default)]
pub(crate) struct ChromosomeSynonyms {
    pairs: HashMap<String, BTreeSet<String>>,
}

impl ChromosomeSynonyms {
    pub(crate) fn read(cache_root: &Path) -> Result<Self> {
        let path = cache_root.join(CHROM_SYNONYMS_FILE);
        let text = match std::fs::read_to_string(&path) {
            Ok(text) => text,
            Err(error) if error.kind() == std::io::ErrorKind::NotFound => {
                return Ok(Self::default());
            }
            Err(error) => {
                return Err(DataFusionError::Execution(format!(
                    "failed to read chromosome synonyms {}: {error}",
                    path.display()
                )));
            }
        };
        let mut pairs: HashMap<String, BTreeSet<String>> = HashMap::new();
        // VEP 116.2 BaseVEP::chromosome_synonyms pairs the first token with
        // each remaining token in both directions. Unlike its database branch,
        // the file reader does not compute a transitive graph closure.
        for line in text
            .lines()
            .filter(|line| !line.trim_start().starts_with('#'))
        {
            let mut names = line.split_whitespace();
            let Some(reference) = names.next() else {
                continue;
            };
            for synonym in names {
                pairs
                    .entry(reference.to_owned())
                    .or_default()
                    .insert(synonym.to_owned());
                pairs
                    .entry(synonym.to_owned())
                    .or_default()
                    .insert(reference.to_owned());
            }
        }
        Ok(Self { pairs })
    }

    pub(crate) fn names(&self) -> impl Iterator<Item = &String> {
        self.pairs.keys()
    }

    pub(crate) fn aliases<'a>(&'a self, name: &str) -> impl Iterator<Item = &'a str> {
        // Stable ordering avoids depending on the hash seed for a malformed or
        // ambiguous synonym file. Exact source names take precedence at callers.
        self.pairs
            .get(name)
            .into_iter()
            .flatten()
            .map(String::as_str)
    }
}

/// Refresh auxiliary metadata independently of whether entity shards need work.
#[cfg(feature = "cache-builder")]
pub(crate) fn preserve_chromosome_synonyms(raw: &Path, output: &Path) -> Result<()> {
    use std::io::Write;
    #[cfg(unix)]
    use std::os::unix::fs::PermissionsExt;
    let source = raw.join(CHROM_SYNONYMS_FILE);
    let destination = output.join(CHROM_SYNONYMS_FILE);
    let bytes = match std::fs::read(&source) {
        Ok(bytes) => bytes,
        Err(error) if error.kind() == std::io::ErrorKind::NotFound => {
            // Absence is authoritative too: never reuse an older source's map.
            return match std::fs::remove_file(&destination) {
                Ok(()) => Ok(()),
                Err(error) if error.kind() == std::io::ErrorKind::NotFound => Ok(()),
                Err(error) => Err(DataFusionError::Execution(format!(
                    "failed to remove obsolete chromosome synonyms {}: {error}",
                    destination.display()
                ))),
            };
        }
        Err(error) => {
            return Err(DataFusionError::Execution(format!(
                "failed to read chromosome synonyms {}: {error}",
                source.display()
            )));
        }
    };
    #[cfg(unix)]
    let source_permissions = std::fs::metadata(&source)?.permissions();
    if std::fs::read(&destination).is_ok_and(|existing| existing == bytes) {
        #[cfg(unix)]
        if std::fs::metadata(&destination)?.permissions().mode() != source_permissions.mode() {
            // Also repair files copied by older builds with tempfile's 0600 mode.
            std::fs::set_permissions(&destination, source_permissions)?;
        }
        return Ok(());
    }
    std::fs::create_dir_all(output)?;
    // Readers must never see a partly written map during a resumed conversion.
    let mut temporary = tempfile::NamedTempFile::new_in(output)?;
    temporary.write_all(&bytes)?;
    // NamedTempFile defaults to 0600 on Unix; shared caches need the source's
    // read permissions before the completed file becomes visible to readers.
    #[cfg(unix)]
    temporary.as_file().set_permissions(source_permissions)?;
    temporary.persist(&destination).map_err(|error| {
        DataFusionError::Execution(format!(
            "failed to preserve chromosome synonyms at {}: {error}",
            destination.display()
        ))
    })?;
    Ok(())
}

#[cfg(all(test, feature = "cache-builder", unix))]
mod tests {
    use super::*;
    use std::os::unix::fs::PermissionsExt;

    #[test]
    fn copied_synonyms_preserve_source_permissions_and_repair_existing_copy() {
        let tmp = tempfile::tempdir().unwrap();
        let raw = tmp.path().join("raw");
        let output = tmp.path().join("output");
        std::fs::create_dir(&raw).unwrap();
        let source = raw.join(CHROM_SYNONYMS_FILE);
        let destination = output.join(CHROM_SYNONYMS_FILE);
        std::fs::write(&source, b"21 NC_000021.9\n").unwrap();
        std::fs::set_permissions(&source, std::fs::Permissions::from_mode(0o644)).unwrap();
        for existing in [false, true] {
            if existing {
                std::fs::set_permissions(&destination, std::fs::Permissions::from_mode(0o600))
                    .unwrap();
            }
            preserve_chromosome_synonyms(&raw, &output).unwrap();
            assert_eq!(
                std::fs::metadata(&destination)
                    .unwrap()
                    .permissions()
                    .mode()
                    & 0o777,
                0o644,
                "existing={existing}"
            );
            assert_eq!(
                std::fs::read(&destination).unwrap(),
                std::fs::read(&source).unwrap()
            );
        }
    }
}
