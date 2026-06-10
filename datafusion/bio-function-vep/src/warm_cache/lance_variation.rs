use std::path::{Path, PathBuf};

use datafusion::arrow::datatypes::Schema;
use datafusion::common::{DataFusionError, Result};

pub const LANCE_VARIATION_DIR: &str = "variation.lance";
pub const WARM_TIER: u8 = 0;
pub const COLD_TIER: u8 = 1;

pub fn lance_variation_dataset_path(cache_root: &Path, chrom: &str) -> PathBuf {
    let bare = chrom.strip_prefix("chr").unwrap_or(chrom);
    cache_root
        .join(LANCE_VARIATION_DIR)
        .join(format!("chr{bare}.lance"))
}

pub async fn read_lance_variation_schema(path: &Path) -> Result<Schema> {
    let dataset = lance::Dataset::open(path.to_string_lossy().as_ref())
        .await
        .map_err(|err| {
            DataFusionError::Execution(format!(
                "failed to open Lance variation dataset '{}': {err}",
                path.display()
            ))
        })?;
    Ok(dataset.schema().into())
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn lance_variation_dataset_path_normalizes_chr_prefix() {
        let root = Path::new("/cache");
        assert_eq!(
            lance_variation_dataset_path(root, "chr1"),
            Path::new("/cache/variation.lance/chr1.lance")
        );
        assert_eq!(
            lance_variation_dataset_path(root, "1"),
            Path::new("/cache/variation.lance/chr1.lance")
        );
    }
}
