use std::collections::{BTreeMap, HashMap};
use std::fs;
use std::path::{Path, PathBuf};

use anyhow::{Context, Result, bail};
use lance_file::version::LanceFileVersion;
use serde::{Deserialize, Serialize};

use crate::payload_sidecar::HeedPayloadBackend;

#[derive(Debug, Clone, Deserialize, Serialize)]
#[serde(deny_unknown_fields)]
pub struct SandboxConfig {
    pub dataset: DatasetConfig,
    pub key: KeyConfig,
    pub indexes: IndexConfig,
    pub benchmark: BenchmarkConfig,
    #[serde(default)]
    pub sidecar: SidecarConfig,
    pub inspect: InspectConfig,
    pub defaults: DefaultsConfig,
    #[serde(default)]
    pub fields: BTreeMap<String, FieldConfig>,
    #[serde(default)]
    pub structs: BTreeMap<String, StructPackConfig>,
}

#[derive(Debug, Clone, Deserialize, Serialize)]
#[serde(deny_unknown_fields)]
pub struct DatasetConfig {
    pub name: String,
    pub cache_root: PathBuf,
    pub chrom: String,
    pub output_root: PathBuf,
    pub lance_version: String,
    pub batch_size: usize,
    #[serde(default)]
    pub overwrite: bool,
}

#[derive(Debug, Clone, Deserialize, Serialize)]
#[serde(deny_unknown_fields)]
pub struct KeyConfig {
    pub mode: KeyMode,
    pub column: String,
    #[serde(rename = "type")]
    pub data_type: KeyDataType,
}

#[derive(Debug, Clone, Copy, Deserialize, Serialize, PartialEq, Eq)]
#[serde(rename_all = "snake_case")]
pub enum KeyMode {
    Position,
    PositionKey,
}

#[derive(Debug, Clone, Copy, Deserialize, Serialize, PartialEq, Eq)]
#[serde(rename_all = "snake_case")]
pub enum KeyDataType {
    Uint32,
    Uint64,
}

#[derive(Debug, Clone, Deserialize, Serialize)]
#[serde(deny_unknown_fields)]
pub struct IndexConfig {
    #[serde(default = "default_true")]
    pub position_btree: bool,
    #[serde(default = "default_true")]
    pub tier_bitmap: bool,
}

#[derive(Debug, Clone, Deserialize, Serialize)]
#[serde(deny_unknown_fields)]
pub struct BenchmarkConfig {
    pub positions_file: PathBuf,
    pub vcf_path: PathBuf,
    pub lookup_batch_sizes: Vec<usize>,
    #[serde(default = "default_true")]
    pub sort_position_keys: bool,
}

#[derive(Debug, Clone, Deserialize, Serialize)]
#[serde(deny_unknown_fields)]
pub struct SidecarConfig {
    #[serde(default)]
    pub payload_backends: Vec<HeedPayloadBackend>,
    #[serde(default = "default_heed_map_size_gib")]
    pub heed_map_size_gib: usize,
}

impl Default for SidecarConfig {
    fn default() -> Self {
        Self {
            payload_backends: Vec::new(),
            heed_map_size_gib: default_heed_map_size_gib(),
        }
    }
}

#[derive(Debug, Clone, Deserialize, Serialize)]
#[serde(deny_unknown_fields)]
pub struct InspectConfig {
    #[serde(default = "default_true")]
    pub include_physical_columns: bool,
    #[serde(default = "default_true")]
    pub include_index_sizes: bool,
}

#[derive(Debug, Clone, Deserialize, Serialize)]
#[serde(deny_unknown_fields)]
pub struct DefaultsConfig {
    pub encoding: EncodingConfig,
}

#[derive(Debug, Clone, Default, Deserialize, Serialize)]
#[serde(deny_unknown_fields)]
pub struct FieldConfig {
    #[serde(default)]
    pub encoding: EncodingOverride,
}

#[derive(Debug, Clone, Deserialize, Serialize)]
#[serde(deny_unknown_fields)]
pub struct EncodingConfig {
    pub structural: String,
    pub compression: String,
    pub compression_level: i32,
    pub dict_values_compression: String,
    pub dict_values_compression_level: i32,
    pub rle_threshold: f64,
    pub dict_size_ratio: f64,
    pub dict_divisor: i32,
    pub minichunk_size: i64,
}

#[derive(Debug, Clone, Default, Deserialize, Serialize)]
#[serde(deny_unknown_fields)]
pub struct EncodingOverride {
    pub structural: Option<String>,
    pub compression: Option<String>,
    pub compression_level: Option<i32>,
    pub dict_values_compression: Option<String>,
    pub dict_values_compression_level: Option<i32>,
    pub rle_threshold: Option<f64>,
    pub dict_size_ratio: Option<f64>,
    pub dict_divisor: Option<i32>,
    pub minichunk_size: Option<i64>,
}

#[derive(Debug, Clone, Deserialize, Serialize)]
#[serde(deny_unknown_fields)]
pub struct StructPackConfig {
    #[serde(default)]
    pub enabled: bool,
    #[serde(default = "default_true")]
    pub packed_metadata: bool,
    pub fields: Vec<String>,
}

fn default_true() -> bool {
    true
}

fn default_heed_map_size_gib() -> usize {
    64
}

impl SandboxConfig {
    pub fn read(path: &Path) -> Result<Self> {
        let raw = fs::read_to_string(path)
            .with_context(|| format!("failed to read config '{}'", path.display()))?;
        let config: Self = toml::from_str(&raw)
            .with_context(|| format!("failed to parse config '{}'", path.display()))?;
        config.validate()?;
        Ok(config)
    }

    pub fn run_root(&self) -> PathBuf {
        self.dataset.output_root.join(&self.dataset.name)
    }

    pub fn dataset_path(&self) -> PathBuf {
        self.run_root()
            .join(format!("{}.lance", self.dataset.chrom))
    }

    pub fn resolved_positions_file(&self, config_path: &Path) -> PathBuf {
        if self.benchmark.positions_file.is_absolute() {
            self.benchmark.positions_file.clone()
        } else {
            config_path
                .parent()
                .and_then(|p| p.parent())
                .unwrap_or_else(|| Path::new("."))
                .join(&self.benchmark.positions_file)
        }
    }

    pub fn lance_file_version(&self) -> Result<LanceFileVersion> {
        match self.dataset.lance_version.as_str() {
            "2.1" => Ok(LanceFileVersion::V2_1),
            "2.2" => Ok(LanceFileVersion::V2_2),
            other => {
                bail!("unsupported Lance version '{other}'; sandbox supports only 2.1 and 2.2")
            }
        }
    }

    pub fn enabled_structs(&self) -> impl Iterator<Item = (&String, &StructPackConfig)> {
        self.structs.iter().filter(|(_, group)| group.enabled)
    }

    pub fn child_to_struct(&self) -> HashMap<String, String> {
        let mut map = HashMap::new();
        for (name, group) in self.enabled_structs() {
            for child in &group.fields {
                map.insert(child.clone(), name.clone());
            }
        }
        map
    }

    pub fn field_encoding(&self, field_name: &str) -> EncodingConfig {
        let mut encoding = self.defaults.encoding.clone();
        if let Some(field) = self.fields.get(field_name) {
            field.encoding.apply_to(&mut encoding);
        }
        encoding
    }

    fn validate(&self) -> Result<()> {
        self.lance_file_version()?;
        if self.key.mode == KeyMode::Position && self.key.data_type != KeyDataType::Uint32 {
            bail!("key.mode=position requires key.type=uint32");
        }
        if self.key.mode == KeyMode::PositionKey && self.key.data_type != KeyDataType::Uint64 {
            bail!("key.mode=position_key requires key.type=uint64");
        }
        if self.benchmark.lookup_batch_sizes.is_empty() {
            bail!("benchmark.lookup_batch_sizes must not be empty");
        }
        if self.enabled_structs().next().is_some() && self.dataset.lance_version != "2.2" {
            bail!("struct packing requires dataset.lance_version=\"2.2\"");
        }
        let mut child_owner = HashMap::<&str, &str>::new();
        for (name, group) in self.enabled_structs() {
            for child in &group.fields {
                if let Some(existing) = child_owner.insert(child.as_str(), name.as_str()) {
                    bail!("field '{child}' appears in both struct '{existing}' and '{name}'");
                }
            }
        }
        Ok(())
    }
}

impl EncodingOverride {
    fn apply_to(&self, encoding: &mut EncodingConfig) {
        if let Some(value) = &self.structural {
            encoding.structural.clone_from(value);
        }
        if let Some(value) = &self.compression {
            encoding.compression.clone_from(value);
        }
        if let Some(value) = self.compression_level {
            encoding.compression_level = value;
        }
        if let Some(value) = &self.dict_values_compression {
            encoding.dict_values_compression.clone_from(value);
        }
        if let Some(value) = self.dict_values_compression_level {
            encoding.dict_values_compression_level = value;
        }
        if let Some(value) = self.rle_threshold {
            encoding.rle_threshold = value;
        }
        if let Some(value) = self.dict_size_ratio {
            encoding.dict_size_ratio = value;
        }
        if let Some(value) = self.dict_divisor {
            encoding.dict_divisor = value;
        }
        if let Some(value) = self.minichunk_size {
            encoding.minichunk_size = value;
        }
    }
}

impl EncodingConfig {
    pub fn to_lance_metadata(&self) -> HashMap<String, String> {
        HashMap::from([
            (
                "lance-encoding:structural-encoding".to_string(),
                self.structural.clone(),
            ),
            (
                "lance-encoding:compression".to_string(),
                self.compression.clone(),
            ),
            (
                "lance-encoding:compression-level".to_string(),
                self.compression_level.to_string(),
            ),
            (
                "lance-encoding:dict-values-compression".to_string(),
                self.dict_values_compression.clone(),
            ),
            (
                "lance-encoding:dict-values-compression-level".to_string(),
                self.dict_values_compression_level.to_string(),
            ),
            (
                "lance-encoding:rle-threshold".to_string(),
                self.rle_threshold.to_string(),
            ),
            (
                "lance-encoding:dict-size-ratio".to_string(),
                self.dict_size_ratio.to_string(),
            ),
            (
                "lance-encoding:dict-divisor".to_string(),
                self.dict_divisor.to_string(),
            ),
            (
                "lance-encoding:minichunk-size".to_string(),
                self.minichunk_size.to_string(),
            ),
        ])
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::{fs, path::Path};

    fn config_toml(dataset_extra: &str, key_extra: &str, benchmark_extra: &str) -> String {
        format!(
            r#"
[dataset]
name = "test"
cache_root = "/tmp/cache"
chrom = "chr1"
output_root = "/tmp/out"
lance_version = "2.1"
batch_size = 65536
overwrite = true
{dataset_extra}

[key]
mode = "position"
column = "position"
type = "uint32"
{key_extra}

[indexes]
position_btree = true
tier_bitmap = true

[benchmark]
positions_file = "inputs/positions.txt"
vcf_path = "/tmp/input.vcf.gz"
lookup_batch_sizes = [238, 512, 1024]
sort_position_keys = true
{benchmark_extra}

[inspect]
include_physical_columns = true
include_index_sizes = true

[defaults.encoding]
structural = "miniblock"
compression = "zstd"
compression_level = 3
dict_values_compression = "zstd"
dict_values_compression_level = 3
rle_threshold = 0.95
dict_size_ratio = 0.99
dict_divisor = 1
minichunk_size = 4096
"#
        )
    }

    #[test]
    fn parses_current_config_shape() {
        let config: SandboxConfig = toml::from_str(&config_toml("", "", "")).unwrap();
        config.validate().unwrap();
    }

    #[test]
    fn parses_heed_payload_sidecar_backends() {
        let raw = format!(
            "{}\n[sidecar]\npayload_backends = [\"heed_payload_raw\", \"heed_payload_zstd\"]\nheed_map_size_gib = 64\n",
            config_toml("", "", "")
        );
        let config: SandboxConfig = toml::from_str(&raw).unwrap();

        assert_eq!(config.sidecar.payload_backends.len(), 2);
        assert_eq!(config.sidecar.heed_map_size_gib, 64);
    }

    #[test]
    fn checked_in_configs_parse() {
        let configs_dir = Path::new(env!("CARGO_MANIFEST_DIR")).join("../../configs");
        let mut checked = 0;
        for entry in fs::read_dir(&configs_dir).unwrap() {
            let path = entry.unwrap().path();
            if path.extension().and_then(|ext| ext.to_str()) != Some("toml") {
                continue;
            }
            SandboxConfig::read(&path)
                .unwrap_or_else(|err| panic!("failed to parse {}: {err:#}", path.display()));
            checked += 1;
        }
        assert!(
            checked > 0,
            "no TOML configs found in {}",
            configs_dir.display()
        );
    }

    #[test]
    fn rejects_removed_dataset_knobs() {
        let err =
            toml::from_str::<SandboxConfig>(&config_toml("cold_fragment_rows = 1024", "", ""))
                .unwrap_err()
                .to_string();
        assert!(err.contains("cold_fragment_rows"), "{err}");

        let err =
            toml::from_str::<SandboxConfig>(&config_toml("warm_fragment_rows = 500000", "", ""))
                .unwrap_err()
                .to_string();
        assert!(err.contains("warm_fragment_rows"), "{err}");
    }

    #[test]
    fn rejects_removed_key_source() {
        let err = toml::from_str::<SandboxConfig>(&config_toml("", "source = \"start\"", ""))
            .unwrap_err()
            .to_string();
        assert!(err.contains("source"), "{err}");
    }

    #[test]
    fn rejects_removed_benchmark_projection() {
        let err =
            toml::from_str::<SandboxConfig>(&config_toml("", "", "projection = \"everything\""))
                .unwrap_err()
                .to_string();
        assert!(err.contains("projection"), "{err}");
    }

    #[test]
    fn rejects_removed_scanner_batch_size_knobs() {
        let err = toml::from_str::<SandboxConfig>(&config_toml("", "", "scan_batch_size = 2000"))
            .unwrap_err()
            .to_string();
        assert!(err.contains("scan_batch_size"), "{err}");

        let err =
            toml::from_str::<SandboxConfig>(&config_toml("", "", "warm_scan_batch_size = 16384"))
                .unwrap_err()
                .to_string();
        assert!(err.contains("warm_scan_batch_size"), "{err}");
    }
}
