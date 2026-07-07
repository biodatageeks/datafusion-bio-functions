#![allow(dead_code)]

#[cfg(test)]
mod tests {
    use super::*;

    fn sample_row(position: u32, variation_name: &str) -> EverythingRow {
        EverythingRow {
            position,
            variant_keys: Some(vec![10, 20]),
            chrom: "1".to_string(),
            start: position as i64,
            end: position as i64 + 1,
            allele_string: "A/T".to_string(),
            variation_name: variation_name.to_string(),
            failed: None,
            somatic: Some(1),
            strand: Some(1),
            minor_allele: Some("T".to_string()),
            minor_allele_freq: Some(0.125),
            clin_sig: Some("benign".to_string()),
            phenotype_or_disease: Some(1),
            clinical_impact: Some("LOW".to_string()),
            clin_sig_allele: Some("T".to_string()),
            af: Some("0.1".to_string()),
            afr: Some("0.2".to_string()),
            amr: None,
            eas: Some("0.3".to_string()),
            eur: Some("0.4".to_string()),
            sas: Some("0.5".to_string()),
            gnomadg: Some("0.6".to_string()),
            gnomadg_afr: None,
            gnomadg_ami: Some("0.7".to_string()),
            gnomadg_amr: Some("0.8".to_string()),
            gnomadg_asj: None,
            gnomadg_eas: Some("0.9".to_string()),
            gnomadg_fin: Some("0.01".to_string()),
            gnomadg_mid: Some("0.02".to_string()),
            gnomadg_nfe: Some("0.03".to_string()),
            gnomadg_sas: Some("0.04".to_string()),
            gnomadg_remaining: Some("0.05".to_string()),
            gnomade: Some("0.11".to_string()),
            gnomade_afr: Some("0.12".to_string()),
            gnomade_amr: Some("0.13".to_string()),
            gnomade_asj: Some("0.14".to_string()),
            gnomade_eas: None,
            gnomade_fin: Some("0.15".to_string()),
            gnomade_nfe: Some("0.16".to_string()),
            gnomade_sas: Some("0.17".to_string()),
            gnomade_mid: Some("0.18".to_string()),
            gnomade_remaining: Some("0.19".to_string()),
            clinvar_ids: Some("RCV0001".to_string()),
            cosmic_ids: None,
            dbsnp_ids: Some("rs1".to_string()),
            pubmed: Some("123".to_string()),
        }
    }

    #[test]
    fn raw_payload_codec_round_trips_all_everything_fields() {
        let payload = EverythingPayload {
            rows: vec![sample_row(100, "var-a"), sample_row(100, "var-b")],
        };

        let encoded = encode_payload_raw(&payload).unwrap();
        let decoded = decode_payload_raw(&encoded).unwrap();

        assert_eq!(decoded, payload);
    }

    #[test]
    fn zstd_payload_codec_round_trips_all_everything_fields() {
        let payload = EverythingPayload {
            rows: vec![sample_row(100, "var-a"), sample_row(100, "var-b")],
        };

        let encoded = encode_payload_zstd(&payload, 3).unwrap();
        let decoded = decode_payload_zstd(&encoded).unwrap();

        assert_eq!(decoded, payload);
    }

    #[test]
    fn position_keys_use_big_endian_order_for_lmdb_range_locality() {
        let small = encode_position_key(2);
        let large = encode_position_key(10);

        assert!(small < large);
        assert_eq!(decode_position_key(&small).unwrap(), 2);
        assert_eq!(decode_position_key(&large).unwrap(), 10);
    }
}
use std::borrow::Cow;
use std::fs;
use std::marker::PhantomData;
use std::path::{Path, PathBuf};
use std::time::Instant;

use anyhow::{Context, Result, anyhow, bail};
use arrow_array::{
    Array, Float64Array, Int8Array, Int64Array, ListArray, RecordBatch, StringArray, StructArray,
    UInt32Array, UInt64Array,
};
use arrow_schema::DataType;
use futures::TryStreamExt;
use heed::types::{Bytes, Str};
use heed::{BoxedError, BytesDecode, BytesEncode, Database, Env, EnvOpenOptions};
use lance::Dataset;
use lance::dataset::scanner::MaterializationStyle;
use serde::{Deserialize, Serialize};

const PAYLOAD_MAGIC: &[u8; 8] = b"LPAYLD01";
const LMDB_DATABASE_NAME: &str = "payloads";
const LMDB_METADATA_DATABASE_NAME: &str = "metadata";
const METADATA_ROWS_KEY: &str = "rows";
const WRITE_COMMIT_POSITIONS: usize = 50_000;
const ZSTD_LEVEL: i32 = 3;

#[derive(Debug, Clone, Copy, PartialEq, Eq, Deserialize, Serialize)]
pub enum HeedPayloadBackend {
    #[serde(rename = "heed_payload_raw")]
    Raw,
    #[serde(rename = "heed_payload_zstd")]
    Zstd,
}

impl HeedPayloadBackend {
    pub fn name(self) -> &'static str {
        match self {
            Self::Raw => "heed_payload_raw",
            Self::Zstd => "heed_payload_zstd",
        }
    }

    pub fn directory_name(self) -> &'static str {
        match self {
            Self::Raw => "payload_heed_raw.lmdb",
            Self::Zstd => "payload_heed_zstd.lmdb",
        }
    }
}

#[derive(Debug, Clone, Serialize)]
pub struct HeedPayloadSidecarReport {
    pub path: String,
    pub backend: HeedPayloadBackend,
    pub built: bool,
    pub build_seconds: f64,
    pub unique_positions: usize,
    pub rows: usize,
    pub size_bytes: u64,
}

#[derive(Debug, Clone, Copy)]
struct PayloadBuildStats {
    unique_positions: usize,
    rows: usize,
}

#[derive(Debug, Clone, PartialEq)]
pub struct EverythingPayload {
    pub rows: Vec<EverythingRow>,
}

#[derive(Debug, Clone, PartialEq)]
pub struct EverythingRow {
    pub position: u32,
    pub variant_keys: Option<Vec<u64>>,
    pub chrom: String,
    pub start: i64,
    pub end: i64,
    pub allele_string: String,
    pub variation_name: String,
    pub failed: Option<i8>,
    pub somatic: Option<i8>,
    pub strand: Option<i8>,
    pub minor_allele: Option<String>,
    pub minor_allele_freq: Option<f64>,
    pub clin_sig: Option<String>,
    pub phenotype_or_disease: Option<i8>,
    pub clinical_impact: Option<String>,
    pub clin_sig_allele: Option<String>,
    pub af: Option<String>,
    pub afr: Option<String>,
    pub amr: Option<String>,
    pub eas: Option<String>,
    pub eur: Option<String>,
    pub sas: Option<String>,
    pub gnomadg: Option<String>,
    pub gnomadg_afr: Option<String>,
    pub gnomadg_ami: Option<String>,
    pub gnomadg_amr: Option<String>,
    pub gnomadg_asj: Option<String>,
    pub gnomadg_eas: Option<String>,
    pub gnomadg_fin: Option<String>,
    pub gnomadg_mid: Option<String>,
    pub gnomadg_nfe: Option<String>,
    pub gnomadg_sas: Option<String>,
    pub gnomadg_remaining: Option<String>,
    pub gnomade: Option<String>,
    pub gnomade_afr: Option<String>,
    pub gnomade_amr: Option<String>,
    pub gnomade_asj: Option<String>,
    pub gnomade_eas: Option<String>,
    pub gnomade_fin: Option<String>,
    pub gnomade_nfe: Option<String>,
    pub gnomade_sas: Option<String>,
    pub gnomade_mid: Option<String>,
    pub gnomade_remaining: Option<String>,
    pub clinvar_ids: Option<String>,
    pub cosmic_ids: Option<String>,
    pub dbsnp_ids: Option<String>,
    pub pubmed: Option<String>,
}

#[derive(Debug, Clone, PartialEq)]
pub struct ResolvedPayloadRows {
    pub requested_positions: usize,
    pub matched_positions: usize,
    pub rows: Vec<EverythingRow>,
}

pub struct U32KeyCodec;
pub struct RawEverythingPayloadCodec;
pub struct ZstdEverythingPayloadCodec;

pub struct HeedPayloadIndex<C> {
    env: Env,
    path: PathBuf,
    backend: HeedPayloadBackend,
    unique_positions: usize,
    rows: usize,
    _codec: PhantomData<C>,
}

pub enum LoadedHeedPayloadIndex {
    Raw(HeedPayloadIndex<RawEverythingPayloadCodec>),
    Zstd(HeedPayloadIndex<ZstdEverythingPayloadCodec>),
}

impl LoadedHeedPayloadIndex {
    pub fn backend(&self) -> HeedPayloadBackend {
        match self {
            Self::Raw(index) => index.backend,
            Self::Zstd(index) => index.backend,
        }
    }

    pub fn unique_positions(&self) -> usize {
        match self {
            Self::Raw(index) => index.unique_positions,
            Self::Zstd(index) => index.unique_positions,
        }
    }

    pub fn rows(&self) -> usize {
        match self {
            Self::Raw(index) => index.rows,
            Self::Zstd(index) => index.rows,
        }
    }

    pub fn path(&self) -> &Path {
        match self {
            Self::Raw(index) => &index.path,
            Self::Zstd(index) => &index.path,
        }
    }

    pub fn resolve_sorted_positions(&self, positions: &[u64]) -> Result<ResolvedPayloadRows> {
        match self {
            Self::Raw(index) => index.resolve_sorted_positions(positions),
            Self::Zstd(index) => index.resolve_sorted_positions(positions),
        }
    }
}

impl<'a> BytesEncode<'a> for U32KeyCodec {
    type EItem = u32;

    fn bytes_encode(value: &'a Self::EItem) -> Result<Cow<'a, [u8]>, BoxedError> {
        Ok(Cow::Owned(encode_position_key(*value).to_vec()))
    }
}

impl<'a> BytesDecode<'a> for U32KeyCodec {
    type DItem = u32;

    fn bytes_decode(bytes: &'a [u8]) -> Result<Self::DItem, BoxedError> {
        decode_position_key(bytes).map_err(Into::into)
    }
}

impl<'a> BytesEncode<'a> for RawEverythingPayloadCodec {
    type EItem = EverythingPayload;

    fn bytes_encode(value: &'a Self::EItem) -> Result<Cow<'a, [u8]>, BoxedError> {
        encode_payload_raw(value)
            .map(Cow::Owned)
            .map_err(Into::into)
    }
}

impl<'a> BytesDecode<'a> for RawEverythingPayloadCodec {
    type DItem = EverythingPayload;

    fn bytes_decode(bytes: &'a [u8]) -> Result<Self::DItem, BoxedError> {
        decode_payload_raw(bytes).map_err(Into::into)
    }
}

impl<'a> BytesEncode<'a> for ZstdEverythingPayloadCodec {
    type EItem = EverythingPayload;

    fn bytes_encode(value: &'a Self::EItem) -> Result<Cow<'a, [u8]>, BoxedError> {
        encode_payload_zstd(value, ZSTD_LEVEL)
            .map(Cow::Owned)
            .map_err(Into::into)
    }
}

impl<'a> BytesDecode<'a> for ZstdEverythingPayloadCodec {
    type DItem = EverythingPayload;

    fn bytes_decode(bytes: &'a [u8]) -> Result<Self::DItem, BoxedError> {
        decode_payload_zstd(bytes).map_err(Into::into)
    }
}

impl<C> HeedPayloadIndex<C>
where
    for<'a> C: BytesDecode<'a, DItem = EverythingPayload>,
    C: 'static,
{
    pub fn resolve_sorted_positions(&self, positions: &[u64]) -> Result<ResolvedPayloadRows> {
        let txn = self
            .env
            .read_txn()
            .context("failed to open heed payload read transaction")?;
        let db: Database<U32KeyCodec, C> = self
            .env
            .open_database(&txn, Some(LMDB_DATABASE_NAME))?
            .context("heed payload database is missing")?;
        let mut matched_positions = 0usize;
        let mut rows = Vec::new();

        for &position in positions {
            let Ok(position) = u32::try_from(position) else {
                continue;
            };
            if let Some(payload) = db
                .get(&txn, &position)
                .with_context(|| format!("failed to read heed payload for position {position}"))?
            {
                matched_positions += 1;
                rows.extend(payload.rows);
            }
        }

        Ok(ResolvedPayloadRows {
            requested_positions: positions.len(),
            matched_positions,
            rows,
        })
    }
}

pub async fn load_or_build_heed_payload_index(
    dataset_path: &Path,
    dataset: &Dataset,
    projection: &[String],
    backend: HeedPayloadBackend,
    map_size_gib: usize,
) -> Result<(LoadedHeedPayloadIndex, HeedPayloadSidecarReport)> {
    let path = payload_sidecar_path(dataset_path, backend);
    let started = Instant::now();

    if path.exists() {
        match open_existing_payload_index(&path, backend, map_size_gib) {
            Ok(index) => {
                let report = HeedPayloadSidecarReport {
                    path: path.display().to_string(),
                    backend,
                    built: false,
                    build_seconds: 0.0,
                    unique_positions: index.unique_positions(),
                    rows: index.rows(),
                    size_bytes: directory_size(&path)?,
                };
                return Ok((index, report));
            }
            Err(error) => {
                eprintln!(
                    "Ignoring unreadable heed payload sidecar '{}': {error:#}",
                    path.display()
                );
                fs::remove_dir_all(&path).with_context(|| {
                    format!("failed to remove stale heed sidecar '{}'", path.display())
                })?;
            }
        }
    }

    let temp_path = path.with_extension("lmdb.tmp");
    if temp_path.exists() {
        fs::remove_dir_all(&temp_path).with_context(|| {
            format!(
                "failed to remove stale temp heed sidecar '{}'",
                temp_path.display()
            )
        })?;
    }
    let stats =
        match build_payload_index(&temp_path, dataset, projection, backend, map_size_gib).await {
            Ok(stats) => stats,
            Err(error) => {
                let _ = fs::remove_dir_all(&temp_path);
                return Err(error);
            }
        };
    fs::rename(&temp_path, &path).with_context(|| {
        format!(
            "failed to move heed sidecar '{}' to '{}'",
            temp_path.display(),
            path.display()
        )
    })?;
    let index = open_payload_index_with_counts(&path, backend, map_size_gib, stats)?;
    let report = HeedPayloadSidecarReport {
        path: path.display().to_string(),
        backend,
        built: true,
        build_seconds: started.elapsed().as_secs_f64(),
        unique_positions: index.unique_positions(),
        rows: index.rows(),
        size_bytes: directory_size(&path)?,
    };
    Ok((index, report))
}

pub fn payload_sidecar_path(dataset_path: &Path, backend: HeedPayloadBackend) -> PathBuf {
    dataset_path
        .parent()
        .unwrap_or_else(|| Path::new("."))
        .join("reports")
        .join(backend.directory_name())
}

pub fn encode_position_key(position: u32) -> [u8; 4] {
    position.to_be_bytes()
}

pub fn decode_position_key(bytes: &[u8]) -> Result<u32> {
    let bytes: [u8; 4] = bytes
        .try_into()
        .map_err(|_| anyhow!("position key must be exactly 4 bytes"))?;
    Ok(u32::from_be_bytes(bytes))
}

pub fn encode_payload_zstd(payload: &EverythingPayload, level: i32) -> Result<Vec<u8>> {
    let raw = encode_payload_raw(payload)?;
    zstd::stream::encode_all(raw.as_slice(), level).context("failed to zstd-compress payload")
}

pub fn decode_payload_zstd(bytes: &[u8]) -> Result<EverythingPayload> {
    let raw = zstd::stream::decode_all(bytes).context("failed to zstd-decompress payload")?;
    decode_payload_raw(&raw)
}

pub fn encode_payload_raw(payload: &EverythingPayload) -> Result<Vec<u8>> {
    let mut out = Vec::new();
    out.extend_from_slice(PAYLOAD_MAGIC);
    write_len(&mut out, payload.rows.len())?;
    for row in &payload.rows {
        encode_row(&mut out, row)?;
    }
    Ok(out)
}

pub fn decode_payload_raw(bytes: &[u8]) -> Result<EverythingPayload> {
    let mut reader = ByteReader::new(bytes);
    let magic = reader.read_exact(PAYLOAD_MAGIC.len())?;
    if magic != PAYLOAD_MAGIC {
        bail!("invalid payload magic");
    }
    let rows_len = reader.read_u32()? as usize;
    let mut rows = Vec::with_capacity(rows_len);
    for _ in 0..rows_len {
        rows.push(decode_row(&mut reader)?);
    }
    if !reader.is_empty() {
        bail!("payload has {} trailing bytes", reader.remaining());
    }
    Ok(EverythingPayload { rows })
}

async fn build_payload_index(
    path: &Path,
    dataset: &Dataset,
    projection: &[String],
    backend: HeedPayloadBackend,
    map_size_gib: usize,
) -> Result<PayloadBuildStats> {
    if path.exists() {
        fs::remove_dir_all(path)
            .with_context(|| format!("failed to remove heed sidecar '{}'", path.display()))?;
    }
    fs::create_dir_all(path)
        .with_context(|| format!("failed to create heed sidecar '{}'", path.display()))?;

    match backend {
        HeedPayloadBackend::Raw => {
            create_and_fill_index::<RawEverythingPayloadCodec>(
                path,
                dataset,
                projection,
                backend,
                map_size_gib,
            )
            .await
        }
        HeedPayloadBackend::Zstd => {
            create_and_fill_index::<ZstdEverythingPayloadCodec>(
                path,
                dataset,
                projection,
                backend,
                map_size_gib,
            )
            .await
        }
    }
}

async fn create_and_fill_index<C>(
    path: &Path,
    dataset: &Dataset,
    projection: &[String],
    backend: HeedPayloadBackend,
    map_size_gib: usize,
) -> Result<PayloadBuildStats>
where
    for<'a> C:
        BytesEncode<'a, EItem = EverythingPayload> + BytesDecode<'a, DItem = EverythingPayload>,
    C: 'static,
{
    let started = Instant::now();
    let env = open_env(path, map_size_gib)?;
    let mut wtxn = env.write_txn()?;
    let db: Database<U32KeyCodec, C> = env.create_database(&mut wtxn, Some(LMDB_DATABASE_NAME))?;
    let metadata: Database<Str, Bytes> =
        env.create_database(&mut wtxn, Some(LMDB_METADATA_DATABASE_NAME))?;

    let projection_refs = projection.iter().map(String::as_str).collect::<Vec<_>>();
    let mut scanner = dataset.scan();
    scanner
        .filter("tier = 1")?
        .project(&projection_refs)?
        .materialization_style(MaterializationStyle::AllLate);

    let mut stream = scanner.try_into_stream().await?;
    let mut current_position = None;
    let mut current_rows = Vec::new();
    let mut unique_positions = 0usize;
    let mut total_rows = 0usize;
    let mut pending_positions = 0usize;

    while let Some(batch) = stream.try_next().await? {
        for row_idx in 0..batch.num_rows() {
            let row = row_from_batch(&batch, row_idx)?;
            if let Some(position) = current_position {
                if row.position < position {
                    bail!(
                        "heed payload sidecar input is not monotonic: previous position {position}, current position {}",
                        row.position
                    );
                }
                if row.position != position {
                    put_payload(&db, &mut wtxn, position, current_rows)?;
                    unique_positions += 1;
                    pending_positions += 1;
                    current_rows = Vec::new();
                    current_position = Some(row.position);
                }
            } else {
                current_position = Some(row.position);
            }

            current_rows.push(row);
            total_rows += 1;

            if pending_positions >= WRITE_COMMIT_POSITIONS {
                wtxn.commit()?;
                eprintln!(
                    "Heed payload build progress: backend={}, unique_positions={}, rows={}, elapsed={:.1}s, temp_size={}",
                    backend.name(),
                    unique_positions,
                    total_rows,
                    started.elapsed().as_secs_f64(),
                    format_bytes(directory_size(path).unwrap_or_default())
                );
                wtxn = env.write_txn()?;
                pending_positions = 0;
            }
        }
    }

    if let Some(position) = current_position {
        put_payload(&db, &mut wtxn, position, current_rows)?;
        unique_positions += 1;
    }
    put_metadata(&metadata, &mut wtxn, total_rows)?;
    wtxn.commit()?;

    eprintln!(
        "Heed payload build complete: backend={}, unique_positions={}, rows={}, elapsed={:.1}s, temp_size={}",
        backend.name(),
        unique_positions,
        total_rows,
        started.elapsed().as_secs_f64(),
        format_bytes(directory_size(path).unwrap_or_default())
    );

    drop(env);

    Ok(PayloadBuildStats {
        unique_positions,
        rows: total_rows,
    })
}

fn open_existing_payload_index(
    path: &Path,
    backend: HeedPayloadBackend,
    map_size_gib: usize,
) -> Result<LoadedHeedPayloadIndex> {
    match backend {
        HeedPayloadBackend::Raw => {
            let index =
                open_existing_index::<RawEverythingPayloadCodec>(path, backend, map_size_gib)?;
            Ok(LoadedHeedPayloadIndex::Raw(index))
        }
        HeedPayloadBackend::Zstd => {
            let index =
                open_existing_index::<ZstdEverythingPayloadCodec>(path, backend, map_size_gib)?;
            Ok(LoadedHeedPayloadIndex::Zstd(index))
        }
    }
}

fn open_payload_index_with_counts(
    path: &Path,
    backend: HeedPayloadBackend,
    map_size_gib: usize,
    stats: PayloadBuildStats,
) -> Result<LoadedHeedPayloadIndex> {
    match backend {
        HeedPayloadBackend::Raw => {
            let index = open_index_with_counts::<RawEverythingPayloadCodec>(
                path,
                backend,
                map_size_gib,
                stats,
            )?;
            Ok(LoadedHeedPayloadIndex::Raw(index))
        }
        HeedPayloadBackend::Zstd => {
            let index = open_index_with_counts::<ZstdEverythingPayloadCodec>(
                path,
                backend,
                map_size_gib,
                stats,
            )?;
            Ok(LoadedHeedPayloadIndex::Zstd(index))
        }
    }
}

fn open_index_with_counts<C>(
    path: &Path,
    backend: HeedPayloadBackend,
    map_size_gib: usize,
    stats: PayloadBuildStats,
) -> Result<HeedPayloadIndex<C>>
where
    for<'a> C: BytesDecode<'a, DItem = EverythingPayload>,
    C: 'static,
{
    let env = open_env(path, map_size_gib)?;
    let rtxn = env.read_txn()?;
    let _db: Database<U32KeyCodec, C> = env
        .open_database(&rtxn, Some(LMDB_DATABASE_NAME))?
        .context("heed payload database is missing")?;
    drop(rtxn);

    Ok(HeedPayloadIndex {
        env,
        path: path.to_path_buf(),
        backend,
        unique_positions: stats.unique_positions,
        rows: stats.rows,
        _codec: PhantomData,
    })
}

fn open_existing_index<C>(
    path: &Path,
    backend: HeedPayloadBackend,
    map_size_gib: usize,
) -> Result<HeedPayloadIndex<C>>
where
    for<'a> C: BytesDecode<'a, DItem = EverythingPayload>,
    C: 'static,
{
    let env = open_env(path, map_size_gib)?;
    let rtxn = env.read_txn()?;
    let db: Database<U32KeyCodec, C> = env
        .open_database(&rtxn, Some(LMDB_DATABASE_NAME))?
        .context("heed payload database is missing")?;
    let unique_positions = db.stat(&rtxn)?.entries;
    let rows = read_metadata_rows(&env, &rtxn)?.unwrap_or_default();
    drop(rtxn);

    Ok(HeedPayloadIndex {
        env,
        path: path.to_path_buf(),
        backend,
        unique_positions,
        rows,
        _codec: PhantomData,
    })
}

fn put_metadata(db: &Database<Str, Bytes>, wtxn: &mut heed::RwTxn<'_>, rows: usize) -> Result<()> {
    let rows = u64::try_from(rows)
        .context("heed payload row count does not fit u64")?
        .to_le_bytes();
    db.put(wtxn, METADATA_ROWS_KEY, &rows)?;
    Ok(())
}

fn read_metadata_rows(env: &Env, rtxn: &heed::RoTxn<'_>) -> Result<Option<usize>> {
    let Some(db): Option<Database<Str, Bytes>> =
        env.open_database(rtxn, Some(LMDB_METADATA_DATABASE_NAME))?
    else {
        return Ok(None);
    };
    let Some(bytes) = db.get(rtxn, METADATA_ROWS_KEY)? else {
        return Ok(None);
    };
    let bytes: [u8; 8] = bytes
        .try_into()
        .context("heed payload metadata rows value must be 8 bytes")?;
    Ok(Some(u64::from_le_bytes(bytes) as usize))
}

fn format_bytes(bytes: u64) -> String {
    const GIB: f64 = 1024.0 * 1024.0 * 1024.0;
    const MIB: f64 = 1024.0 * 1024.0;
    if bytes as f64 >= GIB {
        format!("{:.2} GiB", bytes as f64 / GIB)
    } else {
        format!("{:.2} MiB", bytes as f64 / MIB)
    }
}

fn open_env(path: &Path, map_size_gib: usize) -> Result<Env> {
    let map_size = map_size_gib
        .checked_mul(1024 * 1024 * 1024)
        .context("heed map size overflows usize")?;
    let env = unsafe {
        EnvOpenOptions::new()
            .map_size(map_size)
            .max_dbs(4)
            .open(path)
            .with_context(|| format!("failed to open heed env '{}'", path.display()))?
    };
    Ok(env)
}

fn put_payload<C>(
    db: &Database<U32KeyCodec, C>,
    wtxn: &mut heed::RwTxn<'_>,
    position: u32,
    rows: Vec<EverythingRow>,
) -> Result<()>
where
    for<'a> C: BytesEncode<'a, EItem = EverythingPayload>,
{
    let payload = EverythingPayload { rows };
    db.put(wtxn, &position, &payload)?;
    Ok(())
}

fn directory_size(path: &Path) -> Result<u64> {
    let mut total = 0u64;
    for entry in
        fs::read_dir(path).with_context(|| format!("failed to read dir '{}'", path.display()))?
    {
        let entry = entry?;
        let metadata = entry.metadata()?;
        if metadata.is_dir() {
            total += directory_size(&entry.path())?;
        } else {
            total += metadata.len();
        }
    }
    Ok(total)
}

fn row_from_batch(batch: &RecordBatch, row: usize) -> Result<EverythingRow> {
    let lookup = BatchLookup { batch };
    Ok(EverythingRow {
        position: lookup.required_u32("position", row)?,
        variant_keys: lookup.optional_u64_list("variant_keys", row)?,
        chrom: lookup.required_string("chrom", row)?,
        start: lookup.required_i64("start", row)?,
        end: lookup.required_i64("end", row)?,
        allele_string: lookup.required_string("allele_string", row)?,
        variation_name: lookup.required_string("variation_name", row)?,
        failed: lookup.optional_i8("failed", row)?,
        somatic: lookup.optional_i8("somatic", row)?,
        strand: lookup.optional_i8("strand", row)?,
        minor_allele: lookup.optional_string("minor_allele", row)?,
        minor_allele_freq: lookup.optional_f64("minor_allele_freq", row)?,
        clin_sig: lookup.optional_string("clin_sig", row)?,
        phenotype_or_disease: lookup.optional_i8("phenotype_or_disease", row)?,
        clinical_impact: lookup.optional_string("clinical_impact", row)?,
        clin_sig_allele: lookup.optional_string("clin_sig_allele", row)?,
        af: lookup.optional_string("AF", row)?,
        afr: lookup.optional_string("AFR", row)?,
        amr: lookup.optional_string("AMR", row)?,
        eas: lookup.optional_string("EAS", row)?,
        eur: lookup.optional_string("EUR", row)?,
        sas: lookup.optional_string("SAS", row)?,
        gnomadg: lookup.optional_string("gnomADg", row)?,
        gnomadg_afr: lookup.optional_string("gnomADg_AFR", row)?,
        gnomadg_ami: lookup.optional_string("gnomADg_AMI", row)?,
        gnomadg_amr: lookup.optional_string("gnomADg_AMR", row)?,
        gnomadg_asj: lookup.optional_string("gnomADg_ASJ", row)?,
        gnomadg_eas: lookup.optional_string("gnomADg_EAS", row)?,
        gnomadg_fin: lookup.optional_string("gnomADg_FIN", row)?,
        gnomadg_mid: lookup.optional_string("gnomADg_MID", row)?,
        gnomadg_nfe: lookup.optional_string("gnomADg_NFE", row)?,
        gnomadg_sas: lookup.optional_string("gnomADg_SAS", row)?,
        gnomadg_remaining: lookup.optional_string("gnomADg_REMAINING", row)?,
        gnomade: lookup.optional_string("gnomADe", row)?,
        gnomade_afr: lookup.optional_string("gnomADe_AFR", row)?,
        gnomade_amr: lookup.optional_string("gnomADe_AMR", row)?,
        gnomade_asj: lookup.optional_string("gnomADe_ASJ", row)?,
        gnomade_eas: lookup.optional_string("gnomADe_EAS", row)?,
        gnomade_fin: lookup.optional_string("gnomADe_FIN", row)?,
        gnomade_nfe: lookup.optional_string("gnomADe_NFE", row)?,
        gnomade_sas: lookup.optional_string("gnomADe_SAS", row)?,
        gnomade_mid: lookup.optional_string("gnomADe_MID", row)?,
        gnomade_remaining: lookup.optional_string("gnomADe_REMAINING", row)?,
        clinvar_ids: lookup.optional_string("clinvar_ids", row)?,
        cosmic_ids: lookup.optional_string("cosmic_ids", row)?,
        dbsnp_ids: lookup.optional_string("dbsnp_ids", row)?,
        pubmed: lookup.optional_string("pubmed", row)?,
    })
}

struct BatchLookup<'a> {
    batch: &'a RecordBatch,
}

struct ColumnView<'a> {
    array: &'a dyn Array,
    parent: Option<&'a dyn Array>,
}

impl ColumnView<'_> {
    fn is_null(&self, row: usize) -> bool {
        self.parent.is_some_and(|parent| parent.is_null(row)) || self.array.is_null(row)
    }
}

impl<'a> BatchLookup<'a> {
    fn column(&self, name: &str) -> Result<ColumnView<'a>> {
        if let Some(idx) = self.batch.schema_ref().index_of(name).ok() {
            return Ok(ColumnView {
                array: self.batch.column(idx).as_ref(),
                parent: None,
            });
        }

        for idx in 0..self.batch.num_columns() {
            let field = self.batch.schema_ref().field(idx);
            if let DataType::Struct(fields) = field.data_type() {
                if let Some(child_idx) = fields.iter().position(|child| child.name() == name) {
                    let struct_array = self
                        .batch
                        .column(idx)
                        .as_any()
                        .downcast_ref::<StructArray>()
                        .with_context(|| {
                            format!("column '{}' is not a StructArray", field.name())
                        })?;
                    return Ok(ColumnView {
                        array: struct_array.column(child_idx).as_ref(),
                        parent: Some(struct_array),
                    });
                }
            }
        }

        bail!("missing payload field '{name}'")
    }

    fn required_string(&self, name: &str, row: usize) -> Result<String> {
        self.optional_string(name, row)?
            .with_context(|| format!("required string field '{name}' is null at row {row}"))
    }

    fn optional_string(&self, name: &str, row: usize) -> Result<Option<String>> {
        let column = self.column(name)?;
        if column.is_null(row) {
            return Ok(None);
        }
        let values = column
            .array
            .as_any()
            .downcast_ref::<StringArray>()
            .with_context(|| format!("field '{name}' must be Utf8"))?;
        Ok(Some(values.value(row).to_string()))
    }

    fn required_i64(&self, name: &str, row: usize) -> Result<i64> {
        let column = self.column(name)?;
        if column.is_null(row) {
            bail!("required i64 field '{name}' is null at row {row}");
        }
        let values = column
            .array
            .as_any()
            .downcast_ref::<Int64Array>()
            .with_context(|| format!("field '{name}' must be Int64"))?;
        Ok(values.value(row))
    }

    fn required_u32(&self, name: &str, row: usize) -> Result<u32> {
        let column = self.column(name)?;
        if column.is_null(row) {
            bail!("required u32 field '{name}' is null at row {row}");
        }
        let values = column
            .array
            .as_any()
            .downcast_ref::<UInt32Array>()
            .with_context(|| format!("field '{name}' must be UInt32"))?;
        Ok(values.value(row))
    }

    fn optional_i8(&self, name: &str, row: usize) -> Result<Option<i8>> {
        let column = self.column(name)?;
        if column.is_null(row) {
            return Ok(None);
        }
        let values = column
            .array
            .as_any()
            .downcast_ref::<Int8Array>()
            .with_context(|| format!("field '{name}' must be Int8"))?;
        Ok(Some(values.value(row)))
    }

    fn optional_f64(&self, name: &str, row: usize) -> Result<Option<f64>> {
        let column = self.column(name)?;
        if column.is_null(row) {
            return Ok(None);
        }
        let values = column
            .array
            .as_any()
            .downcast_ref::<Float64Array>()
            .with_context(|| format!("field '{name}' must be Float64"))?;
        Ok(Some(values.value(row)))
    }

    fn optional_u64_list(&self, name: &str, row: usize) -> Result<Option<Vec<u64>>> {
        let column = self.column(name)?;
        if column.is_null(row) {
            return Ok(None);
        }
        let values = column
            .array
            .as_any()
            .downcast_ref::<ListArray>()
            .with_context(|| format!("field '{name}' must be List(UInt64)"))?;
        let offsets = values.value_offsets();
        let start = offsets[row] as usize;
        let end = offsets[row + 1] as usize;
        let value_array = values
            .values()
            .as_any()
            .downcast_ref::<UInt64Array>()
            .with_context(|| format!("field '{name}' values must be UInt64"))?;
        let mut out = Vec::with_capacity(end - start);
        for idx in start..end {
            if !value_array.is_null(idx) {
                out.push(value_array.value(idx));
            }
        }
        Ok(Some(out))
    }
}

fn encode_row(out: &mut Vec<u8>, row: &EverythingRow) -> Result<()> {
    write_u32(out, row.position);
    write_opt_u64_list(out, &row.variant_keys)?;
    write_string(out, &row.chrom)?;
    write_i64(out, row.start);
    write_i64(out, row.end);
    write_string(out, &row.allele_string)?;
    write_string(out, &row.variation_name)?;
    write_opt_i8(out, row.failed);
    write_opt_i8(out, row.somatic);
    write_opt_i8(out, row.strand);
    write_opt_string(out, &row.minor_allele)?;
    write_opt_f64(out, row.minor_allele_freq);
    write_opt_string(out, &row.clin_sig)?;
    write_opt_i8(out, row.phenotype_or_disease);
    write_opt_string(out, &row.clinical_impact)?;
    write_opt_string(out, &row.clin_sig_allele)?;
    write_opt_string(out, &row.af)?;
    write_opt_string(out, &row.afr)?;
    write_opt_string(out, &row.amr)?;
    write_opt_string(out, &row.eas)?;
    write_opt_string(out, &row.eur)?;
    write_opt_string(out, &row.sas)?;
    write_opt_string(out, &row.gnomadg)?;
    write_opt_string(out, &row.gnomadg_afr)?;
    write_opt_string(out, &row.gnomadg_ami)?;
    write_opt_string(out, &row.gnomadg_amr)?;
    write_opt_string(out, &row.gnomadg_asj)?;
    write_opt_string(out, &row.gnomadg_eas)?;
    write_opt_string(out, &row.gnomadg_fin)?;
    write_opt_string(out, &row.gnomadg_mid)?;
    write_opt_string(out, &row.gnomadg_nfe)?;
    write_opt_string(out, &row.gnomadg_sas)?;
    write_opt_string(out, &row.gnomadg_remaining)?;
    write_opt_string(out, &row.gnomade)?;
    write_opt_string(out, &row.gnomade_afr)?;
    write_opt_string(out, &row.gnomade_amr)?;
    write_opt_string(out, &row.gnomade_asj)?;
    write_opt_string(out, &row.gnomade_eas)?;
    write_opt_string(out, &row.gnomade_fin)?;
    write_opt_string(out, &row.gnomade_nfe)?;
    write_opt_string(out, &row.gnomade_sas)?;
    write_opt_string(out, &row.gnomade_mid)?;
    write_opt_string(out, &row.gnomade_remaining)?;
    write_opt_string(out, &row.clinvar_ids)?;
    write_opt_string(out, &row.cosmic_ids)?;
    write_opt_string(out, &row.dbsnp_ids)?;
    write_opt_string(out, &row.pubmed)?;
    Ok(())
}

fn decode_row(reader: &mut ByteReader<'_>) -> Result<EverythingRow> {
    Ok(EverythingRow {
        position: reader.read_u32()?,
        variant_keys: reader.read_opt_u64_list()?,
        chrom: reader.read_string()?,
        start: reader.read_i64()?,
        end: reader.read_i64()?,
        allele_string: reader.read_string()?,
        variation_name: reader.read_string()?,
        failed: reader.read_opt_i8()?,
        somatic: reader.read_opt_i8()?,
        strand: reader.read_opt_i8()?,
        minor_allele: reader.read_opt_string()?,
        minor_allele_freq: reader.read_opt_f64()?,
        clin_sig: reader.read_opt_string()?,
        phenotype_or_disease: reader.read_opt_i8()?,
        clinical_impact: reader.read_opt_string()?,
        clin_sig_allele: reader.read_opt_string()?,
        af: reader.read_opt_string()?,
        afr: reader.read_opt_string()?,
        amr: reader.read_opt_string()?,
        eas: reader.read_opt_string()?,
        eur: reader.read_opt_string()?,
        sas: reader.read_opt_string()?,
        gnomadg: reader.read_opt_string()?,
        gnomadg_afr: reader.read_opt_string()?,
        gnomadg_ami: reader.read_opt_string()?,
        gnomadg_amr: reader.read_opt_string()?,
        gnomadg_asj: reader.read_opt_string()?,
        gnomadg_eas: reader.read_opt_string()?,
        gnomadg_fin: reader.read_opt_string()?,
        gnomadg_mid: reader.read_opt_string()?,
        gnomadg_nfe: reader.read_opt_string()?,
        gnomadg_sas: reader.read_opt_string()?,
        gnomadg_remaining: reader.read_opt_string()?,
        gnomade: reader.read_opt_string()?,
        gnomade_afr: reader.read_opt_string()?,
        gnomade_amr: reader.read_opt_string()?,
        gnomade_asj: reader.read_opt_string()?,
        gnomade_eas: reader.read_opt_string()?,
        gnomade_fin: reader.read_opt_string()?,
        gnomade_nfe: reader.read_opt_string()?,
        gnomade_sas: reader.read_opt_string()?,
        gnomade_mid: reader.read_opt_string()?,
        gnomade_remaining: reader.read_opt_string()?,
        clinvar_ids: reader.read_opt_string()?,
        cosmic_ids: reader.read_opt_string()?,
        dbsnp_ids: reader.read_opt_string()?,
        pubmed: reader.read_opt_string()?,
    })
}

fn write_len(out: &mut Vec<u8>, len: usize) -> Result<()> {
    let len = u32::try_from(len).context("payload length exceeds u32")?;
    write_u32(out, len);
    Ok(())
}

fn write_u32(out: &mut Vec<u8>, value: u32) {
    out.extend_from_slice(&value.to_le_bytes());
}

fn write_i64(out: &mut Vec<u8>, value: i64) {
    out.extend_from_slice(&value.to_le_bytes());
}

fn write_opt_i8(out: &mut Vec<u8>, value: Option<i8>) {
    match value {
        Some(value) => {
            out.push(1);
            out.push(value as u8);
        }
        None => out.push(0),
    }
}

fn write_opt_f64(out: &mut Vec<u8>, value: Option<f64>) {
    match value {
        Some(value) => {
            out.push(1);
            out.extend_from_slice(&value.to_le_bytes());
        }
        None => out.push(0),
    }
}

fn write_string(out: &mut Vec<u8>, value: &str) -> Result<()> {
    write_len(out, value.len())?;
    out.extend_from_slice(value.as_bytes());
    Ok(())
}

fn write_opt_string(out: &mut Vec<u8>, value: &Option<String>) -> Result<()> {
    match value {
        Some(value) => {
            out.push(1);
            write_string(out, value)?;
        }
        None => out.push(0),
    }
    Ok(())
}

fn write_opt_u64_list(out: &mut Vec<u8>, value: &Option<Vec<u64>>) -> Result<()> {
    match value {
        Some(values) => {
            out.push(1);
            write_len(out, values.len())?;
            for value in values {
                out.extend_from_slice(&value.to_le_bytes());
            }
        }
        None => out.push(0),
    }
    Ok(())
}

struct ByteReader<'a> {
    bytes: &'a [u8],
    pos: usize,
}

impl<'a> ByteReader<'a> {
    fn new(bytes: &'a [u8]) -> Self {
        Self { bytes, pos: 0 }
    }

    fn remaining(&self) -> usize {
        self.bytes.len() - self.pos
    }

    fn is_empty(&self) -> bool {
        self.pos == self.bytes.len()
    }

    fn read_exact(&mut self, len: usize) -> Result<&'a [u8]> {
        let end = self
            .pos
            .checked_add(len)
            .context("payload cursor overflow")?;
        if end > self.bytes.len() {
            bail!("payload ended early");
        }
        let out = &self.bytes[self.pos..end];
        self.pos = end;
        Ok(out)
    }

    fn read_u8(&mut self) -> Result<u8> {
        Ok(self.read_exact(1)?[0])
    }

    fn read_u32(&mut self) -> Result<u32> {
        let bytes: [u8; 4] = self.read_exact(4)?.try_into().unwrap();
        Ok(u32::from_le_bytes(bytes))
    }

    fn read_i64(&mut self) -> Result<i64> {
        let bytes: [u8; 8] = self.read_exact(8)?.try_into().unwrap();
        Ok(i64::from_le_bytes(bytes))
    }

    fn read_string(&mut self) -> Result<String> {
        let len = self.read_u32()? as usize;
        let bytes = self.read_exact(len)?;
        String::from_utf8(bytes.to_vec()).context("payload string is not valid UTF-8")
    }

    fn read_opt_i8(&mut self) -> Result<Option<i8>> {
        match self.read_u8()? {
            0 => Ok(None),
            1 => Ok(Some(self.read_u8()? as i8)),
            tag => bail!("invalid optional i8 tag {tag}"),
        }
    }

    fn read_opt_f64(&mut self) -> Result<Option<f64>> {
        match self.read_u8()? {
            0 => Ok(None),
            1 => {
                let bytes: [u8; 8] = self.read_exact(8)?.try_into().unwrap();
                Ok(Some(f64::from_le_bytes(bytes)))
            }
            tag => bail!("invalid optional f64 tag {tag}"),
        }
    }

    fn read_opt_string(&mut self) -> Result<Option<String>> {
        match self.read_u8()? {
            0 => Ok(None),
            1 => Ok(Some(self.read_string()?)),
            tag => bail!("invalid optional string tag {tag}"),
        }
    }

    fn read_opt_u64_list(&mut self) -> Result<Option<Vec<u64>>> {
        match self.read_u8()? {
            0 => Ok(None),
            1 => {
                let len = self.read_u32()? as usize;
                let mut values = Vec::with_capacity(len);
                for _ in 0..len {
                    let bytes: [u8; 8] = self.read_exact(8)?.try_into().unwrap();
                    values.push(u64::from_le_bytes(bytes));
                }
                Ok(Some(values))
            }
            tag => bail!("invalid optional u64 list tag {tag}"),
        }
    }
}
