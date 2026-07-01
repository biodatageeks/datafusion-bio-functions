use std::io::{Read, Write};
use std::mem::size_of;
use std::slice;

use anyhow::{Context, Result, bail};

const MAGIC: &[u8; 8] = b"LRSIDX02";

#[cfg(test)]
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct PositionRowIdPair {
    pub position: u32,
    pub row_id: u64,
}

#[derive(Debug, Clone, PartialEq, Eq)]
pub struct PositionRowIdIndex {
    positions: Vec<u32>,
    row_ids: Vec<u64>,
    unique_positions: usize,
}

#[derive(Debug, Default)]
pub struct PositionRowIdIndexBuilder {
    positions: Vec<u32>,
    row_ids: Vec<u64>,
    unique_positions: usize,
    last_position: Option<u32>,
}

#[derive(Debug, Clone, PartialEq, Eq)]
pub struct ResolvedRowIds {
    pub requested_positions: usize,
    pub matched_positions: usize,
    pub row_ids: Vec<u64>,
}

impl PositionRowIdIndex {
    #[cfg(test)]
    pub fn from_pairs(mut pairs: Vec<PositionRowIdPair>) -> Self {
        pairs.sort_unstable_by_key(|pair| (pair.position, pair.row_id));

        let mut builder = PositionRowIdIndexBuilder::with_capacity(pairs.len());
        for pair in pairs {
            builder
                .push_monotonic(pair.position, pair.row_id)
                .expect("sorted pairs must be monotonic");
        }
        builder
            .finish()
            .expect("sorted pairs must build a valid index")
    }

    #[cfg(test)]
    pub fn positions(&self) -> &[u32] {
        &self.positions
    }

    pub fn unique_positions(&self) -> usize {
        self.unique_positions
    }

    pub fn row_ids_len(&self) -> usize {
        self.row_ids.len()
    }

    pub fn pair_set_mismatch(&self, other: &Self) -> Option<String> {
        let mut left = 0usize;
        let mut right = 0usize;

        loop {
            match (self.positions.get(left), other.positions.get(right)) {
                (None, None) => return None,
                (Some(&position), None) => {
                    return Some(format!(
                        "right side ended before position {position}; left_remaining_rows={}",
                        self.positions.len() - left
                    ));
                }
                (None, Some(&position)) => {
                    return Some(format!(
                        "left side ended before position {position}; right_remaining_rows={}",
                        other.positions.len() - right
                    ));
                }
                (Some(&left_position), Some(&right_position))
                    if left_position != right_position =>
                {
                    return Some(format!(
                        "position mismatch: left position {left_position}, right position {right_position}"
                    ));
                }
                (Some(&position), Some(_)) => {
                    let left_start = left;
                    while left < self.positions.len() && self.positions[left] == position {
                        left += 1;
                    }
                    let right_start = right;
                    while right < other.positions.len() && other.positions[right] == position {
                        right += 1;
                    }

                    let mut left_row_ids = self.row_ids[left_start..left].to_vec();
                    let mut right_row_ids = other.row_ids[right_start..right].to_vec();
                    left_row_ids.sort_unstable();
                    right_row_ids.sort_unstable();

                    if left_row_ids != right_row_ids {
                        return Some(format!(
                            "position {position} differs: left_row_ids={}, right_row_ids={}",
                            format_row_id_sample(&left_row_ids),
                            format_row_id_sample(&right_row_ids)
                        ));
                    }
                }
            }
        }
    }

    pub fn resolve_sorted_positions(&self, query_positions: &[u64]) -> ResolvedRowIds {
        let mut cursor = 0usize;
        self.resolve_sorted_positions_from_cursor(query_positions, &mut cursor)
    }

    pub fn resolve_sorted_positions_from_cursor(
        &self,
        query_positions: &[u64],
        cursor: &mut usize,
    ) -> ResolvedRowIds {
        let mut matched_positions = 0usize;
        let mut row_ids = Vec::new();
        *cursor = (*cursor).min(self.positions.len());

        for &query_position in query_positions {
            let Ok(position) = u32::try_from(query_position) else {
                continue;
            };

            while *cursor < self.positions.len() && self.positions[*cursor] < position {
                *cursor += 1;
            }
            if *cursor == self.positions.len() {
                break;
            }
            if self.positions[*cursor] == position {
                matched_positions += 1;
                while *cursor < self.positions.len() && self.positions[*cursor] == position {
                    row_ids.push(self.row_ids[*cursor]);
                    *cursor += 1;
                }
            }
        }

        ResolvedRowIds {
            requested_positions: query_positions.len(),
            matched_positions,
            row_ids,
        }
    }

    pub fn write_to<W: Write>(&self, mut writer: W) -> Result<()> {
        writer.write_all(MAGIC)?;
        write_u64(&mut writer, self.unique_positions as u64)?;
        write_u64(&mut writer, self.row_ids.len() as u64)?;
        write_u32_slice(&mut writer, &self.positions)?;
        write_u64_slice(&mut writer, &self.row_ids)?;
        Ok(())
    }

    pub fn read_from<R: Read>(mut reader: R) -> Result<Self> {
        let mut magic = [0_u8; 8];
        reader.read_exact(&mut magic)?;
        if &magic != MAGIC {
            bail!("invalid row sidecar magic");
        }

        let unique_positions = read_u64(&mut reader)? as usize;
        let row_id_count = read_u64(&mut reader)? as usize;

        let positions = read_u32_vec(&mut reader, row_id_count)?;
        let row_ids = read_u64_vec(&mut reader, row_id_count)?;

        validate_index(&positions, &row_ids, unique_positions)
            .context("invalid row sidecar contents")?;

        Ok(Self {
            positions,
            row_ids,
            unique_positions,
        })
    }
}

fn format_row_id_sample(row_ids: &[u64]) -> String {
    const LIMIT: usize = 8;
    if row_ids.len() <= LIMIT {
        return format!("{row_ids:?}");
    }

    let sample = row_ids
        .iter()
        .take(LIMIT)
        .map(u64::to_string)
        .collect::<Vec<_>>()
        .join(", ");
    format!("[{sample}, ...; len={}]", row_ids.len())
}

impl PositionRowIdIndexBuilder {
    pub fn new() -> Self {
        Self::default()
    }

    #[cfg(test)]
    pub fn with_capacity(row_id_count: usize) -> Self {
        Self {
            positions: Vec::with_capacity(row_id_count),
            row_ids: Vec::with_capacity(row_id_count),
            ..Self::default()
        }
    }

    pub fn reserve(&mut self, additional: usize) {
        self.positions.reserve(additional);
        self.row_ids.reserve(additional);
    }

    pub fn push_monotonic(&mut self, position: u32, row_id: u64) -> Result<()> {
        if let Some(last_position) = self.last_position {
            if position < last_position {
                bail!(
                    "position sidecar input is not monotonic: previous position {last_position}, current position {position}"
                );
            }
        }

        if self.last_position != Some(position) {
            self.unique_positions += 1;
            self.last_position = Some(position);
        }
        self.positions.push(position);
        self.row_ids.push(row_id);
        Ok(())
    }

    pub fn finish(self) -> Result<PositionRowIdIndex> {
        validate_index(&self.positions, &self.row_ids, self.unique_positions)
            .context("invalid row sidecar builder contents")?;
        Ok(PositionRowIdIndex {
            positions: self.positions,
            row_ids: self.row_ids,
            unique_positions: self.unique_positions,
        })
    }
}

fn validate_index(positions: &[u32], row_ids: &[u64], unique_positions: usize) -> Result<()> {
    if positions.len() != row_ids.len() {
        bail!("positions and row ids must have the same length");
    }
    if count_unique_positions(positions) != unique_positions {
        bail!("unique position count does not match positions array");
    }
    if !positions.windows(2).all(|window| window[0] <= window[1]) {
        bail!("positions must be monotonically sorted");
    }
    Ok(())
}

fn count_unique_positions(positions: &[u32]) -> usize {
    positions
        .iter()
        .copied()
        .fold((0usize, None), |(count, previous), position| {
            let next_count = if previous == Some(position) {
                count
            } else {
                count + 1
            };
            (next_count, Some(position))
        })
        .0
}

#[cfg(target_endian = "big")]
fn write_u32<W: Write>(writer: &mut W, value: u32) -> Result<()> {
    writer.write_all(&value.to_le_bytes())?;
    Ok(())
}

fn write_u64<W: Write>(writer: &mut W, value: u64) -> Result<()> {
    writer.write_all(&value.to_le_bytes())?;
    Ok(())
}

fn write_u32_slice<W: Write>(writer: &mut W, values: &[u32]) -> Result<()> {
    #[cfg(target_endian = "little")]
    {
        writer.write_all(plain_slice_as_bytes(values))?;
        Ok(())
    }

    #[cfg(target_endian = "big")]
    {
        for &value in values {
            write_u32(writer, value)?;
        }
        Ok(())
    }
}

fn write_u64_slice<W: Write>(writer: &mut W, values: &[u64]) -> Result<()> {
    #[cfg(target_endian = "little")]
    {
        writer.write_all(plain_slice_as_bytes(values))?;
        Ok(())
    }

    #[cfg(target_endian = "big")]
    {
        for &value in values {
            write_u64(writer, value)?;
        }
        Ok(())
    }
}

fn read_u32_vec<R: Read>(reader: &mut R, len: usize) -> Result<Vec<u32>> {
    let values = read_plain_vec::<R, u32>(reader, len)?;
    #[cfg(target_endian = "little")]
    {
        Ok(values)
    }

    #[cfg(target_endian = "big")]
    {
        let mut values = values;
        for value in &mut values {
            *value = u32::from_le(*value);
        }
        Ok(values)
    }
}

fn read_u64_vec<R: Read>(reader: &mut R, len: usize) -> Result<Vec<u64>> {
    let values = read_plain_vec::<R, u64>(reader, len)?;
    #[cfg(target_endian = "little")]
    {
        Ok(values)
    }

    #[cfg(target_endian = "big")]
    {
        let mut values = values;
        for value in &mut values {
            *value = u64::from_le(*value);
        }
        Ok(values)
    }
}

fn read_plain_vec<R: Read, T: Copy>(reader: &mut R, len: usize) -> Result<Vec<T>> {
    let byte_len = len
        .checked_mul(size_of::<T>())
        .context("sidecar array byte length overflow")?;
    let mut values = Vec::<T>::with_capacity(len);
    // SAFETY: u32/u64 are plain-old-data values. We read exactly len * size_of::<T>() bytes into
    // spare capacity and only mark the vector initialized after read_exact succeeds.
    unsafe {
        reader.read_exact(slice::from_raw_parts_mut(
            values.spare_capacity_mut().as_mut_ptr().cast::<u8>(),
            byte_len,
        ))?;
        values.set_len(len);
    }
    Ok(values)
}

#[cfg(target_endian = "little")]
fn plain_slice_as_bytes<T: Copy>(values: &[T]) -> &[u8] {
    let byte_len = std::mem::size_of_val(values);
    // SAFETY: little-endian u32/u64 memory representation is the sidecar on-disk format.
    unsafe { slice::from_raw_parts(values.as_ptr().cast::<u8>(), byte_len) }
}

fn read_u64<R: Read>(reader: &mut R) -> Result<u64> {
    let mut bytes = [0_u8; 8];
    reader.read_exact(&mut bytes)?;
    Ok(u64::from_le_bytes(bytes))
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::io::{Cursor, Read};

    struct CountingRead {
        inner: Cursor<Vec<u8>>,
        read_calls: usize,
    }

    impl CountingRead {
        fn new(bytes: Vec<u8>) -> Self {
            Self {
                inner: Cursor::new(bytes),
                read_calls: 0,
            }
        }
    }

    impl Read for CountingRead {
        fn read(&mut self, buf: &mut [u8]) -> std::io::Result<usize> {
            self.read_calls += 1;
            self.inner.read(buf)
        }
    }

    #[test]
    fn sidecar_groups_duplicate_positions_and_resolves_sorted_queries() {
        let index = PositionRowIdIndex::from_pairs(vec![
            PositionRowIdPair {
                position: 20,
                row_id: 200,
            },
            PositionRowIdPair {
                position: 10,
                row_id: 100,
            },
            PositionRowIdPair {
                position: 10,
                row_id: 101,
            },
            PositionRowIdPair {
                position: 30,
                row_id: 300,
            },
        ]);

        assert_eq!(index.unique_positions(), 3);
        assert_eq!(index.row_ids_len(), 4);

        let resolved = index.resolve_sorted_positions(&[5, 10, 20, 25, 30]);
        assert_eq!(resolved.matched_positions, 3);
        assert_eq!(resolved.row_ids, vec![100, 101, 200, 300]);
    }

    #[test]
    fn cursor_resolver_carries_progress_across_sorted_chunks() {
        let index = PositionRowIdIndex::from_pairs(vec![
            PositionRowIdPair {
                position: 20,
                row_id: 200,
            },
            PositionRowIdPair {
                position: 10,
                row_id: 100,
            },
            PositionRowIdPair {
                position: 10,
                row_id: 101,
            },
            PositionRowIdPair {
                position: 30,
                row_id: 300,
            },
            PositionRowIdPair {
                position: 50,
                row_id: 500,
            },
        ]);

        let query = [5, 10, 20, 25, 30, 50];
        let expected = index.resolve_sorted_positions(&query);

        let mut cursor = 0usize;
        let mut matched_positions = 0usize;
        let mut row_ids = Vec::new();
        for chunk in query.chunks(2) {
            let resolved = index.resolve_sorted_positions_from_cursor(chunk, &mut cursor);
            matched_positions += resolved.matched_positions;
            row_ids.extend(resolved.row_ids);
        }

        assert_eq!(matched_positions, expected.matched_positions);
        assert_eq!(row_ids, expected.row_ids);
        assert_eq!(cursor, index.row_ids_len());
    }

    #[test]
    fn sidecar_round_trips_binary_format() {
        let index = PositionRowIdIndex::from_pairs(vec![
            PositionRowIdPair {
                position: 7,
                row_id: 1,
            },
            PositionRowIdPair {
                position: 8,
                row_id: 2,
            },
            PositionRowIdPair {
                position: 8,
                row_id: 3,
            },
        ]);

        let mut bytes = Vec::new();
        index.write_to(&mut bytes).unwrap();
        let loaded = PositionRowIdIndex::read_from(bytes.as_slice()).unwrap();

        assert_eq!(loaded.positions(), &[7, 8, 8]);
        assert_eq!(loaded.unique_positions(), 2);
        assert_eq!(loaded.resolve_sorted_positions(&[8]).row_ids, vec![2, 3]);
    }

    #[test]
    fn sidecar_reader_bulk_reads_primitive_arrays() {
        let index = PositionRowIdIndex::from_pairs(vec![
            PositionRowIdPair {
                position: 7,
                row_id: 1,
            },
            PositionRowIdPair {
                position: 8,
                row_id: 2,
            },
            PositionRowIdPair {
                position: 8,
                row_id: 3,
            },
            PositionRowIdPair {
                position: 9,
                row_id: 4,
            },
        ]);

        let mut bytes = Vec::new();
        index.write_to(&mut bytes).unwrap();
        let mut reader = CountingRead::new(bytes);
        let loaded = PositionRowIdIndex::read_from(&mut reader).unwrap();

        assert_eq!(loaded.positions(), &[7, 8, 8, 9]);
        assert!(
            reader.read_calls <= 5,
            "sidecar reader should bulk-read header and arrays, got {} read calls",
            reader.read_calls
        );
    }

    #[test]
    fn pair_set_mismatch_ignores_row_id_order_within_position() {
        let left = PositionRowIdIndex::from_pairs(vec![
            PositionRowIdPair {
                position: 10,
                row_id: 2,
            },
            PositionRowIdPair {
                position: 10,
                row_id: 1,
            },
            PositionRowIdPair {
                position: 20,
                row_id: 3,
            },
        ]);
        let right = PositionRowIdIndex::from_pairs(vec![
            PositionRowIdPair {
                position: 10,
                row_id: 1,
            },
            PositionRowIdPair {
                position: 10,
                row_id: 2,
            },
            PositionRowIdPair {
                position: 20,
                row_id: 3,
            },
        ]);

        assert_eq!(left.pair_set_mismatch(&right), None);
    }

    #[test]
    fn pair_set_mismatch_reports_first_different_position_group() {
        let left = PositionRowIdIndex::from_pairs(vec![
            PositionRowIdPair {
                position: 10,
                row_id: 1,
            },
            PositionRowIdPair {
                position: 20,
                row_id: 3,
            },
        ]);
        let right = PositionRowIdIndex::from_pairs(vec![
            PositionRowIdPair {
                position: 10,
                row_id: 1,
            },
            PositionRowIdPair {
                position: 20,
                row_id: 4,
            },
        ]);

        let mismatch = left
            .pair_set_mismatch(&right)
            .expect("expected pair mismatch");
        assert!(mismatch.contains("position 20"));
        assert!(mismatch.contains("left_row_ids=[3]"));
        assert!(mismatch.contains("right_row_ids=[4]"));
    }
}
