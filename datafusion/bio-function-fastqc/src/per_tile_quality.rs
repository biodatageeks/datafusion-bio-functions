use std::collections::BTreeMap;

use crate::{QcModule, TidyRow};

#[derive(Debug, Default, Clone)]
struct TileAcc {
    sum: Vec<f64>,
    count: Vec<u64>,
}

impl TileAcc {
    fn add(&mut self, pos: usize, phred: f64) {
        if self.sum.len() <= pos {
            self.sum.resize(pos + 1, 0.0);
            self.count.resize(pos + 1, 0);
        }
        self.sum[pos] += phred;
        self.count[pos] += 1;
    }
}

/// FastQC "Per tile sequence quality": per (tile, position), the mean quality
/// minus the *unweighted* mean of all tiles' means at that position (FastQC
/// divides the summed per-tile means by the total number of tiles).
#[derive(Debug, Default)]
pub struct PerTileQuality {
    tiles: BTreeMap<u32, TileAcc>,
    split_position: Option<usize>,
}

impl PerTileQuality {
    pub fn new() -> Self {
        Self::default()
    }

    /// Parse the tile from the read header, mirroring FastQC's `getTile`:
    /// split the id on ':'; >=7 fields -> index 4, else >=5 -> index 2, else
    /// no tile. The chosen index is cached for subsequent reads.
    fn tile_of(&mut self, name: &[u8]) -> Option<u32> {
        let id = std::str::from_utf8(name).ok()?;
        let fields: Vec<&str> = id.split(':').collect();
        let idx = match self.split_position {
            Some(i) => i,
            None => {
                let i = if fields.len() >= 7 {
                    4
                } else if fields.len() >= 5 {
                    2
                } else {
                    return None;
                };
                self.split_position = Some(i);
                i
            }
        };
        fields.get(idx)?.trim().parse::<u32>().ok()
    }
}

impl QcModule for PerTileQuality {
    fn as_any(&self) -> &dyn std::any::Any {
        self
    }

    fn name(&self) -> &'static str {
        "per_tile_quality"
    }

    fn update(&mut self, name: &[u8], _seq: &[u8], qual: &[u8]) {
        let Some(tile) = self.tile_of(name) else {
            return;
        };
        let acc = self.tiles.entry(tile).or_default();
        for (pos, &b) in qual.iter().enumerate() {
            acc.add(pos, (b as f64) - 33.0);
        }
    }

    fn merge(&mut self, other: &dyn QcModule) {
        let o = other
            .as_any()
            .downcast_ref::<PerTileQuality>()
            .expect("merge type mismatch");
        if self.split_position.is_none() {
            self.split_position = o.split_position;
        }
        for (tile, oacc) in &o.tiles {
            let acc = self.tiles.entry(*tile).or_default();
            let n = oacc.sum.len();
            if acc.sum.len() < n {
                acc.sum.resize(n, 0.0);
                acc.count.resize(n, 0);
            }
            for i in 0..n {
                acc.sum[i] += oacc.sum[i];
                acc.count[i] += oacc.count[i];
            }
        }
    }

    fn finalize(&self, out: &mut Vec<TidyRow>) {
        let m = "per_tile_quality";
        if self.tiles.is_empty() {
            out.push(TidyRow::status(m, "PASS"));
            return;
        }
        let max_len = self.tiles.values().map(|a| a.sum.len()).max().unwrap_or(0);
        let n_tiles = self.tiles.len() as f64;
        // Unweighted mean of per-tile means at each position. Fixtures use
        // uniform read length so every tile has data at every position.
        let mut group_mean = vec![0f64; max_len];
        for acc in self.tiles.values() {
            for (pos, gm) in group_mean.iter_mut().enumerate() {
                let tm = if pos < acc.count.len() && acc.count[pos] > 0 {
                    acc.sum[pos] / acc.count[pos] as f64
                } else {
                    0.0
                };
                *gm += tm;
            }
        }
        for gm in group_mean.iter_mut() {
            *gm /= n_tiles;
        }
        let mut max_abs_dev = 0f64;
        for (tile, acc) in &self.tiles {
            for (pos, (&s, &c)) in acc.sum.iter().zip(acc.count.iter()).enumerate() {
                if c == 0 {
                    continue;
                }
                let tm = s / c as f64;
                let dev = tm - group_mean[pos];
                max_abs_dev = max_abs_dev.max(dev.abs());
                out.push(TidyRow {
                    module: m,
                    label: Some(tile.to_string()),
                    position: Some((pos + 1) as i32),
                    metric: "mean".to_string(),
                    value: Some(dev),
                    value_str: None,
                });
            }
        }
        let status = if max_abs_dev > 10.0 {
            "FAIL"
        } else if max_abs_dev > 5.0 {
            "WARN"
        } else {
            "PASS"
        };
        out.push(TidyRow::status(m, status));
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    // Header format: split on ':' has >=7 fields -> tile is field index 4.
    fn hdr(tile: u32) -> Vec<u8> {
        format!("SRR.1 D0:7:HG:1:{tile}:1330:1935/1").into_bytes()
    }

    #[test]
    fn deviation_is_tile_mean_minus_unweighted_group_mean() {
        let mut m = PerTileQuality::new();
        // tile 1101: two reads, quality 'I'(40) at both positions.
        m.update(&hdr(1101), b"AC", b"II");
        m.update(&hdr(1101), b"AC", b"II");
        // tile 1102: one read, quality '5'(20) at both positions.
        m.update(&hdr(1102), b"AC", b"55");
        let mut rows = Vec::new();
        m.finalize(&mut rows);
        // group mean per position = (40 + 20)/2 = 30 (unweighted over 2 tiles).
        // tile 1101 deviation = 40 - 30 = +10; tile 1102 = 20 - 30 = -10.
        let dev = |tile: &str, pos: i32| {
            rows.iter()
                .find(|r| {
                    r.label.as_deref() == Some(tile)
                        && r.position == Some(pos)
                        && r.metric == "mean"
                })
                .and_then(|r| r.value)
                .unwrap()
        };
        assert!((dev("1101", 1) - 10.0).abs() < 1e-9);
        assert!((dev("1102", 1) + 10.0).abs() < 1e-9);
        // max abs deviation = 10 -> not > 10 -> WARN (>5), per limits tile warn=5.
        let status = rows
            .iter()
            .find(|r| r.metric == "status")
            .and_then(|r| r.value_str.clone())
            .unwrap();
        assert_eq!(status, "WARN");
    }
}
