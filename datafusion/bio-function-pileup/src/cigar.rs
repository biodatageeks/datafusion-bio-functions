// Binary BAM CIGAR op codes (low 4 bits of packed u32).
const BAM_CIGAR_M: u32 = 0;
#[cfg(test)]
const BAM_CIGAR_I: u32 = 1;
const BAM_CIGAR_D: u32 = 2;
const BAM_CIGAR_N: u32 = 3;
#[cfg(test)]
const BAM_CIGAR_S: u32 = 4;
#[cfg(test)]
const BAM_CIGAR_H: u32 = 5;
// BAM_CIGAR_P (6, padding) is part of the BAM spec but not used in coverage computation.
const BAM_CIGAR_EQ: u32 = 7;
const BAM_CIGAR_X: u32 = 8;

/// Represents a CIGAR operation kind.
#[derive(Debug, Clone, Copy, PartialEq)]
pub enum CigarOpKind {
    /// M, X, = — consumes reference and produces coverage events
    AlignmentMatch,
    /// D, N — consumes reference but no coverage events
    ReferenceSkip,
    /// I, S, H, P — does not consume reference
    NoRefConsumption,
}

/// A single CIGAR operation with its kind and length.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct CigarOp {
    pub kind: CigarOpKind,
    pub len: u32,
}

/// Represents a coverage event: a position delta pair.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct CoverageEvent {
    pub position: u32,
    pub delta: i32,
}

/// Parse a CIGAR string into a list of operations.
pub fn parse_cigar(cigar: &str) -> Vec<CigarOp> {
    let mut ops = Vec::new();
    let mut num_start = 0;

    for (i, c) in cigar.char_indices() {
        if c.is_ascii_digit() {
            continue;
        }
        let len: u32 = cigar[num_start..i].parse().unwrap_or(0);
        let kind = match c {
            'M' | 'X' | '=' => CigarOpKind::AlignmentMatch,
            'D' | 'N' => CigarOpKind::ReferenceSkip,
            'I' | 'S' | 'H' | 'P' => CigarOpKind::NoRefConsumption,
            _ => {
                num_start = i + c.len_utf8();
                continue;
            }
        };
        ops.push(CigarOp { kind, len });
        num_start = i + c.len_utf8();
    }
    ops
}

/// Apply CIGAR operations directly to a dense depth array.
///
/// For each alignment match (M/X/=), increments depth[start] by +1 and depth[end] by -1.
/// This is the mosdepth-style approach: O(1) writes with no per-read allocation.
/// Bounds checks prevent writes past the array end.
///
/// Parses the CIGAR string inline to avoid allocating an intermediate `Vec<CigarOp>`.
pub fn apply_cigar_to_depth(start: u32, cigar: &str, depth: &mut [i32]) {
    let bytes = cigar.as_bytes();
    let depth_len = depth.len();
    let mut ref_pos = start as usize;
    let mut num_start = 0;

    for i in 0..bytes.len() {
        let b = bytes[i];
        if b.is_ascii_digit() {
            continue;
        }
        let len: usize = cigar[num_start..i].parse().unwrap_or(0);
        match b {
            b'M' | b'X' | b'=' => {
                let end = ref_pos + len;
                if ref_pos < depth_len {
                    depth[ref_pos] += 1;
                }
                if end < depth_len {
                    depth[end] -= 1;
                }
                ref_pos = end;
            }
            b'D' | b'N' => {
                ref_pos += len;
            }
            _ => {} // I, S, H, P — no ref consumption
        }
        num_start = i + 1;
    }
}

/// Push coverage events directly to a `Vec<(u32, i32)>`, avoiding per-read allocation.
///
/// Same logic as `generate_events` but appends to a caller-owned Vec.
/// This is the preferred hot path for the sparse accumulator.
///
/// Parses the CIGAR string inline to avoid allocating an intermediate `Vec<CigarOp>`.
#[inline]
pub fn apply_cigar_to_event_list(start: u32, cigar: &str, out: &mut Vec<(u32, i32)>) {
    let bytes = cigar.as_bytes();
    let mut ref_pos = start;
    let mut num_start = 0;

    for i in 0..bytes.len() {
        let b = bytes[i];
        if b.is_ascii_digit() {
            continue;
        }
        let len: u32 = cigar[num_start..i].parse().unwrap_or(0);
        match b {
            b'M' | b'X' | b'=' => {
                out.push((ref_pos, 1));
                ref_pos += len;
                out.push((ref_pos, -1));
            }
            b'D' | b'N' => {
                ref_pos += len;
            }
            _ => {} // I, S, H, P — no ref consumption
        }
        num_start = i + 1;
    }
}

/// Apply binary CIGAR operations directly to a dense depth array.
///
/// Same logic as `apply_cigar_to_depth` but reads from packed binary CIGAR bytes
/// (4 bytes per op, little-endian u32: op_len = packed >> 4, op_code = packed & 0xF).
pub fn apply_binary_cigar_to_depth(start: u32, cigar_bytes: &[u8], depth: &mut [i32]) {
    let depth_len = depth.len();
    let mut ref_pos = start as usize;

    let mut offset = 0;
    while offset + 4 <= cigar_bytes.len() {
        let packed = u32::from_le_bytes([
            cigar_bytes[offset],
            cigar_bytes[offset + 1],
            cigar_bytes[offset + 2],
            cigar_bytes[offset + 3],
        ]);
        let op_len = (packed >> 4) as usize;
        let op_code = packed & 0xF;
        offset += 4;

        match op_code {
            BAM_CIGAR_M | BAM_CIGAR_EQ | BAM_CIGAR_X => {
                let end = ref_pos + op_len;
                if ref_pos < depth_len {
                    depth[ref_pos] += 1;
                }
                if end < depth_len {
                    depth[end] -= 1;
                }
                ref_pos = end;
            }
            BAM_CIGAR_D | BAM_CIGAR_N => {
                ref_pos += op_len;
            }
            _ => {} // I, S, H, P — no ref consumption
        }
    }
}

/// Push coverage events from binary CIGAR to a `Vec<(u32, i32)>`.
///
/// Same logic as `apply_cigar_to_event_list` but reads from packed binary CIGAR bytes.
#[inline]
pub fn apply_binary_cigar_to_event_list(start: u32, cigar_bytes: &[u8], out: &mut Vec<(u32, i32)>) {
    let mut ref_pos = start;

    let mut offset = 0;
    while offset + 4 <= cigar_bytes.len() {
        let packed = u32::from_le_bytes([
            cigar_bytes[offset],
            cigar_bytes[offset + 1],
            cigar_bytes[offset + 2],
            cigar_bytes[offset + 3],
        ]);
        let op_len = packed >> 4;
        let op_code = packed & 0xF;
        offset += 4;

        match op_code {
            BAM_CIGAR_M | BAM_CIGAR_EQ | BAM_CIGAR_X => {
                out.push((ref_pos, 1));
                ref_pos += op_len;
                out.push((ref_pos, -1));
            }
            BAM_CIGAR_D | BAM_CIGAR_N => {
                ref_pos += op_len;
            }
            _ => {} // I, S, H, P — no ref consumption
        }
    }
}

/// Generate coverage events from a CIGAR string starting at a reference position.
///
/// For each alignment match (M/X/=), emits +1 at the start and -1 at the end.
/// Deletions and reference skips (D/N) advance the position but don't generate events.
/// Insertions, soft clips, hard clips, and padding (I/S/H/P) are ignored.
pub fn generate_events(start: u32, cigar: &str) -> Vec<CoverageEvent> {
    let ops = parse_cigar(cigar);
    let mut events = Vec::new();
    let mut ref_pos = start;

    for op in &ops {
        match op.kind {
            CigarOpKind::AlignmentMatch => {
                events.push(CoverageEvent {
                    position: ref_pos,
                    delta: 1,
                });
                ref_pos += op.len;
                events.push(CoverageEvent {
                    position: ref_pos,
                    delta: -1,
                });
            }
            CigarOpKind::ReferenceSkip => {
                ref_pos += op.len;
            }
            CigarOpKind::NoRefConsumption => {
                // Does not advance reference position
            }
        }
    }
    events
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_simple_match() {
        let events = generate_events(100, "10M");
        assert_eq!(
            events,
            vec![
                CoverageEvent {
                    position: 100,
                    delta: 1
                },
                CoverageEvent {
                    position: 110,
                    delta: -1
                },
            ]
        );
    }

    #[test]
    fn test_complex_cigar() {
        // 5M2I3M = 5 match, 2 insert, 3 match
        let events = generate_events(0, "5M2I3M");
        assert_eq!(
            events,
            vec![
                CoverageEvent {
                    position: 0,
                    delta: 1
                },
                CoverageEvent {
                    position: 5,
                    delta: -1
                },
                CoverageEvent {
                    position: 5,
                    delta: 1
                },
                CoverageEvent {
                    position: 8,
                    delta: -1
                },
            ]
        );
    }

    #[test]
    fn test_soft_clips() {
        // 3S10M2S
        let events = generate_events(50, "3S10M2S");
        assert_eq!(
            events,
            vec![
                CoverageEvent {
                    position: 50,
                    delta: 1
                },
                CoverageEvent {
                    position: 60,
                    delta: -1
                },
            ]
        );
    }

    #[test]
    fn test_deletion() {
        // 5M3D5M
        let events = generate_events(0, "5M3D5M");
        assert_eq!(
            events,
            vec![
                CoverageEvent {
                    position: 0,
                    delta: 1
                },
                CoverageEvent {
                    position: 5,
                    delta: -1
                },
                // deletion advances ref_pos by 3 (5 -> 8)
                CoverageEvent {
                    position: 8,
                    delta: 1
                },
                CoverageEvent {
                    position: 13,
                    delta: -1
                },
            ]
        );
    }

    #[test]
    fn test_insertion() {
        // 5M3I5M
        let events = generate_events(10, "5M3I5M");
        assert_eq!(
            events,
            vec![
                CoverageEvent {
                    position: 10,
                    delta: 1
                },
                CoverageEvent {
                    position: 15,
                    delta: -1
                },
                // insertion does NOT advance ref_pos
                CoverageEvent {
                    position: 15,
                    delta: 1
                },
                CoverageEvent {
                    position: 20,
                    delta: -1
                },
            ]
        );
    }

    #[test]
    fn test_reference_skip() {
        // 5M100N5M (spliced alignment)
        let events = generate_events(0, "5M100N5M");
        assert_eq!(
            events,
            vec![
                CoverageEvent {
                    position: 0,
                    delta: 1
                },
                CoverageEvent {
                    position: 5,
                    delta: -1
                },
                CoverageEvent {
                    position: 105,
                    delta: 1
                },
                CoverageEvent {
                    position: 110,
                    delta: -1
                },
            ]
        );
    }

    #[test]
    fn test_eq_and_x_operations() {
        // = and X are alignment match variants
        let events = generate_events(0, "3=2X5M");
        assert_eq!(
            events,
            vec![
                CoverageEvent {
                    position: 0,
                    delta: 1
                },
                CoverageEvent {
                    position: 3,
                    delta: -1
                },
                CoverageEvent {
                    position: 3,
                    delta: 1
                },
                CoverageEvent {
                    position: 5,
                    delta: -1
                },
                CoverageEvent {
                    position: 5,
                    delta: 1
                },
                CoverageEvent {
                    position: 10,
                    delta: -1
                },
            ]
        );
    }

    #[test]
    fn test_empty_cigar() {
        let events = generate_events(0, "");
        assert!(events.is_empty());
    }

    #[test]
    fn test_hard_clips() {
        // 5H10M5H
        let events = generate_events(0, "5H10M5H");
        assert_eq!(
            events,
            vec![
                CoverageEvent {
                    position: 0,
                    delta: 1
                },
                CoverageEvent {
                    position: 10,
                    delta: -1
                },
            ]
        );
    }

    // --- Tests for apply_cigar_to_depth ---

    #[test]
    fn test_dense_simple_match() {
        let mut depth = vec![0i32; 20];
        apply_cigar_to_depth(5, "10M", &mut depth);
        assert_eq!(depth[5], 1);
        assert_eq!(depth[15], -1);
        // All other positions should be 0
        assert_eq!(depth[0], 0);
        assert_eq!(depth[10], 0);
    }

    #[test]
    fn test_dense_two_overlapping_reads() {
        let mut depth = vec![0i32; 20];
        apply_cigar_to_depth(0, "10M", &mut depth);
        apply_cigar_to_depth(5, "10M", &mut depth);
        assert_eq!(depth[0], 1);
        assert_eq!(depth[5], 1); // +1 from second read
        assert_eq!(depth[10], -1); // -1 from first read
        assert_eq!(depth[15], -1); // -1 from second read
    }

    #[test]
    fn test_dense_insertion_no_ref_advance() {
        let mut depth = vec![0i32; 20];
        apply_cigar_to_depth(0, "5M3I5M", &mut depth);
        // 5M: +1@0, -1@5; 3I: no ref; 5M: +1@5, -1@10
        assert_eq!(depth[0], 1);
        assert_eq!(depth[5], 0); // -1 + 1 cancel out
        assert_eq!(depth[10], -1);
    }

    #[test]
    fn test_dense_deletion() {
        let mut depth = vec![0i32; 20];
        apply_cigar_to_depth(0, "5M3D5M", &mut depth);
        // 5M: +1@0, -1@5; 3D: skip 3; 5M: +1@8, -1@13
        assert_eq!(depth[0], 1);
        assert_eq!(depth[5], -1);
        assert_eq!(depth[8], 1);
        assert_eq!(depth[13], -1);
    }

    #[test]
    fn test_dense_bounds_check() {
        // Array is smaller than the alignment end position
        let mut depth = vec![0i32; 8];
        apply_cigar_to_depth(5, "10M", &mut depth);
        // +1 at pos 5 (within bounds), -1 at pos 15 (out of bounds, skipped)
        assert_eq!(depth[5], 1);
        // No panic, no out-of-bounds write
    }

    #[test]
    fn test_dense_reference_skip() {
        let mut depth = vec![0i32; 120];
        apply_cigar_to_depth(0, "5M100N5M", &mut depth);
        assert_eq!(depth[0], 1);
        assert_eq!(depth[5], -1);
        assert_eq!(depth[105], 1);
        assert_eq!(depth[110], -1);
    }

    #[test]
    fn test_dense_empty_cigar() {
        let mut depth = vec![0i32; 10];
        apply_cigar_to_depth(0, "", &mut depth);
        assert!(depth.iter().all(|&v| v == 0));
    }

    // --- Tests for binary CIGAR ---

    /// Encode a single CIGAR op as 4 LE bytes: packed = (len << 4) | code.
    fn encode_test_op(len: u32, code: u32) -> [u8; 4] {
        ((len << 4) | code).to_le_bytes()
    }

    #[test]
    fn test_binary_simple_match_events() {
        // 10M
        let bytes = encode_test_op(10, BAM_CIGAR_M);
        let mut out = Vec::new();
        apply_binary_cigar_to_event_list(100, &bytes, &mut out);
        assert_eq!(out, vec![(100, 1), (110, -1)]);
    }

    #[test]
    fn test_binary_complex_multi_op_events() {
        // 5M2I3M
        let mut bytes = Vec::new();
        bytes.extend_from_slice(&encode_test_op(5, BAM_CIGAR_M));
        bytes.extend_from_slice(&encode_test_op(2, BAM_CIGAR_I));
        bytes.extend_from_slice(&encode_test_op(3, BAM_CIGAR_M));
        let mut out = Vec::new();
        apply_binary_cigar_to_event_list(0, &bytes, &mut out);
        assert_eq!(out, vec![(0, 1), (5, -1), (5, 1), (8, -1)]);
    }

    #[test]
    fn test_binary_deletion_events() {
        // 5M3D5M
        let mut bytes = Vec::new();
        bytes.extend_from_slice(&encode_test_op(5, BAM_CIGAR_M));
        bytes.extend_from_slice(&encode_test_op(3, BAM_CIGAR_D));
        bytes.extend_from_slice(&encode_test_op(5, BAM_CIGAR_M));
        let mut out = Vec::new();
        apply_binary_cigar_to_event_list(0, &bytes, &mut out);
        assert_eq!(out, vec![(0, 1), (5, -1), (8, 1), (13, -1)]);
    }

    #[test]
    fn test_binary_reference_skip_events() {
        // 5M100N5M
        let mut bytes = Vec::new();
        bytes.extend_from_slice(&encode_test_op(5, BAM_CIGAR_M));
        bytes.extend_from_slice(&encode_test_op(100, BAM_CIGAR_N));
        bytes.extend_from_slice(&encode_test_op(5, BAM_CIGAR_M));
        let mut out = Vec::new();
        apply_binary_cigar_to_event_list(0, &bytes, &mut out);
        assert_eq!(out, vec![(0, 1), (5, -1), (105, 1), (110, -1)]);
    }

    #[test]
    fn test_binary_soft_hard_clips_events() {
        // 3S10M2S5H
        let mut bytes = Vec::new();
        bytes.extend_from_slice(&encode_test_op(3, BAM_CIGAR_S));
        bytes.extend_from_slice(&encode_test_op(10, BAM_CIGAR_M));
        bytes.extend_from_slice(&encode_test_op(2, BAM_CIGAR_S));
        bytes.extend_from_slice(&encode_test_op(5, BAM_CIGAR_H));
        let mut out = Vec::new();
        apply_binary_cigar_to_event_list(50, &bytes, &mut out);
        assert_eq!(out, vec![(50, 1), (60, -1)]);
    }

    #[test]
    fn test_binary_eq_x_codes_events() {
        // 3=2X5M
        let mut bytes = Vec::new();
        bytes.extend_from_slice(&encode_test_op(3, BAM_CIGAR_EQ));
        bytes.extend_from_slice(&encode_test_op(2, BAM_CIGAR_X));
        bytes.extend_from_slice(&encode_test_op(5, BAM_CIGAR_M));
        let mut out = Vec::new();
        apply_binary_cigar_to_event_list(0, &bytes, &mut out);
        assert_eq!(
            out,
            vec![(0, 1), (3, -1), (3, 1), (5, -1), (5, 1), (10, -1)]
        );
    }

    #[test]
    fn test_binary_empty_bytes_events() {
        let mut out = Vec::new();
        apply_binary_cigar_to_event_list(0, &[], &mut out);
        assert!(out.is_empty());
    }

    #[test]
    fn test_binary_simple_match_depth() {
        let bytes = encode_test_op(10, BAM_CIGAR_M);
        let mut depth = vec![0i32; 20];
        apply_binary_cigar_to_depth(5, &bytes, &mut depth);
        assert_eq!(depth[5], 1);
        assert_eq!(depth[15], -1);
        assert_eq!(depth[0], 0);
        assert_eq!(depth[10], 0);
    }

    #[test]
    fn test_binary_deletion_depth() {
        // 5M3D5M
        let mut bytes = Vec::new();
        bytes.extend_from_slice(&encode_test_op(5, BAM_CIGAR_M));
        bytes.extend_from_slice(&encode_test_op(3, BAM_CIGAR_D));
        bytes.extend_from_slice(&encode_test_op(5, BAM_CIGAR_M));
        let mut depth = vec![0i32; 20];
        apply_binary_cigar_to_depth(0, &bytes, &mut depth);
        assert_eq!(depth[0], 1);
        assert_eq!(depth[5], -1);
        assert_eq!(depth[8], 1);
        assert_eq!(depth[13], -1);
    }

    #[test]
    fn test_binary_bounds_check_depth() {
        let bytes = encode_test_op(10, BAM_CIGAR_M);
        let mut depth = vec![0i32; 8];
        apply_binary_cigar_to_depth(5, &bytes, &mut depth);
        assert_eq!(depth[5], 1);
        // No panic, -1 at pos 15 is out of bounds and skipped
    }

    #[test]
    fn test_binary_empty_bytes_depth() {
        let mut depth = vec![0i32; 10];
        apply_binary_cigar_to_depth(0, &[], &mut depth);
        assert!(depth.iter().all(|&v| v == 0));
    }
}
