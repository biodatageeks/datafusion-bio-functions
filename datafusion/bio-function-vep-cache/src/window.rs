//! Window grouping and boundary handling for genomic variant windows.
//!
//! Groups sorted variants into ~1Mb windows by their start position.
//! Cross-window indels are stored in all windows they overlap.

use crate::key_encoding::window_id_for_position;

/// A genomic window identified by chromosome and window ID.
#[derive(Debug, Clone, PartialEq, Eq, Hash)]
pub struct GenomicWindow {
    pub chrom: String,
    pub window_id: u64,
}

/// Group sorted variants into windows based on their start position.
///
/// Input arrays must be sorted by (chrom, start). Returns a list of
/// (window, row_range) pairs where row_range indexes into the input arrays.
pub fn group_into_windows(
    chroms: &[&str],
    starts: &[i64],
    window_size: u64,
) -> Vec<(GenomicWindow, std::ops::Range<usize>)> {
    if chroms.is_empty() {
        return Vec::new();
    }

    let mut result = Vec::new();
    let mut current_chrom = chroms[0];
    let mut current_wid = window_id_for_position(starts[0], window_size);
    let mut range_start = 0;

    for i in 1..chroms.len() {
        let wid = window_id_for_position(starts[i], window_size);
        if chroms[i] != current_chrom || wid != current_wid {
            result.push((
                GenomicWindow {
                    chrom: current_chrom.to_string(),
                    window_id: current_wid,
                },
                range_start..i,
            ));
            current_chrom = chroms[i];
            current_wid = wid;
            range_start = i;
        }
    }

    // Push the final group.
    result.push((
        GenomicWindow {
            chrom: current_chrom.to_string(),
            window_id: current_wid,
        },
        range_start..chroms.len(),
    ));

    result
}

/// Return all windows that a variant overlaps (for cross-window indels).
///
/// A variant spanning [start, end] may cross window boundaries.
/// This returns the window IDs for every window the variant touches.
pub fn windows_for_variant(
    chrom: &str,
    start: i64,
    end: i64,
    window_size: u64,
) -> Vec<GenomicWindow> {
    let start_wid = window_id_for_position(start, window_size);
    let end_wid = window_id_for_position(end, window_size);

    (start_wid..=end_wid)
        .map(|wid| GenomicWindow {
            chrom: chrom.to_string(),
            window_id: wid,
        })
        .collect()
}

/// Compute the start position of a window.
pub fn window_start(window_id: u64, window_size: u64) -> i64 {
    (window_id * window_size) as i64
}

/// Compute the end position (exclusive) of a window.
pub fn window_end(window_id: u64, window_size: u64) -> i64 {
    ((window_id + 1) * window_size) as i64
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::key_encoding::DEFAULT_WINDOW_SIZE;

    #[test]
    fn test_group_single_window() {
        let chroms = vec!["1", "1", "1"];
        let starts = vec![100, 200, 300];
        let groups = group_into_windows(&chroms, &starts, DEFAULT_WINDOW_SIZE);
        assert_eq!(groups.len(), 1);
        assert_eq!(groups[0].0.chrom, "1");
        assert_eq!(groups[0].0.window_id, 0);
        assert_eq!(groups[0].1, 0..3);
    }

    #[test]
    fn test_group_multiple_windows() {
        let chroms = vec!["1", "1", "1"];
        let starts = vec![500_000, 1_500_000, 2_500_000];
        let groups = group_into_windows(&chroms, &starts, DEFAULT_WINDOW_SIZE);
        assert_eq!(groups.len(), 3);
        assert_eq!(groups[0].0.window_id, 0);
        assert_eq!(groups[0].1, 0..1);
        assert_eq!(groups[1].0.window_id, 1);
        assert_eq!(groups[1].1, 1..2);
        assert_eq!(groups[2].0.window_id, 2);
        assert_eq!(groups[2].1, 2..3);
    }

    #[test]
    fn test_group_chromosome_change() {
        let chroms = vec!["1", "1", "2"];
        let starts = vec![100, 200, 100];
        let groups = group_into_windows(&chroms, &starts, DEFAULT_WINDOW_SIZE);
        assert_eq!(groups.len(), 2);
        assert_eq!(groups[0].0.chrom, "1");
        assert_eq!(groups[1].0.chrom, "2");
    }

    #[test]
    fn test_group_empty() {
        let chroms: Vec<&str> = Vec::new();
        let starts: Vec<i64> = Vec::new();
        let groups = group_into_windows(&chroms, &starts, DEFAULT_WINDOW_SIZE);
        assert!(groups.is_empty());
    }

    #[test]
    fn test_windows_for_variant_single_window() {
        let windows = windows_for_variant("1", 100, 200, DEFAULT_WINDOW_SIZE);
        assert_eq!(windows.len(), 1);
        assert_eq!(windows[0].window_id, 0);
    }

    #[test]
    fn test_windows_for_variant_cross_boundary() {
        // Variant spanning the 1Mb boundary
        let windows = windows_for_variant("1", 999_990, 1_000_010, DEFAULT_WINDOW_SIZE);
        assert_eq!(windows.len(), 2);
        assert_eq!(windows[0].window_id, 0);
        assert_eq!(windows[1].window_id, 1);
    }

    #[test]
    fn test_window_start_end() {
        assert_eq!(window_start(0, DEFAULT_WINDOW_SIZE), 0);
        assert_eq!(window_end(0, DEFAULT_WINDOW_SIZE), 1_000_000);
        assert_eq!(window_start(1, DEFAULT_WINDOW_SIZE), 1_000_000);
        assert_eq!(window_end(1, DEFAULT_WINDOW_SIZE), 2_000_000);
    }
}
