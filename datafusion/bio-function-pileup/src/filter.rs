/// Default SAM flag filter: unmapped(4) | secondary(256) | failed_qc(512) | duplicate(1024)
pub const DEFAULT_FILTER_FLAG: u32 = 1796;

/// Default minimum mapping quality.
pub const DEFAULT_MIN_MAPPING_QUALITY: u32 = 0;

/// Filter for alignment reads based on SAM flags and mapping quality.
#[derive(Debug, Clone)]
pub struct ReadFilter {
    /// Bitwise flag mask. Reads with any of these flags set are filtered out.
    pub filter_flag: u32,
    /// Minimum mapping quality. Reads below this threshold are filtered out.
    pub min_mapping_quality: u32,
}

impl Default for ReadFilter {
    fn default() -> Self {
        Self {
            filter_flag: DEFAULT_FILTER_FLAG,
            min_mapping_quality: DEFAULT_MIN_MAPPING_QUALITY,
        }
    }
}

impl ReadFilter {
    /// Check if a read passes the filter criteria.
    ///
    /// Returns `true` if:
    /// - None of the `filter_flag` bits are set in the read's flags
    /// - The read's mapping quality meets or exceeds `min_mapping_quality`
    pub fn passes(&self, flags: u32, mapping_quality: u32) -> bool {
        (flags & self.filter_flag) == 0 && mapping_quality >= self.min_mapping_quality
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_default_filter() {
        let filter = ReadFilter::default();
        assert_eq!(filter.filter_flag, 1796);
        assert_eq!(filter.min_mapping_quality, 0);
    }

    #[test]
    fn test_passing_read() {
        let filter = ReadFilter::default();
        assert!(filter.passes(0, 60));
    }

    #[test]
    fn test_unmapped_read() {
        let filter = ReadFilter::default();
        assert!(!filter.passes(4, 60));
    }

    #[test]
    fn test_secondary_read() {
        let filter = ReadFilter::default();
        assert!(!filter.passes(256, 60));
    }

    #[test]
    fn test_duplicate_read() {
        let filter = ReadFilter::default();
        assert!(!filter.passes(1024, 60));
    }

    #[test]
    fn test_failed_qc_read() {
        let filter = ReadFilter::default();
        assert!(!filter.passes(512, 60));
    }

    #[test]
    fn test_combined_flags() {
        let filter = ReadFilter::default();
        assert!(!filter.passes(4 | 1024, 60));
    }

    #[test]
    fn test_proper_pair_passes() {
        let filter = ReadFilter::default();
        assert!(filter.passes(2, 30));
    }

    #[test]
    fn test_mapq_threshold() {
        let filter = ReadFilter {
            filter_flag: DEFAULT_FILTER_FLAG,
            min_mapping_quality: 30,
        };
        assert!(filter.passes(0, 30));
        assert!(filter.passes(0, 60));
        assert!(!filter.passes(0, 29));
        assert!(!filter.passes(0, 0));
    }

    #[test]
    fn test_mapq_zero_always_passes_default() {
        let filter = ReadFilter::default();
        assert!(filter.passes(0, 0));
    }
}
