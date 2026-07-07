use std::collections::HashSet;

#[derive(Debug)]
pub struct MissWorklist {
    /// Bare chromosome names (without "chr" prefix), e.g. "1", "X", "MT".
    pub chroms: HashSet<String>,
}

impl MissWorklist {
    /// Build a single-chrom worklist.
    ///
    /// Used by the partitioned path where the per-chrom Parquet dataset IS the
    /// filter, so no variant scanning is needed. The `expanded_chroms()` set
    /// will cover both bare and "chr"-prefixed forms.
    pub fn for_chrom(chrom: &str) -> Self {
        let bare = chrom.strip_prefix("chr").unwrap_or(chrom).to_string();
        let mut chroms = HashSet::new();
        chroms.insert(bare);
        MissWorklist { chroms }
    }

    pub fn is_empty(&self) -> bool {
        self.chroms.is_empty()
    }

    pub fn expanded_chroms(&self) -> HashSet<String> {
        let mut expanded = HashSet::new();
        for c in &self.chroms {
            let escaped = c.replace('\'', "''");
            expanded.insert(escaped.clone());
            if let Some(bare) = escaped.strip_prefix("chr") {
                expanded.insert(bare.to_string());
            } else {
                expanded.insert(format!("chr{escaped}"));
            }
        }
        expanded
    }

    pub fn chrom_filter_clause(&self) -> String {
        if self.chroms.is_empty() {
            return String::new();
        }
        let literals: Vec<String> = self
            .expanded_chroms()
            .iter()
            .map(|c| format!("'{c}'"))
            .collect();
        format!(" WHERE chrom IN ({})", literals.join(", "))
    }

    /// SQL filter for the per-chrom Parquet scan. The partitioned path filters by
    /// chromosome only (the per-chrom dataset is itself the position filter).
    pub fn interval_filter_sql(&self) -> String {
        self.chrom_filter_clause()
    }
}
