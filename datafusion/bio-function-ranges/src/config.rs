use datafusion::common::extensions_options;
use datafusion::config::ConfigExtension;

use crate::algorithms::Algorithm;

extensions_options! {
    pub struct BioConfig {
        pub prefer_interval_join: bool, default = true
        pub interval_join_algorithm: Algorithm, default = Algorithm::default()
        pub interval_join_low_memory: bool, default = false
    }
}

impl ConfigExtension for BioConfig {
    const PREFIX: &'static str = "bio";
}
