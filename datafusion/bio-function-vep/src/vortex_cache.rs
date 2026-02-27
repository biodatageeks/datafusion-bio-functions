//! Helpers for registering Vortex-backed variation cache tables.

use std::sync::Arc;

use datafusion::common::{DataFusionError, Result};
use datafusion::datasource::file_format::FileFormatFactory;
use datafusion::datasource::listing::{
    ListingOptions, ListingTable, ListingTableConfig, ListingTableUrl,
};
use datafusion::prelude::SessionContext;
use vortex_datafusion::VortexFormatFactory;

/// Register a Vortex cache directory as a table for `lookup_variants()`.
///
/// The `cache_path` should point to a directory (or object-store prefix)
/// containing `.vortex` files with the standard variation cache schema.
pub async fn register_vortex_cache(
    ctx: &SessionContext,
    table_name: &str,
    cache_path: &str,
) -> Result<()> {
    if table_name.is_empty() {
        return Err(DataFusionError::Plan(
            "register_vortex_cache(): table_name must not be empty".to_string(),
        ));
    }
    if cache_path.is_empty() {
        return Err(DataFusionError::Plan(
            "register_vortex_cache(): cache_path must not be empty".to_string(),
        ));
    }

    let table_url = ListingTableUrl::parse(cache_path)?;
    let format = VortexFormatFactory::new().default();
    let listing_options =
        ListingOptions::new(format).with_session_config_options(ctx.state().config());

    let config = ListingTableConfig::new(table_url)
        .with_listing_options(listing_options)
        .infer_schema(&ctx.state())
        .await?;

    let table = Arc::new(ListingTable::try_new(config)?);
    ctx.register_table(table_name, table)?;
    Ok(())
}
