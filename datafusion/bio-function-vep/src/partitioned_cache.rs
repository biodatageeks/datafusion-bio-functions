//! Session table helpers shared by the partitioned cache backends.

use datafusion::common::Result;
use datafusion::prelude::SessionContext;

/// Deregister an ephemeral table from the session.
pub async fn deregister_table(session: &SessionContext, name: &str) -> Result<()> {
    // deregister_table returns Option<Arc<dyn TableProvider>>; ignore it.
    let _ = session.deregister_table(name)?;
    Ok(())
}
