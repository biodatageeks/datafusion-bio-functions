#!/usr/bin/env bash
# Parse a vepyr profiling stderr capture into a structured timing breakdown.
# Usage: breakdown.sh /tmp/prof.profile
set -euo pipefail
P="${1:?usage: breakdown.sh <profile-stderr-file>}"

echo "### PIPELINE (VEP_PROFILE) ###"
grep -a "pipeline_profile" "$P" | tail -1 | tr ' ' '\n' \
  | grep -aE "output_rows=|annotate=|engine=|lookup_wait=|hydrate=|context_load=|tx_window=|prepared_ctx=" || true

echo
echo "### VARIATION LOOKUP (VEP_LOOKUP_PROFILE) ###"
grep -a "vep-lookup-profile" "$P" || true

echo
echo "### ENGINE — top level (VEP_ENGINE_PROFILE), additive ###"
grep -a "VEP_ENGINE_PROFILE" "$P" | grep -aoE "[a-z_0-9]+=[0-9]+\.[0-9]+s" \
  | awk -F= '{ v=$2; sub(/s$/,"",v); t[$1]+=v } END { for(k in t) printf "  %-26s %.3fs\n", k, t[k] }' | sort -k2 -rn

echo
echo "### TRANSCRIPT ENGINE — nested in evaluate_prepared (VEP_TX_ENGINE_PROFILE) ###"
grep -a "VEP_TX_ENGINE_PROFILE" "$P" | grep -aoE "[a-z_0-9]+=[0-9]+\.[0-9]+s" \
  | awk -F= '{ v=$2; sub(/s$/,"",v); t[$1]+=v } END { for(k in t) printf "  %-34s %.3fs\n", k, t[k] }' | sort -k2 -rn
