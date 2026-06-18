#!/usr/bin/env bash
# Parse a vepyr profiling stderr capture into a structured timing breakdown.
# Usage: breakdown.sh /tmp/prof.profile
set -euo pipefail
P="${1:?usage: breakdown.sh <profile-stderr-file>}"

echo "### PIPELINE (VEP_PROFILE) ###"
grep -a "pipeline_profile" "$P" | tail -1 | tr ' ' '\n' \
  | grep -aE "output_rows=|annotate=|engine=|lookup_wait=|hydrate=|context_load=|tx_window=|prepared_ctx=" || true

echo
echo "### LANCE LOOKUP (VEP_LANCE_PROFILE) ###"
grep -a "vep-lance-profile" "$P" | awk '{
  ev=$2; s=0;
  for(i=1;i<=NF;i++){ if($i ~ /^(seconds|open_s)=/){ split($i,a,"="); s=a[2] } }
  tot[ev]+=s; c[ev]++
} END { for(e in tot) printf "  %-22s n=%-5d total=%.3fs\n", e, c[e], tot[e] }' | sort -t= -k3 -rn
grep -a "variation_take" "$P" | grep -aoE "seconds=[0-9.]+" | sed 's/seconds=//' | sort -n \
  | awk '{a[NR]=$1; s+=$1} END{ if(NR) printf "  variation_take/batch: min=%.3f median=%.3f mean=%.3f max=%.3f (n=%d)\n", a[1], a[int(NR/2)+1], s/NR, a[NR], NR }'

echo
echo "### ENGINE — top level (VEP_ENGINE_PROFILE), additive ###"
grep -a "VEP_ENGINE_PROFILE" "$P" | grep -aoE "[a-z_0-9]+=[0-9]+\.[0-9]+s" \
  | awk -F= '{ v=$2; sub(/s$/,"",v); t[$1]+=v } END { for(k in t) printf "  %-26s %.3fs\n", k, t[k] }' | sort -k2 -rn

echo
echo "### TRANSCRIPT ENGINE — nested in evaluate_prepared (VEP_TX_ENGINE_PROFILE) ###"
grep -a "VEP_TX_ENGINE_PROFILE" "$P" | grep -aoE "[a-z_0-9]+=[0-9]+\.[0-9]+s" \
  | awk -F= '{ v=$2; sub(/s$/,"",v); t[$1]+=v } END { for(k in t) printf "  %-34s %.3fs\n", k, t[k] }' | sort -k2 -rn
