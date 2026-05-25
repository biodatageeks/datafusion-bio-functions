#!/usr/bin/env bash
# Build the Tier 1 v115 fixture at vep-benchmark/data/port/_cache115/ by
# slicing chr21 (+ MT) from a pre-built whole-genome vepyr cache plus the
# upstream Ensembl v115 native cache.
#
# Why not vepyr.build_cache(chromosomes=["21"]): the chromosomes kwarg
# does not exist in the installed vepyr version (verified 2026-05-26).
# So we drive a whole-genome build once via
#   uv run python -c "import vepyr; vepyr.build_cache(release=115, cache_dir=...)"
# and this script slices chr21 + MT out for the Tier 1 fixture.
#
# Idempotent: skips if _cache115/parquet/ is already populated.
#
# Pre-reqs:
#   - $VEPYR_WHOLE_GENOME_PARQUET points at a whole-genome vepyr parquet
#     output directory containing per-entity subdirs and chr21.parquet.
#     Default: /Users/wojtek/Documents/vepyr/_cache_v115/parquet/115_GRCh38_vep
#   - $VEP_V115_NATIVE_CACHE points at the unpacked Ensembl v115 cache
#     (containing info.txt + per-chromosome subdirs). Default:
#     /Users/wojtek/Documents/vepyr/_cache_v115/homo_sapiens/115_GRCh38
#   - $VEPYR_REFERENCE_FASTA points at a whole-genome reference FASTA with
#     a corresponding .fai. Default:
#     /Users/wojtek/Documents/vepyr/_cache_v115/Homo_sapiens.GRCh38.dna.primary_assembly.fa
#   - samtools on PATH
#
# Outputs (all under _cache115/):
#   parquet/115_GRCh38_vep/{transcript,variation,regulatory,motif_feature}/{21,MT}.parquet
#   parquet/115_GRCh38_vep/{variation,translation_sift}.fjall/   (if present in whole-genome build)
#   native_cache/homo_sapiens/115_GRCh38/{21,MT,info.txt}
#   reference.fa[.fai]   (chr21 + MT slice)
set -euo pipefail

REPO_ROOT="$(cd "$(dirname "$0")/../.." && pwd)"
OUT="$REPO_ROOT/vep-benchmark/data/port/_cache115"

WHOLE_GENOME_PARQUET="${VEPYR_WHOLE_GENOME_PARQUET:-/Users/wojtek/Documents/vepyr/_cache_v115/parquet/115_GRCh38_vep}"
NATIVE_CACHE="${VEP_V115_NATIVE_CACHE:-/Users/wojtek/Documents/vepyr/_cache_v115/homo_sapiens/115_GRCh38}"
REFERENCE_FASTA="${VEPYR_REFERENCE_FASTA:-/Users/wojtek/Documents/vepyr/_cache_v115/Homo_sapiens.GRCh38.dna.primary_assembly.fa}"

if [ -d "$OUT/parquet/115_GRCh38_vep" ] && [ -d "$OUT/native_cache" ] && [ -f "$OUT/reference.fa" ]; then
  echo "Tier 1 v115 fixture already built at $OUT — skipping (delete to rebuild)"
  exit 0
fi

if ! command -v samtools >/dev/null 2>&1; then
  echo "ERROR: samtools not on PATH" >&2
  exit 1
fi
if [ ! -d "$WHOLE_GENOME_PARQUET" ]; then
  echo "ERROR: whole-genome vepyr parquet not found at $WHOLE_GENOME_PARQUET" >&2
  echo "  Build it first via: uv run python -c 'import vepyr; vepyr.build_cache(release=115, cache_dir=\"<root>\")'" >&2
  exit 1
fi
if [ ! -d "$NATIVE_CACHE" ]; then
  echo "ERROR: native v115 cache not found at $NATIVE_CACHE" >&2
  echo "  Either set VEP_V115_NATIVE_CACHE or rebuild the whole-genome cache to populate it." >&2
  exit 1
fi
if [ ! -f "$REFERENCE_FASTA" ]; then
  echo "ERROR: reference FASTA not found at $REFERENCE_FASTA" >&2
  exit 1
fi

CHROMS_PARQUET=(21 MT)
CHROMS_NATIVE=(21 MT)

# ── parquet slice ─────────────────────────────────────────────────────
echo "Slicing parquet chr21+MT from $WHOLE_GENOME_PARQUET"
mkdir -p "$OUT/parquet/115_GRCh38_vep"
for entity_dir in "$WHOLE_GENOME_PARQUET"/*/; do
  entity="$(basename "$entity_dir")"
  # Skip non-entity files like *.fjall (handled separately below)
  case "$entity" in
    *.fjall) continue ;;
  esac
  mkdir -p "$OUT/parquet/115_GRCh38_vep/$entity"
  for chrom in "${CHROMS_PARQUET[@]}"; do
    src="$entity_dir$chrom.parquet"
    if [ -f "$src" ]; then
      cp "$src" "$OUT/parquet/115_GRCh38_vep/$entity/$chrom.parquet"
      echo "  copied $entity/$chrom.parquet ($(du -h "$src" | cut -f1))"
    fi
  done
done

# Copy fjall stores in full (they're per-genome KV, can't be sliced trivially).
# Keep only if the directories exist; otherwise port_common runs parquet-only.
for fjall in "$WHOLE_GENOME_PARQUET"/*.fjall; do
  [ -d "$fjall" ] || continue
  name="$(basename "$fjall")"
  echo "  copying $name (full; KV stores are not per-chrom)"
  cp -r "$fjall" "$OUT/parquet/115_GRCh38_vep/$name"
done

# ── native cache slice (for the Perl-VEP oracle in Docker) ────────────
echo "Slicing native cache chr21+MT from $NATIVE_CACHE"
mkdir -p "$OUT/native_cache/homo_sapiens/115_GRCh38"
cp "$NATIVE_CACHE/info.txt" "$OUT/native_cache/homo_sapiens/115_GRCh38/"
for chrom in "${CHROMS_NATIVE[@]}"; do
  if [ -d "$NATIVE_CACHE/$chrom" ]; then
    cp -r "$NATIVE_CACHE/$chrom" "$OUT/native_cache/homo_sapiens/115_GRCh38/"
    echo "  copied native/$chrom/ ($(du -sh "$NATIVE_CACHE/$chrom" | cut -f1))"
  fi
done

# ── reference FASTA slice ─────────────────────────────────────────────
echo "Slicing reference FASTA (chr21 + MT)"
# samtools faidx -r writes one record per name to stdout, preserving the >chrom header.
samtools faidx "$REFERENCE_FASTA" 21 MT > "$OUT/reference.fa"
samtools faidx "$OUT/reference.fa"
echo "  reference.fa: $(du -h "$OUT/reference.fa" | cut -f1)"

echo "DONE → $OUT"
echo "Sizes:"
du -sh "$OUT"/{parquet,native_cache,reference.fa,reference.fa.fai} 2>/dev/null || true
