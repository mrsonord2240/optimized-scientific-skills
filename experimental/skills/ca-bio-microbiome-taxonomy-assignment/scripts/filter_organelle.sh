#!/usr/bin/env bash
# Checked on QIIME2 2024.10+
# Remove Mitochondria and Chloroplast features after taxonomy assignment, before diversity or differential abundance.
# The exclude is a case-insensitive substring match on the lineage, so it works for SILVA, GTDB, Greengenes2 and RDP.
#
# Usage: filter_organelle.sh <table.qza> <taxonomy.qza> <rep-seqs.qza> <outdir>

set -euo pipefail

TABLE=${1:-}
TAX=${2:-}
SEQS=${3:-}
OUT=${4:-}
if [ -z "$TABLE" ] || [ -z "$TAX" ] || [ -z "$SEQS" ] || [ -z "$OUT" ]; then
    echo "Usage: filter_organelle.sh table.qza taxonomy.qza rep-seqs.qza outdir" >&2
    exit 1
fi
for f in "$TABLE" "$TAX" "$SEQS"; do
    [ -f "$f" ] || { echo "ERROR: not found: $f" >&2; exit 1; }
done
mkdir -p "$OUT"

qiime taxa filter-table --i-table "$TABLE" --i-taxonomy "$TAX" \
    --p-exclude mitochondria,chloroplast \
    --o-filtered-table "$OUT/table-no-organelle.qza"
qiime taxa filter-seqs --i-data "$SEQS" --i-taxonomy "$TAX" \
    --p-exclude mitochondria,chloroplast \
    --o-filtered-data "$OUT/rep-seqs-no-organelle.qza"

echo "Done: $OUT/table-no-organelle.qza $OUT/rep-seqs-no-organelle.qza (check per-sample read retention)"
