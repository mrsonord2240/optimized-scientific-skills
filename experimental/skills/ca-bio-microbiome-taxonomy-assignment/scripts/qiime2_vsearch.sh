#!/usr/bin/env bash
# Checked on QIIME2 2024.10+ (q2-feature-classifier), vsearch 2.22+
# Alignment-consensus classification: no trained model, so no scikit-learn version pinning.
#
# Usage: qiime2_vsearch.sh <rep-seqs.qza> <ref-seqs.qza> <ref-tax.qza> <outdir>
#   For amplicon data pass a region-extracted reference (qiime feature-classifier extract-reads).
# Env:   MAXACCEPTS 10   PERC_IDENTITY 0.8   MIN_CONSENSUS 0.51   THREADS 8
# A rank is reported only if MIN_CONSENSUS of the accepted top hits agree on it.

set -euo pipefail

QUERY=${1:-}
REF_SEQS=${2:-}
REF_TAX=${3:-}
OUT=${4:-}
if [ -z "$QUERY" ] || [ -z "$REF_SEQS" ] || [ -z "$REF_TAX" ] || [ -z "$OUT" ]; then
    echo "Usage: qiime2_vsearch.sh rep-seqs.qza ref-seqs.qza ref-tax.qza outdir" >&2
    exit 1
fi
for f in "$QUERY" "$REF_SEQS" "$REF_TAX"; do
    [ -f "$f" ] || { echo "ERROR: not found: $f" >&2; exit 1; }
done
mkdir -p "$OUT"

qiime feature-classifier classify-consensus-vsearch \
    --i-query "$QUERY" \
    --i-reference-reads "$REF_SEQS" \
    --i-reference-taxonomy "$REF_TAX" \
    --p-maxaccepts "${MAXACCEPTS:-10}" --p-perc-identity "${PERC_IDENTITY:-0.8}" --p-min-consensus "${MIN_CONSENSUS:-0.51}" \
    --p-threads "${THREADS:-8}" \
    --o-classification "$OUT/taxonomy-vsearch.qza" --o-search-results "$OUT/vsearch-hits.qza"

echo "Done: $OUT/taxonomy-vsearch.qza"
