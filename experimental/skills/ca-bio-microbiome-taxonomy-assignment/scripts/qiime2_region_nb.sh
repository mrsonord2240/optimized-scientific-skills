#!/usr/bin/env bash
# Checked on QIIME2 2024.10+ (q2-feature-classifier)
# Region-matched naive-Bayes classification: extract the primer region from a full-length reference,
# train on it, classify. A full-length classifier on V4 reads fabricates and erases calls.
#
# Usage: F_PRIMER=<fwd> R_PRIMER=<rev> qiime2_region_nb.sh <ref-seqs.qza> <ref-tax.qza> <rep-seqs.qza> <outdir>
# Env:   F_PRIMER, R_PRIMER  required; R_PRIMER given 5'->3', NOT reverse-complemented
#                            (515F/806R V4: GTGYCAGCMGCCGCGGTAA / GGACTACNVGGGTWTCTAAT)
#        CONFIDENCE  default 0.7 (below it the lineage is truncated to a shallower rank)
#        N_JOBS      default 1 (each job holds a copy of the classifier)

set -euo pipefail

REF_SEQS=${1:-}
REF_TAX=${2:-}
REP_SEQS=${3:-}
OUT=${4:-}
if [ -z "$REF_SEQS" ] || [ -z "$REF_TAX" ] || [ -z "$REP_SEQS" ] || [ -z "$OUT" ] || [ -z "${F_PRIMER:-}" ] || [ -z "${R_PRIMER:-}" ]; then
    echo "Usage: F_PRIMER=<fwd> R_PRIMER=<rev> qiime2_region_nb.sh ref-seqs.qza ref-tax.qza rep-seqs.qza outdir" >&2
    exit 1
fi
for f in "$REF_SEQS" "$REF_TAX" "$REP_SEQS"; do
    [ -f "$f" ] || { echo "ERROR: not found: $f" >&2; exit 1; }
done
mkdir -p "$OUT"

# 1. In-silico PCR so the reference k-mer composition matches the reads. --p-min-length 50 drops too-short amplicons.
qiime feature-classifier extract-reads \
    --i-sequences "$REF_SEQS" \
    --p-f-primer "$F_PRIMER" --p-r-primer "$R_PRIMER" \
    --p-min-length 50 --p-max-length 0 \
    --o-reads "$OUT/ref-seqs-region.qza"

# 2. Train on the EXTRACTED region. The classifier is a pickled scikit-learn model tied to this QIIME2 release.
#    MEMORY: a full, un-subsampled reference (400K+ sequences) can need tens of GB of RAM and OOM-kill this step;
#    subsample the extracted reference (e.g. ~60,000 sequences) if it does.
qiime feature-classifier fit-classifier-naive-bayes \
    --i-reference-reads "$OUT/ref-seqs-region.qza" \
    --i-reference-taxonomy "$REF_TAX" \
    --o-classifier "$OUT/region-nb-classifier.qza"

# 3. Classify. --p-confidence 0 computes but never truncates; 'disable' skips.
qiime feature-classifier classify-sklearn \
    --i-classifier "$OUT/region-nb-classifier.qza" \
    --i-reads "$REP_SEQS" \
    --p-confidence "${CONFIDENCE:-0.7}" --p-read-orientation auto --p-n-jobs "${N_JOBS:-1}" \
    --o-classification "$OUT/taxonomy.qza"

echo "Done: $OUT/taxonomy.qza (classifier: $OUT/region-nb-classifier.qza)"
