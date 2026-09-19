#!/bin/bash
# Reference: HMMER 3.4+, Pfam 36+ | Verify API if version differs
# Canonical Pfam-A domain annotation with --cut_ga (calibrated gathering thresholds).

set -euo pipefail

QUERY="${1:-proteins.fa}"
PFAM_DIR="${2:-pfam_db}"

mkdir -p "${PFAM_DIR}"

if [ ! -f "${PFAM_DIR}/Pfam-A.hmm.h3i" ]; then
    echo "=== One-time: download and press Pfam-A ==="
    cd "${PFAM_DIR}"
    wget -q https://ftp.ebi.ac.uk/pub/databases/Pfam/current_release/Pfam-A.hmm.gz
    gunzip Pfam-A.hmm.gz
    hmmpress Pfam-A.hmm
    cd -
fi

echo "=== Annotate ${QUERY} against Pfam-A with --cut_ga ==="
# --cut_ga uses each model's pre-calibrated gathering threshold.
# This is the Pfam-recommended cutoff -- preferred over arbitrary E-value.
hmmscan --cut_ga \
        --domtblout query.domtbl \
        --tblout query.tbl \
        --cpu 8 \
        "${PFAM_DIR}/Pfam-A.hmm" "${QUERY}" > query.hmmscan.txt

echo
echo "=== Per-protein Pfam domain summary ==="
# domtbl columns (1-indexed): 1 target_name, 2 target_acc, 3 tlen, 4 query_name, 5 query_acc,
#                 6 qlen, 7 full_evalue, 8 full_score, 9 full_bias, 10 dom#, 11 of, 12 c-Evalue,
#                 13 i-Evalue, 14 dom_score, 15 dom_bias, 16-21 hmm/ali/env coords, 22 acc, 23+ desc
awk '!/^#/ {print $4"\t"$1"\t"$2"\t"$7"\t"$8}' query.domtbl | sort -k1,1 -k4,4g | head -20
echo
echo "Columns: query_name | pfam_name | pfam_acc | full_evalue | full_score"
