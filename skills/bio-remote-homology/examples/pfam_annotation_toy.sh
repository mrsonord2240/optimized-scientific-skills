#!/bin/bash
# Reference: HMMER 3.4+ | Toy fixture, no download.
# Same pattern as pfam_annotation.sh, but against one bundled Pfam family instead of the full
# ~1.7 GB Pfam-A.hmm -- use this to sanity-check hmmscan/awk locally before scaling up.
#
# Fixture: examples/data/PF00069.hmm (Pfam Pkinase, GA 31.7, InterPro API, CC0) and
# examples/data/P17612.fasta (UniProt P17612, human PRKACA, CC BY 4.0). Expected result: a single
# Pkinase hit at full-sequence E~1.4e-79 / score~253.2, far above the GA 31.7 gathering threshold.

set -euo pipefail

DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
DATA_DIR="${DIR}/data"
HMM="${DATA_DIR}/PF00069.hmm"
QUERY="${DATA_DIR}/P17612.fasta"

if [ ! -f "${HMM}.h3i" ]; then
    echo "=== One-time: press the bundled PF00069.hmm ==="
    hmmpress "${HMM}"
fi

echo "=== Annotate ${QUERY} against PF00069 with --cut_ga ==="
hmmscan --cut_ga --domtblout query.domtbl --cpu 4 "${HMM}" "${QUERY}" > query.hmmscan.txt

echo
echo "=== Per-protein Pfam domain summary ==="
awk '!/^#/ {print $4"\t"$1"\t"$2"\t"$7"\t"$8}' query.domtbl
echo
echo "Columns: query_name | pfam_name | pfam_acc | full_evalue | full_score"
echo "Expected: sp|P17612|KAPCA_HUMAN  Pkinase  PF00069.32  ~1.4e-79  ~253.2"
