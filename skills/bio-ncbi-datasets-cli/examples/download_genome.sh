#!/bin/bash
# Reference: NCBI Datasets CLI 18.37.0 (checked 2026-09-19) | Verify API if version differs
# Single-assembly download with zip checksum validation + multiple file types.

set -euo pipefail

ACC="${1:-GCF_000001405.40}"   # default: human GRCh38 RefSeq reference
OUT_ZIP="${2:-genome.zip}"
INCLUDE="${3:-genome,gff3,gtf,protein,cds,seq-report}"

echo "=== Download ${ACC} (${INCLUDE}) ==="
datasets download genome accession "${ACC}" \
    --include "${INCLUDE}" \
    --filename "${OUT_ZIP}" \
    --no-progressbar

echo
echo "=== Unzip ==="
DEST="${OUT_ZIP%.zip}"
unzip -q "${OUT_ZIP}" -d "${DEST}/"
ls -lh "${DEST}/ncbi_dataset/data/${ACC}/"

echo
echo "=== Inspect with dataformat (extract metadata from data report JSON) ==="
JSONL="${DEST}/ncbi_dataset/data/assembly_data_report.jsonl"
if [ -f "${JSONL}" ]; then
    dataformat tsv genome \
        --inputfile "${JSONL}" \
        --fields accession,organism-name,assminfo-level,assmstats-scaffold-n50,assmstats-contig-n50,assmstats-total-sequence-len \
        | column -t
fi

echo
echo "datasets download validates the zip checksum -- no manual md5sum step needed."
