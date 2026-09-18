#!/bin/bash
# Reference: sra-tools 3.0+, ENA portal API 2.0+ | Verify API if version differs
# Batch SRA download via the ENA mirror (preferred default in 2026): direct FASTQ + MD5 verification.

set -euo pipefail

ACC_FILE="${1:-accessions.txt}"
OUT="${2:-./fastq}"

if [ ! -f "${ACC_FILE}" ]; then
    echo "Usage: $0 <accessions.txt> [out_dir]"
    echo "  one SRR/ERR/DRR per line; '#' starts a comment"
    exit 1
fi

mkdir -p "${OUT}"
FAILED="${OUT}/failed.txt"
: > "${FAILED}"

total=$(grep -cv '^[[:space:]]*#\|^[[:space:]]*$' "${ACC_FILE}")
count=0
ok=0

while read -r ACC; do
    [[ -z "${ACC}" || "${ACC}" == \#* ]] && continue
    count=$((count+1))
    echo
    echo "[${count}/${total}] ${ACC}"

    # Query ENA portal API for FASTQ URLs + md5. Locate columns by their documented
    # field name, not a fixed index: ENA's filereport ALWAYS prepends run_accession as
    # column 1 regardless of what `fields=` lists, so with fields=fastq_ftp,fastq_md5,
    # read_count the real data is in columns 2-4, not 1-3 (cut -f1/-f2 silently grabbed
    # the accession and the URL, one column short -- confirmed against the live API;
    # fails 100% of the time as a fixed-index cut).
    RESPONSE=$(curl -fsS "https://www.ebi.ac.uk/ena/portal/api/filereport?accession=${ACC}&result=read_run&fields=fastq_ftp,fastq_md5,read_count&format=tsv") || {
        echo "  ENA portal API failed for ${ACC}"
        echo "${ACC}" >> "${FAILED}"
        continue
    }
    HEADER=$(echo "${RESPONSE}" | head -1)
    ROW=$(echo "${RESPONSE}" | tail -1)
    FTP_COL=$(echo "${HEADER}" | tr '\t' '\n' | grep -nx 'fastq_ftp' | cut -d: -f1)
    MD5_COL=$(echo "${HEADER}" | tr '\t' '\n' | grep -nx 'fastq_md5' | cut -d: -f1)

    URLS=$(echo "${ROW}" | cut -f"${FTP_COL}" | tr ';' '\n')
    MD5S=$(echo "${ROW}" | cut -f"${MD5_COL}" | tr ';' '\n')
    if [ -z "${URLS}" ]; then
        echo "  No FASTQ URLs in ENA mirror for ${ACC}"
        echo "${ACC}" >> "${FAILED}"
        continue
    fi

    success=true
    i=0
    while read -r url; do
        i=$((i+1))
        expected_md5=$(echo "${MD5S}" | sed -n "${i}p")
        fname="${OUT}/$(basename ${url})"
        echo "  Downloading $(basename ${url})"
        if ! curl -fsSL -o "${fname}" "https://${url}"; then
            echo "  curl failed"
            success=false
            break
        fi
        actual_md5=$(md5sum "${fname}" | awk '{print $1}')
        if [ "${actual_md5}" != "${expected_md5}" ]; then
            echo "  MD5 MISMATCH: expected ${expected_md5}, got ${actual_md5}"
            success=false
            break
        fi
        echo "  md5 OK"
    done <<< "${URLS}"

    if [ "${success}" = "true" ]; then
        ok=$((ok+1))
    else
        echo "${ACC}" >> "${FAILED}"
    fi
done < "${ACC_FILE}"

echo
echo "=== Summary ==="
echo "OK:     ${ok}/${total}"
echo "Failed: $((total - ok)) (see ${FAILED})"
