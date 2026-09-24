#!/bin/bash
# Reference: NCBI Datasets CLI 18.37.0 (checked 2026-09-19), aria2c 1.36+ | Verify API if version differs
# Bulk pull via --dehydrated + parallel transfer (the cloud / HPC pattern).

set -euo pipefail

TAXON="${1:-Salmonella enterica}"
OUT_ZIP="${2:-bulk_dehydrated.zip}"
DEST="${3:-bulk_dataset}"
THREADS="${4:-8}"

echo "=== Step 1: dehydrated discovery (small ZIP, just metadata) ==="
datasets download genome taxon "${TAXON}" \
    --reference \
    --annotated \
    --assembly-source RefSeq \
    --include genome,gff3,protein \
    --dehydrated \
    --filename "${OUT_ZIP}" \
    --no-progressbar

unzip -q "${OUT_ZIP}" -d "${DEST}/"

FETCH="${DEST}/ncbi_dataset/fetch.txt"
echo "  Files queued: $(wc -l < ${FETCH})"

echo
echo "=== Step 2: parallel transfer with aria2c ==="
# The fetch.txt format (checked on Datasets CLI 18.37.0) is 3 tab-separated columns:
# <url> [TAB] <byte-size-placeholder, always "0"> [TAB] <path-relative-to-data-dir>
# The real path is the 3rd field, not the 2nd -- using $2 here writes "out=0" for
# every row and silently collides every download onto a file literally named "0".
awk -F'\t' '{print $1"\n  out="$3}' "${FETCH}" > "${DEST}/aria2_input.txt"

# Sanity check: fail loudly instead of silently corrupting every filename.
if grep -q '^  out=0$' "${DEST}/aria2_input.txt"; then
    echo "ERROR: aria2_input.txt has an 'out=0' line -- fetch.txt's column layout" >&2
    echo "       changed again; re-check with 'awk -F\"\\t\" \"{print NF}\" ${FETCH}'" >&2
    exit 1
fi

aria2c \
    --input-file="${DEST}/aria2_input.txt" \
    --dir="${DEST}/ncbi_dataset/" \
    --max-concurrent-downloads="${THREADS}" \
    --max-connection-per-server="${THREADS}" \
    --retry-wait=5 \
    --quiet=true

echo
echo "=== Step 3: verify sizes, re-fetch bad files with datasets rehydrate ==="
# `datasets rehydrate` does NOT verify files already on disk (checked on 18.37.0): a blocked or
# truncated aria2c transfer that wrote a file at the right path is reported "already rehydrated".
# So compare every file with the byte length the CLI recorded in dataset_catalog.json, delete
# mismatches, and let rehydrate fetch just those.
bad_files() {
python3 - "${DEST}/ncbi_dataset/data" <<'PY'
import json, os, sys
root = sys.argv[1]
for asm in json.load(open(os.path.join(root, 'dataset_catalog.json')))['assemblies']:
    for f in asm['files']:
        p = root + '/' + f['filePath']
        if not os.path.exists(p) or os.path.getsize(p) != int(f['uncompressedLengthBytes']):
            print(p)
PY
}

BAD="$(bad_files | tr -d '\r')"
if [ -n "${BAD}" ]; then
    echo "  $(echo "${BAD}" | wc -l) file(s) missing or wrong size; deleting and rehydrating"
    echo "${BAD}" | while IFS= read -r f; do rm -f "$f"; done
    datasets rehydrate --directory "${DEST}/" --max-workers "${THREADS}" --no-progressbar
    BAD="$(bad_files | tr -d '\r')"
fi
if [ -n "${BAD}" ]; then
    echo "ERROR: still wrong after rehydrate:" >&2
    echo "${BAD}" >&2
    exit 1
fi
echo "  all files match dataset_catalog.json sizes"

echo
echo "=== Done ==="
ls "${DEST}/ncbi_dataset/data/" | head
echo
echo "Files: $(find ${DEST}/ncbi_dataset/data/ -type f | wc -l)"
echo "Total size: $(du -sh ${DEST}/ncbi_dataset/data/ | cut -f1)"
