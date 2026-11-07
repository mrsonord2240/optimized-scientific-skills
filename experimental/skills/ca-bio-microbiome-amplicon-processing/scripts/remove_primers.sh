#!/usr/bin/env bash
# Checked on cutadapt 4.6+
# Remove primers from paired-end amplicon FASTQ BEFORE DADA2. Primers left on corrupt the error model,
# shift the truncLen frame and masquerade as chimeras.
#
# Usage: FWD=<fwd primer> REV=<rev primer> remove_primers.sh <raw_dir> <out_dir>
# Env:   FWD, REV   required, IUPAC degenerate bases allowed (V4 example: FWD=GTGYCAGCMGCCGCGGTAA REV=GGACTACNVGGGTWTCTAAT)
#        PATTERN_R1 default _R1_001.fastq.gz   PATTERN_R2 default _R2_001.fastq.gz
# Pairs lacking the primer are discarded (--discard-untrimmed).

set -euo pipefail

RAW=${1:-}
OUT=${2:-}
PAT1=${PATTERN_R1:-_R1_001.fastq.gz}
PAT2=${PATTERN_R2:-_R2_001.fastq.gz}

if [ -z "$RAW" ] || [ -z "$OUT" ] || [ -z "${FWD:-}" ] || [ -z "${REV:-}" ]; then
    echo "Usage: FWD=<fwd primer> REV=<rev primer> remove_primers.sh raw_dir out_dir" >&2
    echo "       Both primers are required; use the ones that made the library." >&2
    exit 1
fi
[ -d "$RAW" ] || { echo "ERROR: not a directory: $RAW" >&2; exit 1; }
mkdir -p "$OUT"

n=0
for r1 in "$RAW"/*"$PAT1"; do
    [ -e "$r1" ] || { echo "ERROR: no files matching *$PAT1 in $RAW" >&2; exit 1; }
    r2="${r1/$PAT1/$PAT2}"
    [ -f "$r2" ] || { echo "ERROR: mate not found for $r1: $r2" >&2; exit 1; }
    base=$(basename "$r1" "$PAT1")
    # -g/-G match the forward/reverse primer as 5' adapters on R1/R2
    cutadapt \
        -g "$FWD" \
        -G "$REV" \
        --discard-untrimmed \
        -o "$OUT/${base}${PAT1}" \
        -p "$OUT/${base}${PAT2}" \
        "$r1" "$r2"
    n=$((n + 1))
done
echo "Done: $n sample(s) trimmed into $OUT"
