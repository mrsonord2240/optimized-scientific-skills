#!/bin/bash
# Checked on samtools 1.24
# Merge BAMs after checking they share one sort order. samtools merge itself does not check.
#
# Usage: merge_bams.sh <output.bam> <in1.bam> <in2.bam> [...]
# Env:   ORDER    coordinate (default) | name-natural | name-ascii; must match how the inputs were sorted
#        THREADS  default 4

set -euo pipefail

OUTPUT=${1:-}
if [ -z "$OUTPUT" ] || [ "$#" -lt 3 ]; then
    echo "Usage: [ORDER=coordinate|name-natural|name-ascii] merge_bams.sh out.bam in1.bam in2.bam [...]" >&2
    exit 1
fi
shift
THREADS=${THREADS:-4}
ORDER=${ORDER:-coordinate}

case "$ORDER" in
    coordinate)   FLAGS=() ;;
    name-natural) FLAGS=(-n) ;;
    name-ascii)   FLAGS=(-N) ;;
    *) echo "ERROR: ORDER must be coordinate, name-natural or name-ascii" >&2; exit 2 ;;
esac

for f in "$@"; do
    [ -f "$f" ] || { echo "ERROR: input not found: $f" >&2; exit 1; }
done

# Exactly one distinct @HD line (name-sorted: -n and -N inputs print different SS: lines; do not mix them)
N_HD=$(for f in "$@"; do samtools view -H "$f" | grep -m1 '^@HD' || echo "@HD missing"; done | sort -u | tee /dev/stderr | wc -l)
if [ "$N_HD" -ne 1 ]; then
    echo "ERROR: inputs disagree on @HD sort order (above); merge would exit 0 and write records out of order. Re-sort them the same way first." >&2
    exit 3
fi

# -c dedups @RG by ID, -p dedups @PG. Two different read groups sharing an ID collapse into one.
samtools merge "${FLAGS[@]}" -c -p -@ "$THREADS" -o "$OUTPUT" "$@"
samtools quickcheck -v "$OUTPUT" || { echo "ERROR: merged output failed quickcheck" >&2; exit 3; }

if [ "$ORDER" = coordinate ]; then
    samtools index "$OUTPUT"
fi
echo "Done: $OUTPUT"
