#!/bin/bash
# Reference: picard 3.1+, pysam 0.22+, samtools 1.19+ | Verify API if version differs
# Checked on samtools 1.24
# Mark duplicates using the samtools collate | fixmate | sort | markdup pipeline.
#
# Usage:   ASSAY=<assay> markdup_pipeline.sh <input.bam> <output.bam> [threads]
# Env:     ASSAY         required; see SKILL.md "When to Mark Duplicates". Refuses the assays where markdup is wrong.
#          OPTICAL_DIST  samtools markdup -d (default 2500; use 100 for HiSeq 2000/2500, MiSeq, NextSeq 500/550)

set -euo pipefail

INPUT=${1:-}
OUTPUT=${2:-}
THREADS=${3:-4}
OPTICAL_DIST=${OPTICAL_DIST:-2500}

if [ -z "$INPUT" ] || [ -z "$OUTPUT" ]; then
    echo "Usage: ASSAY=<assay> markdup_pipeline.sh <input.bam> <output.bam> [threads]" >&2
    exit 1
fi

# Step 0: gate on assay. Standard markdup is wrong for these (SKILL.md decision table).
case "${ASSAY:-}" in
    wgs|wes|capture|somatic|chipseq|cutrun|atac|adna|pacbio-amplicon) ;;
    rnaseq|scrna|umi|amplicon|longread|16s|its)
        echo "ERROR: ASSAY=$ASSAY: samtools markdup is the wrong tool. See SKILL.md 'When to Mark Duplicates' for the right one." >&2
        exit 2 ;;
    *)
        echo "ERROR: set ASSAY to one of: wgs wes capture somatic chipseq cutrun atac adna pacbio-amplicon" >&2
        echo "       (rnaseq scrna umi amplicon longread 16s its are refused; see SKILL.md 'When to Mark Duplicates')" >&2
        exit 2 ;;
esac

# Step 0b: ASSAY is the caller's word, so also check the BAM. A splice-aware aligner in @PG, or N in the CIGARs
# of the first 100000 records, means RNA-seq: refuse. ALLOW_SPLICED=1 overrides (e.g. a DNA BAM realigned by STAR).
if [ "${ALLOW_SPLICED:-0}" != 1 ]; then
    n_pg=$(samtools view -H "$INPUT" | grep -Eic '^@PG.*(ID|PN):(STAR|HISAT2)' || true)
    read -r n_spliced n_sampled < <(samtools view "$INPUT" 2>/dev/null | awk 'NR>100000{exit} $6 ~ /N/{n++} END{print n+0, NR+0}' || true)
    if [ "$n_pg" -gt 0 ] || [ $((n_spliced * 100)) -gt "$n_sampled" ]; then
        echo "ERROR: $INPUT looks like RNA-seq (splice-aware aligner in @PG: $n_pg; spliced CIGARs: $n_spliced of $n_sampled records)." >&2
        echo "       samtools markdup is the wrong tool for it (SKILL.md 'When to Mark Duplicates'). Set ALLOW_SPLICED=1 to override." >&2
        exit 2
    fi
fi

echo "Marking duplicates in $INPUT (assay $ASSAY, $THREADS threads, -d $OPTICAL_DIST)..."

TMP=$(mktemp -d)
trap 'rm -rf "$TMP"' EXIT
STATS="${OUTPUT%.bam}.markdup_stats.txt"

# Stats go to a temp path first: markdup opens -f for writing as soon as the pipe starts, before an
# upstream failure (e.g. a missing input) reaches it, so writing straight to $STATS left an empty file
# behind on a failed run. Move it beside $OUTPUT only once the whole pipeline has actually succeeded.
samtools collate -O -u -@ "$THREADS" "$INPUT" "$TMP/collate" | \
    samtools fixmate -m -u -@ "$THREADS" - - | \
    samtools sort -u -@ "$THREADS" -T "$TMP/sort" - | \
    samtools markdup -@ "$THREADS" -d "$OPTICAL_DIST" --use-read-groups -f "$TMP/markdup_stats.txt" - "$TMP/marked.bam"

# markdup never drops records (no -r): the output must have exactly as many as the input.
n_in=$(samtools view -c "$INPUT")
n_out=$(samtools view -c "$TMP/marked.bam")
if [ "$n_in" -eq 0 ] || [ "$n_in" -ne "$n_out" ]; then
    echo "ERROR: record count changed or empty: input $n_in, output $n_out" >&2
    exit 3
fi

mv "$TMP/marked.bam" "$OUTPUT"
mv "$TMP/markdup_stats.txt" "$STATS"
echo "Indexing..."
samtools index "$OUTPUT"

echo "Done: $OUTPUT (stats: $STATS)"
echo ""
echo "Duplicate statistics:"
samtools flagstat "$OUTPUT" | grep -E "(total|duplicates)"

n_dup=$(samtools view -c -f 1024 -F 2304 "$OUTPUT")
n_pri=$(samtools view -c -F 2304 "$OUTPUT")
awk -v d="$n_dup" -v t="$n_pri" 'BEGIN{printf "Flagged as duplicates: %d of %d primary reads (%.2f%%)\n", d, t, (t ? d * 100 / t : 0)}'
if [ $((n_dup * 2)) -gt "$n_pri" ]; then
    echo "WARNING: over 50% of primary reads flagged as duplicates; confirm ASSAY=$ASSAY is right (SKILL.md 'When to Mark Duplicates')." >&2
fi
