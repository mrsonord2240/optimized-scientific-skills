#!/bin/bash
# Reference: samtools 1.24 (ampliconclip needs 1.11+), pysam 0.24.1 for the residual check
# Standard amplicon primer-clipping workflow:
#   input checks -> ampliconclip -> collate -> fixmate -> sort -> calmd -> index -> residual check
# Inputs: input BAM (coord-sorted, indexed), primer BED (strand in col 6),
#         reference FASTA (indexed with samtools faidx)
# CLIP_OPTS (env) picks the clip mode, default "--both-ends --strand":
#   --both-ends --strand   reads can span the whole amplicon (nanopore, HiFi, read length >= amplicon)
#                          and is harmless for short reads that never reach the opposite primer.
#   --strand               short reads only: clips the 5' primer and leaves any 3' primer in place.
# Every check below stops the script with a non-zero exit; a run that prints the final line succeeded.

set -euo pipefail

BAM=${1:?usage: $0 input.bam primers.bed reference.fa output.bam}
PRIMERS=${2:?primer BED required}
REF=${3:?reference FASTA required}
OUT=${4:?output BAM required}
THREADS=${THREADS:-4}
CLIP_OPTS=${CLIP_OPTS:---both-ends --strand}
HERE=$(cd "$(dirname "$0")" && pwd)

fail() { echo "ERROR: $*" >&2; exit 1; }

WORK=$(mktemp -d)
trap 'rm -rf "$WORK"' EXIT

# 0. Input checks: contig names must agree across BAM, BED and FASTA (ARTIC BEDs are on
#    MN908947.3; BAMs aligned to MT192765.1 never match them and would clip nothing).
[ -f "$REF.fai" ] || fail "index the reference first: samtools faidx $REF"
samtools view -H "$BAM" | awk -F'\t' '$1=="@SQ"{sub("SN:","",$2); print $2}' | sort -u > "$WORK/bam.contigs"
awk -F'\t' '!/^#/ && NF{print $1}' "$PRIMERS" | sort -u > "$WORK/bed.contigs"
cut -f1 "$REF.fai" | sort -u > "$WORK/ref.contigs"
comm -12 "$WORK/bam.contigs" "$WORK/bed.contigs" > "$WORK/shared.contigs"
[ -s "$WORK/shared.contigs" ] || fail "no contig shared by BAM header ($(paste -sd, "$WORK/bam.contigs")) and BED ($(paste -sd, "$WORK/bed.contigs"))"
[ -z "$(comm -23 "$WORK/shared.contigs" "$WORK/ref.contigs")" ] || fail "reference FASTA lacks contig(s) the BAM and BED use: $(comm -23 "$WORK/shared.contigs" "$WORK/ref.contigs" | paste -sd, -)"
case "$CLIP_OPTS" in
    *--strand*) awk -F'\t' '!/^#/ && NF && NF<6{bad=1} END{exit bad}' "$PRIMERS" \
        || fail "--strand needs the strand in BED column 6 (chrom start end name score strand); a 5-column BED is rejected" ;;
esac

# 1. Soft-clip primers from BED (--soft-clip is reversible: CIGAR S, bases retained).
#    Stats go to a file so they can be asserted on; TOTAL CLIPPED counts clip events (one per end).
samtools ampliconclip \
    $CLIP_OPTS \
    --soft-clip \
    -f "$WORK/clip.stats" \
    -b "$PRIMERS" \
    "$BAM" \
    -o "$WORK/clipped.bam"
cat "$WORK/clip.stats"
clipped=$(awk '/^TOTAL CLIPPED/{print $3}' "$WORK/clip.stats")
[ "${clipped:-0}" -gt 0 ] || fail "ampliconclip clipped nothing: BED coordinates do not match the reads"

# 2. Repair mate info (CIGARs changed by clipping invalidate MC, ms, TLEN).
#    collate (faster) -> fixmate -m -> coord sort. The sort also fixes the header: the raw
#    ampliconclip output is SO:unknown and cannot be indexed.
samtools collate -O -u "$WORK/clipped.bam" "$WORK/collate" | \
    samtools fixmate -m -u - - | \
    samtools sort -@ "$THREADS" -o "$WORK/sorted.bam" -

# 3. Restore MD/NM (ampliconclip removes them from clipped reads; bcftools mpileup BAQ
#    and IGV mismatch coloring depend on them). calmd stderr stays visible: a wrong
#    reference prints "fail to find sequence" and exits 0 with no MD written.
samtools calmd -@ "$THREADS" -b "$WORK/sorted.bam" "$REF" > "$OUT"
n_map=$(samtools view -c -F 4 "$OUT")
n_md=$(samtools view "$OUT" | awk 'index($0, "\tMD:Z:"){n++} END{print n+0}')
[ "$n_md" -gt 0 ] || fail "calmd wrote no MD tags (wrong reference?)"
echo "records mapped=$n_map with MD=$n_md"

# 4. Index for region access.
samtools index -@ "$THREADS" "$OUT"

# 5. Independent check: no read may still begin (and, with --both-ends, end) inside a primer.
case "$CLIP_OPTS" in *--both-ends*) T3=--three-prime ;; *) T3= ;; esac
if python3 -c 'import pysam' 2>/dev/null; then
    python3 "$HERE/check_primer_residual.py" "$OUT" "$PRIMERS" $T3 || fail "primer bases remain in $OUT"
else
    echo "pysam not importable: skipped the residual-primer check" >&2
fi

echo "Clipped BAM: $OUT"
