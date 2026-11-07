#!/bin/bash
# Reference: samtools 1.24 (ampliconclip needs 1.11+), pysam 0.24.1 for the residual check
# Standard amplicon primer-clipping workflow:
#   input checks -> ampliconclip -> collate -> fixmate -> sort -> calmd -> index -> residual check
# Inputs: input BAM (coord-sorted, indexed), primer BED (strand in col 6),
#         reference FASTA (indexed with samtools faidx)
# MAX_NOT_CLIPPED_PCT (env, default 1) is the largest accepted share of mapped reads
# that ampliconclip may leave unmatched. Raise it only for a documented sparse/negative control.
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
MAX_NOT_CLIPPED_PCT=${MAX_NOT_CLIPPED_PCT:-1}
HERE=$(cd "$(dirname "$0")" && pwd)

fail() { echo "ERROR: $*" >&2; exit 1; }
[ ! -e "$OUT" ] && [ ! -e "$OUT.bai" ] \
    || fail "output already exists ($OUT or $OUT.bai); choose a new path or remove it explicitly"

WORK=$(mktemp -d)
trap 'rm -rf "$WORK"' EXIT

# 0. Input checks: contig names must agree across BAM, BED and FASTA (ARTIC BEDs are on
#    MN908947.3; BAMs aligned to MT192765.1 never match them and would clip nothing).
[ -f "$REF.fai" ] || fail "index the reference first: samtools faidx $REF"
python3 -c 'import pysam' 2>/dev/null || fail "pysam is required for the residual-primer check"
# Normalize samtools-valid whitespace-delimited BED while ignoring UCSC track/browser headers.
awk 'BEGIN{OFS="\t"}
     {sub(/\r$/, "")}
     NF==0 || $1 ~ /^#/ || $1=="track" || $1=="browser" {next}
     {for(i=1;i<=NF;i++) printf "%s%s", $i, (i==NF ? ORS : OFS)}' \
    "$PRIMERS" > "$WORK/primers.bed"
[ -s "$WORK/primers.bed" ] || fail "primer BED has no data rows after comments/track/browser headers"
samtools view -H "$BAM" | awk -F'\t' '$1=="@SQ"{sub("SN:","",$2); print $2}' | sort -u > "$WORK/bam.contigs"
awk '{print $1}' "$WORK/primers.bed" | sort -u > "$WORK/bed.contigs"
cut -f1 "$REF.fai" | sort -u > "$WORK/ref.contigs"
comm -12 "$WORK/bam.contigs" "$WORK/bed.contigs" > "$WORK/shared.contigs"
[ -s "$WORK/shared.contigs" ] || fail "no contig shared by BAM header ($(paste -sd, "$WORK/bam.contigs")) and BED ($(paste -sd, "$WORK/bed.contigs"))"
[ -z "$(comm -23 "$WORK/shared.contigs" "$WORK/ref.contigs")" ] || fail "reference FASTA lacks contig(s) the BAM and BED use: $(comm -23 "$WORK/shared.contigs" "$WORK/ref.contigs" | paste -sd, -)"
case "$CLIP_OPTS" in
    *--strand*) awk 'NF<6 || $6 !~ /^[+-]$/{bad=1} END{exit bad}' "$WORK/primers.bed" \
        || fail "primer BED problem: --strand needs >=6 columns and + or - in column 6" ;;
esac

# 1. Soft-clip primers from BED (--soft-clip is reversible: CIGAR S, bases retained).
#    Stats go to a file so they can be asserted on; TOTAL CLIPPED counts clip events (one per end).
samtools ampliconclip \
    $CLIP_OPTS \
    --soft-clip \
    -f "$WORK/clip.stats" \
    -b "$WORK/primers.bed" \
    "$BAM" \
    -o "$WORK/clipped.bam"
cat "$WORK/clip.stats"
clipped=$(awk '/^TOTAL CLIPPED/{print $3}' "$WORK/clip.stats")
[ "${clipped:-0}" -gt 0 ] || fail "ampliconclip clipped nothing (wrong BED/build, or no reads overlap primers)"
read -r total_reads not_clipped < <(awk '
  $1=="TOTAL" && $2=="READS:" {total=$3}
  $1=="NOT" && $2=="CLIPPED:" {not=$3}
  END {print total+0, not+0}' "$WORK/clip.stats")
[ "$total_reads" -gt 0 ] || fail "ampliconclip reported no reads (empty/sparse input or wrong BED/build)"
not_clipped_pct=$(awk -v n="$not_clipped" -v t="$total_reads" 'BEGIN{printf "%.3f", 100*n/t}')
awk -v p="$not_clipped_pct" -v m="$MAX_NOT_CLIPPED_PCT" 'BEGIN{exit !(p<=m)}' \
    || fail "NOT CLIPPED share ${not_clipped_pct}% exceeds ${MAX_NOT_CLIPPED_PCT}% (wrong primer scheme/build, or a sparse/negative-control sample; review clip.stats and set MAX_NOT_CLIPPED_PCT only with justification)"

# 2. Repair mate info (CIGARs changed by clipping invalidate MC, ms, TLEN).
#    collate (faster) -> fixmate -m -> coord sort. The sort also fixes the header: the raw
#    ampliconclip output is SO:unknown and cannot be indexed.
samtools collate -O -u "$WORK/clipped.bam" "$WORK/collate" | \
    samtools fixmate -m -u - - | \
    samtools sort -@ "$THREADS" -o "$WORK/sorted.bam" -

# 3. Restore MD/NM (ampliconclip removes them from clipped reads; IGV mismatch coloring
#    and NM/MD-based filters read them). calmd stderr stays visible: a wrong
#    reference prints "fail to find sequence" and exits 0 with no MD written.
samtools calmd -@ "$THREADS" -b "$WORK/sorted.bam" "$REF" > "$WORK/final.bam"
n_map=$(samtools view -c -F 4 "$WORK/final.bam")
n_md=$(samtools view "$WORK/final.bam" | awk 'index($0, "\tMD:Z:"){n++} END{print n+0}')
[ "$n_md" -gt 0 ] || fail "calmd wrote no MD tags (wrong reference?)"
echo "records mapped=$n_map with MD=$n_md"

# 4. Index for region access.
samtools index -@ "$THREADS" "$WORK/final.bam"

# 5. Independent check: no read may still begin (and, with --both-ends, end) inside a primer.
case "$CLIP_OPTS" in *--both-ends*) T3=--three-prime ;; *) T3= ;; esac
set +e
python3 "$HERE/check_primer_residual.py" "$WORK/final.bam" "$WORK/primers.bed" $T3
check_rc=$?
set -e
case "$check_rc" in
    0) ;;
    1) fail "primer bases remain after clipping" ;;
    2) fail "residual-primer check could not run because an input was invalid" ;;
    *) fail "residual-primer check failed unexpectedly (exit $check_rc)" ;;
esac

# Publish only after every assertion passes, so failed runs never leave a finished-looking BAM.
mv -f "$WORK/final.bam" "$OUT"
mv -f "$WORK/final.bam.bai" "$OUT.bai"

echo "Clipped BAM: $OUT"
