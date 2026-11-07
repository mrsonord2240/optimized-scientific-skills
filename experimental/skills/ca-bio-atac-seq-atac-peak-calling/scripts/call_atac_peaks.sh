#!/bin/bash
# ENCODE-style ATAC-seq replicate reproducibility: MACS peak calls, disjoint pseudoreplicates, IDR, rescue and
# self-consistency ratios, blacklist-filtered conservative peak set, optional bigWig.
# Tested: MACS3 3.0.4, samtools 1.24, bedtools 2.31.1, idr 2.0.4.2 (needs numpy<1.24), UCSC bedGraphToBigWig 482.
# Differences from the ENCODE pipeline: reads BAM (not Tn5-shifted tagAlign) and does not use --call-summits.
# `-f BAM --shift -75 --extsize 150` treats each read end as a cut; `-f BAMPE` would ignore --shift/--extsize.
#
# Usage: call_atac_peaks.sh REP1.bam REP2.bam GENOME_SIZE BLACKLIST.bed [OUTDIR] [CHROM.sizes]
#   GENOME_SIZE  effective size matched to read length, e.g. 2.806e9 (hg38, 100 bp; see usage-guide.md)
#   MACS=macs2 or IDR=/path/to/idr may be set in the environment (defaults: macs3, idr).
# The two BAMs must be coordinate-sorted, indexed, deduplicated, MAPQ-filtered and chrM-free.
# Outputs: idr/reproducibility.txt (Nt N1 N2 Np, ratios, verdict); final/conservative.narrowPeak (IDR output cut to
# 10 columns: column 5 is the IDR score, columns 7 and 9 are -1); final/pooled.bw when CHROM.sizes is given.

set -euo pipefail

die() { echo "ERROR: $*" >&2; exit 1; }

REP1=${1:?usage: $0 REP1.bam REP2.bam GENOME_SIZE BLACKLIST.bed [OUTDIR] [CHROM.sizes]}
REP2=${2:?missing REP2.bam}
GENOME=${3:?missing GENOME_SIZE (effective size, e.g. 2.806e9)}
BLACKLIST=${4:?missing BLACKLIST.bed}
OUTDIR=${5:-peaks_out}
CHROMSIZES=${6:-}
MACS=${MACS:-macs3}
IDR=${IDR:-idr}
IDR_T=0.05   # one IDR threshold for every comparison, as in the ENCODE pipeline

for tool in "$MACS" "$IDR" samtools bedtools awk sort; do
    command -v "$tool" >/dev/null || die "$tool not on PATH"
done
[[ -r $BLACKLIST ]] || die "blacklist not readable: $BLACKLIST"
for bam in "$REP1" "$REP2"; do
    [[ -r $bam ]] || die "BAM not readable: $bam"
    [[ -f $bam.bai || -f ${bam%.bam}.bai || -f $bam.csi ]] || die "BAM not indexed (samtools index): $bam"
    samtools view -H "$bam" | grep -q 'SO:coordinate' || die "BAM not coordinate-sorted: $bam"
    n=$(samtools idxstats "$bam" | awk '$1 ~ /^(chrM|chrMT|MT|M)$/ {s += $3} END {print s + 0}')
    [[ $n -eq 0 ]] || die "$bam has $n reads on the mitochondrial contig; remove chrM first"
done
[[ -z $CHROMSIZES || -r $CHROMSIZES ]] || die "chrom sizes not readable: $CHROMSIZES"

mkdir -p "$OUTDIR"/{rep1,rep2,pooled,psr1_1,psr1_2,psr2_1,psr2_2,poolpsr_1,poolpsr_2,idr,final}

call_peaks() {   # call_peaks NAME OUTDIR EXTRA_FLAGS... -- BAM...
    local name=$1 outdir=$2; shift 2
    local extra=() bams=()
    while [[ $1 != -- ]]; do extra+=("$1"); shift; done
    shift; bams=("$@")
    "$MACS" callpeak -t "${bams[@]}" -f BAM -g "$GENOME" -n "$name" --outdir "$outdir" \
        --nomodel --shift -75 --extsize 150 --keep-dup all -p 0.01 "${extra[@]}"
}

# 1. Per-replicate and pooled peaks (pooled also writes the signal track)
call_peaks rep1 "$OUTDIR/rep1" -- "$REP1"
call_peaks rep2 "$OUTDIR/rep2" -- "$REP2"
call_peaks pooled "$OUTDIR/pooled" -B --SPMR -- "$REP1" "$REP2"

# 2. Disjoint pseudoreplicates: -s selects reads by read-name hash, -U writes the complement, so every
#    read (pair) lands in exactly one half.
samtools view -b -s 1.5 -U "$OUTDIR/psr1_2/rep1.psr2.bam" -o "$OUTDIR/psr1_1/rep1.psr1.bam" "$REP1"
samtools view -b -s 1.5 -U "$OUTDIR/psr2_2/rep2.psr2.bam" -o "$OUTDIR/psr2_1/rep2.psr1.bam" "$REP2"
samtools merge -f "$OUTDIR/poolpsr_1/pool.psr1.bam" "$OUTDIR/psr1_1/rep1.psr1.bam" "$OUTDIR/psr2_1/rep2.psr1.bam"
samtools merge -f "$OUTDIR/poolpsr_2/pool.psr2.bam" "$OUTDIR/psr1_2/rep1.psr2.bam" "$OUTDIR/psr2_2/rep2.psr2.bam"

call_peaks rep1_psr1 "$OUTDIR/psr1_1" -- "$OUTDIR/psr1_1/rep1.psr1.bam"
call_peaks rep1_psr2 "$OUTDIR/psr1_2" -- "$OUTDIR/psr1_2/rep1.psr2.bam"
call_peaks rep2_psr1 "$OUTDIR/psr2_1" -- "$OUTDIR/psr2_1/rep2.psr1.bam"
call_peaks rep2_psr2 "$OUTDIR/psr2_2" -- "$OUTDIR/psr2_2/rep2.psr2.bam"
call_peaks pool_psr1 "$OUTDIR/poolpsr_1" -- "$OUTDIR/poolpsr_1/pool.psr1.bam"
call_peaks pool_psr2 "$OUTDIR/poolpsr_2" -- "$OUTDIR/poolpsr_2/pool.psr2.bam"

# 3. IDR (ranked by p-value, column 8) for Nt, N1, N2 and Np
sort_pk() { sort -k8,8nr "$1" > "${1%.narrowPeak}.sorted.narrowPeak"; }
run_idr() {   # run_idr PEAKS_A PEAKS_B LABEL
    sort_pk "$1"; sort_pk "$2"
    "$IDR" --samples "${1%.narrowPeak}.sorted.narrowPeak" "${2%.narrowPeak}.sorted.narrowPeak" \
        --input-file-type narrowPeak --rank p.value \
        --output-file "$OUTDIR/idr/$3.idr" --idr-threshold "$IDR_T" --plot \
        --log-output-file "$OUTDIR/idr/$3.log"
}
run_idr "$OUTDIR/rep1/rep1_peaks.narrowPeak" "$OUTDIR/rep2/rep2_peaks.narrowPeak" true_reps
run_idr "$OUTDIR/psr1_1/rep1_psr1_peaks.narrowPeak" "$OUTDIR/psr1_2/rep1_psr2_peaks.narrowPeak" rep1_pseudoreps
run_idr "$OUTDIR/psr2_1/rep2_psr1_peaks.narrowPeak" "$OUTDIR/psr2_2/rep2_psr2_peaks.narrowPeak" rep2_pseudoreps
run_idr "$OUTDIR/poolpsr_1/pool_psr1_peaks.narrowPeak" "$OUTDIR/poolpsr_2/pool_psr2_peaks.narrowPeak" pooled_pseudoreps

# 4. Rescue and self-consistency ratios (ENCODE: pass if both <= 2, borderline if one > 2, fail if both > 2)
idr_n() { awk -v t="$IDR_T" 'BEGIN {m = int(-125 * log(t) / log(2))} $5 >= m' "$OUTDIR/idr/$1.idr" | wc -l; }
NT=$(idr_n true_reps); N1=$(idr_n rep1_pseudoreps); N2=$(idr_n rep2_pseudoreps); NP=$(idr_n pooled_pseudoreps)
ratio() { awk -v a="$1" -v b="$2" 'BEGIN {m = (a > b) ? a : b; n = (a < b) ? a : b; if (n > 0) printf "%.3f", m / n; else print "inf"}'; }
RESCUE=$(ratio "$NP" "$NT"); SELF=$(ratio "$N1" "$N2")
VERDICT=$(awk -v r="$RESCUE" -v s="$SELF" 'BEGIN {
    br = (r == "inf" || r > 2); bs = (s == "inf" || s > 2)
    print (br && bs) ? "FAIL" : ((br || bs) ? "BORDERLINE" : "PASS")}')
echo "Nt=$NT N1=$N1 N2=$N2 Np=$NP rescue=$RESCUE self_consistency=$SELF verdict=$VERDICT (IDR <= $IDR_T)" \
    | tee "$OUTDIR/idr/reproducibility.txt"

# 5. Conservative set: true-replicate IDR peaks, narrowPeak columns, blacklist removed
awk -v t="$IDR_T" 'BEGIN {m = int(-125 * log(t) / log(2))} $5 >= m' "$OUTDIR/idr/true_reps.idr" \
    | cut -f1-10 \
    | bedtools intersect -v -a stdin -b "$BLACKLIST" > "$OUTDIR/final/conservative.narrowPeak"
echo "conservative peaks after blacklist: $(wc -l < "$OUTDIR/final/conservative.narrowPeak")"

# 6. Optional bigWig from the pooled SPMR pileup (bedGraphToBigWig cannot read a pipe)
if [[ -n $CHROMSIZES ]]; then
    command -v bedGraphToBigWig >/dev/null || die "bedGraphToBigWig not on PATH"
    sort -k1,1 -k2,2n "$OUTDIR/pooled/pooled_treat_pileup.bdg" > "$OUTDIR/final/pooled.sorted.bdg"
    bedGraphToBigWig "$OUTDIR/final/pooled.sorted.bdg" "$CHROMSIZES" "$OUTDIR/final/pooled.bw"
fi
