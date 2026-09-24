#!/bin/bash
# fgbio UMI consensus calling from a BAM with the UMI in the RX tag.
#
# Usage:  fgbio_consensus.sh single in.bam out.bam   # single-strand molecular consensus (RX: one UMI or UMI1-UMI2), adjacency grouping
#         fgbio_consensus.sh duplex in.bam out.bam   # duplex (xGen-Prism, NEBNext duplex): needs RX as UMI1-UMI2, --strategy=paired
# Needs:  fgbio, samtools. Checked on fgbio 4.1.1, samtools 1.24.
# Notes:  If the UMI is in a separate FASTQ instead of the RX tag, annotate first and pass the annotated BAM:
#           fgbio AnnotateBamWithUmis -i raw.bam -f umi.fastq -o annotated.bam
#         The mate MQ tag is added with samtools fixmate -m (alternative: fgbio SetMateInformation -i qn.bam -o mated.bam).
#         Consensus reads are written UNMAPPED: re-align them before variant calling.
#         CallDuplexConsensusReads on adjacency-grouped reads crashes (StringIndexOutOfBoundsException), and
#         --strategy=paired on a single-UMI RX fails (IllegalArgumentException), so the mode picks the strategy.

set -euo pipefail

MODE=${1:-}
IN=${2:-}
OUT=${3:-}
if [ -z "$MODE" ] || [ -z "$IN" ] || [ -z "$OUT" ]; then
    echo "Usage: fgbio_consensus.sh single|duplex in.bam out.bam" >&2
    exit 1
fi

TMP=$(mktemp -d)
trap 'rm -rf "$TMP"' EXIT

samtools sort -n -o "$TMP/qn.bam" "$IN"
samtools fixmate -m "$TMP/qn.bam" "$TMP/mated.bam"        # or: fgbio SetMateInformation -i qn.bam -o mated.bam

case "$MODE" in
    single)
        # Single-strand molecular consensus
        fgbio GroupReadsByUmi -i "$TMP/mated.bam" -o "$TMP/grouped.bam" --strategy=adjacency --edits=1 --raw-tag=RX
        fgbio CallMolecularConsensusReads -i "$TMP/grouped.bam" -o "$OUT" --min-reads=1
        ;;
    duplex)
        # Duplex: --strategy=paired writes MI tags with /A /B strand suffixes
        fgbio GroupReadsByUmi -i "$TMP/mated.bam" -o "$TMP/grouped_duplex.bam" --strategy=paired --edits=1 --raw-tag=RX
        fgbio CallDuplexConsensusReads -i "$TMP/grouped_duplex.bam" -o "$OUT" --min-reads 1 1 0
        ;;
    *)
        echo "ERROR: mode must be single or duplex" >&2
        exit 1 ;;
esac

n_out=$(samtools view -c "$OUT")
echo "consensus records: $n_out"
test "$n_out" -gt 0
