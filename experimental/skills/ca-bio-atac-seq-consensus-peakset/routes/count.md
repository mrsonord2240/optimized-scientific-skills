# Blacklist filter and count matrix

```bash
bedtools intersect -v -a consensus.bed -b hg38-blacklist.v2.bed.gz | sort -k1,1 -k2,2n > consensus_final.bed
awk 'BEGIN{OFS="\t"; print "GeneID","Chr","Start","End","Strand"} {print $1"_"$2"_"$3, $1, $2+1, $3, "+"}' consensus_final.bed > consensus.saf
featureCounts -F SAF -a consensus.saf -o consensus_counts.tsv -p --countReadPairs -T 8 sample1.bam sample2.bam sample3.bam
```

`iterative_overlap.sh` already filters and writes the SAF; use this for any other consensus. BED start is 0-based, SAF start is 1-based, hence `+1`. `-p --countReadPairs` is for paired-end ATAC (Subread 2.0.2 or later). Pick the blacklist for the same build as the peaks.

- Non-unique SAF GeneIDs mean duplicate peaks: the `chr_start_end` ID needs a non-redundant set.
- Negative coordinates after re-centering: clamp with `bedtools slop -g chrom.sizes -b 0` or `awk '$2 >= 0'`.

Done when the count matrix has one column per sample against the same coordinates.
