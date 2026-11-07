## STAR 2-Pass Alignment

**Goal:** Maximize novel-junction sensitivity with one junction reference shared by all samples.

**Approach:** Run STAR once per sample to discover novel junctions (pass 1), merge the novel junctions across the cohort, then re-align every sample with the merged set (pass 2). Per-sample `--twopassMode Basic` is simpler but each sample then inserts its own junctions, so junction sets differ across samples and differential calls do not replicate; use the cohort version for differential splicing. Veeneman 2016 *Bioinformatics* benchmarked per-sample two-pass (>=94% of simulated novel junctions had improved quantification); the cohort-vs-per-sample comparison is the STAR manual's multi-sample 2-pass rationale, not a number from that paper.

```bash
# Pass 1: per sample; only SJ.out.tab is needed
STAR --runMode alignReads \
    --runThreadN 8 \
    --genomeDir genome_index \
    --sjdbGTFfile gencode.v45.basic.gtf \
    --sjdbOverhang 149 \
    --readFilesIn sample_R1.fq.gz sample_R2.fq.gz \
    --readFilesCommand zcat \
    --outSAMtype None \
    --outFileNamePrefix pass1_${sample}_ \
    --outSJtype Standard \
    --outFilterMultimapNmax 20 \
    --alignSJoverhangMin 8 \
    --alignSJDBoverhangMin 3
```

`--sjdbOverhang` = max read length - 1 (149 for 2x150). It must equal the value used at `genomeGenerate`, or STAR stops (Common Errors).

```bash
# Cohort merge. SJ.out.tab columns: 1 chr, 2 start, 3 end, 4 strand (0 undefined, 1 +, 2 -),
# 5 motif (0 non-canonical), 6 annotated (0/1), 7 unique reads, 8 multi-mapped reads, 9 max overhang.
# Keep novel (6==0), canonical-motif junctions with >=3 unique reads; STAR wants only columns 1-4.
cat pass1_*_SJ.out.tab | awk '$6 == 0 && $5 > 0 && $7 >= 3' | cut -f1-4 | sort -u > cohort_novel_SJ.tab

# Pass 2: re-align with the augmented junction set
STAR --runMode alignReads \
    --runThreadN 8 \
    --genomeDir genome_index \
    --sjdbGTFfile gencode.v45.basic.gtf \
    --sjdbFileChrStartEnd cohort_novel_SJ.tab \
    --sjdbOverhang 149 \
    --readFilesIn sample_R1.fq.gz sample_R2.fq.gz \
    --readFilesCommand zcat \
    --outSAMtype BAM SortedByCoordinate \
    --outSAMstrandField intronMotif \
    --outFileNamePrefix pass2_${sample}_ \
    --outSJtype Standard \
    --quantMode GeneCounts \
    --alignSJoverhangMin 8 \
    --alignSJDBoverhangMin 3

samtools index pass2_${sample}_Aligned.sortedByCoord.out.bam   # pysam fetch() needs the index
```

Notes:
- `$5 > 0` alone is not a novelty filter: on the four chrX test samples the merged file had 5,717 lines but 2,184 distinct junctions, and 2,178 of those were already annotated (`sort -u` on whole lines keeps one line per differing read count). `$6 == 0` plus `cut -f1-4` leaves 6 novel junctions.
- `--outSAMstrandField intronMotif` adds the `XS` tag that leafcutter, regtools `-s XS` and Shiba need (0 of 20,000 spliced reads carried it without the flag). It also removes reads with non-canonical unannotated introns, so drop it if you need those reads.
- Pass-2 `SJ.out.tab` marks inserted novel junctions as annotated (col 6 = 1); use pass 1 to tell novel from annotated.
- Raise the unique-read threshold (or require the junction in >=N samples) if the merged file is large enough to hit the `limitSjdbInsertNsj` error.
