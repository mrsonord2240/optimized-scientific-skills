# Duplicate Handling: other markers (samblaster, Picard, biobambam2, sambamba, mapDamage, pbmarkdup)

## Alternative: From Aligner

Some aligners can mark duplicates directly during streaming:

### BWA-MEM2 with samblaster
```bash
# -R sets the @RG line that Picard MarkDuplicates needs (see Common Errors)
bwa-mem2 mem -R '@RG\tID:s1\tSM:s1\tLB:lib1\tPL:ILLUMINA' ref.fa R1.fq R2.fq | \
    samblaster | \
    samtools sort -o marked.bam
```

### Picard MarkDuplicates
```bash
java -jar picard.jar MarkDuplicates \
    I=input.bam \
    O=marked.bam \
    M=metrics.txt \
    OPTICAL_DUPLICATE_PIXEL_DISTANCE=2500
```
Picard reports READ_PAIR_DUPLICATES in pairs (samtools counts reads: 50 pairs = 100 reads). Input does not need fixmate.

### biobambam2, sambamba
```bash
# Coordinate-sorted input; no fixmate step. Metrics use Picard's column names
bammarkduplicates2 I=input.bam O=marked.bam M=metrics.txt
sambamba markdup -t 4 input.bam marked.bam
```

### Ancient DNA: mapDamage rescale after markdup
```bash
# Needs mapDamage 2.2.x (2.1.x rejects paired-end BAMs). The Bayesian rescaling step takes minutes even on small BAMs.
mapDamage -i marked.bam -r ref.fa --rescale -d mapdamage_out   # writes mapdamage_out/marked.rescaled.bam
```

### pbmarkdup (PacBio HiFi amplicons, unaligned BAM/FASTQ)
```bash
# Writes duplicate flags (0x400) into the BAM; --rmdup drops them, --dup-file keeps them in a separate file
pbmarkdup -j 4 hifi.bam marked.bam
```
