# Cohort two-pass alignment

```bash
STAR --runMode alignReads --runThreadN 8 --genomeDir genome_index --sjdbGTFfile annotation.gtf --sjdbOverhang 149 --readFilesIn sample_R1.fq.gz sample_R2.fq.gz --readFilesCommand zcat --outSAMtype None --outFileNamePrefix pass1_sample_
```

Set sjdbOverhang to the index value and maximum read length minus one. Repeat pass 1 by sample, merge supported novel junctions, then run shared pass 2 with `references/star.md`. XS tags for downstream tools remove unannotated noncanonical junction reads.

Done when all samples are aligned with the shared cohort junctions.
