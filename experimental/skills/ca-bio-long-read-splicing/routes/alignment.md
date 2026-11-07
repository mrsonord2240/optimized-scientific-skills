# Alignment and microexons

```bash
minimap2 -ax splice:hq --secondary=no --junc-bed annotation.bed12 -t 16 reference.fa isoseq.fastq.gz | samtools sort -o aligned.bam
```

This command is for HiFi. For ONT, direct RNA, preparing annotation BED and unannotated microexons use `references/alignment.md` and `references/microexons.md`. Use -uf only for reads known to be in transcript orientation. For microexons compare inclusion and skipping chains and inspect dangling junctions before increasing junction bonus.

Done when alignment and the requested junctions are checked.
