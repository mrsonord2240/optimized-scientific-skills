# Alternative clipping

```bash
ivar trim -i input.bam -b primers.bed -p ivar_trimmed -q 0 -m 1
```

iVar expects coordinate-sorted input and drops reads without a primer unless `-e` is supplied. The command disables its quality and minimum-length filters; its default quality threshold is 20. For nanopore use `-q 0 -m 1`, inspect residuals and repair mate fields, sorting, MD/NM and index as in `routes/clip.md`. For irreversible archival hard clipping, use `samtools ampliconclip --both-ends --strand --hard-clip -b primers.bed input.bam -o clipped.bam` and perform the same repairs; preserve original reads. Pre-alignment FASTQ trimming is handled by read QC.

Done when the selected clipped BAM is repaired and verified.
