# Merge BAMs

```bash
bash scripts/merge_bams.sh merged.bam sample1.bam sample2.bam sample3.bam
ORDER=name-ascii bash scripts/merge_bams.sh merged_N.bam a_N.bam b_N.bam
```

`ORDER` is `coordinate` (default), `name-natural` or `name-ascii` and must match how the inputs were sorted; a plain merge does not give name order. The script stops unless all inputs print one identical `@HD` line (`samtools merge` itself exits 0 and writes out-of-order records for mixed inputs; `-n` and `-N` inputs differ in `SS:`), then runs `samtools merge -c -p`, `quickcheck`, and indexes a coordinate merge. `THREADS` defaults to 4.

- `-c` deduplicates `@RG` by ID, so two different read groups (different SM/PU) sharing an ID collapse and their reads are mislabelled; GATK BQSR then recalibrates wrongly. Make colliding IDs unique with `samtools addreplacerg` before merging.
- Region only (every input indexed, else `Could not retrieve index file`): `samtools merge -R chr1:1000000-2000000 out.bam a.bam b.bam`.
- Many files: `samtools merge -b files.txt merged.bam` (one path per line). Overwrite: `-f`.

Done when the merged BAM passed `quickcheck` and, for coordinate order, has an index.
