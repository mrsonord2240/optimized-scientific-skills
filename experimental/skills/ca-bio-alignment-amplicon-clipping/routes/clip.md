# Clip and repair primers

```bash
bash scripts/ampliconclip_workflow.sh input.bam primers.bed reference.fa output.bam
```

Input is a coordinate-sorted indexed BAM, indexed reference FASTA, and BED with 0-based half-open intervals and strand in column 6. Verify that scheme, contigs and reference build match. Default soft clipping uses `--both-ends --strand`; use `CLIP_OPTS=--strand` only when reads cannot reach the opposite primer. It leaves any 3-prime primer in place. `THREADS` defaults to 4. The workflow rejects more than 1% not clipped; raise `MAX_NOT_CLIPPED_PCT` only after review for a documented sparse control. It repairs mate fields, sorts, recalculates MD/NM, indexes, and checks residual primer footprints.

Done when the command succeeds and output BAM and index exist.
