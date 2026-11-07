# Sierra alternative polyadenylation

```bash
regtools junctions extract -a 8 -m 50 -M 500000 -s RF -o junctions.bed possorted_genome_bam.bam
```

Use the GTF, BAM, barcodes and junctions in `references/sierra.md`. This tests APA, not cassette exon splicing. Use donor-defined replicates where available.

Done when requested APA results are produced with their replication stated.
