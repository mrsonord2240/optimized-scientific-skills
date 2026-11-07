# Annotation-free junction outliers

```bash
regtools junctions extract -a 8 -m 50 -s XS sample.bam -o sample.junc
```

Run once per indexed BAM with XS strand tags, then cluster junctions and run LeafcutterMD using `references/leafcutter.md`. Its p-values are raw; apply BH over cluster-by-sample tests as shown there. It is junction-only; intron retention and pseudoexon clusters may lack p-values.

Done when the corrected outlier table is produced and limitations reported.
