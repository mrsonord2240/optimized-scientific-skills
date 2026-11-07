# Junction clusters

```bash
regtools junctions extract -a 8 -m 50 -s XS sample.bam -o sample.junc
```

BAMs need XS tags; check that junction strands are not all unknown. Extract each sample, list files, then cluster with the leafcutter repository script as in `references/leafcutter.md`. Read `references/failure-modes-and-errors.md` for topology and reconciliation.

Done when the requested intron clusters and usage counts are available.
