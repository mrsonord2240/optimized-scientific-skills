# Replicate correlation and fingerprints

```bash
multiBamSummary bins -bs 10000 -p 8 --bamfiles rep1.bam rep2.bam rep3.bam -o multi.npz
plotCorrelation -in multi.npz --corMethod spearman --whatToPlot heatmap --skipZeros -o spearman_heatmap.png
plotFingerprint -p 8 -b rep1.bam rep2.bam rep3.bam --labels rep1 rep2 rep3 --skipZeros --numberOfSamples 50000 -o fingerprint.png --outQualityMetrics fingerprint_metrics.txt
```

Shown, not scripted.

- Spearman, not Pearson (ATAC signal is heavy-tailed); `--skipZeros` is mandatory or empty regions dominate.
- Below 0.85 Spearman is non-trivial divergence; 0.9 or above is a common working target (not an ENCODE ATAC threshold). Check PCA: if replicates separate by batch, add batch as a covariate downstream (`~Batch + Condition`); if an outlier, drop it or rerun. Do not silently merge replicates that correlate badly.
- Poor correlation but sensible PCA: background dominates; correlate on peak counts only.
- The plain Jensen-Shannon distance column in fingerprint metrics exists only with a reference sample (`--JSDsample`); without one, the synthetic JS distance is reported. Larger means stronger enrichment.

Done when each replicate has a correlation and a keep/drop decision.
