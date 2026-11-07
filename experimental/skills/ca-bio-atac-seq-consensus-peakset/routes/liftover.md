# Peaks from another genome build

```bash
liftOver published_peaks.hg19.bed hg19ToHg38.over.chain.gz published_peaks.hg38.bed published_peaks.unmapped.bed
cat published_peaks.hg38.bed my_peaks.hg38.bed | sort -k1,1 -k2,2n | bedtools merge -i - > merged_consensus.bed
```

Shown, not scripted. Convert to the build of your own peaks first, then build the consensus on one build.

- Report the chain file version and the unmapped fraction. Roughly 3-5% of regions fail (complex regions, sex chromosomes, MHC).
- For high-stakes cross-build comparison, re-align instead of lifting over.
- Human to mouse: use synteny-based mapping (`bnMapper` or `halLiftover`), not chain liftOver.
- Lifted peaks have no summit; for fixed width extend from the midpoint and say so.

Done when the answer reports the unmapped fraction and the chain used.
