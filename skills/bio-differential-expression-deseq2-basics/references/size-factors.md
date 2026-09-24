# Size factor alternatives and majority-DE normalization

## Size Factor Alternatives

`estimateSizeFactors()` defaults to `type='ratio'` (median-of-ratios). Edge cases:

| Situation | Use |
|-----------|-----|
| Zero counts in some samples for many genes | `type='poscounts'` -- uses only positive entries per gene |
| Very small libraries, hard to converge | `type='iterate'` |
| Spike-ins (ERCC) or known stable housekeeping genes | `controlGenes = indices` |
| Majority-DE biology (prokaryotic stress, viral host shutoff, MYC amplification) | `controlGenes` with curated stable genes; or spike-in SBN (Jiang 2011 *Genome Res* 21:1543) |

Single-cell pseudobulk: most genes have zeros across donors, so `type='poscounts'` is often required for pseudobulk DESeq2.

### Median-of-ratios fails on prokaryotic stress

**Trigger:** Bacterial RNA-seq under stress where >50% of genes change in one direction; PCA shows huge global shift.

**Mechanism:** Median-of-ratios assumes most genes are not DE. Under massive global perturbation, the reference is dominated by DE genes; size factors absorb the biology.

**Symptom:** MA plot shows the bulk cloud shifted off zero; reported fold changes don't match qPCR.

**Fix:** Use `controlGenes` with curated stable housekeeping genes; or spike-in normalization (Jiang 2011); or RUVg with negative controls.
