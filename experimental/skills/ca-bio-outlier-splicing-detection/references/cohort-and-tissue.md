# Cohort Size and Power


Measured 2026-09-20 on synthetic data (FRASER: planted events in a 30-sample cohort, 2x75 reads; OUTRIDER: `makeExampleOutriderDataSet(n = 2000 genes, freq = 1e-3)`, 3 seeds, padj<0.05).

**FRASER 2** (single patient + n-1 controls, PCA, q from `estimateBestQ`, FRASER 2.6.1): n=8, 12, 16: planted exon-skipping missed (0 calls); n=20, 25: flagged (3 junctions); n=30: all 4 planted events (skipping, cryptic donor, intron retention, pseudoexon) flagged on both FRASER 2.2.0 and 2.6.1; 0 false calls at every n. 4 real chrX samples: runs, 0 calls (n too small).

**OUTRIDER** (fraction of injected outliers detected; false calls 0-2 per 3 seeds):

| Samples | OUTRIDER 1.24.0 | OUTRIDER 1.28.1 |
|---------|-----------------|-----------------|
| 20 | 0/107 | 0/107 |
| 30 | 0/156 | 5/156 (3%) |
| 50 | 111/297 (37%) | 123/297 (41%) |
| 60 | 171/359 (48%) | 201/359 (56%) |
| 100 | 414/583 (71%) | 434/583 (74%) |
| 200 | 857/1203 (71%) | 921/1203 (77%) |

The q that `estimateBestQ` returns matters: on a 3000-gene synthetic matrix with a hidden batch (8 planted outliers) the OUTRIDER block above found 0/8 at n=30 on both versions and, at n=100, 4/8 with the 1.24.0 heuristic q=20 versus 2/8 with the 1.28.1 optimal-hard-threshold q (1, raised to 2; 1 false call). Fixed small q does not rescue n=30 (q=2 or 4: 0/156 on both versions) and q=2 at n=50 gave up to 39 false calls. Practical rule: OUTRIDER **>=50 samples** (useful power from ~100); FRASER **>=20** (more is better; the 4-event result above is n=30). Autoencoder covariance is learned across samples, so batch- and tissue-matched controls matter more than count.

| Cohort size | FRASER 2 | OUTRIDER |
|-------------|----------|----------|
| n < 20 | Missed at n=8-16; do not interpret a null | No power |
| n = 20-49 | Detected at n=20-30 in testing | Effectively none (0-3% at n=20-30) |
| n >= 50 | Recommended | ~40-55% at 50-60 |
| n >= 100 | Best calibration | ~70-75% |

GTEx-derived tissue-matched controls can supplement small in-house cohorts but introduce batch effects (and need dbGaP access); use only when in-house n < 30, keep them a minority, and document the pooling strategy. The tools do not correct for batch.

## Tissue Choice for Mendelian RNA-seq

| Tissue | Pros | Cons |
|--------|------|------|
| Whole blood (PAXgene) | Easy, standard | Globin contamination; many disease genes silent |
| Fibroblast (skin biopsy) | Reasonable expression | Requires culture; senescence variability |
| Muscle biopsy | Best for muscular dystrophy | Invasive |
| iPSC-derived neuron / cardiomyocyte | Disease-relevant tissue | Cost, variability |
| Urine sediment | Non-invasive | Low yield |

For UDN-style cases: blood first, then fibroblast if blood lacks expression of the candidate gene. **A negative blood RNA-seq does not rule out a candidate gene that is silent in blood** - verify gene expression with the GTEx tissue panel before committing to the tissue.

**Detecting a tissue/batch mismatch before interpreting.** Strict tissue matching is required, and a mismatched patient does not always look like "hundreds of outliers": on a synthetic sample with globally shifted splicing, FRASER PCA reported 0 calls (the fit absorbed it), FRASER AE q=10 reported 15, LeafcutterMD 9. So do not read a zero-call patient as tissue-matched. Before fitting FRASER, inspect sample-level count correlations and PCA using a matched gene-count matrix (rows = genes, columns = the same samples); a patient isolated from the control cluster is a stop-and-investigate signal:

```r
# Replace gene_counts.tsv with your pre-FRASER, unnormalised gene x sample matrix.
counts <- as.matrix(read.delim('gene_counts.tsv', row.names = 1, check.names = FALSE))
log_cpm <- log2(sweep(counts + 0.5, 2, colSums(counts + 0.5), '/') * 1e6 + 1)
sample_correlation <- cor(log_cpm, method = 'spearman')
sort(sample_correlation['PATIENT_001', setdiff(colnames(sample_correlation), 'PATIENT_001')])

sample_pca <- prcomp(t(log_cpm), scale. = TRUE)
plot(sample_pca$x[, 1], sample_pca$x[, 2], xlab = 'PC1', ylab = 'PC2')
text(sample_pca$x[, 1], sample_pca$x[, 2], labels = rownames(sample_pca$x), pos = 3)

# Calls are a secondary check only; PCA may absorb a global mismatch.
sort(table(factor(all_results$sampleID, levels = colnames(fds))), decreasing = TRUE)[1:5]
```
