# PCA, batch detection and exclusion rules

Read when checking for batch effects, deciding whether to exclude a sample, or before asking the user about exclusions (the stop conditions are in the closing paragraph).

## PCA and Batch Detection

**Goal:** See whether the dominant variance is biology or batch, and flag outlier samples.

**Approach:** On the normalized survivors, run PCA, color by condition and by batch, and test whether top PCs associate with batch.

```python
import sys; sys.path.insert(0, "scripts")
from pca_batch import pca_batch_check  # scripts/pca_batch.py
coords, evr, tests = pca_batch_check(normalized_log2, sample_info, batch_col="batch")  # sample_info indexed by sample name
```

The third return value `tests` carries the per-PC status (`tested` / `not_testable`); a `not_testable` row is never evidence of no batch effect. A sample isolated from its group is a removal/re-run candidate, but a high-missing sample is judged by `raw_sample_qc`, not by PCA. If batch is PC1, keep batch in the design matrix for the differential test (preferred when batch and condition are balanced), and use `limma::removeBatchEffect` (or ComBat) only on the matrix used for PCA/plots to re-inspect biology; do not test on a batch-corrected matrix and also model batch. If batch is FULLY confounded with condition (every batch level holds exactly one condition) nothing can be corrected: batch and condition are the same variable, and removing one removes the other -- on a fully confounded synthetic set the mean |log2FC| of 104 truly-changed proteins went from 1.55 to 0.00 after batch removal. Report the design as non-identifiable and stop; do not correct, and do not test. Document and justify every exclusion, and re-run the downstream check with and without borderline samples (`references/qc-report-template.md` has the exclusion log and the with/without table). Stop and ask before excluding samples, when n < 5 per group makes PCA unstable, or when no un-normalized column is available for the loading check. Visualization of the projection routes to data-visualization/dimensionality-reduction-plots.
