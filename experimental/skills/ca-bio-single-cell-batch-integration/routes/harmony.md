# Harmony in PCA space
~~~sh
python scripts/harmony_integration.py merged.h5ad integrated.h5ad --batch-key batch --seed 0
~~~
Input is a merged, QC-reviewed object with raw counts in X and batch metadata. The source script normalizes, log-transforms, selects 2000 variable genes, scales, computes PCA, Harmony, neighbors, UMAP and Leiden. Keep original counts and a separate uncorrected object; output contains an uncorrected PCA and corrected embedding, with genes subset to HVGs. Verify installed APIs and inspect [Harmony notes](../references/harmony.md).
**Done when:** the saved object, uncorrected comparison, batch mixing and biological preservation have been reviewed.
