# FRASER splicing outliers

```bash
Rscript scripts/fraser2_rare_disease.R bams PATIENT_001 fraser_workdir
```

Input is indexed coordinate-sorted paired-end BAMs named by sample ID, including the patient and tissue- and batch-matched controls. Inspect sample correlations and PCA before fitting as in `references/cohort-and-tissue.md`; an isolated patient requires investigation. Use at least 20 samples; 50 or more recommended. The script estimates q, fits PCA using Jaccard, filters padj < 0.05 and absolute delta >= 0.1, and writes patient results plus a volcano plot. Check q and tuning in `references/latent-dimension.md`. An empty result from too few or mismatched samples is inconclusive.

Done when the patient output and fit diagnostics are reviewed and limitations reported.
