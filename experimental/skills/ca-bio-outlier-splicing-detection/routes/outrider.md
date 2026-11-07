# Expression outliers

```bash
Rscript scripts/outrider.R counts.tsv PATIENT_001 patient_outliers.tsv
```

Input is a raw gene-by-sample count table including patient and matched controls. Use at least 50 samples; report low power and matching limitations. The script estimates q, floors it at 2, and uses padj 0.05 with zScoreCutoff 0. Check tissue and batch using `references/cohort-and-tissue.md`.

Done when the patient table and cohort limitations are reported.
