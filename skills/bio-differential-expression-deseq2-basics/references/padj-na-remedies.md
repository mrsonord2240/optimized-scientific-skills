# padj = NA remedies: Cook's outliers and independent filtering

### Cook's silences the gene of interest

**Trigger:** Rare-disease cohort or CNV-amplified patient where ONE sample drives the biology; that gene has `padj=NA`.

**Mechanism:** Cook's distance flagged the patient as an outlier and zeroed the gene's p-value (or replaced its counts if n>=7).

**Symptom:** A gene of clear biological interest comes back as NA in the results table even though the count matrix shows the expected pattern.

**Fix:** `results(dds, cooksCutoff = FALSE)`. Optionally cross-validate the result by running a sensitivity analysis with and without the outlier sample.

### Independent filtering kills the master regulator

**Trigger:** A transcription factor expressed at ~10 counts but consistently across all samples shows `padj=NA` despite obvious biology.

**Mechanism:** baseMean is below the data-driven independent filtering threshold; the gene was excluded from FDR adjustment.

**Symptom:** Low-count genes with clean signal end up NA.

**Fix:** `results(dds, independentFiltering = FALSE)`; or `filterFun = ihw` from the IHW package (often less aggressive on low-count genes).
