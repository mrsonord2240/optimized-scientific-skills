# Cohort-complete OncoPrint

```bash
Rscript examples/oncoprint_phd.R cohort.maf clinical.tsv oncoprint.pdf
```

The clinical TSV needs `Tumor_Sample_Barcode`, `Subtype`, and `Stage`; add colors for every observed level. The example makes a top-20 cohort-complete matrix, aligns clinical records by exact ID and orders genes by number of altered samples. It also runs `somaticInteractions` and writes a second PDF, even for small cohorts; review [interaction cautions](interactions.md) and adapt the script if pairwise testing is inappropriate. Its TMB fallback is a raw MAF-row count when `tmb` is absent, not a normalized rate; label it accordingly or supply a valid measure. Output names may be overwritten. Additional CNV/fusion calls are supported by `scripts/maf_to_oncoprint.R` but the example does not pass them; adapt the input pipeline. See [data rules](../references/data.md).

Done when matrix classes, cohort denominator, annotation alignment, ordering, output and statistical limits are checked.
