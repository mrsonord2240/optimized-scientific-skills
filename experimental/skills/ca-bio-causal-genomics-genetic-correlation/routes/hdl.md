# HDL

```bash
Rscript scripts/hdl_rg.R trait1.tsv trait2.tsv UKB_array_SVD_eigen90_extraction hdl_rg.txt 0
```

Use 0 only for known independent cohorts. Supply correctly formatted tables and matched eigen reference as in `references/hdl.md`. The source contains conflicting statements about N0 correction and overlap robustness; do not infer independence from an intercept alone. If overlap is unknown or substantial prefer LDSC while checking the installed HDL interface.

Done when rg, SE, overlap and reference limitations are reported.
