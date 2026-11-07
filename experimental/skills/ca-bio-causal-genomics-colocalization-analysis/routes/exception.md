# Exceptional regions and reporting

```bash
Rscript scripts/check_excluded_region.R 6 30000000 region_check.tsv
```

Supply the actual chromosome and hg38 position. The result checks the lead position against the source's named exclusion regions; it does not establish that the whole locus is outside them. MHC chr6:25-35 Mb requires an HLA-aware method or exclusion; chr8:8.1-11.9 Mb requires inversion conditioning or exclusion. Check the entire analysis window and genome build. For window sensitivity read `references/lead-snp-swap.md`; for low power and LD mismatch read `references/failure-modes.md`; for reporting read `references/reporting.md`.

Done when the exclusion check and whole-locus review are complete and any exception is resolved or explicitly excluded.
