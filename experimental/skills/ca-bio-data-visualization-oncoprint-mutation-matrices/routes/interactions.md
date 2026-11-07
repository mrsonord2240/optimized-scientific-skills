# Pairwise co-occurrence and exclusion

```bash
Rscript examples/oncoprint_phd.R cohort.maf clinical.tsv oncoprint.pdf
```

This full example also performs maftools pairwise interactions and creates an interactions PDF; see [input contract](r-oncoprint.md). It does not save the returned table, so adapt it to export raw 2-by-2 counts, p-values, adjusted p-values, odds ratios and event labels. maftools includes only samples represented in the MAF, not zero-row cohort members. Fisher tests ignore sample-specific mutation burden; the source suggests DISCOVER for heterogeneous backgrounds. Below 50 samples treat results as descriptive and report exact frequency intervals rather than discovery claims; 50–99 emphasize effect size and 100+ is a suggested pair-testing starting point. Do not infer direction solely from a corrected odds ratio in sparse tables.

Done when cohort denominator, raw counts, multiple testing and limitations are reported with any pairwise result.
