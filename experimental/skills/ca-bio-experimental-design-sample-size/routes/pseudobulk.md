# Donor-level pseudobulk sizing
~~~sh
Rscript scripts/pseudobulk_donor_ssize.R cell_counts.rds donor_condition.csv 1.5 0.05 0.80 200 20260918 > donor_summary.txt
~~~
RDS input is a named list of genes-by-cells count matrices, one per donor; CSV has unique donor and condition columns. Source sums counts per donor, fits DESeq2 and calls ssizeRNA_vary. Gene order and identity must match across donors; enough independent donors and an appropriate model are required. The output is donors per group, not cells. Check assumptions and installed APIs; simulation and pilot uncertainty remain.
**Done when:** donor alignment, assumptions, seed, result and limitations are documented.
