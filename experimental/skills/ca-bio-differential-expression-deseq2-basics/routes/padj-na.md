# A gene shows padj = NA

The script prints three counts. Find which one the gene falls under in `results.csv`:

| In the results file | Cause | What to do |
|---------------------|-------|------------|
| `baseMean` is 0 | No reads in any sample | Nothing to test. Report it as not detected. |
| `pvalue` is NA, `baseMean` above 0 | One sample is an extreme outlier for this gene | Look at that gene's counts. If the outlier sample is the biology, rerun with `all_genes_bh=true`. |
| `pvalue` present, `padj` NA | Count too low to be worth testing; set aside to keep power | Report it as below the filter. To adjust over every gene, rerun with `all_genes_bh=true`. |

- `padj = NA` is not "not significant". Report these genes separately.
- A gene with zero counts in one group and real counts in the other is tested normally and gets a very large fold change. Check it by eye.

Done when the gene's status is stated with its cause.

More detail: `references/padj-na-remedies.md`.
