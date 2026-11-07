# Independent groups

```bash
Rscript scripts/annotate_pairwise.R data.csv annotated.png wilcox holm
```

CSV requires nonmissing `group,value`; the script compares every pair, adjusts the full supplied family, and saves a figure and `.results.csv` with raw and adjusted p-values, method, family size and signed effect sizes. Modes are `wilcox`, `dunn`, and `tukey`; use `tukey tukey` for Tukey's own simultaneous adjustment. Optional fifth argument sets positive bracket spacing. Python alternative: `python scripts/annotate_pairwise.py data.csv annotated.png holm` uses Mann–Whitney and supplied adjusted values. Use textual group labels for the Python plotting route. Output files may be overwritten. Neither route defines a different scientific family automatically; adapt if the prespecified family differs.

Done when assumptions, family, exact results, effect direction and rendered labels are checked, including nonsignificant comparisons.
