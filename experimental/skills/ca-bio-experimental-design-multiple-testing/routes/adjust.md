# Adjust a p-value family
~~~sh
Rscript scripts/adjust_pvalues.R tests.csv pvalue BH adjusted.csv
~~~
Choose BH, BY, holm or bonferroni explicitly; CSV needs the named p-value column. Missing values remain missing; values outside 0-1 stop the script. BH has its standard guarantee under independence or suitable positive dependence (PRDS); arbitrary dependence is not guaranteed by a few simulations. BY controls FDR under arbitrary dependence at a power cost. Holm and Bonferroni control family-wise error. Record the family and exclusions.
**Done when:** adjusted values and the selected method are recorded and assumptions reviewed.
