# Ruv
~~~bash
Rscript scripts/ruv.R counts.tsv ntcs.txt factors.tsv adjusted.tsv 2
~~~
Counts are feature-indexed; list one control per line. This uses the source example k=2; evaluate k. Include unwanted-variation factors in downstream models; assess adjusted counts for visualization. See [RUV notes](../references/ruv.md).
Done when factor table and adjusted matrix are saved and k and controls reviewed.
