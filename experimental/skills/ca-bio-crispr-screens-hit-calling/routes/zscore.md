# Custom z-score hit calling
~~~bash
python scripts/custom_zscore_hit_calling.py counts.tsv --ctrl T0_1,T0_2 --treat T18_1,T18_2 --ntc-prefix NonTargeting -o zscore_hits.tsv
~~~
Count TSV has guide ID first, Gene and sample columns. The source uses RPM, per-guide log2 fold change with pseudocount 1, gene mean, a control-derived median and standard deviation, two-sided normal tails and BH correction. Ensure non-targeting genes match the prefix and yield a nonzero finite null standard deviation; omitting the prefix uses all genes and assumes fewer than 40% change.
Done when input columns, null assumptions, outputs and guide-level evidence are reviewed.
