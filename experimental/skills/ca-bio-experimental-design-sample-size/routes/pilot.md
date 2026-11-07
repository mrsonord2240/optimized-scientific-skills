# Pilot-based PROPER simulation
~~~sh
Rscript scripts/proper_power.R pilot_counts.csv 3,6,10,20 20 1.5 0 20260918 > proper_summary.txt
~~~
Input CSV is genes by samples with gene IDs in column one. Optional positional inputs are per-group replicate grid, simulation count, fold change, maximum genes (0 means all) and seed. Source uses p.DE=0.05, DESeq2 and nominal FDR 0.05; check the installed PROPER API and the source's scale conventions before use. Output includes marginal power and actual FDR; do not select a design based on power alone. Monte Carlo uncertainty and pilot uncertainty matter.
**Done when:** summary, simulation settings, seed and realized error-rate assessment are retained and reviewed.
