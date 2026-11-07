# Compare discrete-trait models and sample histories
```bash
Rscript scripts/stochastic_mapping.R --seed=20260928 > discrete-results.txt
```
**Done when:** output metadata is recorded and node probabilities and mapping histories are retained after successful completion. The source script reads `species_tree.nwk` and `traits.csv` from the working directory; the CSV has taxon row names and a `state` column. It requires exact taxon matching, at least two states and an explicit integer seed. Source settings compare ER/SYM/ARD, use 1000 maps and a Fitzjohn root prior, and fit a two-category hidden-rate alternative; assess whether the data support that model. Inspect node and hidden-rate outputs in the script and adapt output persistence before claiming a complete analysis. See [method workflows](../references/method-workflows.md).
