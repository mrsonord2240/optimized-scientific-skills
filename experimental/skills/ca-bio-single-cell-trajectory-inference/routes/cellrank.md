# CellRank fate mapping
~~~sh
python scripts/cellrank_fate.py trajectory.h5ad --time-key dpt_pseudotime --cluster-key leiden --n-states 10 --w-pseudotime 0.8 --out fate.h5ad
~~~
Requires neighbors and a justified pseudotime and cluster column. Source combines pseudotime and connectivity kernels, runs GPCCA and prints state/entropy summaries. Select state count and kernel weights from diagnostics and sensitivity analysis; probabilities inherit the input direction and are not observed fate. Inspect actual saved object to confirm expected probabilities were persisted, and verify the installed API.
**Done when:** terminal states, probability output, kernel assumptions and root sensitivity are reviewed.
