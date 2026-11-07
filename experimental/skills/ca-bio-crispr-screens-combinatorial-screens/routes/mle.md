# Fit MAGeCK MLE interaction
~~~bash
mageck mle --count-table combo_counts.txt --design-matrix combo_design.txt --output-prefix combo_mle
~~~
Prepare a sample-matched four-condition design with baseline, geneA, geneB and interaction columns; interaction equals one only in the double condition. See [source design and limitations](../references/mle.md). MAGeCK permutation-based FDR can vary on reruns; an interaction coefficient alone is not validation of an additive null. This source's claims about formal significance and deterministic z-score significance should be independently evaluated.
Done when output and design assumptions are reviewed against explicit GI scoring.
