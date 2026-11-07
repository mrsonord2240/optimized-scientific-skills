# Covariate
~~~bash
mageck mle --count-table counts.txt --design-matrix design.txt --permutation-round 10 --output-prefix batch_aware_mle
~~~
Prepare a sample-by-covariate design matrix with baseline, batch and treatment columns as in the source example. Check sample names and full rank. Fully confounded treatment and batch are not identifiable. See [design and validation](../references/source-details.md).
Done when MAGeCK results are written and the design and adjustment reviewed.
