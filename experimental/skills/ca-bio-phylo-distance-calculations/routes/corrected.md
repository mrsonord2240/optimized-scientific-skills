# Model-corrected distance tree
~~~bash
Rscript scripts/model_corrected_tree.R alignment.fasta distance_tree.nwk 100
~~~
The source example estimates gamma shape from a JC starting tree with phangorn GTR, computes TN93+gamma distances, FastME balanced tree and 100 bootstrap replicates with seed 42. Confirm alignment, model assumptions, nonfinite distances and installed ape/phangorn versions. Bootstrap support describes sampling precision, not accuracy.
Done when model, gamma, tree and numerical limitations are reviewed.
