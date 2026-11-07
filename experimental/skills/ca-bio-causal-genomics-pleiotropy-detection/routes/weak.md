# Weak-instrument and cluster sensitivity
~~~sh
Rscript scripts/raps_clust.R harmonized.tsv raps_clust.rds
~~~
Uses the source RAPS settings: overdispersion on, Huber loss, shrinkage off. MR-Clust estimates clusters from the same retained harmonized variants. Check installed package APIs and clustering fit. Review [method notes](../references/mr-raps-mr-clust.md) and instrument strength before interpretation.
**Done when:** both results are saved and checked for convergence and assumptions.
