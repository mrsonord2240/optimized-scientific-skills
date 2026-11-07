# SIMEX correction for MR-Egger
~~~sh
Rscript scripts/simex_egger.R harmonized.tsv simex.rds
~~~
The source workflow fits weighted Egger regression and applies quadratic SIMEX with lambda 0.5 to 2 in 0.5 steps and 1000 simulations. The source applies SIMEX only when I-squared is at least 0.6 and below 0.9; the wrapper stops outside that interval. Inspect I-squared and extrapolation; correction can be unstable and does not resolve correlated pleiotropy. Consult [failure modes](../references/failure-modes.md).
**Done when:** the RDS contains the Egger fit and SIMEX result, with extrapolation checked.
