# Storey q-value summary
~~~sh
Rscript scripts/qvalue_summary.R tests.csv pvalue qvalue_summary.rds 0.05
~~~
This runs the source qvalue_safe helper in an isolated process and saves its pi0, discovery count, method and worker status. It does not save per-feature q-values or local FDR; obtain these through a validated qvalue workflow if required. Small families can make pi0 estimation unstable. A failed worker may fall back to BH, and the saved method must be reported. Consult [methods and limitations](../references/methods.md).
**Done when:** summary method, pi0 availability and worker status are reviewed.
