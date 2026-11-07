# Curated phenome scan
~~~sh
Rscript scripts/phewas_curated_endpoints.R pqtl.tsv endpoints.tsv results.tsv 50000 European
~~~
Use a preselected, clumped cis-pQTL table with SNP, BETA, SE, A1, A2, EAF and P, and an endpoint table with OpenGWAS study IDs in id. Requires OpenGWAS access; the source script filters the catalog by sample size and population and uses the full curated-list length as the Bonferroni denominator. It uses IVW and requires at least two outcome variants. Document skipped, missing and failed outcomes; lack of a result is not evidence of no effect. See [phenome notes](../references/phewas.md).
**Done when:** the TSV contains results for successfully queried outcomes and the denominator, exclusions and adverse-effect signals are reviewed.
