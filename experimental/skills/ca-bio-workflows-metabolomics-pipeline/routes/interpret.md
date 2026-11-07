# Interpret the QC handoff
~~~bash
Rscript scripts/qc_cli.R imported.rds qc.rds
~~~
Done when each feature has its acquisition mode and annotation confidence, statistics have BH FDR plus validated multivariate evidence, and pathway background is declared. This command creates the QC handoff only; annotation, guarded OPLS-DA and pathway mapping require component workflows and adaptation, not a complete pipeline CLI. Use only MSI levels 1-2 for identified ORA with assay-coverage background, or the full feature table for mummichog. See [pathway recipe](../references/stage5-pathway-mapping.md); remote services may be involved.
