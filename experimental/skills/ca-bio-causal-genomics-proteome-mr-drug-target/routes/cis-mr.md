# Cis-MR and PAV sensitivity
~~~sh
Rscript scripts/cis_panel.R harmonized.tsv vep_annotations.tsv result.rds
~~~
Provide already selected, ancestry-matched, clumped and TwoSampleMR-harmonized cis instruments, and VEP annotations with SNP and Consequence columns. The source selects variants within +/-500 kb of the gene, P < 5e-8 and exposure F >= 10, clumping at r2 < 0.1 within 500 kb and harmonizing with action=2. Confirm those steps and the genome build first. The wrapper selects Wald for one instrument, IVW and weighted median for 2-9, and adds Egger for at least 10, then reruns excluding the listed PAV consequences. Do not treat an empty PAV-excluded panel as confirmation. See [failure modes](../references/failure-modes.md) and [source thresholds](../references/decision-framework.md).
**Done when:** both estimates and retained counts have been inspected, and colocalization, platform replication, neighbouring genes and overlap have been assessed.
