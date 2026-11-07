# Colocalization and target nomination
The source coloc example assumes fixed study-specific sample sizes, outcome type and case fraction. It is not a general turnkey analysis: prepare full locus-wide matched-allele summary statistics and validate study-specific inputs, then use the colocalization skill. For multiple signals use ancestry-matched LD and coloc.susie per credible set. Use cis-MR, colocalization, cross-platform replication, PAV exclusion, neighbour-gene checks and sample-overlap review together; source thresholds and caveats are in [decision framework](../references/decision-framework.md), [failure modes](../references/failure-modes.md) and [target nomination](../references/target-nomination.md).
PAV annotation requires a valid study VCF and assembly-matched VEP installation:
~~~sh
vep --species homo_sapiens --assembly GRCh38 --canonical --check_existing --input_file cis_pqtls.vcf --output_file pqtl_vep.tsv --tab
~~~
Adapt assembly and output parsing for the study; do not use fabricated variant records. **Done when:** evidence and limitations for each stream are documented. These are population research results, not individual treatment recommendations.
Related shelf skills: bio-causal-genomics-colocalization-analysis, bio-causal-genomics-fine-mapping, bio-causal-genomics-effector-gene-prioritization, bio-causal-genomics-transcriptome-wide-association.
