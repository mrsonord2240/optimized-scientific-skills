# BaalChIP

```bash
Rscript scripts/baalchip_workflow.R --samples samples.tsv --hets hets.tsv --group TUMOR --blacklist blacklist.bed --imprinted imprinted.bed --sex female --assembly GRCh38 --samples-assembly GRCh38 --hets-assembly GRCh38 --blacklist-assembly GRCh38 --imprinted-assembly GRCh38 --correction raf --samtools /path/to/samtools --out results/TUMOR
```

Run from the skill directory. Supply indexed bias-filtered BAMs and peak paths in the five-column sample sheet, a sample-specific heterozygous table with measured per-variant RAF, and assembly-matched blacklist and imprinting BEDs. Alternatively use `--correction gdna --gdna-bam matched.bam`; population AF or detached copy-number segments are not RAF. All declared assemblies and contig styles must agree. The output directory must not exist and its parent must exist. `--preflight-only` validates and writes a preflight output directory without model execution; choose a separate output for the model run. Source runner binds BaalChIP 1.36.0; verify installed API. It uses minimum base quality 10, MAPQ 15, 5,000 iterations and 95% confidence. See [input contracts](../references/methods-and-commands.md).

**Done when:** the new output directory contains status, correction provenance, exclusions and, for a model run, full and significant tables or an explicit empty state. Report the assay, corrections, depth and stage attrition.
