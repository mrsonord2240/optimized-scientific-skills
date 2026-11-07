# Salmon isoform switches

```bash
Rscript scripts/isoform_switch_analysis.R salmon_quant annotation.gtf transcripts.fa sample_metadata.tsv isoswitch_out
```

Metadata needs sample_id and condition; optional columns are covariates joined by sample ID. Exactly two conditions and at least two replicates each are required; two per group is exploratory. The script uses raw NumReads, DEXSeq when all groups have at most five replicates, otherwise satuRn, q < 0.05 and |dIF| > 0.1. Run label permutations and a manual gene-level confirmation before trusting results, as in `references/switch-workflow.md`; compare count routes if necessary. Optional sixth argument supplies annotator files. For visualization, reconciliation and error guidance see `references/interpretation.md`.

Done when switch tables, consequences available from supplied annotations, and null checks are interpreted.
