# Validate co-accessibility

```bash
Rscript scripts/cicero_workflow.R peaks.mtx peak_metadata.tsv cell_metadata.tsv
```

For matching Hi-C or Micro-C loop overlap and a locus-level arc plot, use `references/validation.md`; ensure compatible genome build and coordinate conventions. Overlap percentages are dataset dependent, and the source's numerical concordance bands are unsourced heuristics. Compare candidates with orthogonal experimental evidence where available.

Done when the comparison or visualization is produced with genome build, overlap rule, threshold and limitations reported.
