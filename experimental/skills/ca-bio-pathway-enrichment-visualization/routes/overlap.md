# Show term overlap

```bash
Rscript scripts/plot_result.R result.rds emap figure.pdf 30
```

The wrapper computes `pairwise_termsim` with the source default Jaccard method before the enrichment map. State similarity, edge settings, selected terms and total result count. GO semantic similarity, treeplot, cnetplot, heatplot, upset and compareCluster layouts require input-specific adaptation; inspect installed signatures because enrichplot APIs vary by release. The retained examples use synthetic inputs and may run enrichment before plotting; they are not input-file wrappers. Similarity clustering illustrates overlap, not independent biological replication.

Done when term selection, overlap method, truncation and visual interpretation are recorded.
