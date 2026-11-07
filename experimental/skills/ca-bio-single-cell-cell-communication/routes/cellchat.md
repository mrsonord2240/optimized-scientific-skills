# CellChat pathway summary
~~~sh
Rscript scripts/cellchat_run.R annotated.rds cell_type cellchat.rds interactions.csv
~~~
The source workflow uses human secreted-signaling database, trimean and minimum 10 cells, then aggregates pathways and centrality. Supply a compatible annotated Seurat object with normalized expression; use a species-matched resource and inspect assumptions. The output is model-derived evidence, not measured signaling. See [CellChat notes](../references/cellchat.md).
**Done when:** model, interactions, resource and population sizes are reviewed and retained.
