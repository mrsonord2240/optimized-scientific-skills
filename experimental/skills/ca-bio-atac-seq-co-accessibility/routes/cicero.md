# Cicero peak-pair links

```bash
Rscript scripts/cicero_workflow.R peaks.mtx peak_metadata.tsv cell_metadata.tsv tss.bed
```

Omit the fourth argument to skip gene mapping. Matrix rows are peaks and columns cells; metadata row names must be in matching order. Peak names must be `chr_start_end`; TSS BED column 4 contains gene names. The script computes its own LSI and UMAP, uses k=50, 500 kb window, seed 1 and score >0.25, and writes `cicero_connections.tsv` plus optional `cicero_enhancer_gene_pairs.csv` in the current directory. Candidate pairs flagged `enhancer_is_promoter` are promoter-promoter. It removes NA and duplicate orientations. Assess coordinate conventions before mapping and validate the promoter-window conversion; see `references/cicero.md`.

Run separately by cluster when cell types are heterogeneous; adjust the function arguments in an adapted workflow if a different window or threshold is required. The CLI does not accept an output prefix or those options.

Done when the pair table and optional gene candidates are written and filtering, window, cluster and coordinate conventions are reported.
