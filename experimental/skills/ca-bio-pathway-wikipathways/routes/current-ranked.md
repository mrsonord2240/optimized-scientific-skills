# Current-release ORA and ranked GSEA

```bash
Rscript examples/wikipathways_ora.R
```

This source example requires `de_results.csv` in the current directory with `gene_symbol`, `padj`, and `log2FoldChange`; adapt input paths, organism and annotation package before use. It maps symbols to Entrez, uses tested genes as the ORA universe, and runs GSEA on a decreasing named Entrez ranking with seed 123. It contacts the current WikiPathways release, which can change over time, and writes example outputs to R's temporary directory. For reproducible dated ORA use [pinned ORA](pinned-ora.md). Dated GSEA requires adapting the pinned TERM2GENE/TERM2NAME construction to `GSEA`; no complete file-based dated GSEA interface is supplied. Verify identifier mapping, duplicates and ranking alignment. See `bio-pathway-gsea` for ranking mechanics and `bio-pathway-enrichment-visualization` for plotting.

**Done when:** results, mapping, ranking, universe where applicable, and release provenance are retained and interpreted with the caveats in the reference.
