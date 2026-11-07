# Ranked and topology alternatives

```bash
Rscript examples/kegg_spia_topology.R
```

This source example expects `de_results.csv` with `gene,log2FoldChange,pvalue,padj`, uses human-specific mapping and fixed bootstrap settings, and writes results in a temporary directory. It can require live/reference resources and substantial compute; adapt organism, inputs, universe and output before use. SPIA and graphite can use different topologies and disagree in direction; inspect both and report the actual topology. The source also describes `gseKEGG`, module enrichment, multi-condition comparison and pathview; these require input-specific adaptation and current API checks. Pathview downloads reference images. Preranked GSEA and ORA test different questions, and no topology score alone proves activation.

Done when the method, gene universe, reference, mapping losses, assumptions and uncertainty are documented.
