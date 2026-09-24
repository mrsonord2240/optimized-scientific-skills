# Comparing conditions with compareCluster

Read when comparing KEGG enrichment across several gene lists (up vs down, contrasts, clusters).

## Compare Multiple Conditions

**Goal:** See shared and condition-specific KEGG pathways across groups in one faceted figure.

**Approach:** Pass named gene lists to compareCluster with fun='enrichKEGG'; it fits one model and produces a faceted dotplot. Compare pathway-ID SETS across conditions, never raw p-values (they depend on sample size, DE gene count, and the KEGG release).

```r
clusters <- list(up=up_entrez, down=down_entrez)
ck <- compareCluster(geneClusters=clusters, fun='enrichKEGG', organism='hsa', keyType='ncbi-geneid')
ck <- setReadable(ck, OrgDb=org.Hs.eg.db, keyType='ENTREZID')
# dotplot(ck) -> enrichment-visualization for the plot grammar
```
