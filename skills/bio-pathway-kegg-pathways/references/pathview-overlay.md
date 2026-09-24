# pathview: overlay per-gene data on a KEGG map

Read when the user wants fold-changes drawn on the KEGG map image.

## Overlay Data on the KEGG Map (pathview)

pathview downloads a KEGG pathway's KGML and image, joins per-gene values to the nodes, and writes a colored map PNG/PDF (a KEGG-specific operation owned here; generic dot/cnet/emap plots route to enrichment-visualization). It writes files to the working directory and queries KEGG live.

```r
library(pathview)
vals <- setNames(de$log2FoldChange, de$entrez)
pathview(gene.data=vals, pathway.id='hsa04110', species='hsa', gene.idtype='entrez')   # writes hsa04110.pathview.png
```
