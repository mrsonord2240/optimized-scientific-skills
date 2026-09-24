# DESeq2 with tximport (Salmon / kallisto / RSEM)

## Tximport (Salmon / kallisto / RSEM)

For salmon/kallisto/RSEM input, use `DESeqDataSetFromTximport()` which carries the per-sample length matrix as an offset automatically:

```r
library(tximport)
txi <- tximport(files, type = 'salmon', tx2gene = tx2gene)
dds <- DESeqDataSetFromTximport(txi, colData = samples, design = ~ condition)
dds <- DESeq(dds)
```

The `tximport(..., countsFromAbundance='lengthScaledTPM')` form is for limma-voom (no offset mechanism). Full mechanics, the four `countsFromAbundance` options, RSEM zero-length traps, and tximeta provenance: see `expression-matrix/counts-ingest`.
