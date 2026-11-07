# Salmon, kallisto or RSEM output

```r
library(tximport); library(DESeq2)
source('scripts/deseq2_de.R')
txi <- tximport(files, type = 'salmon', tx2gene = tx2gene)   # files: named vector of quant.sf paths
dds <- DESeqDataSetFromTximport(txi, colData = samples, design = ~ 1)
out <- run_deseq2(dds, condition = 'condition', ref = 'control')
write.csv(out$results, 'results.csv', row.names = FALSE)
```

- Do not round the estimated counts and pass them as a matrix. `DESeqDataSetFromTximport()` carries the transcript-length correction; a rounded matrix loses it.
- `type` is `'salmon'`, `'kallisto'` or `'rsem'`. `tx2gene` is a two-column table of transcript ID and gene ID.
- `run_deseq2()` takes the same options as the command line: `block`, `versus`, `test`, `shrink`.

Done when `results.csv` exists and you have answered from it.

More detail: `references/tximport.md`.
