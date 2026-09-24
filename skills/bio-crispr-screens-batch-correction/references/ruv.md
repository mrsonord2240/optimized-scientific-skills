## RUV (Remove Unwanted Variation)

**Goal:** Identify hidden batch sources via control sgRNAs whose true signal is known.

**Approach:** Designate non-targeting controls as "negative controls" (assumed unchanged); RUV decomposes their variance into unwanted factors, then subtracts these from all data.

```r
library(RUVSeq)
# counts_df: rows = sgRNAs, columns = samples
# RUVg's SeqExpressionSet method takes cIdx as control ROWNAMES (character), not positions;
# which() returns integers and fails S4 dispatch here (the matrix method would accept them).
ntc_rownames <- rownames(counts_df)[rownames(counts_df) %in% ntc_sgrna_names]
stopifnot(length(ntc_rownames) > 0)
seqset <- newSeqExpressionSet(counts = as.matrix(counts_df))
ruv_corrected <- RUVg(seqset, cIdx = ntc_rownames, k = 2)  # k = 2 unwanted factors
# Two outputs. The W factors are what goes into a downstream model (MAGeCK MLE design matrix,
# edgeR/DESeq2 design); normCounts() is the adjusted matrix for PCA and visual checks.
W <- pData(ruv_corrected)          # W_1, W_2: the estimated unwanted factors, one column per k
corrected_counts <- normCounts(ruv_corrected)
```

**When to use:** RUV preferred over ComBat when batches are not annotated (e.g., unknown technical confounders). Worse than ComBat when batch is known and well-annotated; ComBat is more direct.
