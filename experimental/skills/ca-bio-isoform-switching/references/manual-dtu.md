## Manual DTU Pipeline (DRIMSeq + DEXSeq + stageR)

The canonical reference is the *F1000Research* "Swimming downstream" workflow (Love, Soneson, Patro 2018; Bioconductor `rnaseqDTU`). Below: tximport with raw counts (`countsFromAbundance = 'no'`, the route that worked on real data; see Count route).

```r
library(tximport); library(DRIMSeq); library(DEXSeq); library(stageR)

# meta: data.frame(sample_id, condition); files: named vector of quant.sf paths (names = sample_id)
# tx2gene: data.frame(tx = transcript ID, gene = gene ID)
# DRIMSeq normalizes count-column names. Normalize both inputs here so IDs such as
# `1-ctrl` remain aligned; keep the original metadata separately for display/reporting.
stopifnot(identical(names(files), meta$sample_id))
manual_sample_id <- make.names(meta$sample_id, unique = TRUE)
names(files) <- meta$sample_id <- manual_sample_id
txi <- tximport(files, type = 'salmon', txOut = TRUE, countsFromAbundance = 'no')
cts <- txi$counts
txdf <- data.frame(gene_id = tx2gene$gene[match(rownames(cts), tx2gene$tx)], feature_id = rownames(cts), cts, check.names = FALSE)
txdf <- txdf[!is.na(txdf$gene_id), ]

samples <- data.frame(sample_id = colnames(cts), condition = factor(meta$condition[match(colnames(cts), meta$sample_id)]))
stopifnot(!anyNA(samples$condition))

n <- nrow(samples); n_small <- min(table(samples$condition))
d <- dmDSdata(counts = txdf, samples = samples)
d <- dmFilter(d, min_samps_gene_expr = n, min_gene_expr = 10,
              min_samps_feature_expr = n_small, min_feature_expr = 10,
              min_samps_feature_prop = n_small, min_feature_prop = 0.1)

# Qualify with DRIMSeq:: -- once DEXSeq is attached, Biobase::samples masks DRIMSeq::samples() ("unable to find an
# inherited method for function 'samples' for signature 'object = "dmDSdata"'"); counts() is also exported by
# DEXSeq, DESeq2 and BiocGenerics, so qualify it too.
cnt <- DRIMSeq::counts(d)
dxd <- DEXSeqDataSet(
    countData = round(as.matrix(cnt[, samples$sample_id])),
    sampleData = DRIMSeq::samples(d),
    design = ~ sample + exon + condition:exon,
    featureID = cnt$feature_id,
    groupID = cnt$gene_id
)
dxd <- estimateSizeFactors(dxd)
dxd <- estimateDispersions(dxd, quiet = TRUE)
dxd <- testForDEU(dxd, reducedModel = ~ sample + exon)
dxr <- DEXSeqResults(dxd, independentFiltering = FALSE)
qval <- perGeneQValue(dxr)

pConfirmation <- matrix(dxr$pvalue, ncol = 1, dimnames = list(dxr$featureID, NULL))
pConfirmation[is.na(pConfirmation)] <- 1
tx2gene_d <- data.frame(transcript = dxr$featureID, gene = dxr$groupID)

stageRObj <- stageRTx(pScreen = qval, pConfirmation = pConfirmation, pScreenAdjusted = TRUE, tx2gene = tx2gene_d)
stageRObj <- stageWiseAdjustment(stageRObj, method = 'dtu', alpha = 0.05)

results <- getAdjustedPValues(stageRObj, order = FALSE, onlySignificantGenes = FALSE)   # geneID, txID, gene, transcript
```

Confirm IsoformSwitchAnalyzeR calls (`aSwitchList` from the workflow) with the gene-level test above; see "Count route" for why:

```r
f <- aSwitchList$isoformFeatures
hit <- unique(f$isoform_id[!is.na(f$isoform_switch_q_value) & f$isoform_switch_q_value < 0.05 & abs(f$dIF) > 0.1])
called <- unique(tx2gene$gene[match(hit, tx2gene$tx)])            # gene IDs of the tx2gene table, not gene names
confirmed <- intersect(called, names(qval)[!is.na(qval) & qval < 0.05])
```

Checked (synthetic 6 v 6, planted truth): `perGeneQValue` < 0.05 in 27 genes = 20/20 planted + 1 null-gene false positive + 5 planted 0.15-shift and 1 planted 0.06-shift genes; DRIMSeq (`dmPrecision`/`dmFit`/`dmTest`) called 27 genes, Jaccard 0.93 with DEXSeq; stageR confirmed the truly switching transcript (poison B or skip C) in 20/20. Control-versus-control split: 2 gene calls of 291. Real chrX 2 v 2 (same filter; raw / `lengthScaledTPM` / `scaledTPM` counts): DEXSeq 13 / 13 / 2 genes; DRIMSeq 11 / 11 / 7; DEXSeq and DRIMSeq overlap 7 genes with raw counts. Label-shuffled splits: see "Null check".

**Importing with tximeta instead.** `tximeta(coldata)` needs a linked transcriptome that matches the Salmon index; on Salmon 2.7.0 with tximeta 1.24.0 it did not match ("couldn't find matching transcriptome") and returned an empty `rowData`, so `rowData(se)$gene_id` / `tx_id` could **not** be run here. Simulated only: `gene_id` is a list-like column, so `unlist()` it before building the `data.frame` above. `tximeta(coldata, skipMeta = TRUE)` works (used for swish below) but carries no gene mapping.

**stageR semantics:**
- **Stage 1 (screening)**: gene-level p-value (`perGeneQValue` from DEXSeq, or DRIMSeq's gene-level adjusted p) is filtered at the desired Overall FDR.
- **Stage 2 (confirmation)**: only within significant genes, individual transcripts are tested at a within-gene FWER computed to maintain global OFDR.
- **Net effect**: gene-level FDR is properly controlled, AND the transcript that drove the call is known.
- Without stageR: naive transcript-level BH overcounts because the gene-level multiple-testing burden is ignored.
