#!/usr/bin/env Rscript
# Purpose: GO over-representation of an RNA-seq DE list corrected for gene-length bias (GOseq, Wallenius), hg38,
#   fully offline: lengths from the local TxDb, GO map from the local OrgDb.
# Usage:   Rscript goseq_ora.R de_results.csv results.csv [gene_col=gene_id] [padj=0.05] [lfc=1]
# Input:   DE table with Ensembl gene IDs (version suffixes such as .15 are stripped) and columns padj,
#          log2FoldChange, pvalue. Universe = rows with a non-NA pvalue.
# Output:  results.csv: GOseq table plus padj (BH; goseq returns unadjusted p-values); results.csv.log holds
#          the mapping rates.
# Needs:   goseq, GenomicFeatures, TxDb.Hsapiens.UCSC.hg38.knownGene, org.Hs.eg.db, clusterProfiler
# Checked: Bioconductor 3.20 pins goseq 1.58.0, GenomicFeatures 1.58.0, TxDb 3.20.0, org.Hs.eg.db 3.20.0
#          (not yet executed after the recut).
suppressPackageStartupMessages({library(clusterProfiler); library(org.Hs.eg.db); library(goseq)})

args <- commandArgs(trailingOnly = TRUE)
opts <- list(gene_col = 'gene_id', padj = '0.05', lfc = '1')
positional <- character(0)
for (a in args) {
  kv <- strsplit(a, '=', fixed = TRUE)[[1]]
  if (length(kv) >= 2 && kv[1] %in% names(opts)) opts[[kv[1]]] <- paste(kv[-1], collapse = '=')
  else if (grepl('=', a, fixed = TRUE)) stop('unknown argument ', a)
  else positional <- c(positional, a)
}
if (length(positional) != 2) stop('usage: goseq_ora.R de_results.csv results.csv [gene_col=gene_id] [padj=0.05] [lfc=1]')
out <- positional[2]

required <- c('goseq', 'GenomicFeatures', 'TxDb.Hsapiens.UCSC.hg38.knownGene', 'org.Hs.eg.db')
missing <- required[!vapply(required, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing)) stop('Install the pinned GOseq prerequisites first: ', paste(missing, collapse = ', '))

de <- read.csv(positional[1], check.names = FALSE)
for (column in c(opts$gene_col, 'padj', 'log2FoldChange', 'pvalue'))
  if (!column %in% names(de)) stop('DE table has no column "', column, '"')
ids <- sub('\\..*$', '', as.character(de[[opts$gene_col]]))
if (!any(grepl('^ENSG', ids))) stop('goseq_ora.R needs Ensembl gene IDs (ENSG...) in ', opts$gene_col)
all_genes <- unique(ids[!is.na(de$pvalue)])
hit <- !is.na(de$padj) & de$padj < as.numeric(opts$padj) & abs(de$log2FoldChange) > as.numeric(opts$lfc)
sig_genes <- unique(ids[hit])

id_map <- bitr(all_genes, fromType = 'ENSEMBL', toType = 'ENTREZID', OrgDb = org.Hs.eg.db)
id_map <- id_map[!duplicated(id_map$ENSEMBL) & !duplicated(id_map$ENTREZID), ]
conversion_rate <- nrow(id_map) / length(all_genes)
log_lines <- sprintf('ENSEMBL -> ENTREZID: %d/%d (%.1f%%)', nrow(id_map), length(all_genes), 100 * conversion_rate)
message(log_lines)
if (conversion_rate < 0.85) warning('More than 15% of tested genes failed ID conversion; resolve the annotation/build mismatch.')

de_genes <- as.integer(id_map$ENSEMBL %in% sig_genes)
names(de_genes) <- id_map$ENTREZID

txdb <- getExportedValue('TxDb.Hsapiens.UCSC.hg38.knownGene', 'TxDb.Hsapiens.UCSC.hg38.knownGene')
exons_by_gene <- GenomicFeatures::exonsBy(txdb, by = 'gene')
bias_data <- vapply(exons_by_gene, function(exons) sum(IRanges::width(GenomicRanges::reduce(exons))), numeric(1))
keep <- names(de_genes) %in% names(bias_data)
de_genes <- de_genes[keep]
bias_data <- bias_data[names(de_genes)]
if (length(de_genes) < 10L || anyNA(bias_data) || any(bias_data <= 0))
  stop('Too few valid hg38 TxDb length records after mapping; check that IDs and genome build are hg38.')

go_map <- AnnotationDbi::select(org.Hs.eg.db, keys = names(de_genes), keytype = 'ENTREZID', columns = 'GO')
go_map <- go_map[!is.na(go_map$GO), c('ENTREZID', 'GO')]
go_genes <- intersect(names(de_genes), unique(go_map$ENTREZID))
go_rate <- length(go_genes) / length(de_genes)
log_lines <- c(log_lines, sprintf('TxDb-mapped genes with GO annotations: %d/%d (%.1f%%)', length(go_genes), length(de_genes), 100 * go_rate))
message(log_lines[2])
if (length(go_genes) < 10L) stop('Too few genes have local GO annotations after mapping; check the ID namespace.')
if (go_rate < 0.85) warning('More than 15% of TxDb-mapped genes lack GO annotations; resolve the annotation mismatch before interpreting enrichment.')
de_genes <- de_genes[go_genes]
bias_data <- bias_data[names(de_genes)]
gene2cat <- split(go_map$ENTREZID, go_map$GO)

pwf <- nullp(de_genes, bias.data = bias_data, plot.fit = FALSE)
go <- goseq(pwf, gene2cat = gene2cat, method = 'Wallenius')
go$padj <- p.adjust(go$over_represented_pvalue, method = 'BH')
write.csv(go, out, row.names = FALSE)
writeLines(c(log_lines, sprintf('foreground %d of %d tested genes', sum(de_genes), length(de_genes)),
             paste('org.Hs.eg.db', as.character(packageVersion('org.Hs.eg.db')))), paste0(out, '.log'))
message('wrote ', out, ' (', nrow(go), ' terms)')
