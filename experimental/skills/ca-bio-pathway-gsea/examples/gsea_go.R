# Verified 2026-09-23: R 4.5.3 / Bioconductor 3.22, clusterProfiler 4.18.4, org.Hs.eg.db 3.20.0 | Require exit status 0; verify API if version differs.
# Preranked GSEA on GO biological processes from a DESeq2-style ranked vector.
# Self-contained: builds a synthetic signed ranking over real Entrez IDs so it runs offline (no network, no DE fixture).

library(clusterProfiler)
library(org.Hs.eg.db)

set.seed(123)                                   # fixes the synthetic ranking AND the multilevel Monte Carlo

entrez_ids <- head(keys(org.Hs.eg.db, keytype = 'ENTREZID'), 4000)

# In real use this vector is de$stat (DESeq2 Wald) named by de$entrez_id: signed + variance-calibrated.
gene_list <- rnorm(length(entrez_ids))
names(gene_list) <- entrez_ids

# Plant a coordinated shift in a real GO BP term's own member genes (DNA replication,
# GO:0006260) so this demo exercises its own reporting branch instead of returning zero
# terms by construction - a pure rnorm() ranking is a null by definition and cannot enrich.
term_genes <- unique(AnnotationDbi::select(org.Hs.eg.db, keys = 'GO:0006260',
                                            keytype = 'GOALL', columns = 'ENTREZID')$ENTREZID)
planted <- intersect(names(gene_list), term_genes)
gene_list[planted] <- gene_list[planted] + 1.1   # ~1 SD shift, same order as a real effect size

gene_list <- gene_list[!is.na(gene_list)]
gene_list <- gene_list[!duplicated(names(gene_list))]   # one statistic per gene; duplicates double-count hits
gene_list <- sort(gene_list, decreasing = TRUE)         # REQUIRED: unsorted input silently mis-ranks

cat('Ranked gene list:', length(gene_list), 'genes\n')

gse_bp <- gseGO(geneList = gene_list, OrgDb = org.Hs.eg.db, keyType = 'ENTREZID',
                ont = 'BP', exponent = 1,       # p=1 weights hits by |statistic| (Subramanian 2005 default)
                minGSSize = 10, maxGSSize = 500,  # drop tiny (outlier-driven) and huge (always-enriched) sets
                eps = 0,                          # exact tiny p-values; nPerm still exists but silently downgrades to fgseaSimple - do not pass it
                pvalueCutoff = 0.05, pAdjustMethod = 'BH',
                seed = TRUE, by = 'fgsea', verbose = FALSE)

gse_bp <- setReadable(gse_bp, OrgDb = org.Hs.eg.db, keyType = 'ENTREZID')

results <- as.data.frame(gse_bp)
cat('Enriched GO terms:', nrow(results), '\n')
if (nrow(results) > 0) {
    cat('Positive NES (top of ranking):', sum(results$NES > 0), '\n')
    cat('Negative NES (bottom of ranking):', sum(results$NES < 0), '\n')
    print(head(results[, c('Description', 'NES', 'p.adjust')]))
    leading_edge <- strsplit(results$core_enrichment[1], '/')[[1]]   # the interpretable core
    cat('Leading-edge genes for top term:', length(leading_edge), '\n')
}

write.csv(results, file.path(tempdir(), 'gsea_go_results.csv'), row.names = FALSE)
