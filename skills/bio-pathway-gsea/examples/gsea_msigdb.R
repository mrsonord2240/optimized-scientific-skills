# Reference: clusterProfiler 4.18.4+, org.Hs.eg.db 3.22+, msigdbr 26+ | Verify API if version differs
# GSEA against the MSigDB Hallmark collection using the generic GSEA(TERM2GENE) function.
# Self-contained: builds a synthetic signed ranking over real Entrez IDs so it runs offline.

library(clusterProfiler)
library(org.Hs.eg.db)
library(msigdbr)

set.seed(123)

entrez_ids <- head(keys(org.Hs.eg.db, keytype = 'ENTREZID'), 4000)

# In real use this vector is de$stat (DESeq2 Wald) named by Entrez IDs.
gene_list <- rnorm(length(entrez_ids))
names(gene_list) <- entrez_ids

# msigdbr 26.x: collection= (was category=); the Entrez column is ncbi_gene (older releases used entrez_gene).
hallmarks <- msigdbr(species = 'Homo sapiens', collection = 'H')
hallmarks_t2g <- hallmarks[, c('gs_name', 'ncbi_gene')]

# Plant a coordinated shift in the real HALLMARK_OXIDATIVE_PHOSPHORYLATION member genes so
# this demo exercises its own reporting branch instead of returning zero terms by construction
# - a pure rnorm() ranking is a null by definition and cannot enrich.
oxphos_genes <- as.character(hallmarks_t2g$ncbi_gene[hallmarks_t2g$gs_name == 'HALLMARK_OXIDATIVE_PHOSPHORYLATION'])
planted <- intersect(names(gene_list), oxphos_genes)
gene_list[planted] <- gene_list[planted] + 1.1   # ~1 SD shift, same order as a real effect size

gene_list <- gene_list[!is.na(gene_list)]
gene_list <- gene_list[!duplicated(names(gene_list))]
gene_list <- sort(gene_list, decreasing = TRUE)

gse_hallmark <- GSEA(geneList = gene_list, TERM2GENE = hallmarks_t2g, exponent = 1,
                     minGSSize = 10, maxGSSize = 500, eps = 0,
                     pvalueCutoff = 0.05, pAdjustMethod = 'BH',
                     seed = TRUE, verbose = FALSE)

results <- as.data.frame(gse_hallmark)
cat('Enriched Hallmarks:', nrow(results), '\n')
if (nrow(results) > 0) {
    up_hallmarks <- results[results$NES > 0, c('ID', 'NES', 'p.adjust')]
    down_hallmarks <- results[results$NES < 0, c('ID', 'NES', 'p.adjust')]
    cat('\nTop up (positive NES):\n')
    print(head(up_hallmarks[order(-up_hallmarks$NES), ], 5))
    cat('\nTop down (negative NES):\n')
    print(head(down_hallmarks[order(down_hallmarks$NES), ], 5))
}

write.csv(results, file.path(tempdir(), 'gsea_hallmark_results.csv'), row.names = FALSE)
