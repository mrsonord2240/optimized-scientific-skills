# Reference: ReactomePA 1.54+, clusterProfiler 4.18+ | Verify API if version differs
# reactome.db is LOCAL, so this runs offline. gsePathway needs a NAMED numeric vector with ENTREZ
# names, sorted DECREASING; there is no keyType argument, so map non-Entrez ids first.
library(ReactomePA)
library(clusterProfiler)
library(org.Hs.eg.db)
library(reactome.db)
library(enrichplot)

set.seed(42)

# Plant a real coordinated shift on one Reactome pathway's own member genes so this demo
# exercises its own reporting/plotting branch instead of returning 0 terms by construction -
# a pure rnorm() ranking over every gene is a null by definition and cannot enrich.
planted_id <- 'R-HSA-877300'   # Interferon gamma signaling
planted_entrez <- unique(na.omit(AnnotationDbi::select(reactome.db, keys = planted_id,
                          columns = 'ENTREZID', keytype = 'PATHID')$ENTREZID))
planted_symbols <- unique(na.omit(AnnotationDbi::select(org.Hs.eg.db, keys = planted_entrez,
                          columns = 'SYMBOL', keytype = 'ENTREZID')$SYMBOL))

# Background: a ~3000-gene sample that includes the planted pathway plus random noise genes -
# a realistic "genes measured" size, not the whole ~20,000-gene genome.
all_symbols <- keys(org.Hs.eg.db, keytype = 'SYMBOL')
noise_bg <- sample(setdiff(all_symbols, planted_symbols), 3000 - length(planted_symbols))
background_symbols <- unique(c(planted_symbols, noise_bg))

# stand-in ranked statistic: a per-gene value (t-stat / signed -log10 p / shrunken log2FC) for every gene.
de <- data.frame(symbol = background_symbols, stringsAsFactors = FALSE)
de$stat <- rnorm(nrow(de))
de$stat[de$symbol %in% planted_symbols] <- de$stat[de$symbol %in% planted_symbols] + 1.5   # planted shift

mapped <- bitr(de$symbol, fromType = 'SYMBOL', toType = 'ENTREZID', OrgDb = org.Hs.eg.db)
de <- merge(de, mapped, by.x = 'symbol', by.y = 'SYMBOL')

gene_list <- de$stat
names(gene_list) <- de$ENTREZID                 # names MUST be ENTREZ
gene_list <- sort(gene_list, decreasing = TRUE)

pvalue_cutoff <- 0.05   # filters on p.adjust (BH) by default

set.seed(123)           # gsePathway permutes; fix the seed so p-values reproduce across runs
gse <- gsePathway(geneList = gene_list, organism = 'human',
                  pvalueCutoff = pvalue_cutoff, pAdjustMethod = 'BH', verbose = FALSE,
                  BPPARAM = BiocParallel::SerialParam()) # bounded, reproducible: avoid implicit worker fan-out

results_df <- as.data.frame(gse)
cat('Rows returned:', nrow(results_df), '\n')
cat('Is planted pathway', planted_id, 'recovered?', planted_id %in% results_df$ID, '\n')
results_df

# gseaResult plots are owned by enrichment-visualization; shown here for a Reactome GSEA result.
# Route to a tempdir device so nothing lands in the working directory (ridgeplot needs ggridges).
if (nrow(results_df) > 0) {
    pdf(file.path(tempdir(), 'reactome_gsea_plots.pdf'))
    print(gseaplot2(gse, geneSetID = 1:min(3, nrow(results_df)), title = 'Reactome GSEA'))
    print(ridgeplot(gse, showCategory = 15))
    dev.off()
}

write.csv(results_df, file.path(tempdir(), 'reactome_gsea_results.csv'), row.names = FALSE)
