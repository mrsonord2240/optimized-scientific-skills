# Reference: SPIA 2.50+, graphite 1.56+, clusterProfiler 4.18.4+ | Verify API if version differs
# Signed-topology perturbation analysis (the third pathway-analysis generation).
# SPIA propagates DE fold-changes through KEGG's signed signaling wiring (KGML), combining
# over-representation (pNDE) with perturbation (pPERT) into a global pG.
# NOTE: SPIA and graphite query / depend on the live KEGG release (graphite ships versioned
#       topology data, SPIA's bundled hsaSPIA is an OLDER snapshot). Pin the release before
#       publishing. SPIA is SIGNALING-ONLY: it is undefined for metabolic maps (e.g. glycolysis).

library(SPIA)
library(graphite)
library(clusterProfiler)
library(org.Hs.eg.db)

n_boot   <- 2000   # SPIA default; bootstrap replicates for the pPERT null
padj_cut <- 0.05   # DESeq2 adjusted-p gate for selecting DE genes

de <- read.csv('de_results.csv')

# SPIA needs a NAMED vector of log2 fold-changes (DE genes only) plus the universe, in Entrez space.
# bitr drops unmapped symbols and can map many-to-one, so it must be MERGED back by symbol -- never
# assigned as names directly (that recycles and silently attaches the wrong fold-change to each Entrez).
sig <- de[de$padj < padj_cut, ]
sig_map <- bitr(sig$gene, fromType = 'SYMBOL', toType = 'ENTREZID', OrgDb = org.Hs.eg.db)
de_vec  <- setNames(sig$log2FoldChange[match(sig_map$SYMBOL, sig$gene)], sig_map$ENTREZID)
de_vec  <- de_vec[!duplicated(names(de_vec))]

# universe = all measured genes (same ID space); SPIA aborts if >1% of DE IDs are absent from it
universe <- bitr(de$gene[!is.na(de$pvalue)], fromType = 'SYMBOL', toType = 'ENTREZID', OrgDb = org.Hs.eg.db)$ENTREZID

# Direct SPIA against KEGG (organism code; signaling maps only)
set.seed(123)   # SPIA's pPERT is a stochastic bootstrap; fix the seed for reproducibility
res <- spia(de = de_vec, all = universe, organism = 'hsa', nB = n_boot, plots = FALSE)
# output cols: Name, ID, pSize, NDE, pNDE, tA, pPERT, pG, pGFdr, pGFWER, Status, KEGGLINK
cat('SPIA scored', nrow(res), 'pathways;', sum(res$pGFdr < 0.05), 'significant after FDR\n')

# graphite route: harmonizes node IDs, resolves complexes/families, removes compounds,
# and runs SPIA over the cleaned graphs (also works on Reactome topology).
# runSPIA checks `datasetName(pathwaySetName) %in% dir()`, and bare dir() lists only the
# CURRENT WORKING DIRECTORY's filenames -- an absolute/tempdir() pathwaySetName can never
# match, so prepareSPIA/runSPIA must both run with a RELATIVE name from a matching setwd().
# convertIdentifiers() also prefixes graphite's node IDs ('ENTREZID:1017'), so de_vec/all
# need the same prefix or every ID join returns 0 rows even once the path bug is worked
# around. Confirmed against installed graphite 1.52.0 and current Bioconductor-release
# graphite 1.56.0 source.
db <- pathways('hsapiens', 'kegg')
db <- convertIdentifiers(db, 'ENTREZID')
de_vec_gr   <- setNames(de_vec, paste0('ENTREZID:', names(de_vec)))
universe_gr <- paste0('ENTREZID:', universe)
owd <- getwd(); setwd(tempdir())
prepareSPIA(db, 'kegg_hsa_spia')              # writes kegg_hsa_spiaSPIA.RData into tempdir()
set.seed(123)   # graphite's runSPIA bootstraps pPERT the same way spia() does
gr <- runSPIA(de = de_vec_gr, all = universe_gr, 'kegg_hsa_spia')
setwd(owd)
cat('graphite runSPIA scored', nrow(gr), 'pathways;', sum(gr$pGFdr < 0.05), 'significant after FDR\n')

write.csv(res, file.path(tempdir(), 'spia_results.csv'), row.names = FALSE)
