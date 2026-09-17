# Reference: clusterProfiler 4.18.4+, org.Hs.eg.db 3.22+ | Verify API if version differs
# GO ORA across all three ontologies, simplified PER ONTOLOGY.
# simplify() operates on one ontology (GOSemSim similarity is defined within a single DAG),
# so an ont='ALL' object must be split into BP/MF/CC and each simplified separately.
# Self-contained, like go_enrichment_basic.R: draws foreground/universe straight from
# org.Hs.eg.db so it runs offline instead of requiring an unshipped de_results.csv. In a real
# analysis, gene_list/universe_ids instead come from a DE table's hits and tested genes (see
# "Build the Foreground and Universe from DE Results" in SKILL.md for the bitr-from-SYMBOL route).

library(clusterProfiler)
library(org.Hs.eg.db)

simplify_cut <- 0.7   # semantic-similarity redundancy cutoff; lower keeps more terms

all_entrez   <- keys(org.Hs.eg.db, keytype = 'ENTREZID')
universe_ids <- head(all_entrez, 3000)   # the genes TESTED in this assay, NOT the whole genome
gene_list    <- head(universe_ids, 200)  # foreground = the flagged hits

simplified <- lapply(c('BP', 'MF', 'CC'), function(ont) {
    ego <- enrichGO(gene = gene_list, universe = universe_ids, OrgDb = org.Hs.eg.db, keyType = 'ENTREZID',
                    ont = ont, pAdjustMethod = 'BH', pvalueCutoff = 0.05, qvalueCutoff = 0.2, readable = TRUE)
    ego <- simplify(ego, cutoff = simplify_cut, by = 'p.adjust', select_fun = min, measure = 'Wang')
    df <- as.data.frame(ego)
    cat(ont, ':', nrow(df), 'terms after simplify\n')
    df
})

combined <- do.call(rbind, simplified)
write.csv(combined, file.path(tempdir(), 'go_all_simplified.csv'), row.names = FALSE)
