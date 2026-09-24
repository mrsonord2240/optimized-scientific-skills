## Downstream Ligand-Activity (NicheNet)

**Goal:** Identify which sender ligand best explains the receiver's observed transcriptional response - the distinct, better-grounded question.

**Approach:** Define a receiver gene set of interest (DE genes from a condition contrast), restrict to ligands expressed in senders with receptors expressed in the receiver, and rank ligands by how well their predicted regulatory targets recover that gene set (AUPR).

```r
library(nichenetr)
library(Seurat)
library(tidyverse)

ligand_target_matrix <- readRDS('ligand_target_matrix.rds')
lr_network <- readRDS('lr_network.rds')

# Receiver gene set: garbage in -> garbage out; a noisy/batch-confounded DE list invalidates the ranking
geneset_oi <- FindMarkers(seurat_obj, ident.1 = 'activated_T', ident.2 = 'naive_T') %>%
    filter(p_val_adj < 0.05, avg_log2FC > 0.5) %>% rownames()
background <- get_expressed_genes('T_cell', seurat_obj, pct = 0.10)

expressed_ligands <- intersect(unique(lr_network$from), get_expressed_genes(c('Macrophage', 'Dendritic'), seurat_obj, 0.10))
expressed_receptors <- intersect(unique(lr_network$to), background)
potential_ligands <- lr_network %>% filter(from %in% expressed_ligands, to %in% expressed_receptors) %>% pull(from) %>% unique()

ligand_activities <- predict_ligand_activities(
    geneset = geneset_oi, background_expressed_genes = background,
    ligand_target_matrix = ligand_target_matrix, potential_ligands = potential_ligands)

# Current model ranks by aupr_corrected (AUPR is the headline metric; v1 used pearson)
best_ligands <- ligand_activities %>% top_n(30, aupr_corrected) %>% arrange(-aupr_corrected) %>% pull(test_ligand)
```
