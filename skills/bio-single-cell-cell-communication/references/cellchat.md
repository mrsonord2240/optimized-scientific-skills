## Pathway Probability (CellChat v2)

**Goal:** Summarize communication at the signaling-pathway level with sender/receiver roles.

**Approach:** Build the object, pick a database subset, identify over-expressed interactions, compute the mass-action probability with trimean, filter tiny populations, aggregate to pathways, then compute centrality for role analysis. Order matters.

```r
library(CellChat)

cellchat <- createCellChat(object = seurat_obj, group.by = 'cell_type')
cellchat@DB <- CellChatDB.human   # or CellChatDB.mouse; subsetDB(..., search='Secreted Signaling') to restrict
cellchat <- subsetData(cellchat)
cellchat <- identifyOverExpressedGenes(cellchat)
cellchat <- identifyOverExpressedInteractions(cellchat)
cellchat <- computeCommunProb(cellchat, type = 'triMean')   # trimean ~25% truncated mean: conservative
cellchat <- filterCommunication(cellchat, min.cells = 10)   # drop populations under 10 cells
cellchat <- computeCommunProbPathway(cellchat)
cellchat <- aggregateNet(cellchat)
cellchat <- netAnalysis_computeCentrality(cellchat, slot.name = 'netP')   # sender/receiver/mediator roles
# Viz: netVisual_aggregate(signaling='WNT'), netVisual_bubble(), netAnalysis_signalingRole_heatmap()
```
