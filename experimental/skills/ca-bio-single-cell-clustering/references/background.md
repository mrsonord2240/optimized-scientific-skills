# Background

Not needed to run an analysis. This explains choices the routes and scripts already make.

## Governing principle

Resolution is not a truth knob, and clusters are not discoveries. Graph-based clustering partitions a kNN graph built in PCA space; there is no ground truth and no correct number of clusters. Clusters are hypotheses to validate by markers, stability and, where a claim is made, a significance test. Over-clustering is the default failure: any homogeneous blob can be bisected. UMAP and tSNE distances, sizes and gaps are artifacts of a neighbour-preserving objective. Clustering chooses labels to maximize between-group separation, so a marker test on those clusters tests a hypothesis built from the same data (double-dipping): the p-values are invalid, not merely inflated.

## Leiden, Louvain, SLM

| Method | Use when | Fails when |
|--------|----------|-----------|
| Leiden | Default; refinement guarantees connected communities; faster on large graphs | Backend or iteration count unpinned gives different labels across versions |
| Louvain | Legacy; Seurat's default (`algorithm=1`) | Can return internally disconnected communities |
| SLM | Seurat `algorithm=3`; tighter modularity optima | Slower; rarely needed |

Scanpy `leiden` historically defaulted to the `leidenalg` backend with a FutureWarning that the default would switch to igraph; the flip version was unconfirmed in the source, so the scripts pin `flavor='igraph', n_iterations=2, directed=False`. Seurat 5 `algorithm=4` uses `leidenbase` by default.

## Parameters

| Parameter | Range | Note |
|-----------|-------|------|
| `n_pcs` | max(elbow, about 30), commonly 30-50 | The elbow is a lower bound; do not go below about 30 because the elbow appears early |
| `n_neighbors` | 10-30 (15 default) | Secondary lever; vary only after `n_pcs` is set |
| resolution | 0.2-2.0, swept | Unidentifiable from the data; tuning it to match known cell types is confirmation bias |
| UMAP `min_dist` | 0.1-0.5 | Display only |

## Why stable does not mean real

A perfectly reproducible split can be technical: cell cycle, dissociation stress, mitochondrial fraction, ambient RNA or batch. Stability (Jaccard at least about 0.6-0.7) is necessary only. Significance is a separate question (scSHC, CHOIR). UMAP beliefs such as "far apart means more different" are folklore; PCA initialization and high perplexity make tSNE less misleading but not metric.

## Open points for audit

- Installed Scanpy is 1.12.4, newer than the source's tested 1.10+; whether `flavor='igraph'` remains the recommended pin in 1.12 was not checked against the release notes.
- The source's Seurat sweep example fed clustree directly from `RNA_snn_res.<r>`; with SCTransform the assay prefix is `SCT_snn_res.`. The script derives the prefix from the default assay.
