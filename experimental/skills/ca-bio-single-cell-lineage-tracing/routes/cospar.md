# Clone-state transition map
~~~sh
python scripts/cospar_run.py lineage_traced.h5ad transition.h5ad unique_dataset
~~~
Input needs raw counts in X, X_clone in obsm, and time_info in obs. The script computes PCA and UMAP on a normalized copy and uses source smoothing [15,10,5] and sparsity threshold 0.1. The dataset label must be unique to avoid cache reuse. CoSpar needs multiple independent informative clones and timepoints; output is a transition model, not a phylogenetic tree. Verify version-specific API.
**Done when:** the map, clone coverage, timepoints, settings and model uncertainty are reviewed.
