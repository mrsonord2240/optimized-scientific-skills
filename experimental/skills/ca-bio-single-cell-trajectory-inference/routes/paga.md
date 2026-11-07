# PAGA and rooted diffusion pseudotime
~~~sh
python scripts/paga_dpt.py clustered.h5ad root_cell_barcode trajectory.h5ad --groupby leiden
~~~
Input needs X_pca and at least two existing clusters. Choose a root cell using independent biological or temporal evidence before running. The source uses 15 neighbors, diffusion map with 15 components and DPT with 10 diffusion components. PAGA connectivity supports topology review but no single threshold proves a continuum. Sweep choices, inspect disconnected components and compare roots.
**Done when:** output, root, graph settings, uncertainty and biological support are documented.
