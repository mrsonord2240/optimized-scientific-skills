# Psix: regulated AS along trajectories

Psix tests whether observed PSI is smooth over a cell-state neighbourhood. It takes PSI and mRNA tables (events x cells) plus a cells x dimensions latent space; it builds its own neighbour metric and does not read `adata.obsp`.

```python
import pandas as pd
import scanpy as sc
import psix

adata = sc.read_h5ad('cells.h5ad')
latent = pd.DataFrame(adata.obsm['X_pca'][:, :20], index=adata.obs_names)
psix_obj = psix.Psix(psi_table='psi_matrix.tsv', mrna_table='mrna_matrix.tsv')
psix_obj.run_psix(latent=latent, n_jobs=4)
regulated = psix_obj.psix_results.query('qvals < 0.05')
```

Results are `psix_score`, `pvals`, and `qvals`, not `pvalue`. `junctions2psi(...)` can build the tables from per-cell STAR junction files but was not run for this guide; verify its current signature. Defaults (`n_random_exons=2000`, `n_neighbors=100`) can be expensive. Reducing both makes exploratory runs faster but increases null calls; tune with held-out or null controls and state the chosen values rather than quoting a universal runtime.
