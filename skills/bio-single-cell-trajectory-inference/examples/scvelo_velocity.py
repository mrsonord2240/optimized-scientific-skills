# Reference: scVelo 0.3+, scanpy 1.10+ | Verify API if version differs
#
# Compatibility (checked on scvelo 0.3.4 + numpy 2.5.3 + pandas 3.0.5): mode='dynamical'
# (via recover_dynamics) and mode='stochastic' both crash inside scvelo's own internals
# on this stack (pandas>=3 / numpy>=2 incompatibilities not fixable from the call site --
# see SKILL.md's RNA Velocity compatibility note). This example uses mode='deterministic',
# the one mode confirmed to run end-to-end here, and substitutes velocity_pseudotime for
# latent_time (which needs recover_dynamics). If your installed versions are older
# (numpy<2, pandas<3) or a newer scvelo has fixed these internals, mode='dynamical' plus
# scv.tl.latent_time is the richer, preferred path -- verify it runs before switching.
import scvelo as scv
import scanpy as sc
import numpy as np
import pandas as pd

scv.settings.verbosity = 3
scv.settings.set_figure_params('scvelo')

adata = sc.read_h5ad('adata_clustered.h5ad')
ldata = sc.read_loom('velocyto_output.loom')        # scv.read does not exist in scvelo 0.3+; use scanpy's reader
adata = scv.utils.merge(adata, ldata)

scv.pp.filter_and_normalize(adata, min_shared_counts=20)   # n_top_genes was removed from this call in scvelo 0.3+
adata.layers['normalized_X'] = adata.X.copy()               # do HVG selection as a separate step (checked on scvelo 0.3.4),
sc.pp.log1p(adata)                                           # then restore the non-log normalized X moments() expects
sc.pp.highly_variable_genes(adata, n_top_genes=2000)
adata = adata[:, adata.var['highly_variable']].copy()
adata.X = adata.layers.pop('normalized_X')
scv.pp.moments(adata, n_pcs=30, n_neighbors=30)

scv.tl.velocity(adata, mode='deterministic')
scv.tl.velocity_graph(adata, n_jobs=1, show_progress_bar=False)   # show_progress_bar=False avoids a Windows
                                                                    # multiprocessing.Manager() crash in a plain .py script

scv.pl.velocity_embedding_stream(adata, basis='umap', color='clusters',
                                  save='velocity_stream.png', dpi=150)

scv.tl.velocity_pseudotime(adata)                   # ordering proxy in place of latent_time (needs recover_dynamics)
sc.pl.umap(adata, color='velocity_pseudotime', cmap='gnuplot', save='_velocity_pseudotime.png')   # scv.pl.scatter raises KeyError: 0 on numeric .obs colors
                                                                                                        # under pandas>=3 (scvelo 0.3.4); scanpy's plot is unaffected

scv.tl.velocity_confidence(adata)
sc.pl.umap(adata, color=['velocity_confidence', 'velocity_length'], save='_velocity_confidence.png')

scv.tl.rank_velocity_genes(adata, groupby='clusters', min_corr=0.3)
velocity_genes = pd.DataFrame(adata.uns['rank_velocity_genes']['names']).head(10)   # ranks x groups (a recarray does not flatten to gene names)
# scv.pl.velocity (and scv.pl.scatter with a list of genes) raise the pandas>=3 `unique requires a Series...` TypeError in scvelo 0.3.4;
# per-gene spliced/unspliced phase portraits still work one gene at a time with explicit x/y
for gene in list(dict.fromkeys(velocity_genes.to_numpy().ravel()))[:6]:
    scv.pl.scatter(adata, gene, x='spliced', y='unspliced', color='clusters', save=f'phase_{gene}.png')

adata.write('adata_with_velocity.h5ad')
