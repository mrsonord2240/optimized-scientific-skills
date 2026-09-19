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
scv.pl.scatter(adata, color='velocity_pseudotime', cmap='gnuplot', save='latent_time.png')

scv.tl.velocity_confidence(adata)
scv.pl.scatter(adata, color=['velocity_confidence', 'velocity_length'],
               save='velocity_confidence.png')

scv.tl.rank_velocity_genes(adata, groupby='clusters', min_corr=0.3)
velocity_genes = adata.uns['rank_velocity_genes']['names'][:10]
scv.pl.velocity(adata, var_names=list(velocity_genes.flatten()[:6]),
                basis='umap', save='top_velocity_genes.png')

adata.write('adata_with_velocity.h5ad')
