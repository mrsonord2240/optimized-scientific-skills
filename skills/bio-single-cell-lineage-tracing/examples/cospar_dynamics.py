'''Clonal dynamics analysis with CoSpar'''
# Reference: cospar 0.5.0, scanpy 1.12.4, numpy 2.5.3 | Verify API if version differs
# Run in a dedicated CoSpar env (numpy>=2) -- see SKILL.md "Installation and Version
# Compatibility" for why Cassiopeia and CoSpar/scanpy cannot share one env.
import cospar as cs
import scanpy as sc

# Load AnnData with lineage information
# Requires clone_id or barcode in obs
adata = sc.read_h5ad('lineage_traced.h5ad')

# CoSpar needs a precomputed PCA and 2-D embedding (initialize_adata_object only warns when they are
# missing; the crash comes later as KeyError: 'X_emb'). Compute them on a log-normalized copy so
# adata.X stays raw counts, which CoSpar expects.
emb = adata.copy()
sc.pp.normalize_total(emb)
sc.pp.log1p(emb)
sc.pp.pca(emb)
sc.pp.neighbors(emb)
sc.tl.umap(emb)

# CoSpar expects time_info, state_info, and a clonal matrix X_clone on the AnnData
# initialize_adata_object wires those fields into the object CoSpar operates on.
# data_des must be unique per dataset: CoSpar caches similarity matrices on disk under it.
adata = cs.pp.initialize_adata_object(adata, X_clone=adata.obsm['X_clone'], time_info=adata.obs['time_info'],
                                      X_pca=emb.obsm['X_pca'], X_emb=emb.obsm['X_umap'], data_des='my_dataset')

# Infer the transition map jointly from clones at multiple timepoints and state similarity
# smooth_array applies multi-scale smoothing; results land in adata.uns['transition_map']
adata = cs.tmap.infer_Tmap_from_multitime_clones(adata, smooth_array=[15, 10, 5], sparsity_threshold=0.1)

# Fate map: probability of reaching each terminal state, propagated onto cells lacking clones
cs.tl.fate_map(adata, selected_fates=['Monocyte', 'Neutrophil'], source='transition_map')
cs.pl.fate_map(adata, selected_fates=['Monocyte', 'Neutrophil'], source='transition_map')

# Fate bias recovers the early, transcriptomically-hidden bias Weinreb 2020 showed state cannot predict
cs.tl.fate_bias(adata, selected_fates=['Monocyte', 'Neutrophil'], source='transition_map')
cs.pl.fate_bias(adata, selected_fates=['Monocyte', 'Neutrophil'], source='transition_map')
