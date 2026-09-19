'''Clonal dynamics analysis with CoSpar'''
# Reference: cospar 0.5.0, scanpy 1.12.4, numpy 2.5.3 | Verify API if version differs
# Run in a dedicated CoSpar env (numpy>=2) -- see SKILL.md "Installation and Version
# Compatibility" for why Cassiopeia and CoSpar/scanpy cannot share one env.
import cospar as cs
import scanpy as sc

# Load AnnData with lineage information
# Requires clone_id or barcode in obs
adata = sc.read_h5ad('lineage_traced.h5ad')

# Standard preprocessing
sc.pp.normalize_total(adata)
sc.pp.log1p(adata)
sc.pp.pca(adata)
sc.pp.neighbors(adata)
sc.tl.umap(adata)

# CoSpar expects time_info, state_info, and a clonal matrix X_clone on the AnnData
# initialize_adata_object wires those fields into the object CoSpar operates on
adata = cs.pp.initialize_adata_object(adata, X_clone=adata.obsm['X_clone'], time_info=adata.obs['time_info'])

# Infer the transition map jointly from clones at multiple timepoints and state similarity
# smooth_array applies multi-scale smoothing; results land in adata.uns['transition_map']
adata = cs.tmap.infer_Tmap_from_multitime_clones(adata, smooth_array=[15, 10, 5], sparsity_threshold=0.1)

# Fate map: probability of reaching each terminal state, propagated onto cells lacking clones
cs.tl.fate_map(adata, selected_fates=['Monocyte', 'Neutrophil'], source='transition_map')
cs.pl.fate_map(adata, selected_fates=['Monocyte', 'Neutrophil'], source='transition_map')

# Fate bias recovers the early, transcriptomically-hidden bias Weinreb 2020 showed state cannot predict
cs.tl.fate_bias(adata, selected_fates=['Monocyte', 'Neutrophil'], source='transition_map')
cs.pl.fate_bias(adata, selected_fates=['Monocyte', 'Neutrophil'], source='transition_map')
