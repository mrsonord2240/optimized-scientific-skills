'''Cluster-based compositional testing with scCODA.'''
# Reference: scCODA 0.1.9+ | Verify API if version differs
# Expects an AnnData 'adata' with obs columns 'sample' (replicate),
# 'condition' (group) and 'cell_type' (discrete clusters).
import pandas as pd
import scanpy as sc
import tensorflow as tf
from sccoda.util import cell_composition_data as dat
from sccoda.util import comp_ana as mod

adata = sc.read_h5ad('annotated.h5ad')

counts = pd.crosstab(adata.obs['sample'], adata.obs['cell_type']).reset_index()
meta = adata.obs[['sample', 'condition']].drop_duplicates()
counts = counts.merge(meta, on='sample')

data = dat.from_pandas(counts, covariate_columns=['sample', 'condition'])

# reference_cell_type must be a stable, abundant type; 'automatic' picks a
# low-dispersion type present in all samples. The reference is assumed unchanged,
# so every other effect is reported relative to it.
# sample_hmc has no seed argument and the initial state is drawn when the model is built,
# so set the TensorFlow global seed BEFORE CompositionalAnalysis().
tf.random.set_seed(42)
model = mod.CompositionalAnalysis(data, formula='condition', reference_cell_type='automatic')
# Defaults num_results=20000, num_burnin=5000 (~70 s); the argument is num_burnin, not n_burnin.
result = model.sample_hmc()

# est_fdr=0.1 sets the spike-and-slab threshold for ~10% expected FDR.
result.set_fdr(est_fdr=0.1)
result.summary()
print(result.credible_effects())
