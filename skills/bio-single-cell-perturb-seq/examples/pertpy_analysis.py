'''Perturb-seq analysis with pertpy: mixture guide assignment, Mixscape escaper removal, E-distance.'''
# Reference: pertpy 1.3+, scanpy 1.12+ | Verify API if version differs
# Verified end to end against pt.dt.papalexi_2021() on 2026-09-19 (pertpy 1.3.0): see fix log
# fixes/bio-single-cell-perturb-seq.md in the records repo for the run that produced these numbers.
import numpy as np
import scanpy as sc
import pertpy as pt

mdata = pt.dt.papalexi_2021()
adata = mdata.mod['rna']

sc.pp.normalize_total(adata, target_sum=1e4)
sc.pp.log1p(adata)
sc.pp.highly_variable_genes(adata, n_top_genes=2000)
sc.pp.pca(adata, n_comps=50)

# Guide assignment is a mixture problem, not a flat threshold: ambient contamination is biased to abundant guides.
# assign_mixture_model needs the optional JAX backend (pip install 'pertpy[jax]'); without it, fall back to
# assign_by_threshold, which takes a required assignment_threshold and writes a binary layer, not an obs column.
gdo = mdata.mod['gdo']
gdo.X = gdo.X.tocsr()  # ships as CSC from this loader; both assignment methods need CSR or dense
gdo.layers['counts'] = gdo.X.copy()
ga = pt.pp.GuideAssignment()
try:
    ga.assign_mixture_model(gdo, assigned_guides_key='assigned_guide')
    print(gdo.obs['assigned_guide'].value_counts().head())
except ImportError as e:
    print(f"assign_mixture_model unavailable ({e}); using assign_by_threshold instead")
    ga.assign_by_threshold(gdo, assignment_threshold=5, output_layer='assigned_guides')

# Mixscape removes non-perturbed escapers before any DE (assignment != effective perturbation).
# The perturbation/gene_target/replicate columns live on the joined mdata.obs, not on adata.obs
# (adata = mdata.mod['rna'] only carries nCount_RNA/nFeature_RNA/percent.mito) -- pull them across first.
mdata.push_obs(columns=['perturbation', 'gene_target', 'replicate'], mods=['rna'])
# Memory: the default (unbatched, sparse) path grew past 20 GB on these 20,729 cells and was killed in an audit. Densify X and
# pass batch_size: peak ~7 GB, ~2 min, KO/NP/NT counts within a few cells (checked pertpy 1.3.0). batch_size on sparse X raises ValueError.
adata.X = adata.X.toarray()
ms = pt.tl.Mixscape()
ms.perturbation_signature(adata, pert_key='perturbation', control='NT', n_neighbors=20, batch_size=1000)
ms.mixscape(adata, pert_key='gene_target', control='NT', layer='X_pert')   # pert_key renamed from labels
# An all-NP target confounds no-phenotype with no-editing: report the perturbed fraction, do not call the gene dead
print(adata.obs['mixscape_class_global'].value_counts())

# E-distance is the modern effect-size; pin the embedding and metric (sqeuclidean vs euclidean default changed)
dist = pt.tl.Distance(metric='edistance', obsm_key='X_pca')
pairwise = dist.pairwise(adata, groupby='gene_target')
print(pairwise.iloc[:5, :5])

# E-test: smallest p ~ 1/(n_perms+1), crushed by multiple testing across many perturbations.
# DistanceTest exposes no seed/random_state of its own; seed the global RNG for reproducible permutations.
np.random.seed(0)
etest = pt.tl.DistanceTest('edistance', n_perms=1000)
results = etest(adata, groupby='gene_target', contrast='NT')
print(results.sort_values('pvalue').head())
