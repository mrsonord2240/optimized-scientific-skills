'''scArches scANVI label transfer with out-of-distribution gating (real-data template).

Input contract (both files in the working directory):
  reference_labeled.h5ad  raw integer counts in layers['counts']; obs['cell_type'] (reference labels),
                          obs['batch']; genes = the HVGs to use.
  query.h5ad              raw counts in layers['counts'] (same gene naming); obs['batch'] present (as in the reference).
Trains the reference scVI/scANVI itself (about 90 s on CPU for ~1.2k + 2.6k cells), maps the query by
surgery, and writes query_annotated.h5ad (predicted_label, transfer_uncertainty, X_scANVI). Seeded:
identical input gives identical output on the same CPU and package versions.
The non-obvious point: gate transferred labels on a weighted-kNN transfer uncertainty (does the cell
belong?), NOT on the classifier softmax max (which only says which label). See
ood_gating_demo.py for a synthetic, runnable version of the gating logic.
'''
# Reference: anndata 0.13, scanpy 1.12, scvi-tools 1.5, scikit-learn 1.9 | Verify API if version differs

import scvi
import scanpy as sc
import numpy as np
from sklearn.neighbors import KNeighborsClassifier

scvi.settings.seed = 0   # reproducibility: same input -> same labels (CPU, same package versions)

adata_ref = sc.read_h5ad('reference_labeled.h5ad')

# Reference scANVI head from a trained scVI model. unlabeled_category is REQUIRED.
scvi.model.SCVI.setup_anndata(adata_ref, layer='counts', batch_key='batch')
ref_vae = scvi.model.SCVI(adata_ref, n_latent=30, n_layers=2)
ref_vae.train(max_epochs=100, early_stopping=True)
ref_scanvi = scvi.model.SCANVI.from_scvi_model(ref_vae, unlabeled_category='Unknown', labels_key='cell_type')
ref_scanvi.train(max_epochs=20, n_samples_per_label=100)

adata_query = sc.read_h5ad('query.h5ad')
# Mandatory gene alignment: zero-pads missing, reorders. Silent corruption if skipped.
scvi.model.SCANVI.prepare_query_anndata(adata_query, ref_scanvi)
query_scanvi = scvi.model.SCANVI.load_query_data(adata_query, ref_scanvi)
# weight_decay=0.0 keeps the shared latent fixed so queries stay cross-comparable.
query_scanvi.train(max_epochs=100, plan_kwargs={'weight_decay': 0.0})

adata_query.obs['predicted_label'] = query_scanvi.predict()
adata_query.obsm['X_scANVI'] = query_scanvi.get_latent_representation()

# OOD gate on the shared latent. HLCA (Sikkema 2023) sets uncertainty > 0.2 to Unknown.
ref_latent = ref_scanvi.get_latent_representation()
knn = KNeighborsClassifier(n_neighbors=15, weights='distance').fit(ref_latent, adata_ref.obs['cell_type'])
uncertainty = 1.0 - knn.predict_proba(adata_query.obsm['X_scANVI']).max(axis=1)
adata_query.obs['transfer_uncertainty'] = uncertainty
adata_query.obs.loc[uncertainty > 0.2, 'predicted_label'] = 'Unknown'   # 0.2: HLCA default

print(f'Flagged Unknown (OOD / novel): {(uncertainty > 0.2).mean():.1%}')
print(adata_query.obs['predicted_label'].value_counts())

adata_query.write_h5ad('query_annotated.h5ad')
print('Wrote query_annotated.h5ad. Next: python scripts/label_marker_check.py --markers <curated.json> (the gate cannot see a missing type).')
