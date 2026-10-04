'''Why an out-of-distribution signal (not a class-probability) catches a novel population FAR from every reference type.

Synthetic latent: three reference cell types and a query with a fourth population far from all of them.
A kNN class-vote probability (stand-in for a classifier softmax) gives every novel cell a confident label,
while a distance gate (mean distance to the nearest reference cells, calibrated on the reference itself;
the idea behind Symphony's Mahalanobis gate) flags them. LIMIT: this demonstrates the easy case only. On real
data a missing type that lands next to a reference type is NOT flagged by this gate or by the kNN
uncertainty gate; see SKILL.md and scripts/label_marker_check.py.
'''
# Reference: scikit-learn 1.9, anndata 0.13 | Verify API if version differs

import numpy as np
import anndata as ad
from sklearn.neighbors import KNeighborsClassifier

rng = np.random.default_rng(0)

ref_centers = {'Tcell': [0, 0], 'Bcell': [6, 0], 'Myeloid': [0, 6]}
ref_latent = np.vstack([rng.normal(c, 0.6, (300, 2)) for c in ref_centers.values()])
ref_labels = np.repeat(list(ref_centers), 300)

query_known = np.vstack([rng.normal(c, 0.6, (80, 2)) for c in ref_centers.values()])
query_novel = rng.normal([12, 12], 0.6, (80, 2))      # hepatocyte-like: belongs to nothing
query_latent = np.vstack([query_known, query_novel])
novel_mask = np.array([False] * query_known.shape[0] + [True] * query_novel.shape[0])

knn = KNeighborsClassifier(n_neighbors=15).fit(ref_latent, ref_labels)
predicted = knn.predict(query_latent)
vote_conf = knn.predict_proba(query_latent).max(axis=1)    # WRONG signal: which label, not whether

# Distance-based OOD gate: mean distance to the k nearest REFERENCE cells. Calibrate the
# threshold on the reference itself (99th percentile of its own nearest-neighbor distances).
ref_self_dist = knn.kneighbors(ref_latent)[0].mean(axis=1)
threshold = np.percentile(ref_self_dist, 99)
query_dist = knn.kneighbors(query_latent)[0].mean(axis=1)
gated = predicted.copy()
gated[query_dist > threshold] = 'Unknown'

print(f'Novel cells passing a 0.5 kNN-probability filter (WRONG): {(vote_conf[novel_mask] >= 0.5).mean():.0%}')
print(f'Novel cells caught by distance OOD gate (RIGHT):  {(gated[novel_mask] == "Unknown").mean():.0%}')
print(f'Known cells wrongly flagged Unknown:              {(gated[~novel_mask] == "Unknown").mean():.0%}')

adata_query = ad.AnnData(query_latent)
adata_query.obs['predicted_label'] = gated
adata_query.obs['ood_distance'] = query_dist
print(adata_query.obs['predicted_label'].value_counts())
