# Neighborhood and batch validation

A 2D embedding optimizes a neighborhood objective; it does not guarantee that every high-dimensional neighbor survives. Measure rather than assert preservation.

For each point, compute its `k` nearest neighbors in the high-dimensional input and in the 2D embedding, then report the mean overlap fraction:

```python
from sklearn.neighbors import NearestNeighbors

def neighbor_retention(high_dim, embedding, k=15):
    hi_model = NearestNeighbors(n_neighbors=k + 1).fit(high_dim)
    lo_model = NearestNeighbors(n_neighbors=k + 1).fit(embedding)
    hi = hi_model.kneighbors(high_dim, return_distance=False)[:, 1:]
    lo = lo_model.kneighbors(embedding, return_distance=False)[:, 1:]
    return sum(len(set(a) & set(b)) / k for a, b in zip(hi, lo)) / len(hi)
```

Report the value, `k`, preprocessing, and representation used. Values well below 0.5 are common and do not by themselves invalidate a display; they prevent claims that local neighborhoods were preserved wholesale.

For batch assessment, a picture is insufficient. Screen batch association across at least PC1-PC5 (for example, per-PC R2 for categorical batch indicators) and compute same-batch kNN fractions in the analysis representation. Compare against the expected fraction under mixing and stratify by biological group so composition does not masquerade as batch separation.
