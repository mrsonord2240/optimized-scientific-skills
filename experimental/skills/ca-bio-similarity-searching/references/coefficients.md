## Similarity Coefficient Taxonomy

| Coefficient | Formula | Range | Symmetric | Use case | Fails when |
|-------------|---------|-------|-----------|----------|------------|
| Tanimoto | c / (a + b - c) | 0-1 | Yes | Default for ECFP4 similarity, ranking analogs | Saturates at low similarity (drug vs natural product) |
| Dice | 2c / (a + b) | 0-1 | Yes | Bit or nonnegative sparse-count vectors when Dice semantics are intended | Thresholds depend on vector type; analog choice subjective |
| Cosine (Ochiai) | c / sqrt(a*b) | 0-1 | Yes | Count vectors, weighted similarity | Not standard for bit vectors |
| Tversky alpha,beta | c / (alpha*(a-c) + beta*(b-c) + c) | 0-1 | No when alpha != beta | Asymmetric "is A a substructure of B" queries | Parameter choice subjective; alpha=1,beta=0 = substructure-like |
| Hamming | (a + b - 2c) / nBits | 0-1 | Yes | Binary bit vectors when bit-wise disagreement matters | Does not preserve count magnitude |
| Russell-Rao | c / nBits | 0-1 | Yes | Sparse fingerprints | Biased by fingerprint density |
| Kulczynski | (c/a + c/b) / 2 | 0-1 | Yes | When fingerprints have very different bit-density | Less standard |

Where a = set bits in fp1, b = set bits in fp2, c = bits in common.

## When to Use Which Coefficient

| Scenario | Coefficient | Why |
|----------|-------------|-----|
| Standard analog search (drug-like, ECFP4) | Tanimoto, start near 0.7 | Repository starting heuristic; calibrate against project actives and analog judgments |
| Sensitive search at lower similarity | Dice, threshold 0.45 | Dice is roughly 2*Tanimoto/(1+Tanimoto); more sensitive in middle range |
| Substructure-like ranking | Tversky alpha=1, beta=0 + SMARTS confirmation | Asymmetric ranking rewards query-feature overlap; it does not prove containment |
| Count fingerprints (neural, atom-environment) | Cosine | Bit-vector Tanimoto loses information |
| Activity-cliff diagnosis | Tanimoto + property difference | Detect ECFP4>=0.85 but |delta(activity)|>=2 log units |
| Cross-target / scaffold-hopping | FCFP4 Tanimoto OR AtomPair Tanimoto | Pharmacophore-equivalent matches different scaffolds |
| Metabolomics / natural products | MHFP6 Jaccard | ECFP4 saturates near 0.2 across diverse classes |
| 3D shape | Tanimoto on shape volume overlap | See shape-similarity skill |

## Tanimoto Thresholds (Repository Starting Heuristics)

| Threshold | Interpretation | Caveat |
|-----------|----------------|--------|
| >=0.85 | Likely same scaffold + close analog | Activity cliffs still possible |
| 0.70-0.85 | Same series, R-group variation | Standard "similar" threshold |
| 0.55-0.70 | Related chemotype, different decoration | Useful for series expansion |
| 0.35-0.55 | Distant analog, possible scaffold hop | Many false positives |
| <0.35 | Mostly noise; use 3D shape or pharmacophore instead | ECFP4 not informative |

These bands are working defaults for ECFP4-like fingerprints, not transferable calibration. Inspect the target dataset's similarity distribution and known series before setting a cutoff. To make thresholds comparable across fingerprint families, calculate the actual library's pairwise-score distribution and choose the score at a retained-pair percentile (for example, the 99.5th percentile retains roughly the top 0.5% of pairs); sample pairs rather than materializing all pairs for a large library. Maggiora's similarity principle states "similar molecules tend to have similar activity" -- but **activity cliffs** (Stumpfe & Bajorath 2012) violate this. Treat high ECFP4 similarity as a prioritization signal, not evidence that activity will be preserved.

```python
import numpy as np

threshold = float(np.quantile(pairwise_similarities, 0.995))
```

## Decision Tree by Scenario

| Goal | Workflow | Tools |
|------|----------|-------|
| Find analogs of a hit (lead opt) | ECFP4 Tanimoto >=0.7 search | RDKit `BulkTanimotoSimilarity` |
| Find scaffold hops | FCFP4 OR AtomPair Tanimoto >=0.5 + filter MCS | RDKit + rdFMCS |
| Cluster library by chemotype | Butina clustering at Tanimoto 0.6 cutoff | RDKit `Butina.ClusterData` |
| Diversity sampling | MaxMin selection on Tanimoto | RDKit `rdSimDivPickers.MaxMinPicker` |
| Nearest neighbors in >1M library | LSH forest with MHFP6 | `mhfp.lsh_forest.LSHForestHelper` |
| Activity cliff diagnosis | Tanimoto + pIC50 delta scatter | Custom analysis |
| 3D similarity (shape) | USRCAT / Open3DAlign / ROCS | shape-similarity skill |
