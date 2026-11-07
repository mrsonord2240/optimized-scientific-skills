## Tanimoto Similarity (single query, large library)

**Goal:** Rank a library by ECFP4 Tanimoto similarity to a query molecule, returning hits above a threshold.

**Approach:** Generate ECFP4 fingerprints for all molecules once, then use `BulkTanimotoSimilarity` for O(N) lookup.

```python
from rdkit import Chem, DataStructs
from rdkit.Chem import rdFingerprintGenerator

def precompute_fps(smiles_list, radius=2, nBits=2048, include_chirality=False):
    generator = rdFingerprintGenerator.GetMorganGenerator(
        radius=radius, fpSize=nBits, includeChirality=include_chirality)
    fps = []
    for smi in smiles_list:
        mol = Chem.MolFromSmiles(smi)
        if mol is None:
            fps.append(None)
        else:
            fps.append(generator.GetFingerprint(mol))
    return fps

def search(query_smi, library_fps, threshold=0.7, include_chirality=False):
    qmol = Chem.MolFromSmiles(query_smi)
    if qmol is None:
        raise ValueError('invalid query SMILES')
    generator = rdFingerprintGenerator.GetMorganGenerator(
        radius=2, fpSize=2048, includeChirality=include_chirality)
    qfp = generator.GetFingerprint(qmol)
    valid = [(i, fp) for i, fp in enumerate(library_fps) if fp is not None]
    sims = DataStructs.BulkTanimotoSimilarity(qfp, [fp for _, fp in valid])
    return [(source_i, sim) for (source_i, _), sim in zip(valid, sims)
            if sim >= threshold]
```

## Tversky for Asymmetric Substructure-Like Search

**Goal:** Rank a library by how much each compound "contains" the features of the query (asymmetric).

**Approach:** Tversky with alpha=1, beta=0 rewards compounds containing query bits (substructure-like) while ignoring extra bits in the compound.

```python
from rdkit import DataStructs

def tversky_substructure_like(qfp, lib_fps, alpha=1.0, beta=0.0):
    return [DataStructs.TverskySimilarity(qfp, f, alpha, beta) for f in lib_fps if f]
```

Use case: identifying analogs that extend a pharmacophore vs. exact-similarity ranking. This is a feature-overlap ranker, not a substructure test: confirm every retained hit with `Chem.MolFromSmarts(query_smarts)` and `candidate.HasSubstructMatch(query_pattern)` (see `chemoinformatics/substructure-search`) before calling it a containing compound. Its precision depends on query size, fingerprint, and threshold.

## Butina Clustering

**Goal:** Group a library around Taylor-Butina centroids whose assigned neighbors are within the selected distance cutoff.

**Approach:** Compute upper-triangle distance matrix, apply Taylor-Butina with chosen distance cutoff.

```python
from rdkit.ML.Cluster import Butina

def cluster(mols, cutoff=0.4):
    generator = rdFingerprintGenerator.GetMorganGenerator(radius=2, fpSize=2048)
    fps = [generator.GetFingerprint(m) for m in mols]
    n = len(fps)
    dists = []
    for i in range(1, n):
        sims = DataStructs.BulkTanimotoSimilarity(fps[i], fps[:i])
        dists.extend([1 - s for s in sims])
    return Butina.ClusterData(dists, n, cutoff, isDistData=True)
```

`cutoff=0.4` means each assigned member was a neighbor of its selected centroid at Tanimoto >= 0.6. It does **not** guarantee that every pair of non-centroid members has Tanimoto >= 0.6. The first molecule in each returned cluster is the cluster centroid.

**Trade-off:** Butina materializes O(N^2) pairwise distances. Benchmark memory and runtime on the actual library; for much larger collections, use an approximate method such as an MHFP6 LSH forest.

## Diversity Selection (MaxMin)

**Goal:** Select N diverse compounds from a library by maximizing the minimum pairwise distance.

```python
from rdkit.SimDivFilters import rdSimDivPickers

picker = rdSimDivPickers.MaxMinPicker()
n_pick = 100
n_lib = len(fps)
selected = picker.LazyBitVectorPick(fps, n_lib, n_pick, seed=42)
```

`LazyBitVectorPick` is memory-efficient (does not materialize full distance matrix).
