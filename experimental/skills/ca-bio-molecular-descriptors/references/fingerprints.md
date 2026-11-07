## Fingerprint Taxonomy

| Fingerprint | Type | Radius/Path | Bits | Use case | Fails when |
|-------------|------|-------------|------|----------|------------|
| Morgan (ECFP) | Circular | r=2 (ECFP4), r=3 (ECFP6) | 2048 typical | Drug-like similarity, ML default | Loses long-range topology; bit collisions at low nBits |
| FCFP | Functional Morgan | r=2 default | 2048 | Pharmacophore-aware similarity | Same caveats as ECFP; less specific |
| MACCS | Substructure key | 166 fixed bits | 167 | Quick fingerprint, drug-likeness | Too sparse for large diverse libraries |
| RDKit FP | Path/subgraph-based | paths and branched subgraphs up to 7 bonds by default | 2048 | RDKit-native ECFP alternative | Drug-like only; not optimal for scaffold hopping |
| AtomPair | Pair + topological distance | All atom pairs | 2048 | Long-range topological similarity | Slower than ECFP; harder to interpret |
| TopologicalTorsion | 4-atom torsion | All TT | 2048 | Path-pattern similarity | Like AP, slower than ECFP |
| Avalon | Substructure + atom pairs | Mixed | 512/1024 | Fast similarity | Less standard; older |
| MAP4 (MinHashed atom-pair) | MinHash atom-pair | r=1,2 | 1024/2048 | Biological + metabolite diversity | No supported PyPI `map4` distribution; see the MAP4 installation boundary below |
| MHFP6 (MinHash) | MinHash ECFP-like | r=3 (diam 6) | 2048 | Large-library nearest-neighbor with a compatible MinHash/LSH index | Different distance semantics from folded-bit Tanimoto |
| Pharm2D | 2D pharmacophore | feature pairs/triplets | sparse | Pharmacophore search | Sparse, slower |

**Decision:** For drug-like similarity ranking, start with **ECFP4 2048 bit** because it is fast and well characterized. MHFP6 outperformed ECFP4 for analog recovery in the benchmark reported by Probst and Reymond (2018), making it a candidate for large, diverse libraries. For scaffold hopping, benchmark ECFP4, AtomPair, TopologicalTorsion, and pharmacophore fingerprints on target-relevant actives and decoys; published comparisons do not support a universal AtomPair advantage (Gardiner et al. 2011; Riniker & Landrum 2013).

## Bit vs Count Vectors

| Form | Use | Library impact |
|------|-----|----------------|
| Bit (0/1) | Tanimoto similarity, BulkTanimotoSimilarity, RDKit fingerprint folding | Standard for similarity |
| Count (integer) | Some ML methods, RF on counts, neural fingerprints | Loses bit-level fast operations; richer signal |
| Sparse (dict) | Direct chemical interpretation (which fragments at which atoms) | Use for SHAP / atomic attribution |

```python
from rdkit import Chem
from rdkit.Chem import rdFingerprintGenerator

mol = Chem.MolFromSmiles('CCO')

morgan = rdFingerprintGenerator.GetMorganGenerator(radius=2, fpSize=2048)
ecfp4_bit = morgan.GetFingerprint(mol)
ecfp4_count = morgan.GetCountFingerprint(mol)
ecfp4_sparse = morgan.GetSparseCountFingerprint(mol)
```

## Morgan / ECFP Radius Math

ECFP-X notation: X is the **diameter** in bonds. RDKit's `radius` parameter is half of X.

| Notation | RDKit radius | Diameter | Captures |
|----------|--------------|----------|----------|
| ECFP0 | 0 | 0 | Atom identity only |
| ECFP2 | 1 | 2 | Atom + immediate neighbors |
| ECFP4 | 2 | 4 | Atom + 2-bond environment |
| ECFP6 | 3 | 6 | Atom + 3-bond environment |

**Trade-off:** Larger radius captures more specific local environments but increases collisions at fixed `nBits`. **ECFP4 2048** is a common baseline (Rogers & Hahn 2010; Wu et al. 2018). O'Boyle and Sayle (2016) showed that increasing folded-fingerprint length can improve virtual-screening performance, but they did not establish a universal 4096-bit setting or a 1-5% collision rate. Measure collision occupancy and model performance for the dataset; increase `nBits` or use an unhashed sparse representation when needed.

## FCFP vs ECFP

FCFP (Functional-Class) uses RDKit's Morgan feature invariants (donor, acceptor, aromatic, halogen, basic, and acidic) instead of atom identity. Hydrophobe is a family in `BaseFeatures.fdef`, but it is not one of the default Morgan feature-invariant classes. FCFP trades atom-specificity for functional-equivalence.

```python
ecfp_generator = rdFingerprintGenerator.GetMorganGenerator(radius=2, fpSize=2048)
feature_invariants = rdFingerprintGenerator.GetMorganFeatureAtomInvGen()
fcfp_generator = rdFingerprintGenerator.GetMorganGenerator(
    radius=2, fpSize=2048, atomInvariantsGenerator=feature_invariants)
ecfp4 = ecfp_generator.GetFingerprint(mol)
fcfp4 = fcfp_generator.GetFingerprint(mol)
```

**When to use FCFP4:** Scaffold-hopping campaigns, pharmacophore-driven similarity, cross-target activity prediction.

**When to use ECFP4:** Within-series QSAR, lead optimization, when chemotype identity matters.
