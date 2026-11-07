## Scaffold Hopping

**Goal:** Find compounds with different scaffold but similar 3D shape / pharmacophore / activity.

| Method | Approach | Tools |
|--------|----------|-------|
| 2D similarity with FCFP4 | Functional-class fingerprint Tanimoto | similarity-searching skill |
| 3D shape (ROCS) | Tanimoto on shape + color volumes | shape-similarity skill |
| Pharmacophore | Common pharmacophore features | pharmacophore-modeling skill |
| Maximum Common Substructure (MCS) | Largest shared substructure | similarity-searching skill (rdFMCS) |
| Deep scaffold hopping | Conditional molecular generation | DeepHop (Zheng et al. 2021) |

For systematic scaffold-hop discovery, combine:
1. Find target's bioactive series
2. Compute 3D pharmacophore from bound conformer
3. ROCS / pharmacophore search against vendor catalogs
4. Filter to compounds with Bemis-Murcko scaffold NOT in training data

## Series Detection

**Goal:** Identify "analog series" within a library -- compounds sharing a scaffold + co-varying R-groups.

```python
def detect_series(smiles_list, min_size=3):
    clusters = scaffold_clusters(smiles_list)
    series = {scaff: cmpds for scaff, cmpds in clusters.items()
              if len(cmpds) >= min_size}
    return series
```

Series counts depend on library provenance, standardization, scaffold definition, and minimum size. Report the observed distribution and use series as one possible unit for SAR analysis.

## Per-Tool Failure Modes

### Bemis-Murcko -- linear molecule yields empty

**Trigger:** Compound has no rings (e.g., fatty acid, simple amine).

**Mechanism:** Bemis-Murcko strips R-groups; no rings = nothing remains.

**Symptom:** Scaffold is empty string; molecules cluster together as "no scaffold".

**Fix:** For linear-rich libraries, augment with linear chain length / functional group features.

### Bemis-Murcko -- spiro / bridged ring confusion

**Trigger:** Compound has spiro or bridged ring system.

**Mechanism:** All ring atoms included; result is the entire ring system without R-groups.

**Symptom:** Apparently different drugs share a "scaffold" because of common spiro center.

**Fix:** Validate visually; use generic framework for topology-only comparison.

### Generic framework -- loses heteroatom info

**Trigger:** Distinguishing pyridine vs benzene scaffolds.

**Mechanism:** `MakeScaffoldGeneric` sets all atoms to C.

**Symptom:** Pyridine and benzene scaffolds reported as identical.

**Fix:** Use Bemis-Murcko (heteroatoms preserved); generic framework for topology only.

### Scaffold split -- imbalanced classes

**Trigger:** Library has many singletons + few large scaffolds.

**Mechanism:** A largest-first greedy assignment can put all multi-member scaffolds in train.

**Symptom:** Test set is mostly singleton scaffolds; metrics misleading.

**Fix:** Use an allocation that spreads large scaffold groups across partitions (the bundled splitter does this and reports singleton fractions); or use scaffold-balanced cross-validation. Audit both scaffold overlap and label balance.

### MMPA -- low pair count for novel transformations

**Trigger:** Transformation rare in dataset.

**Mechanism:** Need enough pairs to estimate delta(activity).

**Symptom:** Transformation reports N=2 with very large delta.

**Fix:** Report uncertainty and context diversity, avoid overinterpreting sparse transformations, and seek additional matched evidence where appropriate.

### R-group decomposition -- ambiguous match

**Trigger:** Multiple positions in scaffold could match same R-group.

**Mechanism:** Multiple core embeddings, symmetry, and unlabeled attachment choices can yield assignments that differ from the medicinal-chemistry convention.

**Symptom:** R1/R2 columns mixed up.

**Fix:** Specify labeled attachment points, inspect the returned rows and unmatched indices, and use `RGroupDecompositionParameters` for the intended matching/alignment behavior.

## Reconciliation: Scaffold Definition Disagreements

| Concept | Definition A | Definition B | Pick which |
|---------|--------------|--------------|------------|
| Bemis-Murcko scaffold | Atoms in rings + linkers | Same | RDKit default |
| Generic framework | All C, all single bonds | All C, original bonds | `MakeScaffoldGeneric` implements the first; preserve bond orders with an explicit custom transformation |
| Cyclic skeleton | Only ring atoms | Only ring atoms, generic | Implement explicitly; it is not an RDKit Murcko flag |
| "Series" | Same Bemis-Murcko | Tanimoto > 0.8 + same MW | Bemis-Murcko for SAR; Tanimoto for screening |

For ML splits: Bemis-Murcko. For library diversity: Bemis-Murcko + cluster size. For series detection: Bemis-Murcko + R-group decomposition.

## Common Errors

| Symptom | Cause | Fix |
|---------|-------|-----|
| Murcko scaffold includes unexpected linker atoms | Bemis-Murcko linkers connect ring systems by definition | Inspect the definition; for hierarchical networks use `rdScaffoldNetwork.ScaffoldNetworkParams` with `CreateScaffoldNetwork` |
| Singleton scaffolds dominate library | Aggressive standardization | Check for tautomer-induced scaffold variation; canonicalize first |
| R-group decomposition empty | Mol doesn't match scaffold | Use FMCS to find actual shared core |
| mmpdb missing transformations | Cores too restrictive | Try smaller core requirement |
| Scaffold split test is all singletons | Largest-first greedy allocation | Use the bundled balanced allocation and inspect `test_singleton_compound_fraction`; consider scaffold-balanced cross-validation |
| Generic framework same for different drugs | Stripped heteroatom info | Use Bemis-Murcko (preserves heteroatoms) |
| MakeScaffoldGeneric error | RDKit version issue | RDKit 2024.09+ uses `Chem.Scaffolds.MurckoScaffold` |
