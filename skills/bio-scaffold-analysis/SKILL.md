---
name: bio-scaffold-analysis
description: Analyzes chemical libraries by scaffold using Bemis-Murcko scaffolds, generic frameworks, cyclic skeletons, matched molecular pair (MMP) analysis via mmpdb, R-group decomposition, Free-Wilson analysis, scaffold hopping, and chemotype-aware ML train/test splits. Use when identifying chemotype clusters in a library, deriving SAR transformation rules, decomposing series into R-groups, performing scaffold-balanced QSAR splits, or planning analog campaigns.
tool_type: python
primary_tool: RDKit
license: MIT
author: GPTomics
---

## Version Compatibility

Reference examples tested with: RDKit 2024.09+, mmpdb 3.1+, scikit-learn 1.4+, datamol 0.12+.

Before using code patterns, verify installed versions match. If versions differ:
- Python: `pip show <package>` then `help(module.function)` to check signatures

If code throws ImportError, AttributeError, or TypeError, introspect the installed
package and adapt the example to match the actual API rather than retrying.

# Scaffold Analysis

Analyze chemical libraries by their underlying scaffolds. Bemis-Murcko (1996) is the canonical scaffold decomposition: ring systems + linkers, with all R-groups stripped. Generic framework + cyclic skeleton are progressively-more-abstract views. Scaffold analysis underpins QSAR train/test splits (preventing data leakage), library diversity assessment, chemotype clustering, R-group decomposition for SAR modeling, and matched molecular pair analysis (MMPA). The choice of scaffold representation determines whether two compounds are "the same series" -- a critical decision for medicinal chemistry workflows.

For reaction-based enumeration and Free-Wilson, see `chemoinformatics/reaction-enumeration`. For scaffold-hopping via fingerprints, see `chemoinformatics/similarity-searching`. For 3D shape-based scaffold hopping, see `chemoinformatics/shape-similarity`.

## Scaffold Representation Taxonomy

| Representation | Origin | Definition | Use case | Fails when |
|----------------|--------|------------|----------|------------|
| Bemis-Murcko scaffold | Bemis & Murcko 1996 | Ring systems + linkers, R-groups stripped | Default chemotype identifier | Linear molecules (no rings) -> empty scaffold |
| Generic framework | Bemis & Murcko 1996 | Bemis-Murcko with all atoms set to C, all bonds single | Topology comparison | Loses heteroatom info; may merge distinct chemotypes |
| Cyclic skeleton (CSK) | Custom RDKit transformation | Ring atoms only, all C, all single | Pure ring-topology view | Loses linker info; not a built-in Murcko option |
| Murcko atom indices | Derived by matching the scaffold to the parent | Parent-molecule atom indices | Programmatic operations | Symmetry can yield multiple equivalent matches |

```python
from rdkit import Chem
from rdkit.Chem.Scaffolds import MurckoScaffold

def all_scaffold_views(smi):
    mol = Chem.MolFromSmiles(smi)
    bm = MurckoScaffold.GetScaffoldForMol(mol)
    bm_smi = Chem.MolToSmiles(bm)

    generic = MurckoScaffold.MakeScaffoldGeneric(bm)
    generic_smi = Chem.MolToSmiles(generic)

    return {
        'bemis_murcko': bm_smi,
        'generic_framework': generic_smi,
    }
```

Example (RDKit 2026.03.6 canonical SMILES): `Cc1ccc(C(=O)NCC2CCCC2)cc1` -> Bemis-Murcko `O=C(NCC1CCCC1)c1ccccc1`; generic `CC(CCC1CCCC1)C1CCCCC1`.

The amount of collapse is library-dependent. For scale only, 1,500 ChEMBL hERG compounds in the reference run produced 857 Bemis-Murcko scaffolds but 583 generic frameworks (1.47-fold fewer groups); its largest generic framework merged four Bemis-Murcko scaffolds. Measure this ratio on the library at hand before using the generic view for a decision.

## Library Chemotype Clustering

**Goal:** Group compounds by shared Bemis-Murcko scaffold.

**Approach:** Compute scaffold for each compound; group by scaffold SMILES.

```python
from collections import defaultdict

def scaffold_clusters(smiles_list):
    clusters = defaultdict(list)
    for smi in smiles_list:
        mol = Chem.MolFromSmiles(smi)
        if mol is None:
            continue
        scaffold = MurckoScaffold.GetScaffoldForMol(mol)
        scaffold_smi = Chem.MolToSmiles(scaffold)
        clusters[scaffold_smi].append(smi)
    return clusters
```

Output: dict {scaffold_smiles: [compound_smiles, ...]}. Cluster sizes inform library diversity.

## Bemis-Murcko Scaffold Split (ML)

For QSAR / ML, random train/test split causes data leakage: compounds from the same chemotype (analogs in same series) end up in both. Bemis-Murcko split puts entire scaffolds in train or test, never both.

```python
from rdkit.Chem.Scaffolds import MurckoScaffold

def scaffold_split(df, smiles_col='smiles', train_frac=0.8, seed=42,
                   return_diagnostics=False):
    import random
    rng = random.Random(seed)

    scaffolds = defaultdict(list)
    invalid_positions = []
    for pos, smi in enumerate(df[smiles_col].tolist()):
        mol = Chem.MolFromSmiles(smi)
        if mol is None:
            invalid_positions.append(pos)
            continue
        scaff = Chem.MolToSmiles(MurckoScaffold.GetScaffoldForMol(mol))
        scaffolds[scaff].append(pos)

    if invalid_positions:
        raise ValueError(f'Invalid SMILES at row positions: {invalid_positions}')

    if not 0 < train_frac < 1:
        raise ValueError('train_frac must be strictly between 0 and 1')
    scaffold_sets = list(scaffolds.values())
    # Shuffle before size sorting so equal-size groups do not get a positional bias.
    rng.shuffle(scaffold_sets)
    scaffold_sets.sort(key=lambda x: len(x), reverse=True)

    n_total = sum(len(s) for s in scaffold_sets)
    n_train = int(n_total * train_frac)
    if len(scaffold_sets) < 2:
        raise ValueError('A scaffold split requires at least two scaffolds')
    n_test = n_total - n_train
    if n_train == 0 or n_test == 0:
        raise ValueError('train_frac produces an empty partition; use a larger dataset or a less extreme fraction')
    train_idx, test_idx = [], []
    train_group_sizes, test_group_sizes = [], []
    # Allocate large groups to the currently less-filled partition. This preserves
    # whole-scaffold isolation while spreading multi-member groups across splits;
    # singleton groups then make the requested fraction as close as possible.
    for scaff_set in scaffold_sets:
        train_fill = len(train_idx) / n_train
        test_fill = len(test_idx) / n_test
        if train_fill <= test_fill:
            train_idx.extend(scaff_set)
            train_group_sizes.append(len(scaff_set))
        else:
            test_idx.extend(scaff_set)
            test_group_sizes.append(len(scaff_set))

    # A last singleton-only adjustment can hit an exact target without moving a
    # multi-member scaffold. Otherwise report the closest whole-group split.
    if len(train_idx) != n_train:
        # Move singleton rows only; every singleton is its own scaffold group.
        train_singleton_groups = [s for s in scaffold_sets if len(s) == 1 and s[0] in train_idx]
        test_singleton_groups = [s for s in scaffold_sets if len(s) == 1 and s[0] in test_idx]
        while len(train_idx) < n_train and test_singleton_groups:
            group = test_singleton_groups.pop()
            test_idx.remove(group[0]); train_idx.extend(group)
            test_group_sizes.remove(1); train_group_sizes.append(1)
        while len(train_idx) > n_train and train_singleton_groups:
            group = train_singleton_groups.pop()
            train_idx.remove(group[0]); test_idx.extend(group)
            train_group_sizes.remove(1); test_group_sizes.append(1)

    diagnostics = {
        'requested_train_fraction': train_frac,
        'achieved_train_fraction': len(train_idx) / n_total,
        'scaffold_overlap': 0,
        'train_singleton_compound_fraction': sum(size == 1 for size in train_group_sizes) / len(train_idx),
        'test_singleton_compound_fraction': sum(size == 1 for size in test_group_sizes) / len(test_idx),
    }
    if any(len(group) > 1 for group in scaffold_sets):
        assert diagnostics['test_singleton_compound_fraction'] < 1, (
            'All test compounds are singleton scaffolds; inspect the split diagnostics')
    train, test = df.iloc[train_idx], df.iloc[test_idx]
    if return_diagnostics:
        return train, test, diagnostics
    return train, test
```

**Effect on benchmark metrics:** A scaffold split often produces different performance from a random split because it tests transfer across scaffold groups. The size and meaning of the gap are dataset- and deployment-dependent; it is not a direct universal measure of memorization.

**Caveat:** Bemis-Murcko split is *one* scaffold-split; for production ML, consider time split (newer compounds in test) or activity-cliff-balanced split.

**Class-imbalanced datasets:** The bundled splitter spreads large scaffold groups across the two partitions and returns `test_singleton_compound_fraction` with `return_diagnostics=True`; audit it with scaffold overlap and endpoint balance. Chemprop's `scaffold_balanced` split balances scaffold-group sizes; it is not label-stratified. If both group isolation and label balance are required, use a validated group-aware stratification procedure such as `StratifiedGroupKFold` where its assumptions fit, then audit every fold for scaffold overlap and endpoint balance.

## R-Group Decomposition

**Goal:** Given a defined scaffold and a set of analog compounds, extract the R-group at each numbered attachment point into a tabular SAR matrix.

```python
from rdkit.Chem import rdRGroupDecomposition as rgd

def decompose_series(compounds, scaffold_smiles_with_R):
    scaffold = Chem.MolFromSmiles(scaffold_smiles_with_R)
    if scaffold is None:
        raise ValueError('Invalid scaffold SMARTS/SMILES')
    parsed = [(i, Chem.MolFromSmiles(s)) for i, s in enumerate(compounds)]
    invalid = [i for i, mol in parsed if mol is None]
    if invalid:
        raise ValueError(f'Invalid compound SMILES at positions: {invalid}')
    mols = [mol for _, mol in parsed]
    decomp, unmatched = rgd.RGroupDecompose([scaffold], mols, asSmiles=True)
    unmatched_set = set(unmatched)
    matched_positions = [i for i in range(len(mols)) if i not in unmatched_set]
    return decomp, matched_positions, list(unmatched)

scaffold = 'c1ccc(C(=O)N[*:1])cc1-[*:2]'
compounds = ['c1ccc(C(=O)NCC)cc1F', 'c1ccc(C(=O)NCCC)cc1Cl']
table = decompose_series(compounds, scaffold)
```

Output: list of {'Core': scaffold, 'R1': r1_smiles, 'R2': r2_smiles} dicts. Used for Free-Wilson analysis (see reaction-enumeration skill).

## Matched Molecular Pair Analysis (MMPA) via mmpdb

**Goal:** Mine a SAR dataset for substructure transformations and their associated activity changes.

**Approach:** Fragment all compounds into core + variable side; index pairs differing by one transformation; report delta(activity) per transformation.

```bash
mmpdb fragment data.smi -o data.fragments
mmpdb index data.fragments -o data.mmpdb
mmpdb loadprops -p props.tsv data.mmpdb
mmpdb transform --smiles 'COc1ccccc1' --property pIC50 data.mmpdb
```

`props.tsv` must be a tab-separated file with an `ID` column matching the compound IDs in `data.smi` and a numeric `pIC50` column, for example:

```text
ID	pIC50
cmpd-001	6.4
cmpd-002	7.1
```

`mmpdb transform` emits TSV rows in tool order, not a ranked confidence table. Sort the output explicitly for the question at hand; useful reported columns are property-prefixed (for example, `pIC50_count`, `pIC50_avg`, `pIC50_std`, `pIC50_paired_t`, and `pIC50_p_value`). There is no mmpdb `confidence` column.

Interpret transformation effects from pair count, chemical-context diversity, dependence among pairs, uncertainty intervals, and prospective validation. Do not convert a universal pair-count/effect-size table into reliability labels.

## Context-Based MMPA

Classical MMPA: "Me -> F always +0.5 log units."
Context-based MMPA: "Me -> F adjacent to amide is +0.5; Me -> F adjacent to ester is -0.1."

Matched-pair effects can depend strongly on the local chemical environment, so report the transformation together with its attachment-point context rather than treating a global mean as universal (Raut & Dixit 2025). Use mmpdb's stored environments or a custom stratified analysis to compare context-specific effects.

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

## References

- Bemis GW, Murcko MA. *J. Med. Chem.* 39:2887-2893 (1996) -- original scaffold framework (DOI 10.1021/jm9602928).
- Hu, Stumpfe & Bajorath, *J. Med. Chem.* 60:1238-1246 (2017), DOI 10.1021/acs.jmedchem.6b01437 -- modern scaffold hopping review.
- Hussain J, Rea C. *J. Chem. Inf. Model.* 50:339-348 (2010) -- MMPA core method (DOI 10.1021/ci900450m).
- Raut & Dixit, *RSC Med. Chem.* 16:3281-3290 (2025), DOI 10.1039/D4MD01012D -- local-environment effects in matched molecular pairs.
- Zheng et al., *J. Cheminformatics* 13:87 (2021), DOI 10.1186/s13321-021-00565-5 -- DeepHop conditional scaffold hopping.
- Yang K et al., *J. Chem. Inf. Model.* 59:3370-3388 (2019) -- Chemprop molecular-property prediction (DOI 10.1021/acs.jcim.9b00237).
- RDKit R-group decomposition API: https://www.rdkit.org/docs/source/rdkit.Chem.rdRGroupDecomposition.html
- Chemprop splitting documentation: https://chemprop.readthedocs.io/en/main/tutorial/python/data/splitting.html

## Related Skills

- chemoinformatics/molecular-io - Parse compounds
- chemoinformatics/molecular-standardization - Standardize before scaffold extraction
- chemoinformatics/reaction-enumeration - Free-Wilson analysis on R-decomposition
- chemoinformatics/similarity-searching - 2D scaffold-hopping (FCFP4, AtomPair)
- chemoinformatics/shape-similarity - 3D scaffold-hopping
- chemoinformatics/qsar-modeling - Scaffold-aware splitting for QSAR
- chemoinformatics/generative-design - Scaffold-decoration generative tasks
