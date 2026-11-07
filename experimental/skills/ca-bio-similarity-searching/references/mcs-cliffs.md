## Maximum Common Substructure

**Goal:** Find the largest substructure shared across a set of molecules.

**Approach:** `rdFMCS.FindMCS` with parameters controlling atom/bond equivalence.

```python
from rdkit.Chem import rdFMCS

def mcs_smarts(mols, timeout=60, ring_match='strict', atom_match='elements'):
    params = rdFMCS.MCSParameters()
    params.Timeout = timeout
    if ring_match == 'strict':
        params.BondCompareParameters.MatchFusedRings = True
        params.BondCompareParameters.MatchFusedRingsStrict = True
        params.BondCompareParameters.RingMatchesRingOnly = True
    if atom_match == 'elements':
        params.AtomCompareParameters.MatchValences = False
    result = rdFMCS.FindMCS(mols, params)
    return result.smartsString, result.numAtoms, result.numBonds, result.canceled
```

Use cases: identify scaffold across a series, build scaffold hopping queries, generate consensus pharmacophore.

**Limit:** MCS search can become combinatorial as input count, size, and structural divergence grow. Set a finite timeout and inspect `result.canceled`: if `True`, reduce the comparison set or raise the timeout only when the additional runtime is acceptable. If it is `False` but the MCS is tiny, the supplied molecules genuinely share only a trivial fragment; pre-cluster by similarity and run MCS within a coherent cluster instead. No molecule-count or atom-count boundary guarantees tractability.

## Activity Cliff Diagnosis

**Goal:** Detect pairs of similar molecules with dissimilar activities (cliffs).

**Approach:** Compute pairwise ECFP4 Tanimoto + pIC50 delta. Flag pairs with high similarity and large activity gap.

```python
def activity_cliffs(df, sim_threshold=0.85, activity_gap=2.0, activity_col='pIC50',
                    include_chirality=True):
    generator = rdFingerprintGenerator.GetMorganGenerator(
        radius=2, fpSize=2048, includeChirality=include_chirality)
    mols = [Chem.MolFromSmiles(s) for s in df['smiles']]
    if any(mol is None for mol in mols):
        raise ValueError('activity-cliff input contains invalid SMILES')
    fps = [generator.GetFingerprint(mol) for mol in mols]
    cliffs = []
    for i in range(len(fps)):
        sims = DataStructs.BulkTanimotoSimilarity(fps[i], fps[i+1:])
        for j_off, sim in enumerate(sims):
            j = i + 1 + j_off
            if sim >= sim_threshold:
                gap = abs(df[activity_col].iloc[i] - df[activity_col].iloc[j])
                if gap >= activity_gap:
                    cliffs.append((i, j, sim, gap))
    return cliffs
```

Activity cliffs flag (a) measurement noise, (b) cryptic SAR (e.g. ring-flip changing dihedral), (c) protein conformational selection, or (d) actually informative SAR. Cliffs are an opportunity for medchem investigation, not necessarily an error. Standardize and strip salts before comparison (`chemoinformatics/molecular-standardization`); leave `include_chirality=True` when stereoisomers are relevant, and record any deliberate achiral setting.
