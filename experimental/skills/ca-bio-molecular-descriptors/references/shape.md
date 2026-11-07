## 3D Descriptors and Conformer Dependence

Conformer-dependent descriptors (asphericity, eccentricity, principal moments of inertia, RDF) require a generated 3D structure. A single conformer may be unrepresentative when the molecule is flexible; measure descriptor variation across a conformer ensemble when the downstream conclusion depends on 3D shape.

**Goal:** Compute 3D shape descriptors over a conformer ensemble rather than from a single (possibly unrepresentative) conformer.

**Approach:** Add explicit hydrogens, embed N conformers with ETKDGv3, MMFF-optimize them all, then evaluate the descriptor across each conformer for downstream averaging.

```python
from rdkit.Chem import AllChem, Descriptors3D

mol = Chem.MolFromSmiles('CCCCO')
mol = Chem.AddHs(mol)

params = AllChem.ETKDGv3()
params.randomSeed = 42
conf_ids = AllChem.EmbedMultipleConfs(mol, numConfs=20, params=params)
if not conf_ids:
    raise RuntimeError('ETKDGv3 failed to generate any conformers')
if not AllChem.MMFFHasAllMoleculeParams(mol):
    raise ValueError('MMFF94 parameters are unavailable for this molecule')
optimization_results = AllChem.MMFFOptimizeMoleculeConfs(mol, maxIters=2000)
converged_ids = [
    conf_id for conf_id, (status, _) in zip(conf_ids, optimization_results)
    if status == 0
]
if not converged_ids:
    raise RuntimeError('MMFF94 did not converge for any generated conformer')
if len(converged_ids) != len(conf_ids):
    import warnings
    warnings.warn(
        f'MMFF94 did not converge for {len(conf_ids) - len(converged_ids)} '
        'conformer(s); excluding them from the descriptor aggregate',
        RuntimeWarning,
    )

asphericities = [Descriptors3D.Asphericity(mol, confId=c) for c in converged_ids]
```

**Decision:** For QSAR / ML, choose and document the conformer count using a convergence check on representative molecules. RDKit MMFF status `1` means that the optimizer needs more iterations, not that embedding or force-field setup failed. Use a higher iteration budget, retain only converged conformers, record how many were excluded, and fail only if none converge. Report the aggregation rule, such as a simple mean or a Boltzmann-weighted average, and the energy model used for any weights.
