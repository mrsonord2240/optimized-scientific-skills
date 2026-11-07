## ETKDGv3 (Modern Default)

ETKDGv3 (Wang et al. 2020), building on the original ETKDG method (Riniker & Landrum 2015), incorporates experimental torsion preferences and updated macrocycle handling into distance geometry.

**Goal:** Generate an ensemble of 3D conformers from a SMILES with the modern default embedding algorithm.

**Approach:** Add explicit hydrogens, configure ETKDGv3 parameters (random seed, maximum embedding iterations, random coordinates), and embed multiple conformers via `EmbedMultipleConfs`.

```python
from rdkit import Chem
from rdkit.Chem import AllChem

def gen_conformers(smiles, n_conf=20, seed=42):
    mol = Chem.MolFromSmiles(smiles)
    if mol is None or mol.GetNumAtoms() == 0:
        raise ValueError(f'Invalid SMILES: {smiles!r}')
    mol = Chem.AddHs(mol)
    params = AllChem.ETKDGv3()
    params.randomSeed = seed
    params.useRandomCoords = True
    params.maxIterations = 1000
    ids = AllChem.EmbedMultipleConfs(mol, numConfs=n_conf, params=params)
    return mol, list(ids)
```

`useRandomCoords=True` can improve convergence for macrocycles and highly flexible molecules. On an `EmbedParameters` object the documented limit is `maxIterations`; `maxAttempts` is only present in legacy positional overloads.

A runnable reference implementation of this pipeline ships in `scripts/gen_conformers.py`; there the equivalent functions are named `gen_conformer_ensemble` (this `gen_conformers`), `prune_rmsd` (this `prune_conformers_rmsd`, below), and `filter_energy_window` (this `filter_by_energy`, below) -- same operations, shorter names.

## Force-Field Optimization

After embedding, minimize each conformer to a local minimum.

**Goal:** Reduce strain in each embedded conformer to a stable local minimum and record the resulting energies.

**Approach:** Build MMFF94s force-field parameters, minimize each conformer in place, and collect energies; fall back to UFF when MMFF94 cannot parameterize the molecule.

```python
def optimize_conformers(mol, conf_ids, force_field='mmff94'):
    results = []
    if force_field == 'mmff94':
        mmff_props = AllChem.MMFFGetMoleculeProperties(mol, mmffVariant='MMFF94s')
        if mmff_props is not None:
            for cid in conf_ids:
                ff = AllChem.MMFFGetMoleculeForceField(mol, mmff_props, confId=cid)
                if ff is None:
                    results.append({'conf_id': cid, 'status': 'force_field_failed'})
                    continue
                status = ff.Minimize(maxIts=1000)
                results.append({'conf_id': cid, 'energy': ff.CalcEnergy(),
                                'converged': status == 0, 'force_field': 'MMFF94s'})
            return results
        force_field = 'uff'
    if force_field == 'uff':
        if not AllChem.UFFHasAllMoleculeParams(mol):
            raise ValueError('Neither MMFF94s nor UFF covers this molecule')
        for cid in conf_ids:
            ff = AllChem.UFFGetMoleculeForceField(mol, confId=cid)
            if ff is None:
                results.append({'conf_id': cid, 'status': 'force_field_failed'})
                continue
            status = ff.Minimize(maxIts=1000)
            results.append({'conf_id': cid, 'energy': ff.CalcEnergy(),
                            'converged': status == 0, 'force_field': 'UFF'})
        return results
    raise ValueError(f'Unsupported force field: {force_field}')
```

**MMFF94 vs MMFF94s:** MMFF94s is the static-structure variant, with modified out-of-plane and torsional terms intended to preserve planarity in selected functional groups. It is useful for geometry optimization, but it is not a universally preferred replacement for MMFF94; record which variant was used and validate it for the downstream task.

**UFF (Universal Force Field):** UFF has different and often broader parameter coverage than MMFF94, but RDKit does not guarantee coverage for every molecule and generic UFF does not validate metal coordination chemistry. Call `UFFHasAllMoleculeParams` before use and use a chemistry-appropriate method for metal complexes.

## RMSD Pruning

Remove near-duplicate conformers within a chosen RMSD cutoff to keep the ensemble diverse:

```python
import numpy as np

def prune_conformers_rmsd(mol, conf_ids, rmsd_cutoff=0.5):
    n = len(conf_ids)
    keep = []
    for i, cid in enumerate(conf_ids):
        is_unique = True
        for kept_cid in keep:
            rmsd = AllChem.GetBestRMS(mol, mol, cid, kept_cid)
            if rmsd < rmsd_cutoff:
                is_unique = False
                break
        if is_unique:
            keep.append(cid)
    return keep
```

**Typical RMSD cutoff (Source / Rationale):**

| Cutoff | Use case | Source |
|--------|----------|--------|
| 0.5 Å | Drug-like ensemble for descriptors / docking | Repository clustering heuristic; validate for the downstream task |
| 1.0 Å | Drug-like ensemble for pharmacophore | Standard ROCS / pharmacophore practice |
| 1.5-2.0 Å | Macrocycles / peptides | Repository clustering heuristic for higher conformational freedom |
| 2.0+ Å | Cluster-centroid representative ensembles | Coarse representative sampling |

## Energy Window Filtering

Remove conformers above a project-justified energy cutoff only when the energy model and downstream purpose support that choice. Bound conformers can be strained relative to an isolated-molecule minimum.

```python
def filter_by_energy(mol, conf_ids, energies, window_kcal=10.0):
    if len(conf_ids) != len(energies):
        raise ValueError('conf_ids and energies must have the same length')
    if not conf_ids:
        return []
    min_e = min(energies)
    keep = []
    for cid, e in zip(conf_ids, energies):
        if e - min_e <= window_kcal:
            keep.append(cid)
    return keep
```

Treat any numerical energy window as a starting parameter. Calibrate it against conformer recovery or downstream metric convergence and record the energy method, solvent treatment, protonation state, and temperature assumptions.
