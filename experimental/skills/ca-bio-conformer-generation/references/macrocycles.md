## Macrocycle Handling

Macrocycles (>=12 atom rings) have distinct conformational issues. In current RDKit, `ETKDGv3()` enables
`useMacrocycleTorsions` by default; preserve that setting and test whether the ensemble is adequate for the
specific molecule and downstream metric. Do not claim that this flag universally improves sampling: use the
packaged default-versus-opt-out comparison below to check a representative macrocycle before escalating.

```python
from rdkit.Chem import AllChem

def macrocycle_conformers(smiles, n_conf=200, seed=42):
    mol = Chem.MolFromSmiles(smiles)
    if mol is None or mol.GetNumAtoms() == 0:
        raise ValueError(f'Invalid SMILES: {smiles!r}')
    mol = Chem.AddHs(mol)
    params = AllChem.ETKDGv3()
    params.randomSeed = seed
    params.useRandomCoords = True
    params.useMacrocycleTorsions = True
    params.useSmallRingTorsions = True
    params.maxIterations = 5000
    ids = list(AllChem.EmbedMultipleConfs(mol, numConfs=n_conf, params=params))
    if not ids:
        raise RuntimeError(f'No macrocycle conformers embedded for {smiles!r}')
    return mol, ids
```

For difficult macrocycles, CREST + GFN2-xTB is a useful higher-cost option; validate coverage against experimental or downstream evidence rather than treating one method as universally definitive.

### Macrocycle setting check

Run the packaged comparison on a representative macrocycle before treating the default ensemble as adequate:

```bash
python scripts/compare_macrocycle_embedding.py --smiles 'O=C1CCCCCCCCCCCNC1' --n-conf 40 --seed 42
```

It reports embedding counts and MMFF94s energy ranges for the ETKDGv3 default and an explicit
`useMacrocycleTorsions=False` opt-out. A difference is molecule-specific evidence, not a general benchmark;
equal results mean only that this test did not show a benefit. Increase sampling or use an orthogonal method
when the downstream metric has not converged.
