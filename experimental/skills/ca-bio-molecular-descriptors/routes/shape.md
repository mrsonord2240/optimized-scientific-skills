# 3D shape descriptors

```python
from rdkit import Chem
from scripts.calculate_descriptors import calculate_3d_descriptors
import json
mol = Chem.MolFromSmiles(open('input.smi').readline().strip())
if mol is None:
    raise ValueError('Invalid SMILES')
json.dump(calculate_3d_descriptors(mol, random_seed=42), open('shape.json', 'w'), indent=2)
```

This helper produces one conformer and warns if optimization needs more iterations. Use an ensemble when shape-sensitive conclusions matter; `references/shape.md` shows converged-conformer aggregation and `references/errors.md` covers failures.

Done when the chosen descriptors and conformer assumptions are reported.
