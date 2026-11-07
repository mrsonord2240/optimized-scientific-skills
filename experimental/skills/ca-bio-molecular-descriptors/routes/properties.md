# Physicochemical properties

```python
from rdkit import Chem
from scripts.calculate_descriptors import batch_calculate_descriptors
smiles = [line.strip() for line in open('input.smi') if line.strip()]
molecules = [Chem.MolFromSmiles(s) for s in smiles]
batch_calculate_descriptors(molecules).to_csv('descriptors.csv', index=False)
```

Invalid inputs produce missing descriptor rows. The helper also offers QED and Lipinski; see `references/properties.md` for supported rules and caveats. QED and rule cutoffs are project-specific aids, not universal exclusion rules.

Done when results, missing rows and descriptor definitions are recorded.
