# Molecular fingerprints

```python
from rdkit import Chem
from scripts.calculate_descriptors import batch_fingerprints
import numpy as np
smiles = [line.strip() for line in open('input.smi') if line.strip()]
molecules = [Chem.MolFromSmiles(s) for s in smiles]
if any(m is None for m in molecules):
    raise ValueError('Invalid SMILES: correct or report inputs before fingerprinting')
np.save('fingerprints.npy', batch_fingerprints(molecules, fp_type='ecfp4', n_bits=2048))
```

The helper is `scripts/calculate_descriptors.py`. It supports ECFP4, ECFP6 and MACCS; ECFP notation is diameter and radius is half. See `references/fingerprints.md` for feature invariants, other fingerprints, bit versus count and similarity semantics.

Done when the fingerprint matrix and input order and settings are saved.
