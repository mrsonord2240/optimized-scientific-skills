# MinHash fingerprints

```python
import numpy as np
from rdkit import Chem
from mhfp.encoder import MHFPEncoder
smiles = open('input.smi').readline().strip()
mol = Chem.MolFromSmiles(smiles)
if mol is None or mol.GetNumAtoms() == 0: raise ValueError('Invalid molecule')
np.save('mhfp6.npy', MHFPEncoder(2048).encode_mol(mol, radius=3))
```

For MHFP6 distance semantics see `references/minhash.md`. The source reports an mhfp/NumPy 2 incompatibility and no supported MAP4 installation in its environment; MAP4 is unsupported here until a maintained implementation is verified.

Done when the fingerprint and its representation and similarity conventions are recorded.
