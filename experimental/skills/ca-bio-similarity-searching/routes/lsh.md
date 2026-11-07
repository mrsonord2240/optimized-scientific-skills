# Approximate MHFP6 neighbors

```python
import json
from rdkit import Chem
from mhfp.encoder import MHFPEncoder
from mhfp.lsh_forest import LSHForestHelper
smiles = [x.strip() for x in open('library.smi') if x.strip()]
encoder = MHFPEncoder(2048); forest = LSHForestHelper()
fps = [encoder.encode(s, radius=3) for s in smiles]
for i, fp in enumerate(fps): forest.add(i, fp)
forest.index()
query = Chem.MolFromSmiles(open('query.smi').readline().strip())
if query is None: raise ValueError('Invalid query')
neighbors = forest.query(encoder.encode_mol(query, radius=3), k=10, data=fps)
json.dump([str(i) for i in neighbors], open('neighbors.json', 'w'))
```

The source uses approximate MHFP6 retrieval; verify its version-specific return format and benchmark recall against exact neighbors. See `references/lsh.md`.

Done when neighbors and recall/latency evaluation are recorded.
