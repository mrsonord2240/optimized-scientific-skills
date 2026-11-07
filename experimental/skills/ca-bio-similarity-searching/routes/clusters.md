# Butina clustering

```python
import json
from rdkit import Chem
from scripts.similarity_search import cluster_molecules
smiles = [x.strip() for x in open('library.smi') if x.strip()]
molecules = [Chem.MolFromSmiles(s) for s in smiles]
json.dump(cluster_molecules(molecules, cutoff=0.4), open('clusters.json', 'w'), indent=2)
```

Indices refer to input positions; inspect invalid records. Cutoff 0.4 gives centroid-neighbor similarity >=0.6, not all-pairs similarity. For MaxMin diversity selection and memory scaling see `references/clustering.md`.

Done when clusters, original order and thresholds are saved.
