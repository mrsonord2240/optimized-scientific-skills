# Scaffold clusters

```python
import json, pandas as pd
from scripts.scaffold_split import scaffold_clusters
df = pd.read_csv('library.csv')
clusters = scaffold_clusters(df, smiles_col='smiles')
json.dump(clusters, open('clusters.json', 'w'), indent=2)
```

The helper is `scripts/scaffold_split.py`; cluster values are zero-based row positions. For generic frameworks and series definitions see `references/representations.md` and `references/interpretation.md`.

Done when cluster assignments and scaffold definitions are recorded.
